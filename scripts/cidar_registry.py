#!/usr/bin/env python3
# coding: utf-8
"""Automatic CIDAR sequences annotation pipeline.

Downloads depositor full sequences from each CIDAR plasmid page,
rotate the plasmids so that they have the same orientation using an ubiquitous
feature as a reference position. The reference position makes sure no feature
overlaps the $0$ reference. Missing sequences (`DVK_AE` and `DVK_AF`) are
derived from `DVK_EF`.

Edits:
- Recolor all AmpR with the same color as YTK parts
- Add AmpR terminator feature with standard color

"""

import copy
import io
import itertools
import json
import re
import os
import warnings
import sys

import bs4 as bs
import six
import tqdm
import requests
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq, translate
from Bio.SeqFeature import (
    SeqFeature,
    FeatureLocation,
    CompoundLocation,
    Reference,
)
from Bio.SeqIO import read, write
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import BsaI
from fs.zipfs import ReadZipFS

from moclo.record import CircularRecord
from moclo.regex import DNARegex

ZIP_URL = "https://media.addgene.org/cms/filer_public/00/cc/00cc5b92-36b4-48ed-a49d-b84e343f149d/iverson-2015-genbank-files.zip"
URL = "https://www.addgene.org/cloning/moclo/densmore/"
UA = "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"

DVK_AF = None
DVLA_RX = None
BB_PREFIX = "gaattcgcggccgcttctagag".upper()
AMPR_TERM = "gattatcaaaaaggatctt".upper()  # Reverse 3' of AmpR terminator
KANR_PROM = "aacaccccttgtattactgtttatgtaagcagacagtttt".upper()
KANR_TERM = "gtgttacaaccaattaaccaattctga".upper()

NAME_REGEX = re.compile(r"([^ ]*) \(([^\)]*)\)(_[A-Z]{2})")
COLOR_REGEX = re.compile(r"color: (#[0-9a-fA-F]{6})")


def translate_color(feature):
    notes = feature.qualifiers.get("note", [])
    color_note = next((n for n in notes if n.startswith("color: #")), None)

    if color_note is None:
        return

    hex_color = COLOR_REGEX.match(color_note).group(1).lower()
    feature.qualifiers["note"].remove(color_note)
    feature.qualifiers.update(
        {
            "ApEinfo_fwdcolor": [hex_color],
            "ApEinfo_revcolor": [hex_color],
            "ApEinfo_graphicformat": [
                "arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0"
            ],
        }
    )


if __name__ == "__main__":

    warnings.simplefilter("ignore")
    session = requests.Session()

    # load the kit inventory page
    with session.get(URL) as res:
        soup = bs.BeautifulSoup(res.text, "html.parser")

    # load the zip archive
    with session.get(ZIP_URL) as res:
        archive = ReadZipFS(six.BytesIO(res.content)).opendir(
            "/Iverson 2015 - Addgene Final Genbank Files 15-3-30 (1)/"
        )

    # load inventory
    inventory = soup.find("table", class_="kit-inventory-table")
    it = tqdm.tqdm(inventory.find_all("tr")[1:])
    for row in it:

        # extract each row
        row_text = row.find("a").text

        # get antibiotics resistances
        resistance = row.find("span", class_="resistance-spacing").text.strip()

        # Find a name / ID
        match = NAME_REGEX.match(row_text)
        if match is not None:
            id_ = match.group(1) + match.group(3)
            name = match.group(2)
        # Find the corresponding zip sequence archive for the ID
        elif row_text.startswith(("pJ02B2Gm", "pJ02B2Rm")):
            if resistance == "Ampicillin" and not ':' in row_text:
                name = id_ = row_text + "(A)"
            else:
                name = id_ = row_text
        else:
            id_ = name = row_text

        # Update the progress bar
        it.set_description(id_)

        # extract info
        info = {
            "resistance": resistance,
            # "name": id_,
            "id": id_,
            # "type": type_,
            "location": row.find("b").text.strip().replace(" / ", ""),
            "addgene_id": row.find("a").get("href").strip("/"),
        }

        # get the ZIP sequence
        with archive.open("{}.gb".format(id_)) as f:
            gb_archive = CircularRecord(read(f, "gb"))
            pref = gb_archive.seq.find(BB_PREFIX)
            gb_archive <<= pref - 1

        # get the AddGene sequences page
        url = "https://www.addgene.org/{}/sequences/".format(info["addgene_id"])
        with session.get(url) as res:
            soup = bs.BeautifulSoup(res.text, "html.parser")

        # No full sequence available, but it's actually almost the same
        # as DVK_EF, so we can simply patch the sequence and it will be
        # fine
        if id_ in ("DVK_AE", "DVK_AF"):

            # DVK_EF
            gb_url = "https://media.addgene.org/snapgene-media/v1.5.26-0-g82d7ad5/sequences/44/95/114495/addgene-plasmid-66069-sequence-114495.gbk"
            with requests.get(gb_url) as res:
                gb = info["gb"] = CircularRecord(
                    read(io.StringIO(res.text), "gb")
                )

            if id_ == "DVK_AE":
                gb.seq = Seq(  # replace E:lacz:F with A:lacz:E
                    str(gb.seq)
                    .lower()
                    .replace("gcttagagacc", "ggagagagacc")
                    .replace("ggtctctcgct", "ggtctctgctt")
                    .upper()
                )
                gb.description = gb.description.replace(
                    "[E:LacZa:F]", "[A:LacZa:E]"
                )

            elif id_ == "DVK_AF":
                gb.seq = Seq(  # replace E:lacz:F with A:lacz:F
                    str(gb.seq)
                    .lower()
                    .replace("gcttagagacc", "ggagagagacc")
                    .upper()
                )
                gb.description = gb.description.replace(
                    "[E:LacZa:F]", "[A:LacZa:F]"
                )

        # For the other plasmids a full sequence is available so simply
        # download it from the Depositor full sequence
        else:

            # get the addgene full sequence
            section = soup.find("section", id="depositor-full")
            gb_url = section.find("a", class_="genbank-file-download").get(
                "href"
            )
            with requests.get(gb_url) as res:
                gb = info["gb"] = CircularRecord(
                    read(io.StringIO(res.text), "gb")
                )

        # Copy well documented information from one record to the other
        gb.seq, gb_archive.seq = gb.seq.upper(), gb_archive.seq.upper()
        gb.seq.alphabet = gb_archive.seq.alphabet = IUPAC.unambiguous_dna
        gb.id = gb_archive.id = id_
        gb.name = gb_archive.name = name
        gb_archive.description = gb.description
        gb.annotations["keywords"] = gb_archive.annotations["keywords"]
        gb_archive.annotations = copy.deepcopy(gb.annotations)

        # quick feature accessor
        def get_features(label):
            return (
                f for f in gb.features if label in f.qualifiers.get("label", [])
            )

        def get_features_from_note(note):
            return (
                f for f in gb.features if note in f.qualifiers.get("note", [])
            )

        # Correct overlapping features by setting the origin just before the
        # biobrick prefix
        pref = next(get_features("BioBrick prefix"))
        gb <<= pref.location.start - 1

        # Check sequences are the same, or warn the user
        if len(gb) != len(gb_archive):
            print(
                "lengths differ for", id_, ":", len(gb), "VS", len(gb_archive)
            )
        elif gb.seq != gb_archive.seq:
            print("seqs differ for", id_)

        # Fix bad position of lacz-alpha
        if id_.startswith(("DVA_", "DVK_")):
            lacz = next(get_features("lacZ-alpha"))
            start = gb.seq.lower().find("ctagagactagtgggtctca")
            end = start + 324
            lacz.location = FeatureLocation(start, end, -1)

        # Fix bad position of AmpR
        if id_.endswith("m_AF") or "A_" in id_:
            # write(gb, '/tmp/{}.gb'.format(id_), 'gb')
            ampr = next(get_features("AmpR"))
            start = gb.seq.lower().find("ttaccaatgcttaatcagtg")
            end = start + 861
            ampr.location = FeatureLocation(start, end, -1)

        # Patch duplicated fusion sites in AddGene sequences
        if "pJ02B2Gm" in id_ and "AGGTAGGT" in gb.seq:
            d_site = gb.seq.find("AGGTAGGT")
            new_gb = gb[:d_site] + gb[d_site + 4 :]
            new_gb.annotations = gb.annotations
            new_gb.id = gb.id
            new_gb.name = gb.name
            new_gb.description = gb.description
            gb = new_gb

        # Add the LuxR activator CDS
        if id_.startswith("C0062"):
            start = gb.seq.find("ATGAAAAACATAAATGCCGACGACACATACAGAATAATT")
            end = start + gb.seq[start:].translate().find("*") * 3
            loc = FeatureLocation(start, end, 1)
            lux = SeqFeature(location=loc, type="CDS")
            lux.qualifiers.update(
                {
                    "codon_start": "1",
                    "label": "LuxR repressor",
                    "product": "transcription factor LuxR",
                    "function": ["represses Lux pL promoter"],
                    "gene": "luxR",
                    "operon": "lux",
                    "db_xref": [
                        "UniProtKB/Swiss-Prot:P12746",
                        "GO:0008218",
                        "GO:0045893",
                        "GO:0009371",
                        "GO:0006351",
                        "InterPro:IPR016032",
                        "InterPro:IPR005143",
                        "InterPro:IPR036693",
                        "InterPro:IPR000792",
                        "InterPro:IPR036388",
                        "PFAM:PF03472",
                        "PFAM:PF00196",
                    ],
                    "note": "iGEM Part: BBa_R0040",
                    "translation": lux.extract(gb.seq).translate(),
                    "inference": [
                        "DESCRIPTION:alignment:blastx:UniProtKB/Swiss-Prot:P12746"
                    ],
                }
            )
            gb.features.append(lux)

        # Add missing Lux pL promoter
        elif id_.startswith("R0063"):
            start = gb.seq.lower().find(
                "cctgtacgatcctacaggtgcttatgttaagtaattgt"
            )
            plux = SeqFeature(
                location=FeatureLocation(start, start + 150, 1),
                type="promoter",
                qualifiers={
                    "label": ["Lux pL promoter"],
                    "gene": "luxR",
                    "operon": "lux",
                    "function": ["LuxR repressed weak promoter"],
                    "note": ["color: #00a1ee", "iGEM Part: BBa_R0063"],
                },
            )
            gb.features.append(plux)

        # Add pTetR promoter
        elif id_.startswith("R0040"):
            start = gb.seq.find(
                "TCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATACTGAGCAC"
            )
            ptet = SeqFeature(
                location=FeatureLocation(start, start + 54, 1),
                type="promoter",
                qualifiers={
                    "label": ["Tet promoter"],
                    "gene": "tetR",
                    "operon": "tet",
                    "function": ["TetR repressed medium-strength promoter"],
                    "note": ["color: #00a1ee", "iGEM Part: BBa_R0040"],
                },
            )
            gb.features.append(ptet)

        # Add LacI regulator
        elif id_.startswith("R0010"):
            plac = next(get_features("lac promoter"))
            start = gb.seq.find("CAATACGCAAACCGCCTCTCCCCGCG")
            plac.location = FeatureLocation(start, start + 200, 1)
            # plac = SeqFeature(
            # location=
            # type='promoter',
            plac.qualifiers.update(
                {
                    "label": "Lac regulatory sequence",
                    "gene": "lacZ",
                    "operon": "lac",
                    "function": "LacI repressed promoter",
                    "note": ["color: #00a1ee", "iGEM Part: BBa_R0010"],
                }
            )
            # gb.features.append(plac)

        # AmpR recolor and annotations
        ampr = next(get_features("AmpR"), None)
        if ampr is not None:
            ampr.qualifiers = {
                "label": "AmpR",
                "codon_start": 1,
                "gene": "bla",
                "product": "beta-lactamase",
                "function": "ampicilin and caribenicillin resistance",
                "translation": ampr.extract(gb.seq).translate(),
                "note": ["color: #9F4240"],
                "db_xref": [
                    "GO:0005515",
                    "GO:0008800",
                    "GO:0016787",
                    "GO:0030655",
                    "GO:0046677",
                    "InterPro:IPR000871",
                    "InterPro:IPR023650",
                    "InterPro:IPR012338",
                    "PDB:1ZG4",
                    "UniProtKB/Swiss-Prot:P62593",
                ],
                "EC_number": "3.5.2.6",
            }
        ampr_prom = next(get_features("AmpR promoter"), None)
        if ampr_prom is not None:
            ampr_prom.qualifiers["label"] = ["AmpR Promoter"]
            ampr_prom.qualifiers["note"] = ["color: #ff6666"]
        ampr_term_start = gb.seq.find(AMPR_TERM)
        if ampr is not None and ampr_term_start >= 0:
            ampr_term = SeqFeature(
                location=FeatureLocation(
                    ampr_term_start, ampr_term_start + 94, -1
                ),
                type="terminator",
                qualifiers={
                    "label": "AmpR Terminator",
                    "note": ["color: #ff6666"],
                },
            )
            gb.features.append(ampr_term)

        # KanR recolor and annotations
        kanr = next(get_features("KanR"), None)
        if kanr is not None:
            kanr.qualifiers.update(
                {
                    "gene": "aphA1",
                    "product": "aminoglycoside phosphotransferase",
                    "EC_number": "2.7.1.95",
                    "label": "KanR",
                    "function": "kanamicyn resistance",
                    "db_xref": [
                        "CDD:cd05150",
                        "GO:0000166",
                        "GO:0005524",
                        "GO:0008910",
                        "GO:0016301",
                        "GO:0016740",
                        "GO:0016773",
                        "GO:0016310",
                        "GO:0046677",
                        "InterPro:IPR024165",
                        "InterPro:IPR011009",
                        "InterPro:IPR002575",
                        "PFAM:PF01636",
                        "UniProtKB/Swiss-Prot:P00551",
                    ],
                    "note": ["color: #008000"],
                }
            )
        kanr_term_start = gb.seq.find(KANR_TERM)
        if kanr is not None and kanr_term_start >= 0:
            kanr_term = SeqFeature(
                location=FeatureLocation(
                    kanr_term_start, kanr_term_start + 27, -1
                ),
                type="terminator",
                qualifiers={
                    "label": ["KanR Terminator"],
                    "note": ["color: #93ff35"],
                },
            )
            gb.features.append(kanr_term)
        kanr_prom_start = gb.seq.find(KANR_PROM)
        if kanr is not None and kanr_prom_start >= 0:
            kanr_prom = SeqFeature(
                location=FeatureLocation(
                    kanr_prom_start, kanr_prom_start + 148, -1
                ),
                type="terminator",
                qualifiers={
                    "label": ["KanR Promoter"],
                    "note": ["color: #93ff35"],
                },
            )
            gb.features.append(kanr_prom)

        # GFP recolor and annotations
        gfp = next(get_features("GFP"), None)
        if gfp is not None:
            gfp.qualifiers.update(
                {
                    "label": "GFP",
                    "note": ["color: #34ff03"],
                    "product": ["green fluorescent protein"],
                    "gene": ["GFP"],
                    "db_xref": [
                        "PDB:1H6R",
                        "InterPro:IPR009017",
                        "InterPro:IPR011584",
                        "InterPro:IPR000786",
                        "PFAM:PF01353",
                        "GO:0008218",
                        "GO:0006091",
                        "GO:0018298",
                        "UniProtKB/Swiss-Prot:P42212",
                    ],
                    "inference": [
                        "DESCRIPTION:alignment:blastx:UniProtKB/Swiss-Prot:P42212"
                    ],
                }
            )

        # mRFP1 recolor and annotations
        rfp = next(get_features("mRFP1"), None)
        if rfp is not None:
            rfp.qualifiers.update(
                {
                    "label": "mRFP",
                    "product": "mRFP1",
                    "note": [
                        "monomeric derivative of DsRed (Campbell et al., 2002)",
                        "iGEM Part: BBa_E1010",
                        "color: #c16969",
                    ],
                    "db_xref": [
                        "UniProtKB/Swiss-Prot:Q9U6Y8",
                        "GO:0008218",
                        "GO:0006091",
                        "GO:0018298",
                        "PDB:2H5R",
                    ],
                }
            )

        # Copy AddGene annotations to the archive record
        for feature in gb.features:
            # get the feature sequence
            seq = feature.extract(gb.seq + gb.seq)
            if feature.location.strand == -1:
                seq = seq.reverse_complement()
            # search the feature in the
            match = DNARegex(seq).search(gb_archive)
            if match is not None:
                start, end = match.span()
                loc = FeatureLocation(start, end, feature.location.strand)
                # remove possibly duplicate annotation
                if any(f.location == loc for f in gb_archive.features):
                    other = next(
                        f for f in gb_archive.features if f.location == loc
                    )
                    gb_archive.features.remove(other)
                # add the annotation to the archive record
                new_feature = copy.deepcopy(feature)
                new_feature.location = loc
                gb_archive.features.append(new_feature)

        # quick feature accessor for amalgamated record
        def get_features(label):
            return (
                f
                for f in gb_archive.features
                if label in f.qualifiers.get("label", [])
            )

        def get_features_from_note(note):
            return (
                f
                for f in gb_archive.features
                if note in f.qualifiers.get("note", [])
            )

        # Extract the DVLA sequence to relocate bad DVLA locations
        if DVLA_RX is None and any(get_features("DVLA")):
            dvla = next(get_features("DVLA"))
            if dvla.location is not None:
                DVLA_RX = DNARegex(
                    next(get_features("DVLA")).extract(gb_archive)
                )

        # relocate and improve DVLA annotation
        dvla = next(get_features("DVLA"), None)
        if dvla is not None:
            if dvla.location is None:
                match = DVLA_RX.search(gb_archive.seq)
                dvla.location = FeatureLocation(*match.span(), strand=None)
            # TODO: set to source / synthetic

        # remove bla annotation since we have a better AmpR
        bla = next(get_features("bla"), None)
        if bla is not None:
            gb_archive.features.remove(bla)

        for f in gb_archive.features:
            if f.location is None:
                print(f)

        # sort features by start location, source always first
        gb_archive.features.sort(
            key=lambda f: (-len(gb.seq)) * (f.type == "source")
            + f.location.start
        )

        # translate color from notes to ApEinfo
        for feature in gb_archive.features:
            translate_color(feature)

        # Fix the direct submission reference
        if (
            gb_archive.annotations["references"][-1].title
            == "Direct Submission"
        ):
            ref = gb_archive.annotations["references"][-1]
        else:
            ref = Reference()
            ref.title = "Direct Submission"
            gb_archive.annotations.append(ref)
        ref.authors = "Larralde M"
        ref.journal = "Distributed with the MoClo Python library\nhttps://github.com/althonos/moclo"

        # write the final record
        dst_dir = os.path.abspath(
            os.path.join(
                __file__, "..", "..", "moclo-cidar", "registry", "cidar"
            )
        )
        dst_file = os.path.join(dst_dir, "{}.gb").format(info["id"])
        write(gb_archive, dst_file, "gb")

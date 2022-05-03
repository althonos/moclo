#!/usr/bin/env python3
# coding: utf-8
"""Automatic Icon Genetics sequences annotation pipeline.
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
from Bio.Seq import Seq, translate
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, Reference
from Bio.SeqIO import read, write
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import BsaI
from fs.zipfs import ReadZipFS

from moclo.record import CircularRecord
from moclo.regex import DNARegex

ZIP_URL = "https://media.addgene.org/cms/filer_public/b6/f8/b6f82f82-4604-4444-9886-f8577018aee4/moclo_tool_kit_genbank_files_2_1.zip"
URL = "https://www.addgene.org/cloning/moclo/marillonnet/"
UA = "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"

# Part sequence for automatic annotation / annotation relocation
BSAI = DNARegex("ggtctc")
BPII = DNARegex("gaagac")

# DVK_AF = None
# DVLA_RX = None
# BB_PREFIX = DNARegex("gaattcgcggccgcttctagag")
# SSRA_TAG = DNARegex("gctgcaaacgacgaaaactacgctttagtagct")
AMPR_TERM = "gattatcaaaaaggatctt".upper()  # Reverse 3' of AmpR terminator
KANR_PROM = "aacaccccttgtattactgtttatgtaagcagacagtttt".upper()
KANR_TERM = "gtgttacaaccaattaaccaattctga".upper()

# Some descriptions (e.g. promoter) for device / cassettes annotation
DESCS = {}

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
            "/MoClo Tool Kit genbank files 2/"
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
        id_ = name = row_text
        if name == "pAGM1311" or name == "pAGM9121":
            continue

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
        for path in ("{} cor.gbk", "{}.gbk", "{}.gb"):
            if archive.exists(path.format(id_)):
                break
        with archive.open(path.format(id_), encoding="latin-1") as f:
            rec = f.read().replace("Exported File", "Exported      ")
            gb_archive = CircularRecord(read(six.StringIO(rec), "gb"))

        # get the AddGene sequences page
        url = "https://www.addgene.org/{}/sequences/".format(info["addgene_id"])
        with session.get(url) as res:
            soup = bs.BeautifulSoup(res.text, "html.parser")

        # get the addgene full sequence
        gb_url = soup.find("a", class_="genbank-file-download").get("href")
        with requests.get(gb_url) as res:
            gb = info["gb"] = CircularRecord(read(io.StringIO(res.text), "gb"))

        if id_ == "PAGM1276":
            gb = gb.reverse_complement(True, True, True, True, True, True, True)

        # Copy well documented information from one record to the other
        gb.seq, gb_archive.seq = (gb.seq.upper(), gb_archive.seq.upper())
        gb.id = gb_archive.id = id_
        gb.name = gb_archive.name = name
        gb_archive.description = gb.description
        gb_archive.annotations = copy.deepcopy(gb.annotations)

        # quick feature accessor
        def get_features(label):
            return (f for f in gb.features if label in f.qualifiers.get("label", []))

        def get_features_from_note(note):
            return (f for f in gb.features if note in f.qualifiers.get("note", []))

        # Check sequences are the same, or skip
        if len(gb) != len(gb_archive):
            # FIXME: pAGM1276 sequences differ
            print("lengths differ for", id_, ":", len(gb), "VS", len(gb_archive))
            # gb = gb.reverse_complement(True, True, True, True, True, True, True)
            continue
        elif gb.seq != gb_archive.seq:
            print("sequences differ for", id_, ":", len(gb))
            continue

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
                    other = next(f for f in gb_archive.features if f.location == loc)
                    gb_archive.features.remove(other)
                # add the annotation to the archive record
                new_feature = copy.deepcopy(feature)
                new_feature.location = loc
                gb_archive.features.append(new_feature)

        # quick feature accessor for amalgamated record
        def get_features(label):
            return (
                f for f in gb_archive.features if label in f.qualifiers.get("label", [])
            )

        def get_features_from_note(note):
            return (
                f for f in gb_archive.features if note in f.qualifiers.get("note", [])
            )

        # AmpR recolor and annotations
        ampr = next(get_features("AmpR"), None)
        if ampr is not None:
            ampr.qualifiers = {
                "label": ["AmpR"],
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
        ampr_prom = ampr_prom or next(get_features("AmpR Promoter"), None)
        if ampr_prom is not None:
            ampr_prom.qualifiers["label"] = ["AmpR Promoter"]
            ampr_prom.qualifiers["note"] = ["color: #ff6666"]
        ampr_term_start = gb.seq.find(AMPR_TERM)
        if ampr is not None and ampr_term_start >= 0:
            ampr_term = SeqFeature(
                location=FeatureLocation(ampr_term_start, ampr_term_start + 94, -1),
                type="terminator",
                qualifiers={"label": "AmpR Terminator", "note": ["color: #ff6666"]},
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
                    "label": ["KanR"],
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
                location=FeatureLocation(kanr_term_start, kanr_term_start + 27, -1),
                type="terminator",
                qualifiers={"label": ["KanR Terminator"], "note": ["color: #93ff35"]},
            )
            gb.features.append(kanr_term)
        kanr_prom = next(get_features("KanR Promoter"), None)
        kanr_prom_start = gb.seq.find(KANR_PROM)
        if kanr is not None and kanr_prom_start >= 0:
            kanr_prom = SeqFeature(
                location=FeatureLocation(kanr_prom_start, kanr_prom_start + 148, -1),
                type="terminator",
                qualifiers={"label": ["KanR Promoter"], "note": ["color: #93ff35"]},
            )
            gb.features.append(kanr_prom)
        if kanr_prom is None and ampr_prom is not None and kanr is not None:
            kanr_prom = ampr_prom
            kanr_prom.qualifiers["label"] = ["KanR Promoter"]
            kanr_prom.qualifiers["note"] = ["color: #93ff35"]


        # SpecR recolor and annotations
        smr = next(get_features('SmR'), None)
        if smr is not None:
            smr.qualifiers.update({
                'label': ['SmR'],
                'gene': 'aadA',
                'product': 'aminoglycoside adenylyltransferase',
                'function': 'spectinomycin and streptomycin resistance',
                'note': ['color: #ffff00'],
            })

        # Remove duplicate crtI
        if len(list(get_features('crtI'))) == 2:
            x, y = get_features('crtI')
            gb_archive.features.remove(x if x.type == 'misc_feature' else y)

        # Remove duplicate crtB
        if len(list(get_features('crtB'))) == 2:
            x, y = get_features('crtB')
            gb_archive.features.remove(x if x.type == 'misc_feature' else y)

        # Remove duplicate oriV
        if len(list(get_features('oriV'))) == 2:
            gb_archive.features.remove(min(get_features('oriV'), key=lambda f: len(f.qualifiers)))

        # GFP recolor and annotations
        # gfp = next(get_features("GFP"), None)
        # if gfp is not None:
        #     gfp.qualifiers.update(
        #         {
        #             "label": "GFP",
        #             "note": ["color: #34ff03"],
        #             "product": ["green fluorescent protein"],
        #             "gene": ["GFP"],
        #             "db_xref": [
        #                 "PDB:1H6R",
        #                 "InterPro:IPR009017",
        #                 "InterPro:IPR011584",
        #                 "InterPro:IPR000786",
        #                 "PFAM:PF01353",
        #                 "GO:0008218",
        #                 "GO:0006091",
        #                 "GO:0018298",
        #                 "UniProtKB/Swiss-Prot:P42212",
        #             ],
        #             "inference": [
        #                 "DESCRIPTION:alignment:blastx:UniProtKB/Swiss-Prot:P42212"
        #             ],
        #         }
        #     )

        # mRFP1 recolor and annotations
        # rfp = next(get_features("mRFP1"), None)
        # if rfp is not None:
        #     rfp.qualifiers.update(
        #         {
        #             "label": "mRFP",
        #             "product": "mRFP1",
        #             "note": [
        #                 "monomeric derivative of DsRed (Campbell et al., 2002)",
        #                 "iGEM Part: BBa_E1010",
        #                 "color: #c16969",
        #             ],
        #             "db_xref": [
        #                 "UniProtKB/Swiss-Prot:Q9U6Y8",
        #                 "GO:0008218",
        #                 "GO:0006091",
        #                 "GO:0018298",
        #                 "PDB:2H5R",
        #             ],
        #         }
        #     )

        # remove bla annotation since we have a better AmpR
        # bla = next(get_features("bla"), None)
        # if bla is not None:
        #     gb_archive.features.remove(bla)

        # Remove bad Sm/Sp
        rm_feats = [
            "Sm/Sp",
            r"Sm/Sp\no\DraIII",
            r"spec orf?",
            r"spec\orf?",
            r"Kan\(no\BpiI)",
            r"spec",
            r"rep\(pMB1)",
            r'rep - pMB1',
            "NPTII",
            "AP(R)",
            r"AP\r",
            "ALPHA",
            "Kan",
        ]
        for label in rm_feats:
            for f in get_features(label):
                gb_archive.features.remove(f)

        # FIXME
        # have all the plasmids in the same direction, i.e. so that the
        # antibiotics resistance cassette is always on the reverse strand
        # antibio = next(x for y in ('AmpR', 'SmR', 'KanR') for x in get_features(y))
        # if antibio.location.strand != -1:
        #     gb_archive = gb_archive.reverse_complement(True, True, True, True, True, True, True,)

        # Remove all "/vntifkey" feature qualifier
        for feature in gb_archive.features:
            feature.qualifiers.pop("vntifkey", None)

        # sort features by start location, source always first
        gb_archive.features.sort(
            key=lambda f: (-len(gb.seq)) * (f.type == "source") + f.location.start
        )

        # translate color from notes to ApEinfo
        for feature in gb_archive.features:
            translate_color(feature)

        # Fix the direct submission reference
        if gb_archive.annotations["references"][-1].title == "Direct Submission":
            ref = gb_archive.annotations["references"][-1]
        else:
            ref = Reference()
            ref.title = "Direct Submission"
            gb_archive.annotations.append(ref)
        ref.authors = "Larralde M"
        ref.journal = "Distributed with the MoClo Python library\nhttps://github.com/althonos/moclo"

        # write the final record
        dst_dir = os.path.abspath(
            os.path.join(__file__, "..", "..", "moclo-moclo", "registry", "moclo")
        )
        dst_file = os.path.join(dst_dir, "{}.gb").format(info["id"])
        write(gb_archive, dst_file, "gb")

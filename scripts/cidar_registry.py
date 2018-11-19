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
import re
import os
import warnings

import bs4 as bs
import six
import tqdm
import requests
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, Reference
from Bio.SeqIO import read, write
from fs.zipfs import ReadZipFS

from moclo.record import CircularRecord
from moclo.regex import DNARegex

from features import annotate, translate_color


ZIP_URL = "https://media.addgene.org/cms/filer_public/00/cc/00cc5b92-36b4-48ed-a49d-b84e343f149d/iverson-2015-genbank-files.zip"
URL = "https://www.addgene.org/cloning/moclo/densmore/"
UA = "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"

# Part sequence for automatic annotation / annotation relocation
DVK_AF = None
DVLA_RX = None
BB_PREFIX = DNARegex("gaattcgcggccgcttctagag")
SSRA_TAG = DNARegex("gctgcaaacgacgaaaactacgctttagtagct")
AMPR_TERM = "gattatcaaaaaggatctt".upper()  # Reverse 3' of AmpR terminator
KANR_PROM = "aacaccccttgtattactgtttatgtaagcagacagtttt".upper()
KANR_TERM = "gtgttacaaccaattaaccaattctga".upper()

# Some descriptions (e.g. promoter) for device / cassettes annotation
DESCS = {}

# Regexes
NAME_REGEX = re.compile(r"([^ ]*) \(([^\)]*)\)(_[A-Z]{2})")


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
            if resistance == "Ampicillin" and not ":" in row_text:
                name = id_ = row_text + "(A)"
            else:
                name = id_ = row_text
        else:
            id_ = name = row_text

        # DVA and DVK are empty vectors, they are not valid MoClo elements
        if id_ in ('DVA', 'DVK'):
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
        with archive.open("{}.gb".format(id_)) as f:
            gb_archive = CircularRecord(read(f, "gb"))
            pref, _ = BB_PREFIX.search(gb_archive).span()
            gb_archive <<= pref - 1

        # get the AddGene sequences page
        url = "https://www.addgene.org/{}/sequences/".format(info["addgene_id"])
        with session.get(url) as res:
            soup = bs.BeautifulSoup(res.text, "html.parser")

        # No full sequence available, but it's actually almost the same
        # as DVK_EF, so we can simply patch the sequence and it will be
        # fine
        if id_ in ("DVK_AE", "DVK_AF"):

            # Load the DVK_EF page and find the GenBank file link
            with requests.get("https://www.addgene.org/66069/sequences/") as res:
                soup = bs.BeautifulSoup(res.text, "html.parser")
                section = soup.find("section", id="depositor-full")
                gb_url = soup.find("a", class_="genbank-file-download").get('href')

            # DVK_EF
            with requests.get(gb_url) as res:
                gb = info["gb"] = CircularRecord(read(io.StringIO(res.text), "gb"))

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
            gb_url = section.find("a", class_="genbank-file-download").get("href")
            with requests.get(gb_url) as res:
                gb = info["gb"] = CircularRecord(read(io.StringIO(res.text), "gb"))

        # Copy well documented information from one record to the other
        gb.seq, gb_archive.seq = (gb.seq.upper(), gb_archive.seq.upper())
        gb.seq.alphabet = gb_archive.seq.alphabet = IUPAC.unambiguous_dna
        gb.id = gb_archive.id = id_
        gb.name = gb_archive.name = name
        gb_archive.description = gb.description
        gb.annotations["keywords"] = gb_archive.annotations["keywords"]
        gb_archive.annotations = copy.deepcopy(gb.annotations)

        # quick feature accessor
        def get_features(label):
            return (f for f in gb.features if label in f.qualifiers.get("label", []))

        def get_features_from_note(note):
            return (f for f in gb.features if note in f.qualifiers.get("note", []))

        # Correct overlapping features by setting the origin just before the
        # biobrick prefix
        pref = next(get_features("BioBrick prefix"))
        gb <<= pref.location.start - 1

        # Check sequences are the same, or warn the user
        # if len(gb) != len(gb_archive):
        #     print(
        #         "lengths differ for", id_, ":", len(gb), "VS", len(gb_archive)
        #     )
        # elif gb.seq != gb_archive.seq:
        #     print("seqs differ for", id_)

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
            new_gb = gb[:d_site] + gb[d_site + 4:]
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
            annotate("luxr", lux, gb.seq)
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
            annotate("ampr", ampr, gb.seq)

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
            annotate("kanr", kanr, gb.seq)

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
            annotate("gfp", gfp, gb.seq)

        # mRFP1 recolor and annotations
        rfp = next(get_features("mRFP1"), None)
        if rfp is not None:
            annotate("mrfp", rfp, gb.seq)

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

        # improve B0015 annotation
        b0015 = next(get_features('B0015'), None)
        if b0015 is not None:
            b0015.type = 'terminator'
            b0015.qualifiers = {
                'label': 'B0015 Double Terminator',
                'note': ['color: #ff8eff', '"iGEM Part: BBa_B0015"']
            }

        # merge Lux pL promoter
        r0063 = next(get_features('R0063'), None)
        if r0063 is not None:
            luxpl = next(get_features('Lux pL promoter'))
            luxpl.location = r0063.location
            gb_archive.features.remove(r0063)

        # add LVA ssrA tag
        ssra_match = SSRA_TAG.search(gb_archive.seq)
        if ssra_match is not None:
            ssra = SeqFeature(type="CDS")
            ssra.location = FeatureLocation(*ssra_match.span(), strand=1)
            ssra.qualifiers = {
                "label": ["ssrA tag (LVA)"],
                "product": [
                    "C-terminal peptide that mediates degradation in bacteria through the ClpXP and ClpAP proteases (McGinness et al., 2006)"
                ],
                "translation": "AANDENYALVA",
                "note": [
                    "mutant LVA variant that confers accelerated degradation under some conditions (Andersen et al., 1998)",
                    "color: #cc99b2",
                ],
            }
            gb_archive.features.append(ssra)

        # Replace E0040m with well annotated GFP
        e0040m = next(get_features("E0040m"), None)
        if e0040m is not None:
            if any(get_features("GFP")):
                gb_archive.features.remove(next(get_features("GFP")))
            e0040m.qualifiers.update(gfp.qualifiers)

        # Replace E1010m with well annotated mRFP
        e1010m = next(get_features("E1010m"), None)
        if e1010m is not None:
            gb_archive.features.remove(e1010m)

        # Replace E0030m with well annotated YFP (TODO)
        e0030 = next(get_features("E0030"), None)
        if e0030 is not None:
            gb_archive.features.remove(e0030)

        # Replace C0080 with well annotated araC (TODO)
        c0080 = next(get_features("C0080"), None)
        if c0080 is not None:
            gb_archive.features.remove(c0080)

        # Replace C0040 with well annotated TetR (TODO)
        c0040 = next(get_features("C0040"), None)
        if c0040 is not None:
            gb_archive.features.remove(c0040)

        # Replace C0012 with well annotated LacI (TODO)
        c0012 = next(get_features("C0012"), None)
        if c0012 is not None:
            gb_archive.features.remove(c0012)

        # Replace C0062 with well annotated LuxR repressor (TODO)
        # and fix LuxR repressor CDS to include STOP codon
        c0062 = next(get_features('C0062'), None)
        if c0062 is not None:
            gb_archive.features.remove(c0062)
            luxr_rep = next(get_features('LuxR repressor'))
            luxr_rep.location = FeatureLocation(
                luxr_rep.location.start,
                luxr_rep.location.end + 3,
                luxr_rep.location.strand
            )

        # Remove duplicate LacO
        lacO = next(get_features("LacO"), None)
        if lacO is not None:
            gb_archive.features.remove(lacO)

        # Annotate promoters
        promoters = [
            "J23102",
            "J23100",
            "J23103",
            "J23106",
            "J23107",
            "J23116",
            "I13453",
        ]
        j231 = next(iter(itertools.chain(*map(get_features, promoters))), None)
        if j231 is not None:
            name = j231.qualifiers["label"]
            if isinstance(name, list):
                name = name.pop()
            if name not in DESCS:
                DESCS[name] = re.search(r": (.*) \[", gb_archive.description).group(1)
            j231.type = "promoter"
            j231.qualifiers = {
                "label": "{} promoter".format(name),
                "note": [DESCS[name], "color: #00a1ee"],
            }

        # Annotate RBS
        rbses = ["B0032m", "B0033m", "B0034m", "BCD2", "BCD8", "BCD12"]
        rbs = next(iter(itertools.chain(*map(get_features, promoters))), None)
        if rbs is not None:
            name = rbs.qualifiers["label"]
            if isinstance(name, list):
                name = name.pop()
            if name not in DESCS:
                DESCS[name] = re.search(
                    r": (.*) \[", gb_archive.description
                ).group(1)
            rbs.type = "RBS"
            rbs.qualifiers["note"] = [DESCS[name], "color: #f58a5e"]

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

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
import fs
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

ZIP_URL = "https://media.addgene.org/cms/filer_public/39/6d/396d2c07-f428-4658-b723-e1fb765b2df5/plant_parts_genbank_updated_mar2015.zip"
URL = "https://www.addgene.org/kits/patron-moclo/"
UA = "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"

# Part sequence for automatic annotation / annotation relocation
BSAI = DNARegex("ggtctc")
BPII = DNARegex("gaagac")
BSAI_R = DNARegex("gagacc")

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
            "Plant Parts genbank updated Mar2015"
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

        if id_ == "pICSL30008":
            gb = gb.reverse_complement(True, True, True, True, True, True, True)
        elif id_ == "pICSL50004":
            gb_archive = gb_archive.reverse_complement(True, True, True, True, True, True, True)

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

        # put the BsaI on 3..13 bases via plasmid rotation
        gb_bsa = (gb.seq + gb.seq).find(BsaI.site)
        gba_bsa = (gb_archive.seq + gb_archive.seq).find(BsaI.site)
        if gb_bsa < 0:
            raise RuntimeError("no BsaI site found in gb !")
        if gba_bsa < 0:
            raise RuntimeError("no BsaI site found in gb_archive !")
        gb <<= gb_bsa - 2
        gb_archive <<= gba_bsa - 2

        # Check sequences are the same, or skip
        if len(gb) != len(gb_archive):
            it.write(f"lengths differ for {id_}: {len(gb)} VS {len(gb_archive)}")
        elif gb.seq != gb_archive.seq:
            for i, (x, y) in enumerate(zip(gb.seq, gb_archive.seq)):
                if x != y:
                    it.write(f"sequences differ for {id_}: at position {i}, direct={x}, archive={y}")

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

        # remove duplicate Sm/Sp and SmR
        sm_sp = next(get_features("Sm/Sp no DraIII"), None) or next(get_features("SpecR"), None)
        smr = next(get_features("SmR"), None)
        if sm_sp is not None and smr is not None:
            gb_archive.features.remove(smr)
            sm_sp.qualifiers = {
                'label': ['SmR'],
                'gene': 'aadA',
                'product': 'aminoglycoside adenylyltransferase',
                'function': 'spectinomycin and streptomycin resistance',
                'note': ['color: #ffc400', 'no DraIII restriction site'],
            }

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

        # SmR recolor and annotations
        smr = next(get_features('SmR'), None) or next(get_features('SpecR'), None) or next(get_features("Sm/Sp no DraIII"), None)
        if smr is not None:
            smr.qualifiers.update({
                'label': ['SmR'],
                'gene': 'aadA',
                'product': 'aminoglycoside adenylyltransferase',
                'function': 'spectinomycin and streptomycin resistance',
                'note': ['color: #ffc400'],
            })
        smr_prom = next(get_features("SmR Promoter"), None)
        if smr_prom is None and ampr_prom is not None and smr is not None:
            smr_prom = ampr_prom
            smr_prom.qualifiers["label"] = ["SmR Promoter"]
            smr_prom.qualifiers["note"] = ["color: #ffe080"]

        # # Remove duplicate crtI
        # if len(list(get_features('crtI'))) == 2:
        #     x, y = get_features('crtI')
        #     gb_archive.features.remove(x if x.type == 'misc_feature' else y)
        #
        # # Remove duplicate crtB
        # if len(list(get_features('crtB'))) == 2:
        #     x, y = get_features('crtB')
        #     gb_archive.features.remove(x if x.type == 'misc_feature' else y)
        #
        # # Remove duplicate oriV
        # if len(list(get_features('oriV'))) == 2:
        #     gb_archive.features.remove(min(get_features('oriV'), key=lambda f: len(f.qualifiers)))

        # eGFP recolor and annotations
        egfp = next(get_features("eGFP"), None) or next(get_features("EGFP"), None)
        if egfp is not None:
            egfp.qualifiers.update(
                {
                    "label": "eGFP",
                    "note": [
                        "mammalian codon-optimized",
                        "color: #34ff03",
                    ],
                    "product": ["enhanced green fluorescent protein"],
                    "gene": ["GFP"],
                    "db_xref": [
                        "PDB:1EMA",
                        "InterPro:IPR009017",
                        "InterPro:IPR011584",
                        "InterPro:IPR000786",
                        "PFAM:PF01353",
                        "GO:0008218",
                        "GO:0006091",
                        "GO:0018298",
                        "UniProtKB/Swiss-Prot:P42212",
                    ],
                }
            )

        # eYFP
        eyfp = next(get_features("eYFP"), None)
        if eyfp is not None:
            eyfp.qualifiers.update(
                {
                    "label": "eYFP",
                    "note": [
                        "mammalian codon-optimized",
                        "color: #ccff00"
                    ],
                    "product": ["enhanced yellow fluorescent protein"],
                    "gene": ["YFP"],
                    "db_xref": [
                        "GO:0008218",
                        "GO:0006091",
                    ]
                }
            )

        # mCherry recolor and annotations
        mcherry = next(get_features("mCherry"), None)
        if mcherry is not None:
            mcherry.qualifiers.update(
                {
                    "label": "mCherry",
                    "product": "mCherry",
                    "note": [
                        "monomeric derivative of DsRed (Shaner et al., 2004)",
                        "color: #c16969",
                    ],
                    "db_xref": [
                        "UniProtKB/Swiss-Prot:A0A4D6FVK6",
                        "GO:0008218",
                        "GO:0006091",
                        "PDB:4ZIN",
                    ],
                }
            )

        # pMB1 ORI
        pmb1 = next(get_features("rep (pMB1)"), None)
        if pmb1 is not None:
            pmb1.type = "rep_origin"
            pmb1.qualifiers["label"] = ["pMB1"]
            pmb1.qualifiers["note"] = ["color: #7f7f7f"]

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

        # Remove all "/vntifkey" feature qualifier
        for feature in gb_archive.features:
            feature.qualifiers.pop("vntifkey", None)

        # recolor promoter features with #00a1ee
        end_site = BSAI_R.search(gb_archive.seq).span(0)[0]
        promoter_features = [
            f
            for f in gb_archive.features
            if any ("promoter" in label.lower() for label in f.qualifiers.get("label", []))
            and f.location.start < end_site
        ]
        if promoter_features:
            for f in promoter_features:
                f.qualifiers.setdefault("label", [])
                color = next((label for label in f.qualifiers["label"] if label.startswith("color:")), None)
                if color is not None:
                    f.qualifiers["label"].remove(color)
                f.qualifiers["label"].append("color: #00a1ee")

        # recolor terminator features with #ff8eff
        end_site = BSAI_R.search(gb_archive.seq).span(0)[0]
        terminator_features = [
            f
            for f in gb_archive.features
            if any ("terminator" in label.lower() for label in f.qualifiers.get("label", []))
            and f.location.start < end_site
        ]
        if terminator_features:
            for f in terminator_features:
                f.qualifiers.setdefault("label", [])
                color = next((label for label in f.qualifiers["label"] if label.startswith("color:")), None)
                if color is not None:
                    f.qualifiers["label"].remove(color)
                f.qualifiers["label"].append("color: #ff8eff")

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
            os.path.join(__file__, "..", "..", "moclo-plant", "registry", "plant")
        )
        with fs.open_fs(os.path.join(__file__, "..", ".."), create=True) as dst_fs:
            dir_fs = dst_fs.makedirs(fs.path.join("moclo-plant", "registry", "plant"), recreate=True)
            with dir_fs.open("{}.gb".format(info["id"]), "w") as dst_file:
                write(gb_archive, dst_file, "gb")

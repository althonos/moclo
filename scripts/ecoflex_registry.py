#!/usr/bin/env python3
# coding: utf-8
"""Automatic CIDAR sequences annotation pipeline.

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

ZIP_URL = "https://media.addgene.org/cms/filer_public/1a/00/1a00a9f1-608f-453a-937a-7f46cf872dfc/ecoflex-kit-genbank-files.zip"
URL = "https://www.addgene.org/cloning/moclo/freemont-ecoflex/"
UA = "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"

# Part sequence for automatic annotation / annotation relocation

AMPR_TERM = DNARegex("gattatcaaaaaggatctt")  # Reverse 3' of AmpR terminator
BB_PREFIX = DNARegex("gaattcgcggccgcttctag")
CMR_PROMOTER = DNARegex('tttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaa')
CMR_TERMINATOR = DNARegex('accaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaa')
AMPR_PROMOTER = DNARegex('actcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccacctg')
AMPR_TERMINATOR = DNARegex('gattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacag')


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
            "/EcoFlex - GenBank/"
        )

    # load inventory
    inventory = soup.find("table", class_="kit-inventory-table")
    it = tqdm.tqdm(inventory.find_all("tr")[1:])
    for row in it:

        # extract each row
        row_text = row.find("a").text

        # get antibiotics resistances
        resistance = row.find("span", class_="resistance-spacing").text.strip()
        name = id_ = row_text.strip()

        # Update the progress bar
        it.set_description(id_)

        # TODO: entry vector not supported
        if id_ in ('pBP', 'pBP-ORF', 'pBP-lacZ'):
            continue

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
        path = next(archive.walk.files('/', filter=['{}.gb'.format(id_)]), None)
        if path is None:
            print("COULD NOT FIND", id_)
            continue
        with archive.open(path) as f:
            gb = CircularRecord(read(f, "gb"))



        # Copy well documented information from one record to the other
        gb.seq = gb.seq.upper()
        gb.seq.alphabet = IUPAC.unambiguous_dna
        gb.id = id_
        gb.name = name

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
        pref = next(get_features_from_note("BioBrick prefix"))
        if pref.location is None:
            match = BB_PREFIX.search(gb)
            pref.location = FeatureLocation(
                start=match.start(),
                end=match.end(),
                strand=1,
            )
        gb <<= pref.location.start - 1

        # AmpR recolor and annotations
        ampr = next(get_features_from_note("AmpR"), None)
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

            old_prom = next(get_features_from_note('AmpR promoter'), None)
            if old_prom is not None:
                gb.features.remove(old_prom)

            ampr_prom = next(get_features("AmpR promoter"), None)
            if ampr_prom is None:
                start, end = AMPR_PROMOTER.search(gb.seq).span()
                ampr_prom = SeqFeature(FeatureLocation(start, end, -1))
                gb.features.append(ampr_prom)
            ampr_prom.type = "promoter"
            ampr_prom.qualifiers["label"] = ["AmpR Promoter"]
            ampr_prom.qualifiers["note"] = ["color: #ff6666"]

            ampr_term = next(get_features("AmpR terminator"), None)
            if ampr_term is None:
                start, end = AMPR_TERMINATOR.search(gb.seq).span()
                ampr_term = SeqFeature(FeatureLocation(start, end, -1))
                gb.features.append(ampr_term)
            ampr_term.type = 'terminator'
            ampr_term.qualifiers['label'] = 'AmpR Terminator'
            ampr_term.qualifiers['note'] = ['color: #ff6666']

        # CmR recolor and annotations
        cmr = next(get_features_from_note('CmR'), None)
        if cmr is not None:
            cmr.qualifiers.update(
                {
                    "codon_start": [1],
                    "gene": ["cat"],
                    "product": ["chloramphenicol acetyltransferase"],
                    "label": ["CmR"],
                    "function": ["chloramphenicol resistance"],
                    "note": ["color: #0000ff; direction: LEFT"],
                    "EC_number": ["2.3.1.28"],
                    "db_xref": [
                        "UniProtKB/Swiss-Prot:P62577",
                        "GO:0008811",
                        "GO:0016740",
                        "GO:0016746",
                        "GO:0046677",
                        "PFAM:PF00302",
                    ],
                }
            )

            cmr_prom = next(get_features("CamR Promoter"), None)
            if cmr_prom is None:
                start, end = CMR_PROMOTER.search(gb.seq).span()
                cmr_prom = SeqFeature(location=FeatureLocation(start, end, -1))
                gb.features.append(cmr_prom)
            cmr_prom.type = "promoter"
            cmr_prom.qualifiers.update(
                {
                    "label": ["CmR Promoter"],
                    "note": ["color: #66ccff; direction: LEFT"],
                }
            )

            cmr_term = next(get_features("CamR Terminator"), None)
            if cmr_term is None:
                start, end = CMR_TERMINATOR.search(gb.seq).span()
                cmr_term = SeqFeature(location=FeatureLocation(start, end, -1))
                gb.features.append(cmr_term)
            cmr_term.type = "terminator"
            cmr_term.qualifiers.update(
                {
                    "label": ["CmR Terminator"],
                    "note": ["color: #66ccff; direction: LEFT"],
                }
            )

            old_term = next(get_features_from_note('lambda t0 terminator'), None)
            if old_term is not None:
                gb.features.remove(old_term)


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


        # if any(f.location is None for f in gb.features):
        #     continue
        for f in gb.features:
            if f.location is None:
                print(gb, f)

        # sort features by start location, source always first
        gb.features.sort(
            key=lambda f: (-len(gb.seq)) * (f.type == "source")
            + f.location.start
        )

        # translate color from notes to ApEinfo
        for feature in gb.features:
            translate_color(feature)

        # Add an EcoFlex article reference
        ref = Reference()
        ref.authors = 'Moore SJ, Lai HE, Kelwick RJ, Chee SM, Bell DJ, Polizzi KM, Freemont PS.'
        ref.title = 'EcoFlex: A Multifunctional MoClo Kit for E. coli Synthetic Biology.'
        ref.journal = 'ACS Synth Biol 2016;5:1059-1069.'
        ref.pubmed_id = '27096716'
        gb.annotations['references'].insert(0, ref)

        # Fix the direct submission reference
        ref = gb.annotations["references"][-1]
        ref.authors = "Larralde M"
        ref.journal = "Distributed with the MoClo Python library\nhttps://github.com/althonos/moclo"

        # write the final record
        dst_dir = os.path.abspath(
            os.path.join(
                __file__, "..", "..", "moclo-ecoflex", "registry", "ecoflex"
            )
        )
        dst_file = os.path.join(dst_dir, "{}.gb").format(info["id"])
        write(gb, dst_file, "gb")

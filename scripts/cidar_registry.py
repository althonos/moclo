#!/usr/bin/env python3
# coding: utf-8
"""Automatic PTK sequences annotation pipeline.

Downloads depositor and AddGene full sequences from each PTK plasmid page,
rotate the plasmids so that they have the same orientation as YTK official
sequences, merge the two records, remove duplicated annotations, add BsaI
sites as annotations, and then:

Edits:
- rename all chloramphenicol resistance associated sequences (promoter, CDS,
  terminator) to use the same nomenclature and color scheme as the YTK records
- improve annotations of pPTK001, pPKT002, pPTK003, pPTK004 (promoters) using
  the same format, using '/label', '/function' and '/gene' and '/note' in a
  standardized way.
- change type of 'pENO1' to 'promoter' and use YTK color
- set type of 'RFP' to CDS, remove duplicate annotaiton, add '/db_xref'
  qualifiers with UniProt/GO/PFAM/PDB/InterPro, set '/gene' to 'mCherry'
  (accepted name), use YTK color scheme
-

"""

import copy
import io
import json
import re
import os
import warnings
import sys

import tqdm
import requests
import bs4 as bs
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

from moclo.record import CircularRecord
from moclo.regex import DNARegex


URL = "https://www.addgene.org/cloning/moclo/densmore/"
UA = "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"

DVK_AF = None
AMPR_TERM = "gattatcaaaaaggatctt".upper()  # Reverse 3' of AmpR terminator

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

    # load inventory
    inventory = soup.find("table", class_="kit-inventory-table")
    for row in tqdm.tqdm(inventory.find_all("tr")[1:]):

        # extract each row
        row_text = row.find("a").text

        match = NAME_REGEX.match(row_text)
        if match is not None:
            id_ = match.group(1) + match.group(3)
            name = match.group(2)
        else:
            id_ = name = row_text

        # id_, type_, name = map(str.strip, row_text.split("-", 2))
        info = {
            "resistance": row.find("span", class_="resistance-spacing").text,
            # "name": id_,
            "id": id_,
            # "type": type_,
            "location": row.find("b").text.strip().replace(" / ", ""),
            "addgene_id": row.find("a").get("href").strip("/"),
        }

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

        # Make sure the name and id are not empty
        gb.seq.alphabet = IUPAC.unambiguous_dna
        gb.seq = gb.seq.upper()
        gb.id = id_
        gb.name = name

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
        pref = next(
            f
            for f in gb.features
            if "BioBrick prefix" in f.qualifiers.get("label", [])
        )
        gb <<= pref.location.start - 1

        # Fix bad position of lacz-alpha
        if id_.startswith(("DVA_", "DVK_")):
            lacz = next(
                f
                for f in gb.features
                if "lacZ-alpha" in f.qualifiers.get("label", [])
            )
            start = gb.seq.lower().find("ctagagactagtgggtctca")
            end = start + 324
            lacz.location = FeatureLocation(start, end, -1)

        # Fix bad position of AmpR
        elif id_.endswith("m_AF"):
            ampr = next(
                f
                for f in gb.features
                if "AmpR" in f.qualifiers.get("label", [])
            )
            start = gb.seq.lower().find("taccaatgcttaatcagtg")
            end = start + 861
            ampr.location = FeatureLocation(start, end, -1)

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

        ampr_prom = next(get_features("AmpR Promoter"), None)
        if ampr_prom is not None:
            ampr_prom.qualifiers["note"] = ["color: #ff6666"]

        ampr_term_start = gb.seq.find(AMPR_TERM)
        if ampr_term_start >= 0:
            ampr_term = SeqFeature(
                location=FeatureLocation(
                    ampr_term_start, ampr_term_start + 94, -1
                ),
                type="terminator",
                qualifiers={
                    "label": "AmpR Terminator",
                    "note": "color: #ff6666",
                },
            )
            gb.features.append(ampr_term)

        # sort features by start location, source first
        gb.features.sort(
            key=lambda f: (-len(gb.seq)) * (f.type == "source")
            + f.location.start
        )

        # translate color from notes to ApEinfo
        for feature in gb.features:
            translate_color(feature)

        # Fix the direct submission annotation
        if gb.annotations["references"][-1].title == "Direct Submission":
            ref = gb.annotations["references"][-1]
        else:
            ref = Reference()
            ref.title = "Direct Submission"
            gb.annotations.append(ref)
        ref.authors = "Larralde M"
        ref.journal = "Distributed with the MoClo Python library\nhttps://github.com/althonos/moclo"

        # write the final record
        dst_dir = os.path.abspath(
            os.path.join(
                __file__, "..", "..", "moclo-cidar", "registry", "cidar"
            )
        )
        dst_file = os.path.join(dst_dir, "{}.gb").format(info["id"])
        write(gb, dst_file, "gb")

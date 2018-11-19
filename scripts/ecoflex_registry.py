#!/usr/bin/env python3
# coding: utf-8
"""Automatic EcoFlex sequences annotation pipeline.

Edits:
- Recolor all AmpR with the same color as YTK parts
- Add AmpR terminator feature with standard color

"""

import io
import itertools
import re
import os
import warnings

import bs4 as bs
import fs.path
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

from features import annotate, translate_color

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


FULL_SEQUENCES = {
    "pBP-BBa_B0034": "https://www.addgene.org/72980/sequences/",
    "pBP-SJM901": "https://www.addgene.org/72966/sequences/",
}

# Partial sequences from the reference EcoFlex paper
PROMOTERS = {
    "pBP-SJM901": "CTATTTTACAGCTAGCTCAGTCCTAGGTATAATGCTAGCGTAC",
    "pBP-SJM902": "CTATTTTACAGCTAGCTCAGTCCTAGGGATTATGCTAGCGTAC",
    "pBP-SJM903": "CTATCTTATAGCTAGCTCAGTCCTTGGGATTATGCTAGCGTAC",
    "pBP-SJM905": "CTATTTTATAGCTAGCTCAGTCCTTGGGATTATGCTAGCGTAC",
    "pBP-SJM906": "CTATTTGATGGCTAGCTCAGTCCTAGGGATTGTGCTAGCGTAC",
    "pBP-SJM908": "CTATTTTATAGCTAGCTCAGCCCTTGGTATTATGCTAGCGTAC",
    "pBP-SJM910": "CTATTTGATGGCTAGCTCAGTCCTTGGTATTATGCTAGCGTAC",
    "pBP-SJM911": "CTATTTGACAGCTAGCTCAGTCCTTGGTACTGTGCTAGCGTAC",
    "pBP-SJM912": "CTATTTGATAGCTAGCTCAGTCCTAGGTACTATGCTAGCGTAC",
    "pBP-SJM914": "CTATTTGATGGCTAGCTCAGTCCTAGGGATTGTGCTAGCGTAC",
    "pBP-SJM915": "CTATTTTATGGCTAGCTCAGTCCTTGGTATTATGCTAGCGTAC",
}


if __name__ == "__main__":

    warnings.simplefilter("ignore")
    session = requests.Session()

    # load the kit inventory page
    with session.get(URL) as res:
        soup = bs.BeautifulSoup(res.text, "html.parser")

    # load the zip archive
    with session.get(ZIP_URL) as res:
        archive = ReadZipFS(six.BytesIO(res.content)).opendir("/EcoFlex - GenBank/")

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
        elif id_ == "pBP-T7_RBS-His6-Thrombin":
            name = id_ = "pBP-T7-RBS-His6"
        elif id_.startswith("pBP-T7_"):
            name = id_ = id_.replace("_", "-")
        elif id_.startswith("pBP-ORF-"):
            name = id_ = id_.replace("pBP-ORF-", "pBP-")
        elif id_ == "pBP-HexHis":
            name = id_ = "pBP-His6"
        elif id_.startswith("pBP_BBa"):
            name = id_ = id_.replace("pBP_BBa", "pBP-BBa")

        # extract info
        info = {
            "resistance": resistance,
            # "name": id_,
            "id": id_,
            # "type": type_,
            "location": row.find("b").text.strip().replace(" / ", ""),
            "addgene_id": row.find("a").get("href").strip("/"),
        }

        # get the online full sequence
        if id_ in FULL_SEQUENCES:
            # Load the AddGene sequences page and get the full sequence
            with requests.get(FULL_SEQUENCES[id_]) as res:
                soup = bs.BeautifulSoup(res.text, "html.parser")
                section = soup.find("section", id="depositor-full")
                gb_url = soup.find("a", class_="genbank-file-download").get('href')
            # Get the Genbank file
            with requests.get(gb_url) as res:
                gb = CircularRecord(read(io.StringIO(res.text), "gb"))

        # get the pBP-SJM901 sequence and patch it
        elif id_.startswith("pBP-SJM"):
            # get pBP-SJM
            # Load the AddGene sequences page and get the full sequence
            with requests.get(FULL_SEQUENCES["pBP-SJM901"]) as res:
                soup = bs.BeautifulSoup(res.text, "html.parser")
                section = soup.find("section", id="depositor-full")
                gb_url = soup.find("a", class_="genbank-file-download").get('href')
            # Get the Genbank file
            with requests.get(gb_url) as res:
                gb = CircularRecord(read(io.StringIO(res.text), "gb"))
            # replace the target sequence
            gb.seq = Seq(
                str(gb.seq.upper()).replace(PROMOTERS["pBP-SJM901"], PROMOTERS[id_])
            )
            gb.description = gb.description.replace("SJM901", id_[4:])
            gb.keywords = [id_[4:]]

        # get the ZIP sequence
        else:
            path = next(
                (
                    f
                    for f in archive.walk.files('/')
                    if fs.path.basename(f).lower() == '{}.gb'.format(id_).lower()
                ),
                None,
            )
            if id_ == "pBP-His6":
                path = "/Level 0/Tags/pBP-His6_tag.gb"
            elif id_ == "pBP-T7-RBS-His6":
                path = "/Level 0/T7 parts/pBP-T7_RBS_His6.gb"
            elif id_ == "pBP-T7-RBS":
                path = "/Level 0/T7 parts/pBP-T7_RBS.gb"
            elif id_ == "pBP-Strep(II)":
                path = "/Level 0/Tags/pBP-StrepII_tag.gb"
            elif id_ == "pBP-pET-RBS":
                path = "/Level 0/RBS/pBP-PET_RBS.gb"
            elif id_ == "pBP-BBa_B0034":
                path = "/Level 0/Promoters/pBP_BBa_B0034.gb"
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
        gb.annotations['references'].clear()  # FIXME ?

        # quick feature accessor
        def get_features_from_label(label):
            return (
                f for f in gb.features if label in f.qualifiers.get("label", [])
            )

        def get_features_from_note(note):
            return (
                f for f in gb.features if note in f.qualifiers.get("note", [])
            )

        def get_features(name):
            return itertools.chain(
                get_features_from_label(name),
                get_features_from_note(name),
            )

        # Correct overlapping features by setting the origin just before the
        # biobrick prefix
        pref = next(itertools.chain(
            get_features("BioBrick prefix"),
            get_features_from_note("BioBrick prefix")
        ))
        if pref.location is None:
            match = BB_PREFIX.search(gb)
            pref.location = FeatureLocation(
                start=match.start(),
                end=match.end(),
                strand=1,
            )
        gb <<= pref.location.start - 1

        # AmpR recolor and annotations
        ampr = next(get_features("AmpR"), None)
        if ampr is not None:
            annotate("ampr", ampr, gb.seq)

            old_prom = next(get_features_from_note('AmpR promoter'), None)
            if old_prom is not None:
                gb.features.remove(old_prom)

            ampr_prom = next(get_features_from_label("AmpR promoter"), None)
            if ampr_prom is None:
                start, end = AMPR_PROMOTER.search(gb.seq).span()
                ampr_prom = SeqFeature(FeatureLocation(start, end, -1))
                gb.features.append(ampr_prom)
            ampr_prom.type = "promoter"
            ampr_prom.qualifiers["label"] = ["AmpR Promoter"]
            ampr_prom.qualifiers["note"] = ["color: #ff6666"]

            ampr_term = next(get_features_from_label("AmpR terminator"), None)
            if ampr_term is None:
                start, end = AMPR_TERMINATOR.search(gb.seq).span()
                ampr_term = SeqFeature(FeatureLocation(start, end, -1))
                gb.features.append(ampr_term)
            ampr_term.type = 'terminator'
            ampr_term.qualifiers['label'] = 'AmpR Terminator'
            ampr_term.qualifiers['note'] = ['color: #ff6666']

        # CmR recolor and annotations
        cmr = next(get_features('CmR'), None)
        if cmr is not None:
            annotate("cmr", cmr, gb.seq)

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

            cmr_term = next(get_features_from_label("CamR Terminator"), None)
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
        gfp = next(get_features_from_label("GFP"), None)
        if gfp is not None:
            annotate("gfp", gfp, gb.seq)

        # mRFP1 recolor and annotations
        rfp = next(get_features_from_label("mRFP1"), None)
        if rfp is not None:
            annotate("mrfp", rfp, gb.seq)

        # patch pBP-SJM promoters
        if id_.startswith("pBP-SJM"):
            promoter = next(get_features_from_label("J23119 promoter"))
            promoter.type = "promoter"
            promoter.qualifiers.update({
                "function": ["strong constitutive promoter"],
                "note": ["color: #00a1ee; direction: RIGHT"],
            })
            if id_ == "pBP-SJM901":
                promoter.qualifiers['label'] = "J23119 Promoter"
                promoter.qualifiers['note'].insert(0, "Anderson series consensus promoter")
            else:
                promoter.qualifiers['label'] = "{} Promoter".format(id_[4:])
                promoter.qualifiers['note'].insert(0, "derived from pBP-SJM901 (BBa_J23119)")

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
        gb.annotations['references'].append(ref)

        # Fix the direct submission reference
        ref = Reference()
        ref.authors = "Larralde M"
        ref.title = "Direct Submission"
        ref.journal = "Distributed with the MoClo Python library\nhttps://github.com/althonos/moclo"
        gb.annotations['references'].append(ref)

        # write the final record
        dst_dir = os.path.abspath(
            os.path.join(
                __file__, "..", "..", "moclo-ecoflex", "registry", "ecoflex"
            )
        )
        dst_file = os.path.join(dst_dir, "{}.gb").format(info["id"])
        write(gb, dst_file, "gb")

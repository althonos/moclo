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
from Bio.Seq import Seq, translate
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqIO import read, write
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import BsaI

from moclo.record import CircularRecord
from moclo.regex import DNARegex

from features import annotate, translate_color


URL = "https://www.addgene.org/kits/sieber-moclo-pichia-toolkit/#protocols-and-resources"
UA = "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.99 Safari/537.36"

# minimal Sequence from Bxb1 attB site from 10.1371/journal.pgen.1003490
ATTB = DNARegex("GGCC(GGCTTGTCGACGACGGCGGTCNCCGTCGTCAGGATCAT)CCGG")
ALPHA_MF_DELTA = None

COLOR_REGEX = re.compile(r"color: (#[0-9a-fA-F]{6})")


if __name__ == "__main__":

    warnings.simplefilter("ignore")
    session = requests.Session()

    # load the kit inventory page
    with session.get(URL) as res:
        soup = bs.BeautifulSoup(res.text, "html.parser")

    # load inventory
    inventory = soup.find("table", class_="kit-inventory-table")
    it = tqdm.tqdm(inventory.find_all("tr")[1:])

    for row in it:

        # extract each row
        row_text = row.find("a").text
        id_, type_, name = map(str.strip, row_text.split("-", 2))
        it.set_description(id_)
        info = {
            "resistance": row.find("span", class_="resistance-spacing").text,
            "name": name,
            "id": id_,
            "type": type_,
            "location": row.find("b").text.strip().replace(" / ", ""),
            "addgene_id": row.find("a").get("href").strip("/"),
        }

        # get the AddGene sequences page
        url = "https://www.addgene.org/{}/sequences/".format(info["addgene_id"])
        with session.get(url) as res:
            soup = bs.BeautifulSoup(res.text, "html.parser")

        # get the addgene full sequence
        section = soup.find("section", id="addgene-full")
        gb_url = section.find("a", class_="genbank-file-download").get("href")
        with session.get(gb_url) as res:
            gbd = info["gb_depositor"] = CircularRecord(read(io.StringIO(res.text), "gb"))

        # get the AddGene plasmid page
        url = "https://www.addgene.org/{}/".format(info["addgene_id"])
        with session.get(url) as res:
            soup = bs.BeautifulSoup(res.text, "html.parser")

        # get the deposited record
        section = soup.find("ul", class_="addgene-document-list")
        gb_url = section.find("a").get("href")
        with session.get(gb_url) as res:
            gba = info["gb_addgene"] = CircularRecord(read(io.StringIO(res.text), "gb"))

        # Sanity check
        if len(gba) != len(gbd):
            raise RuntimeError("sequences do not have the same length")
        if gba.seq not in (gbd.seq + gbd.seq):
            gbd = gbd.reverse_complement(id=True, name=True, annotations=True)
        if gba.seq not in (gbd.seq + gbd.seq):
            raise RuntimeError("sequences differ")

        # put the BsaI on 3..13 bases via plasmid rotation
        gba_bsa = (gba.seq + gba.seq).find(BsaI.site)
        gbd_bsa = (gbd.seq + gbd.seq).find(BsaI.site)
        if gba_bsa < 0:
            raise RuntimeError("no BsaI site found in gba !")
        if gbd_bsa < 0:
            raise RuntimeError("no BsaI site found in gbd !")
        gba <<= gba_bsa - 2
        gbd <<= gbd_bsa - 2

        # Check sequences are exactly the same so we don't fuck up features
        if not gba.seq == gbd.seq:
            raise RuntimeError("alignment failed")

        # Fix duplicated features and add additional annotations
        features = gba.features + gbd.features

        def get_features(label):
            return (f for f in features if label in f.qualifiers.get("label", []))

        def get_features_from_note(note):
            return (f for f in features if note in f.qualifiers.get("note", []))

        # Add direct BsaI site
        features.append(
            SeqFeature(
                type="protein_bind",
                qualifiers={
                    # "label": ["BsaI"],
                    "bound_moiety": ["BsaI"],
                    "note": ["color: #ff0000; direction: RIGHT"],
                    # "note": ["This forward directional feature has 2 segments:\n"
                    #          "1: 3 ..  8 / #ff0000\n"
                    #          "2: 10 .. 13 / #ff0000\n"]
                },
                location=CompoundLocation(
                    [FeatureLocation(2, 8, strand=1), FeatureLocation(9, 13, strand=1)]
                ),
            )
        )

        # Add reversed BsaI site
        BsaI_site_r = Seq(BsaI.site).reverse_complement()
        pos = (gba.seq + gba.seq).find(BsaI_site_r)
        features.append(
            SeqFeature(
                type="protein_bind",
                qualifiers={
                    # "label": ["BsaI(1)"],
                    "bound_moiety": ["BsaI"],
                    "note": ["color: #ff0000; direction: LEFT"],
                    # "note": [("""This forward directional feature has 2 segments:
                    #                  1: {} .. {} / #ff0000
                    #                  2: {} .. {} / #ff0000
                    #           """).format(pos+1, pos + 6, pos - 4, pos - 1)]
                },
                location=CompoundLocation(
                    [
                        FeatureLocation(pos, pos + 6, strand=-1),
                        FeatureLocation(pos - 5, pos - 1, strand=-1),
                    ]
                ),
            )
        )

        # Fix annotations of Chloramphenicol resistance cassette
        camr, cmr = next(get_features("CamR")), next(get_features("CmR"))
        features.remove(camr)
        cmr.qualifiers.update({k: v for k, v in camr.qualifiers.items() if k not in cmr.qualifiers})
        annotate("cmr", cmr, gba.seq)

        cmr_prom = next(get_features("CamR Promoter"))
        annotate("cmr-prom", cmr_prom, gba.seq)
        cmr_term = next(get_features("CamR Terminator"))
        annotate("cmr-term", cmr_term, gba.seq)

        # Make sure ColE1 is grayed
        cole1 = next(get_features("ColE1"))
        cole1.type = "rep_origin"
        cole1.qualifiers["note"] = ["color: #7f7f7f"]

        # plasmid-specific annotations
        if id_ == "pPTK001":
            p1, p2 = [
                f
                for f in features
                if any("AOX1" in l for l in f.qualifiers.get("gene", []))
            ]
            features.remove(p2)
            p1.type = "promoter"
            p1.qualifiers.update(
                {
                    "label": ["PpAOX1 Promoter"],
                    "gene": ["Pichia pastoris AOX1"],
                    "function": ["methanol inducible promoter"],
                    "note": [
                        "promoter for Pichia pastoris alcohol oxydase",
                        "color: #00a1ee; direction: RIGHT",
                    ],
                }
            )

        elif id_ == "pPTK002":
            p1 = next(get_features_from_note("GAP promoter"))
            p2 = next(get_features("GAP promoter"))
            features.remove(p2)
            p1.type = "promoter"
            p1.qualifiers.update(
                {
                    "label": ["PpGAP Promoter"],
                    "gene": ["Pichia pastoris GAP"],
                    "function": ["glucose inducible promoter"],
                    "note": [
                        "promoter for Pichia pastoris glyceraldehyde-3-phosphate dehydrogenase",
                        "color: #00a1ee; direction: RIGHT",
                    ],
                }
            )

            mod = next(get_features_from_note("Modification for no digestion"))
            mod.type = "misc_difference"
            mod.qualifiers = {"replace": ["c"], "note": ["remove BsmBI cutting site"]}

        elif id_ == "pPTK003":
            pENO = next(get_features_from_note("pENO"))
            pENO.type = "promoter"
            pENO.qualifiers = {
                "label": ["PpENO1 Promoter"],
                "gene": ["Pichia pastoris ENO1"],
                "function": ["hypoxia inducible promoter"],
                "note": [
                    "promoter for Pichia pastoris alpha-enolase",
                    "color: #00a1ee; direction: RIGHT",
                ],
            }

        elif id_ == "pPTK004":
            pTP1 = next(get_features_from_note("TPI1"))
            pTP1.type = "promoter"
            pTP1.qualifiers = {
                "label": ["PpTPI1 Promoter"],
                "gene": ["Pichia pastoris TPI1"],
                "function": ["strong constitutive promoter"],
                "note": [
                    "promoter for Pichia pastoris triose-phosphate isomerase",
                    "color: #00a1ee; direction: RIGHT",
                ],
            }

        elif id_ == "pPTK005":
            a1 = next(get_features_from_note("alpha-factor secretion signal"))
            a2 = next(get_features("-alpha-factor secretion signal"))
            features.remove(a1)
            start = a2.location.start
            a2.location = CompoundLocation(
                [
                    FeatureLocation(start, start + 57, 1),
                    FeatureLocation(start + 57, start + 57 + 198, 1),
                    FeatureLocation(start + 57 + 198, start + 57 + 198 + 12, 1),
                ]
            )
            # QUESTION: CDS or sig_peptide ?
            a2.type = "sig_peptide"
            a2.qualifiers.update(
                {
                    "codon_start": ["1"],
                    "label": ["MF-alpha-1 signal"],
                    "gene": ["S. cerevisiae MF(ALPHA)1"],
                    "product": ["mating factor alpha-1 secretion signal"],
                    "note": [
                        "Cleavage by the Kex2 protease occurs after the "
                        "dibasic KR sequence. The EA dipeptides are then "
                        "removed by dipeptidyl aminopeptidase A.",
                        "color: #ffcbbf; direction: RIGHT",
                    ],
                    "db_xref": [
                        "UniProtKB/Swiss-Prot:P01149",
                        "InterPro:IPR008675",
                        "PFAM:PF05436",
                        "SGD:S000006108",
                    ],
                }
            )

        elif id_ == "pPTK006":
            a1 = next(get_features_from_note("alpha-factor secretion signal"))
            start = a1.location.start
            a1.location = CompoundLocation(
                [
                    FeatureLocation(start, start + 57, 1),
                    FeatureLocation(start + 57, start + 57 + 198, 1),
                ]
            )
            # QUESTION: CDS or sig_peptide ?
            a1.type = "sig_peptide"
            a1.qualifiers.update(
                {
                    "codon_start": ["1"],
                    # 'direction': ['RIGHT'],
                    "label": ["MF-alpha-1 signal (no EAEA)"],
                    "gene": ["S. cerevisiae MF(ALPHA)1"],
                    "product": ["mating factor alpha-1 secretion signal, no EAEA"],
                    "note": [
                        "Cleavage by the Kex2 protease occurs after the "
                        "dibasic KR sequence.",
                        "color: #ffcbbf; direction: RIGHT",
                    ],
                    "db_xref": [
                        "UniProtKB/Swiss-Prot:P01149",
                        "InterPro:IPR008675",
                        "PFAM:PF05436",
                        "SGD:S000006108",
                    ],
                }
            )

        elif id_ == "pPTK007":
            ad = next(get_features_from_note("Alpha"))
            # QUESTION: CDS or sig_peptide ?
            ad.type = "sig_peptide"
            ad.qualifiers.update(
                {
                    "codon_start": ["1"],
                    # 'direction': ['RIGHT'],
                    "label": ["MF-alpha-1-delta signal"],
                    "gene": ["S. cerevisiae MF(ALPHA)1"],
                    "product": ["mating factor alpha-1 secretion signal, shortened"],
                    "note": ["color: #ffcbbf; direction: RIGHT"],
                    "db_xref": [
                        "UniProtKB/Swiss-Prot:P01149",
                        "InterPro:IPR008675",
                        "PFAM:PF05436",
                        "SGD:S000006108",
                    ],
                }
            )
            ALPHA_MF_DELTA = DNARegex(
                str(ad.extract(gba.seq)).replace("GCCGCTA", "GCC[GT]CTA")
            )

        elif id_ == "pPTK008":
            adk = next(get_features_from_note("AlphaT*"))
            adk.type = "sig_peptide"
            adk.location = FeatureLocation(
                adk.location.start, adk.location.start + 57, 1
            )
            adk.qualifiers = {
                "label": ["MF-alpha-1 signal pre-sequence"],
                "gene": ["S. cerevisiae MF(ALPHA1)"],
                "product": ["mating factor alpha-1 secretion signal pre-sequence"],
                "note": ["color: #ff6600; direction: RIGHT"],
                "db_xref": [
                    "UniProtKB/Swiss-Prot:P01149",
                    "InterPro:IPR008675",
                    "PFAM:PF05436",
                    "SGD:S000006108",
                ],
            }

            moda = next(get_features_from_note("AlphaT"))
            moda.type = "misc_difference"
            moda.qualifiers = {"replace": ["G"]}

            ank = next(get_features_from_note("AlphaT 3"))
            ank.type = "sig_peptide"
            ank.location = FeatureLocation(adk.location.start + 58, ank.location.end, 1)
            ank.qualifiers = {
                "codon_start": ["1"],
                # 'direction': ['RIGHT'],
                "label": ["alpha-MF-delta (no Kex)"],
                "gene": ["S. cerevisiae MF(ALPHA)1"],
                "product": [
                    "mating factor alpha-1 secretion signal, shortened, no Kex2p site"
                ],
                "note": ["color: #ffcbbf; direction: RIGHT"],
                "db_xref": [
                    "UniProtKB/Swiss-Prot:P01149",
                    "InterPro:IPR008675",
                    "PFAM:PF05436",
                    "SGD:S000006108",
                ],
            }

        elif id_ == "pPTK009":
            am = next(get_features_from_note("alpha-amylase"))
            am.type = "sig_peptide"
            am.location = FeatureLocation(am.location.start, am.location.start + 60, +1)
            am.qualifiers = {
                "gene": ["amy"],
                "product": ["alpha-amylase signal"],
                "note": ["color: #ff6600; direction: RIGHT"],
                "label": ["Alpha-amylase signal"],
                "db_xref": ["UniProtKB/Swiss-Prot:P30292"],
            }

        elif id_ == "pPTK010":
            gam = next(get_features_from_note("Glucoamylase"))
            gam.type = "sig_peptide"
            gam.location = FeatureLocation(
                gam.location.start, gam.location.start + 54, +1
            )
            gam.qualifiers = {
                "gene": ["gaI"],
                "product": ["glucoamylase I signal"],
                "note": ["color: #ff6600; direction: RIGHT"],
                "label": ["Glucoamylase signal"],
                "db_xref": ["UniProtKB/Swiss-Prot:P23176"],
            }

        elif id_ == "pPTK011":
            alb = next(get_features_from_note("hSA"))
            alb.type = "sig_peptide"
            alb.location = FeatureLocation(
                alb.location.start, alb.location.start + 54, +1
            )
            alb.qualifiers = {
                "gene": ["ALB"],
                "product": ["serum albumin signal"],
                "note": ["color: #ff6600; direction: RIGHT"],
                "label": ["Serum Albumin signal"],
                "db_xref": ["UniProtKB/Swiss-Prot:P02768"],
            }

        elif id_ == "pPTK012":
            inu = next(get_features_from_note("Inulinase"))
            inu.type = "sig_peptide"
            inu.location = FeatureLocation(
                inu.location.start, inu.location.start + 48, +1
            )
            inu.qualifiers = {
                "gene": ["INU1"],
                "product": ["inulinase signal peptide"],
                "note": ["color: #ff6600; direction: RIGHT"],
                "label": ["Inulinase signal"],
                "db_xref": ["UniProtKB/Swiss-Prot:P28999"],
            }

        elif id_ == "pPTK013":
            inv = next(get_features_from_note("Invertase"))
            inv.type = "sig_peptide"
            inv.location = FeatureLocation(
                inv.location.start, inv.location.start + 57, +1
            )
            inv.qualifiers = {
                "gene": ["SUC1"],
                "product": ["invertase 1 signal peptide"],
                "note": ["color: #ff6600; direction: RIGHT"],
                "label": ["Invertase signal"],
                "db_xref": ["UniProtKB/Swiss-Prot:P10594"],
            }

        elif id_ == "pPTK014":
            kil = next(get_features_from_note("Killer"))
            kil.type = "sig_peptide"
            kil.location = FeatureLocation(
                kil.location.start, kil.location.start + 78, 1
            )
            kil.qualifiers = {
                "product": ["M1-1 killer toxin signal peptide"],
                "note": ["color: #ff6600; direction: RIGHT"],
                "label": ["Killer toxin signal"],
                "db_xref": ["UniProtKB/Swiss-Prot:P01546"],
            }

        elif id_ == "pPTK020":
            attb1 = next(get_features_from_note("BxbI attB"))
            attb2 = next(get_features_from_note("attB"))
            features.remove(attb1)
            start, end = ATTB.search(gba.seq).span(1)
            attb2.type = "protein_bind"
            attb2.location = FeatureLocation(start, end, 1)
            attb2.qualifiers = {
                "label": ["attB"],
                "bound_moiety": ["Bxb1 integrase"],
                "function": ["Bxb1 integrase attB recognition site"],
                "note": ["color: #99cc00"],
            }

            uas = next(get_features_from_note("UASrpg"))
            uas.type = "misc_feature"
            uas.qualifiers = {
                "label": ["TKL2 3' Homology"],
                "function": ["upstream activation site for ribosomal protein genes"],
                "note": [
                    "color: #6685ca",
                    "TEF2 UASrpg from S. cerevisiae chromosome II",
                ],
            }

        egfp2 = next(get_features("EGFP"), None)
        if egfp2 is not None:
            egfp1 = next(get_features_from_note("EGFP"))
            features.remove(egfp2)
            annotate("egfp", egfp1, gba.seq)

        rfp = next(get_features_from_note("RFP"), None)
        if rfp is not None:
            start = rfp.location.start + rfp.extract(gba.seq).find("ATG", 1)
            rfp.location = FeatureLocation(start, start + 708, 1)
            annotate("mcherry", rfp, gba.seq)
            alert = next(get_features_from_note("Start with second codon!!! ATG deleted!"))
            features.remove(alert)

        pars = next(get_features_from_note("PARS-1"), None)
        if pars is not None:
            pars.type = "rep_origin"
            pars.qualifiers["label"] = ["PARS-1"]
            pars.qualifiers["note"] = ["color: #9A969B"]

        t1 = next(get_features_from_note("AOX1 terminator"), None)
        t2 = next(get_features("AOX1 terminator"), None)
        if t1 is not None and t2 is not None:
            features.remove(t2)
            t1.qualifiers.update(
                {
                    "note": ["transcription terminator for AOX1", "color: #ff8eff"],
                    "label": ["PpAOX1 Terminator"],
                }
            )

        match = ALPHA_MF_DELTA and ALPHA_MF_DELTA.search(gba.seq)
        if id_ != "pPTK007" and match is not None:
            start, end = match.span(0)
            features.append(
                SeqFeature(
                    location=FeatureLocation(start, end, 1),
                    type="sig_peptide",
                    id="alpha-MF-delta",
                    qualifiers={
                        "label": ["alpha-MF-delta"],
                        "gene": ["S. cerevisiae MF-alpha-1"],
                        "product": ["alpha-factor shortened secretion signal"],
                        "note": ["color: #ffcbbf; direction: RIGHT"],
                    },
                )
            )

        # sort features
        features.sort(
            key=lambda feat: (-len(gba.seq)) * (feat.type == "source")
            + feat.location.start
        )

        # translate color from notes to ApEinfo
        for feature in features:
            translate_color(feature)

        # merge annotations
        annotations = copy.deepcopy(gba.annotations)
        annotations.update(gbd.annotations)

        # Fix the direct submission annotation
        ref = annotations["references"][-1]
        ref.authors = "Larralde M"
        ref.journal = "Distributed with the MoClo Python library\nhttps://github.com/althonos/moclo"

        # Add the YTK type to the record comments
        annotations["comment"] = ["YTK:{}".format(type_)]

        # create the final record
        final = CircularRecord(
            seq=gba.seq,
            id=info["id"],
            name=info["id"],
            description=info["name"],
            dbxrefs=gba.dbxrefs + gbd.dbxrefs,
            features=features,
            annotations=annotations,
        )

        # write the final record
        dst_dir = os.path.abspath(
            os.path.join(__file__, "..", "..", "moclo-ytk", "registry", "ptk")
        )
        dst_file = os.path.join(dst_dir, "{}.gb").format(info["id"])
        write(final, dst_file, "gb")

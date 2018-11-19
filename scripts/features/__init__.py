# coding: utf-8

import re

import pkg_resources
import yaml


COLOR_REGEX = re.compile(r"color: (#[0-9a-fA-F]{6})")


def annotate(name, feature, seq):
    # open YAML file with feature metadata
    with pkg_resources.resource_stream(__name__, "{}.yml".format(name)) as f:
        data = yaml.load(f)
    # patch feature
    feature.type = type_ = data.pop("type")
    feature.qualifiers.update(data)
    # add feature translation in case of a CDS
    if type_ == "CDS":
        feature.qualifiers['translation'] = feature.extract(seq).translate()


def translate_color(feature):
    # find the note with the color information
    notes = feature.qualifiers.get("note", [])
    color_note = next((n for n in notes if n.startswith("color: #")), None)
    if color_note is not None:
        # translate color to ApEinfo tags
        hex_color = COLOR_REGEX.match(color_note).group(1).lower()
        feature.qualifiers["note"].remove(color_note)
        feature.qualifiers["ApEinfo_fwdcolor"] = [hex_color]
        feature.qualifiers["ApEinfo_revcolor"] = [hex_color]
        feature.qualifiers["ApEinfo_graphicformat"] = [
            "arrow_data {{0 1 2 0 0 -1} {} 0} width 5 offset 0"
        ]

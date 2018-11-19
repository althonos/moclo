# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Imports -----------------------------------------------------------------

import configparser
import datetime
import os
import re
import sys
import shutil
import semantic_version
import sphinx.util.logging
import sphinx_bootstrap_theme
from sphinx.util.console import bold, darkgreen, reset

# -- Globals -----------------------------------------------------------------

KITS = ["cidar", "ecoflex", "gb3", "ig", "ytk"]
LIBS = ["moclo"] + ["moclo-{}".format(kit) for kit in KITS]
LOGGER = sphinx.util.logging.getLogger(__name__)

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

docssrc_dir = os.path.abspath(os.path.join(__file__, ".."))
project_dir = os.path.dirname(os.path.dirname(docssrc_dir))

# Setup path to main moclo directory
sys.path.insert(0, os.path.join(project_dir, "moclo"))
sys.path.append(os.path.join(docssrc_dir))

# Add path to moclo kits
import moclo.kits

for kit in KITS:
    name = "moclo-{}".format(kit)
    moclo.kits.__path__.append(os.path.join(project_dir, name, "moclo", "kits"))

# -- Sphinx setup ------------------------------------------------------------

from _scripts import ytk_parts, registries


def setup(app):
    # Add custom math admonition classes
    app.add_stylesheet("bootstrap-math.css")

    # Generate SVG files from template SVG
    ytk_parts.generate_svg()

    # Force `build_ext --inplace` in moclo kits
    n = len(darkgreen(max(KITS, key=len)))
    msg = "{} [{{:4.0%}}] {{:{n}}}\r".format(bold("building registries..."), n=n)
    for index, kit in enumerate(KITS):
        percent = index / (len(KITS) - 1)
        LOGGER.info(msg.format(percent, darkgreen(kit)), nonl=index != len(KITS) - 1)
        registries.build_registries(kit)

    # Copy CHANGELOG files to the doc source directory
    n = len(darkgreen(max(KITS, key=len)))
    msg = "{} [{{:4.0%}}] {{:{n}}}\r".format(bold("collecting changelogs..."), n=n)
    for index, lib in enumerate(LIBS):
        LOGGER.info(msg.format(percent, darkgreen(kit)), nonl=index != len(LIBS) - 1)
        changelog_src = os.path.join(project_dir, lib, "CHANGELOG")
        changelog_dst = os.path.join(docssrc_dir, "changes", "{}.rst".format(lib))
        if os.path.exists(changelog_src):
            with open(changelog_src, 'rb') as src:
                src.readline()  # remove title
                src.readline()  # and underline
                with open(changelog_dst, 'wb') as dst:
                    dst.write(b":tocdepth: 2\n\n")
                    dst.write("``{}``\n".format(lib).encode('utf-8'))
                    dst.write("{}\n".format("=" * (len(lib) + 4)).encode('utf-8'))
                    shutil.copyfileobj(src, dst)

# -- Project information -----------------------------------------------------


# General information
project = moclo.__name__
author = re.match('(.*) <.*>', moclo.__author__).group(1)
year = datetime.date.today().year
copyright = "{}, {}".format(
    "2018-{}".format(year) if year > 2018 else "2018",
    author,
)

# The parsed semantic version
semver = semantic_version.Version.coerce(moclo.__version__)
# The short X.Y version
version = "{v.major}.{v.minor}.{v.patch}".format(v=semver)
# The full version, including alpha/beta/rc tags
release = str(semver)

# Project URLs
_parser = configparser.ConfigParser()
_parser.read(os.path.join(project_dir, "moclo", "setup.cfg"))
project_urls = dict(
    map(str.strip, line.split(" = ", 1))
    for line in _parser.get("metadata", "project_urls").splitlines()
    if line.strip()
)


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    # "sphinx.ext.imgmath",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx_bootstrap_theme",
    "nbsphinx",
    "IPython.sphinxext.ipython_console_highlighting",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# The name of the default role for inline references
default_role = "py:obj"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "bootstrap"

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    # Bootswatch (http://bootswatch.com/) theme.
    "bootswatch_theme": "simplex",
    # Choose Bootstrap version.
    "bootstrap_version": "3",
    # Tab name for entire site. (Default: "Site")
    "navbar_site_name": "Documentation",
    # HTML navbar class (Default: "navbar") to attach to <div> element.
    # For black navbar, do "navbar navbar-inverse"
    "navbar_class": "navbar navbar-inverse",
    # Render the next and previous page links in navbar. (Default: true)
    "navbar_sidebarrel": True,
    # Render the current pages TOC in the navbar. (Default: true)
    "navbar_pagenav": False,
    # A list of tuples containing pages or urls to link to.
    "navbar_links": [
        ("GitHub", _parser.get("metadata", "home-page").strip(), True)
    ] + [
        (k, v, True) for k, v in project_urls.items() if k != "Documentation"
    ],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
html_sidebars = {
    "*": ["localtoc.html"],
    os.path.join("api", "*"): ["localtoc.html"],
    os.path.join("concepts", "*"): ["localtoc.html"],
    os.path.join("examples", "*"): ["localtoc.html"],
    os.path.join("kits", "*", "*"): ["localtoc.html"],
    os.path.join("theory", "*"): ["localtoc.html"],
    os.path.join("changes", "*"): ["localtoc.html"],
}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "moclo"


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "moclo.tex", "moclo Documentation", "Martin Larralde", "manual")
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "moclo", "moclo Documentation", [author], 1)]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "moclo",
        "moclo Documentation",
        author,
        "moclo",
        "One line description of project.",
        "Miscellaneous",
    )
]


# -- Extension configuration -------------------------------------------------

# -- Options for imgmath extension -------------------------------------------

imgmath_image_format = "svg"

# -- Options for napoleon extension ------------------------------------------

napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = True
napoleon_include_private_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True

# -- Options for autodoc extension -------------------------------------------

autoclass_content = "class"

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "six": ("http://six.readthedocs.io/", None),
}

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

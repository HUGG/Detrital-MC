# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Detrital MC'
copyright = '2020-2026, David Whipp, University of Helsinki'
author = 'David Whipp'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_nb",
    "sphinxcontrib.bibtex",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"

# HTML theme options
html_theme_options = {
    # "external_links": [],
    "repository_url": "https://github.com/HUGG/Detrital-MC",
    "repository_branch": "main",
    "path_to_docs": "docs/source/",
    "use_edit_page_button": True,
    "use_repository_button": True,
}

# Allow myst admonition style
myst_admonition_enable = True

# Define level for myst heading implicit anchors
myst_heading_anchors = 3

# Enable math config options
myst_enable_extensions = ["dollarmath"]

# MathJax config
mathjax3_config = {
    "loader": {"load": ["[tex]/upgreek"]},
    "tex": {"packages": {"[+]": ["upgreek"]}},
}

# Use bibtex for citations
bibtex_bibfiles = ["refs.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "author_year"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

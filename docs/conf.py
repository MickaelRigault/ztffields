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

project = 'ztffields'
copyright = '2022, Mickael Rigault'
author = 'Mickael Rigault'

import os, sys
sys.path.insert(0, os.path.abspath('..'))
for x in os.walk(f'../{project}'):
  sys.path.insert(0, x[0])



# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'matplotlib.sphinxext.plot_directive',
    # extra
    "numpydoc",
    'myst_nb',
    "nbsphinx",
    'sphinx_copybutton'
    ]


intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Autodoc configuration
autodoc_default_options = {
    'members': True,            # Document all members
    'undoc-members': True,      # ... including undocumented ones
    'ignore-module-all': True,  # do not stick to __all__
}
autoclass_content = "both"              # Insert class and __init__ docstrings
autodoc_member_order = "bysource"       # Keep source order

source_suffix = ['.rst', '.ipynb', '.md']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_logo = '_static/logo_ztffields.png'

html_theme = 'sphinx_book_theme'

html_theme_options = {
    'logo_only': True,
    'show_toc_level': 2,
    'repository_url': f'https://github.com/MickaelRigault/{project}',
    'use_repository_button': True,     # add a "link to repository" button
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

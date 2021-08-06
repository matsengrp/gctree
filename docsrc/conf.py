# Configuration file for the Sphinx documentation builder.
#
# Taken from https://github.com/JamesALeedham/Sphinx-Autosummary-Recursion

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
import os
import sys

sys.path.insert(0, os.path.abspath(".."))  # Source code dir relative to this file

# -- Project information -----------------------------------------------------

project = "gctree"
author = "William DeWitt"
copyright = '2020, William DeWitt'

# No version in docs, doesn't play nice with versioneer
# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # Core Sphinx library for auto html doc generation from docstrings
    "sphinx.ext.autodoc",
    # Create neat summary tables for modules/classes/methods etc
    "sphinx.ext.autosummary",
    "sphinx.ext.githubpages",
    # Link to other project's documentation (see mapping below)
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    # Add a link to the Python source code for classes, functions etc.
    "sphinx.ext.viewcode",
    # support NumPy and Google style docstrings
    "sphinx.ext.napoleon",
    # Automatically document param types (less noise in class signature)
    "sphinx_autodoc_typehints",
    # track to do list items
    "sphinx.ext.todo",
    "sphinxarg.ext",
    # render command line output
    "sphinxcontrib.programoutput",
]

# show todos in output
todo_include_todos = True

# Mappings for sphinx.ext.intersphinx. Projects have to have Sphinx-generated doc! (.inv file)
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
    'ete3': ('http://etetoolkit.org/docs/latest/', None),
}

autosummary_generate = True  # Turn on sphinx.ext.autosummary
autoclass_content = "both"  # Add __init__ doc (ie. params) to class summaries
html_show_sourcelink = (
    False
)  # Remove 'view source code' from top of page (for html, not python)
autodoc_inherit_docstrings = True  # If no class summary, inherit base class summary

autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    'logo_only': True,
}

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = '_static/logo.png'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

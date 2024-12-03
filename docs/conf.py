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
import os
import subprocess
import toml
# import sys
# sys.path.insert(0, os.path.abspath('.'))



# -- Project information -----------------------------------------------------
# obtain project information from pyproject.toml
# Path to your pyproject.toml file
pyproject_path = os.path.join(os.path.dirname(__file__), '..', 'pyproject.toml')

# Load the pyproject.toml file
with open(pyproject_path, 'r') as f:
    pyproject_data = toml.load(f)

# Extract the relevant fields from the pyproject.toml
project = pyproject_data['project']['name']
authors = pyproject_data['project']['authors']
author = ', '.join([author['name'] for author in authors])  # Comma-separated list of author names
release = pyproject_data['project']['version']
copyright = f"2024, {author}"  # Use author and release for copyright

# -- General configuration ---------------------------------------------------

# Settings to determine if we are building on readthedocs
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
if read_the_docs_build:
    subprocess.call('pwd', shell=True)
    subprocess.call('ls', shell=True)
    subprocess.call('cmake -B ../build -S .. -DBUILD_DOCS=ON', shell=True)
    subprocess.call('cmake --build ../build  --target Doxygen', shell=True)
    xml_dir = os.path.abspath("../build/docs/xml")

    

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # For Google/NumPy-style docstrings
    'sphinx.ext.viewcode',  # Links to source code
    'breathe',              # For C++ API
]

# Breathe Configuration
breathe_default_project = "polytopewalk"
if read_the_docs_build:
    breathe_projects = {"polytopewalk": xml_dir}


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autoclass_content = 'both'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
"""
InflateORG
Inflate or shrink the membrane to resolve clash between membrane and protein.
"""

# Add imports here
from .inflateorg import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

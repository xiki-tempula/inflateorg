"""
Unit and regression test for the inflateorg package.
"""

# Import package, test suite, and other packages as needed
import inflateorg
import pytest
import sys

def test_inflateorg_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "inflateorg" in sys.modules

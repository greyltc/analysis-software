from setuptools_scm import get_version
import setuptools_scm_git_archive
from pkg_resources import get_distribution, DistributionNotFound

# check if we're running from git or a github archive
try:
    __version__ = get_version(root='..', relative_to=__file__)
except LookupError:
    __version__ = None

# check if we're running from a proper install
if __version__ == None:
    try:
        __version__ = get_distribution(__name__).version
    except DistributionNotFound:
        __version__ = '0.0.0'
        pass

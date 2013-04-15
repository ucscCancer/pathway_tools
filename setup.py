import sys
import os

if sys.version_info[:2] < (2, 6):
    print "pathway_tools requires Python 2.6 or better.  Python %d.%d detected" % \
          sys.version_info[:2]
    sys.exit(-1)


from distutils.core import setup
from distutils.core import Command
from distutils.command.install import install
from distutils.command.build_py import build_py
from distutils.command.build_ext import build_ext
from distutils.extension import Extension

class test_pathwayTools(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass


    def run(self):
        os.chdir("tests")
        sys.path.insert(0, '')
        import runTests 
        runTests.main([])


PACKAGES = [
    'pathway_tools'
]

__version__ = "0.1a"

setup(
    name='pathway_tools',
    version=__version__,
    author='Kyle Ellrott',
    author_email='',
    url='https://github.com/ucscCancer/pathway_tools',
    description='UCSC Pathway Toolkit',
    download_url='https://github.com/ucscCancer/pathway_tools',
    cmdclass = {
        "test" : test_pathwayTools,
    },
    
    packages=PACKAGES,
)

import os, sys
import setuptools

"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import pathlib
from codecs import open

here = os.getenv('SRC_DIR')

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
    
if ((sys.version_info < (3, 9)) and (sys.version_info >= (3, 10))):
    sys.exit('Python>=3.9 is required by vvv2_display.')

setuptools.setup(
    name="vvv2_display",  # Required
    version="0.1.1",  # Required
    description="Viral Variant Visualizer 2 display",  # Optional
    long_description=long_description,
    long_description_content_type="text/plain",  # Optional (see note above)
    url="https://github.com/ANSES-Ploufragan/vvv2_display",  # Optional
    author="FTouzain",  # Optional
    author_email="fabrice.touzain@anses.fr",  # Optional
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 3 - Alpha",
        # Indicate who your project is intended for
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: GPL3",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ],
    keywords="display, variant, virus, viral",  # Optional
    packages=setuptools.find_packages(),
    # packages=find_packages(where='src'),  # Required
    python_requires='>=3.9',
    include_package_data=True,
    install_requires=["r-ggplot2==3.6.6",
                      "pysam==0.19.1",
                      "numpy==1.23.1"],  # Optional
    scripts=[
        os.path.join(here, "vvv2_display.py"),
        os.path.join(here, "PYTHON_SCRIPTS/convert_tbl2json.py"),
        os.path.join(here, "PYTHON_SCRIPTS/convert_vcffile_to_readablefile2.py"),
        os.path.join(here, "PYTHON_SCRIPTS/correct_multicontig_vardict_vcf.py")
        ],
    project_urls={  # Optional
        "Bug Reports": "https://github.com/ANSES-Ploufragan/vvv2_display/issues",
        "Source": "https://github.com/ANSES-Ploufragan/vvv2_display/",
    },
)

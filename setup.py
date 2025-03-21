#!/usr/bin/env python3
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
# src_dir = os.path.join(here, 'src/')

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
    
if (sys.version_info < (3, 9)):
    sys.exit('Python>=3.9 is required by vvv2_display.')

setuptools.setup(
    name="vvv2_display",  # Required
    version="0.2.3.7",  # Required
    description="Viral Variant Visualizer 2 display",  # Optional
    long_description=long_description,
    long_description_content_type="text/markdown",  # Optional (see note above)
    url="https://github.com/ANSES-Ploufragan/vvv2_display",  # Optional
    author="Fabrice Touzain",  # Optional
    author_email="fabrice.touzain@anses.fr",  # Optional
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 5 - Production/Stable",        
#        "Development Status :: 4 - Beta",
#        "Development Status :: 3 - Alpha",        
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
    packages=setuptools.find_packages(
        where='src',
         include=["vvv2_display",
                  "convert_tbl2json.py",
                  "convert_vcffile_to_readablefile2.py",
                  "correct_multicontig_vardict_vcf.py",
                  "correct_covdepth_f.py",
                  "visualize_coverage_depth.R",
                  "visualize_snp_v4.R"
         ]

    ),
    python_requires=">=3.9",
    include_package_data=True,
    install_requires=[
        "python>=3.9",
        "r-ggplot2>=3.4.4",
        "r-gridextra>=2.3",
        "r-cowplot>=1.1.1",
        "r-stringr>=1.5.1",
        "r-jsonlite>=1.8.8",
        "pysam==0.19.1",
        "numpy>=1.23.1"],  # Optional
    scripts=[
        "src/vvv2_display.py",
        "src/convert_tbl2json.py",
        "src/convert_vcffile_to_readablefile2.py",
        "src/correct_multicontig_vardict_vcf.py",
        "src/correct_covdepth_f.py",
        "src/visualize_coverage_depth.R",
        "src/visualize_snp_v4.R"
        ],
    zip_safe=False,
    project_urls={  # Optional
        "Bug Reports": "https://github.com/ANSES-Ploufragan/vvv2_display/issues",
        "Source": "https://github.com/ANSES-Ploufragan/vvv2_display/",
    }
)

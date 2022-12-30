import pathlib
import os
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '1.2.2'
PACKAGE_NAME = 'SCellBOW'
AUTHOR = 'NAMRATA BHATTACHARYA, SAM KOSHY THOMAS, DEBARKA SENGUPTA'
AUTHOR_EMAIL = 'namrata.bhattacharya@hdr.qut.edu.au, samkoshy.thomas@student.adelaide.edu.au, debarka@iiitd.ac.in'
URL = 'https://github.com/cellsemantics/SCellBOW'

LICENSE = 'Apache License 2.0'
DESCRIPTION = 'SCellBOW is an unsupervised transfer learning algorithm for clustering scRNA-seq data and performing phenotype algebra analysis on the clusters'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy',
      'pandas',
      'gensim',
      'scanpy',
      'nltk',
      'scikit-learn',
      'scikit-survival',
      'imbalanced-learn',
      'xgbse',
      'tqdm',
      'matplotlib']

# with open('./reqtxt') as f:
#     required=f.read.splitlines()

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      packages=find_packages()
      )
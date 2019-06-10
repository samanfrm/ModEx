#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'ModEx',
    version='0.1.dev1',
    author="Saman Farahmand",
    author_email="saman.farahmand001@umb.edu",
    description=(
        'A Text Mining Approach For Exploring Transcription Factor-Target Gene Direction Extraction Through PubMed Abstracts'),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/umbibio/ChIPMiner',
    packages=setuptools.find_packages(),
    license='MIT',
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent"
    ],
    install_requires=['numpy', 'pandas', 'biopython', 'nltk','networkx','jellyfish',
                      'becas','requests','string','json','re'],
    zip_safe=False,
)

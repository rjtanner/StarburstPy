# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="StarburstPy",
    version="0.9.1",
    author="Ryan Tanner",
    author_email="ryan.tanner@nasa.gov",
    description="Python Wrapper for Starburst99",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rjtanner/StarburstPy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
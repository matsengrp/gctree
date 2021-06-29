import setuptools
import versioneer


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gctree",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="William DeWitt",
    author_email="wsdewitt@gmail.com",
    description="phylogenetic inference of genotype-collapsed trees",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/matsengrp/gctree",
    entry_points={"console_scripts": ["gctree=gctree.cli:main"]},
    # packages=setuptools.find_packages(exclude=['tests', 'docs', 'docsrc']),
    # packages=['gctree'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        "PyQt5",
        "ete3",
        "biopython",
        "matplotlib",
        "pandas",
        "scipy",
        "seaborn",
        "jellyfish",
    ],
)

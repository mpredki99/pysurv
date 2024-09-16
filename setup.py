# Coding: UTF-8

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pysurv",
    version="1.0.0",
    author="Michal Predki",
    description="Ordinary, weighted, and robust least squares adjustment of surveying control networks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mpredki99/pysurv.git",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "pandas",
        "pytest"
    ],
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU General Public License v3.0",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)

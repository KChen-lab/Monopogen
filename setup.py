import setuptools
from pathlib import Path

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scPopGene",
    version="v0.0",
    author="Jinzhuang Dou",
    author_email="",
    description="",
    long_description=long_description,
    long_description_content_type="",
    url="",
    #packages=['scPopGene'], #setuptools.find_packages(),
    install_requires=[
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)


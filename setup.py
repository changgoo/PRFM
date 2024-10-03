from setuptools import setup, find_packages

with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setup(
    name="prfm",
    version="1.0",
    description="Tools for application of the PRFM theory",
    author="Chang-Goo Kim",
    author_email="changgoo@princeton.edu",
    packages=find_packages(),
    install_requires=install_requires,
)

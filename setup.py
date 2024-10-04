from setuptools import setup

with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setup(
    name="prfm",
    version="1.0",
    description="Tools for application of the PRFM theory",
    author="Chang-Goo Kim",
    author_email="changgoo@princeton.edu",
    packages=["prfm"],
    install_requires=install_requires,
    package_data={"prfm": ["prfm_ring.csv", "tigress_ncr_K24.nc"]},
)

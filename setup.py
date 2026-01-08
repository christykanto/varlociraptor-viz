from setuptools import setup, find_packages

setup(
    name="varlociraptor-viz",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "altair>=5.0.0",
        "pandas>=2.0.0",
        "pysam>=0.22.0",
    ],
)

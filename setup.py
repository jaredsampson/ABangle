from setuptools import find_packages, setup

setup(
    name="abangle",
    version="0.1.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "ABangle=abangle.cli:main",
        ],
    },
)

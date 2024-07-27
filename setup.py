from setuptools import setup, find_packages

setup(
    name="synenrich",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "gseapy",
    ],
    entry_points={
        'console_scripts': [
            'synenrich=synenrich.run:main',
        ],
    },
)

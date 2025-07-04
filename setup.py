from setuptools import setup, find_packages

setup(
    name="ktypr",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        "pyhmmer==0.10.11",
        "pandas==2.2.3",
    ],
    entry_points={
        'console_scripts': [
            'ktypr = ktypr.ktypr_cli:main'
        ]
    },
)
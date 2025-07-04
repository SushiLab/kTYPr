from setuptools import setup, find_packages

setup(
    name="ktypr",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        "pandas", "joblib",  # add other deps here
    ],
    entry_points={
        'console_scripts': [
            'ktypr = ktypr.ktypr_cli:main'
        ]
    },
)

from setuptools import setup

# Dependencies and other variables
install_requires = [
    "pyhmmer==0.10.11",
    "pandas==2.2.3",
    "pyrodigal==3.4.1",
    "clinker"
]

long_description = read('README.md')

setup(
    name="ktypr",
    version="1.0",
    author = "Samuel Miravet-Verde",
    author_email = "smiravet@ethz.ch",
    description = ("kTYPr: predicting K-antigen classifications using HMM-based annotation profiling"),
    license = "GPLv3",
    include_package_data=True,
    packages=['kTYPr'],
    install_requires=install_requires,
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords = "bioinformatics genome serotype k ktype antigen escherichia coli",
    url = "https://github.com/SushiLab/kTYPr",
    download_url = "https://github.com/SushiLab/kTYPr/archive/refs/heads/master.zip",
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    entry_points={
        'console_scripts': ['ktypr = ktypr.ktypr_cli:main']
    }
)
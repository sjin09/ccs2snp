# -*- coding: utf-8 -*-

from setuptools import setup
setup(
        name="ccs2snp", 
        version='0.0.1',
        project_urls={
            "homepage": "https://github.com/sjin09/ccs2snp",
            "repository": "https://github.com/sjin09/ccs2snp"
        },
        author='Sangjin Lee',
        author_email='sl17@sanger.ac.uk',
        license='MIT',
        classifiers=[
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9'
        ],
        entry_points={"console_scripts": ["ccs2snp = ccs2snp.__main__:main"]},
        packages=['ccs2snp'],
        package_dir={"": "src"},
        package_data={"ccs2snp": ["*.typed"]},
        install_requires=[
            'argparse==1.*,>=1.4.0', 'biopython==1.*,>=1.78.0',
            'click==7.*,>=7.0.0', 'natsort==8.*,>=8.0.0', 'numpy==1.*,>=1.20.2',
            'psutil==5.*,>=5.8.0', 'pysam==0.*,>=0.16.0'
        ],
)

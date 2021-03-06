#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=7.0', 'numpy', 'pandas', 'mygene', 'statsmodels',
                'missingpy', 'seaborn', 'matplotlib']

setup_requirements = ['pytest-runner']

test_requirements = ['pytest>=3']

setup(
    author="Xu Ren, Pei Fen Kuan",
    author_email='xuren2120@gmail.com, peifen.kuan@stonybrook.edu',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="RNASeq Age Calculation in Python",
    entry_points={
        'console_scripts': [
            'racpy=racpy.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='racpy',
    name='racpy',
    packages=find_packages(include=['racpy', 'racpy.*']),
    package_data={"racpy": ["internal_data/*.csv",
                            "internal_data/all/*.csv",
                            "internal_data/Caucasian/*.csv"]},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/reese3928/racpy',
    version='0.1.6',
    zip_safe=False,
)

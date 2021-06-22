#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=7.0', 'numpy', 'tqdm', 'ase']

test_requirements = ['pytest>=3', ]

setup(
    author="Christopher Braxton Owens",
    author_email='cbraxtonowens@gmail.com',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="TBD",
    entry_points={
        'console_scripts': [
            'pyrelate=pyrelate.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pyrelate',
    name='pyrelate',
    packages=find_packages(include=['pyrelate', 'pyrelate.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/braxtonowens/pyrelate',
    version='0.1.0',
    zip_safe=False,
)

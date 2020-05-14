# `relate`
<span style="font-size:23px;"> __RE__presenting and __L__earning on __AT__omic __E__nvironments __</span>__

[![Build Status](https://travis-ci.org/jayspendlove/relate.svg?branch=master)](https://travis-ci.org/jayspendlove/relate)


Unlike machine learning in some applications that are in all senses of the phrase 'big data', with millions if not billions of pieces of training data, in materials applications we are often limited to much smaller datasets to train a model (thousands, or even perhaps a few hundred). In addition, the accuracy requirements for materials applications are extremely high. To work with these data limitations there need to be expressive ways to _represent_ this data to allow effective machine learning. This package aims to implement some of the 'top' descriptors, and allow implementation of custom descriptors in such a way to be a user-friendly tool in materials research to facilitate effective Machine Learning.

`relate` is implemented with the use of ASE Atoms objects, with several built in descriptors, as well as the added capability of allowing the user to utilize their own descriptor function.

To see the documentation:
https://msg-byu.github.io/relate/

## Installation
To pip install from github run the following command:

`pip install -e git+git://github.com/msg-byu/relate.git#egg=relate`

#### Common installation errors
You may get an error 'No pandoc found' when installing and trying to read the long description on a setup.py. This error happens because pandoc is not installed on your system, but pypandoc (which was pip installed with the package) does not check to make sure it exists before it tries to use it in your setup.py. There are 2 solutions:

###### Solution 1
`pip uninstall pypandoc` then pip install the package again. This will activate the except statement in your setup.py that will not use the convert() pandoc function (which is not on your system), and everything will install fine.

###### Solution 2
This solution involves getting pandoc on your system.  Pypandoc will automatically install pandoc on your system for you _if you install [from the github page](https://github.com/bebraw/pypandoc/blob/master/setup.py)_. You can then proceed with pip installing `relate` on your system.

#### Optional dependencies
Depending what descriptors you plan to use, you will not need to install all packages that `relate` uses. Here are some of the dependencies for particular built in descriptors.

###### SOAP:
`pip install pycsoap`

###### ASR
SOAP

###### LER:
SOAP

`pip install annoy`

## How to use `relate`
See the 'tutorials' folder for tutorials to get started with using `relate`.

For some data to get started with, [here](https://data.mendeley.com/datasets/bz46sgxm8g/1) is a publicly available dataset which includes a variety of 100 disorientation axis grain boundaries created as part of two different publications.

## Citing

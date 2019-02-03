#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
	name='batch-iv-analysis',
	description='GUI to analyze solar cell data',
	author='Grey Christoforo',
	author_email='grey@mutovis.com',
	url='https://github.com/mutovis/analysis-software',
	entry_points={'gui_scripts': ['batch-iv-analysis = batch_iv_analysis:main', ],},
	packages=find_packages(),
)

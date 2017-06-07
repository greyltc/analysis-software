#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
	name='batch-iv-analysis',
	description='GUI to analyze solar cell data',
	author='Grey Christoforo',
	author_email='grey@christoforo.net',
	url='https://github.com/greysAcademicCode/batch-iv-analysis',
	entry_points={'gui_scripts': ['batch-iv-analysis = batch_iv_analysis.__main__:main', ],},
	packages=find_packages(),
)

from setuptools import setup, find_packages

PACKAGE = 'peaksCaller'
desc_lines = open('README.rst', 'r').readlines()

setup(name=PACKAGE, #"dfp%s" % PACKAGE,
	version = "0.0.7",
	author = "Douglas Porter",
	author_email = "dfporter@wisc.edu",
	description = ("Experimental CLIP-seq peak caller program."),
	long_description = ''.join(desc_lines),
	license = "MIT",
	include_package_data = True,
	exclude_package_data = { '': ['*.pyc']},
	url = "https://github.com/dfporter/peakCaller",
	packages = ['call_peaks', 'call_peaks.lib'],
	package_data={'call_peaks.lib':['cerevisiae_genes.bed',
		'ef4.fa', 'ef4.fa.fai'],
			'call_peaks':['extractAllRegionsAboveCutoff.pl']},
	classifiers=['Development Status :: 2 - Pre-Alpha',
	'License :: OSI Approved :: MIT License',
 	'Programming Language :: Python :: 3.0',
	],
	keywords='clip_seq HITS',
)



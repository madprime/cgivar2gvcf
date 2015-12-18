from setuptools import setup
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='cgivar2gvcf',
    version='0.1dev3',
    description='Lossy conversion of Complete Genomics var file to VCF',
    url='https://github.com/madprime/cgivar2vcf',
    author_email='mpball@gmail.com',
    license='MIT',

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
    ],

    packages=['cgivar2vcf'],
    package_data={'cgivar2vcf': ['*.pyx'],
                  '': ['LICENSE.TXT', 'requirements.txt']},

    # Core dependencies should be listed here (will be installed by pip).
    install_requires=[
        'argparse>=1.2.1', 'twobitreader==3.1.0', 'wsgiref>=0.1.2'],

)

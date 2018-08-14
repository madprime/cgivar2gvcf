from setuptools import setup
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='cgivar2gvcf',
    version='0.1.8',
    description='Conversion of Complete Genomics var file to gVCF',
    url='https://github.com/madprime/cgivar2gvcf',
    author='Mad Price Ball',
    author_email='mpball@gmail.com',
    license='MIT',

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
    ],

    packages=['cgivar2gvcf'],
    package_data={'': ['LICENSE.TXT', 'requirements.txt']},

    # Core dependencies should be listed here (will be installed by pip).
    install_requires=['twobitreader==3.1.0'],

)

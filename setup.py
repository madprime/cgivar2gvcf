from setuptools import setup
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='cgivar2vcf',
    version='0.1dev1',
    description='Lossy conversion of Complete Genomics var file to VCF',
    url='https://github.com/madprime/cgivar2vcf',
    author_email='mpball@gmail.com',
    license='MIT',

    classifiers=[
        'Development Status :: 1 - Planning',
    ],

    packages=['cgivar2vcf'],
    package_data={'cgivar2vcf': ['*.pyx'],
                  '': ['LICENSE.TXT', 'requirements.txt']},

    # Core dependencies should be listed here (will be installed by pip).
    install_requires=['Cython==0.22'],

)

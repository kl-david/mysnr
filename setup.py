from setuptools import setup, find_packages

 

setup(
    name='mySNR',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[],
    author='Kilian David',
    description='A library for calculating SNR and SNDR from power spectral density',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/username/my_library',
)
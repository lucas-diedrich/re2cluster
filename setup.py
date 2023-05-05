from setuptools import setup, find_packages


with open('README.md', 'r') as file: 
    long_description = file.read()

setup(
    name='re2cluster',
    version='0.1',
    author='Lucas Diedrich',
    author_email='ldiedric@broadinstitute.org',
    description='Repetitive unsupervised clustering based on ARBOLpy package.',
    long_description=long_description,
    packages=find_packages(),
    install_requires=[
        'scanpy', 
        'pandas',
        'sklearn',
        'matplotlib'
    ]
)

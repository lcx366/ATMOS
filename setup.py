import setuptools
from setuptools import setup 

setup(
    name='pyatmos',
    version='1.0.1',
    long_description_content_type='text/markdown',
    description='A package to estimate the atmosphere parameters',
    long_description=open('README.md', 'rb').read().decode('utf-8'),
    license='MIT',
    author='Chunxiao Li',
    author_email='lcx366@126.com',
    url='https://github.com/lcx366/ATMOS',
    classifiers=[
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data = {
        'pyatmos.data': ['*.npz'],
        },
    install_requires=[
        'scipy',
        'numpy',
        'pyshtools',
        'astropy'
        ],
)

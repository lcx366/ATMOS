from setuptools import setup,find_packages 

setup(
    name = 'pyatmos',
    version = '1.2.6',
    description = 'A package to estimate the vertical structure of atmosphere with various atmospheric density models',
    author = 'Chunxiao Li',
    author_email = 'lcx366@126.com',
    url = 'https://github.com/lcx366/ATMOS',
    license = 'MIT',
    long_description_content_type = 'text/markdown',
    long_description = open('README.md', 'rb').read().decode('utf-8'),
    keywords = ['coesa76','nrlmsise00','jb2008'],
    python_requires = '>=3.10',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        ],
    packages = find_packages(),
    include_package_data = True,
    package_data = {'pyatmos.data': ['*.npz'],},
    install_requires=[
        'scipy',
        'numpy',
        'numba',
        'pandas',
        'astropy',
        'wget'
        ]
    )

from setuptools import setup,find_packages 

setup(
    name = 'pyatmos',
    version = '1.2.1',
    description = 'A package to estimate the atmosphere parameters',
    author = 'Chunxiao Li',
    author_email = 'lcx366@126.com',
    url = 'https://github.com/lcx366/ATMOS',
    license = 'MIT',
    long_description_content_type = 'text/markdown',
    long_description = open('README.md', 'rb').read().decode('utf-8'),
    keywords = ['atmosphere models','coesa76','nrlmsise00'],
    python_requires = '>=3.6',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: MIT License',
        ],
    packages = find_packages(),
    include_package_data = True,
    package_data = {'pyatmos.data': ['*.npz'],},
    install_requires=[
        'scipy',
        'numpy',
        'astropy',
        ],
    extras_require={
        'colored-progress': ['tqdm','colorama'],
        'legendre': ['pyshtools'],
        },    
    )

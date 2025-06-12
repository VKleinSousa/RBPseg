from setuptools import setup, find_packages

setup(
    name='RBPseg',
    version='1.1.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'pandas',
        'numpy<2',  # Ensure that NumPy version is below 2
        'biopython',
        'scikit-learn',
        'umap-learn',
        'hdbscan',
        'matplotlib',
        'joblib',
    ],
    extras_require={
        'pdbfixer': ['pdbfixer'],# Recommend users install via conda
        'foldseek' : ['foldseek'],  #required for classify module
        'usalign': ['usalign'], #required for merging module
    },
    entry_points={
        'console_scripts': [
            'rbpseg-sdp=rbpseg.sdp.main:main',
            'rbpseg-merge=rbpseg.merge.merge_main:main',
            'rbpseg-classify=rbpseg.classify:main',
            'rbpseg-strclust=rbpseg.strclust.strclust:main',
            'rbpseg-prepare=rbpseg.merge.prepare_files_for_merge:main',
            'rbpseg-validate=rbpseg.validate.validate_main:main'
        ],
    },
    python_requires='==3.9.21',  # Define Python version compatibility
    author='Victor Klein-Sousa',
    description='A package for tail fiber structure analysis.',
    long_description=open('README.md').read(),  # Make sure this file exists in the root directory
    long_description_content_type='text/markdown',
    url='https://github.com/VKleinSousa/RBPseg',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
)

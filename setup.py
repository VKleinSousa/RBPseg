from setuptools import setup, find_packages

setup(
    name='RBPseg',
    version='0.1.2',
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
        'pdbfixer': ['pdbfixer'],  # Recommend users install pdbfixer via conda
    },
    entry_points={
        'console_scripts': [
            'rbpseg-sdp=rbpseg.sdp.main:main',
            'rbpseg-merge=rbpseg.merge.merge_main:main',
            'rbpseg-strclust=rbpseg.strclust.strclust:main'
        ],
    },
    python_requires='>=3.8',  
    author='Victor Klein-Sousa',
    description='A package for tail fiber structure merging and SDP analysis.',
    long_description=open('README.md').read(), 
    long_description_content_type='text/markdown',
    url='https://github.com/VKleinSousa/RBPseg',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
)

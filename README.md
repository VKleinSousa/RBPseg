# RBPseg: A Tool for Tail Fiber Structure Prediction

**RBPseg** is a pipeline designed to predict and analyze phage tail fiber proteins. It has three major modules. First, it uses structural information (ESMfold/ColabFold/Alphafold) monomeric prediction to find pseudo-domains in the fiber and fractionate its sequence to '.FASTA' files (using the sDp approach) that can be further predicted using AlphaFold-multimer as trimers. The fraction modules can be merged together into a full fiber structure. RBPseg also has a built structural clustering metric (SM/pSM) that estimate the optimal number of clusters giving a TM-score matrix. 

## Installation

To get started with `RBPseg`, follow the steps below.

### Requirements

- Python 3.8+
- Conda 
- CUDA-enabled GPU (optional for faster computations)


Clone this repository to your local machine:
```bash
git clone https://github.com/VKleinSousa/rbpseg.git
cd rbpseg
```

```bash
conda create -n rbpseg_env python=3.10
conda activate rbpseg_env
```

Install RBPseg using pip:

```bash
pip install .
```

Conda install the remaining dependencies

```
conda install -c conda-forge pdbfixer openmm
```


Optionally, verify that the installation was successful by running:

```bash
rbpseg-sdp --help
```


### Usage

Once installed, you can run RBPseg from the command line to perform structural segmentation and merging.

The main entry points are rbpseg-sdp for domain segmentation and rbpseg-merge for merging chains.
Basic Command

```bash
rbpseg-sdp -p <pdb_file> -s <min_domain_size> -pd <pair_distance_constant>
```


## Examples
#### Example 1: Creating fraction fastas using an ESMfold model

To run the segmentation on a PDB file:

```bash
rbpseg-sdp -p Examples/Example1-sDp/rbp_11.pdb -c HDBSCAN -so -np 
```

outputs: FASTA files; overlaps file; sDp plot:
![sdpExample1](./Examples/Example1-sDp/rbp_11_combined_plots.png)


#### Example 2: Merging fractions

To run the merge process:

```bash
rbpseg-merge -d Examples/Example2-Merge -of Examples/Example2-Merge/overlaps.csv -n rbp_11.pdb
```
outputs: merged file (rbp_11.pdb). 

Add ```-r ``` to run relaxation. The relaxation step needs CUDA to run. You can skip relaxation by adding the flag ``` -r False ```

#### Example 3: Finding pseudo-domains

```bash
rbpseg-sdp -p Examples/Example3-PseudoDomain/rbp_11.pdb -k 20 -sv
```

![sdpspectral](./Examples/Example3-PseudoDomain/rbp_11_combined_plots_spectral.png)


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
git clone https://github.com/VKleinSousa/RBPseg.git
cd RBPseg

```

```bash
conda create -n rbpseg_env python=3.9
conda activate rbpseg_env
```

Install RBPseg using pip:

```bash
pip install .
```

Conda install the remaining dependencies

```
conda install -c conda-forge -c bioconda foldseek pdbfixer openmm usalign
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


#### Example 3: Classifying your own RBP into TC or D-Classes

```bash
rbpseg-classify -p input_protein.pdb -o output_dir -db 0
````

```bash
  -p PDB, --pdb PDB     Path to the input PDB file.
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Directory to store output results.
  -db {0,1}, --target_db {0,1}
                        Classification target: 0 for TC classes, 1 for domain classification.

```

## Reference ##

If you applied any of these codes in your work, please consider citing:

***Towards a complete phage tail fiber structure atlas.
Victor Klein-Sousa, Aritz Roa-Eguiara, Claudia Sybille Kielkopf, Nicholas Sofos, Nicholas M. I. Taylor
bioRxiv 2024.10.28.620165; doi: https://doi.org/10.1101/2024.10.28.620165***

```
@article {Klein-Sousa2024.10.28.620165,
	author = {Klein-Sousa, Victor and Roa-Eguiara, Aritz and Kielkopf, Claudia Sybille and Sofos, Nicholas and Taylor, Nicholas M. I.},
	title = {Towards a complete phage tail fiber structure atlas.},
	elocation-id = {2024.10.28.620165},
	year = {2024},
	doi = {10.1101/2024.10.28.620165},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Bacteriophages use receptor-binding proteins (RBPs) to adhere to bacterial hosts. Understanding the structure of these RBPs can provide insights into their target interactions. Tail fibers, a prominent type of RBP, are typically elongated, flexible, and trimeric proteins, making it challenging to obtain high-resolution experimental data of their full-length structures. Recent advancements in deep learning-based protein structure prediction, such as AlphaFold2-multimer (AF2M) and ESMfold, allow for the generation of high-confidence predicted models of complete tail fibers. In this paper, we introduce RBPseg, a method that combines monomeric ESMfold predictions with a novel sigmoid distance pair (sDp) protein segmentation technique. This method segments the tail fiber sequences into smaller fractions, preserving domain boundaries. These segments are then predicted in parallel using AF2M and assembled into a full fiber model. We demonstrate that RBPseg significantly improves AF2M v2.3.2 in terms of model confidence, running time, and memory usage. To validate our approach, we used single-particle cryo-electron microscopy to analyze five tail fibers from three phages of the BASEL collection. Additionally, we conducted a structural classification of 67 fibers and their domains, which identified 16 well-defined tail fiber classes and 89 domains. Our findings suggest the existence of modular fibers as well as fibers with different sequences and shared structure, indicating possible sequence convergence, divergence, and domain swapping. We further demonstrate that these structural classes account for at least 24\% of the known tail fiber universe.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2024/10/28/2024.10.28.620165},
	eprint = {https://www.biorxiv.org/content/early/2024/10/28/2024.10.28.620165.full.pdf},
	journal = {bioRxiv}
}

```

---

## FAQ

### **Which clustering mode should I use to find my fractions?**
The choice of clustering mode depends on the type of tail fiber you are analyzing and the available computing resources. Here’s a breakdown:
- **HDBSCAN:**  
   - Used in the original study, it is reliable for capturing elongated domains.  
   - Tends to generate larger FASTA files, often grouping consecutive globular domains into a single domain.  
   - **Best for:** Overall assembly of the fiber.  
- **Spectral/Kmeans Clustering:**  
   - Produces smaller, more detailed domain fragments.  
   - Suitable for studies focused on domain organization and finer structural details.  
   - However, these smaller fragments may lead to challenges in assembling the final fiber correctly, with risks of missed assignments or poor superposition.  
   - **Best for:** Analyzing detailed domain organization.

**Recommendation:**  
If your goal is to assemble the fiber structure, HDBSCAN is the preferred option. For detailed domain studies, Spectral/Kmeans may be better suited.

---

### **I have my AlphaFold predictions of my fractions. How do I prepare my files to run the merging module?**
RBPseg-merge works best when your AlphaFold models are organized into a specific directory structure. The directory should look like this:

- * I have my AF predictions of my fractions, how do I prepare my files to run the merging module?*
    RBPseg-merge often works better if you have an specific directiory with all your AF models in this manner:
    ```
	AFmodels/
	    {fraction_number}_fibername_{model_number}
	Example:
	    0_coolfiber_0.pdb
	    0_coolfiber_1.pdb
	    ...
	    1_coolfiber_0.pdb
	    ...
	    5_coolfiber_5.pdb

    ```
If your AlphaFold results are in a single folder, you can use the following script to automatically prepare your files:

``` rbpseg/merge/prepare_files_for_merge.py``` 

---

### **How many models per fraction should I have?**
Currently, **RBPseg works best with 5 models per fraction.** We are working to generalize this in future updates.

---

### **Have more questions?**
For additional inquiries, please don't hesitate to contact us. We’re happy to help!

---

# Frequent Problems

## Superposition and Chain Pairing

### **Superposition**
In some cases—most often with beta-sandwich domains or helix-coil-helix domains—the consecutive domain may be superimposed with the wrong directionality, leading to a misaligned and incorrect assembly.

#### **Recommendations to resolve this issue:**
- **Adjust the segmentation parameters:**  
   Try running `rbpseg-sdp` with modified parameters. Increasing the values for `--min_domain_size` and `--min_ov_size` often helps resolve this issue by producing fewer but larger FASTA files. Larger models are typically easier to superimpose correctly.
- **Manual assembly:**  
   If automatic assembly fails, consider manually assembling your fractions using software like ChimeraX or PyMOL.

---

### **Chain Pairing**
Determining the correct chain pairing between fractions can be a complex challenge. For a random assignment of two fractions, there is a 1/6 chance of finding the correct pair. This probability further decreases exponentially with the number of fractions (n), following the formula \( (1/6)^{(n-1)} \).

RBPseg provides two methods to find the optimal pairing, both of which aim to minimize the distance between chains:
1. **Spherical Constraint Method (Default):**  
   This method applies a spherical constraint that prevents chains from crossing the inner chain space. However, if the algorithm cannot converge, it may randomly assign a chain.
2. **Hierarchical Clustering Method:**  
   This method uses `scipy.cluster.hierarchy.linkage` to find the chain pair with the minimal distance. Unlike the default method, it does not apply any geometric constraints.

#### **Common Problems and Causes:**
- Issues often arise in regions with loops or disordered areas where symmetry breaks subtly. In such cases, the correct chain may not always be the closest neighbor, leading to incorrect assignments.

#### **How to avoid this issue:**
- **Reduce the number of fractions:**  
   Creating fewer fractions can simplify the pairing process and reduce errors.

Version 1.1

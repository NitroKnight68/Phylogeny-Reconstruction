# Phylogeny - Reconstruction
A repository in C++/Python to reconstruct Phylogenic Trees using K-Mismatch Average Common Substrings to generate a distance matrix and UPGMA (Unweighted Pair Group Method with Arithmetic Mean) to construct a tree from the matrix.
## Installation
Clone the repository from Github.
```bash
git clone https://github.com/NitroKnight68/Phylogeny-Reconstruction.git
```
## KMACS - Distance Matrix Algorithm (C++)
Uses Enhanced Suffix Arrays and ACS with K-Mismatches to generate a Distance Matrix.
```bash
make #Only works on Linux machines
```
```bash
./kmacs example2.fa 1 #<FASTA_file> <mismatch_count>
./kmacs example3.fa 1
```
## UPGMA - Reconstruction Algorithm (Python)
Uses a library called Biopython to construct and visualise the Phylogenic Tree formed from the KMACS distance matrix.
```bash
pip install biopython
python3 UPGMA-Reconstruction.py example2.fa
python3 UPGMA-Reconstruction.py example3.fa
```
## Resources
- [Project Brief - Powerpoint Presentation](CSPE43-ADSA.pptx)
- [Research Paper - ACS Approach to Phylogenetic Reconstruction](https://sci-hub.se/http://doi.org/10.1089/cmb.2006.13.336)
- [Research Paper - Induced Sorting of Suffix & LPS Array](https://arxiv.org/pdf/1101.3448.pdf)

# Phylogeny - Reconstruction
A library in C++/Python to reconstruct Phylogenic Trees using K-Mismatch Average Common Substrings to generate a distance matrix and Neighbour-Joining Clustering Method to construct a tree from the matrix.
## KMACS - Distance Matrix Algorithm (C++)
Uses Enhanced Suffix Arrays and ACS with K-Mismatches to generate a Distance Matrix.
```bash
make
```
```bash
./kmacs randomInput.fa 1
```
## NJ - Reconstruction Algorithm (Python)
Uses a library called Biopython to construct and visualise the Phylogenic Tree formed from the KMACS distance matrix.
```bash
python NJ-Reconstruction.py randomInput.fa
```

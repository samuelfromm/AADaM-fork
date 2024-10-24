# AADaM - Antibody-Antigen Dataset Maker

This is a *slightly modified* version of the Antibody-Antigen Dataset Maker (AADaM) [1] (see also https://github.com/kaymccoy/AADaM). It is a **Python** script that takes a larger dataset downloaded from SAbDab [2] and uses it to create benchmark/testing datasets intended for **ML methods**. It may also create complementary training datasets for ML methods that use antibody-antigen structures for training.

Compared to the original Antibody-Antigen Dataset Maker, **this version includes some minor bug fixes and mild refactoring, as well as a change** in how the "number of missing residues" is calculated. For details, see the description below.

### Dataset Creation

The dataset is created from a SAbDab dataset in two steps:

1. The SAbDab dataset is split into two datasets according to the split date (-d flag). From the "before" split date dataset, only structures for which (one of) the antigen chain(s) is a peptide or protein are kept. From the "after" split date dataset, only structures for which (one of) the antigen chain(s) is a protein are kept. The "after" split date dataset is further filtered by method (-m flag), resolution (-r flag), and non-natural residues (-nx flag).

   If the -cd (--cutoffDate) flag is provided, only structures after the cutoff date are considered (both for the "before" and "after" split date datasets). If the flag --minAtomSeqresFraction is provided, structures with "too many missing residues"<sup>+</sup> **are discarded** from the "after" split date dataset.

2. Structures in the "after" split date dataset are first filtered by sequence similarity to the "before" split date dataset. Sequence identity is calculated separately for the H, L, and antigen chain(s) using either a local or global alignment (-g, --globalSeqID flag). If the -cs (--cutoffStrict) flag is provided, structures are discarded if any of the sequence identity percentages are above the provided threshold (--abCompSeqCut). Otherwise, structures are discarded if the (maximum of the) sequence identity percentage(s) of the antigen(s) **and** the maximum of the sequence identity percentages of the H and L chains **are** above the provided threshold.

   The resulting dataset is further filtered by sequence similarity within the dataset (using the same procedure as above) with the --withinDatasetCut threshold. When one structure "knocks out" another from consideration due to high sequence similarity, the structure with the "fewest missing residues"<sup>+</sup> within the H, L, and antigen chain(s) is preferred. If both structures share the same number of "missing residues"<sup>+</sup>, the structure with the shorter antigen sequence is selected.

<sup>+</sup>: We use the minimum fraction of atom sequence length to SEQRES sequence length (as provided by the **PDB** file) as a measure of the number of missing residues. In particular, if this fraction is less than --minAtomSeqresFraction, the structure is discarded.

### Setup

To set up AADaM, please first set up Mosaist, **available on** GitHub at Grigoryanlab/Mosaist. Then, provide the path to Mosaist's `lib` directory in the second line of `src/utils.py`, replacing "/path/to/Mosaist/lib" with your path. The script also uses the `pandas` library.

If you want to use `conda`, follow these steps:

1. Download the Mosaist repository: https://github.com/Grigoryanlab/Mosaist
2. Create a conda environment named `AADaB_env` using `conda create --name AADaB_env`
3. Activate the conda environment using `conda activate AADaB_env`
4. Run `conda install conda-forge::boost`
5. Run `make` in the Mosaist directory
6. Run `make libs python` in the Mosaist directory
7. Run `conda install anaconda::pandas`

### Helper Scripts

Also included in the repo are helper scripts to IMGT number and otherwise clean up antibody-antigen structures, search antibody-antigen interfaces for structural motifs, check the Neff of both paired and single-chain MSAs, and calculate interfacial pLDDT. These scripts may be useful for analyzing antibody-antigen models in your future projects. Those scripts relying on Mosaist similarly require updating the Mosaist lib path.

### License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for more details.

### References

[1] McCoy KM, Ackerman ME, Grigoryan G. A comparison of antibody-antigen complex sequence-to-structure prediction methods and their systematic biases. Protein Sci. 2024 Sep;33(9):e5127. doi: 10.1002/pro.5127. PMID: 39167052; PMCID: PMC11337930.

[2] Constantin Schneider, Matthew I J Raybould, Charlotte M Deane. SAbDab in the age of biotherapeutics: updates including SAbDab-nano, the nanobody structure tracker. Nucleic Acids Research, Volume 50, Issue D1, 7 January 2022, Pages D1368–D1372, https://doi.org/10.1093/nar/gkab1050.

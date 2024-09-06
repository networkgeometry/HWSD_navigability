## Optimal navigability of weighted human brain connectomes in physical space

Laia Barjuan, Jordi Soriano and M. Ángeles Serrano<br>
NeuroImage (2024)<br>
DOI: [10.1016/j.neuroimage.2024.120703](https://doi.org/10.1016/j.neuroimage.2024.120703)<br>
arXiv: [arXiv:2311.10669](https://arxiv.org/abs/2311.10669)


Here we upload our source code to run the High Weight Short Distance navigation protocol for a given value of the time-out


### Code files in repository
* `grw_probab.f90`    High Weight Short Distance code
* `mtfort90_without_main.f90 & mtmod.mod`      Random number generation files

### Input files format

The code extracts the network information from two files with the following formatting
* `*_edges.txt`
```
[node i] [node j] [weight of connection between node i & node j]
```
* `*_coord_eucl.txt`
```
[node] [x coord] [y coord] [z coord]
```

### Running the code

The code requires a fortran 90 (or newer) compliant compiler

```
# Command line
gfortran grw_probab.f90 mtfort90_without_main.f90 -o grw_probab.out
./grw_probab.out
```

### Output files

The program outputs a single file containing the list of lambdas, the success rate, the average stretch, and the information cost of Euclidean and weight distances for each value of lambda


### Zenodo repository
_Optimal navigability of weighted human brain connectomes in physical space._<br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11198346.svg)](https://doi.org/10.5281/zenodo.11198346)

The repository comprises the weighted connectomes (edge list with fiber density) of two different datasets. The first dataset contains a subsample of 44 subjects extracted from the test-retest dataset included in the Human Connectome Project (HCP). The second dataset is compound by 40 healty controls selected from Lausanne psychosis cohort (UL). The age of the subjects belonging to the HCP dataset ranges between 22 and 35 years old while the mean age for the subjects included in the UL dataset is around 25 years old.

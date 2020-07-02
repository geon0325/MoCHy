# Hypergraph Motifs: Concepts, Algorithms, and Discoveries
Source code for the paper [Hypergraph Motifs: Concepts, Algorithms, and Discoveries](https://arxiv.org/abs/2003.01853), Geon Lee, Jihoon Ko, Kijung Shin, [VLDB 2020](https://vldb2020.org/).

We propose **Hypergraph Motifs (h-motifs)**, whose occurrences capture local structural patterns of real-world hypergraphs.

**H-motifs** describe connectivity patterns of three connected hyperedges with the following properties:
* *Exhaustive*: h-motifs capture connectivity patterns of all possible three connected hyperedges
* *Unique*: connectivity pattern of any three connected hyperedges is captured by exactly one h-motif
* *Size Independent*: h-motifs capture connectivity patterns independently of the sizes of hyperedges

**MoCHy** (**Mo**tif **C**ounting in **Hy**pergraphs) is a family of parallel algorithms for counting hypergraph motifs' instances.
* *MoCHy-E (MoCHy Exact)* exactly counts the instances of each h-motif.
* *MoCHy-A (MoCHy Approximate)*: approximately counts the instances of each h-motif.
* The advanced approximated version *MoCHy-A+* is up to 25X more accurate than *MoCHy-A*, and is up to 32X faster than *MoCHy-E*.

## Datasets
* The sample dataset is available from [here](https://gist.github.com/pszufe/02666497d2c138d1b2de5b7f67784d2b#sec_dblp).
* The real-world datasets used in the paper are from [here](https://www.cs.cornell.edu/~arb/data/).

## Running MoCHy
1. To run **MoCHy-E**, type 'run_exact.sh'.
2. To run *parallelized* **MoCHy-E**, type 'run_exact_par.sh'.
3. To run **MoCHy-A**, type 'run_approx_ver1.sh'.
4. To run **MoCHy-A+**, type 'run_approx_ver2.sh'.
5. To run *parallelized* **MoCHy-A+**, type 'run_approx_ver2_par.sh'.
6. To run *memory-bounded* **MoCHy-A+**, type 'run_approx_ver2_memory.sh'.

## Terms and Conditions
If you use this code as part of any published research, please acknowledge our VLDB 2020 paper.

## Contact Information
If you have any questions, please contact [Geon Lee](geonlee0325@kaist.ac.kr).

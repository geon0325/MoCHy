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
* The advanced approximated version *MoCHy-A+* is up to 25X more accurate than *MoCHy-A*, and it is up to 32X faster than *MoCHy-E*.

## Datasets
* The sample dataset is available [here](https://gist.github.com/pszufe/02666497d2c138d1b2de5b7f67784d2b#sec_dblp).
* The real-world datasets used in the paper are available [here](https://www.cs.cornell.edu/~arb/data/) or [here](http://dmlab.kaist.ac.kr/hmotif/).
* In the paper, we used datasets with unique hyperedges, where duplicated hyperedges are removed. 

## Input & Output Format
* The input format should be lines of hyperedges, where each line represents the nodes contained in each hyperedge.
* The index of the nodes should start from 0.
* For example, with 3 hyperedges: {0, 1, 2}, {2, 3}, and {1, 3, 4, 5}, the input file should be:
```
0,1,2
2,3
1,3,4,5
```
* The output of the code will be:
```
motif 1: 123
motif 2: 22
...
motif 26: 31
```

## Running Demo
You can run demo with the sample dataset (dblp_graph.txt).
1. To run **MoCHy-E**, type 'run_exact.sh'.
2. To run *parallelized* **MoCHy-E**, type 'run_exact_par.sh'.
3. To run **MoCHy-A**, type 'run_approx_ver1.sh'.
4. To run **MoCHy-A+**, type 'run_approx_ver2.sh'.
5. To run *parallelized* **MoCHy-A+**, type 'run_approx_ver2_par.sh'.
6. To run *memory-bounded* **MoCHy-A+**, type 'run_approx_ver2_memory.sh'.

## Terms and Conditions
If you use this code as part of any published research, please acknowledge our VLDB 2020 paper.
```
@article{lee2020hypergraph,
  title={Hypergraph Motifs: Concepts, Algorithms, and Discoveries},
  author={Lee, Geon and Ko, Jihoon and Shin, Kijung},
  journal={Proceedings of the VLDB Endowment},
  year={2020},
  publisher={VLDB Endowment}
}
```

## Contact Information
If you have any questions, please contact [Geon Lee (geonlee0325@kaist.ac.kr)](mailto:geonlee0325@kaist.ac.kr).

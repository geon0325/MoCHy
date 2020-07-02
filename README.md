# Hypergraph Motifs: Concepts, Algorithms, and Discoveries
Source code for the paper [Hypergraph Motifs: Concepts, Algorithms, and Discoveries](https://arxiv.org/abs/2003.01853), Geon Lee, Jihoon Ko, Kijung Shin, VLDB 2020.

We proposer **Hypergraph Motifs (h-motifs)**, whose occurrences capture local structural patterns of real-world hypergraphs.

**H-motifs** describe connectivity patterns of three connected hyperedges with the following properties:
* *Exhaustive*: h-motifs capture connectivity patterns of all possible three connected hyperedges
* *Unique*: connectivity pattern of any three connected hyperedges is captured by exactly one h-motif
* *Size Independent*: h-motifs capture connectivity patterns independently of the sizes of hyperedges

**MoCHy** (**Mo**tif **C**ounting in **Hy**pergraphs) is a family of parallel algorithms for counting hypergraph motifs' instances.
* *MoCHy-E (MoCHy Exact)* exactly counts the instances of each h-motif.
* *MoCHy-A (MoCHy Approximate)*: approximately counts the instances of each h-motif.
* The advanced approximated version *MoCHy-A^{+}*

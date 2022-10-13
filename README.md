Parallel-SSSP 
====================== 

This repository includes the implmentation of $\rho$-stepping, $\Delta$\*-stepping, and Bellman-Ford. 

## Developing 

### Prerequisites 
* g++ &gt;= 7 with support for Cilk Plus and C++17 (It is tested with g++ 7.5.0) 

### Setting up 
Clone the library with submodule 
```shell
git clone --recurse-submodules https://github.com/ucrparlay/Parallel-SSSP.git 
cd Parallel-SSSP/ 
```

Alternatively, you can first clone it and add the submodule 
```shell
git clone https://github.com/ucrparlay/Parallel-SSSP.git 
git submodule update --init --recursive 
cd Parallel-SSSP/ 
```

### Building
A makefile is given in the repository, you can compile the code by: 
```shell
make 
```

## Usage
```shell
./sssp [-i input_file] [-p parameter] [-w] [-s] [-v] [-a algorithm] 
```
Options: 
* -i input file path 
* -p parameter(e.g. delta, rho) 
* -w weighted input graph 
* -s symmetrized input graph 
* -v verify result 
* -a algorithm: [rho-stepping] [delta-stepping] [bellman-ford] 

For example, if you want to run $\rho$-stepping on a symmetrized weighted graph INPUT_NAME, set
$\rho$=2000000, and use Dijkstra's algorithm to verify the result after the test, you can run: 
```shell
./sssp -i INPUT_NAME -p 2000000 -w -s -v -a rho-stepping
```
## Graph Formats
The application can auto-detect the format of the input graph based on the suffix of the filename. Here is a list of supported graph formats: 
+ `.bin` The binary graph format from [GBBS](https://github.com/ParAlg/gbbs). 
+ `.adj` The adjacency graph format from [Problem Based Benchmark suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html). 
+ `.wsg` The weighted serialized pre-built graph format from [GAPBS](https://github.com/sbeamer/gapbs). 
+ `.gr` The galois graph file from [Galois](https://github.com/IntelligentSoftwareSystems/Galois). 

Some unweighted binary graphs can be found in our [Google Drive](https://drive.google.com/drive/u/3/folders/1ZuhfaLmdL-EyOiWYqZGD1rOy_oSFRWe4). For storage limit, we don't provide the large graphs used in our paper. They can be found in [Stanford Network Analysis Project](http://snap.stanford.edu/) and [Web Data Commons](http://webdatacommons.org/hyperlinkgraph/). 

## Reference 
X. Dong, Y. Gu, Y. Sun, and Y. Zhang. Efficient Stepping Algorithms and Implementations for Parallel Shortest Paths. In *ACM Symposium on Parallelism in Algorithms and Architectures (SPAA)*, 2021. 


X. Dong, Y. Gu, Y. Sun, and Y. Zhang. [Efficient Stepping Algorithms and Implementations for Parallel Shortest Paths.](https://arxiv.org/abs/2105.06145) *arXiv preprint 2105.06145*, 2021. 

Parallel-SSSP 
====================== 

This repository includes the implmentation of $\rho$-stepping, $\Delta$\*-stepping, and Bellman-Ford. 

## Developing 

### Prerequisites 
* g++ &gt;= 5.3.0 with support for Cilk Plus 

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

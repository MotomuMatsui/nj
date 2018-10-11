# nj

Reference: Matsui and Iwasaki, ???, 2018  
Online tool: [GS analysis server](http://gs.bs.s.u-tokyo.ac.jp/)  
Our Laboratory: [Iwasaki Lab](http://iwasakilab.bs.s.u-tokyo.ac.jp/eindex.html)  

[![Build Status](https://travis-ci.org/MotomuMatsui/nj.svg?branch=master)](https://travis-ci.org/MotomuMatsui/nj)

## History
version 1.0 (2018/10/12)   
  - Implemented in C++    

## Installation
### 0. Requirements

- [GNU GCC compiler](https://gcc.gnu.org/) (5.0+) is required to compile `nj`

### 1. Compile from source code:

````
    $ git clone https://github.com/MotomuMatsui/nj
    $ cd nj
    $ make
````

## Usage
To get on-line help:
```
    $ ./nj -h
```

The following command enables you to calculate NJ tree (phylogenetic tree reconstructed by Neighbor Joining method):
```
    $ ./nj [arguments] input > output
```

Arguments:

|Option| Description                                                                                         |
|:----:|:----------------------------------------------------------------------------------------------------|
|  -e  |<strong>[integer(>=0)]</strong> <em>The number of replicates for EP method. Default: 0</em>          |
|  -r  |<strong>[integer(>=1)]</strong> <em>The random seed number for EP method. Default: random number</em>|
|  -s  |<em>Silent mode: do not report progress. Default: Off</em>                                           |
|  -h  |<em>Show help messages. Default: Off</em>                                                            |
|  -v  |<em>Show the version. Default: Off</em>                                                              |


## License
This software is distributed under the GNU GPL, see [LICENSE](LICENSE)   
Copyright &copy; 2018, Motomu Matsui

## Author
[Motomu Matsui](https://sites.google.com/site/motomumatsui/)

## Reference
Naruya Saitou and Masatoshi Nei, The neighbor-joining method: a new method for reconstructing phylogenetic trees, Mol. Biol. Evol., 1987
Motomu Matsui and Wataru Iwasaki, ???, 2018

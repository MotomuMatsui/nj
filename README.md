# nj
`nj` is a software to conduct the naiive Neighbor Joining (NJ) method (Saito and Nei, MBE, 1987) with Edge Perturbation (EP) method (Matsui and Iwasaki, xxx, xxxx).    
`nj` is open-source software (GPL v3.0) implemented in C++ for <strong>Linux</strong>, <strong>Mac (macOS)</strong> and <strong>Windows (Cygwin)</strong>.    

Reference: Matsui and Iwasaki, ???, 2018  
Online tool: [GS analysis server](http://gs.bs.s.u-tokyo.ac.jp/)  
Our Laboratory: [Iwasaki Lab](http://iwasakilab.bs.s.u-tokyo.ac.jp/eindex.html)  

[![Build Status](https://travis-ci.org/MotomuMatsui/nj.svg?branch=master)](https://travis-ci.org/MotomuMatsui/nj)
[![Ubuntu](https://img.shields.io/badge/Linux-Ubuntu-green.svg)](https://www.ubuntu.com/)
[![CentOS](https://img.shields.io/badge/Linux-CentOS-green.svg)](https://www.centos.org/)
[![Mac](https://img.shields.io/badge/Mac-macOS-green.svg)](https://www.apple.com/macos/)
[![Windows](https://img.shields.io/badge/Windows-Cygwin-green.svg)](https://www.cygwin.com/)
[![Language](https://img.shields.io/badge/C%2B%2B-5.0%2B-green.svg)](https://gcc.gnu.org/)
[![GPL License](https://img.shields.io/badge/license-GPL-blue.svg)](LICENSE)

## History
version 2.0 (2018/11/16)   
  - Add Transfer Bootstrap Expectation algorithm (F. Lemoine, et al., Nature, 2018)    

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
|  -b  |<strong>[string(tbe/fbs)]</strong> <em>The bootstrap method. Default: tbe</em>                       |
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
Frederic Lemoine, Jean-Baka Domelevo Entfellner, Eduan Wilkinson, Damien Correia, Miraine Davila Felipe, Tulio De Oliveira, and Olivier Gascuel, Renewing Felsensteins phylogenetic bootstrap in the era of big data, Nature, 2018   
Motomu Matsui and Wataru Iwasaki, ???, 2018

/********************************************\
| Edge Perturbation Method v2.0 (2018/11/16) |
|                                            |
|  Copyright (c) 2015-2018 Motomu Matsui     |
|      Distributed under the GNU GPL         |
|                                            |
|      Matsui M and Iwasaki W (2018)         |
|      Systematic Biology, xx:xx-xx.         |
|                                            |
|      http://gs.bs.s.u-tokyo.ac.jp/         |
\********************************************/

#ifndef NJ_H
#define NJ_H

#include <iostream>
#include <string>

using namespace std;

int NJ(double* const (&oW), int* (&nj), int const& size);

#endif /*NJ_H*/

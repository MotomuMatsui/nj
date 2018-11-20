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

#ifndef EP_H
#define EP_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "ep_function.h"
#include "format.h"

int EP_fbs(double* const (&oW), unordered_map<string, double>& ep, function<double()>& R, int const& size);
int EP_tbe(double* const (&oW), int* const (&list_ori), unordered_map<string, double>& ep, function<double()>& R, int const& size);

#endif /*EP_H*/

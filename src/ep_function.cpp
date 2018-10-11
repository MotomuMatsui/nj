/******************************************\
| Graph Splitting Method v2.0 (2018/06/01) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <algorithm>
#include <cmath> 

using namespace std;

//Generalized Extreme Value function (inverse function)                                                                        
double gev(double const& x, double const& mu){
  double theta = mu*(1-mu)/3;
  double gamma = exp(-3*mu)-1;

  if(gamma == 0){ //Gummbel distribution                                                                                        
    return mu - theta*log(-log(x));
  }
  else{
    return mu +( pow(-log(x),-gamma)-1 )*theta/gamma;
  }
}

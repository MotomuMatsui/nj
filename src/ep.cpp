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

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <cmath>

using namespace std;

//ep_functions.cpp
extern double gev(double const&, double const&);

//format.cpp
extern void sc2list(int* const&, int*&, int const&);

int EP(double* const (&oW), unordered_map<string, double>& ep, function<double()>& R, int const& size){

  // Edge perturbation
  auto N = size*size;
  auto W = new double[N]();
  for(int n = 0; n < N; n++){
    auto a = oW[n];       // distance
    auto b = exp(-a);     // distance -> similarity    
    auto c = gev(R(), b); // Edge perturbation
    auto d = -1*log(c);   // similarity -> distance
    
    W[n] = (d<0)? 0: d;
    // d: Similarity scores perturved according to Generalized Extreme Value distribution
    // W[n]: Perturved distance matrix (0 <= W[n])
  }

  // Recursive Neighbor Joining (N-2 times)
  auto r    = new double[size](); // r_i = (1/(n-2))*SIGMA_{k}(d_ik)
  auto step = new int[N]();       // Result (e.g. n=3)
  for(int n = size; n>2; n--){

    auto m   = size-n;
    double t = 0.0;

    // r_i = (1/(n-2))*SIGMA_{k}(d_ik)
    for(int i=0; i<size; i++){
      r[i] = 0.0;
      for(int j=0; j<size; j++){
	if(W[i*size+j]>=0){
	  r[i] += W[i*size+j];
	  t    += W[i*size+j];
	}
      }
      r[i] /= size-2;
    }    

    // Search minimum D_ij = d_ij-r_i-r_j    
    auto Dm = 2.0*t;
    int im  = 0;
    int jm  = 0;

    for(int i=0; i<size-1; i++){
      for(int j=i+1; j<size; j++){
	if(W[i*size+j]>=0){

	  auto D = W[i*size+j]-r[i]-r[j];
	  
	  if(D<Dm){
	    Dm = D;
	    im = i;
	    jm = j;
	  }
	}
      }
    }

    // Result
    for(int i=0; i<size; i++){
      step[i*size+m+1] = step[i*size+m];
    }
    
    int gi = step[im*size+m];
    int gj = step[jm*size+m];

    if(gi==0){
      step[im*size+m+1] = m+1;
    }
    else{
      for(int i=0; i<size; i++){
	if(step[i*size+m]==gi){
	  step[i*size+m+1] = m+1;
	}
      }
    }
    
    if(gj==0){
      step[jm*size+m+1] = m+1;
    }
    else{
      for(int i=0; i<size; i++){
	if(step[i*size+m]==gj){
	  step[i*size+m+1] = m+1;
	}
      }
    }
  
    if(n > 2){
      // Re-calculate d_km = (d_im+d_jm-d_ij)/2
      // Replace im with k
      for(int i=0; i<size; i++){
	if(i != im && W[i*size+im] >= 0 && W[i*size+jm] >= 0){
	  W[i*size+im] = (W[i*size+im] + W[i*size+jm] - W[im*size+jm])/2;
	  W[i+size*im] = W[i*size+im];
	}
      }

      // Delete jm
      for(int i=0; i<size; i++){
	if(i != jm){
	  W[i*size+jm] = -1;
	  W[i+size*jm] = -1;
	}
      }
    }
  }

  for(int i=0; i<size; i++){
    step[(i+1)*size-1] = size-1;
  }

  // Transform
  int* nj = new int[size*size](); // Result (e.g. n=3)                                                                                 
  for(int i=0; i<size; i++){ // | 1, 1, 1 |    _|-1
    nj[i*size] = 1;          // | 1, 1, 3 |  _| --2
  }                          // | 1, 2, 2 |   L___3

  for(int i=1; i<size; i++){
    int m = size-i;
    int flag = -1;
    for(int j=0; j<size; j++){

      // Copy
      nj[j*size+i] = nj[j*size+i-1];

      if(step[(j+1)*size-i] == m){
        if(flag == -1){ // cluster (1st time)
          flag = step[(j+1)*size-i-1];
          nj[j*size+i] = i+1;
        }
        else if(flag == 0){}
        else if(flag == step[(j+1)*size-i-1]){ // cluster
          nj[j*size+i] = i+1;
        }
      }
    }
  }

  int* list_ep;
  sc2list(nj, list_ep, size);

  stringstream ss;
  for(int x=0; x<2*(size-3); x++){
    ss.str("");
    for(int z=0; z<size; z++){
      if(list_ep[x*size+z]>0){
        ss << z+1 << "|";
      }
    }
    ep[ss.str()] ++;
  }
  
  delete[] W;
  delete[] r;
  delete[] step;
  delete[] list_ep;
  return 1;
}

int EP2(double* const (&oW), int* const (&list), unordered_map<string, double>& ep, function<double()>& R, int const& size){

  // Edge perturbation
  auto N = size*size;
  auto W = new double[N]();
  for(int n = 0; n < N; n++){
    auto a = oW[n];       // distance
    auto b = exp(-a);     // distance -> similarity    
    auto c = gev(R(), b); // Edge perturbation
    auto d = -1*log(c);   // similarity -> distance
    
    W[n] = (d<=0)? 0: d;
      // d: Similarity scores perturved according to Generalized Extreme Value distribution
      // W[n]: Perturved distance matrix (0 <= W[n])
  }
    
  // Recursive Neighbor Joining (N-2 times)
  auto r    = new double[size]();   // r_i = (1/(n-2))*SIGMA_{k}(d_ik)
  auto step = new int[size*size](); // Result (e.g. n=3)
  for(int n = size; n>2; n--){

    int m    = size-n;
    double t = 0.0;

    // r_i = (1/(n-2))*SIGMA_{k}(d_ik)
    for(int i=0; i<size; i++){
      r[i] = 0.0;
      for(int j=0; j<size; j++){
	if(W[i*size+j]>=0){
	  r[i] += W[i*size+j];
	  t    += W[i*size+j];
	}
      }
      r[i] /= size-2;
    }    

    // Search minimum D_ij = d_ij-r_i-r_j    
    auto Dm = 2.0*t;
    int im  = 0;
    int jm  = 0;

    for(int i=0; i<size-1; i++){
      for(int j=i+1; j<size; j++){
	if(W[i*size+j]>=0){

	  auto D = W[i*size+j]-r[i]-r[j];
	  
	  if(D<Dm){
	    Dm = D;
	    im = i;
	    jm = j;
	  }
	}
      }
    }

    // Result
    for(int i=0; i<size; i++){
      step[i*size+m+1] = step[i*size+m];
    }
    
    int gi = step[im*size+m];
    int gj = step[jm*size+m];

    if(gi==0){
      step[im*size+m+1] = m+1;
    }
    else{
      for(int i=0; i<size; i++){
	if(step[i*size+m]==gi){
	  step[i*size+m+1] = m+1;
	}
      }
    }
    
    if(gj==0){
      step[jm*size+m+1] = m+1;
    }
    else{
      for(int i=0; i<size; i++){
	if(step[i*size+m]==gj){
	  step[i*size+m+1] = m+1;
	}
      }
    }
  
    if(n > 2){
      // Re-calculate d_km = (d_im+d_jm-d_ij)/2
      // Replace im with k
      for(int i=0; i<size; i++){
	if(i != im && W[i*size+im] >= 0 && W[i*size+jm] >= 0){
	  W[i*size+im] = (W[i*size+im] + W[i*size+jm] - W[im*size+jm])/2;
	  W[i+size*im] = W[i*size+im];
	}
      }

      // Delete jm
      for(int i=0; i<size; i++){
	if(i != jm){
	  W[i*size+jm] = -1;
	  W[i+size*jm] = -1;
	}
      }
    }
  }

  for(int i=0; i<size; i++){
    step[(i+1)*size-1] = size-1;
  }

  // Transform
  int* nj = new int[size*size](); // Result (e.g. n=3)
  for(int i=0; i<size; i++){ // | 1, 1, 1 |    _|-1
    nj[i*size] = 1;          // | 1, 1, 3 |  _| --2
  }                          // | 1, 2, 2 |   L___3

  for(int i=1; i<size; i++){
    int m = size-i;
    int flag = -1;
    for(int j=0; j<size; j++){
      
      // Copy
      nj[j*size+i] = nj[j*size+i-1];

      if(step[(j+1)*size-i] == m){
        if(flag == -1){ // cluster (1st time)
          flag = step[(j+1)*size-i-1];
          nj[j*size+i] = i+1;
        }
        else if(flag == 0){}
        else if(flag == step[(j+1)*size-i-1]){ // cluster
          nj[j*size+i] = i+1;
        }
      }
    }
  }
  
  int* list_ep;
  sc2list(nj, list_ep, size);

  for(int x=0; x<size-3; x++){
    double min_share = size;

    // calculate subset sizes
    double num_a1 = 0;
    double num_a2 = 0;
    for(int z=0; z<size; z++){
      num_a1 += list[2*x*size+z];
      num_a2 += list[(2*x+1)*size+z];      
    }

    auto num_a = (num_a1 <= num_a2)? num_a1: num_a2;

    for(int y=0; y<2*(size-3); y++){

      double num_b = 0;
      for(int z=0; z<size; z++){
        num_b += list_ep[y*size+z];
      }

      double share1 = 0;
      double share2 = 0;
      for(int z=0; z<size; z++){
	share1 += list[2*x*size+z]^list_ep[y*size+z];
	share2 += list[(2*x+1)*size+z]^list_ep[y*size+z];
      }
      min_share = (share1 < min_share)? share1: min_share;
      min_share = (share2 < min_share)? share2: min_share;
    }

    stringstream ss1;
    stringstream ss2;
    for(int z=0; z<size; z++){
      if(list[2*x*size+z]>0){
	ss1 << z+1 << "|";	
      }
      else{
	ss2 << z+1 << "|";	
      }
    }
    if(min_share < num_a - 1){
      ep[ss1.str()] += 1 - min_share / (num_a - 1);
      ep[ss2.str()] += 1 - min_share / (num_a - 1);
    }
  }
  
  delete[] W;
  delete[] r;
  delete[] step;
  delete[] nj;
  delete[] list_ep;
  return 1;
}

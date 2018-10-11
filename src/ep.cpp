/********************************************\
| Edge Perturbation Method v1.0 (2018/10/12) |
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

using namespace std;

//ep_functions.cpp
extern double gev(double const&, double const&);

int EP(double* const (&oW), unordered_map<string, double>& ep, function<double()>& R, int const& size){

  // Edge perturbation
  int N = size*size;
  double* W = new double[N]();
  for(int n = 0; n < N; n++){
    double a = oW[n];
    if(a==0){
      W[n] = 0; // Graph topology is not changed
    }
    else{
      auto b = gev(R(), a);
      W[n] = 
        (b>1)? 1: // b: Similarity scores perturved according to Generalized Extreme Value distribution
        (b<0)? 0: // E[n]: Perturved sequence similarity score
	b; // 0 <= E[n] <= 1
    }
  }
  
  auto r    = new double[size]();      // r_i = (1/(n-2))*SIGMA_{k}(d_ik)
  auto step = new int[size*size]();    // Result (e.g. n=3)
  for(int i=0; i<size; i++){           // | 0, 1, 2 |    _|-1
    step[i*size] = 0;                  // | 0, 1, 1 |  _| --2
  }                                    // | 0, 0, 0 |   L___3

  // Recursive Neighbor Joining (N-2 times)
  for(int n = size; n>2; n--){

    int m  = size-n;
    int im = 0;
    int jm = 0;
    double t  = 0.0;

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
    double Dm = 2.0*t;
    double D  = 0.0;
    for(int i=0; i<size-1; i++){
      for(int j=i+1; j<size; j++){
	if(W[i*size+j]>=0){

	  D = W[i*size+j]-r[i]-r[j];
	  
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
  stringstream ss;
  
  for(int i=1; i<size; i++){
    int m = size-i;
    int flag = -1;

    vector<int> a;
    vector<int> b;

    for(int j=0; j<size; j++){
      if(step[(j+1)*size-i] == m){
	if(flag == -1){ // cluster (1st time)
	  flag = step[(j+1)*size-i-1];
	  a.push_back(j);
	}
	else if(flag == 0){
	  b.push_back(j);
	}
	else if(flag == step[(j+1)*size-i-1]){ // cluster
	  a.push_back(j);
	}
	else{
	  b.push_back(j);	  
	}
      }
    }
    
    ss.str("");
    sort(a.begin(), a.end());
    for(int n : a){
      ss << n+1 << "|";
    }
    ep[ss.str()]++;
    
    ss.str("");
    sort(b.begin(), b.end());
    for(int n : b){
      ss << n+1 << "|";
    }
    ep[ss.str()]++;
  }
  
  delete[] W;
  delete[] r;
  delete[] step;
  return 1;
}

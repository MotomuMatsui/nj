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

#include "nj.h"

using namespace std;

int NJ(double* const (&oW), int* (&nj), int const& size){
    
  // Copy (oW -> W) & Check
  if(size < 2){
    return -1;
  }
  else if(size == 2){
    nj = new int[4](); // Result (e.g. n=2)
    nj[0] = 1; nj[1] = 1;
    nj[2] = 1; nj[3] = 2;

    return 1;
  }

  auto W = new double[size*size](); // Distance matrix
  for(int i=0; i<size-1; i++){
    for(int j=i+1; j<size; j++){
      auto v = oW[i*size+j];
      if(v < 0){
	delete[] W;
	return -1;
      }
      W[i*size+j] = v;
      W[i+size*j] = v;
    }
    W[i*size+i] = 0.0;
  }
  
  // Recursive Neighbor Joining (N-2 times)
  auto r    = new double[size]();     // r_i = (1/(n-2))*SIGMA_{k}(d_ik)
  auto L    = new double[size*2-3](); // branch lengths of the NJ tree
  auto step = new int[size*size]();   // Result (e.g. n=3)
  for(int n = size; n>2; n--){

    int m  = size-n;
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
  
    // Bramch length
    L[2*m]   = (W[im*size+jm]+r[im]-r[jm])/2; // d_ik
    L[2*m+1] = (W[im*size+jm]-r[im]+r[jm])/2; // d_jk
    
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

  // Hope (last pair)
  int li = 0;
  int lj = 0;
  for(int i=0; i<size-1; i++){
    for(int j=i+1; j<size; j++){
      if(W[i*size+j]>=0){
	li = i;
	lj = j;

	//break
	i = size;
	j = size;
      }
    }
  }
  
  L[2*(size-2)] = W[li*size+lj];

  for(int i=0; i<size; i++){
    step[(i+1)*size-1] = size-1;
  }

  // Transform
  nj = new int[size*size](); // Result (e.g. n=3)
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
  
  delete[] W;
  delete[] r;
  delete[] L;
  delete[] step;
  return 1;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

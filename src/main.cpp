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

#include <random>   
#include <regex>
#include <unistd.h>

#include "ep.h"
#include "format.h"
#include "messages.h"
#include "nj.h"

using namespace std;

int main(int argc, char* argv[]){

  /*/Getopt/*/
  int    silence   = 0;         // -s
  int    ep_num    = 0;         // -e
  int    seed      = 0;         // -r
  string bs_method = "";        // -b [fbp (Felsenstein's bootstrap proportion) or tbe (transfer bootstrap expectation)]

  opterr = 0; // default error messages -> OFF
  int opt;
  regex renum(R"(^[\d\.]+$)"); // -e/-r option requires an integer/flout number
  while ((opt = getopt(argc, argv, "shve:r:b:")) != -1){
    if(opt == 'e'){ // OK! (./gs -e 100 IN.fst)
      if(regex_match(optarg, renum)){
        ep_num = atoi(optarg);
      }
      else{ // NG! (./gs -e hundred IN.fst)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
    else if(opt == 'r'){ // OK! (./gs -r 12345 IN.fst)
      if(regex_match(optarg, renum)){
        seed = atoi(optarg);
      }
      else{ // NG! (./gs -r one_two_three IN.fst)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -r requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
    else if(opt == 'h'){ // HELP message (./gs -h)
      /*PRINT*/ print_banner();
      /*PRINT*/ print_usage(argv[0]);
      return 0;
    }
    else if(opt == 'v'){ // Version (./gs -v)
      /*PRINT*/ print_banner();
      return 0;
    }
    else if(opt == 's'){ // SILENT mode (./gs -s -e 100 IN.fst)
      silence = 1;
    }
    else if(opt == 'b'){ // Statistical method to assess the robustness of inffered branches (./gs -e 100 -b tbe IN.fst)
      bs_method = optarg;
      if(bs_method != "fbs" && bs_method != "tbe"){
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -b requires a string that equals 'fbs' or 'tbe'.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
    else if (opt == '?'){
      if(optopt == 'e'){ // NG! (./gs IN.fst -e)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
      else if(optopt == 'r'){ // NG! (./gs IN.fst -r)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -r requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
      else{ // NG! (./gs -Z)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << argv[0] << ": invalid option\n" <<  endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
  }  
  
  /*/Input file/*/
  string input = "";
  if(optind < argc){ // OK!
    input = argv[optind];
  }
  else{ // NG!
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << argv[0] << " requires an input file (matrix).\n" << endl;
    /*PRINT*/ print_usage(argv[0]);
    return -1;
  }

  /*/Variables/*/
  double* W;     // Distance matrix
  int size;      // Row size of W (W is a synmetry matrix)
  int* nj;       // Result of NJ method
  string newick; // NJ tree (without EP values)

  /*/File I/O/*/
  auto original_file = string(input);
  ifstream ifs1(original_file); // Matrix file (original)

  if(ifs1.fail()){
    /*PRINT*/ cerr << "\nCannot access " << original_file << "!" << endl;
    return -1;
  }

  /*/Parsing matrix file/*/
  auto same_sequence = readMAT(ifs1, W, size); 
    // ifs1: INPUT (original matrix file)
    // W:    OUTPUT (matrix)
    // size: # of sequence = row size of sequence similarity matrix
  
  /*/Parameters/*/
  if(!silence){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "Settings:" << endl;
    /*PRINT*/ cerr << "-Input" << endl;
    /*PRINT*/ cerr << "  file = " << input << endl;
    /*PRINT*/ cerr << "  # of sequences = " << size << endl << endl;
    /*PRINT*/ if(same_sequence>0) cerr << "  <WARNING> This dataset has " << same_sequence << " duplicated sequence pair(s)" << endl << endl;
    /*PRINT*/ cerr << "-EP method" << endl;

    if(ep_num>0){
      if(bs_method == "fbs"){
	/*PRINT*/ cerr << "  method = Felsenstein's bootstrap proportion" << endl;      
      }
      else{
	/*PRINT*/ cerr << "  method = Transfer bootstrap expectation" << endl;            
      }
    }

    if(seed>0){
      /*PRINT*/ cerr << "  random seed = " << seed << endl;
    }
    else{
      /*PRINT*/ cerr << "  random seed = " << "a random number (default)" << endl;
    }

    /*PRINT*/ cerr << "  # of iterations = " << ep_num << endl << endl;
    /*PRINT*/ cerr << "Progress:" << endl;
  }

  /*/NJ method/*/
  /*PRINT*/ if(!silence) cerr << "-NJ method\n" << "  executing...\r" << flush;
  NJ(W, nj, size);
    // W: INPUT (sequence similarity matrix)
    // nj: OUTPUT (result of stepwise neighbor joining method)

  /*PRINT*/ if(!silence) cerr << "  done.         " << endl << endl;

  /*/Generating NJ tree Newick/*/
  sc2nwk(nj, newick, size);
    // nj: INPUT (result of stepwise spectral clustering)
    // newick: OUTPUT (GS tree [newick format])

  /*/EP method/*/
  if(ep_num>0){
    /*PRINT*/ if(!silence) cerr << "-EP method" << endl;

    string newick_EP; // GS+EP tree
    unordered_map<string, double> ep;

    if(bs_method == "fbs"){
      // Random number generator (Uniform distribution->Mersenne Twister)
      function<double()> R;
      uniform_real_distribution<double> urd(0,1);    // uniform distributed random number
      
      if(seed>0){
	mt19937 mt(static_cast<unsigned int>(seed)); // mersenne twister
	R = bind(urd, ref(mt));                      // random number generator    
      }
      else{
	random_device rd;                            // random seed
	mt19937 mt(rd());                            // mersenne twister
	R = bind(urd, ref(mt));                      // random number generator        
      }    

      for(int n=1; n<=ep_num; n++){
	/*PRINT*/ if(!silence) cerr << "  " << n << "/" << ep_num << " iterations" << "\r"<< flush;

	EP_fbs(W, ep, R, size);
          // W: INPUT (sequence similarity matrix)
          // ep: OUTPUT (result of Edge Perturbation method)
          // R: random number generator
      }      
    }
    else{
      int* list_ori; 
      sc2list(nj, list_ori, size);
        // nj: INPUT (result of stepwise spectral clustering)
        // list_ori: OUTPUT (NJ tree [leaves])      
      
      // Random number generator (Uniform distribution->Mersenne Twister)
      function<double()> R;
      uniform_real_distribution<double> urd(0,1);    // uniform distributed random number
      
      if(seed>0){
	mt19937 mt(static_cast<unsigned int>(seed)); // mersenne twister
	R = bind(urd, ref(mt));                      // random number generator    
      }
      else{
	random_device rd;                            // random seed
	mt19937 mt(rd());                            // mersenne twister
	R = bind(urd, ref(mt));                      // random number generator        
      }    

      for(int n=1; n<=ep_num; n++){
	/*PRINT*/ if(!silence) cerr << "  " << n << "/" << ep_num << " iterations" << "\r"<< flush;
	
	EP_tbe(W, list_ori, ep, R, size);
          // W: INPUT (sequence similarity matrix)
          // ep: OUTPUT (result of Edge Perturbation method)
          // R: random number generator
      }
      delete[] list_ori;
    }

    /*PRINT*/ if(!silence) cerr << "\n  done." << endl << endl;
    /*PRINT*/ if(!silence) cerr << "------------------------------------------\n" << endl;
    
    addEP(newick, newick_EP, ep, ep_num, size);
      // newick: INPUT (GS tree [newick format])
      // newick_EP: OUTPUT (GS+EP tree [newick format])
      // ep: INPUT (result of Edge Perturbation method)
      // ep_num: INPUT (# of Edge Perturbation method)
    
    /*/ GS tree WITH EP values ->STDOUT /*/
    cout << newick_EP << endl;
  }
  else{ // skip the EP method
    /*PRINT*/ if(!silence) cerr << "------------------------------------------\n" << endl;

    /*/ GS tree WITHOUT EP values ->STDOUT /*/
    cout << newick << endl;
  }
  
  return 0;
}

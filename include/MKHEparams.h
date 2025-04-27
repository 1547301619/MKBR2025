#ifndef MKPARAMS
#define MKPARAMS

#include <NTL/ZZX.h>
#include <iostream>
#include <vector>
#include <cmath> 

using namespace std;

struct MKHEParams{
        const static int q_static=32749;   
        const static int Q_static=133919213;  
        const static int N_static=2048;
        const static int N2p1_static=N_static/2+1;

        int q;  //LWE modulus   
        int Q;  //RLWE,NTRU modulus
        int n=500;   //LWE dimension
        int N=2048;   //dimension of RQ
        int N2p1=N/2+1;
        int parties=4;  //number of parties

        int B_n=512;  //gadget base
        int d_n=3;    //gadget length
        int shift_n=9;  //log2B

        int B_r=32;
        int d_r=6;
        int shift_r=5;

        int B_l=32;
        int d_l=4;
        int shift_l=5;

        double stdev_LWEerr=1.9;  
        double stdev_RLWEkey=0.25;
        double stdev_RLWEerr=0.25;
        double stdev_NTRUerr=0.25;
         

        int t=4;
        int half_delta_base;
        int delta_base;
        int nand_const;
        int half_q_base;
        
        MKHEParams();
        MKHEParams(int parties,int q,int Q,int n,int N,int B_n,int B_r,int B_l);


        vector<long> Q_list;
};

enum BINFHE_PARAMSET {
      MKHE2party_v1,
      MKHE2party_v2,
      MKHE4party_v2,
      MKHE8party_v2,
      MKHE16party_v2,
};

//ParamsSet (k,q,Q,n,N,B_n,B_r,B_l)
const MKHEParams PARAM_MKHE2party_v1(2,32749,133919213,500,2048,512,4,16);  
const MKHEParams PARAM_MKHE2party_v2(2,32749,133919213,500,2048,512,32,16);  

const MKHEParams PARAM_MKHE4party_v2(4,32749,133919213,500,2048,512,32,8); 

const MKHEParams PARAM_MKHE8party_v2(8,32749,133919213,500,2048,128,32,8);  

const MKHEParams PARAM_MKHE16party_v2(16,32749,133919213,500,2048,128,32,8);  

const int Q=MKHEParams::Q_static;
const int Q_half=Q/2;
const long Q_long=MKHEParams::Q_static;
const long Q_long_half=Q_long/2;


#endif
#include "MKHEparams.h"

MKHEParams::MKHEParams(){};
MKHEParams::MKHEParams(int parties,int q,int Q,int n,int N,int B_n,int B_r,int B_l) {
    this->q = q;
    this->Q = Q;
    this->n = n;
    this->N = N;
    this->parties = parties;
    this->B_n = B_n;
    this->B_r = B_r;
    this->B_l = B_l;
    this->shift_n=log2(B_n);
    this->shift_r=log2(B_r);
    this->shift_l=log2(B_l);
    this->d_n=int(ceil(log(double(Q))/log(double(B_n))));
    this->d_r=int(ceil(log(double(Q))/log(double(B_r))));
    this->d_l=int(ceil(log(double(q))/log(double(B_l))));
    this->half_delta_base = q/(2*t);
    this->nand_const = 5*half_delta_base;
    this->delta_base=2*half_delta_base;
    this->half_q_base=q/2;
}





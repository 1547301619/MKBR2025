#ifndef POLY
#define POLY

#include <NTL/ZZX.h>
#include "FINAL.h"
#include <iostream>
#include <cmath> 
#include "MKHEparams.h"


void external_product(vector<long>& res, const vector<int>& poly, const vector<FFTPoly>& poly_vector, int b, int shift, int l,int N,int N2p1);

void mult_fft_poly_by_int(FFTPoly& a, const int b);

template <typename T>
inline void printpoly(vector<T> p){
    for (int i = 0; i < p.size(); ++i) {
        cout << p[i] << " ";
    }
    cout<<endl<<endl;
}

inline void modulo_switch_a(vector<int>& a, int old_q, int new_q,int length)
{
    for (size_t i = 0; i < length; i++)
        a[i] = static_cast<int>((static_cast<long long>(a[i]) * new_q) / old_q);
}

void modulo_switch_poly(ModQPoly& v,int old_q,int new_q);

void polymul(vector<int>& res,const vector<int>& p1,const vector<int>& p2,long q);


void get_uniform_vector(vector<int>& vec,int Q);

void get_gaussian_vector(vector<int>& vec, double st_dev);

void get_ternary_vector(vector<int>& vec);

void get_binary_vector(vector<int>& vec);

void get_hwt_vector(vector<int>& vec, int h);

void gadget_decomp(vector<FFTPoly>& res_vector,const vector<int>& poly,int b, int shift, int l);

void get_gadget(vector<FFTPoly> &g,int N,int B,int d);

int mod_q_poly(const int input);

int mod_q_poly(const int input,int Q);

void mod_q_poly(ModQPoly& v,int Q);

void mod_q_poly(ModQPoly& v);

void mod_q_poly(ModQPoly& v,vector<long> tmp_long,int Q);

void mod_q_poly(ModQPoly& v,vector<long> tmp_long);

int mod_q_lwe(long input,int q,int q_half);

int mod_q_lwe(int input,int q,int q_half);


void get_invertible_vector(vector<int>& vec, vector<int>& vec_inv, int scale, int shift);

void get_invertible_vector(vector<int>& vec, vector<int>& vec_inv, int scale,int N,int Q);

void get_invertible_vector_gaussic(vector<int>& vec, vector<int>& vec_inv,double st_dev);

void get_invertible_poly(vector<int>& vec_inv, vector<int>& vec,int N,int Q);


int errestimator(const ModQPoly& v,int std);

void print_hamin(vector<int> vec);








#endif
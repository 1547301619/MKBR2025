#ifndef MKNTRUHE
#define MKNTRUHE
#include "poly.h"

typedef struct{
    ModQPoly sk;
    ModQPoly sk_inv;
}SKey_boot;

void NTRU_scarlar_enc(ModQPoly &res,const ModQPoly& m,const SKey_boot& sk,int Q,int N);

void NTRU_scarlar_dec(ModQPoly &m_dec,const ModQPoly& c,const SKey_boot& sk,int Q,int N);

void NTRU_vector_enc(NGSFFTctxt& ct, const ModQPoly& m, int l, int B, const SKey_boot& sk_boot,MKHEParams *parmk);

void NTRU_vector_enc(NGSFFTctxt& ct, int m, int l, int B, const SKey_boot& sk_boot,MKHEParams *parmk);




#endif
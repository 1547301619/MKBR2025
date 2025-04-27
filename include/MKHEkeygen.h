#ifndef MKKEYGEN
#define MKKEYGEN

#include "MKHEparams.h"
#include "poly.h"
#include "NTRUhe.h"
#include <iostream>

typedef struct{
    ModQMatrix A;
    std::vector<int> b;
}  KSKey_LWE;

typedef std::vector<int> SKey_base_LWE;
typedef std::vector<int> MKLweKey;
typedef vector<vector<FFTPoly>> MKRLwePK;
typedef vector<ModQPoly> ACCList;
typedef vector<vector<NGSFFTctxt>> MKBRKey;

typedef struct MKRLweKey{
    vector<FFTPoly> MKRLweKey_SK; // (z0,...zk) z0=1
    vector<ModQPoly> MKRLweKey_SK_poly;   //used in extractkey
    vector<vector<FFTPoly>> MKRLweKey_PK; // (b0,b1,...,bk,a) bj=-zj*a +e
}MKRLweKey;

typedef struct HPK{
    vector<vector<FFTPoly>> hpk_FFT;//d0,d1,d2
    int party; //party index in 1,..k
}HPK;


typedef vector<vector<KSKey_LWE>> MKlksk;

typedef vector<vector<vector<FFTPoly>>> RGSWKey;

typedef struct MKnrksk{
    vector<vector<FFTPoly>> nrksk0;
    vector<RGSWKey> nrksk1;
}MKnrksk;


void HybridProduct_rksk(vector<ModQPoly>& c,const HPK& hpk,const MKRLweKey& mkrlwe_key,MKHEParams *parmk);

void MKLweKeyGen(MKLweKey &result,MKHEParams *parmk);

void MKRLweKeyGen(MKRLweKey &result,MKHEParams *parmk);

void HPKGen(HPK& hpk,const MKRLweKey &mkrlwe_key,const ModQPoly& f,int party,MKHEParams *parmk);


void MKrkskGen(vector<vector<FFTPoly>> &mkrksk,const vector<HPK>& mkhpk,const MKRLweKey &mkrlwe_key,MKHEParams *parmk);

void MKHPKGen(vector<HPK>& mkhpk,const MKRLweKey &mkrlwe_key,const vector<SKey_boot>& vec_sk_boot,MKHEParams *parmk);

void ReKeyGen(vector<FFTPoly> &ReKey,const SKey_boot& sk_boot,MKHEParams *parmk);

void MKReKeyGen(vector<vector<FFTPoly>> &mkReKey,const vector<SKey_boot>& vec_sk_boot,MKHEParams *parmk);

void NRKGen(vector<vector<FFTPoly>> &nrksk0,vector<RGSWKey> &nrksk1,const vector<SKey_boot>& vec_sk_boot,const MKRLweKey &mkrlwe_key,MKHEParams *parmk);

void BRKGen(vector<NGSFFTctxt> &BRK,const SKey_base_LWE& sk_base, const SKey_boot& sk_boot,MKHEParams *parmk);

void MKBRKGen(MKBRKey &mkBRK,const MKLweKey& mklwe_sk,const vector<SKey_boot>& vec_sk_boot,MKHEParams *parmk);
 
void MKBOOTSKGen(vector<SKey_boot> &vec_sk_boot,MKHEParams *parmk);

void LKSKGen(vector<KSKey_LWE> &lksk,const vector<int>& sk_in,const vector<int>& sk_out,MKHEParams *parmk);

void MKLKSKGen(MKlksk &mklksk,const MKLweKey& sk_in,const MKLweKey& sk_out,MKHEParams *parmk);

#endif
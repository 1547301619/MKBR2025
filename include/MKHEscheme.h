#ifndef MKHESCHEME
#define MKHESCHEME

#include "MKLwe.h"
#include "NTRUhe.h"
#include "RLWEhe.h"


class MKHEscheme{
    public:
    MKHEParams parmk;
    MKLweSample MKLwe;
    MKLweKey mklwe_sk;
    MKRLweKey mkrlwe_key;
    vector<SKey_boot> vec_sk_boot;
    MKBRKey mkBRK;
    vector<HPK> mkhpk;
    vector<vector<FFTPoly>> mkReKey; //brk_n
    vector<vector<FFTPoly>> mkrksk;
    vector<vector<FFTPoly>> nrk0; //v1
    vector<RGSWKey> nrk1;  //V2...Vk
    MKLweKey mklwe_sk_z;
    MKlksk mklksk;
    

    MKHEscheme(BINFHE_PARAMSET set,int v);

    void MKHEencrypt(MKLweSample &result,const MKLweKey& key,int m,MKHEParams *parmk);

    int MKHEDecrypt(const MKLweSample &ctxt,const MKLweKey& key,MKHEParams *parmk);


};



void MKBootstrap_v1(MKLweSample &result,const MKLweSample &ct,const MKBRKey& mkBRK,const vector<vector<FFTPoly>> &mkReKey,const vector<vector<FFTPoly>> &mkrksk,const MKlksk& mklksk,MKHEParams *parmk);

void MKBootstrap_v2(MKLweSample &result,const MKLweSample &ct,const MKBRKey& mkBRK,const vector<vector<FFTPoly>> &mkreKey,const vector<vector<FFTPoly>>& nrk0,const vector<RGSWKey> &nrk1,const MKlksk& mklksk,MKHEParams *parmk);

void SingleKeyReBR(ModQPoly &acc, int *al,const ModQPoly &c,const NGSFFTctxt &ReKey,const vector<NGSFFTctxt>& BRK,MKHEParams *parmk);

void MKkeyswitch(MKLweSample &result,const MKLweSample &input,const MKlksk& mklksk,MKHEParams *parmk);

void MKkeyswitch_v2(MKLweSample &result,const MKLweSample &input,const MKlksk& mklksk,MKHEParams *parmk);

void MKExtract(MKLweSample &result,const vector<ModQPoly>& mkrlwe_c,MKHEParams *parmk);

void RlweExtract(MKLweSample& result,const vector<ModQPoly>& rlwe_c,MKHEParams *parmk);

void MKExtractKey(MKLweKey& mklwe_sk_z,const MKRLweKey& mkrlwe_key,MKHEParams *parmk);

void MKModSwitch(MKLweSample &result,int oldQ,int newQ,MKHEParams *parmk);

void MKModSwitch(MKLweSample &result,MKHEParams *parmk);

void LWEModSwitch(MKLweSample &result,int oldQ,int newQ,MKHEParams *parmk);

void MKRlweDec(ModQPoly &res,const vector<ModQPoly>& c,const vector<FFTPoly>& z,MKHEParams *parmk,int Q);



#endif
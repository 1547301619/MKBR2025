#ifndef MKRLWEHE
#define MKRLWEHE
#include "MKHEkeygen.h"


void RLWE_Enc(vector<ModQPoly> &rlwe_c,const ModQPoly &m,const ModQPoly &z);

void RLWE_Dec(ModQPoly &res,const vector<ModQPoly> &rlwe_c,const ModQPoly &z);

void RLWE_RGSWproduct(vector<ModQPoly>& res,const vector<ModQPoly> &rlwe_c,const RGSWKey &RGSW_c,MKHEParams *parmk);









#endif
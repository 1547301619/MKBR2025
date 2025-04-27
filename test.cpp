#include <iostream>
#include <cassert>
#include <time.h>
#include <cstdint>
#include <stdexcept>
#include <chrono>
#include <limits.h>

#include <NTL/ZZX.h>
#include <NTL/ZZ_pE.h>

#include "FINAL.h"

#include <cmath> 
#include "MKHEscheme.h"


using namespace std;
using namespace NTL;


enum GateType {NAND, AND, OR, XOR, NOT};




void test_nandbootstrap_v1(BINFHE_PARAMSET set){
    MKHEscheme MKhe(set,1);
    int n=MKhe.parmk.n;
    int k=MKhe.parmk.parties;
    int q=MKhe.parmk.q;
    int nand_const=MKhe.parmk.nand_const;

    int m1=binary_sampler(rand_engine);
    int m2=binary_sampler(rand_engine);
    MKLweSample c1(q,n,k);
    MKLweSample c2(q,n,k);

    MKhe.MKHEencrypt(c1,MKhe.mklwe_sk,m1,&MKhe.parmk);
    MKhe.MKHEencrypt(c2,MKhe.mklwe_sk,m2,&MKhe.parmk);

    MKLweSample c_nand(q,n,k);
    MKLweSample c_nandboot(q,n,k);


 

    clock_t start = clock();
    c_nand=nand_const-(c1+c2);
    MKBootstrap_v1(c_nandboot,c_nand,MKhe.mkBRK,MKhe.mkReKey,MKhe.mkrksk,MKhe.mklksk,&MKhe.parmk); 
    cout<<"NAND_Bootstrap_v1:"<<float(clock()-start)* 1000/CLOCKS_PER_SEC<<"ms"<<endl;
    

    int m_nand=MKhe.MKHEDecrypt(c_nand,MKhe.mklwe_sk,&MKhe.parmk);
    int m_boot=MKhe.MKHEDecrypt(c_nandboot,MKhe.mklwe_sk,&MKhe.parmk);

    //correct test
    if(m_boot == !(m1 & m2)){
        cout<<"correct"<<endl;
    }
    else{
        cout<<"error"<<endl;
    }
}

void test_nandbootstrap_v2(BINFHE_PARAMSET set){

    MKHEscheme MKhe(set,2);
    int n=MKhe.parmk.n;
    int k=MKhe.parmk.parties;
    int q=MKhe.parmk.q;
    int m1=1;
    int m2=0;
    MKLweSample c1(q,n,k);
    MKLweSample c2(q,n,k);

    MKhe.MKHEencrypt(c1,MKhe.mklwe_sk,m1,&MKhe.parmk);
    MKhe.MKHEencrypt(c2,MKhe.mklwe_sk,m2,&MKhe.parmk);

    MKLweSample c_nand(q,n,k);
    MKLweSample c_nandboot(q,n,k);

    clock_t start = clock();
    c_nand=parLWE.nand_const-(c1+c2);
    MKBootstrap_v2(c_nandboot,c_nand,MKhe.mkBRK,MKhe.mkReKey,MKhe.nrk0,MKhe.nrk1,MKhe.mklksk,&MKhe.parmk); 
    cout<<"NAND_Bootstrap_v2:"<<float(clock()-start)* 1000/CLOCKS_PER_SEC<<"ms"<<endl;
    

    int m_nand=MKhe.MKHEDecrypt(c_nand,MKhe.mklwe_sk,&MKhe.parmk);
    int m_boot=MKhe.MKHEDecrypt(c_nandboot,MKhe.mklwe_sk,&MKhe.parmk);

    //correct test
    if(m_boot == !(m1 & m2)){
        cout<<"correct"<<endl;
    }
    else{
        cout<<"error"<<endl;
    }
}







int main()
{   


    test_nandbootstrap_v1(MKHE2party_v1);

    //PARAM_MKHE4party_v2/  PARAM_MKHE8party_v2/ PARAM_MKHE16party_v2
    test_nandbootstrap_v2(MKHE2party_v2);
    


 



    
    return 0;
    
    
}

#include "MKHEkeygen.h"
#include "MKHEparams.h"
#include "sampler.h"
#include "poly.h"
#include <cmath>
#include <iostream>

//s1,...,sk
void MKLweKeyGen(MKLweKey &result,MKHEParams *parmk){
     Sampler s(parLWE);
     int n=parmk->n;
     int parties=parmk->parties;
     result.resize(n*parties);
     s.get_binary_vector(result);
}

void MKRLweKeyGen(MKRLweKey &result,MKHEParams *parmk){
     int n=parmk->n;
     int parties=parmk->parties;
     int N=parmk->N;
     int N2p1=N/2+1;
     int d=parmk->d_r;
     int Q=parmk->Q;
     ModQPoly tmp(N);
     FFTPoly tmp_FFT(N2p1);
     //generate CRS a in RQ^d
     vector<FFTPoly> a;
     for(int i=0;i<d;i++){
          get_uniform_vector(tmp,Q);
          fftN.to_fft(tmp_FFT,tmp);
          a.push_back(tmp_FFT);
     }
  
     FFTPoly b_FFT(N2p1);
     ModQPoly e(N);
     FFTPoly e_FFT(N2p1);
     FFTPoly z_FFT(N2p1);

     //generate mkrlwe_sk=(z0,...,zk) and mkrlwe_pk=(b0,b1,..bk,a)
     for(int i=0;i<=parties;i++){     
          //z_0=1;
          if(i==0){
               tmp.clear();
               tmp.resize(N,0);
               tmp[0]=1;
               result.MKRLweKey_SK_poly.push_back(tmp);
               fftN.to_fft(z_FFT,tmp);
               result.MKRLweKey_SK.push_back(z_FFT);

          }
          else{ //zi
               get_gaussian_vector(tmp,parmk->stdev_RLWEkey);  
               result.MKRLweKey_SK_poly.push_back(tmp);
               fftN.to_fft(z_FFT,tmp);
               result.MKRLweKey_SK.push_back(z_FFT);
          }

          vector<FFTPoly> b_vec(d);
          for(int j=0;j<d;j++){
               get_gaussian_vector(e,parmk->stdev_RLWEerr);
               fftN.to_fft(e_FFT,e);
               b_FFT=e_FFT-a[j]*z_FFT;
               b_vec[j]=b_FFT;
          }
          result.MKRLweKey_PK.push_back(b_vec);
     }
     result.MKRLweKey_PK.push_back(a);
}

void HybridProduct_rksk(vector<ModQPoly>& c,const HPK& hpk,const MKRLweKey& mkrlwe_key,MKHEParams *parmk){
    int k=parmk->parties;
    assert(c.size()==k+1);  //c0,..,ck
    int N=parmk->N;
    int N2p1=parmk->N2p1;
    int d=parmk->d_r;
    int B=parmk->B_r;
    int shift=parmk->shift_r;
    int Q=parmk->Q;
    int party=hpk.party;
    vector<int> tmp_poly(N);
    vector<long> tmp_poly_long(N);

    //vi= g^-1(ci)*bi  
    vector<ModQPoly> v(k+1);
    for(int i=0;i<=k;i++){
        external_product(tmp_poly_long,c[i],mkrlwe_key.MKRLweKey_PK[i],B,shift,d,N,N2p1);  //g^-1(ci)*bi
        mod_q_boot(tmp_poly,tmp_poly_long);
        v[i]=tmp_poly;
    }

    ModQPoly u(N);
    for(int i=0;i<=k;i++){
        external_product(tmp_poly_long,c[i],hpk.hpk_FFT[0],B,shift,d,N,N2p1);  //g^-1(ci)*d0       
        mod_q_boot(u,tmp_poly_long);
        c[i]=u;
    }

    ModQPoly w0(N,0);
    for(int i=0;i<=k;i++){
        external_product(tmp_poly_long,v[i],hpk.hpk_FFT[2],B,shift,d,N,N2p1);  //g^-1(vi)*d2
        mod_q_boot(tmp_poly,tmp_poly_long);
        for(int j=0;j<N;j++){
            w0[j]+=tmp_poly[j];
        }
    }
    mod_q_boot(w0);

    ModQPoly w1(N,0);
    for(int i=0;i<=k;i++){
        external_product(tmp_poly_long,v[i],hpk.hpk_FFT[1],B,shift,d,N,N2p1);  //g^-1(vi)*d1
        mod_q_boot(tmp_poly,tmp_poly_long);
        for(int j=0;j<N;j++){
            w1[j]+=tmp_poly[j];
        }
    } 
    mod_q_boot(w1);

    for(int i=0;i<N;i++){
        c[0][i]+=w0[i];
        c[party][i]+=w1[i];
    }
    mod_q_boot(c[0]);
    mod_q_boot(c[party]);
}

//Only z_i and a in mkrlwekey is used
void HPKGen(HPK& hpk,const MKRLweKey &mkrlwe_key,const ModQPoly& f,int party,MKHEParams *parmk){
     int parties=parmk->parties;
     int Q=parmk->Q; 
     int N=parmk->N;
     int N2p1=parmk->N2p1;
     int B=parmk->B_r;
     int d=parmk->d_r;

     ModQPoly tmp(N);
     FFTPoly tmp_FFT(N2p1);
     vector<FFTPoly> d0(d);
     vector<FFTPoly> d1(d);
     vector<FFTPoly> d2(d);
     FFTPoly z_FFT=mkrlwe_key.MKRLweKey_SK[party];
     vector<FFTPoly> a=mkrlwe_key.MKRLweKey_PK[parties+1];
     ModQPoly r(N);
     ModQPoly e(N);
     FFTPoly r_FFT(N2p1);
     FFTPoly e_FFT(N2p1);
     FFTPoly f_FFT(N2p1);
     fftN.to_fft(f_FFT,f);
     FFTPoly fPOWB_FFT=f_FFT;
     get_gaussian_vector(r,parmk->stdev_RLWEkey); 
     fftN.to_fft(r_FFT,r);
     FFTPoly rPOWB_FFT=r_FFT;

     for(int i=0;i<d;i++){
          //d0
          get_gaussian_vector(e,parmk->stdev_RLWEerr);
          fftN.to_fft(e_FFT,e);
          tmp_FFT=r_FFT*a[i]+fPOWB_FFT+e_FFT;
          d0[i]=tmp_FFT;
          mult_fft_poly_by_int(fPOWB_FFT,B);

          //d1
          get_uniform_vector(tmp,Q);
          fftN.to_fft(tmp_FFT,tmp);
          d1[i]=tmp_FFT;

          //d2
          fftN.to_fft(f_FFT,f);
          get_gaussian_vector(e,parmk->stdev_RLWEerr);
          fftN.to_fft(e_FFT,e);
          tmp_FFT=rPOWB_FFT+e_FFT-z_FFT*d1[i];
          d2[i]=tmp_FFT;
          mult_fft_poly_by_int(rPOWB_FFT,B);
     }
     hpk.party=party;
     hpk.hpk_FFT={d0,d1,d2};
}

void MKrkskGen(vector<vector<FFTPoly>> &mkrksk,const vector<HPK>& mkhpk,const MKRLweKey &mkrlwe_key,MKHEParams *parmk){
     int Q=parmk->Q;
     int N=parmk->N;
     int k=parmk->parties;
     int B=parmk->B_r;
     int d=parmk->d_r;

     ModQPoly tmp_poly(N,0);
     vector<long> tmp_long(N);
     vector<ModQPoly> mkrlwe_c(k+1);
     vector<FFTPoly> mkrlwe_c_fft(k+1);
     FFTPoly tmp_fft(N/2+1);
     mkrksk=vector<vector<FFTPoly>>(k+1,vector<FFTPoly>(d,FFTPoly(N/2+1)));

     int powB=1;
     for(int i=0;i<d;i++){
          tmp_poly.assign(N,0);
          for(int j=0;j<=k;j++){
               mkrlwe_c[j]=tmp_poly;
          }
          mkrlwe_c[0][0]=powB;

          for(int j=0;j<k;j++){
             HybridProduct_rksk(mkrlwe_c,mkhpk[j],mkrlwe_key,parmk);
          }

          for(int j=0;j<=k;j++){
               fftN.to_fft(tmp_fft,mkrlwe_c[j]);
               mkrksk[j][i]=tmp_fft;
          }
          powB*=B;
     }
}

void MKHPKGen(vector<HPK>& mkhpk,const MKRLweKey &mkrlwe_key,const vector<SKey_boot>& vec_sk_boot,MKHEParams *parmk){
     int parties=parmk->parties;
     HPK hpk;
     for(int i=0;i<parties;i++){
          // party index in 1...parties ,so party=i+1
          HPKGen(hpk,mkrlwe_key,vec_sk_boot[i].sk,i+1,parmk);  //vec_sk_boot[i]=f_i+1
          mkhpk.push_back(hpk);
     }
}




void ReKeyGen(vector<FFTPoly> &ReKey,const SKey_boot& sk_boot,MKHEParams *parmk){
     int Q=parmk->Q;
     int N=parmk->N;
     int B=parmk->B_n;
     int d=parmk->d_n;
     int shift=parmk->shift_n;
     NTRU_vector_enc(ReKey,sk_boot.sk_inv,d,B,sk_boot,parmk);
}

void MKReKeyGen(vector<vector<FFTPoly>> &mkReKey,const vector<SKey_boot>& vec_sk_boot,MKHEParams *parmk){
     int parties=parmk->parties;
     int N=parmk->N;
     int d=parmk->d_n;
     vector<FFTPoly> Rekey(d);
     for(int i=0;i<parties;i++){
          ReKeyGen(Rekey,vec_sk_boot[i],parmk);
          mkReKey.push_back(Rekey);
     }
}



void NRKGen(vector<vector<FFTPoly>> &nrksk0,vector<RGSWKey> &nrksk1,const vector<SKey_boot>& vec_sk_boot,const MKRLweKey &mkrlwe_key,MKHEParams *parmk){
     int Q=parmk->Q;
     int N=parmk->N;
     int B=parmk->B_r;
     int d=parmk->d_r;
     int shift=parmk->shift_r;
     int N2p1=parmk->N2p1;
     int k=parmk->parties;

     ModQPoly tmp_poly(N);
     ModQPoly e(N);
     ModQPoly r(N);
     FFTPoly tmp_FFT(N2p1);
     FFTPoly tmp_FFT2(N2p1);
     FFTPoly e_FFT(N2p1);
     FFTPoly r_FFT(N2p1);
     ModQPoly f(N);
     FFTPoly f_FFT(N2p1);
     FFTPoly fPOWB_FFT(N2p1);

     vector<FFTPoly> b_sum(d,FFTPoly(N2p1));
     for(int i=0;i<d;i++){
          //b1=party1
          b_sum[i]=mkrlwe_key.MKRLweKey_PK[1][i];
          for(int j=2;j<=k;j++){
               b_sum[i]+=mkrlwe_key.MKRLweKey_PK[j][i];
          }
     }
     vector<FFTPoly> a=mkrlwe_key.MKRLweKey_PK[k+1];
     //nrksk0=(br+f1g+e,ar+e) party1
     nrksk0=vector<vector<FFTPoly>>(2,vector<FFTPoly>(d));
     get_gaussian_vector(r,parmk->stdev_RLWEkey);
     fftN.to_fft(r_FFT,r);
     f=vec_sk_boot[0].sk;
     fftN.to_fft(f_FFT,f);
     fPOWB_FFT=f_FFT;
     for(int i=0;i<d;i++){
          get_gaussian_vector(e,parmk->stdev_RLWEerr);
          fftN.to_fft(e_FFT,e);
          tmp_FFT=r_FFT*b_sum[i]+fPOWB_FFT+e_FFT;
          get_gaussian_vector(e,parmk->stdev_RLWEerr);
          fftN.to_fft(e_FFT,e);
          tmp_FFT2=r_FFT*a[i]+e_FFT;
          nrksk0[0][i]=tmp_FFT;
          nrksk0[1][i]=tmp_FFT2;
          mult_fft_poly_by_int(fPOWB_FFT,B);
     }
     // nrksk1 party 2-k
     vector<vector<FFTPoly>> RGSW_row1(2,vector<FFTPoly>(d));
     vector<vector<FFTPoly>> RGSW_row2(2,vector<FFTPoly>(d));
     vector<vector<vector<FFTPoly>>> RGSWkey(2);
     for(int i=1;i<k;i++){
          f=vec_sk_boot[i].sk;
          fftN.to_fft(f_FFT,f);
          //row1
          get_gaussian_vector(r,parmk->stdev_RLWEkey);
          fftN.to_fft(r_FFT,r);
          fPOWB_FFT=f_FFT;
          for(int j=0;j<d;j++){
               get_gaussian_vector(e,parmk->stdev_RLWEerr);
               fftN.to_fft(e_FFT,e);
               tmp_FFT=r_FFT*b_sum[j]+fPOWB_FFT+e_FFT;
               get_gaussian_vector(e,parmk->stdev_RLWEerr);
               fftN.to_fft(e_FFT,e);
               tmp_FFT2=r_FFT*a[j]+e_FFT;
               RGSW_row1[0][j]=tmp_FFT;
               RGSW_row1[1][j]=tmp_FFT2;
               mult_fft_poly_by_int(fPOWB_FFT,B);
          }
          //row2
          get_gaussian_vector(r,parmk->stdev_RLWEkey);
          fftN.to_fft(r_FFT,r);
          fPOWB_FFT=f_FFT;
          for(int j=0;j<d;j++){
               get_gaussian_vector(e,parmk->stdev_RLWEerr);
               fftN.to_fft(e_FFT,e);
               tmp_FFT=r_FFT*b_sum[j]+e_FFT;
               get_gaussian_vector(e,parmk->stdev_RLWEerr);
               fftN.to_fft(e_FFT,e);
               tmp_FFT2=r_FFT*a[j]+fPOWB_FFT+e_FFT;
               RGSW_row2[0][j]=tmp_FFT;
               RGSW_row2[1][j]=tmp_FFT2;
               mult_fft_poly_by_int(fPOWB_FFT,B);
          }
          RGSWkey[0]=RGSW_row1;
          RGSWkey[1]=RGSW_row2;
          nrksk1.push_back(RGSWkey);
     }
}

void BRKGen(vector<NGSFFTctxt>& BRK,const SKey_base_LWE& sk_base, const SKey_boot& sk_boot,MKHEParams *parmk){
     int n=parmk->n;
     int B=parmk->B_n;
     int d=parmk->d_n;
     int N2p1=parmk->N2p1;
     BRK=vector<NGSFFTctxt> (n);
     for(int i=0;i<n;i++){
          NGSFFTctxt brk(d);
          NTRU_vector_enc(brk, sk_base[i], d, B, sk_boot,parmk);
          BRK[i]=brk;
     }
}

void MKBRKGen(MKBRKey& mkBRK,const MKLweKey& mklwe_sk,const vector<SKey_boot>& vec_sk_boot,MKHEParams *parmk){
     int n=parmk->n;
     int N=parmk->N;
     int parties=parmk->parties;
     vector<NGSFFTctxt> BRK;
     SKey_boot sk_boot;
     for(int i=0;i<parties;i++){
          SKey_base_LWE sk_base(mklwe_sk.begin() + i*n, mklwe_sk.begin() + (i+1)*n);
          sk_boot=vec_sk_boot[i];
          BRKGen(BRK,sk_base,sk_boot,parmk);
          mkBRK.push_back(BRK);
     }
}

void MKBOOTSKGen(vector<SKey_boot> &vec_sk_boot,MKHEParams *parmk){
     int N=parmk->N;
     int parties=parmk->parties;
     SKey_boot sk_boot;
     for(int i=0;i<parties;i++){
          sk_boot.sk = ModQPoly(N,0);
          sk_boot.sk_inv = ModQPoly(N,0);
          get_invertible_vector(sk_boot.sk, sk_boot.sk_inv,1,0);
          vec_sk_boot.push_back(sk_boot);
     }
}

void LKSKGen(vector<KSKey_LWE> &lksk,const vector<int>& sk_in,const vector<int>& sk_out,MKHEParams *parmk){
     int B=parmk->B_l;
     int d=parmk->d_l;
     int q=parmk->q;
     int q_half=q/2;
     int n=sk_out.size();
     int N=sk_in.size();
     assert(lksk.size()==N);
     //generate Ai
     for(int i=0;i<N;i++){
          lksk[i].A.clear();
          vector<int> row(n,0);
          for(int j=0;j<d;j++){
          get_uniform_vector(row,q);
          lksk[i].A.push_back(row);
          }
     }

     vector<int> e(d,0);
     //generate bi
     for(int i=0;i<N;i++){
          lksk[i].b=vector<int>(d,0);
          get_gaussian_vector(e,parmk->stdev_LWEerr);
          int powB=1;
          for(int j=0;j<d;j++){
             //-As
             long tmp=0;
             for(int k=0;k<n;k++){
               tmp-=long(lksk[i].A[j][k])*long(sk_out[k]);
             }
             tmp+=e[j]+long(sk_in[i])*long(powB);
             lksk[i].b[j]=mod_q_lwe(tmp,q,q_half);
             powB*=B;
          }
     }
}

void MKLKSKGen(MKlksk &mklksk,const MKLweKey& sk_in,const MKLweKey& sk_out,MKHEParams *parmk){
     int parties=parmk->parties;
     int N=parmk->N;
     int n=parmk->n;
     for(int i=0;i<parties;i++){
          vector<KSKey_LWE> lksk(N);
          vector<int> z(sk_in.begin() + i*N, sk_in.begin() + (i+1)*N);
          vector<int> s(sk_out.begin() + i*n, sk_out.begin() + (i+1)*n);
          LKSKGen(lksk,z,s,parmk);
          mklksk.push_back(lksk);
     }
}








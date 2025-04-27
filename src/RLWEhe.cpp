#include "RLWEhe.h"




void RLWE_Enc(vector<ModQPoly> &rlwe_c,const ModQPoly &m,const ModQPoly &z){
    int N=MKHEParams::N_static;
    int Q=MKHEParams::Q_static;
    int N2p1=MKHEParams::N2p1_static;
    assert(rlwe_c.size()==2);
    ModQPoly tmp_poly(N);
    vector<long> tmp_long(N);
    FFTPoly tmp_fft(N2p1);

    ModQPoly a(N);
    ModQPoly e(N);
    FFTPoly e_fft(N2p1);
    FFTPoly a_fft(N2p1);
    FFTPoly z_fft(N2p1);
    FFTPoly m_fft(N2p1);
    fftN.to_fft(m_fft,m);
    fftN.to_fft(z_fft,z);
    get_uniform_vector(a,Q);
    fftN.to_fft(a_fft,a);
    get_gaussian_vector(e,0.01);
    fftN.to_fft(e_fft,e);
    tmp_fft=e_fft-a_fft*z_fft+m_fft;
    fftN.from_fft(tmp_long,tmp_fft);
    mod_q_poly(tmp_poly,tmp_long,Q);

    rlwe_c[0]=tmp_poly;
    rlwe_c[1]=a;
    
}


void RLWE_Dec(ModQPoly &res,const vector<ModQPoly> &rlwe_c,const ModQPoly &z){
    int N=MKHEParams::N_static;
    int Q=MKHEParams::Q_static;
    int N2p1=MKHEParams::N2p1_static;
    assert(rlwe_c.size()==2);

    ModQPoly tmp_poly(N);
    vector<long> tmp_long(N);
    FFTPoly tmp_fft(N2p1);

    ModQPoly c0=rlwe_c[0];
    ModQPoly c1=rlwe_c[1];
    FFTPoly c_fft(N2p1);
    FFTPoly z_fft(N2p1);
    fftN.to_fft(c_fft,c1);
    fftN.to_fft(z_fft,z);

    tmp_fft=c_fft*z_fft;
    fftN.from_fft(tmp_long,tmp_fft);
    mod_q_poly(tmp_poly,tmp_long,Q);

    for(int i=0;i<N;i++){
        res[i]=c0[i]+tmp_poly[i];
    }
    mod_q_poly(res,Q);
}


void RLWE_RGSWproduct(vector<ModQPoly>& res,const vector<ModQPoly> &rlwe_c,const RGSWKey &RGSW_c,MKHEParams *parmk){
    int N=parmk->N;
    int k=parmk->parties;
    int Q=parmk->Q;
    int N2p1=parmk->N2p1;
    int B=parmk->B_r;
    int shift=parmk->shift_r;
    int d=parmk->d_r;

    ModQPoly tmp_poly(N);
    ModQPoly tmp_poly2(N);
    vector<long> tmp_long(N);
    FFTPoly tmp_fft(N2p1);

    res=vector<ModQPoly>(2,ModQPoly(N));


    external_product(tmp_long,rlwe_c[0],RGSW_c[0][0],B,shift,d,N,N2p1);
    mod_q_boot(tmp_poly,tmp_long);

    external_product(tmp_long,rlwe_c[1],RGSW_c[1][0],B,shift,d,N,N2p1);
    mod_q_boot(tmp_poly2,tmp_long);

    for(int i=0;i<N;i++){
        res[0][i]=tmp_poly[i]+tmp_poly2[i];
    }


    external_product(tmp_long,rlwe_c[0],RGSW_c[0][1],B,shift,d,N,N2p1);
    mod_q_boot(tmp_poly,tmp_long);

    external_product(tmp_long,rlwe_c[1],RGSW_c[1][1],B,shift,d,N,N2p1);
    mod_q_boot(tmp_poly2,tmp_long);

    for(int i=0;i<N;i++){
        res[1][i]=tmp_poly[i]+tmp_poly2[i];
    }

    mod_q_boot(res[0]);
    mod_q_boot(res[1]);

}


#include "NTRUhe.h"

void NTRU_scarlar_enc(ModQPoly &res,const ModQPoly& m,const SKey_boot& sk,int Q,int N){
    ModQPoly g(N);
    get_gaussian_vector(g,0.25);
    for(int i=0;i<N;i++){
        res[i]=m[i]+g[i];//m+g
    }
    FFTPoly res_FFT(N/2+1);
    fftN.to_fft(res_FFT,res);
    ModQPoly f_inv=sk.sk_inv;
    FFTPoly f_inv_FFT(N/2+1);
    fftN.to_fft(f_inv_FFT,f_inv);
    res_FFT*=f_inv_FFT; // m/f+ g/f
    vector<long> tmp_poly_long(N);
    fftN.from_fft(tmp_poly_long,res_FFT);
    mod_q_poly(res,tmp_poly_long,Q);
}

void NTRU_scarlar_dec(ModQPoly& m_dec,const ModQPoly& c,const SKey_boot& sk,int Q,int N){
    Sampler s(parLWE);
    m_dec.clear();
    m_dec.resize(N,0);
    FFTPoly c_FFT(N/2+1);
    FFTPoly f_FFT(N/2+1);
    FFTPoly tmp_FFT(N/2+1);
    vector<long> tmp_poly_long(N);
    ModQPoly f=sk.sk;
    fftN.to_fft(f_FFT,f);
    fftN.to_fft(c_FFT,c);
    tmp_FFT=f_FFT*c_FFT;  //c*f
    fftN.from_fft(tmp_poly_long,tmp_FFT);
    mod_q_poly(m_dec,tmp_poly_long,Q);
}

void NTRU_vector_enc(NGSFFTctxt& ct, int m, int l, int B, const SKey_boot& sk_boot,MKHEParams *parmk)
{
    ModQPoly msg(parmk->N,0L);
    msg[0] = m; // msg = m (degree-0 polynomial)
    NTRU_vector_enc(ct, msg, l, B, sk_boot,parmk);
}


void NTRU_vector_enc(NGSFFTctxt& ct, const ModQPoly& m, int l, int B, const SKey_boot& sk_boot,MKHEParams *parmk)
{
    if(ct.size() != l)
        ct = NGSFFTctxt(l);
    int N=parmk->N;
    int N2p1=parmk->N2p1;
    double stdev_NTRUerr=0.09;
    FFTPoly sk_boot_inv_fft(N2p1); // f^-1 in FFT form
    fftN.to_fft(sk_boot_inv_fft, sk_boot.sk_inv);
    FFTPoly g_fft(N2p1);
    ModQPoly msg(m); // at each iteration i, msg will be equal to m * B^i
    FFTPoly msg_powB(N2p1);
    fftN.to_fft(msg_powB, msg); // FFT of m * B^i
    FFTPoly tmp_ct(N2p1);
    vector<long> tmp_ct_long(N);
    vector<int> tmp_ct_int(N);
    for (int i = 0; i < l; i++)
    {
        // sample random ternary vector
        ModQPoly g(N,0L);
        get_gaussian_vector(g,parmk->stdev_NTRUerr); 
        //  get_gaussian_vector(g,0.09); 
        // FFT transform it
        fftN.to_fft(g_fft, g);
        // compute g * sk_boot^(-1)
        tmp_ct = g_fft * sk_boot_inv_fft;
        // compute g * sk_boot^(-1) + B^i * m
        tmp_ct += msg_powB;
        // inverse FFT of the above result
        fftN.from_fft(tmp_ct_long, tmp_ct);
        // reduction modulo q_boot
        mod_q_boot(tmp_ct_int, tmp_ct_long);
        // FFT transform for further use
        fftN.to_fft(tmp_ct, tmp_ct_int);
        ct[i] = tmp_ct;
        mult_fft_poly_by_int(msg_powB, B); // msg_powB = msg * B^i
    }
}


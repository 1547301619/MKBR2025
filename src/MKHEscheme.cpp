#include "MKHEscheme.h"
using namespace std;


MKHEscheme::MKHEscheme(BINFHE_PARAMSET set,int v){

    static const std::map<BINFHE_PARAMSET, MKHEParams> paramMap = {
        {MKHE2party_v1, PARAM_MKHE2party_v1},
        {MKHE2party_v2, PARAM_MKHE2party_v2},
        {MKHE4party_v2, PARAM_MKHE4party_v2},
        {MKHE8party_v2, PARAM_MKHE8party_v2},
        {MKHE16party_v2, PARAM_MKHE16party_v2}
    };

   
    auto it = paramMap.find(set);
    if (it != paramMap.end()) {
        this->parmk = it->second; 
    } else {
        throw std::invalid_argument("Invalid BINFHE_PARAMSET value");
    }

    cout<<"----------------------KEYGEN---------------------"<<endl;

    MKLweKeyGen(mklwe_sk,&parmk);
    MKBOOTSKGen(vec_sk_boot,&parmk);
    MKRLweKeyGen(mkrlwe_key,&parmk);
    MKHPKGen(mkhpk,mkrlwe_key,vec_sk_boot,&parmk);
    MKExtractKey(mklwe_sk_z,mkrlwe_key,&parmk);
    MKLKSKGen(mklksk,mklwe_sk_z,mklwe_sk,&parmk);
    MKReKeyGen(mkReKey,vec_sk_boot,&parmk);
    MKBRKGen(mkBRK,mklwe_sk,vec_sk_boot,&parmk);

    if(v==1){
       MKrkskGen(mkrksk,mkhpk,mkrlwe_key,&parmk);
    }
    if(v==2){
       NRKGen(nrk0,nrk1,vec_sk_boot,mkrlwe_key,&parmk);
    }

    cout<<"------------------KEYGEN SUCCESS------------------"<<endl;
}


void MKHEscheme::MKHEencrypt(MKLweSample &result,const MKLweKey& key,int m,MKHEParams *parmk){
    int n=result.n;
    int parties=result.parties;

    Sampler s(parLWE);
    normal_distribution<double> gaussian_sampler(0.0, parmk->stdev_LWEerr);
    int b=parLWE.delta_base*m+static_cast<int>(round(gaussian_sampler(rand_engine)));
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result.a[i*n +j] =s.mod_q_base_sampler(rand_engine);
            b -= result.a[i*n +j]*key[i*n+j];
        } 
    }
    result.b=parLWE.mod_q_base(b);
}

int MKHEscheme::MKHEDecrypt(const MKLweSample &ctxt,const MKLweKey& key,MKHEParams *parmk){
    int n=ctxt.n;
    int parties=ctxt.parties;
    int q=parmk->q;
    int t=parmk->t;

    long output_long=ctxt.b;
        for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            output_long += ctxt.a[i*n +j]*key[i*n+j];
        } 
    }

    int output = parLWE.mod_q_base(output_long);
    output = int(round(double(output*t)/double(q)));
    return output;
}



void MKBootstrap_v1(MKLweSample &result,const MKLweSample &ct,const MKBRKey& mkBRK,const vector<vector<FFTPoly>> &mkReKey,const vector<vector<FFTPoly>> &mkrksk,const MKlksk& mklksk,MKHEParams *parmk){
    int N=parmk->N;
    int k=parmk->parties;
    int q=parmk->q;
    int n=parmk->n;
    int Q=parmk->Q;
    int N2=N*2;
    int N2p1=parmk->N2p1;
    int B=parmk->B_r;
    int shift=parmk->shift_r;
    int d=parmk->d_r;

    
    int b=int((ct.b*N2)/q);
    int *a=ct.a;

    ModQPoly tmp_poly(N,0);
    vector<long> tmp_long(N,0L);
    vector<long> acc_long(N);
    vector<int> acc(N, Q/8);
    int b_pow = (N/2 + b)%N2;
    if (b_pow < 0)
        b_pow += N2;
    int b_sign = 1;
    if (b_pow >= N)
    {
        b_pow -= N;
        b_sign = -1;
    }
    for (int i = 0; i < b_pow; ++i)
        acc[i] = (b_sign == 1) ? -acc[i]: acc[i];
    for (int i = b_pow; i < N; ++i)
        acc[i] = (b_sign == 1) ? acc[i]: -acc[i];

    ModQPoly c_ntru=acc;
    for(int i=0;i<k;i++){
        SingleKeyReBR(acc,a + n * i,c_ntru,mkReKey[i],mkBRK[i],parmk);
        c_ntru=acc;
    }

 
    vector<ModQPoly> mkrlwe_c(k+1);

    //g^(c) * RKSK
    for(int i=0;i<=k;i++){
        external_product(tmp_long,acc,mkrksk[i],B,shift,d,N,N2p1);
        mod_q_poly(tmp_poly,tmp_long,Q);
        mkrlwe_c[i]=tmp_poly;
    }

    mkrlwe_c[0][0]+=Q/8; 

    MKLweSample c_z(Q,N,k);
    MKExtract(c_z,mkrlwe_c,parmk);   //in ZQ
    
    MKModSwitch(c_z,Q,q,parmk);  //in zq

    MKkeyswitch(result,c_z,mklksk,parmk);

    a = nullptr; 
}



void MKBootstrap_v2(MKLweSample &result,const MKLweSample &ct,const MKBRKey &mkBRK,const vector<vector<FFTPoly>> &mkreKey,const vector<vector<FFTPoly>> &nrk0,const vector<RGSWKey> &nrk1,const MKlksk& mklksk,MKHEParams *parmk){
     
    int N=parmk->N;
    int k=parmk->parties;
    int q=parmk->q;
    int n=parmk->n;
    int Q=parmk->Q;
    int N2=N*2;
    int N2p1=parmk->N2p1;
    int B=parmk->B_r;
    int shift=parmk->shift_r;
    int d=parmk->d_r;
    
    int b=int((ct.b*N2)/q);
    int *a=ct.a;
 
    ModQPoly tmp_poly(N,0);
    vector<long> tmp_long(N,0L);
    vector<long> acc_long(N);
    vector<int> acc(N, Q/8);
    int b_pow = (N/2 + b)%N2;
    if (b_pow < 0)
        b_pow += N2;
    int b_sign = 1;
    if (b_pow >= N)
    {
        b_pow -= N;
        b_sign = -1;
    }
    for (int i = 0; i < b_pow; ++i)
        acc[i] = (b_sign == 1) ? -acc[i]: acc[i];
    for (int i = b_pow; i < N; ++i)
        acc[i] = (b_sign == 1) ? acc[i]: -acc[i];


    ModQPoly c_ntru=acc;
    SingleKeyReBR(acc,a,c_ntru,mkreKey[0],mkBRK[0],parmk);


    vector<ModQPoly> rlwe_c(2);
    for(int i=0;i<2;i++){
        external_product(tmp_long,acc,nrk0[i],B,shift,d,N,N2p1);
        mod_q_poly(tmp_poly,tmp_long,Q);
        rlwe_c[i]=tmp_poly;
    }
   
 
    vector<ModQPoly> rlwe_res(2);
    for(int i=1;i<k;i++){
        SingleKeyReBR(acc,a+n*i,rlwe_c[0],mkreKey[i],mkBRK[i],parmk);
        rlwe_c[0]=acc;

        SingleKeyReBR(acc,a+n*i,rlwe_c[1],mkreKey[i],mkBRK[i],parmk);
        rlwe_c[1]=acc;

        //nrk1[i-1]=partyi Vi
        RLWE_RGSWproduct(rlwe_res,rlwe_c,nrk1[i-1],parmk);
        rlwe_c=rlwe_res;
    }



    rlwe_c[0][0]+=Q/8; 

    MKLweSample c_z(Q,N,k);
    RlweExtract(c_z,rlwe_c,parmk);   //in ZQ
    LWEModSwitch(c_z,Q,q,parmk);  //in zq
    MKkeyswitch_v2(result,c_z,mklksk,parmk);
    a = nullptr; 
}



void SingleKeyReBR(ModQPoly &acc, int *al,const ModQPoly& c,const NGSFFTctxt& ReKey,const vector<NGSFFTctxt>& BRK,MKHEParams *parmk){
    int N = parmk->N;
    int N2=N*2;
    int N2p1 = parmk->N2p1;
    int q=parmk->q;
    int n=parmk->n;
    int B=parmk->B_n;
    int shift=parmk->shift_n;
    int d=parmk->d_n;


    vector<int> tmp_poly(N);
    vector<long> tmp_poly_long(N);

    acc.resize(N,0L);


    //acc=c * brk_n , Rekey=brk_n
    external_product(tmp_poly_long, c, ReKey, B, shift, d,N,N2p1);   
    mod_q_boot(acc, tmp_poly_long);


    vector<int> a(al,al+n);
    // switch to modulus 2*N
    modulo_switch_a(a, q, N2,n);


    int coef, coef_sign;

    for(int i=0;i<n;i++){
            coef = a[i];
            if (coef == 0) continue;
            coef_sign = 1;
            if (coef < 0) coef += N2;
            if (coef >= N)
            {
                coef -= N;
                coef_sign = -1;
            }

            // acc * (X^coef - 1)
            if (coef_sign == 1)//coef>0,right rotation
            {
                for (int i = 0; i<coef; ++i)
                    tmp_poly[i] = mod_q_boot(-acc[i-coef+N] - acc[i]);
                for (int i = coef; i < N; ++i)
                    tmp_poly[i] = mod_q_boot(acc[i-coef] - acc[i]);
            }
            else 
            {   
                for (int i = 0; i<coef; ++i)
                    tmp_poly[i] = mod_q_boot(acc[i-coef+N] - acc[i]);
                for (int i = coef; i < N; ++i)
                    tmp_poly[i] = mod_q_boot(-acc[i-coef] - acc[i]);
            }

            // acc * (X^coef - 1) x brk[i]
            external_product(tmp_poly_long, tmp_poly, BRK[i], B, shift, d,N,N2p1);
            mod_q_boot(tmp_poly, tmp_poly_long);
            // acc * (X^coef - 1) x brk[i] + acc
            for (int i = 0; i<N; ++i)
                acc[i] += tmp_poly[i];

            mod_q_boot(acc); //if Q=2^27
        }

}




void MKkeyswitch(MKLweSample &result,const MKLweSample &input,const MKlksk& mklksk,MKHEParams *parmk){
    int n=parmk->n;
    assert(result.n==n);
    int q=parmk->q;
    int q_half=q/2;
    int Q=q;
    int Q_half=q_half;
    int B=parmk->B_l;
    int d=parmk->d_l;
    int shift=parmk->shift_l;
    int parties=parmk->parties;
    int N=parmk->N;
    int bound=B>>1;
    int digit;
    int tmp;
    int sgn;
    result.b=input.b;
    int *a_in=input.a;
    vector<long> a_long(n,0L);
    vector<int> a_bar(n,0);
    vector<int> poly_decomp(d);

    for(int k=0;k<parties;k++){
        long b_long=0L;
        int b_bar=0;
        fill(a_long.begin(), a_long.end(), 0L);
        fill(a_bar.begin(), a_bar.end(), 0);
      
        for(int index=0;index<N;index++){
            //g^-1(aij)
            tmp=abs(a_in[index+N*k]);
            sgn=(a_in[index+N*k] < 0)? -1 : 1;
            for (int j = 0; j < d; ++j){
                digit = tmp % B;
                if (digit > bound)
                {
                poly_decomp[j] = (sgn == 1) ? (digit - B): (B - digit);
                tmp >>= shift;
                ++tmp;
                }
                else
                {
                poly_decomp[j] = (sgn == 1) ? digit:  - digit;
                tmp >>= shift;
                } 
            }

            // g^-1(aij)*b
            for(int j=0;j<d;j++){
            b_long+=long(poly_decomp[j])*mklksk[k][index].b[j];
            }
            //g^-1(aij)*A
            for(int j=0;j<d;j++){
            const vector<int>& ksk_row = mklksk[k][index].A[j];
                for(int t=0;t<n;t++){
                 a_long[t] += long(ksk_row[t]) * long(poly_decomp[j]);
                }
            }
        }
        mod_q_poly(a_bar,a_long,Q);
        b_bar=mod_q_lwe(b_long,Q,Q_half);
        result.b+=b_bar;
        for(int i=0;i<n;i++){
            result.a[i+n*k]=a_bar[i];
        }
    }

    a_in=nullptr;

}

//used in MKboot_v9
void MKkeyswitch_v2(MKLweSample &result,const MKLweSample &input,const MKlksk& mklksk,MKHEParams *parmk){
    int n=parmk->n;
    assert(result.n==n);
    int q=parmk->q;
    int q_half=q/2;
    int Q=q;
    int Q_half=q_half;
    int B=parmk->B_l;
    int d=parmk->d_l;
    int shift=parmk->shift_l;
    int parties=parmk->parties;
    int N=parmk->N;
    int bound=B>>1;
    int digit;
    int tmp;
    int sgn;
    result.b=input.b;
    int *a_in=input.a;
    vector<long> a_long(n,0L);
    vector<int> a_bar(n,0);
    vector<int> poly_decomp(d);

    for(int k=0;k<parties;k++){
        long b_long=0L;
        int b_bar=0;
        fill(a_long.begin(), a_long.end(), 0L);
        fill(a_bar.begin(), a_bar.end(), 0);
      
        for(int index=0;index<N;index++){
            //g^-1(aj)
            tmp=abs(a_in[index]);
            sgn=(a_in[index] < 0)? -1 : 1;
            for (int j = 0; j < d; ++j){
                digit = tmp % B;
                if (digit > bound)
                {
                poly_decomp[j] = (sgn == 1) ? (digit - B): (B - digit);
                tmp >>= shift;
                ++tmp;
                }
                else
                {
                poly_decomp[j] = (sgn == 1) ? digit:  - digit;
                tmp >>= shift;
                } 
            }

            // g^-1(aj)*b
            for(int j=0;j<d;j++){
            b_long+=long(poly_decomp[j])*mklksk[k][index].b[j];
            }
            //g^-1(aj)*A
            for(int j=0;j<d;j++){
            const vector<int>& ksk_row = mklksk[k][index].A[j];
                for(int t=0;t<n;t++){
                 a_long[t] += long(ksk_row[t]) * long(poly_decomp[j]);
                }
            }
        }
        mod_q_poly(a_bar,a_long,Q);
        b_bar=mod_q_lwe(b_long,Q,Q_half);
        result.b+=b_bar;
        for(int i=0;i<n;i++){
            result.a[i+n*k]=a_bar[i];
        }
    }

    a_in=nullptr;

}

void MKExtract(MKLweSample& result,const vector<ModQPoly>& mkrlwe_c,MKHEParams *parmk){
    int N=parmk->N;
    assert(result.n==N);
    int parties=parmk->parties;
    result.b=mkrlwe_c[0][0];

    for(int i=1;i<=parties;i++){
        for(int j=0;j<N;j++){
            result.a[(i-1)*N+j]=mkrlwe_c[i][j];
        }
    }
}



void RlweExtract(MKLweSample& result,const vector<ModQPoly>& rlwe_c,MKHEParams *parmk){
    int N=parmk->N;
    assert(result.n==N);
    result.b=rlwe_c[0][0];

    for(int j=0;j<N;j++){
        result.a[j]=rlwe_c[1][j];
    }

}



void MKExtractKey(MKLweKey& mklwe_sk_z,const MKRLweKey& mkrlwe_key,MKHEParams *parmk){
    int N=parmk->N;
    int k=parmk->parties;
    mklwe_sk_z.clear();
    mklwe_sk_z.resize(N*k,0);
    vector<ModQPoly> Z_vec=mkrlwe_key.MKRLweKey_SK_poly;
    //z1,..,zk without z0=1
    for(int i=1;i<=k;i++){
        mklwe_sk_z[(i-1)*N+0]=Z_vec[i][0];
        for(int j=1;j<N;j++){
            mklwe_sk_z[(i-1)*N+j]=-Z_vec[i][N-j];
        }
    }
}

void MKModSwitch(MKLweSample &result,MKHEParams *parmk){
    int Q=parmk->Q;
    int q=parmk->q;
    int N=parmk->N;
    int k=parmk->parties;
    // result.b*=q/Q;
    double radio = double(q)/double(Q);
    result.b=int(round(double(result.b)*radio));
    for (int i=0;i<N*k;i++)
        result.a[i] = int(round(double(result.a[i])*radio));
}

void MKModSwitch(MKLweSample &result,int oldQ,int newQ,MKHEParams *parmk){
    int N=parmk->N;
    int k=parmk->parties;
    double radio = double(newQ)/double(oldQ);
    result.b=int(round(double(result.b)*radio));
    for (int i=0;i<N*k;i++)
        result.a[i] = int(round(double(result.a[i])*radio));
}


void LWEModSwitch(MKLweSample &result,int oldQ,int newQ,MKHEParams *parmk){
    int N=parmk->N;
    int k=parmk->parties;
    result.q=newQ;
    double radio = double(newQ)/double(oldQ);
    result.b=int(round(double(result.b)*radio));
    for (int i=0;i<N;i++)
        result.a[i] = int(round(double(result.a[i])*radio));
}


void MKRlweDec(ModQPoly &res,const vector<ModQPoly>& c,const vector<FFTPoly>& z,MKHEParams *parmk,int Q){
    int N=parmk->N;
    int k=parmk->parties;
    FFTPoly c_FFT(N/2+1);
    res.clear();
    res.resize(N,0);
    FFTPoly res_FFT(N/2+1);
    fftN.to_fft(res_FFT,res);
    vector<long> tmp_poly_long(N);
    for(int i=0;i<=k;i++){
        fftN.to_fft(c_FFT,c[i]);
        res_FFT=res_FFT+c_FFT*z[i];
    }
    fftN.from_fft(tmp_poly_long,res_FFT);
    // mod_q_boot(res,tmp_poly_long);
    mod_q_poly(res,tmp_poly_long,Q);
}

#include "poly.h"
#include <algorithm>
using namespace std;



void get_uniform_vector(vector<int>& vec,int Q)
{   long q=Q;
    uniform_int_distribution<int> mod_q_base_sampler=uniform_int_distribution<int>(-q/2, q/2);
    for(int i=0; i<vec.size(); i++)
        vec[i] = mod_q_base_sampler(rand_engine);
}

void get_gaussian_vector(vector<int>& vec, double st_dev)
{
    normal_distribution<double> gaussian_sampler(0.0, st_dev);
    for(size_t i=0; i<vec.size(); i++)
        vec[i] = static_cast<int>(round(gaussian_sampler(rand_engine)));
}

void get_hwt_vector(vector<int>& vec, int h) {
    int size = vec.size();
    if (h > size) {
        std::cerr << "Error: Hamming weight cannot be greater than vector size." << std::endl;
        return;
    }
    
    std::vector<int> indices(size);
    for (int i = 0; i < size; ++i) indices[i] = i;
    
    random_device rd;
    mt19937 g(rd());
    shuffle(indices.begin(), indices.end(), g);
    
    std::vector<int> signs(h, 1);
    for (int i = h / 2; i < h; ++i) signs[i] = -1;
    shuffle(signs.begin(), signs.end(), g);
    
    for (int i = 0; i < h; ++i) {
        vec[indices[i]] = signs[i];
    }
}


void get_gadget(vector<FFTPoly> &g,int N,int B,int d){
    FFT_engine fftN(N);
    ModQPoly g_poly(N,0);
    g_poly[0]=1;
    FFTPoly g_FFT(N/2+1);
    fftN.to_fft(g_FFT,g_poly);
    for(int i=0;i<d;i++){
        g[i]=g_FFT;
        mult_fft_poly_by_int(g_FFT,B);
    }
}




void polymul(vector<int>& res,const vector<int>& a,const vector<int>& b,long q){
    assert(a.size()==b.size());
    int N=a.size();
    int N2p1=N/2+1;
    FFTPoly a_fft(N2p1);
    fftN.to_fft(a_fft,a);
    FFTPoly b_fft(N2p1);
    fftN.to_fft(b_fft,b);
    FFTPoly c_ftt=a_fft*b_fft;
    vector<long> c(N);
    fftN.from_fft(c,c_ftt);
    printpoly(c);
    lazy_mod_q(res,c,q,q/2);
}

// void polyvector_product();


void gadget_decomp(vector<FFTPoly>& res_vector,const vector<int>& poly,int b, int shift, int l){
    int N = poly.size();
    int N2p1 = N/2+1;
    FFT_engine fftN(N);

    ModQPoly poly_sign(N);
    ModQPoly poly_abs(N);
    vector<int> poly_decomp(N);
    
    FFTPoly tmp_fft(N2p1);

        for (int i = 0; i < N; ++i)
    {
        const int& polyi = poly[i];
        poly_abs[i] = abs(polyi);
        poly_sign[i] = (polyi < 0)? -1 : 1;
    }
    int mask = b-1;
    int bound = b >> 1;
    int digit, sgn;

    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            int& abs_val = poly_abs[i];
            //digit =low b bit of abs_poly[i]
            digit = abs_val & mask;
            if (digit > bound)
            {
                poly_decomp[i] = (poly_sign[i] == 1) ? (digit - b): (b - digit);
                abs_val >>= shift;
                ++abs_val;
            }
            else
            {
                poly_decomp[i] = (poly_sign[i] == 1) ? digit: -digit;
                abs_val >>= shift;
            }
        }
        fftN.to_fft(tmp_fft, poly_decomp);
        res_vector.push_back(tmp_fft);
    }
}



void mod_q_poly(ModQPoly& v,int Q){
    int Q_half=Q/2;
    for(int i=0;i<v.size();i++){
    int coef = v[i]%Q;
    if (coef > Q_half)
        coef = coef - Q;
    if (coef < -Q_half)
        coef = coef + Q;
    v[i]=coef;
    }
}



void mod_q_poly(ModQPoly& v){
    for(int i=0;i<v.size();i++){
    int coef = v[i]%Q;
    if (coef > Q_half)
        coef = coef - Q;
    if (coef < -Q_half)
        coef = coef + Q;
    v[i]=coef;
    }
}

int mod_q_poly(const int input){
    int coef = input%Q;
    if (coef > Q_half)
        coef = coef - Q;
    if (coef < -Q_half)
        coef = coef + Q;
    return coef;
}

int mod_q_poly(const int input,int Q){
    int Q_half=Q/2;
    int coef = input%Q;
    if (coef > Q_half)
        coef = coef - Q;
    if (coef < -Q_half)
        coef = coef + Q;
    return coef;
}




void mod_q_poly(ModQPoly& v,vector<long> tmp_long,int Q){
    long Q_long=long(Q);
    long Q_long_half=long(Q/2);
    for(int i=0;i<v.size();i++){
        int coef = tmp_long[i]%Q_long;
        if (coef > Q_long_half)
            coef = static_cast<int>(coef - Q_long);
        if (coef < -Q_long_half)
            coef = static_cast<int>(coef + Q_long);
        v[i]=coef;
    }
}





void mod_q_poly(ModQPoly& v,vector<long> tmp_long){
    for(int i=0;i<v.size();i++){
        int coef = tmp_long[i]%Q_long;
        if (coef > Q_long_half)
            coef = static_cast<int>(coef - Q_long);
        if (coef < -Q_long_half)
            coef = static_cast<int>(coef + Q_long);
        v[i]=coef;
    }
}


int mod_q_lwe(long input,int q,int q_half){
    int coef = input%q;
    if (coef > q_half)
        return coef - q;
    if (coef < -q_half)
        return coef + q;
    return coef;
}


int mod_q_lwe(int input,int q,int q_half){
    int coef = input%q;
    if (coef > q_half)
        return coef - q;
    if (coef < -q_half)
        return coef + q;
    return coef;
}


void modulo_switch_poly(ModQPoly& v,int old_q,int new_q){
    double ratio = double(new_q)/double(old_q);
    for (int i=0;i<v.size();i++)
        v[i] = int(round(double(v[i])*ratio));
}





void get_invertible_vector(vector<int>& vec, vector<int>& vec_inv, int scale, int shift)
{
    //polynomial with the coefficient vector vec (will be generated later)
    ZZ_pX poly;
    //element of Z_(q_boot)
    ZZ_p coef;
    coef.init(ZZ(q_boot));
    //the inverse of poly modulo poly_mod (will be generated later)
    ZZ_pX inv_poly;
    //random sampling
    while (true)
    {
        //create the polynomial with the coefficient vector of the desired form
        SetCoeff(poly, 0, ternary_sampler(rand_engine)*scale + shift);
        for (size_t i = 1; i < vec.size(); i++)
        {
            coef = ternary_sampler(rand_engine)*scale;
            SetCoeff(poly, i, coef);
        }
        //test invertibility
        try
        {
            InvMod(inv_poly, poly, Param::get_def_poly());
            break;
        }
        catch(...)
        {
            cout << "Polynomial " << poly << " isn't a unit" << endl;
            continue;
        }
    }

    //extract the coefficient vector of poly
    int tmp_coef;
    for (int i = 0; i <= deg(poly); i++)
    {
        tmp_coef = conv<long>(poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec[i] = tmp_coef;
    }

    for (int i = 0; i <= deg(inv_poly); i++)
    {
        tmp_coef = conv<long>(inv_poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec_inv[i] = tmp_coef;
    }
}

void get_invertible_vector_gaussic(vector<int>& vec, vector<int>& vec_inv,double st_dev)
{
    //polynomial with the coefficient vector vec (will be generated later)
    ZZ_pX poly;
    //element of Z_(q_boot)
    ZZ_p coef;
    coef.init(ZZ(q_boot));
    //the inverse of poly modulo poly_mod (will be generated later)
    ZZ_pX inv_poly;
    //random sampling

    normal_distribution<double> gaussian_sampler(0.0, st_dev);

    while (true)
    {
        //create the polynomial with the coefficient vector of the desired form
        // SetCoeff(poly, 0, int(gaussian_sampler(rand_engine))*scale + shift);
        for (size_t i = 0; i < vec.size(); i++)
        {
            coef = static_cast<int>(round(gaussian_sampler(rand_engine)));
            SetCoeff(poly, i, coef);
        }
        //test invertibility
        try
        {
            InvMod(inv_poly, poly, Param::get_def_poly());
            break;
        }
        catch(...)
        {
            cout << "Polynomial " << poly << " isn't a unit" << endl;
            continue;
        }
    }

    //extract the coefficient vector of poly
    int tmp_coef;
    for (int i = 0; i <= deg(poly); i++)
    {
        tmp_coef = conv<long>(poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec[i] = tmp_coef;
    }

    for (int i = 0; i <= deg(inv_poly); i++)
    {
        tmp_coef = conv<long>(inv_poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec_inv[i] = tmp_coef;
    }
}



void get_invertible_vector(vector<int>& vec, vector<int>& vec_inv, int scale,int N,int Q)
{
    //polynomial with the coefficient vector vec (will be generated later)
    ZZ_pX poly;
    ZZ_p coef;
    coef.init(ZZ(Q));
    //the inverse of poly modulo poly_mod (will be generated later)
    ZZ_pX inv_poly;
    //random sampling

    // X^N + 1
    ZZ_pX modulus;
    SetCoeff(modulus, 0, 1); 
    SetCoeff(modulus, N, 1); 
    int count=1;
    while (true)
    {
        //create the polynomial with the coefficient vector of the desired form
        SetCoeff(poly, 0, ternary_sampler(rand_engine)*scale+1);
        for (size_t i = 1; i < vec.size(); i++)
        {
            coef = ternary_sampler(rand_engine)*scale;
            SetCoeff(poly, i, coef);
        }
        //test invertibility
        try
        {
            InvMod(inv_poly, poly, modulus);
            break;
        }
        catch(...)
        {
            // cout << "Polynomial " << poly << " isn't a unit" << endl;
            count++;
            cout<<count<<endl;
            continue;
        }
    }

    // cout<<poly<<endl;

    //extract the coefficient vector of poly
    int tmp_coef;
    for (int i = 0; i <= deg(poly); i++)
    {
        tmp_coef = conv<long>(poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec[i] = tmp_coef;
    }

    for (int i = 0; i <= deg(inv_poly); i++)
    {
        tmp_coef = conv<long>(inv_poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec_inv[i] = tmp_coef;
    }
}


void get_invertible_poly(vector<int>& vec_inv, vector<int>& vec,int N,int Q)
{
    //polynomial with the coefficient vector vec (will be generated later)
    ZZ_pX poly;
    ZZ_p coef;
    coef.init(ZZ(Q));
    //the inverse of poly modulo poly_mod (will be generated later)
    ZZ_pX inv_poly;
    //random sampling

    // X^N + 1
    ZZ_pX modulus;
    SetCoeff(modulus, 0, 1); 
    SetCoeff(modulus, N, 1); 
    int count=1;
    for (size_t i = 0; i < vec.size(); i++)
        {
            coef = vec[i];
            SetCoeff(poly, i, coef);
        }
    InvMod(inv_poly, poly, modulus);
    //extract the coefficient vector of poly
    int tmp_coef;
    for (int i = 0; i <= deg(inv_poly); i++)
    {
        tmp_coef = conv<long>(inv_poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec_inv[i] = tmp_coef;
    }
}



void get_ternary_vector(vector<int>& vec)
{
    for(int i=0; i<vec.size(); i++)
        vec[i] = ternary_sampler(rand_engine);
}

void get_binary_vector(vector<int>& vec){
    for(int i=0; i<vec.size(); i++)
    vec[i] = binary_sampler(rand_engine);
}




int errestimator(const ModQPoly& v,int std){
    int erracc=0;
    for(int i=0;i<v.size();i++){
        int dev=abs(v[i]-std);
        erracc+=dev;
    }
    int errdev=double(erracc)/double(v.size());
    return errdev;
}



void print_hamin(vector<int> vec){
    int count=0;
    for(int i=0;i<vec.size();i++){
        if(vec[i]!=0){
            count++;
        }
    }
    cout<<"hamin weight="<<count<<endl;
}


void external_product(vector<long>& res, const vector<int>& poly, const vector<FFTPoly>& poly_vector, int b, int shift, int l,int N,int N2p1)
{ 
    ModQPoly poly_sign(N);
    ModQPoly poly_abs(N);
    vector<int> poly_decomp(N);

    for (int i = 0; i < N; ++i)
    {
        const int& polyi = poly[i];
        poly_abs[i] = abs(polyi);
        poly_sign[i] = (polyi < 0)? -1 : 1;
    }
    FFTPoly res_fft(N2p1);
    FFTPoly tmp_fft(N2p1);
    int mask = b-1;
    int bound = b >> 1;
    int digit, sgn;
    for (int j = 0; j < l; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            int& abs_val = poly_abs[i];
            //digit =low b bit of abs_poly[i]
            digit = abs_val & mask;
            if (digit > bound)
            {
                poly_decomp[i] = (poly_sign[i] == 1) ? (digit - b): (b - digit);
                abs_val >>= shift;
                ++abs_val;
            }
            else
            {
                poly_decomp[i] = (poly_sign[i] == 1) ? digit: -digit;
                abs_val >>= shift;
            }
        }
        fftN.to_fft(tmp_fft, poly_decomp);
        tmp_fft *= poly_vector[j];
        res_fft += tmp_fft;
    }
    fftN.from_fft(res, res_fft);
}

void mult_fft_poly_by_int(FFTPoly& a, const int b){
    for(int i = 0; i < a.size(); i++)
        a[i] *= b;
}





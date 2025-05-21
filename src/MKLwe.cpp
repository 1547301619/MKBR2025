#include "MKLwe.h"

MKLweSample::MKLweSample(int q,int n,int parties):
q(q),n(n),parties(parties)
{
    this->a=new int[n*parties];
    this->b=0;
}

MKLweSample::MKLweSample(){};

MKLweSample::~MKLweSample() {}

MKLweSample& MKLweSample::operator=(const MKLweSample& other){
        if (this == &other) { 
            return *this;
        }
        if (a != nullptr) {
            delete[] a;
        }
        n = other.n;
        parties = other.parties;
        b = other.b;

        a = new int[n * parties]; 
        for (int i = 0; i < n * parties; ++i) {
            a[i] = other.a[i]; 
        }
        return *this; 
}

MKLweSample MKLweSample::operator+(const MKLweSample& other){
    int q_half=q/2;
    MKLweSample res(q,n,parties);

    for(int i=0;i<n*parties;i++){
        res.a[i]=mod_q_lwe(a[i]+other.a[i],q,q_half);
    }
    res.b=mod_q_lwe(b+other.b,q,q_half);
    return res;
}

MKLweSample MKLweSample::operator-(const MKLweSample& other){
    int q_half=q/2;
    MKLweSample res(q,n,parties);
    for(int i=0;i<n*parties;i++){
        res.a[i]=mod_q_lwe(a[i]-other.a[i],q,q_half);
    }
    res.b=mod_q_lwe(b+other.b,q,q_half);
    return res;
}

MKLweSample operator-(const int c, const MKLweSample& other){
    int q=other.q;
    int n=other.n;
    int parties=other.parties;
    int q_half=q/2;
    MKLweSample res(q,n,parties);
    for(int i=0;i<n*parties;i++){
        res.a[i]=mod_q_lwe(-other.a[i],q,q_half);
    }
    res.b=mod_q_lwe(c-other.b,q,q_half);
    return res;   
}




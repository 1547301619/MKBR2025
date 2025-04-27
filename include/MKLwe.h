#ifndef MKLWE
#define MKLWE

#include <iostream>
#include "poly.h"

using namespace std;

class MKLweSample{
    public:
    int q;
    int b;
    int *a;
    int n;
    int parties;
    

    MKLweSample();
    MKLweSample(int q,int n,int parties);
    ~MKLweSample();

    MKLweSample& operator=(const MKLweSample& other);
    MKLweSample operator+(const MKLweSample& other);
    MKLweSample operator-(const MKLweSample& other);

};

MKLweSample operator-(const int c, const MKLweSample& other);




#endif
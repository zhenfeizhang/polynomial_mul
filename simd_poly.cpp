/*
 ============================================================================
 Name        : simd_poly.c
 Author      : zhenfei
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <math.h>       /* sqrt */
#include <NTL/ZZ.h>
#include <immintrin.h>

//#include "poly.h"
using namespace std;
NTL_CLIENT;


void print256_num(__m256i var)
{
    uint16_t *val = (uint16_t*) &var;
    printf("%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
           val[0], val[1], val[2], val[3],
           val[4], val[5], val[6], val[7],
           val[8], val[9], val[10], val[11],
           val[12], val[13], val[14], val[15]);
}


void
grade_school_mul(
    uint16_t        *res1,  /* out - a * b in Z[x], must be length 2N */
    uint16_t const  *a,     /*  in - polynomial */
    uint16_t const  *b,     /*  in - polynomial */
    uint16_t const   N)     /*  in - number of coefficients in a and b */
{
    uint16_t i;
    uint16_t j;

    for(j=0; j<N; j++)
    {
        res1[j] = a[0]*b[j];
    }
    for(i=1; i<N; i++)
    {
        res1[i+N-1] = 0;
        for(j=0; j<N; j++)
        {
            res1[i+j] += a[i]*b[j];
        }
    }
    res1[2*N-1] = 0;

    return;
}

void
__m256i_grade_school_mul_16(
    uint16_t        *res1,  /* out - a * b in Z[x], must be length 2N */
    uint16_t        *buf,   /* buf size >= 32 bytes */
    uint16_t const  *a,     /*  in - polynomial */
    uint16_t const  *b,     /*  in - polynomial */
    uint16_t const   N)     /*  in - number of coefficients in a and b <= 16*/
{
    uint16_t i;
    uint16_t j;
    uint16_t n;
    uint16_t *buf1;
    memset(buf, 0, 32*sizeof(uint16_t));
    memcpy(buf, a, N*sizeof(uint16_t));
    buf1 = buf + 16;
    __m256i a256low;
    __m256i a256high;
    __m256i r256low;
    __m256i r256high;
    __m256i tmp;

    a256low = _mm256_loadu_si256((__m256i *)(buf));
    r256high = _mm256_set1_epi16(0);
    tmp = _mm256_set1_epi16(b[0]);
    r256low = _mm256_mullo_epi16 (tmp,a256low);


    for (i=1;i<N;i++)
    {
        tmp = _mm256_set1_epi16(b[i]);
        memcpy(buf+i, a, N*sizeof(uint16_t));       // shifting
        memset(buf, 0, i*sizeof(uint16_t));
        a256low = _mm256_loadu_si256((__m256i *)(buf));
        a256high = _mm256_loadu_si256((__m256i *)(buf1));
        r256low =_mm256_add_epi16 (r256low, _mm256_mullo_epi16 (tmp,a256low));
        r256high = _mm256_add_epi16 (r256high,_mm256_mullo_epi16 (tmp,a256high));
   }
    _mm256_storeu_si256((__m256i *)(res1+16), r256high);
    _mm256_storeu_si256((__m256i *)res1, (r256low));

    return;
}

void
__m256i_grade_school_mul_32(
    uint16_t        *res1,  /* out - a * b in Z[x], must be length 2N */
    uint16_t        *buf,   /* buf size >= 64 bytes */
    uint16_t const  *a,     /*  in - polynomial */
    uint16_t const  *b,     /*  in - polynomial */
    uint16_t const   N)     /*  in - number of coefficients in a and b <= 32*/
{
    uint16_t i;
    uint16_t j;
    uint16_t n;
    uint16_t *buf1, *buf2, *buf3;
    memset(buf, 0, 64*sizeof(uint16_t));
    memcpy(buf, a, N*sizeof(uint16_t));
    buf1 = buf + 16;
    buf2 = buf + 32;
    buf3 = buf + 48;

    __m256i a256[4];
    __m256i r256[4];
    __m256i tmp;

    a256[0] = _mm256_loadu_si256((__m256i *)(buf));
    a256[1] = _mm256_loadu_si256((__m256i *)(buf1));

    tmp = _mm256_set1_epi16(b[0]);
    r256[0] = _mm256_mullo_epi16 (tmp,a256[0]);
    r256[1] = _mm256_mullo_epi16 (tmp,a256[1]);
    r256[2] = _mm256_set1_epi16(0);
    r256[3] = _mm256_set1_epi16(0);

    for (i=1;i<N;i++)
    {
        tmp = _mm256_set1_epi16(b[i]);
        memcpy(buf+i, a, N*sizeof(uint16_t));       // shifting
        memset(buf, 0, i*sizeof(uint16_t));
        a256[0] = _mm256_loadu_si256((__m256i *)(buf));
        a256[1] = _mm256_loadu_si256((__m256i *)(buf1));
        a256[2] = _mm256_loadu_si256((__m256i *)(buf2));
        a256[3] = _mm256_loadu_si256((__m256i *)(buf3));
        r256[0] = _mm256_add_epi16 (r256[0], _mm256_mullo_epi16 (tmp,a256[0]));
        r256[1] = _mm256_add_epi16 (r256[1], _mm256_mullo_epi16 (tmp,a256[1]));
        r256[2] = _mm256_add_epi16 (r256[2], _mm256_mullo_epi16 (tmp,a256[2]));
        r256[3] = _mm256_add_epi16 (r256[3], _mm256_mullo_epi16 (tmp,a256[3]));
   }

    _mm256_storeu_si256((__m256i *)res1, r256[0]);
    _mm256_storeu_si256((__m256i *)(res1+16), r256[1]);
    _mm256_storeu_si256((__m256i *)(res1+32), r256[2]);
    _mm256_storeu_si256((__m256i *)(res1+48), r256[3]);
    return;
}

static void
karatsuba(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n, /*  in - number of coefficients in a and b */
    uint16_t const   k) /*  in - degree of schoolbook multiplication */
{
    if (n <= k)
    {
        grade_school_mul(r, a, b, n);
        return;
    }
    uint16_t i;
    uint16_t s = n/2;
    uint16_t const *a1 = a+s;
    uint16_t const *b1 = b+s;
    uint16_t *t1 = t+s;
    uint16_t *r1 = r+s, *r2 = r+2*s, *r3 = r+3*s;

    __m256i m[6];
    for (i=0; i<s; i+=16) {
      m[0] = _mm256_loadu_si256((__m256i *)(a+i));
      m[1] = _mm256_loadu_si256((__m256i *)(a1+i));
      m[2] = _mm256_loadu_si256((__m256i *)(b1+i));
      m[3] = _mm256_loadu_si256((__m256i *)(b+i));
      m[0] = _mm256_sub_epi16(m[0], m[1]);
      m[1] = _mm256_sub_epi16(m[2], m[3]);
      _mm256_storeu_si256((__m256i *)(r+i), m[0]);
      _mm256_storeu_si256((__m256i *)(r1+i), m[1]);
    }

    karatsuba(t, r2, r, r1, s, k);
    karatsuba(r2, r, a1, b1, s, k);

    for (i=0; i<s; i+=16) {
      m[0] = _mm256_loadu_si256((__m256i *)(r2+i));
      m[1] = _mm256_loadu_si256((__m256i *)(t+i));
      m[2] = _mm256_loadu_si256((__m256i *)(r3+i));
      m[3] = _mm256_loadu_si256((__m256i *)(t1+i));
      m[1] = _mm256_add_epi16(m[1], m[0]);
      m[2] = _mm256_add_epi16(m[2], m[0]);
      m[2] = _mm256_add_epi16(m[2], m[3]);
      _mm256_storeu_si256((__m256i *)(r1+i), m[1]);
      _mm256_storeu_si256((__m256i *)(r2+i), m[2]);
    }

    karatsuba(t, r, a, b, s, k);

    for (i=0; i<s; i++) {
      r[i] = t[i];
    }
    for (i=0; i<s; i+=16) {
      m[0] = _mm256_loadu_si256((__m256i *)(r1+i));
      m[1] = _mm256_loadu_si256((__m256i *)(t+i));
      m[2] = _mm256_loadu_si256((__m256i *)(t1+i));
      m[3] = _mm256_loadu_si256((__m256i *)(r2+i));
      m[0] = _mm256_add_epi16(m[0], m[1]);
      m[0] = _mm256_add_epi16(m[0], m[2]);
      m[3] = _mm256_add_epi16(m[3], m[2]);
      _mm256_storeu_si256((__m256i *)(r1+i), m[0]);
      _mm256_storeu_si256((__m256i *)(r2+i), m[3]);
    }

    return;
}

static void
karatsuba_32(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n, /*  in - number of coefficients in a and b */
    uint16_t const   k) /*  in - degree of schoolbook multiplication */
{
    if (n <= 32)
    {
        __m256i_grade_school_mul_32(r, t, a, b, n);
        return;
    }
    uint16_t i;
    uint16_t s = n/2;
    uint16_t const *a1 = a+s;
    uint16_t const *b1 = b+s;
    uint16_t *t1 = t+s;
    uint16_t *r1 = r+s, *r2 = r+2*s, *r3 = r+3*s;

    __m256i m[6];
    for (i=0; i<s; i+=16) {
      m[0] = _mm256_loadu_si256((__m256i *)(a+i));
      m[1] = _mm256_loadu_si256((__m256i *)(a1+i));
      m[2] = _mm256_loadu_si256((__m256i *)(b1+i));
      m[3] = _mm256_loadu_si256((__m256i *)(b+i));
      m[0] = _mm256_sub_epi16(m[0], m[1]);
      m[1] = _mm256_sub_epi16(m[2], m[3]);
      _mm256_storeu_si256((__m256i *)(r+i), m[0]);
      _mm256_storeu_si256((__m256i *)(r1+i), m[1]);
    }

    karatsuba(t, r2, r, r1, s, k);
    karatsuba(r2, r, a1, b1, s, k);

    for (i=0; i<s; i+=16) {
      m[0] = _mm256_loadu_si256((__m256i *)(r2+i));
      m[1] = _mm256_loadu_si256((__m256i *)(t+i));
      m[2] = _mm256_loadu_si256((__m256i *)(r3+i));
      m[3] = _mm256_loadu_si256((__m256i *)(t1+i));
      m[1] = _mm256_add_epi16(m[1], m[0]);
      m[2] = _mm256_add_epi16(m[2], m[0]);
      m[2] = _mm256_add_epi16(m[2], m[3]);
      _mm256_storeu_si256((__m256i *)(r1+i), m[1]);
      _mm256_storeu_si256((__m256i *)(r2+i), m[2]);
    }

    karatsuba(t, r, a, b, s, k);

    for (i=0; i<s; i++) {
      r[i] = t[i];
    }
    for (i=0; i<s; i+=16) {
      m[0] = _mm256_loadu_si256((__m256i *)(r1+i));
      m[1] = _mm256_loadu_si256((__m256i *)(t+i));
      m[2] = _mm256_loadu_si256((__m256i *)(t1+i));
      m[3] = _mm256_loadu_si256((__m256i *)(r2+i));
      m[0] = _mm256_add_epi16(m[0], m[1]);
      m[0] = _mm256_add_epi16(m[0], m[2]);
      m[3] = _mm256_add_epi16(m[3], m[2]);
      _mm256_storeu_si256((__m256i *)(r1+i), m[0]);
      _mm256_storeu_si256((__m256i *)(r2+i), m[3]);
    }

    return;
}


int main()
{

    uint16_t N;     // dimension
    uint16_t q;     // q = 2048
    uint16_t *a;    // first polynomial
    uint16_t *b;    // second polynomial
    uint64_t *a64;    // first polynomial
    uint64_t *b64;    // second polynomial
    uint64_t *buf64;  // buffer
    uint16_t *buf;    // result
    uint16_t *r;    // result
    uint64_t *r64;    // result
    uint16_t i,j;
    uint16_t test_dim;
    __m256   *a256;
    __m256   *b256;
    __m256   *r256;

    clock_t start, end;
    double karat, scht, toomt,toomtrec;
    ofstream fout;
    float ss1,ss2;
    fout.open("result.txt");
    N = 256;
    a = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    b = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    buf = (uint16_t*) malloc (4*N*sizeof(uint16_t));
    r = (uint16_t*) malloc (4*N*sizeof(uint64_t));

    r64 = (uint64_t*)r;
    fout<<test_dim<<" ";

    test_dim = 1000;



    for (test_dim=8;test_dim<33;test_dim++)
    {
        ss1 = 0;
        ss2 = 0;
        cout<<"dimension: "<<test_dim<<" ";
        for (j=0;j<1000;j++)
        {
            for(i=0; i< test_dim;i++)
            {
                a[i] = rand()&0x07FF;
                b[i] = rand()&0x07FF;
            }


            start = clock();
            grade_school_mul(r, a, b, test_dim);
            end = clock();
            ss1 += (float)(end-start);


            start = clock();
            __m256i_grade_school_mul_32(r, buf,a, b, test_dim);
            end = clock();
            ss2 += (float)(end-start);
        }
        cout<<ss1<<" "<<ss2<<" "<<endl;
    }



    for (test_dim=33;test_dim<500;test_dim++)
    {
        ss1 = 0;
        ss2 = 0;
        cout<<"dimension: "<<test_dim<<" ";
        for (j=0;j<1000;j++)
        {
            for(i=0; i< test_dim;i++)
            {
                a[i] = rand()&0x07FF;
                b[i] = rand()&0x07FF;
            }


            start = clock();
            karatsuba(r, buf, a, b, test_dim, 32);
            end = clock();
            ss1 += (float)(end-start);


            start = clock();
            karatsuba_32(r, buf, a, b, test_dim, 32);
            end = clock();
            ss2 += (float)(end-start);
        }
        cout<<ss1<<" "<<ss2<<" "<<endl;
    }


    return 0;

}

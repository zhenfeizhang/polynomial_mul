#ifndef _SIMD_POLY_H_
#define _SIMD_POLY_H_
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
using namespace std;

void print256_num(__m256i var);


/* ******************** school book area ************************ */

void
grade_school_mul(
    uint16_t        *res1,  /* out - a * b in Z[x], must be length 2N */
    uint16_t const  *a,     /*  in - polynomial */
    uint16_t const  *b,     /*  in - polynomial */
    uint16_t const   N);    /*  in - number of coefficients in a and b */

void
__m256i_grade_school_mul_16(
    uint16_t        *res1,  /* out - a * b in Z[x], must be length 2N */
    uint16_t        *buf,   /* buf size >= 32 bytes */
    uint16_t const  *a,     /*  in - polynomial */
    uint16_t const  *b,     /*  in - polynomial */
    uint16_t const   N);    /*  in - number of coefficients in a and b <= 16*/


void
__m256i_grade_school_mul_32(
    uint16_t        *res1,  /* out - a * b in Z[x], must be length 2N */
    uint16_t        *buf,   /* buf size >= 64 bytes */
    uint16_t const  *a,     /*  in - polynomial */
    uint16_t const  *b,     /*  in - polynomial */
    uint16_t const   N);    /*  in - number of coefficients in a and b <= 32*/


int test_SB_32();


/* ******************** karatsuba area ************************ */

static void
__mm256i_karatsuba_SB(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n);/*  in - number of coefficients in a and b */


static void
__mm256i_karatsuba__mm256_SB(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n);/*  in - number of coefficients in a and b */
int test_karatsuba();

/* ******************** toom 3 area ************************ */

/*
 * toom3 multiplication with uint16_t coefficients,
 * only usable for degree less than 96,
 * uses uint16_t based school book multiplications when
 * degree drops below 32
 * okay speed
 */
int
toom3(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n);/*  in - number of coefficients in a and b */



/*
 * toom3 multiplication with uint16_t coefficients,
 * only usable for degree less than 96,
 * uses vectorized school book multiplications (avx2)
 * when degree drops below 32
 * good speed
 */
int
toom3__mm256i_SB(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n);/*  in - number of coefficients in a and b */


/*
 * vectorized toom3 multiplication with avx2 coefficients,
 * only usable for degree less than 96,
 * uses vectorized school book multiplications (avx2)
 * when degree drops below 32
 * fastest
 */
int
__mm256i_toom3__mm256i_SB(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n);/*  in - number of coefficients in a and b */


int test_toom3();
#endif /*SIMD_POLY_H_*/

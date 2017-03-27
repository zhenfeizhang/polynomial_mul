#include "simd_poly.h"



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
    if (n <= 31)
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

    karatsuba_32(t, r2, r, r1, s, k);
    karatsuba_32(r2, r, a1, b1, s, k);

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

    karatsuba_32(t, r, a, b, s, k);

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


void
__mm256_karatsuba_768(
        uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
        uint16_t        *t, /*  in - n coefficients of scratch space */
        uint16_t const  *a, /*  in - polynomial */
        uint16_t const  *b, /*  in - polynomial */
        uint16_t const   n) /*  in - number of coefficients in a and b */
{
    __m256i ax[48], bx[48], rx[96],tmp;

    uint16_t i,j,k;

    for(i=0;i<48;i++)
    {
        ax[i] = _mm256_loadu_si256((__m256i *)(a+i*16));
    }

    for (i=0;i<48;i++)
    {
        print256_num(ax[i]);
        printf("\n");
        for (j=0;j<16;j++)
            printf ("%i ", a[i*16+j]);
        printf("\n");
    }

}

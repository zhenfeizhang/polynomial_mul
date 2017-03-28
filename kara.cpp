#include "simd_poly.h"



static void
__mm256i_karatsuba_SB(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n) /*  in - number of coefficients in a and b */
{
    if (n < 32)
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

    __mm256i_karatsuba_SB(t, r2, r, r1, s);
    __mm256i_karatsuba_SB(r2, r, a1, b1, s);

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

    __mm256i_karatsuba_SB(t, r, a, b, s);

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
__mm256i_karatsuba__mm256_SB(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n) /*  in - number of coefficients in a and b */
{
    if (n < 32)
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

    __mm256i_karatsuba__mm256_SB(t, r2, r, r1, s);
    __mm256i_karatsuba__mm256_SB(r2, r, a1, b1, s);

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

    __mm256i_karatsuba__mm256_SB(t, r, a, b, s);

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

int test_karatsuba()
{
    uint16_t N;     // dimension
    uint16_t *a;    // first polynomial
    uint16_t *b;    // second polynomial
    uint16_t *buf;    // buffer
    uint16_t *r;    // result
    uint16_t *r2;    // result
    uint16_t *r3;    // result
    uint16_t i,j;
    uint16_t test_dim;
    float ss1,ss2, ss3;
    clock_t start, end;

    N = 96;
    a = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    b = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    buf = (uint16_t*) malloc (4*N*sizeof(uint16_t));
    r = (uint16_t*) malloc (4*N*sizeof(uint64_t));
    r2 = (uint16_t*) malloc (4*N*sizeof(uint64_t));
    r3 = (uint16_t*) malloc (4*N*sizeof(uint64_t));
    cout<<"testing karatsuba 3"<<endl;
    for (test_dim=65;test_dim<96;test_dim++)
    {
        ss1 = 0;
        ss2 = 0;
        ss3 = 0;
        cout<<"dimension: "<<test_dim<<" ";
        for (j=0;j<1000;j++)
        {
            memset(a+test_dim, 0, 2*N*sizeof(uint16_t));
            memset(b+test_dim, 0, 2*N*sizeof(uint16_t));
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
            __mm256i_karatsuba_SB(r2, buf, a, b, test_dim);
            end = clock();
            ss2 += (float)(end-start);


            start = clock();
            __mm256i_karatsuba__mm256_SB(r3, buf, a, b, test_dim);
            end = clock();
            ss3 += (float)(end-start);

            for (i=0;i<test_dim*2-1;i++)
            {
                if ((r[i]%2048)!=(r2[i]%2048))
                {
                    printf("error\n");
                    for (j=0;j<test_dim*2-1;j++)
                    {
                        printf("%d %d %d %d\n", j, r[j],r2[j], r[j]-r2[j]);

                    }
                    return 1;
                }
            }
        }
        cout<<ss1<<" "<<ss2<<" "<<ss3<<endl;
    }
    return 0;
}

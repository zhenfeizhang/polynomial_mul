#include "simd_poly.h"
int
toom3(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n) /*  in - number of coefficients in a and b */
{
    if (n>96)
    {
        printf("degree exceeds the maximum (96) allowed\n");
        return -1;
    }
    if (n<=32)
    {
        grade_school_mul(r, a, b, n);
        return 0;
    }

    uint16_t    i;
    uint16_t    s = 32;
    uint16_t    s2 = 64;
    uint16_t  const  *a1 = a+s, *a2 = a+s2;
    uint16_t  const  *b1 = b+s, *b2 = b+s2;
    uint16_t    *r1 = r+s, *r2 = r+s2, *r3 = r+s*3, *r4 = r+s*4, *r5 = r+s*5;
    uint16_t    *t2 = t+s2, *t4 = t+s2*2, *t6 = t+s2*3, *t8 = t+s2*4;

    /*
     * t  = w0   = a(0) * b(0)
     * t8 = w4   = a(inf) * b(inf)
     *           = a2     * b2
     */

    grade_school_mul(t, a, b, s);
    grade_school_mul(t8, a2, b2, s);

#ifdef DEBUG
    printf("t:");
    print32poly(t);
    print32poly(t+32);
    printf("t8:");
    print32poly(t8);
    print32poly(t8+32);
#endif
    /*
     * t2 = a(1) *b(1)
     *    = (a0+a1+a2)*(b0+b1+b2)
     *    = r         * r2
     * t4 = a(-1) * b(-1)
     *    = (a0-a1+a2)*(b0-b1+b2)
     *    = r1        * r3
     */
    for (i=0;i<s;i++)
    {
        r[i]  = a[i]  + a2[i];
        r1[i] = r[i]  - a1[i];
        r[i]  = r[i]  + a1[i];
        r2[i] = b[i]  + b2[i];
        r3[i] = r2[i] - b1[i];
        r2[i] = r2[i] + b1[i];
    }
    grade_school_mul(t2, r, r2, s);
    grade_school_mul(t4, r1, r3, s);

#ifdef DEBUG
    printf("t2:");
    print32poly(t2);
    print32poly(t2+32);
    printf("t4:");
    print32poly(t4);
    print32poly(t4+32);
#endif


    /*
     * r  = (t2+t4)/2 - w0 - w4
     *    = w2
     * r2 = (t2-t4)/2
     *    = w1 + w3
     */
    for (i=0;i<s2;i++)
    {
        r[i]  = (t2[i]  + t4[i])/2 - t[i] - t8[i];
        r2[i] = (t2[i]  - t4[i])/2;
    }

#ifdef DEBUG
    printf("r:");
    print32poly(r);
    print32poly(r+32);
    printf("r2:");
    print32poly(r2);
    print32poly(r2+32);
#endif

    /*
     * t6 = a(2) *b(2)
     *    = (a0+2a1+4a2)*(b0+2b1+4b2)
     *    = r4          * r5
     */
    for (i=0;i<s;i++)
    {
        r4[i]  = a[i]  + 2*a1[i] + 4*a2[i];
        r5[i]  = b[i]  + 2*b1[i] + 4*b2[i];
    }
    grade_school_mul(t6, r4, r5, s);

    /*
     * t6 = w1 + 4*w3
     *    = (t6 - w0 - 4*w2 - 16*w4)/2
     *    = (t6 - t0 - 4*r  - 16*t8)/2
     */

    for (i=0;i<s2;i++)
    {
        t6[i]  = (t6[i] - t[i]-4*r[i]-16*t8[i])/2;
    }


#ifdef DEBUG
    printf("t6:");
    print32poly(t6);
    print32poly(t6+32);
#endif

    /*
     * t2 = w3 = (t6 - r2)/3
     *    = (t6 - t4) * 43691
     * t4 = w1 = (4*r2-t6)/3
     *    = (4*r4-t6) * 43691
     */
    for (i=0;i<s2;i++)
    {
        t2[i]  = (t6[i] - r2[i])*43691;
        t4[i]  = ((4*r2[i] - t6[i])*43691);
    }
#ifdef DEBUG
    printf("t2:");
    print32poly(t2);
    print32poly(t2+32);
    printf("t4:");
    print32poly(t4);
    print32poly(t4+32);
#endif

    /*
     * now we have
     *  t0 = w0
     *  t4 = w1
     *  r  = w2
     *  t2 = w3
     *  t8 = w4
     * putting them back
     */

    memcpy(r2, r,  sizeof(uint16_t)*s2);
    memcpy(r,  t,  sizeof(uint16_t)*s2);
    memcpy(r4, t8, sizeof(uint16_t)*s2);


    for (i=0;i<s2;i++)
    {
        r1[i] += t4[i];
        r3[i] += t2[i];
    }

    return 0;
}



int
toom3__mm256i_SB(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n) /*  in - number of coefficients in a and b */
{
    if (n>96)
    {
        printf("degree exceeds the maximum (96) allowed\n");
        return -1;
    }
    if (n<=32)
    {
        grade_school_mul(r, a, b, n);
        return 0;
    }

    uint16_t    i;
    uint16_t    s = 32, s2 = 64;
    uint16_t  const  *a1 = a+s, *a2 = a+s2;
    uint16_t  const  *b1 = b+s, *b2 = b+s2;
    uint16_t    *r1 = r+s,  *r2 = r+s2,   *r3 = r+s*3,  *r4 = r+s*4, *r5 = r+s*5;
    uint16_t    *t2 = t+s2, *t4 = t+s2*2, *t6 = t+s2*3, *t8 = t+s2*4;
    uint16_t    *buf = t+s2*10;


    /*
     * t  = w0   = a(0) * b(0)
     * t8 = w4   = a(inf) * b(inf)
     *           = a2     * b2
     */

//    grade_school_mul(t, a, b, s);
//    grade_school_mul(t8, a2, b2, s);

    __m256i_grade_school_mul_32(t, buf, a, b, s);
    __m256i_grade_school_mul_32(t8, buf, a2, b2, s);
#ifdef DEBUG
    printf("t:");
    print32poly(t);
    print32poly(t+32);
    printf("t8:");
    print32poly(t8);
    print32poly(t8+32);
#endif
    /*
     * t2 = a(1) *b(1)
     *    = (a0+a1+a2)*(b0+b1+b2)
     *    = r         * r2
     * t4 = a(-1) * b(-1)
     *    = (a0-a1+a2)*(b0-b1+b2)
     *    = r1        * r3
     */
    for (i=0;i<s;i++)
    {
        r[i]  = a[i]  + a2[i];
        r1[i] = r[i]  - a1[i];
        r[i]  = r[i]  + a1[i];
        r2[i] = b[i]  + b2[i];
        r3[i] = r2[i] - b1[i];
        r2[i] = r2[i] + b1[i];
    }
//    grade_school_mul(t2, r, r2, s);
//    grade_school_mul(t4, r1, r3, s);

    __m256i_grade_school_mul_32(t2, buf,r, r2, s);
    __m256i_grade_school_mul_32(t4, buf,r1, r3, s);

#ifdef DEBUG
    printf("t2:");
    print32poly(t2);
    print32poly(t2+32);
    printf("t4:");
    print32poly(t4);
    print32poly(t4+32);
#endif


    /*
     * r  = (t2+t4)/2 - w0 - w4
     *    = w2
     * r2 = (t2-t4)/2
     *    = w1 + w3
     */
    for (i=0;i<s2;i++)
    {
        r[i]  = (t2[i]  + t4[i])/2 - t[i] - t8[i];
        r2[i] = (t2[i]  - t4[i])/2;
    }

#ifdef DEBUG
    printf("r:");
    print32poly(r);
    print32poly(r+32);
    printf("r2:");
    print32poly(r2);
    print32poly(r2+32);
#endif

    /*
     * t6 = a(2) *b(2)
     *    = (a0+2a1+4a2)*(b0+2b1+4b2)
     *    = r4          * r5
     */
    for (i=0;i<s;i++)
    {
        r4[i]  = a[i]  + 2*a1[i] + 4*a2[i];
        r5[i]  = b[i]  + 2*b1[i] + 4*b2[i];
    }
//    grade_school_mul(t6, r4, r5, s);
    __m256i_grade_school_mul_32(t6, buf, r4, r5, s);;

    /*
     * t6 = w1 + 4*w3
     *    = (t6 - w0 - 4*w2 - 16*w4)/2
     *    = (t6 - t0 - 4*r  - 16*t8)/2
     */

    for (i=0;i<s2;i++)
    {
        t6[i]  = (t6[i] - t[i]-4*r[i]-16*t8[i])/2;
    }


#ifdef DEBUG
    printf("t6:");
    print32poly(t6);
    print32poly(t6+32);
#endif

    /*
     * t2 = w3 = (t6 - r2)/3
     *    = (t6 - t4) * 43691
     * t4 = w1 = (4*r2-t6)/3
     *    = (4*r4-t6) * 43691
     */
    for (i=0;i<s2;i++)
    {
        t2[i]  = (t6[i] - r2[i])*43691;
        t4[i]  = ((4*r2[i] - t6[i])*43691);
    }
#ifdef DEBUG
    printf("t2:");
    print32poly(t2);
    print32poly(t2+32);
    printf("t4:");
    print32poly(t4);
    print32poly(t4+32);
#endif

    /*
     * now we have
     *  t0 = w0
     *  t4 = w1
     *  r  = w2
     *  t2 = w3
     *  t8 = w4
     * putting them back
     */

    memcpy(r2, r,  sizeof(uint16_t)*s2);
    memcpy(r,  t,  sizeof(uint16_t)*s2);
    memcpy(r4, t8, sizeof(uint16_t)*s2);

    for (i=0;i<s2;i++)
    {
        r1[i] += t4[i];
        r3[i] += t2[i];
    }

    return 0;
}




int
__mm256i_toom3__mm256i_SB(
    uint16_t        *r, /* out - a * b in Z[x], must be length 2n */
    uint16_t        *t, /*  in - n coefficients of scratch space */
    uint16_t const  *a, /*  in - polynomial */
    uint16_t const  *b, /*  in - polynomial */
    uint16_t const   n) /*  in - number of coefficients in a and b */
{
    if (n>96)
    {
        printf("degree %d exceeds the maximum (96) allowed\n", n);
        return -1;
    }
    if (n<=32)
    {
        grade_school_mul(r, a, b, n);
        return 0;
    }

    uint16_t    i;
    uint16_t    s = 32, s2 = 64;
    uint16_t  const  *a1 = a+s, *a2 = a+s2;
    uint16_t  const  *b1 = b+s, *b2 = b+s2;
    uint16_t    *r1 = r+s,  *r2 = r+s2,   *r3 = r+s*3,  *r4 = r+s*4, *r5 = r+s*5;
    uint16_t    *t2 = t+s2, *t4 = t+s2*2, *t6 = t+s2*3, *t8 = t+s2*4;
    uint16_t    *buf = t+s2*10;
    __m256i     mr0a, mr0b, mr1a, mr1b, mr2a, mr2b,
                mr3a, mr3b, mr4a, mr4b, mr5a, mr5b;
    __m256i     mt0a, mt0b, mt0c, mt0d,
                mt2a, mt2b, mt2c, mt2d,
                mt4a, mt4b, mt4c, mt4d,
                mt6a, mt6b, mt6c, mt6d,
                mt8a, mt8b, mt8c, mt8d;
    __m256i     tmp;


    /*
     * t  = w0   = a(0) * b(0)
     * t8 = w4   = a(inf) * b(inf)
     *           = a2     * b2
     */

    __m256i_grade_school_mul_32(t, buf, a, b, s);
    __m256i_grade_school_mul_32(t8, buf, a2, b2, s);

    /*
     * t2 = a(1) *b(1)
     *    = (a0+a1+a2)*(b0+b1+b2)
     *    = r         * r2
     * t4 = a(-1) * b(-1)
     *    = (a0-a1+a2)*(b0-b1+b2)
     *    = r1        * r3
     */

    /*
     * r  = (a0+a1+a2)
     * r1 = (a0-a1+a2)
     */

    mr0a = _mm256_loadu_si256((__m256i *)a);
    tmp  = _mm256_loadu_si256((__m256i *)a2);
    mr0a = _mm256_add_epi16 (mr0a, tmp);
    mr0b = _mm256_loadu_si256((__m256i *)a+1);
    tmp  = _mm256_loadu_si256((__m256i *)a2+1);
    mr0b = _mm256_add_epi16 (mr0b, tmp);

    tmp  = _mm256_loadu_si256((__m256i *)a1);
    mr1a = _mm256_sub_epi16 (mr0a, tmp);
    mr0a = _mm256_add_epi16 (mr0a, tmp);
    tmp  = _mm256_loadu_si256((__m256i *)a1+1);
    mr1b = _mm256_sub_epi16 (mr0b, tmp);
    mr0b = _mm256_add_epi16 (mr0b, tmp);

    _mm256_storeu_si256((__m256i *)r, mr0a);
    _mm256_storeu_si256((__m256i *)r+1, mr0b);
    _mm256_storeu_si256((__m256i *)r1, mr1a);
    _mm256_storeu_si256((__m256i *)r1+1, mr1b);

    /*
     * r2 = (b0+b1+b2)
     * r3 = (b0-b1+b2)
     */

    mr2a = _mm256_loadu_si256((__m256i *)b);
    tmp  = _mm256_loadu_si256((__m256i *)b2);
    mr2a = _mm256_add_epi16 (mr2a, tmp);
    mr2b = _mm256_loadu_si256((__m256i *)b+1);
    tmp  = _mm256_loadu_si256((__m256i *)b2+1);
    mr2b = _mm256_add_epi16 (mr2b, tmp);

    tmp  = _mm256_loadu_si256((__m256i *)b1);
    mr3a = _mm256_sub_epi16 (mr2a, tmp);
    mr2a = _mm256_add_epi16 (mr2a, tmp);
    tmp  = _mm256_loadu_si256((__m256i *)b1+1);
    mr3b = _mm256_sub_epi16 (mr2b, tmp);
    mr2b = _mm256_add_epi16 (mr2b, tmp);

    _mm256_storeu_si256((__m256i *)r2, mr2a);
    _mm256_storeu_si256((__m256i *)r2+1, mr2b);
    _mm256_storeu_si256((__m256i *)r3, mr3a);
    _mm256_storeu_si256((__m256i *)r3+1, mr3b);


    __m256i_grade_school_mul_32(t2, buf,r, r2, s);
    __m256i_grade_school_mul_32(t4, buf,r1, r3, s);


    /*
     * r  = (t2+t4)/2 - w0 - w4
     *    = w2
     * r2 = (t2-t4)/2
     *    = w1 + w3
     */

    mt0a = _mm256_loadu_si256((__m256i *)t);
    mt0b = _mm256_loadu_si256((__m256i *)t+1);
    mt0c = _mm256_loadu_si256((__m256i *)t+2);
    mt0d = _mm256_loadu_si256((__m256i *)t+3);

    mt2a = _mm256_loadu_si256((__m256i *)t2);
    mt2b = _mm256_loadu_si256((__m256i *)t2+1);
    mt2c = _mm256_loadu_si256((__m256i *)t2+2);
    mt2d = _mm256_loadu_si256((__m256i *)t2+3);

    mt4a = _mm256_loadu_si256((__m256i *)t4);
    mt4b = _mm256_loadu_si256((__m256i *)t4+1);
    mt4c = _mm256_loadu_si256((__m256i *)t4+2);
    mt4d = _mm256_loadu_si256((__m256i *)t4+3);

    mt8a = _mm256_loadu_si256((__m256i *)t8);
    mt8b = _mm256_loadu_si256((__m256i *)t8+1);
    mt8c = _mm256_loadu_si256((__m256i *)t8+2);
    mt8d = _mm256_loadu_si256((__m256i *)t8+3);

    /*
     * r  = (t2+t4)/2 - w0 - w4
     *    = w2
     */

    mr0a = _mm256_add_epi16 (mt2a, mt4a);
    mr0a = _mm256_srai_epi16(mr0a, 1);
    mr0a = _mm256_sub_epi16 (mr0a, mt0a);
    mr0a = _mm256_sub_epi16 (mr0a, mt8a);

    mr0b = _mm256_add_epi16 (mt2b, mt4b);
    mr0b = _mm256_srai_epi16(mr0b, 1);
    mr0b = _mm256_sub_epi16 (mr0b, mt0b);
    mr0b = _mm256_sub_epi16 (mr0b, mt8b);


    mr1a = _mm256_add_epi16 (mt2c, mt4c);
    mr1a = _mm256_srai_epi16(mr1a, 1);
    mr1a = _mm256_sub_epi16 (mr1a, mt0c);
    mr1a = _mm256_sub_epi16 (mr1a, mt8c);

    mr1b = _mm256_add_epi16 (mt2d, mt4d);
    mr1b = _mm256_srai_epi16(mr1b, 1);
    mr1b = _mm256_sub_epi16 (mr1b, mt0d);
    mr1b = _mm256_sub_epi16 (mr1b, mt8d);


    /*
     * r2 = (t2-t4)/2
     *    = w1 + w3
     */

    mr2a = _mm256_sub_epi16 (mt2a, mt4a);
    mr2a = _mm256_srai_epi16(mr2a, 1);
    mr2b = _mm256_sub_epi16 (mt2b, mt4b);
    mr2b = _mm256_srai_epi16(mr2b, 1);
    mr3a = _mm256_sub_epi16 (mt2c, mt4c);
    mr3a = _mm256_srai_epi16(mr3a, 1);
    mr3b = _mm256_sub_epi16 (mt2d, mt4d);
    mr3b = _mm256_srai_epi16(mr3b, 1);

    /*
     * t6 = a(2) *b(2)
     *    = (a0+2a1+4a2)*(b0+2b1+4b2)
     *    = r4          * r5
     */

    mr4a = _mm256_loadu_si256((__m256i *)a);
    tmp  = _mm256_loadu_si256((__m256i *)a1);
    tmp  = _mm256_slli_epi16(tmp, 1);
    mr4a = _mm256_add_epi16 (mr4a,tmp);
    tmp  = _mm256_loadu_si256((__m256i *)a2);
    tmp  = _mm256_slli_epi16(tmp, 2);
    mr4a = _mm256_add_epi16 (mr4a,tmp);

    mr4b = _mm256_loadu_si256((__m256i *)a+1);
    tmp  = _mm256_loadu_si256((__m256i *)a1+1);
    tmp  = _mm256_slli_epi16(tmp, 1);
    mr4b = _mm256_add_epi16 (mr4b,tmp);
    tmp  = _mm256_loadu_si256((__m256i *)a2+1);
    tmp  = _mm256_slli_epi16(tmp, 2);
    mr4b = _mm256_add_epi16 (mr4b,tmp);

    mr5a = _mm256_loadu_si256((__m256i *)b);
    tmp  = _mm256_loadu_si256((__m256i *)b1);
    tmp  = _mm256_slli_epi16(tmp, 1);
    mr5a = _mm256_add_epi16 (mr5a,tmp);
    tmp  = _mm256_loadu_si256((__m256i *)b2);
    tmp  = _mm256_slli_epi16(tmp, 2);
    mr5a = _mm256_add_epi16 (mr5a,tmp);

    mr5b = _mm256_loadu_si256((__m256i *)b+1);
    tmp  = _mm256_loadu_si256((__m256i *)b1+1);
    tmp  = _mm256_slli_epi16(tmp, 1);
    mr5b = _mm256_add_epi16 (mr5b,tmp);
    tmp  = _mm256_loadu_si256((__m256i *)b2+1);
    tmp  = _mm256_slli_epi16(tmp, 2);
    mr5b = _mm256_add_epi16 (mr5b,tmp);



    _mm256_storeu_si256((__m256i *)r4, mr4a);
    _mm256_storeu_si256((__m256i *)r4+1, mr4b);
    _mm256_storeu_si256((__m256i *)r5, mr5a);
    _mm256_storeu_si256((__m256i *)r5+1, mr5b);


    __m256i_grade_school_mul_32(t6, buf, r4, r5, s);;

    /*
     * t6 = w1 + 4*w3
     *    = (t6 - w0 - 4*w2 - 16*w4)/2
     *    = (t6 - t0 - 4*r  - 16*t8)/2
     */

    mt6a = _mm256_loadu_si256((__m256i *)t6);
    mt6b = _mm256_loadu_si256((__m256i *)t6+1);
    mt6c = _mm256_loadu_si256((__m256i *)t6+2);
    mt6d = _mm256_loadu_si256((__m256i *)t6+3);

    /*  t6 = t6 - t0 */
    mt6a = _mm256_sub_epi16(mt6a, mt0a);
    mt6b = _mm256_sub_epi16(mt6b, mt0b);
    mt6c = _mm256_sub_epi16(mt6c, mt0c);
    mt6d = _mm256_sub_epi16(mt6d, mt0d);

    /*  t6 = t6 - 4r0 */
    mt6a = _mm256_sub_epi16(mt6a, _mm256_slli_epi16(mr0a,2));
    mt6b = _mm256_sub_epi16(mt6b, _mm256_slli_epi16(mr0b,2));
    mt6c = _mm256_sub_epi16(mt6c, _mm256_slli_epi16(mr1a,2));
    mt6d = _mm256_sub_epi16(mt6d, _mm256_slli_epi16(mr1b,2));

    /*  t6 = t6 - 16t8 */
    mt6a = _mm256_sub_epi16(mt6a, _mm256_slli_epi16(mt8a,4));
    mt6b = _mm256_sub_epi16(mt6b, _mm256_slli_epi16(mt8b,4));
    mt6c = _mm256_sub_epi16(mt6c, _mm256_slli_epi16(mt8c,4));
    mt6d = _mm256_sub_epi16(mt6d, _mm256_slli_epi16(mt8d,4));

    /*  t6 = t6/2 */
    mt6a = _mm256_srai_epi16(mt6a, 1);
    mt6b = _mm256_srai_epi16(mt6b, 1);
    mt6c = _mm256_srai_epi16(mt6c, 1);
    mt6d = _mm256_srai_epi16(mt6d, 1);

    /*
     * t2 = w3 = (t6 - r2)/3
     *    = (t6 - t4) * 43691
     * t4 = w1 = (4*r2-t6)/3
     *    = (4*r4-t6) * 43691
     */
    tmp = _mm256_set1_epi16 (43691);

    /* t2 = t6 - r2 */
    mt2a = _mm256_sub_epi16 (mt6a, mr2a);
    mt2b = _mm256_sub_epi16 (mt6b, mr2b);
    mt2c = _mm256_sub_epi16 (mt6c, mr3a);
    mt2d = _mm256_sub_epi16 (mt6d, mr3b);

    /* t2 = t2/3 */
    mt2a = _mm256_mullo_epi16 (mt2a, tmp);
    mt2b = _mm256_mullo_epi16 (mt2b, tmp);
    mt2c = _mm256_mullo_epi16 (mt2c, tmp);
    mt2d = _mm256_mullo_epi16 (mt2d, tmp);

    /* t4 = 4*r2 - t6 */
    mt4a = _mm256_sub_epi16(_mm256_slli_epi16(mr2a,2), mt6a);
    mt4b = _mm256_sub_epi16(_mm256_slli_epi16(mr2b,2), mt6b);
    mt4c = _mm256_sub_epi16(_mm256_slli_epi16(mr3a,2), mt6c);
    mt4d = _mm256_sub_epi16(_mm256_slli_epi16(mr3b,2), mt6d);

    /* t4 = t4/3 */
    mt4a = _mm256_mullo_epi16 (mt4a, tmp);
    mt4b = _mm256_mullo_epi16 (mt4b, tmp);
    mt4c = _mm256_mullo_epi16 (mt4c, tmp);
    mt4d = _mm256_mullo_epi16 (mt4d, tmp);

    /*
     * now we have
     *  t0 = w0
     *  t4 = w1
     *  r  = w2
     *  t2 = w3
     *  t8 = w4
     * putting them back
     */

    mt0c = _mm256_add_epi16 (mt0c, mt4a);
    mt0d = _mm256_add_epi16 (mt0d, mt4b);
    mr0a = _mm256_add_epi16 (mr0a, mt4c);
    mr0b = _mm256_add_epi16 (mr0b, mt4d);
    mr1a = _mm256_add_epi16 (mr1a, mt2a);
    mr1b = _mm256_add_epi16 (mr1b, mt2b);
    mt8a = _mm256_add_epi16 (mt8a, mt2c);
    mt8b = _mm256_add_epi16 (mt8b, mt2d);

    _mm256_storeu_si256((__m256i *)r, mt0a);
    _mm256_storeu_si256((__m256i *)r+1, mt0b);
    _mm256_storeu_si256((__m256i *)r1, mt0c);
    _mm256_storeu_si256((__m256i *)r1+1, mt0d);
    _mm256_storeu_si256((__m256i *)r2, mr0a);
    _mm256_storeu_si256((__m256i *)r2+1, mr0b);
    _mm256_storeu_si256((__m256i *)r3, mr1a);
    _mm256_storeu_si256((__m256i *)r3+1, mr1b);
    _mm256_storeu_si256((__m256i *)r4, mt8a);
    _mm256_storeu_si256((__m256i *)r4+1, mt8b);
    _mm256_storeu_si256((__m256i *)r5, mt8c);
    _mm256_storeu_si256((__m256i *)r5+1, mt8d);

    return 0;
}




int test_toom3()
{
    uint16_t N;     // dimension
    uint16_t *a;    // first polynomial
    uint16_t *b;    // second polynomial
    uint16_t *buf;    // buffer
    uint16_t *r;    // result
    uint16_t *r2;    // result
    uint16_t i,j;
    uint16_t test_dim;
    float ss0, ss1,ss2, ss3;
    clock_t start, end;

    N = 96;
    a = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    b = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    buf = (uint16_t*) malloc (4*N*sizeof(uint16_t));
    r = (uint16_t*) malloc (4*N*sizeof(uint64_t));
    r2 = (uint16_t*) malloc (4*N*sizeof(uint64_t));

    cout<<"testing toom 3"<<endl;
    for (test_dim=65;test_dim<96;test_dim++)
    {
        ss0 = 0;
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
            grade_school_mul(r,  a, b, test_dim);
            end = clock();
            ss0 += (float)(end-start);


            start = clock();
            toom3(r2, buf, a, b, test_dim);
            end = clock();
            ss1 += (float)(end-start);


            start = clock();
            toom3__mm256i_SB(r2, buf, a, b, test_dim);
            end = clock();
            ss2 += (float)(end-start);


            start = clock();
            __mm256i_toom3__mm256i_SB(r2, buf, a, b, test_dim);
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
        cout<<ss0<<" "<<ss1<<" "<<ss2<<" "<<ss3<<endl;
    }
    return 0;
}


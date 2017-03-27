#include "simd_poly.h"

void print256_num(__m256i var)
{
    uint16_t *val = (uint16_t*) &var;
    printf("%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
           val[0], val[1], val[2], val[3],
           val[4], val[5], val[6], val[7],
           val[8], val[9], val[10], val[11],
           val[12], val[13], val[14], val[15]);
}

void print32poly(uint16_t* poly)
{
    uint16_t i;
    for (i=0;i<8;i++)
        printf("%d, ", poly[i]);
    printf("\n");
    for (i=8;i<16;i++)
        printf("%d, ", poly[i]);
    printf("\n");
    for (i=16;i<24;i++)
        printf("%d, ", poly[i]);
    printf("\n");
    for (i=24;i<32;i++)
        printf("%d, ", poly[i]);
    printf("\n\n");
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
    uint16_t *r2;
    uint64_t *r64;    // result
    uint16_t i,j;
    uint16_t test_dim;
    __m256i   *a256;
    __m256i   *b256;
    __m256i   *r256;
    __m256i   *buf256;

    test_SB_32();

    clock_t start, end;
    double karat, scht, toomt,toomtrec;
    ofstream fout;
    float ss1,ss2;
    fout.open("result.txt");
    N = 256*4;
    a = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    b = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    buf = (uint16_t*) malloc (4*N*sizeof(uint16_t));
    r = (uint16_t*) malloc (4*N*sizeof(uint64_t));
    r2 = (uint16_t*) malloc (4*N*sizeof(uint64_t));
    r64 = (uint64_t*)r;
    fout<<test_dim<<" ";

    test_dim = 93;



    srand(3);


    test_toom3();




    a256 = (__m256i *)a;
    b256 = (__m256i *)b;
    r256 = (__m256i *)r;
    buf256 = (__m256i *)buf;
    for(i=0; i< test_dim;i++)
    {
        a[i] = rand()&0x07FF;
        b[i] = rand()&0x07FF;

    }
    printf("a: ");
    for(i=0; i< test_dim;i++)
    {
        printf("%d, ", a[i]);
    }
    printf("\n");

    printf("b: ");
    for(i=0; i< test_dim;i++)
    {
        printf("%d, ", b[i]);
    }
    printf("\n");

    __mm256i_toom3(r, buf,a, b, test_dim);

    printf("final:");
    print32poly(r);
    print32poly(r+32);
    print32poly(r+64);
    print32poly(r+96);
    print32poly(r+128);
    print32poly(r+160);

  /*
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



    for (test_dim=33;test_dim<34;test_dim++)
    {
        ss1 = 0;
        ss2 = 0;
        cout<<"dimension: "<<test_dim<<" ";
        for (j=0;j<1000;j++)
        {
            memset(a, 0, sizeof(uint16_t)*400);
            memset(b, 0, sizeof(uint16_t)*400);
            for(i=0; i< test_dim;i++)
            {
                a[i] = rand()&0x07FF;
                b[i] = rand()&0x07FF;
            }

//            __mm256_karatsuba_768(r, buf, a, b, test_dim);
            start = clock();
            karatsuba(r, buf, a, b, test_dim, 32);
            end = clock();
            ss1 += (float)(end-start);


            start = clock();
            karatsuba_32(r2, buf, a, b, test_dim, 32);
            end = clock();
            ss2 += (float)(end-start);

            for (i=0;i<test_dim*2-1;i++)
            {
                if (r[i]!=r2[i])
                {
                    printf("error\n");

                for (j=0;j<test_dim*2-1;j++)
                            {
                    printf("%d %d %d %d\n", j, r[j],r2[j], r[j]-r2[j]);
                            }
                    return 0;}
            }
        }
        cout<<ss1<<" "<<ss2<<" "<<endl;
    }
*/

    return 0;

}


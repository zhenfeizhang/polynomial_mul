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
void print256_num_mod_q(__m256i var)
{
    uint16_t *val = (uint16_t*) &var;
    printf("%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
            val[0]%(2048), val[1]%(2048), val[2]%(2048), val[3]%(2048),
            val[4]%(2048), val[5]%(2048), val[6]%(2048), val[7]%(2048),
            val[8]%(2048), val[9]%(2048), val[10]%(2048), val[11]%(2048),
            val[12]%(2048), val[13]%(2048), val[14]%(2048), val[15]%(2048));
}
void print32poly(uint16_t* poly)
{
    uint16_t i;
    for (i=0;i<8;i++)
        printf("%d, ", poly[i]%2048);
    printf("\n");
    for (i=8;i<16;i++)
        printf("%d, ", poly[i]%2048);
    printf("\n");
    for (i=16;i<24;i++)
        printf("%d, ", poly[i]%2048);
    printf("\n");
    for (i=24;i<32;i++)
        printf("%d, ", poly[i]%2048);
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

    test_dim = 380;



    srand(3);


    for(i=0; i< test_dim;i++)
    {
        a[i] = rand()&0x07FF;
        b[i] = rand()&0x07FF;

    }


    toom4__mm256i_toom3(r, buf, a, b, test_dim);
    __mm256i_toom4__mm256i_toom3(r, buf, a, b, test_dim);
    printf("\n");
    grade_school_mul(r2, a, b, test_dim);
/*    for(i=0; i< test_dim*2;i++)
    {
        printf("%d %d, %d, %d \n", i, r[i]%2048, r2[i]%2048, (r[i]-r2[i])%2048);
    }
    printf("\n");
*/
    test_toom4();
//    test_toom3();
//    test_karatsuba();


    return 0;

}


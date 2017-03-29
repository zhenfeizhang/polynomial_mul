#include "simd_poly.h"
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

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

void test_ntru()
{
    uint16_t N;     // dimension
    uint16_t *a;    // first polynomial
    uint16_t *b;    // second polynomial
    uint16_t *buf;    // buffer
    uint16_t *r;    // result
    uint16_t i,j;
    float ss0, ss1,ss2, ss3;
    clock_t start, end;
    uint64_t startc, endc;
    N = 768;
    a = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    b = (uint16_t*) malloc (2*N*sizeof(uint16_t));
    buf = (uint16_t*) malloc (4*N*sizeof(uint16_t));
    r = (uint16_t*) malloc (4*N*sizeof(uint64_t));

    ss2 = 0;
    ss3 = 0;
    for (j=0;j<1000;j++)
    {
        memset(a, 0, 2*N*sizeof(uint16_t));
        memset(b, 0, 2*N*sizeof(uint16_t));
        memset(buf, 0, 4*N*sizeof(uint16_t));
        memset(r, 0,4*N*sizeof(uint16_t));
        for(i=0; i< N;i++)
        {
            a[i] = rand()&0x07FF;
            b[i] = rand()&0x07FF;
        }

        start = clock();startc = rdtsc();
        __mm256i_karatsuba__mm256_toom4(r, buf, a, b, N);
        endc = rdtsc();
        end = clock();
        ss3 += (float)(end-start);
        ss0 += (endc-startc);
        cout<<(float)(end-start)<<" "<< (endc-startc)<<" ";


        start = clock();startc = rdtsc();
        karatsuba_old(r, buf, a, b, N);
        endc = rdtsc();
        end = clock();
        ss2 += (float)(end-start);
        cout<<(float)(end-start)<<" "<<(endc-startc)<<endl;
        ss1 +=(endc-startc);
    }
    cout<<endl;
    cout<<ss3<<" "<<ss0<<" "<<ss2<<" "<<ss1<<endl;
}


int main()
{

    srand(3);
/*
    test_SB_32();
    test_toom3();
    test_toom4();
    test_karatsuba();
*/
    test_ntru();

    return 0;

}


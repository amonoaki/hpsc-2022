#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

void print(__m256 vec)
{
	float temp[8];
	_mm256_store_ps(temp,vec);
	for(int i=0;i<8;i++)
	{
		printf("%g ",temp[i]);
	}
	printf("\n");
}

int main() {
  const int N = 8;
  float temp[N],x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }

  for(int i=0; i<N; i++) {
    
    //printf("loop : %d \n",i);
    //init
    __m256 mj = _mm256_load_ps(m);
    __m256 xj = _mm256_load_ps(x);
    __m256 yj = _mm256_load_ps(y);
    
    //print(mj);
    //print(xj);    
        
    // calculate rx
    // float rx = x[i] - x[j];
    __m256 xi = _mm256_set1_ps(x[i]);
    __m256 rx = _mm256_sub_ps(xj, xi);

    
    // calc ry
    // float ry = y[i] - y[j];
    __m256 yi = _mm256_set1_ps(y[i]);
    __m256 ry = _mm256_sub_ps(yj, yi);

    //calculate r
    //float r = std::sqrt(rx * rx + ry * ry);
    __m256 squared_r = _mm256_add_ps(_mm256_mul_ps(rx,rx),_mm256_mul_ps(ry,ry));
    
    //calculate 1/r
    __m256 inverse_r = _mm256_rsqrt_ps(squared_r);
    //print(inverse_r);   
   
   //process i!=j with mask
    __m256 zero = _mm256_setzero_ps();
    __m256 mask = _mm256_cmp_ps(squared_r, zero, _CMP_GT_OQ);
    
    inverse_r = _mm256_blendv_ps(zero, inverse_r, mask);
    
    //calculate 1/r^3
    __m256 inverse_r3 = _mm256_mul_ps(_mm256_mul_ps(inverse_r, inverse_r), inverse_r);
    //o_mm256_store_ps(fx,inverse_r3);
    //printf("inverse_r:%f",fx[i]); 
    //  fx[i] -= rx * m[j] / (r * r * r);
    //  fy[i] -= ry * m[j] / (r * r * r);
    __m256 fxj = _mm256_sub_ps(fxj,_mm256_mul_ps(_mm256_mul_ps(rx, mj), inverse_r3));
    __m256 fyj = _mm256_sub_ps(fyj,_mm256_mul_ps(_mm256_mul_ps(ry, mj), inverse_r3));
    _mm256_store_ps(fx, fxj);
    _mm256_store_ps(fy, fyj);
    

  }
 for(int i =0;i<8;i++)	
    printf("%d %g %g\n",i,fx[i],fy[i]);
}

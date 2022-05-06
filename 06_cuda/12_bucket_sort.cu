#include <cstdio>
#include <cstdlib>

__global__ void init_bucket(int *bucket) {
  bucket[threadIdx.x] = 0;
}

__global__ void load_bucket(int *key,int *bucket) {
  int i =threadIdx.x;
  atomicAdd(&bucket[key[i]], 1);
}

__global__ void load_key(int *key,int *bucket) {
  int pos = 0;
  int i = threadIdx.x;
  for (int j=0; j<i; j++) {
    pos += bucket[j];
  }
  for (; bucket[i]>0; bucket[i]--) {
      key[pos++] = i;
  }  
}

int main() {
  int n = 50;
  int range = 5;
  int *key;
  int *bucket;
  cudaMallocManaged(&key, n*sizeof(int));
  cudaMallocManaged(&bucket, range*sizeof(int));
  
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");
  
  init_bucket<<<1,range>>>(bucket);
  load_bucket<<<1,n>>>(key,bucket);
  load_key<<<1,range>>>(key,bucket);
  
  cudaDeviceSynchronize();
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  cudaFree(key);
  cudaFree(bucket);

}


#include <stdio.h>

__global__ void saxpy(int n, float a, float *x, float *y, char *ad, char *bd)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n){ y[i] = a*x[i] + y[i];
              ad[0] = 'C';
            }
}

int main(void)
{
  int N = 1<<20;
  float *x, *y, *d_x, *d_y;
  x = (float*)malloc(N*sizeof(float));
  y = (float*)malloc(N*sizeof(float));

  char *a, *b, *a_d, *b_d;
  a = (char*)malloc(sizeof(char));
  b = (char*)malloc(sizeof(char));

  cudaMalloc(&d_x, N*sizeof(float)); 
  cudaMalloc(&d_y, N*sizeof(float));

  cudaMalloc(&a_d, sizeof(char));
  cudaMalloc(&b_d, sizeof(char));
  printf("Size of x:%d, d_x:%d\n",sizeof(x), sizeof(d_x));

  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

  a[0] = 'A';
  b[0] = 'B';
  cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, N*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(a_d, a, sizeof(char), cudaMemcpyHostToDevice);
  cudaMemcpy(b_d, b, sizeof(char), cudaMemcpyHostToDevice);

  // Perform SAXPY on 1M elements
  saxpy<<<(N+255)/256, 256>>>(N, 2.0f, d_x, d_y, a_d, b_d);

  cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(a, a_d, sizeof(char), cudaMemcpyDeviceToHost);

  float maxError = 0.0f;
  for (int i = 0; i < N; i++)
    maxError = max(maxError, abs(y[i]-4.0f));
  printf("Max error: %f\n", maxError);

  printf("char a:%c\n", a[0]);

  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(a_d);
  cudaFree(b_d);
  free(x);
  free(y);
  free(a);
  free(b);
}

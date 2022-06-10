#include <math.h>
#include <cstdio>
#include <cstdlib>

const int nit = 50;
const int N = 41;
const float dx = 2.0 / (N - 1);
const float dy = 2.0 / (N - 1);
const float rho = 1;
const float nu = .1;
const float F = 1.0;
const float dt = .01;


__global__ void build_up_b(float* b,float* u, float* v)
{
	
		int j = threadIdx.x + 1;
		for (int i = 1; i < N - 1; i++) {
			b[j*N+i] = rho * (1.0 / dt * ((u[j*N+i+1] - u[j*N+i-1]) / (2.0 * dx) + (v[(j+1)*N+i] - v[(j-1)*N+i]) / (2.0 * dy))
				- pow((u[j*N+i+1] - u[j*N+i-1]) / (2.0 * dx), 2)
				- 2.0 * ((u[(j+1)*N+i] - u[(j-1)*N+i]) / (2.0 * dy) * (v[j*N+i+1] - v[j*N+i-1]) / (2.0 * dx))
				- pow((v[(j+1)*N+i] - v[(j-1)*N+i]) / (2.0 * dy), 2));
		}
	
	
		b[j*N+N-1] = rho * (1.0 / dt * ((u[j*N] - u[j*N+N-2]) / (2.0 * dx) + (v[(j+1)*N+N-1] - v[(j-1)*N+N-1]) / (2.0 * dy))
			- pow((u[j*N] - u[j*N+N-2]) / (2.0 * dx), 2)
			- 2.0 * ((u[(j+1)*N+N-1] - u[(j-1)*N+N-1]) / (2.0 * dy) * (v[j*N] - v[j*N+N-2]) / (2.0 * dx))
			- pow((v[(j+1)*N+N-1] - v[(j-1)*N+N-1]) / (2.0 * dy), 2));
	
	
		b[j*N] = rho * (1.0 / dt * ((u[j*N+1] - u[j*N+N-1]) / (2.0 * dx) + (v[(j+1)*N] - v[(j-1)*N]) / (2.0 * dy))
			- pow((u[j*N+1] - u[j*N+N-1]) / (2.0 * dx), 2)
			- 2.0 * ((u[(j+1)*N] - u[(j-1)*N]) / (2.0 * dy) * (v[j*N+1] - v[j*N+N-1]) / (2.0 * dx))
			- pow((v[(j+1)*N] - v[(j-1)*N]) / (2.0 * dy), 2));
	
}

__global__ void pressure_poisson_periodic(float* p, float* b,float*pn)
{
	int j = threadIdx.x + 1;
		
		
			for (int i = 1; i < N - 1; i++) {
				p[j*N+i] = (((pn[j*N+i+1] - pn[j*N+i-1]) * pow(dy, 2) +
					(pn[(j+1)*N+N-1] - pn[(j-1)*N+N-1]) * pow(dx, 2)) /
					(2.0 * (pow(dx, 2) + pow(dy, 2))) -
					pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j*N+i]);

			}
			p[j*N+N-1] = (((pn[j*N] - pn[j*N+N-2]) * pow(dy, 2) +
				(pn[(j+1)*N+N-1] - pn[(j-1)*N+N-1]) * pow(dx, 2)) /
				(2.0 * (pow(dx, 2) + pow(dy, 2))) -
				pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j*N+N-1]);

			p[j * N] = (((pn[j * N + 1] - pn[j * N + N - 1]) * pow(dy, 2) +
				(pn[(j + 1) * N] - pn[(j - 1) * N]) * pow(dx, 2)) /
				(2.0 * (pow(dx, 2) + pow(dy, 2))) -
				pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j * N]);
			for (int i = 0; i < N; i++)
			{
				p[(N-1)*N+i] = p[(N - 2) * N + i];
				p[i] = p[N+i];
			}
	
}

__global__ void init_vector(float* vec) {
	int j = threadIdx.x ;
	for (int i = 0; i < N; i++)
	{
		vec[j*N+i] = 0;
	}
}

__global__ void init_vector_one(float* vec) {
	int j = threadIdx.x;
	for (int i = 0; i < N; i++)
	{
		vec[j * N + i] = 1.0;
	}
}

__global__ void build_up_v_u(float* v, float* vn, float* u, float* un, float* p)
{
	int j = threadIdx.x;
	for (int i = 1; i < N - 1; i++) {
		u[j * N + i] = (un[j * N + i] -
			un[j * N + i] * dt / dx *
			(un[j * N + i] - un[j * N + i - 1]) -
			vn[j * N + i] * dt / dy *
			(un[j * N + i] - un[(j - 1) * N + i]) -
			dt / (2 * rho * dx) *
			(p[j * N + i + 1] - p[j * N + i - 1]) +
			nu * (dt / pow(dx, 2) *
				(un[j * N + i + 1] - 2 * un[j * N + i] + un[j * N + i - 1]) +
				dt / pow(dy, 2) *
				(un[(j + 1) * N + i] - 2 * un[j * N + i] + un[(j - 1) * N + i])) +
			F * dt);

		v[j * N + i] = (vn[j * N + i] -
			un[j * N + i] * dt / dx *
			(vn[j * N + i] - vn[j * N + i - 1]) -
			vn[j * N + i] * dt / dy *
			(vn[j * N + i] - vn[(j - 1) * N + i]) -
			dt / (2 * rho * dy) *
			(p[(j + 1) * N + i] - p[(j - 1) * N + i]) +
			nu * (dt / pow(dx, 2) *
				(vn[j * N + i + 1] - 2 * vn[j * N + i] + vn[j * N + i - 1]) +
				dt / pow(dy, 2) *
				(vn[(j + 1) * N + i] - 2 * vn[j * N + i] + vn[(j - 1) * N + i])));
	}

	// Periodic BC u @ x = 2
	u[j * N + N - 1] = (un[j * N + N - 1] - un[j * N + N - 1] * dt / dx *
		(un[j * N + N - 1] - un[j * N + N - 2]) -
		vn[j * N + N - 1] * dt / dy *
		(un[j * N + N - 1] - un[(j - 1) * N + N - 1]) -
		dt / (2 * rho * dx) *
		(p[j * N] - p[j * N + N - 2]) +
		nu * (dt / pow(dx, 2) *
			(un[j * N] - 2 * un[j * N + N - 1] + un[j * N + N - 2]) +
			dt / pow(dy, 2) *
			(un[(j + 1) * N + N - 1] - 2 * un[j * N + N - 1] + un[(j - 1) * N + N - 1])) + F * dt);

	// Periodic BC u @ x = 0
	u[j * N] = (un[j * N] - un[j * N] * dt / dx *
		(un[j * N] - un[j * N + N - 1]) -
		vn[j * N] * dt / dy *
		(un[j * N] - un[(j - 1) * N]) -
		dt / (2 * rho * dx) *
		(p[j * N + 1] - p[j * N + N - 1]) +
		nu * (dt / pow(dx, 2) *
			(un[j * N + 1] - 2 * un[j * N] + un[j * N + N - 1]) +
			dt / pow(dy, 2) *
			(un[(j + 1) * N] - 2 * un[j * N] + un[(j - 1) * N])) + F * dt);

	// Periodic BC v @ x = 2
	v[j * N + N - 1] = (vn[j * N + N - 1] - un[j * N + N - 1] * dt / dx *
		(vn[j * N + N - 1] - vn[j * N + N - 2]) -
		vn[j * N + N - 1] * dt / dy *
		(vn[j * N + N - 1] - vn[(j - 1) * N + N - 1]) -
		dt / (2 * rho * dy) *
		(p[(j + 1) * N + N - 1] - p[(j - 1) * N + N - 1]) +
		nu * (dt / pow(dx, 2) *
			(vn[j * N] - 2 * vn[j * N + N - 1] + vn[j * N + N - 2]) +
			dt / pow(dy, 2) *
			(vn[(j + 1) * N + N - 1] - 2 * vn[j * N + N - 1] + vn[(j - 1) * N + N - 1])));

	// Periodic BC v @ x = 0
	v[j * N] = (vn[j * N] - un[j * N] * dt / dx *
		(vn[j * N] - vn[j * N + N - 1]) -
		vn[j * N] * dt / dy *
		(vn[j * N] - vn[(j - 1) * N]) -
		dt / (2 * rho * dy) *
		(p[(j + 1) * N] - p[(j - 1) * N]) +
		nu * (dt / pow(dx, 2) *
			(vn[j * N + 1] - 2 * vn[j * N] + vn[j * N + N - 1]) +
			dt / pow(dy, 2) *
			(vn[(j + 1) * N] - 2 * vn[j * N] + vn[(j - 1) * N])));

}

__global__ void handle_boundary_conditions(float* u, float* v)
{	
	
	int i = threadIdx.x;
	u[i] = 0;
	u[(N - 1) * N + i] = 0;
	v[i] = 0;
	v[(N - 1) * N + i] = 0;
}

__global__ void sum(float &sum,float* vec)
{
	int j = threadIdx.x;
	for (int i = 0; i < N; i++)
	{
		sum += vec[j * N + i];
	}
}

int main()
{
	//cout << "start" << endl;
	clock_t start,finish;
	float udiff = 1.0;
	int	stepcount = 0;
	float* u,*un;
	float* v, *vn;
	float* p, *pn;
	float* b;
	float sumu, sumun;

	cudaMallocManaged(&v, N*N * sizeof(float));
	cudaMallocManaged(&vn, N * N * sizeof(float));
	cudaMallocManaged(&u, N * N * sizeof(float));
	cudaMallocManaged(&un, N * N * sizeof(float));
	cudaMallocManaged(&p, N * N * sizeof(float));
	cudaMallocManaged(&pn, N * N * sizeof(float));
	cudaMallocManaged(&b, N * N * sizeof(float));
	init_vector <<< 1, N >>> (v);
	init_vector << < 1, N >> > (vn);
	init_vector << < 1, N >> > (u);
	init_vector << < 1, N >> > (un);
	init_vector << < 1, N >> > (b);
	init_vector_one<<<1,N>>>(p);
	init_vector_one<<<1,N>>>(pn);
	start = clock();
//	auto tic = chrono::steady_clock::now();
	while (udiff > 0.001)
	{
		//cout << "step:" << stepcount << endl;
		for(int i =0;i<N*N;i++)
		{
			un[i] = u[i];
			vn[i] = v[i];
		}	

		build_up_b<<<1,N-2>>>(b, u, v);
		for (int q = 0; q < nit; q++)
		{
			 for(int i =0;i<N*N;i++)
                {
                        pn[i] = p[i];
                }
			pressure_poisson_periodic <<<1, N - 2>>> (p, b,pn);
		}
		build_up_v_u <<< 1, N - 2 >>> (v, vn, u, un, p);
		//cout << "create b and u" << endl;

		handle_boundary_conditions << <1, N >> > (u, v);
		//cout << "4" << endl;
		cudaDeviceSynchronize();
		sumu = 0.0;
		sumun = 0.0;
		for(int i=0;i<N*N;i++)
		{
			sumu+=u[i];
			sumun+=un[i];	
		}	
		//cudaDeviceSynchronize();
		udiff = (sumu - sumun) / sumu;
		//printf("sumu = %lf, sumun = %lf \n",sumu,sumun);
		//printf("udiff = %lf\n", udiff);	
		stepcount++;
	}
//	auto toc = chrono::steady_clock::now();
//	float time = chrono::duration<float>(toc - tic).count();
	finish = clock();
	float time = (float)(finish-start)/CLOCKS_PER_SEC;
	printf("udiff = %lf\n", udiff);
	printf("step = %d\n", stepcount);
	printf("%1.3lf sec\n", time);
	cudaFree(u);
	cudaFree(v);
	cudaFree(p);
	cudaFree(b);
	cudaFree(un);
	cudaFree(vn);
	cudaFree(pn);
}



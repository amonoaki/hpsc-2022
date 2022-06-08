//#include <cstdlib>
#include<iostream>
#include <cstdio>
#include <vector>
#include<cstring>
#include <chrono>
#include <math.h>
#include<openacc.h>


using namespace std;
typedef vector<vector<double> > matrix;

const int N = 41;
int nt = 10;
int nit = 50;
int c = 1;
double dx = 2.0 / (N - 1);
double dy = 2.0 / (N - 1);
double rho = 1;
double nu = .1;
double F = 1.0;
double dt = .01;


vector<double> linespace(int start, int end, int points)
{
	double step = (end - start) / (double)points;
	vector<double> ans;
	while (start <= end)
	{
		ans.push_back(start);
		start += step;
	}
	return ans;
}

matrix zeros(int x, int y)
{
	return vector<vector<double> >(x, vector<double>(y, 0.0));
}

matrix ones(int x, int y)
{
	return vector<vector<double> >(x, vector<double>(y, 1.0));
}

double* build_up_b(double u[], double v[])
{
	double b[N * N] = { 0 };
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			b[j*N+i] = rho * (1.0 / dt * ((u[j*N+i+1] - u[j*N+i-1]) / (2.0 * dx) + (v[(j+1)*N+i] - v[(j-1)*N+i]) / (2.0 * dy))
				- pow((u[j*N+i+1] - u[j*N+i-1]) / (2.0 * dx), 2)
				- 2.0 * ((u[(j+1)*N+i] - u[(j-1)*N+i]) / (2.0 * dy) * (v[j*N+i+1] - v[j*N+i-1]) / (2.0 * dx))
				- pow((v[(j+1)*N+i] - v[(j-1)*N+i]) / (2.0 * dy), 2));
		}
	}
	for (int j = 1; j < N - 1; j++)
	{
		b[j*N+N-1] = rho * (1.0 / dt * ((u[j*N] - u[j*N+N-2]) / (2.0 * dx) + (v[(j+1)*N+N-1] - v[(j-1)*N+N-1]) / (2.0 * dy))
			- pow((u[j*N] - u[j*N+N-2]) / (2.0 * dx), 2)
			- 2.0 * ((u[(j+1)*N+N-1] - u[(j-1)*N+N-1]) / (2.0 * dy) * (v[j*N] - v[j*N+N-2]) / (2.0 * dx))
			- pow((v[(j+1)*N+N-1] - v[(j-1)*N+N-1]) / (2.0 * dy), 2));
	}
	for (int j = 1; j < N - 1; j++)
	{
		b[j*N] = rho * (1.0 / dt * ((u[j*N+1] - u[j*N+N-1]) / (2.0 * dx) + (v[(j+1)*N] - v[(j-1)*N]) / (2.0 * dy))
			- pow((u[j*N+1] - u[j*N+N-1]) / (2.0 * dx), 2)
			- 2.0 * ((u[(j+1)*N] - u[(j-1)*N]) / (2.0 * dy) * (v[j*N+1] - v[j*N+N-1]) / (2.0 * dx))
			- pow((v[(j+1)*N] - v[(j-1)*N]) / (2.0 * dy), 2));
	}

	return b;
}

double* pressure_poisson_periodic(double p[], double b[])
{
	
	for (int q = 0; q < nit; q++)
	{
		double pn[N*N];
		memcpy(pn, p, sizeof(p));
		
		for (int j = 1; j < N - 1; j++) {
			for (int i = 1; i < N - 1; i++) {
				p[j*N+i] = (((pn[j*N+i+1] - pn[j*N+i-1]) * pow(dy, 2) +
					(pn[(j+1)*N+N-1] - pn[(j-1)*N+N-1]) * pow(dx, 2)) /
					(2.0 * (pow(dx, 2) + pow(dy, 2))) -
					pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j*N+i]);

			}
		}
		//Periodic BC Pressure @ x = 2
		for (int j = 1; j < N - 1; j++)
		{
			p[j*N+N-1] = (((pn[j*N] - pn[j*N+N-2]) * pow(dy, 2) +
				(pn[(j+1)*N+N-1] - pn[(j-1)*N+N-1]) * pow(dx, 2)) /
				(2.0 * (pow(dx, 2) + pow(dy, 2))) -
				pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j*N+N-1]);
		}

		//Periodic BC Pressure @ x = 0
		for (int j = 1; j < N - 1; j++)
		{
			p[j*N] = (((pn[j*N+1] - pn[j*N+N-1]) * pow(dy, 2) +
				(pn[(j+1)*N] - pn[(j-1)*N]) * pow(dx, 2)) /
				(2.0 * (pow(dx, 2) + pow(dy, 2))) -
				pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j*N]);
		}
		//Wall boundary conditions, pressure
		for (int i = 0; i < N; i++)
		{
			p[(N-1)*N+i] = p[(N - 2) * N + i];
			p[i] = p[N+i];
		}
	}
	return p;
}

double sum(double m[])
{
	double ans = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			ans += m[j*N+i];
		}
	}
	return ans;
}

int main()
{
	cout << "OPENACC" << endl;
	double udiff = 1.0;
	int	stepcount = 0;
	double u[N*N] = { 0 };
	double un[N * N] = { 0 };
	double v[N * N] = { 0 };
	double vn[N * N] = { 0 };
	double b[N * N] = { 0 };
	double p[N * N];
	double pn[N * N];
	std::fill(std::begin(p), std::end(p), 1.0);
	std::fill(std::begin(pn), std::end(pn), 1.0);
	auto tic = chrono::steady_clock::now();
#pragma acc parallel
	while (udiff>0.001)
	{
		/*cout << "step:" << stepcount << endl;*/
		memcpy(un, u, sizeof(u));
		memcpy(vn, v, sizeof(u));

		double b[N * N] = { 0 };
		{
#pragma acc loop
			for (int j = 1; j < N - 1; j++) {
				for (int i = 1; i < N - 1; i++) {
					b[j * N + i] = rho * (1.0 / dt * ((u[j * N + i + 1] - u[j * N + i - 1]) / (2.0 * dx) + (v[(j + 1) * N + i] - v[(j - 1) * N + i]) / (2.0 * dy))
						- pow((u[j * N + i + 1] - u[j * N + i - 1]) / (2.0 * dx), 2)
						- 2.0 * ((u[(j + 1) * N + i] - u[(j - 1) * N + i]) / (2.0 * dy) * (v[j * N + i + 1] - v[j * N + i - 1]) / (2.0 * dx))
						- pow((v[(j + 1) * N + i] - v[(j - 1) * N + i]) / (2.0 * dy), 2));
				}
			}
#pragma acc loop
			for (int j = 1; j < N - 1; j++)
			{
				b[j * N + N - 1] = rho * (1.0 / dt * ((u[j * N] - u[j * N + N - 2]) / (2.0 * dx) + (v[(j + 1) * N + N - 1] - v[(j - 1) * N + N - 1]) / (2.0 * dy))
					- pow((u[j * N] - u[j * N + N - 2]) / (2.0 * dx), 2)
					- 2.0 * ((u[(j + 1) * N + N - 1] - u[(j - 1) * N + N - 1]) / (2.0 * dy) * (v[j * N] - v[j * N + N - 2]) / (2.0 * dx))
					- pow((v[(j + 1) * N + N - 1] - v[(j - 1) * N + N - 1]) / (2.0 * dy), 2));
			}
#pragma acc loop
			for (int j = 1; j < N - 1; j++)
			{
				b[j * N] = rho * (1.0 / dt * ((u[j * N + 1] - u[j * N + N - 1]) / (2.0 * dx) + (v[(j + 1) * N] - v[(j - 1) * N]) / (2.0 * dy))
					- pow((u[j * N + 1] - u[j * N + N - 1]) / (2.0 * dx), 2)
					- 2.0 * ((u[(j + 1) * N] - u[(j - 1) * N]) / (2.0 * dy) * (v[j * N + 1] - v[j * N + N - 1]) / (2.0 * dx))
					- pow((v[(j + 1) * N] - v[(j - 1) * N]) / (2.0 * dy), 2));
			}
		}
		{

			for (int q = 0; q < nit; q++)
			{
				double pn[N * N];
				memcpy(pn, p, sizeof(p));
#pragma acc loop
				for (int j = 1; j < N - 1; j++) {
					for (int i = 1; i < N - 1; i++) {
						p[j * N + i] = (((pn[j * N + i + 1] - pn[j * N + i - 1]) * pow(dy, 2) +
							(pn[(j + 1) * N + N - 1] - pn[(j - 1) * N + N - 1]) * pow(dx, 2)) /
							(2.0 * (pow(dx, 2) + pow(dy, 2))) -
							pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j * N + i]);

					}
				}
				//Periodic BC Pressure @ x = 2
#pragma acc loop
				for (int j = 1; j < N - 1; j++)
				{
					p[j * N + N - 1] = (((pn[j * N] - pn[j * N + N - 2]) * pow(dy, 2) +
						(pn[(j + 1) * N + N - 1] - pn[(j - 1) * N + N - 1]) * pow(dx, 2)) /
						(2.0 * (pow(dx, 2) + pow(dy, 2))) -
						pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j * N + N - 1]);
				}

				//Periodic BC Pressure @ x = 0
#pragma acc loop
				for (int j = 1; j < N - 1; j++)
				{
					p[j * N] = (((pn[j * N + 1] - pn[j * N + N - 1]) * pow(dy, 2) +
						(pn[(j + 1) * N] - pn[(j - 1) * N]) * pow(dx, 2)) /
						(2.0 * (pow(dx, 2) + pow(dy, 2))) -
						pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j * N]);
				}
				//Wall boundary conditions, pressure
#pragma acc loop
				for (int i = 0; i < N; i++)
				{
					p[(N - 1) * N + i] = p[(N - 2) * N + i];
					p[i] = p[N + i];
				}
			}
		}

		/*cout << "create b and u" << endl;*/
#pragma acc loop
		for (int j = 1; j < N - 1; j++) {
			for (int i = 1; i < N - 1; i++) {
				u[j*N+i] = (un[j*N+i] -
					un[j*N+i] * dt / dx *
					(un[j*N+i] - un[j*N+i-1]) -
					vn[j*N+i] * dt / dy *
					(un[j*N+i] - un[(j-1)*N+i]) -
					dt / (2 * rho * dx) *
					(p[j*N+i+1] - p[j*N+i-1]) +
					nu * (dt / pow(dx, 2) *
						(un[j*N+i+1] - 2 * un[j*N+i] + un[j*N+i-1]) +
						dt / pow(dy, 2) *
						(un[(j+1)*N+i] - 2 * un[j*N+i] + un[(j-1)*N+i])) +
					F * dt);
			}
		}
		/*cout << "1" << endl;*/
#pragma acc loop
		for (int j = 1; j < N - 1; j++)
		{
			for (int i = 1; i < N - 1; i++)
			{
				v[j*N+i] = (vn[j*N+i] -
					un[j*N+i] * dt / dx *
					(vn[j*N+i] - vn[j*N+i-1]) -
					vn[j*N+i] * dt / dy *
					(vn[j*N+i] - vn[(j-1)*N+i]) -
					dt / (2 * rho * dy) *
					(p[(j+1)*N+i] - p[(j-1)*N+i]) +
					nu * (dt / pow(dx, 2) *
						(vn[j*N+i+1] - 2 * vn[j*N+i] + vn[j*N+i-1]) +
						dt / pow(dy, 2) *
						(vn[(j+1)*N+i] - 2 * vn[j*N+i] + vn[(j-1)*N+i])));
			}

		}
		/*cout << "2" << endl;*/
#pragma acc loop
		for (int j = 1; j < N - 1; j++)
		{
			// Periodic BC u @ x = 2
			u[j*N+N-1] = (un[j*N+N-1] - un[j*N+N-1] * dt / dx *
				(un[j*N+N-1] - un[j*N+N-2]) -
				vn[j*N+N-1] * dt / dy *
				(un[j*N+N-1] - un[(j-1)*N+N-1]) -
				dt / (2 * rho * dx) *
				(p[j*N] - p[j*N+N-2]) +
				nu * (dt / pow(dx, 2) *
					(un[j*N] - 2 * un[j*N+N-1] + un[j*N+N-2]) +
					dt / pow(dy, 2) *
					(un[(j+1)*N+N-1] - 2 * un[j*N+N-1] + un[(j-1)*N+N-1])) + F * dt);
		}
#pragma acc loop
		for (int j = 1; j < N - 1; j++)
		{
			// Periodic BC u @ x = 0
			u[j*N] = (un[j*N] - un[j*N] * dt / dx *
				(un[j*N] - un[j*N+N-1]) -
				vn[j*N] * dt / dy *
				(un[j*N] - un[(j-1)*N]) -
				dt / (2 * rho * dx) *
				(p[j*N+1] - p[j*N+N-1]) +
				nu * (dt / pow(dx, 2) *
					(un[j*N+1] - 2 * un[j*N] + un[j*N+N-1]) +
					dt / pow(dy, 2) *
					(un[(j+1)*N] - 2 * un[j*N] + un[(j-1)*N])) + F * dt);
		}
#pragma acc loop
		for (int j = 1; j < N - 1; j++)
		{
			// Periodic BC v @ x = 2
			v[j*N+N-1] = (vn[j*N+N-1] - un[j*N+N-1] * dt / dx *
				(vn[j*N+N-1] - vn[j*N+N-2]) -
				vn[j*N+N-1] * dt / dy *
				(vn[j*N+N-1] - vn[(j-1)*N+N-1]) -
				dt / (2 * rho * dy) *
				(p[(j+1)*N+N-1] - p[(j-1)*N+N-1]) +
				nu * (dt / pow(dx, 2) *
					(vn[j*N] - 2 * vn[j*N+N-1] + vn[j*N+N-2]) +
					dt / pow(dy, 2) *
					(vn[(j+1)*N+N-1] - 2 * vn[j*N+N-1] + vn[(j-1)*N+N-1])));
		}
#pragma acc loop
		for (int j = 1; j < N - 1; j++)
		{
			// Periodic BC v @ x = 0
			v[j*N] = (vn[j*N] - un[j*N] * dt / dx *
				(vn[j*N] - vn[j*N+N-1]) -
				vn[j*N] * dt / dy *
				(vn[j*N] - vn[(j-1)*N]) -
				dt / (2 * rho * dy) *
				(p[(j+1)*N] - p[(j-1)*N]) +
				nu * (dt / pow(dx, 2) *
					(vn[j*N+1] - 2 * vn[j*N] + vn[j*N+N-1]) +
					dt / pow(dy, 2) *
					(vn[(j+1)*N] - 2 * vn[j*N] + vn[(j-1)*N])));
		}
#pragma acc loop
		for (int i = 0; i < N; i++)
		{
			u[i] = 0;
			u[(N-1)*N+i] = 0;
			v[i] = 0;
			v[(N - 1) * N + i] = 0;
		}
		//cout << "3" << endl;
		//cout << "udiff= " << udiff << endl;
		udiff = (sum(u) - sum(un)) / sum(u);
		stepcount++;
	}
	auto toc = chrono::steady_clock::now();
	double time = chrono::duration<double>(toc - tic).count();
	printf("udiff = %lf\n",udiff);
	printf("step = %d\n", stepcount);
	printf("%1.3lf sec\n", time);
}

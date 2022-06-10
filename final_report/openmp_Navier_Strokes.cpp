#include <cstdlib>
#include <cstdio>
#include <vector>
#include <chrono>
#include <math.h>
#include<iostream>
using namespace std;
typedef vector<vector<double> > matrix;

int nx = 41;
int ny = 41;
int nt = 10;
int nit = 50;
int c = 1;
double dx = 2.0 / (nx - 1);
double dy = 2.0 / (ny - 1);
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

//matrix build_up_b(matrix& u, matrix& v)
//{
//	matrix b = zeros(ny, nx);
//	#pragma omp for
//	for (int j = 1; j < ny - 1; j++) {
//		for (int i = 1; i < nx - 1; i++) {
//			b[j][i] = rho * (1.0 / dt * ((u[j][i + 1] - u[j][i - 1]) / (2.0 * dx) + (v[j + 1][i] - v[j - 1][i]) / (2.0 * dy))
//				- pow((u[j][i + 1] - u[j][i - 1]) / (2.0 * dx), 2)
//				- 2.0 * ((u[j + 1][i] - u[j - 1][i]) / (2.0 * dy) * (v[j][i + 1] - v[j][i - 1]) / (2.0 * dx))
//				- pow((v[j + 1][i] - v[j - 1][i]) / (2.0 * dy), 2));
//		}
//	}
//	#pragma omp for
//	for (int j = 1; j < ny - 1; j++)
//	{
//		b[j][nx - 1] = rho * (1.0 / dt * ((u[j][0] - u[j][nx - 2]) / (2.0 * dx) + (v[j + 1][nx - 1] - v[j - 1][nx - 1]) / (2.0 * dy))
//			- pow((u[j][0] - u[j][nx - 2]) / (2.0 * dx), 2)
//			- 2.0 * ((u[j + 1][nx - 1] - u[j - 1][nx - 1]) / (2.0 * dy) * (v[j][0] - v[j][nx - 2]) / (2.0 * dx))
//			- pow((v[j + 1][nx - 1] - v[j - 1][nx - 1]) / (2.0 * dy), 2));
//	}
//	#pragma omp for
//	for (int j = 1; j < ny - 1; j++)
//	{
//		b[j][0] = rho * (1.0 / dt * ((u[j][1] - u[j][nx - 1]) / (2.0 * dx) + (v[j + 1][0] - v[j - 1][0]) / (2.0 * dy))
//			- pow((u[j][1] - u[j][nx - 1]) / (2.0 * dx), 2)
//			- 2.0 * ((u[j + 1][0] - u[j - 1][0]) / (2.0 * dy) * (v[j][1] - v[j][nx - 1]) / (2.0 * dx))
//			- pow((v[j + 1][0] - v[j - 1][0]) / (2.0 * dy), 2));
//	}
//
//	return b;
//}
//
//matrix pressure_poisson_periodic(matrix& p, matrix& b)
//{
//	
//	for (int q = 0; q < nit; q++)
//	{
//		matrix pn(p);
//#pragma omp for
//		for (int j = 1; j < ny - 1; j++) {
//			for (int i = 1; i < nx - 1; i++) {
//				p[j][i] = (((pn[j][i + 1] - pn[j][i - 1]) * pow(dy, 2) +
//					(pn[j + 1][nx - 1] - pn[j - 1][nx - 1]) * pow(dx, 2)) /
//					(2.0 * (pow(dx, 2) + pow(dy, 2))) -
//					pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j][i]);
//
//			}
//		}
//		//Periodic BC Pressure @ x = 2
//#pragma omp for
//		for (int j = 1; j < ny - 1; j++)
//		{
//			p[j][nx - 1] = (((pn[j][0] - pn[j][nx - 2]) * pow(dy, 2) +
//				(pn[j + 1][nx - 1] - pn[j - 1][nx - 1]) * pow(dx, 2)) /
//				(2.0 * (pow(dx, 2) + pow(dy, 2))) -
//				pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j][nx - 1]);
//		}
//
//		//Periodic BC Pressure @ x = 0
//#pragma omp for
//		for (int j = 1; j < ny - 1; j++)
//		{
//			p[j][0] = (((pn[j][1] - pn[j][nx - 1]) * pow(dy, 2) +
//				(pn[j + 1][0] - pn[j - 1][0]) * pow(dx, 2)) /
//				(2.0 * (pow(dx, 2) + pow(dy, 2))) -
//				pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j][0]);
//		}
//		//Wall boundary conditions, pressure
//#pragma omp for
//		for (int i = 0; i < nx; i++)
//		{
//			p[ny - 1][i] = p[ny - 2][i];
//			p[0][i] = p[1][i];
//		}
//	}
//	return p;
//}

double sum(const matrix& m)
{
	double ans = 0;
	for (int i = 0; i < ny; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			ans += m[j][i];
		}
	}
	return ans;
}

int main()
{
	cout << "OPENMP" << endl;
	double udiff = 1.0;
	int	stepcount = 0;
	matrix u = zeros(ny, nx);
	matrix	un = zeros(ny, nx);

	matrix	v = zeros(ny, nx);
	matrix	vn = zeros(ny, nx);

	matrix	p(ny, vector<double>(nx, 1.0f));
	matrix	pn(ny, vector<double>(nx, 1.0f));

	matrix	b = zeros(ny, nx);
	auto tic = chrono::steady_clock::now();

	while (udiff > 0.001)
	{
		//cout << "step:" << stepcount << endl;
		un = u;
		vn = v;
		matrix b = zeros(ny, nx);
#pragma omp parallel for collapse(2)
		for (int j = 1; j < ny - 1; j++) {
			for (int i = 1; i < nx - 1; i++) {
				b[j][i] = rho * (1.0 / dt * ((u[j][i + 1] - u[j][i - 1]) / (2.0 * dx) + (v[j + 1][i] - v[j - 1][i]) / (2.0 * dy))
					- pow((u[j][i + 1] - u[j][i - 1]) / (2.0 * dx), 2)
					- 2.0 * ((u[j + 1][i] - u[j - 1][i]) / (2.0 * dy) * (v[j][i + 1] - v[j][i - 1]) / (2.0 * dx))
					- pow((v[j + 1][i] - v[j - 1][i]) / (2.0 * dy), 2));
			}
		}
#pragma omp parallel for
		for (int j = 1; j < ny - 1; j++)
		{
			b[j][nx - 1] = rho * (1.0 / dt * ((u[j][0] - u[j][nx - 2]) / (2.0 * dx) + (v[j + 1][nx - 1] - v[j - 1][nx - 1]) / (2.0 * dy))
				- pow((u[j][0] - u[j][nx - 2]) / (2.0 * dx), 2)
				- 2.0 * ((u[j + 1][nx - 1] - u[j - 1][nx - 1]) / (2.0 * dy) * (v[j][0] - v[j][nx - 2]) / (2.0 * dx))
				- pow((v[j + 1][nx - 1] - v[j - 1][nx - 1]) / (2.0 * dy), 2));
	
			b[j][0] = rho * (1.0 / dt * ((u[j][1] - u[j][nx - 1]) / (2.0 * dx) + (v[j + 1][0] - v[j - 1][0]) / (2.0 * dy))
				- pow((u[j][1] - u[j][nx - 1]) / (2.0 * dx), 2)
				- 2.0 * ((u[j + 1][0] - u[j - 1][0]) / (2.0 * dy) * (v[j][1] - v[j][nx - 1]) / (2.0 * dx))
				- pow((v[j + 1][0] - v[j - 1][0]) / (2.0 * dy), 2));
		}

		for (int q = 0; q < nit; q++)
		{
			matrix pn(p);
#pragma omp parallel for collapse(2)
			for (int j = 1; j < ny - 1; j++) {
				for (int i = 1; i < nx - 1; i++) {
					p[j][i] = (((pn[j][i + 1] - pn[j][i - 1]) * pow(dy, 2) +
						(pn[j + 1][nx - 1] - pn[j - 1][nx - 1]) * pow(dx, 2)) /
						(2.0 * (pow(dx, 2) + pow(dy, 2))) -
						pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j][i]);

				}
			}
			//Periodic BC Pressure @ x = 2
#pragma omp parallel for
			for (int j = 1; j < ny - 1; j++)
			{
				p[j][nx - 1] = (((pn[j][0] - pn[j][nx - 2]) * pow(dy, 2) +
					(pn[j + 1][nx - 1] - pn[j - 1][nx - 1]) * pow(dx, 2)) /
					(2.0 * (pow(dx, 2) + pow(dy, 2))) -
					pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j][nx - 1]);
			
				p[j][0] = (((pn[j][1] - pn[j][nx - 1]) * pow(dy, 2) +
					(pn[j + 1][0] - pn[j - 1][0]) * pow(dx, 2)) /
					(2.0 * (pow(dx, 2) + pow(dy, 2))) -
					pow(dx, 2) * pow(dy, 2) / (2.0 * (pow(dx, 2) + pow(dy, 2))) * b[j][0]);
			}
			//Wall boundary conditions, pressure
#pragma omp parallel for
			for (int i = 0; i < nx; i++)
			{
				p[ny - 1][i] = p[ny - 2][i];
				p[0][i] = p[1][i];
			}
		}
		//cout << "create b and u" << endl;
#pragma omp parallel for collapse(2)
		for (int j = 1; j < ny - 1; j++) {
			for (int i = 1; i < nx - 1; i++) {
				u[j][i] = (un[j][i] -
					un[j][i] * dt / dx *
					(un[j][i] - un[j][i - 1]) -
					vn[j][i] * dt / dy *
					(un[j][i] - un[j - 1][i]) -
					dt / (2 * rho * dx) *
					(p[j][i + 1] - p[j][i - 1]) +
					nu * (dt / pow(dx, 2) *
						(un[j][i + 1] - 2 * un[j][i] + un[j][i - 1]) +
						dt / pow(dy, 2) *
						(un[j + 1][i] - 2 * un[j][i] + un[j - 1][i])) +
					F * dt);
			
				v[j][i] = (vn[j][i] -
					un[j][i] * dt / dx *
					(vn[j][i] - vn[j][i - 1]) -
					vn[j][i] * dt / dy *
					(vn[j][i] - vn[j - 1][i]) -
					dt / (2 * rho * dy) *
					(p[j + 1][i] - p[j - 1][i]) +
					nu * (dt / pow(dx, 2) *
						(vn[j][i + 1] - 2 * vn[j][i] + vn[j][i - 1]) +
						dt / pow(dy, 2) *
						(vn[j + 1][i] - 2 * vn[j][i] + vn[j - 1][i])));
			}

		}
		//cout << "2" << endl;
#pragma omp parallel for
		for (int j = 1; j < ny - 1; j++)
		{
			// Periodic BC u @ x = 2
			u[j][nx - 1] = (un[j][nx - 1] - un[j][nx - 1] * dt / dx *
				(un[j][nx - 1] - un[j][nx - 2]) -
				vn[j][nx - 1] * dt / dy *
				(un[j][nx - 1] - un[j - 1][nx - 1]) -
				dt / (2 * rho * dx) *
				(p[j][0] - p[j][nx - 2]) +
				nu * (dt / pow(dx, 2) *
					(un[j][0] - 2 * un[j][nx - 1] + un[j][nx - 2]) +
					dt / pow(dy, 2) *
					(un[j + 1][nx - 1] - 2 * un[j][nx - 1] + un[j - 1][nx - 1])) + F * dt);

			// Periodic BC u @ x = 0
			u[j][0] = (un[j][0] - un[j][0] * dt / dx *
				(un[j][0] - un[j][nx - 1]) -
				vn[j][0] * dt / dy *
				(un[j][0] - un[j - 1][0]) -
				dt / (2 * rho * dx) *
				(p[j][1] - p[j][nx - 1]) +
				nu * (dt / pow(dx, 2) *
					(un[j][1] - 2 * un[j][0] + un[j][nx - 1]) +
					dt / pow(dy, 2) *
					(un[j + 1][0] - 2 * un[j][0] + un[j - 1][0])) + F * dt);

			// Periodic BC v @ x = 2
			v[j][nx - 1] = (vn[j][nx - 1] - un[j][nx - 1] * dt / dx *
				(vn[j][nx - 1] - vn[j][nx - 2]) -
				vn[j][nx - 1] * dt / dy *
				(vn[j][nx - 1] - vn[j - 1][nx - 1]) -
				dt / (2 * rho * dy) *
				(p[j + 1][nx - 1] - p[j - 1][nx - 1]) +
				nu * (dt / pow(dx, 2) *
					(vn[j][0] - 2 * vn[j][nx - 1] + vn[j][nx - 2]) +
					dt / pow(dy, 2) *
					(vn[j + 1][nx - 1] - 2 * vn[j][nx - 1] + vn[j - 1][nx - 1])));

			// Periodic BC v @ x = 0
			v[j][0] = (vn[j][0] - un[j][0] * dt / dx *
				(vn[j][0] - vn[j][nx - 1]) -
				vn[j][0] * dt / dy *
				(vn[j][0] - vn[j - 1][0]) -
				dt / (2 * rho * dy) *
				(p[j + 1][0] - p[j - 1][0]) +
				nu * (dt / pow(dx, 2) *
					(vn[j][1] - 2 * vn[j][0] + vn[j][nx - 1]) +
					dt / pow(dy, 2) *
					(vn[j + 1][0] - 2 * vn[j][0] + vn[j - 1][0])));
		}
#pragma omp parallel for
		for (int i = 0; i < nx; i++)
		{
			u[0][i] = 0;
			u[ny - 1][i] = 0;
			v[0][i] = 0;
			v[ny - 1][i] = 0;
		}
//		printf("udiff = %lf\n", udiff);
		udiff = (sum(u) - sum(un)) / sum(u);
		stepcount++;
	}
	auto toc = chrono::steady_clock::now();
	double time = chrono::duration<double>(toc - tic).count();
	printf("udiff = %lf\n", udiff);
	printf("step = %d\n", stepcount);
	printf("%1.3lf sec\n", time);
}

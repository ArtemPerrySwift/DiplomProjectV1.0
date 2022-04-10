#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>
#include "matrix.h"

using namespace std;

typedef double real;
double nev;
int itarations;

enum method { ORDIN_LOS, LOS_LU, LOS_LU_inv, LOS_LUsq, LOS_DD, LOS_DI, LOS_ID };

const int PREC = 14;
/*
struct matrix {
	int* ig; // массив указателей начала строк (столбцов) в массивах ggl и ggu
	int* jg; // массив номеров строк(столбцов) заданного треугольника матрицы
	real* ggl; // массив внедиагональных элементов нижнего треугольника матрицы
	real* ggu; // массив внедиагональных элементов верхнего треугольника матрицы
	real* di; // массив диагональных элементов матрицы
	int n; // размерность матрицы
};
*/
int read_vec(real* vec, int n, ifstream& in);
int read_vec(int* vec, int n, ifstream& in);
int decomp_mat_LU(matrix & A, matrix & LU);
int decomp_mat_LU_inv(matrix & A, matrix& LU);
int decomp_mat_LUsq(matrix& A, matrix& LU);
int decomp_mat_DI(matrix& A, matrix& LU, bool is_DD);
int count_LOS(matrix& A, real* x, real* r, real* z, real* p, real* f, real* buf_v, int maxiter, real eps);
int count_LOS(matrix& A, matrix& LU, real* x, real* r, real* z, real* p, real* f, real* buf_v, real* buf_v1, int maxiter, real eps, bool is_L_diag_1, bool is_U_diag_1);
int mult_mat_vec(matrix& A, real* vec, real* res);
int calc_Lx(matrix& LU, real* x, real* f, bool is_diag_1);
int calc_Ux(matrix& LU, real* x, real* f, bool is_diag_1);
real skal(real* vec1, real* vec2, int n);
int write(real* x, int n, ofstream& out);


int read_vec(real* vec, int n, ifstream& in)
{
	for (int i = 0; i < n; i++)
		in >> vec[i];
	return 0;
}

int read_vec(int* vec, int n, ifstream& in)
{
	for (int i = 0; i < n; i++)
		in >> vec[i];
	return 0;
}

int decomp_mat_LU(matrix& A, matrix& LU)
{
	int n = A.n;

	real* A_di = A.di;
	real* A_ggl = A.ggl;
	real* A_ggu = A.ggu;
	int* A_ig = A.ig;
	int* A_jg = A.jg;

	real* LU_di = LU.di;
	real* LU_ggl = LU.ggl;
	real* LU_ggu = LU.ggu;
	int* LU_ig = LU.ig;
	int* LU_jg = LU.jg;

	int i, j, ij, ik, kj;
	int i_beg, i_end, j_beg, j_end;
	for (i = 0; i < n; i++)
	{
		i_beg = A_ig[i];
		i_end = A_ig[i + 1];
		LU_di[i] = A_di[i];
		for (ij = i_beg; ij < i_end; ij++)
		{
			j = A_jg[ij];
			j_beg = A_ig[j];
			j_end = A_ig[j + 1];
			LU_ggl[ij] = A_ggl[ij];
			LU_ggu[ij] = A_ggu[ij];
			for (ik = i_beg, kj = j_beg; kj < j_end && ik < ij;)
			{
				if (A_jg[ik] == A_jg[kj])
				{
					LU_ggl[ij] -= LU_ggl[ik] * LU_ggu[kj];
					LU_ggu[ij] -= LU_ggl[kj] * LU_ggu[ik];
					ik++;
					kj++;
				}
				else if (A_jg[ik] > A_jg[kj]) kj++;
				else ik++;
			}
			LU_ggu[ij] /= LU_di[j];
			LU_di[i] -= LU_ggl[ij] * LU_ggu[ij];
		}
	}
	return 0;
}

int decomp_mat_LU_inv(matrix& A, matrix& LU)
{
	int n = A.n;

	real* A_di = A.di;
	real* A_ggl = A.ggl;
	real* A_ggu = A.ggu;
	int* A_ig = A.ig;
	int* A_jg = A.jg;

	real* LU_di = LU.di;
	real* LU_ggl = LU.ggl;
	real* LU_ggu = LU.ggu;
	int* LU_ig = LU.ig;
	int* LU_jg = LU.jg;

	int i, j, ij, ik, kj;
	int i_beg, i_end, j_beg, j_end;
	for (i = 0; i < n; i++)
	{
		i_beg = A_ig[i];
		i_end = A_ig[i + 1];
		LU_di[i] = A_di[i];
		for (ij = i_beg; ij < i_end; ij++)
		{
			j = A_jg[ij];
			j_beg = A_ig[j];
			j_end = A_ig[j + 1];
			LU_ggl[ij] = A_ggl[ij];
			LU_ggu[ij] = A_ggu[ij];
			for (ik = i_beg, kj = j_beg; kj < j_end && ik < ij;)
			{
				if (A_jg[ik] == A_jg[kj])
				{
					LU_ggl[ij] -= LU_ggl[ik] * LU_ggu[kj]; // ik kj
					LU_ggu[ij] -= LU_ggl[kj] * LU_ggu[ik];
					ik++;
					kj++;
				}
				else if (A_jg[ik] > A_jg[kj]) kj++;
				else ik++;
			}
			LU_ggl[ij] /= LU_di[j];
			LU_di[i] -= LU_ggl[ij] * LU_ggu[ij];
		}
	}
	return 0;
}

int decomp_mat_LUsq(matrix& A, matrix& LU)
{
	int n = A.n;

	real* A_di = A.di;
	real* A_ggl = A.ggl;
	real* A_ggu = A.ggu;
	int* A_ig = A.ig;
	int* A_jg = A.jg;

	real* LU_di = LU.di;
	real* LU_ggl = LU.ggl;
	real* LU_ggu = LU.ggu;
	int* LU_ig = LU.ig;
	int* LU_jg = LU.jg;

	int i, j, ij, ik, kj;
	int i_beg, i_end, j_beg, j_end;
	for (i = 0; i < n; i++)
	{
		i_beg = A_ig[i];
		i_end = A_ig[i + 1];
		LU_di[i] = A_di[i];
		for (ij = i_beg; ij < i_end; ij++)
		{
			j = A_jg[ij];
			j_beg = A_ig[j];
			j_end = A_ig[j + 1];
			LU_ggl[ij] = A_ggl[ij];
			LU_ggu[ij] = A_ggu[ij];
			for (ik = i_beg, kj = j_beg; kj < j_end && ik < ij;)
			{
				if (A_jg[ik] == A_jg[kj])
				{
					LU_ggl[ij] -= LU_ggl[ik] * LU_ggu[kj];
					LU_ggu[ij] -= LU_ggl[kj] * LU_ggu[ik];
					ik++;
					kj++;
				}
				else if (A_jg[ik] > A_jg[kj]) kj++;
				else ik++;
			}
			LU_ggl[ij] /= LU_di[j];
			LU_ggu[ij] /= LU_di[j];
			LU_di[i] -= LU_ggl[ij] * LU_ggu[ij];
		}
		LU_di[i] = sqrt(LU_di[i]);
	}
	return 0;
}


int count_LOS(matrix& A, real* x, real* r, real* z, real* p, real* f, real* buf_v, int maxiter, real eps)
{
	int n = A.n;
	int i, j, k;
	//cout << endl;
	real err;
	real norm_ff = sqrt(skal(f, f, n));
	bool fl = true;
	i = 0;
	while (fl)
	{
		//cout << endl << "Begin again" << endl;
		mult_mat_vec(A, x, r);
		for (j = 0; j < n; j++)
		{
			r[j] = f[j] - r[j];
			z[j] = r[j];
		}
		mult_mat_vec(A, z, p);

		real a_k, b_k;
		real skal_pp;

		err = sqrt(skal(r, r, n)) / norm_ff;

		fl = i < maxiter&& err > eps;
		for (k = 0; fl && k < n; i++, k++)
		{
			skal_pp = skal(p, p, n);
			a_k = skal(p, r, n) / skal_pp;
			for (j = 0; j < n; j++)
			{
				x[j] += a_k * z[j];
				r[j] -= a_k * p[j];
			}
			mult_mat_vec(A, r, buf_v);
			b_k = -skal(p, buf_v, n) / skal_pp;
			for (j = 0; j < n; j++)
			{
				z[j] = r[j] + b_k * z[j];
				p[j] = buf_v[j] + b_k * p[j];
			}
			err = sqrt(skal(r, r, n)) / norm_ff;
			fl = i < maxiter&& err > eps;
			//cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
		}
		itarations = i;
		nev = err;
	}
	cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
	if (i == maxiter) return 1;
	return 0;
}

int count_BSGSTAB(matrix& A, real* x, real* r, real* r1, real* z, real* p, real* f, real* buf, real* r0, int maxiter, real eps)
{
	int n = A.n;
	int i, j, k;
	//cout << endl;
	real err;
	real norm_ff = sqrt(skal(f, f, n));
	bool fl = true;
	i = 0;
	while (fl)
	{
		//cout << endl << "Begin again" << endl;
		mult_mat_vec(A, x, r);
		for (j = 0; j < n; j++)
		{
			r[j] = r0[j] = f[j] - r[j];
			z[j] = r[j];
		}

		real a_k, b_k, y_k;
		real skal_rr0, skal_rr1;

		err = sqrt(skal(r, r, n)) / norm_ff;

		fl = i < maxiter&& err > eps;
		skal_rr1 = skal(r, r0, n);
		for (k = 0; fl && k < n; i++, k++)
		{
			mult_mat_vec(A, z, buf);
			a_k = skal_rr1 / skal(r0, buf, n);
			for (j = 0; j < n; j++)
			{
				r1[j] = r[j];
				p[j] = r[j] - a_k * buf[j];
			}
			mult_mat_vec(A, p, buf);
			y_k = skal(p, buf, n) / skal(buf, buf, n);
			for (j = 0; j < n; j++)
			{
				x[j] += a_k * z[j] + y_k * p[j];
				r[j] = p[j] - y_k * buf[j];
			}
			skal_rr0 = skal(r, r0, n);
			b_k = a_k * skal_rr0 / (y_k * skal_rr1);
			skal_rr1 = skal_rr0;
			mult_mat_vec(A, z, buf);
			for (j = 0; j < n; j++)
			{
				z[j] = r[j] + b_k * (r1[j] - y_k * buf[j]);
			}
			err = sqrt(skal(r, r, n)) / norm_ff;
			fl = i < maxiter&& err > eps;
			//cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
		}
		itarations = i;
		nev = err;
	}
	if (i == maxiter) return 1;
	return 0;
}



int mult_mat_vec(matrix& A, real* vec, real* res)
{
	int n = A.n;

	real* A_di = A.di;
	real* A_ggl = A.ggl;
	real* A_ggu = A.ggu;
	int* A_ig = A.ig;
	int* A_jg = A.jg;

	int i, j, ij;
	int i_beg, i_end;

	for (i = 0; i < n; i++)
	{
		res[i] = 0;
		i_beg = A_ig[i];
		i_end = A_ig[i + 1];
		for (ij = i_beg; ij < i_end; ij++)
		{
			j = A_jg[ij];
			res[i] += A_ggl[ij] * vec[j];
			res[j] += A_ggu[ij] * vec[i];
		}
		res[i] += A_di[i] * vec[i];
	}


	return 0;
}

real skal(real* vec1, real* vec2, int n)
{
	real sum = 0;
	for (int i = 0; i < n; i++)
		sum += vec1[i] * vec2[i];
	return sum;
}

int calc_Lx(matrix& LU, real* x, real* f, bool is_diag_1)
{
	int n = LU.n;

	real* L_di = LU.di;
	real* L_ggl = LU.ggl;
	int* L_ig = LU.ig;
	int* L_jg = LU.jg;

	int i, j;
	int i_beg, i_end;
	int ij;
	if (is_diag_1)
		for (i = 0; i < n; i++)
		{
			x[i] = f[i];
			i_beg = L_ig[i];
			i_end = L_ig[i + 1];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = L_jg[ij];
				x[i] -= L_ggl[ij] * x[j];
			}
		}
	else
		for (i = 0; i < n; i++)
		{
			x[i] = f[i];
			i_beg = L_ig[i];
			i_end = L_ig[i + 1];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = L_jg[ij];
				x[i] -= L_ggl[ij] * x[j];
			}
			x[i] /= L_di[i];
		}
	return 0;
}

int calc_Ux(matrix& LU, real* x, real* f, bool is_diag_1)
{
	int n = LU.n;

	real* U_di = LU.di;
	real* U_ggu = LU.ggu;
	int* U_ig = LU.ig;
	int* U_jg = LU.jg;

	int i, j;
	int i_beg, i_end;
	int ij;
	for (i = 0; i < n; i++)
		x[i] = f[i];

	if (is_diag_1)
		for (i = n - 1; i > -1; i--)
		{
			i_beg = U_ig[i];
			i_end = U_ig[i + 1];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = U_jg[ij];
				x[j] -= U_ggu[ij] * x[i];
			}
		}
	else
		for (i = n - 1; i > -1; i--)
		{
			i_beg = U_ig[i];
			i_end = U_ig[i + 1];
			x[i] /= U_di[i];
			for (ij = i_beg; ij < i_end; ij++)
			{
				j = U_jg[ij];
				x[j] -= U_ggu[ij] * x[i];
			}
		}

	return 0;
}

int count_LOS(matrix& A, matrix& LU, real* x, real* r, real* z, real* p, real* f, real* buf_v, real* buf_v1, int maxiter, real eps, bool is_L_diag_1, bool is_U_diag_1)
{
	int n = A.n;
	int i, j, k;
	real norm_ff = sqrt(skal(f, f, n));
	real err;

	bool fl = true;
	i = 0;
	while (fl)
	{
		mult_mat_vec(A, x, buf_v);
		for (j = 0; j < n; j++)
			buf_v[j] = f[j] - buf_v[j];

		calc_Lx(LU, r, buf_v, is_L_diag_1);
		calc_Ux(LU, z, r, is_U_diag_1);

		mult_mat_vec(A, z, buf_v);
		calc_Lx(LU, p, buf_v, is_L_diag_1);

		real a_k, b_k;
		real skal_pp;

		err = sqrt(skal(r, r, n)) / norm_ff;
		fl = i < maxiter&& err > eps;
		for (k = 0; fl && k < n; i++, k++)
		{

			skal_pp = skal(p, p, n);
			a_k = skal(p, r, n) / skal_pp;

			for (j = 0; j < n; j++)
			{
				x[j] += a_k * z[j];
				r[j] -= a_k * p[j];
			}

			calc_Ux(LU, buf_v, r, is_U_diag_1);

			mult_mat_vec(A, buf_v, buf_v1);
			calc_Lx(LU, buf_v, buf_v1, is_L_diag_1);

			b_k = -skal(p, buf_v, n) / skal_pp;
			calc_Ux(LU, buf_v1, r, is_U_diag_1);
			for (j = 0; j < n; j++)
			{
				z[j] = buf_v1[j] + b_k * z[j];
				p[j] = buf_v[j] + b_k * p[j];
			}
			err = sqrt(skal(r, r, n)) / norm_ff;
			fl = i < maxiter&& err > eps;
			cout << setprecision(14) << "iteration: " << i << "; err: " << err << endl;
			itarations = i;
			nev = err;
		}
	}
	return 0;
}

int write(real* x, int n, ofstream& out)
{

	for (int i = 0; i < n; i++)
	{
		cout << x[i] << endl;
	}
	return 0;
}

int decomp_mat_DI(matrix& A, matrix& LU, bool is_DD)
{
	int n = A.n;

	real* A_di = A.di;

	real* LU_di = LU.di;
	real* LU_ggl = LU.ggl;
	real* LU_ggu = LU.ggu;
	int* LU_ig = LU.ig;
	int* LU_jg = LU.jg;

	int i = 0;
	if (is_DD)
	{
		for (i = 0; i < n; i++)
		{
			LU_di[i] = sqrt(A_di[i]);
			LU_ig[i] = 0;
		}

	}
	else
	{
		for (i = 0; i < n; i++)
		{
			LU_di[i] = A_di[i];
			LU_ig[i] = 0;
		}

	}
	LU_ig[n] = 0;


	return 0;
}



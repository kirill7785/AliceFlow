#ifndef GMRES_CPP
#define GMRES_CPP 1


// Скалярное произведение двух векторов
/*
template <typename doublerealT>
doublerealT Scal(doublerealT* &v1, doublerealT* &v2, integer n) {
	doublerealT s = 0.0;
	integer i; // счётчик цикла for


			   // Возможно не имеет смысла так часто постоянно задавать количество потоков.
			   //#ifdef _OPENMP
			   //omp_set_num_threads(inumcore);
			   //#endif


#pragma omp parallel for shared(v1, v2, n) private(i) reduction (+: s) schedule (guided)
	for (i = 0; i<n; i++)
	{
		s += v1[i] * v2[i];
	}
	return s;
} // Scal 
*/
/*
inline float Scal(float* &v1, float* &v2, integer n) {
	return Scal<float>(v1, v2, n);
}

inline double Scal(double* &v1, double* &v2, integer n) {
	return Scal<double>(v1, v2, n);
}
*/

double Scal(double *v1, double *v2, integer n) {
double s = 0.0;
integer i; // счётчик цикла for


// Возможно не имеет смысла так часто постоянно задавать количество потоков.
//#ifdef _OPENMP
//omp_set_num_threads(inumcore);
//#endif


#pragma omp parallel for shared(v1, v2, n) private(i) reduction (+: s) schedule (guided)
for (i = 0; i<n; i++)
{
s += v1[i] * v2[i];
}
return s;
} // Scal

float Scal(float *v1, float *v2, integer n) {
float s = 0.0;
integer i; // счётчик цикла for


// Возможно не имеет смысла так часто постоянно задавать количество потоков.
//#ifdef _OPENMP
//omp_set_num_threads(inumcore);
//#endif


#pragma omp parallel for shared(v1, v2, n) private(i) reduction (+: s) schedule (guided)
for (i = 0; i<n; i++)
{
s += v1[i] * v2[i];
}
return s;
} // Scal

/* Описание стандарта хранения CRS:
*  1. val - ненулевые значения элементов матрицы отсортированные
*  по номерам строк (нумерация начинается с нуля).
*  2. col_ind - соответствующие элементам из val номера столбцов.
*  3. row_ptr - используется для определения начала следующей строки.
*  Пример:
*
*  9.0   0.0   0.0   3.0   1.0   0.0   1.0
*  0.0   11.0   2.0   1.0   0.0   0.0   2.0
*  0.0   2.0   10.0   2.0   0.0   0.0   0.0
*  3.0   1.0   2.0   9.0   1.0   0.0   0.0
*  1.0   0.0   0.0   1.0   12.0   0.0   1.0
*  0.0   0.0   0.0   0.0   0.0   8.0   0.0
*  1.0   2.0   0.0   0.0   1.0   0.0   8.0
*
*------------- Разреженная матрица ------------
* Формат хранения: CRS
* val:      9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
* col_ind:  0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
* row_ptr:  0 4 8 10 12 14 15 16
*------------------------------------------------------
*/


// умножение матрицы на вектор
// используется формат хранения CRS
// Разреженная матрица A (val, col_ind, row_ptr) квадратная размером nxn.
// Число уравнений равно числу неизвестных и равно n.
// Уданной функции три эквивалентных представления отличающихся лишь типом аргументов. 13.января.2018
/*
template <typename doublerealT1, typename doublerealT2>
void MatrixCRSByVector(doublerealT1* val, integer* col_ind, integer* row_ptr, doublerealT2* V, doublerealT2* &tmp, integer n)
{


	// вектор tmp индексируется начиная с нуля так же как и вектор V
#pragma omp parallel for
	for (integer i = 0; i<n; i++) tmp[i] = 0.0;

	// В целях увеличения быстродействия
	// вся необходимая память выделяется заранее.
	//if (tmp == NULL)
	//{
	//printf("malloc: out of memory for vector tmp in MatrixCRSByVector\n"); // нехватка памяти
	//getchar();
	//exit(0);  // завершение программы
	//}



	//omp_set_num_threads(inumcore);

#pragma omp parallel for  schedule (guided)
	for (integer i = 0; i<n; i++) {
		doublerealT1 sum;
		integer rowend, rowbeg;

		sum = 0.0;
		rowend = row_ptr[i + 1];
		rowbeg = row_ptr[i];
		for (integer j = rowbeg; j<rowend; j++)
		{
			sum += val[j] * V[col_ind[j]];
		}
		tmp[i] = sum;
	}

	//return tmp;
} // MatrixCRSByVector
*/

// умножение матрицы на вектор
// используется формат хранения CRS
// Разреженная матрица A (val, col_ind, row_ptr) квадратная размером nxn.
// Число уравнений равно числу неизвестных и равно n.
// Уданной функции три эквивалентных представления отличающихся лишь типом аргументов. 13.января.2018
void MatrixCRSByVector(double* val, integer* col_ind, integer* row_ptr, double* V, double* &tmp, integer n)
{


	// вектор tmp индексируется начиная с нуля так же как и вектор V
#pragma omp parallel for
	for (integer i = 0; i<n; i++) tmp[i] = 0.0;

	// В целях увеличения быстродействия
	// вся необходимая память выделяется заранее.
	//if (tmp == NULL)
	//{
	//printf("malloc: out of memory for vector tmp in MatrixCRSByVector\n"); // нехватка памяти
	//getchar();
	//exit(0);  // завершение программы
	//}



	//omp_set_num_threads(inumcore);

#pragma omp parallel for  schedule (guided)
	for (integer i = 0; i<n; i++) {
		double sum;
		integer rowend, rowbeg;

		sum = 0.0;
		rowend = row_ptr[i + 1];
		rowbeg = row_ptr[i];
		for (integer j = rowbeg; j<rowend; j++)
		{
			sum += val[j] * V[col_ind[j]];
		}
		tmp[i] = sum;
	}

	//return tmp;
} // MatrixCRSByVector

  // Уданной функции три эквивалентных представления отличающихся лишь типом аргументов. 13.января.2018
void MatrixCRSByVector(float* val, integer* col_ind, integer* row_ptr, float* V, float* &tmp, integer n)
{


	// вектор tmp индексируется начиная с нуля так же как и вектор V
#pragma omp parallel for
	for (integer i = 0; i<n; i++) tmp[i] = 0.0;

	// В целях увеличения быстродействия
	// вся необходимая память выделяется заранее.
	//if (tmp == NULL)
	//{
	//printf("malloc: out of memory for vector tmp in MatrixCRSByVector\n"); // нехватка памяти
	//getchar();
	//exit(0);  // завершение программы
	//}



	//omp_set_num_threads(inumcore);

#pragma omp parallel for  schedule (guided)
	for (integer i = 0; i<n; i++) {
		float sum;
		integer rowend, rowbeg;

		sum = 0.0;
		rowend = row_ptr[i + 1];
		rowbeg = row_ptr[i];
		for (integer j = rowbeg; j<rowend; j++)
		{
			sum += val[j] * V[col_ind[j]];
		}
		tmp[i] = sum;
	}

	//return tmp;
} // MatrixCRSByVector

  // Уданной функции три эквивалентных представления отличающихся лишь типом аргументов. 13.января.2018
void MatrixCRSByVector(double* val, integer* col_ind, integer* row_ptr, float* V, float* &tmp, integer n)
{


	// вектор tmp индексируется начиная с нуля так же как и вектор V
#pragma omp parallel for
	for (integer i = 0; i<n; i++) tmp[i] = 0.0;

	// В целях увеличения быстродействия
	// вся необходимая память выделяется заранее.
	//if (tmp == NULL)
	//{
	//printf("malloc: out of memory for vector tmp in MatrixCRSByVector\n"); // нехватка памяти
	//getchar();
	//exit(0);  // завершение программы
	//}



	//omp_set_num_threads(inumcore);

#pragma omp parallel for  schedule (guided)
	for (integer i = 0; i<n; i++) {
		double sum;
		integer rowend, rowbeg;

		sum = 0.0;
		rowend = row_ptr[i + 1];
		rowbeg = row_ptr[i];
		for (integer j = rowbeg; j<rowend; j++)
		{
			sum += val[j] * V[col_ind[j]];
		}
		tmp[i] = sum;
	}

	//return tmp;
} // MatrixCRSByVector


// 13 сентября 2017 года. Алгоритм Саада и Шульца.
// Заимствован из iterative template library. 
// Выделение и уничтожение памяти происходит внутри алгоритма.
//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//    dxret  --  approximate solution to Ax = b
//    maxit  --  the number of iterations performed before the
//               tolerance was reached
//      dterminatedTResudual  --  the residual after the final iteration
//  
//*****************************************************************

template <typename doublerealT>
void mult_givens(doublerealT c, doublerealT s, integer k, doublerealT* &g)

//****************************************************************************80
//
//  Purpose:
//
//    MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 August 2006
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
//    Charles Romine, Henk van der Vorst,
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Methods,
//    SIAM, 1994,
//    ISBN: 0898714710,
//    LC: QA297.8.T45.
//
//    Tim Kelley,
//    Iterative Methods for Linear and Nonlinear Equations,
//    SIAM, 2004,
//    ISBN: 0898713528,
//    LC: QA297.8.K45.
//
//    Yousef Saad,
//    Iterative Methods for Sparse Linear Systems,
//    Second Edition,
//    SIAM, 2003,
//    ISBN: 0898715342,
//    LC: QA188.S17.
//
//  Parameters:
//
//    Input, double C, S, the cosine and sine of a Givens
//    rotation.
//
//    Input, int K, indicates the location of the first vector entry.
//
//    Input/output, double G[K+2], the vector to be modified.  On output,
//    the Givens rotation has been applied to entries G(K) and G(K+1).
//
{
	doublerealT g1;
	doublerealT g2;

	g1 = c * g[k] - s * g[k + 1];
	g2 = s * g[k] + c * g[k + 1];

	g[k] = g1;
	g[k + 1] = g2;

	return;
}

// Норма вектора
// как корень квадратный из суммы квадратов.
double NormaV_for_gmres( double *dV, integer isize)
{
	integer i; // Счетчик цикла
	double dnorma, dsum;

	// инициализация переменных
	dsum = 0.0;
	for (i = 0; i <= (isize - 1); i++)
	{
		dsum += dV[i] * dV[i];
	}
	dnorma = sqrt(dsum); // норма вектора
	return dnorma;
}// NormaV_for_gmres

 // Норма вектора
 // как корень квадратный из суммы квадратов.
float NormaV_for_gmres(float *dV, integer isize)
{
	integer i; // Счетчик цикла
	float dnorma, dsum;

	// инициализация переменных
	dsum = 0.0;
	for (i = 0; i <= (isize - 1); i++)
	{
		dsum += dV[i] * dV[i];
	}
	dnorma = sqrt(dsum); // норма вектора
	return dnorma;
}// NormaV_for_gmres


void Update(double* &x, integer k, integer n, double** &h, double* &s, double** &v)
{//ok
	//Vector y(s);
	double* y = new double[k + 1];
	for (integer i_1 = 0; i_1 <= k; i_1++) y[i_1] = s[i_1];

	// Backsolve:  
	for (integer i = k; i >= 0; i--) {
		y[i] /= h[i][i];
		for (integer j = i - 1; j >= 0; j--)
			y[j] -= h[j][i] * y[i];
	}

	for (integer j = 0; j <= k; j++) {
		for (integer j_1 = 0; j_1 < n; j_1++) {
			x[j_1] += v[j][j_1] * y[j];
		}
	}


	delete[] y;
	y = NULL;
}

void Update(float* &x, integer k, integer n, float** &h, float* &s, float** &v)
{//ok
 //Vector y(s);
	float* y = new float[k + 1];
	for (integer i_1 = 0; i_1 <= k; i_1++) y[i_1] = s[i_1];

	// Backsolve:  
	for (integer i = k; i >= 0; i--) {
		y[i] /= h[i][i];
		for (integer j = i - 1; j >= 0; j--)
			y[j] -= h[j][i] * y[i];
	}

	for (integer j = 0; j <= k; j++) {
		for (integer j_1 = 0; j_1 < n; j_1++) {
			x[j_1] += v[j][j_1] * y[j];
		}
	}


	delete[] y;
	y = NULL;
}


template <typename doublerealT>
void print_Hessenberg(integer k, doublerealT** &h, doublerealT beta) {
	printf("%e\n",beta);
	for (integer i_1 = 0; i_1 <= k; i_1++) {
		for (integer i_2 = 0; i_2 <= k; i_2++) {
			printf("%1.2e ", h[i_1][i_2]);
		}
		printf("\n");
	}
	getchar();
}

void Update_flexible(doublereal* &x, integer k, integer n, doublereal** &h, doublereal* &s, doublereal** &v, doublereal beta, doublereal** &Z)
{//ok
	//Vector y(s);
	doublereal* y = new doublereal[k + 1];
	for (integer i_1 = 0; i_1 <= k; i_1++) y[i_1] = beta*s[i_1];
	//for (integer i_1 = 0; i_1 <= k; i_1++) y[i_1] = 0.0;
	//y[0] = beta;

	//debug
	//print_Hessenberg(k, h, beta);

	// Backsolve:  	
	for (integer i = k; i >= 0; i--) {
		y[i] /= h[i][i];
		for (integer j = i - 1; j >= 0; j--)
			y[j] -= h[j][i] * y[i];
		
	}
	
	/*
	// Backsolve:  
	for (integer i = 0; i <= k; i++) {
		// моя вставка.
		y[i] /= h[i][i];
		for (integer j = i - 1; j >= 0; j--) {
			y[i] -= (h[j][i] / h[i][i]) * y[j];
			//y[j] -= h[j][i] * y[i];
		}		
	}
	*/
	//printf("y=\n");
	//for (integer i = 0; i <= k; i++) {
		//printf("%e ",y[i]);
	//}
	//getchar();

	//for (integer i_1 = 0; i_1 <= k; i_1++) y[i_1] *= beta;

	for (integer j = 0; j <= k; j++) {
		for (integer j_1 = 0; j_1 < n; j_1++) {
			x[j_1] += Z[j][j_1] * y[j];//y[j]
		}
	}


	delete[] y;
	y = NULL;
}


void GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn)
{//ok
	if (fabs(dy) < 1.0e-30) {
		cs = 1.0;
		sn = 0.0;
	}
	else if (fabs(dy) > fabs(dx)) {
		double temp = dx / dy;
		sn = 1.0 / sqrt(1.0 + temp*temp);
		cs = temp * sn;
	}
	else {
		double temp = dy / dx;
		cs = 1.0 / sqrt(1.0 + temp*temp);
		sn = temp * cs;
	}
}

void GeneratePlaneRotation(float &dx, float &dy, float &cs, float &sn)
{//ok
	if (fabs(dy) < 1.0e-30) {
		cs = 1.0;
		sn = 0.0;
	}
	else if (fabs(dy) > fabs(dx)) {
		float temp = dx / dy;
		sn = 1.0 / sqrt(1.0 + temp*temp);
		cs = temp * sn;
	}
	else {
		float temp = dy / dx;
		cs = 1.0 / sqrt(1.0 + temp*temp);
		sn = temp * cs;
	}
}


void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{//ok
	double temp = cs * dx + sn * dy;
	dy = -sn * dx + cs * dy;
	dx = temp;
}

void ApplyPlaneRotation(float &dx, float &dy, float &cs, float &sn)
{//ok
	float temp = cs * dx + sn * dy;
	dy = -sn * dx + cs * dy;
	dx = temp;
}

// Если использовать gmres с предобуславливателем AINV Bridson то нужно брать m_restart=32. 
// Так рекомендуют в статье для решения уравнений теплопередачи с числом неизвестных до 1.5М.
// Если использовать GMRES в качестве предобуславливателя то рекомендуют брать m_restart=2.
// Даный метод существенно замедляется при увеличении числа итераций. Например он не сходится на Диффузионно дрейфовой модели для 
// материала GaN в то время как BiCGStab сходится.
template <typename doublerealT>
integer  gmres(integer n, doublerealT *val, integer* col_ind, integer* row_ptr, doublerealT *dV, doublerealT* &dX0,
	integer maxit, integer &m_restart) {

	dterminatedTResudual = 1.0e-12;

	doublerealT resid;
	integer i, j = 1, k, iter_number=1;
	//Vector s(m + 1), cs(m + 1), sn(m + 1), w;
	doublerealT* w = new doublerealT[n];
	doublerealT* s = new doublerealT[m_restart + 2];
	doublerealT* cs = new doublerealT[m_restart + 2];
	doublerealT* sn = new doublerealT[m_restart + 2];

	

	doublerealT *dx = new doublerealT[n];

	//integer i_1 = 0;

	// начальное приближение
	// X0 ==
	// под X0 понимается вектор поля температур к примеру.
	if (dX0 == NULL) {
		dX0 = new doublerealT[n];
		for (integer i_1 = 0; i_1<n; i_1++) {
			dx[i_1] = 0.0;
			dX0[i_1] = 0.0;
		}
	}
	else {
		for (integer  i_1 = 0; i_1<n; i_1++) dx[i_1] = dX0[i_1];
	}

	//doublerealT normb = norm(M.solve(b));
	doublerealT normb = 0.0;
	// здесь реализованы все три нормы
	// вообще говоря они все эквивалентны

	normb = NormaV_for_gmres(dV, n);// ERROR На правую dV часть надо подействовать предобуславливателем

	//Vector r = M.solve(dV - A * x);
	doublerealT *r = new doublerealT[n];
	MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // результат занесён в  r
	for (integer i_1 = 0; i_1 < n; i_1++) r[i_1] = dV[i_1] - r[i_1]; // ERRoR предобуславливатель.

	//doublerealT beta = norm(r);
	doublerealT beta = 0.0;

	beta = NormaV_for_gmres(r, n);

	if (fabs(normb) < 1.0e-30)
		normb = 1;

	doublerealT norm_r = 0.0;


	norm_r = NormaV_for_gmres(r, n);

	if ((resid = norm_r / normb) <= dterminatedTResudual) {
		//tol = resid;
		maxit = 0;
		return 0;
	}

	doublerealT** H = new doublerealT*[m_restart + 2]; // Hessenberg
	for (integer  i_1 = 0; i_1 < m_restart + 2; i_1++) H[i_1] = new doublerealT[m_restart + 2];
	for (integer i_1 = 0; i_1 < m_restart + 2; i_1++)
	{
		for (integer j_1 = 0; j_1 < m_restart + 2; j_1++)
		{
			H[i_1][j_1] = 0.0;
		}
	}

	//Vector *v = new Vector[m_restart + 1];
	doublerealT** v = new doublerealT*[m_restart + 2];
	for (integer i_1 = 0; i_1 <= m_restart + 1; i_1++) v[i_1] = new doublerealT[n];
	for (integer i_1 = 0; i_1 <= m_restart + 1; i_1++) {
		for (integer j_1 = 0; j_1 < n; j_1++)
		{
			v[i_1][j_1] = 0.0;
		}
	}

	while (j <= maxit) {

		//v[0] = r * (1.0 / beta);    // ??? r / beta
		for (integer j_1 = 0; j_1 < n; j_1++)
		{
			v[0][j_1] = r[j_1] * (1.0 / beta);
		}
		//s = 0.0;
		for (integer i_1 = 0; i_1 <= m_restart + 1; i_1++) {
				s[i_1] = 0.0;
		}
		s[0] = beta;

		for (i = 0; i < m_restart && j <= maxit; i++, j++) {
			//w = M.solve(A * v[i]);
			MatrixCRSByVector(val, col_ind, row_ptr, v[i], w, n); // результат занесён в  w
			for (k = 0; k <= i; k++) {
				H[k][i] = Scal(w, v[k], n);
				for (integer j_1 = 0; j_1 < n; j_1++)
				{
					w[j_1] -= H[k][i] * v[k][j_1];
				}
			}
			//H[i + 1][i] = norm(w);
			H[i + 1][i] = NormaV_for_gmres(w, n);// евклидова

			for (integer j_1 = 0; j_1 < n; j_1++)
			{
				v[i + 1][j_1] = w[j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			}

			for (k = 0; k < i; k++)
				ApplyPlaneRotation(H[k][i], H[k + 1][i], cs[k], sn[k]);

			GeneratePlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);

			//printf("%d %e\n", iter_number++, fabs(s[i + 1]));

			if ((resid = fabs(s[i + 1]) / normb) < dterminatedTResudual) {
				Update(dx, i, n, H, s, v);
				//tol = resid;
				//maxit = j;
				for (integer i_1 = 0; i_1<n; i_1++) {
					dX0[i_1] = dx[i_1];

				}
				for (integer i_1 = 0; i_1 <= m_restart + 1; i_1++) delete[] v[i_1];
				delete[] v;
				delete[] dx;
				delete[] r;
				delete[] w;
			

				delete[] s;
				delete[] cs;
				delete[] sn;
				for (integer i_1 = 0; i_1 < m_restart + 2; i_1++) delete[] H[i_1];
				delete[] H;
				return 0;
			}
		}
		Update(dx, i - 1, n, H, s, v);
		//r = M.solve(b - A * x);
		MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // Результат занесён в r
		for (integer i_1 = 0; i_1 < n; i_1++) r[i_1] = dV[i_1] - r[i_1];
		//beta = norm(r);
		beta = NormaV_for_gmres(r, n); // евклидова.

		if ((resid = beta / normb) < dterminatedTResudual) {
			//tol = resid;
			//maxit = j;
			for (integer i_1 = 0; i_1<n; i_1++) {
				dX0[i_1] = dx[i_1];

			}
			for (integer i_1 = 0; i_1 <= m_restart + 1; i_1++) delete[] v[i_1];
			delete[] v;
			delete[] dx;
			delete[] r;
			delete[] w;
			

			delete[] s;
			delete[] cs;
			delete[] sn;
			for (integer i_1 = 0; i_1 < m_restart + 2; i_1++) delete[] H[i_1];
			delete[] H;
			return 0;
		}
	}

	//tol = resid;
	for (integer i_1 = 0; i_1<n; i_1++) {
		dX0[i_1] = dx[i_1];

	}
	for (integer i_1 = 0; i_1 <= m_restart + 1; i_1++) delete[] v[i_1];
	delete[] v;
	delete[] dx;
	delete[] r;
	delete[] w;
	

	delete[] s;
	delete[] cs;
	delete[] sn;
	for (integer i_1 = 0; i_1 < m_restart + 2; i_1++) delete[] H[i_1];
	delete[] H;
	return 1;

}

 

  /*     Last change:  ERB  22 Aug 2000   10:31 am */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  /*     AMG1R5                                        MAIN SUBROUTINE */

  /*     RELEASE 1.5, OCTOBER 1990 */

  /*     CHANGES AGAINST VERSION 1.1, JULY 1985: */

  /* 1.  A BUG WAS DETECTED WHICH UNDER CERTAIN CIRCUMSTANCES INFLUENCED */
  /*     SLIGHTLY THE CONVERGENCE RATE OF AMG1R1. FOR THAT REASON, THE */
  /*     FOLLOWING LINE IN SUBROUTINE RESC: */
  /*     IW(IMAXW(KC-1)+1) = IA(IMIN(KC)) */
  /*     HAS BEEN CHANGED TO: */
  /*     IW(IMAXW(KC-1)+1) = IAUX */

  /* 2.  A BUG WAS DETECTED IN SUBROUTINE PWINT. UNDER CERTAIN CIRCUM- */
  /*     STANCES AN UNDEFINED VARIABLE WAS USED. ALTHOUGH THIS DID NOT */
  /*     AFFECT THE NUMERICAL RESULTS, PROBLEMS CAN OCCUR IF CHECKING */
  /*     FOR UNDEFINED VARIABLES IS USED. TO FIX THIS ERROR, IN PWINT */
  /*     THE LABEL 1000 WAS MOVED TO THE STATEMENT */
  /*     IBLCK1 = IMINW(K). */

  /* 3.  A PARAMETER LRATIO HAS BEEN INTRODUCED, DENOTING THE RATIO */
  /*     OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
  /*     THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
  /*     TO 2. CHANGE THIS VALUE IF NECESSARY. (IF, FOR EXAMPLE, YOU */
  /*     WANT TO CHANGE THE DOUBLE PRECISION VECTORS TO SINGLE PRE- */
  /*     CISION, LRATIO HAS TO BE SET TO 1. IN THE YALE SMP - ROUTINE */
  /*     NDRV THERE IS A PARAMETER LRATIO, TOO. */

  /* 4.  TYPE DECLARATIONS REAL*4 AND REAL*8 HAVE BEEN CHANGED TO THE */
  /*     STANDARD-CONFORMING KEYWORDS REAL AND DOUBLE PRECISION, RESPEC- */
  /*     TIVELY. */

  /* 5.  CALLS TO THE FOLLOWING INTRINSIC FUNCTIONS HAVE BEEN REPLACED BY */
  /*     CALLS USING GENERIC NAMES: DSQRT, MIN0, MAX0, IABS, DABS, FLOAT, */
  /*     DFLOAT, DMAX1, ISIGN, IDINT, DLOG10. */

  /* 6.  A SAVE STATEMENT HAS BEEN INSERTED IN ALL SUBROUTINES. */

  /* 7.  EXTERNAL DECLARATION STATEMENTS HAVE BEEN INSERTED IN ALL SUB- */
  /*     ROUTINES FOR ALL EXTERNAL REFERENCES. */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.3, APRIL 1986: */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD FAIL UNDER CERTAIN CIRCUMSTANCES. FOR A FIX, THE FOLLOWING */
  /*     STATEMENTS IN SUBROUTINE CHECK HAVE BEEN CHANGED: */

  /*     DO 450 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 450 J=IA(I)+1,ICG(I)-1 */

  /*     DO 430 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 430 J1=IA(I1)+1,ICG(I1)-1 */

  /*     DO 550 J=IA(I)+1,IA(I+1)-1 WAS CHANGED TO */
  /*     DO 550 J=IA(I)+1,ICG(I)-1 */

  /*     DO 530 J1=IA(I1)+1,IA(I1+1)-1 WAS CHANGED TO */
  /*     DO 530 J1=IA(I1)+1,ICG(I1)-1 */

  /* 2.  THE EXPLANATORY PART IN SUBROUTINE AMG1R5 HAS BEEN ENLARGED TO */
  /*     AVOID MISUNDERSTANDINGS IN THE DEFINITION OF THE ARGUMENT LIST. */

  /* ----------------------------------------------------------------------- */

  /*     CHANGE AGAINST VERSION 1.4, OCTOBER, 1990 (BY JOHN W. RUGE) */

  /* 1.  A BUG IN SUBROUTINE CHECK HAS BEEN REMOVED. IF THE ORIGINAL MATRIX */
  /*     WAS STORED IN AN UNSYMMETRIC WAY, THE SYMMETRIZATION BY AMG1R3 */
  /*     COULD STILL FAIL UNDER CERTAIN CIRCUMSTANCES, AND WAS NOT FIXED */
  /*     IN THE PREVIOUS VERSION. IN ADDITION, THE ROUTINE WAS CHANGED */
  /*     IN ORDER TO AVOID SOME UNNECESSARY ROW SEARCHES FOR TRANSOSE */
  /*     ENTRIES. */

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/* Subroutine */ integer amg1r5_fgmres_version(doublereal *a, integer *ia, integer *ja,
	doublereal *u, doublereal *f, integer *ig, integer *nda, integer *
	ndia, integer *ndja, integer *ndu, integer *ndf, integer *ndig,
	integer *nnu, integer *matrix, integer *iswtch, integer *iout,
	integer *iprint, integer *levelx, integer *ifirst, integer *ncyc,
	doublereal *eps, integer *madapt, integer *nrd, integer *nsolco,
	integer *nru, doublereal *ecg1, doublereal *ecg2, doublereal *ewt2,
	integer *nwt, integer *ntr, integer *ierr, integer iVar,
	equation3D* &sl, equation3D_bon* &slb, integer maxelm, integer maxbound,
	bool &bOkfgmres_amg1r5)
{
	// информирует внешний вызывающий код сощелся алгоритм или нет.
	bOkfgmres_amg1r5 = false;

	// 31 декабря 2017.
	// FGMRES - гибкий вариант обобщённого метода минимальной невязки. (Если не BiCGStab то FGMRES).
	// В данной функции amg1r5 алгоритм является предобуславливателем для алгоритма 
	// Ю. Саада и Шульца FGMRES[1986].

	// sl, slb - Матрица СЛАУ для алгоритма FGMRES [1986] Ю.Саада и Шульца.
	// Матрица СЛАУ (a, ia, ja) модифицируется в процессе работы алгебраического многосеточного метода.
	// априори: nnu[0]==maxelm+maxbound
	// ndu[0]==ndf[0] ...


	/* Format strings */


	/* Builtin functions */


	/* Local variables */
	integer mda = 0, mdf = 0, mdu = 0;
	doublereal res = 0.0;
	integer ium = 0, iup = 0;
	doublereal res0 = 0.0;
	extern /* Subroutine */ integer idec_(integer *, integer *, integer *,
		integer *);
	integer mdia = 0, mdja = 0, mdig = 0, iarr[25] = { 0 };
	//static real time[20];
	unsigned int time[20] = { 0 };
	integer imin[25] = { 0 }, imax[25] = { 0 };
	doublereal resi[25] = { 0.0 };
	integer kout = 0, ncyc0 = 0, irow0 = 0, ndicg = 0, icgst = 0, iminw[25] = { 0 }, imaxw[25] = { 0 };
	extern /* Subroutine */ integer first_(integer *, doublereal *, integer *,
		integer *, integer *, integer *), solve_(integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, doublereal
			*, doublereal *, doublereal *, integer *, integer *, integer *,
			doublereal *, integer *, integer *, integer *, integer *, integer
			*, integer *, integer *, integer *, /*real*/ unsigned int *, integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, integer *,
			integer *, integer *, integer *, integer *, integer *, doublereal
			*, doublereal *, doublereal *), setup_(integer *, integer *,
				integer *, doublereal *, doublereal *, doublereal *, integer *,
				integer *, integer *, doublereal *, doublereal *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *, integer *, integer *,
				integer *, integer *, integer *, integer *);
	integer ndigit = 0, levels = 0, kevelx = 0, nstcol[25] = { 0 }, kswtch = 0;
	extern /* Subroutine */ integer wrkcnt_(integer *, integer *, integer *,
		integer *, integer *, integer *, integer *, /*real*/ unsigned int *, integer *,
		integer *, integer *, integer *, integer *, integer *, integer *,
		integer *, doublereal *, doublereal *);

	/* Fortran I/O blocks */

	integer nnz;
	// нумерация векторов начинается с нуля.
	integer n75 = maxelm + maxbound; // число неизвестных на подробном уровне.
	doublereal* val75 = NULL;
	//val75 = new doublereal[nnz];
	integer* col_ind75 = NULL;
	//col_ind75 = new integer[nnz];
	integer* row_ptr75 = NULL;
	//row_ptr75 = new integer[n75 + 1];
	bool bnorelax = true; // Для уравнения теплопроводности не используется релаксация.
	integer m_restart = my_amg_manager.m_restart;

	doublereal resid;
	integer i, j = 1, k, iter_number;
	doublereal* w;
	doublereal* s;
	doublereal* cs;
	doublereal* sn;
	doublereal *dx;
	doublereal *buffer;
	doublereal *Zcopy;
	doublereal *vCopy;

	doublereal normb = 0.0;
	doublereal *r;
	doublereal beta = 0.0;
	doublereal norm_r = 0.0;
	integer maxit = 2000;
	integer i_1 = 0; // счётчик цикла for

	doublereal** v;
	doublereal** Z;
	integer i_copy;
	doublereal** H;




	/*         ----------------------------------------------- */
	/*         | AMG-MODULE FOR SOLVING LINEAR SYSTEMS L*U=F | */
	/*         ----------------------------------------------- */

	/*         -------------------------------------------------------------- */

	/*     ASSUMPTIONS ON L: */

	/*         THE PROGRAM REQUIRES: */ // Требования к входным данным:

										/*             - DIAGONAL ENTRIES ARE ALWAYS POSITIVE (ON ALL GRIDS); */
										/*             - L IS A SQUARE MATRIX WHICH IS EITHER REGULAR OR SINGULAR */
										/*               WITH ROWSUMS=0. */
										// Матрица L это квадратная матрица с всегда положительными диагональными элементами с нулевой суммой коэффициентов в каждой строке (положительная определённость).

										/*         FOR THEORETICAL REASONS THE FOLLOWING SHOULD HOLD: */

										/*             - L POSITIVE DEFINITE (OR SEMI-DEFINITE WITH ROWSUM=0) */
										/*             - L "ESSENTIALLY" POSITIVE TYPE, I.E., */

										/*                  -- DIAGONAL ENTRIES MUST BE > 0 ; */
										/*                  -- MOST OF THE OFF-DIAGONAL ENTRIES <= 0 ; */
										/*                  -- ROWSUMS SHOULD BE >= 0 . */

										// Матрица L положительно определённая : диагональные элементы  строго больше нуля, большинство внедиагональных элементов <= 0.
										// Сумма коэффициентов в строке больше либо равна нулю - диагональное преобладание.


										/*     THE USER HAS TO PROVIDE THE MATRIX L, THE RIGHT HAND SIDE F AND */
										/*     CERTAIN POINTER VECTORS IA AND JA. */

										// Пользователь задаёт матрицу коэффициентов L, правую часть F, а также информацию о связях между элементами матрицы в специальных столбцах IA и JA.
										// IA - позиция первого элемента в строке. JA - номер столбца для каждого элемента, первым записывается диагональный элемент.
										// В фортране нумерация начинается с единицы.

										/*         -------------------------------------------------------------- */

										/*     STORAGE OF L: */ // Требования к хранению матрицы L.

																/*         THE NON-ZERO ENTRIES OF THE MATRIX L ARE STORED IN */
																/*         "COMPRESSED" SKY-LINE FASHION IN A 1-D VECTOR A, I.E., ROW */
																/*         AFTER ROW, EACH ROW STARTING WITH ITS DIAGONAL ELEMENT. THE */
																/*         OTHER NON-ZERO ROW ENTRIES FOLLOW THEIR DIAGONAL ENTRY IN ANY */
																/*         ORDER. */


																/*         IN ORDER TO IDENTIFY EACH ELEMENT IN A, THE USER HAS TO */
																/*         PROVIDE TWO POINTER ARRAYS IA AND JA. IF NNU DENOTES THE TOTAL */
																/*         NUMBER OF UNKNOWNS, THE NON-ZERO ENTRIES OF ANY ROW I OF L */
																/*         (1.LE.I.LE.NNU) ARE STORED IN A(J) WHERE THE RANGE OF J */
																/*         IS GIVEN BY */

																/*                     IA(I) .LE. J .LE. IA(I+1)-1. */

																/*         THUS, IA(I) POINTS TO THE POSITION OF THE DIAGONAL ENTRY OF */
																/*         ROW I WITHIN THE VECTOR A. IN PARTICULAR, */

																/*                     IA(1) = 1 ,  IA(NNU+1) = 1 + NNA */

																/*         WHERE NNA DENOTES THE TOTAL NUMBER OF MATRIX ENTRIES STORED. */
																/*         THE POINTER VECTOR JA HAS TO BE DEFINED SUCH THAT */
																/*         ANY ENTRY A(J) CORRESPONDS TO THE UNKNOWN U(JA(J)), I.E., */
																/*         JA(J) POINTS TO THE COLUMN INDEX OF A(J). */
																/*         IN PARTICULAR, A(IA(I)) IS THE DIAGONAL ENTRY OF ROW I */
																/*         AND CORRESPONDS TO THE UNKNOWN U(I): JA(IA(I))=I. */

																/*         IN THIS TERMINOLOGY, THE I-TH EQUATION READS AS FOLLOWS */
																/*         (FOR ANY I WITH  1.LE.I.LE.NNU): */

																/*                  F(I) =        SUM      A(J) * U(JA(J)) */
																/*                           J1.LE.J.LE.J2 */

																/*         WHERE F(I) DENOTES THE I-TH COMPONENT OF THE RIGHT HAND */
																/*         SIDE AND */

																/*                     J1 = IA(I) ,  J2 = IA(I+1)-1. */

																/*         NOTES: THE ENTRY IA(NNU+1) HAS TO TOCHKA TO THE FIRST FREE */
																/*                ENTRY IN VECTORS A AND JA, RESPECTIVELY. OTHERWISE, */
																/*                AMG CANNOT KNOW THE LENGTH OF THE LAST MATRIX ROW. */

																/*                THE INPUT VECTORS A, IA AND JA ARE CHANGED BY AMG1R5. */
																/*                SO, AFTER RETURN FROM AMG1R5, THE PACKAGE MUST NOT */
																/*                BE CALLED A SECOND TIME WITHOUT HAVING NEWLY DEFINED */
																/*                THE INPUT VECTORS AND USING ISWTCH=4. OTHERWISE, THE */
																/*                SETUP PHASE WILL FAIL. */
																/*                  ON THE OTHER HAND, RUNNING AMG A SECOND TIME ON THE */
																/*                SAME INPUT DATA WITH ISWTCH=4 HAS NO SENSE, BECAUSE */
																/*                THE RESULTS OF THE FIRST SETUP PHASE ARE STILL STORED */
																/*                AND THUS THIS PHASE CAN BE SKIPPED IN A SECOND CALL. */
																/*                IN ORDER TO DO THIS, SET ISWTCH TO 1, 2 OR 3. */

																/* ----------------------------------------------------------------------- */

																/*         THE FORM OF THE CALLING PROGRAM HAS TO BE AS FOLLOWS: */

																/*               PROGRAM DRIVER */
																/*         C */
																/*               DOUBLE PRECISION A(#NDA),U(#NDU),F(#NDF) */
																/*               INTEGER IA(#NDIA),JA(#NDJA),IG(#NDIG) */
																/*         C */
																/*               NDA  = #NDA */
																/*               NDU  = #NDU */
																/*               NDF  = #NDF */
																/*               NDIA = #NDIA */
																/*               NDJA = #NDJA */
																/*               NDIG = #NDIG */
																/*         C */
																/*         C     SET UP A, F, IA, JA AND SPECIFY NECESSARY PARAMETERS */
																/*         C */
																/*               .... */
																/*               .... */
																/*         C */
																/*               CALL AMG1R5(A,IA,JA,U,F,IG, */
																/*        +                  NDA,NDIA,NDJA,NDU,NDF,NDIG,NNU,MATRIX, */
																/*        +                  ISWTCH,IOUT,IPRINT, */
																/*        +                LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU, */
																/*        +                  ECG1,ECG2,EWT2,NWT,NTR, */
																/*        +                  IERR) */
																/*         C */
																/*               .... */
																/*               .... */
																/*         C */
																/*               STOP */
																/*               END */

																/* ----------------------------------------------------------------------- */

																/*     INPUT VIA ARRAYS (SEE ABOVE): */

																/*     A        -   MATRIX L */

																/*     IA       -   POINTER VECTOR */

																/*     JA       -   POINTER VECTOR */

																/*     U        -   FIRST APPROXIMATION TO SOLUTION */

																/*     F        -   RIGHT HAND SIDE */


																/* ----------------------------------------------------------------------- */


																/*     SCALAR INPUT PARAMETERS OF AMG1R5: */

																/*     THE INPUT PARAMETERS OF AMG1R5 IN THE LIST BELOW ARE ARRANGED */
																/*     ACCORDING TO THEIR IMPORTANCE TO THE GENERAL USER. THE PARAMETERS */
																/*     PRECEEDED BY A * MUST BE SPECIFIED EXPLICITELY. ALL THE OTHER */
																/*     PARAMETERS ARE SET TO STANDARD VALUES IF ZERO ON INPUT. */

																/*     THERE ARE FOUR CLASSES OF INPUT PARAMETERS WITH DECREASING PRI- */
																/*     ORITY: */

																/*     1. PARAMETERS DESCRIBING THE USER-DEFINED PROBLEM AND DIMENSIONING */
																/*        OF VECTORS IN THE CALLING PROGRAM */

																/*     2. PARAMETERS SPECIFYING SOME GENERAL ALGORITHMIC ALTERNATIVES AND */
																/*        THE AMOUNT OF OUTPUT DURING SOLUTION */

																/*     3. PARAMETERS CONTROLLING THE MULTIGRID CYCLING DURING THE SOLU- */
																/*        TION PHASE */

																/*     4. PARAMETERS CONTROLLING THE CREATION OF COARSER GRIDS AND INTER- */
																/*        POLATION FORMULAS. */

																/*     ONLY THE CLASS 1 - PARAMETERS MUST BE SPECIFIED EXPLICITELY BY */
																/*     THE USER. CLASS 2 - PARAMETERS CONTROL THE GENERAL PERFORMANCE OF */
																/*     AMG1R5. CHANGING THEM DOESN'T REQUIRE UNDERSTANDING THE AMG - */
																/*     ALGORITHM. SPECIFYING NON-STANDARD-VALUES FOR CLASS 3 - PARAMETERS */
																/*     PRESUPPOSES A GENERAL KNOWLEDGE OF MULTIGRID METHODS, WHEREAS THE */
																/*     FUNCTION OF CLASS 4 - PARAMETERS IS ONLY UNDERSTANDABLE AFTER */
																/*     STUDYING THE AMG-ALGORITHM IN DETAIL. FORTUNATELY IN MOST CASES */
																/*     THE CHOICE OF CLASS 3 AND 4 - PARAMETERS ISN'T CRITICAL AND USING */
																/*     THE AMG1R5 - SUPPLIED STANDARD VALUES SHOULD GIVE SATISFACTORY */
																/*     RESULTS. */

																/*         -------------------------------------------------------------- */

																/*     CLASS 1 - PARAMETERS: */

																/*  *  NDA      -   DIMENSIONING OF VECTOR A IN CALLING PROGRAM */

																/*  *  NDIA     -   DIMENSIONING OF VECTOR IA IN CALLING PROGRAM */

																/*  *  NDJA     -   DIMENSIONING OF VECTOR JA IN CALLING PROGRAM */

																/*  *  NDU      -   DIMENSIONING OF VECTOR U IN CALLING PROGRAM */

																/*  *  NDF      -   DIMENSIONING OF VECTOR F IN CALLING PROGRAM */

																/*  *  NDIG     -   DIMENSIONING OF VECTOR IG IN CALLING PROGRAM */

																/*  *  NNU      -   NUMBER OF UNKNOWNS */

																/*  *  MATRIX   -   INTEGER VALUE CONTAINING INFO ABOUT THE MATRIX L. */

																/*                  1ST DIGIT OF MATRIX  --  ISYM: */
																/*                    =1: L IS SYMMETRIC; */
																/*                    =2: L IS NOT SYMMETRIC. */

																/*                  2ND DIGIT OF MATRIX  --  IROW0: */
																/*                    =1: L HAS ROWSUM ZERO; */
																/*                    =2: L DOES NOT HAVE ROWSUM ZERO. */

																/*         -------------------------------------------------------------- */

																/*     CLASS 2 - PARAMETERS: */

																/*     ISWTCH   -   PARAMETER CONTROLLING WHICH MODULES OF AMG1R5 ARE TO */
																/*                  BE USED. */
																/*                    =1:   CALL FOR -----, -----, -----, WRKCNT. */
																/*                    =2:   CALL FOR -----, -----, SOLVE, WRKCNT. */
																/*                    =3:   CALL FOR -----, FIRST, SOLVE, WRKCNT. */
																/*                    =4:   CALL FOR SETUP, FIRST, SOLVE, WRKCNT. */
																/*                  SETUP DEFINES THE OPERATORS NEEDED IN THE SOLUTION */
																/*                         PHASE. */
																/*                  FIRST INITIALIZES THE SOLUTION VECTOR (SEE PARAMETER */
																/*                         IFIRST). */
																/*                  SOLVE COMPUTES THE SOLUTION BY AMG CYCLING (SEE */
																/*                         PARAMETER NCYC). */
																/*                  WRKCNT PROVIDES THE USER WITH INFORMATION ABOUT */
																/*                         RESIDUALS, STORAGE REQUIREMENTS AND CP-TIMES */
																/*                         (SEE PARAMETER IOUT). */
																/*                  IF AMG1R5 IS CALLED THE FIRST TIME, ISWTCH HAS TO */
																/*                  BE =4. INDEPENDENT OF ISWTCH, SINGLE MODULES CAN BE */
																/*                  BYPASSED BY A PROPER CHOICE OF THE CORRESPONDING */
																/*                  PARAMETER. */

																/*     IOUT     -   PARAMETER CONTROLLING THE AMOUNT OF OUTPUT DURING */
																/*                  SOLUTION PHASE: */

																/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO. */

																/*                  2ND DIGIT: */
																/*                    =0: NO OUTPUT (EXCEPT FOR MESSAGES) */
																/*                    =1: RESIDUAL BEFORE AND AFTER SOLUTION PROCESS */
																/*                    =2: ADD.: STATISTICS ON CP-TIMES AND STORAGE REQUI- */
																/*                        REMENTS */
																/*                    =3: ADD.: RESIDUAL AFTER EACH AMG-CYCLE */

																/*     IPRINT   -   PARAMETER SPECIFYING THE FORTRAN UNIT NUMBERS FOR */
																/*                  OUTPUT: */

																/*                  1ST DIGIT: NOT USED; HAS TO BE NON-ZERO */

																/*                  2ND AND 3RD DIGIT  --  IUP: UNIT NUMBER FOR RESULTS */

																/*                  4TH AND 5TH DIGIT  --  IUM: UNIT NUMBER FOR MESSAGES */

																/*         -------------------------------------------------------------- */

																/*     CLASS 3 - PARAMETERS: */

																/*     LEVELX   -   MAXIMUM NUMBER OF MG-LEVELS TO BE CREATED (>=1). */

																/*     IFIRST   -   PARAMETER FOR FIRST APPROXIMATION. */

																/*                  1ST DIGIT OF IFIRST: NOT USED; HAS TO BE NON-ZERO. */

																/*                  2ND DIGIT OF IFIRST  --  ITYPU: */
																/*                    =0: NO SETTING OF FIRST APPROXIMATION, */
																/*                    =1: FIRST APPROXIMATION CONSTANT TO ZERO, */
																/*                    =2: FIRST APPROXIMATION CONSTANT TO ONE, */
																/*                    =3: FIRST APPROXIMATION IS RANDOM FUNCTION WITH */
																/*                        THE CONCRETE RANDOM SEQUENCE BEING DETERMINED */
																/*                        BY THE FOLLWING DIGITS. */

																/*                  REST OF IFIRST  --  RNDU: */
																/*                    DETERMINES THE CONCRETE RANDOM SEQUENCE USED IN */
																/*                    THE CASE ITYPU=3. (IFIRST=13 IS EQUIVALENT TO */
																/*                    IFIRST=1372815) */

																/*     NCYC     -   INTEGER PARAMETER DESCRIBING THE TYPE OF CYCLE TO BE */
																/*                  USED AND THE NUMBER OF CYCLES TO BE PERFORMED. */

																/*                  1ST DIGIT OF NCYC  --  IGAM: */
																/*                    =1: V -CYCLE, */
																/*                    =2: V*-CYCLE, */
																/*                    =3: F -CYCLE, */
																/*                    =4: W -CYCLE. */
																/*                  IF NCYC IS NEGATIV, THEN THE APPROXIMATION OF THE */
																/*                  PROBLEM ON THE SECOND FINEST GRID IS COMPUTED BY */
																/*                  IGAM V-CYCLES ON THAT PARTICULAR GRID. */

																/*                  2ND DIGIT OF NCYC  --  ICGR: */
																/*                    =0: NO CONJUGATE GRADIENT, */
																/*                    =1: CONJUGATE GRADIENT (ONLY FIRST STEP OF CG), */
																/*                    =2: CONJUGATE GRADIENT (FULL CG). */

																/*                  3RD DIGIT OF NCYC  --  ICONV: */
																/*                    CONVERGENCE CRITERION FOR THE USER-DEFINED PROBLEM */
																/*                    (FINEST GRID): */
																/*                    =1: PERFORM A FIXED NUMBER OF CYCLES AS GIVEN BY */
																/*                        NCYCLE (SEE BELOW) */
																/*                    =2: STOP, IF  ||RES|| < EPS */
																/*                    =3: STOP, IF  ||RES|| < EPS * |F| */
																/*                    =4: STOP, IF  ||RES|| < EPS * |U| * |DIAG| */
																/*                    WITH ||RES|| = L2-NORM OF RESIDUAL, */
																/*                           EPS     (SEE INPUT PARAMETER EPS) */
																/*                           |F|   = SUPREMUM NORM OF RIGHT HAND SIDE */
																/*                           |U|   = SUPREMUM NORM OF SOLUTION */
																/*                         |DIAG|  = MAXIMAL DIAGONAL ENTRY IN MATRIX L */
																/*                    NOTE THAT IN ANY CASE THE SOLUTION PROCESS STOPS */
																/*                    AFTER AT MOST NCYCLE CYCLES. */

																/*                  REST OF NCYC  --  NCYCLE: */
																/*                    MAXIMAL NUMBER OF CYCLES TO BE PERFORMED (>0) OR */
																/*                    NCYCLE=0: NO CYCLING. */

																/*     EPS      -   CONVERGENCE CRITERION FOR SOLUTION PROCESS: (SEE */
																/*                  PARAMETER NCYC). NOTE THAT NO MORE THAN NCYCLE CYCLES */
																/*                  ARE PERFORMED, REGARDLESS OF EPS. */

																/*     MADAPT   -   INTEGER VALUE SPECIFYING THE CHOICE OF COARSEST */
																/*                  GRID IN CYCLING: */

																/*                  1ST DIGIT OF MADAPT  --  MSEL: */
																/*                    =1: IN CYCLING, ALL GRIDS CONSTRUCTED IN THE SETUP */
																/*                        PHASE ARE USED WITHOUT CHECK. */
																/*                    =2: THE NUMBER OF GRIDS IS AUTOMATICALLY REDUCED */
																/*                        IF THE CONVERGENCE FACTOR ON THE COARSER GRIDS */
																/*                        IS FOUND TO BE LARGER THAN A GIVEN VALUE FAC */
																/*                        (SEE BELOW). */

																/*                  REST OF MADAPT  --  FAC */
																/*                        THE REST OF MADAPT DEFINES THE FRACTIONAL PART */
																/*                        OF A REAL NUMBER FAC BETWEEN 0.1 AND 0.99, E.G. */
																/*                        MADAPT=258 MEANS MSEL=2 AND FAC=0.58. IF MADAPT */
																/*                        CONSISTS OF ONLY ONE DIGIT, FAC IS SET TO 0.7 */
																/*                        BY DEFAULT. */


																/*     NRD      -   PARAMETER DESCRIBING RELAXATION (DOWNWARDS): */

																/*                  1ST DIGIT OF NRD: NOT USED; HAS TO BE NON-ZERO. */

																/*                  2ND DIGIT OF NRD  --  NRDX: */
																/*                    ACTUAL NUMBER OF SMOOTHING STEPS TO BE PERFORMED */
																/*                    THE TYPE OF WHICH IS GIVEN BY THE FOLLOWING DIGITS */

																/*                  FOLLOWING DIGITS  --  ARRAY NRDTYP: */
																/*                    =1: RELAXATION OVER THE F-POINTS ONLY */
																/*                    =2: FULL GS SWEEP */
																/*                    =3: RELAXATION OVER THE C-POINTS ONLY */
																/*                    =4: FULL MORE COLOR SWEEP, HIGHEST COLOR FIRST */

																/*     NSOLCO   -   PARAMETER CONTROLLING THE SOLUTION ON COARSEST GRID: */

																/*                  1ST DIGIT  --  NSC: */
																/*                    =1: GAUSS-SEIDEL METHOD */
																/*                    =2: DIRECT SOLVER (YALE SMP) */

																/*                  REST OF NSOLCO  --  NRCX: (ONLY IF NSC=1) */
																/*                  NUMBER OF GS SWEEPS ON COARSEST GRID (>=0). */
																/*                  IF NRCX=0, THEN AS MANY GS SWEEPS ARE PERFORMED */
																/*                  AS ARE NEEDED TO REDUCE THE RESIDUAL BY TWO ORDERS */
																/*                  OF MAGNITUDE. (MAXIMAL 100 RELAXATION SWEEPS) */

																/*     NRU      -   PARAMETER FOR RELAXATION (UPWARDS), ANALOGOUS TO NRD. */

																/*         -------------------------------------------------------------- */

																/*     CLASS 4 - PARAMETERS: */

																/*     ECG1,ECG2-   REAL PARAMETERS AFFECTING THE CREATION OF COARSER */
																/*     EWT2     -   GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
																/*                  THE CHOICE OF THESE PARAMETERS DEPENDS ON */
																/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

																/*     NWT      -   INTEGER PARAMETER AFFECTING THE CREATION OF COARSER */
																/*                  GRIDS AND/OR THE DEFINITION OF THE INTERPOLATION. */
																/*                  THE CHOICE OF THIS PARAMETER DEPENDS ON */
																/*                  THE ACTUAL AMG VERSION (SEE SUBROUTINE CRSNG) */

																/*     NTR      -   PARAMETER CONTROLLING COARSE-GRID OPERATOR TRUNCATION */
																/*                    =0: PAIRS OF ZEROES ARE REMOVED FROM COARSE GRID */
																/*                        OPERATORS */
																/*                    =1: NO COARSE-GRID OPERATOR TRUNCATION */


																/* ----------------------------------------------------------------------- */

																/*     OUTPUT: */

																/*     U        -   CONTAINS THE COMPUTED SOLUTION */


																/*     IERR     -   ERROR PARAMETER: */

																/*                    >0: FATAL ERROR (ABNORMAL TERMINATION OF AMG1R5) */
																/*                    <0: NON-FATAL ERROR (EXECUTION OF AMG1R5 CONTINUES) */

																/*                  ERROR CODES IN DETAIL: */

																/*                  1. DIMENSIONING TOO SMALL FOR VECTOR */
																/*                        A      (IERR = 1) */
																/*                        IA     (IERR = 2) */
																/*                        JA     (IERR = 3) */
																/*                        U      (IERR = 4) */
																/*                        F      (IERR = 5) */
																/*                        IG     (IERR = 6) */

																/*                     NO YALE-SMP BECAUSE OF STORAGE (NDA TOO SMALL): */
																/*                               (IERR = -1) */
																/*                     NO YALE-SMP BECAUSE OF STORAGE (NDJA TOO SMALL): */
																/*                               (IERR = -3) */
																/*                     NO CG BECAUSE OF STORAGE (NDU TOO SMALL): */
																/*                               (IERR = -4) */
																/*                     NO SPACE FOR TRANSPOSE OF INTERPOLATION (NDA OR */
																/*                                                     NDJA TOO SMALL): */
																/*                               (IERR = -1) */

																/*                  2. INPUT DATA ERRONEOUS: */

																/*                     A-ENTRY MISSING, ISYM = 1:           (IERR = -11) */
																/*                     PARAMETER MATRIX MAY BE ERRONEOUS:   (IERR = -12) */
																/*                     DIAGONAL ELEMENT NOT STORED FIRST:   (IERR =  13) */
																/*                     DIAGONAL ELEMENT NOT POSITIV:        (IERR =  14) */
																/*                     POINTER IA ERRONEOUS:                (IERR =  15) */
																/*                     POINTER JA ERRONEOUS:                (IERR =  16) */
																/*                     PARAMETER ISWTCH ERRONEOUS:          (IERR =  17) */
																/*                     PARAMETER LEVELX ERRONEOUS:          (IERR =  18) */

																/*                  3. ERRORS OF THE AMG1R5-SYSTEM (SHOULD NOT OCCUR): */

																/*                     TRANSPOSE A-ENTRY MISSING:           (IERR =  21) */
																/*                     INTERPOLATION ENTRY MISSING:         (IERR =  22) */

																/*                  4. ALGORITHMIC ERRORS: */

																/*                     CG-CORRECTION NOT DEFINED:           (IERR =  31) */
																/*                     NO YALE-SMP BECAUSE OF ERROR IN */
																/*                     FACTORIZATION:                       (IERR = -32) */

																/* ----------------------------------------------------------------------- */

																/*     WORK SPACE: */

																/*     THE INTEGER VECTOR IG HAS TO BE PASSED TO AMG1R5 AS WORK SPACE. */

																/* ----------------------------------------------------------------------- */

																/*     DIMENSIONING OF INPUT VECTORS AND WORK SPACE: */

																/*     IT'S IMPOSSIBLE TO TELL IN ADVANCE THE EXACT STORAGE REQUIREMENTS */
																/*     OF AMG. THUS, THE FOLLOWING FORMULAS GIVE ONLY REASONABLE GUESSES */
																/*     FOR THE VECTOR LENGTHS WHICH HAVE TO BE DECLARED IN THE CALLING */
																/*     PROGRAM. IN THESE FORMULAS NNA DENOTES THE NUMBER OF NON-ZERO */
																/*     ENTRIES IN THE INPUT-MATRIX L AND NNU IS THE NUMBER OF UNKNOWNS. */

																/*     VECTOR         NEEDED LENGTH (GUESS) */
																/*       A               3*NNA + 5*NNU */
																/*       JA              3*NNA + 5*NNU */
																/*       IA              2.2*NNU */
																/*       U               2.2*NNU */
																/*       F               2.2*NNU */
																/*       IG              5.4*NNU */

																/* ----------------------------------------------------------------------- */


																/*     STANDARD CHOICES OF PARAMETERS (AS FAR AS MEANINGFUL): */

																/*          ISWTCH = 4 */
																/*          IOUT   = 12 */
																/*          IPRINT = 10606 */

																/*          LEVELX = 25 */
																/*          IFIRST = 13 */
																/*          NCYC   = 10110 */
																/*          EPS    = 1.D-12 */
																/*          MADAPT = 27 */
																/*          NRD    = 1131 */
																/*          NSOLCO = 110 */
																/*          NRU    = 1131 */

																/*          ECG1   = 0. */
																/*          ECG2   = 0.25 */
																/*          EWT2   = 0.35 */
																/*          NWT    = 2 */
																/*          NTR    = 0 */

																/*     IF ANY ONE OF THESE PARAMETERS IS 0 ON INPUT, ITS CORRESPONDING */
																/*     STANDARD VALUE IS USED BY AMG1R5. */

																/* ----------------------------------------------------------------------- */

																/*     PORTABILITY RESTRICTIONS: */

																/*     1. ROUTINE CTIME IS MACHINE DEPENDENT AND HAS TO BE ADAPTED TO */
																/*        YOUR COMPUTER INSTALLATION OR REPLACED BY A DUMMY ROUTINE. */

																/*     2. MOST INPUT PARAMETERS ARE COMPOSED OF SEVERAL DIGITS, THEIR */
																/*        SIGNIFICANCE HAVING BEEN DESCRIBED ABOVE. BE SURE NOT TO ENTER */
																/*        MORE DIGITS THAN YOUR COMPUTER CAN STORE ON AN INTEGER VARI- */
																/*        ABLE. */

																/*     3. APART FROM FORTRAN INTRINSIC FUNCTIONS AND SERVICE ROUTINES, */
																/*        THERE IS ONLY ONE EXTERNAL REFERENCE TO A PROGRAM NOT CONTAINED */
																/*        IN THE AMG1R5 - SYSTEM, I.E. THE LINEAR SYSTEM SOLVER NDRV OF */
																/*        THE YALE SPARSE MATRIX PACKAGE. IF YOU HAVN'T ACCESS TO THIS */
																/*        PACKAGE, ENTER A DUMMY ROUTINE NDRV AND AVOID CHOOSING NSC=2 */
																/*        (SUBPARAMETER OF NSOLCO). THEN NDRV ISN'T CALLED BY AMG1R5. */
																/*        IN THIS CASE, HOWEVER, INDEFINITE PROBLEMS WILL NOT BE SOLV- */
																/*        ABLE. */
																/*          THE YALE SPARSE MATRIX PACKAGE IS FREELY AVAILABLE FOR NON- */
																/*        PROFIT PURPOSES. CONTACT THE DEPARTMENT OF COMPUTER SCIENCE, */
																/*        YALE UNITVERSITY. */

																/*     4. IN AMG1R5 THERE IS THE PARAMETER LRATIO, DENOTING THE RATIO */
																/*        OF SPACE OCCUPIED BY A DOUBLE PRECISION REAL VARIABLE AND */
																/*        THAT OF AN INTEGER. FOR THE IBM-VERSION LRATIO HAS BEEN SET */
																/*        TO 2. CHANGE THIS VALUE IF NECESSARY. (THE SAME HAS TO BE */
																/*        DONE WITH THE YALE SMP-ROUTINE NDRV.) */


																/* ----------------------------------------------------------------------- */

																/*     AUTHORS: */

																/*          JOHN RUGE, FORT COLLINS (USA), */
																/*              INSTITUTE FOR COMPUTATIONAL STUDIES AT CSU; */

																/*          KLAUS STUEBEN, D-5205 ST. AUGUSTIN (W.-GERMANY), */
																/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

																/*          ROLF HEMPEL, D-5205 ST. AUGUSTIN (W.-GERMANY), */
																/*              GESELLSCHAFT FUER MATHEMATIK UND DATENVERARBEITUNG (GMD). */

																/* ----------------------------------------------------------------------- */


																/* ===> LRATIO HAS TO BE SET TO THE NUMBER OF INTEGERS OCCUPYING THE SAME */
																/*     AMOUNT OF STORAGE AS ONE DOUBLE PRECISION REAL. */


																/* ===> MAXGR IS THE MAXIMAL NUMBER OF GRIDS. CHANGING THIS UPPER LIMIT */
																/*     JUST REQUIRES CHANGING THE PARAMETER STATEMENT. */


																/* Parameter adjustments */
	--ig;
	--f;
	--u;
	--ja;
	--ia;
	--a;

	/* Function Body */
	*ierr = 0;

	/* ===> SET PARAMETERS TO STANDARD VALUES, IF NECCESSARY */


	if (*iout != 0) {
		idec_(iout, &c__2, &ndigit, iarr);
		kout = iarr[1];
	}
	else {
		kout = 2;
	}

	if (*iswtch != 0) {
		kswtch = *iswtch;
	}
	else {
		kswtch = 4;
	}

	if (*levelx > 0) {
		kevelx = my_imin(*levelx, 25);
	}
	else if (*levelx < 0) {
		goto L70;
	}
	else {
		kevelx = 25;
	}

	if (*iprint != 0) {
		idec_(iprint, &c__4, &ndigit, iarr);
		iup = iarr[1] * 10 + iarr[2];
		ium = iarr[3];
	}
	else {
		iup = 6;
		ium = 6;
	}
	icgst = *nnu + 3;
	ndicg = (*ndig - icgst + 1) / 2;
	if (ndicg <= 0) {
		goto L60;
	}

	switch (kswtch) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	}
	//io___10.ciunit = ium;
	//s_wsfe(&io___10);
	//e_wsfe();


	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER I.  SWTCH ***\n");

	*ierr = 17;
	return 0;

L40:
	setup_(nnu, matrix, &kevelx, ecg1, ecg2, ewt2, nwt, ntr, ierr, &a[1], &u[
		1], &ia[1], &ja[1], &ig[1], imin, imax, iminw, imaxw, &ig[icgst],
			&ig[icgst + ndicg], nstcol, iarr, time, &levels, &irow0, nda,
			ndia, ndja, ndu, ndf, &ndicg, &ium, &mda, &mdia, &mdja, &mdu, &
			mdf, &mdig, &c__25, &c__2);
	if (*ierr > 0) {
		return 0;
	}
L30:
	first_(ifirst, &u[1], imin, imax, iarr, &irow0);
L20:

	// Алгебраический Многосеточный Метод как предобуславливатель
	// к алгоритму Крыловского типа Ю.Саала и Шульца FGMRES [1986]
	// FGMRES - гибкий вариант обобщенного метода минимальных невязок
	// в котором в качестве предобуславливателя используется amg1r5 алгоритм.
	// Требует ещё одну память под матрицу А на самом подробном уровне.
	// 5.01.2017 Алгоритм FGMRES изобретён в 1986 году.

	/*
	// Для корректности работы надо передавать информацию о модифицированном размере в вызывающие функции вверх по коду.
	if ((a != NULL) && (ja != NULL)) {

		// На данном этапе доступен истинный размер матрицы СЛАУ - хранимое число ненулевых элементов.
		// Вычислим реальное число ненулевых элементов матрицы a. Если номер столбца равен нулю то число ненулевых элементов
		// в матрице заведомо меньше чем позиция этого нуля.
		// С помощью низкоуровневой операции realloc ужимаем размер матрицы a.
		// Последующие выделения оперативной памяти для внешнего Крыловского итерационного процесса BiCGStab не приведут к ещё большему расходу
		// оперативной памяти, т.к. мы её освободили.

		printf("nda=%lld\n", nda[0]);
		integer isize97 = 0;
		for (integer i_95 = 1; i_95 < ndja[0]; i_95++) {
			if (ja[i_95] == 0) {
				isize97 = i_95;
				if (i_95 + 2 < ndja[0] && ja[i_95 + 1] == 0 && ja[i_95 + 2] == 0) {
					break;
				}
			}
		}
		printf("nda_new=%lld\n", isize97);
		++a;
		++ja;
		a = (doublereal*)realloc(a, ((integer)(isize97)+2) * sizeof(doublereal));
		ja = (integer*)realloc(ja, ((integer)(isize97)+2) * sizeof(integer));
		--a;
		--ja;
		//getchar();
	}
	*/

	ncyc[0] = 1011; // Предобуславливание в один V цикл.

					/*
					// для отладки.
					for (integer i_numberV_cycle = 0; i_numberV_cycle < 10; i_numberV_cycle++) {
					ncyc[0] = 1011;
					solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &u[1], &f[1], &
					ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
					&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
					nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
					res0, &res);
					getchar();
					}
					*/
	printf("sizeof  ndu=%lld nnu=%lld ndf=%lld\n", ndu[0], nnu[0], ndf[0]);
	
	/*
	integer n75 = -1 , nnz;

	for (integer i_83 = 1; i_83 <= nda[0]; i_83++) {
	if (ja[i_83] > n75) n75 = ja[i_83];
	}
	nnz = ia[n75 + 1];
	printf("n_75=%d nnz=%d\n",n75,nnz);
	getchar();
	*/
	

	// Разреженная матрица СЛАУ
	// в CRS формате.


	if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
		/*
		if (ibackregulationgl != NULL) {
		// nested desection версия алгоритма.
		integer ierr = equation3DtoCRSnd(sl, slb, val75, col_ind75, row_ptr75, maxelm, maxbound, 1.0, true, NULL, NULL);
		if (ierr > 0) {
		switch (iVar) {
		case VX: printf("VX equation problem.\n"); break;
		case VY: printf("VY equation problem.\n"); break;
		case VZ: printf("VZ equation problem.\n"); break;
		case PAM: printf("PAM equation problem.\n"); break;
		}
		}
		}
		else {
		*/
		integer ierr = equation3DtoCRS(sl, slb, val75, col_ind75, row_ptr75, maxelm, maxbound, 1.0, true);
		if (ierr > 0) {
			switch (iVar) {
			case VX: printf("VX equation problem.\n"); break;
			case VY: printf("VY equation problem.\n"); break;
			case VZ: printf("VZ equation problem.\n"); break;
			case PAM: printf("PAM equation problem.\n"); break;
			}
		}
		//}
	}
	if (iVar == TEMP) {
		integer ierr = equation3DtoCRS(sl, slb, val75, col_ind75, row_ptr75, maxelm, maxbound, 1.0, true);
		if (ierr > 0) {
			printf("Temperature equation problem.\n");
		}
	}

	if ((val75 == NULL) || (col_ind75 == NULL) || (row_ptr75 == NULL)) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for val, col_ind or row_ptr: bicgStab + camg...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	nnz = row_ptr75[n75];
	printf("n=%lld nnz=%lld ndu=%lld nda=%lld\n", n75, nnz, ndu[0], nda[0]);



	/*
	// Так делать нельзя, т.к. матрица СЛАУ (a,ia,ja)
	// была модифицирована в процессе работы amg1r5.f алгоритма.

	// инициализация матрицы.
	#pragma omp parallel for
	for (integer i_91 = 1; i_91 <= n75; i_91++) {

	for (integer i_92 = ia[i_91]; i_92 < ia[i_91 + 1]; i_92++) {

	val75[i_92 - 1] = a[i_92];
	col_ind75[i_92 - 1] = ja[i_92]- 1;

	}
	row_ptr75[i_91 - 1] = ia[i_91] - 1;
	}
	row_ptr75[n75] = ia[n75 + 1] - 1;
	*/

	//system("PAUSE");

	
	//Vector s(m + 1), cs(m + 1), sn(m + 1), w;
	w = new doublereal[n75];
	s = new doublereal[m_restart + 2];
	cs = new doublereal[m_restart + 2];
	sn = new doublereal[m_restart + 2];

	dx = new doublereal[n75];
	buffer = new doublereal[n75];
	Zcopy = new doublereal[ndu[0] +1];
	vCopy = new doublereal[ndu[0] +1];

	

	--Zcopy;
	--vCopy;

	// начальное приближение
	for (i = 0; i<n75; i++) dx[i] = u[i + 1];


	//doublereal normb = norm(M.solve(b));
	
	// здесь реализованы все три нормы
	// вообще говоря они все эквивалентны



	normb = NormaV_for_gmres(&f[1], n75);
	//normb = NormaV(buffer, n75);

	//Vector r = &f[1] - A * x;
	

	r = new doublereal[n75];
	MatrixCRSByVector(val75, col_ind75, row_ptr75, dx, r, n75); // результат занесён в  r
	for (i = 0; i < n75; i++) r[i] = f[i + 1] - r[i];

	//  calculate residual precontidioning;


	//doublereal beta = norm(r);
	


	beta = NormaV_for_gmres(r, n75);

	if (fabs(normb) < 1.0e-30)
		normb = 1;

	

	norm_r = NormaV_for_gmres(r, n75);

	


	if ((resid = norm_r / normb) <= dterminatedTResudual) {
		//tol = resid;
		maxit = 0;
		return 0;
	}

	

	H = new doublereal*[m_restart + 2]; // Hessenberg
	for (i_1 = 0; i_1 < m_restart + 2; i_1++) H[i_1] = new doublereal[m_restart + 2];


	for (i_1 = 0; i_1 < m_restart + 2; i_1++)
	{
		for (integer j_1 = 0; j_1 < m_restart + 2; j_1++)
		{
			H[i_1][j_1] = 0.0;
		}
	}

	

	//Vector *v = new Vector[m_restart + 1];
	
		v = new doublereal*[m_restart + 2];
	for (i_1 = 0; i_1 <= m_restart + 1; i_1++) v[i_1] = new doublereal[n75];


	for (i_1 = 0; i_1 <= m_restart + 1; i_1++) {
		for (integer j_1 = 0; j_1 < n75; j_1++)
		{
			v[i_1][j_1] = 0.0;
		}
	}

	 Z = new doublereal*[m_restart + 2];
	 
	for (i_1 = 0; i_1 <= m_restart + 1; i_1++) Z[i_1] = new doublereal[n75];

	for (i_1 = 0; i_1 <= m_restart + 1; i_1++) {
		for (integer j_1 = 0; j_1 < n75; j_1++)
		{
			Z[i_1][j_1] = 0.0;
		}
	}

	j = 1; // номер первой итерации
		   //doublereal delta = 1.0e-3;// DOPOLNENIE

	

	while (j <= maxit) {

		//v[0] = r * (1.0 / beta);    // ??? r / beta
		for (integer j_1 = 0; j_1 < n75; j_1++)
		{
			v[0][j_1] = r[j_1] * (1.0 / beta);
		}

		//s = 0.0;
		for (i_1 = 0; i_1 <= m_restart + 1; i_1++) s[i_1] = 0.0;
		s[0] = beta;
		//s[0] = 1.0;


		for (i_1 = 0; i_1 < m_restart + 2; i_1++)
		{ // DOPOLNENIE
			for (integer j_1 = 0; j_1 < m_restart + 2; j_1++)
			{
				H[i_1][j_1] = 0.0;
			}
		}


		// Ортогонализация Арнольди.
		for (i = 0; i < m_restart && j <= maxit; i++, j++) {

			i_copy = i;


			// KZ[i]=v[i]

			// (LU)Z[i]=v[i];

			// multigrid Ruge and Stuben preconditioning [1986].
			// достаточно одного V цикла.
			// K*Z = v;
			for (i_1 = 0; i_1 < n75; i_1++) {
				Zcopy[i_1+1] = 0.0;
				vCopy[i_1+1] = v[i][i_1];
			}
			
			ifirst[0] = 11; // Нулевое начальное приближение
			for (integer i_numberV_cycle = 0; i_numberV_cycle < 1; i_numberV_cycle++) {
				ncyc[0] = 1011; // достаточно одного V цикла.
				solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &Zcopy[1], &vCopy[1], &
					ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
					&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
					nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
					res0, &res);
				//getchar();
			}
		
			for (i_1 = 0; i_1 < n75; i_1++) {
				Z[i][i_1] = Zcopy[i_1+1];
			}

			/*
			if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
			// Очень важно начинать с нуля иначе не будет сходимости.
			#pragma omp parallel for shared(m) private(i_1) schedule (guided)
			for (i_1 = 0; i_1<n75; i_1++) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

			//  9 августа 2015 при внедрении перенумерации узлов nested desection
			if (bpam_gsp && (iVar == PAM)) {
			if (ibackregulationgl != NULL) {
			PAMGSPnd(sl, slb, m.y, v[i], maxelm, maxbound, ifrontregulationgl);
			}
			else {
			PAMGSP(sl, slb, m.y, v[i], maxelm, maxbound);
			}
			}
			else {

			if (brc) {
			for (integer i7 = 0; i7<n75; i7++) m.vec[i7] = m.pi[i7];
			for (integer i7 = 0; i7<m.iwk + 2; i7++) {
			m.alurc[i7] = m.alu[i7];
			m.jlurc[i7] = m.jlu[i7];
			}
			for (integer i7 = 0; i7<n75 + 2; i7++) m.jurc[i7] = m.ju[i7];
			}

			if (ibackregulationgl != NULL) {
			//lusol_2(n75, v[i], m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=v[i];
			lusol_3(n75, v[i], m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=v[i];
			}
			else {
			lusol_(n75, v[i], m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=v[i];

			}

			if (brc) {
			for (integer i7 = 0; i7<n75; i7++) m.pi[i7] = m.vec[i7];
			for (integer i7 = 0; i7<m.iwk + 2; i7++) {
			m.alu[i7] = m.alurc[i7];
			m.jlu[i7] = m.jlurc[i7];
			}
			for (integer i7 = 0; i7<n75 + 2; i7++) m.ju[i7] = m.jurc[i7];
			}

			}
			//for (i_1 = 0; i_1 < n75; i_1++) w[i_1] = m.y[i_1];
			for (i_1 = 0; i_1 < n75; i_1++) Z[i][i_1] = m.y[i_1];
			//for (i_1 = 0; i_1 < n75; i_1++)  v[i + 1][i_1] = m.y[i_1];


			}


			if (iVar == TEMP) {
			// Очень важно начинать с нуля иначе не будет сходимости.
			//#pragma omp parallel for shared(m) private(i) schedule (guided)
			for (i_1 = 0; i_1<n75; i_1++) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

			lusol_(n, v[i], m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=v[i];
			for (i_1 = 0; i_1 < n75; i_1++) Z[i][i_1] = m.ty[i_1];
			//for (i_1 = 0; i_1 < n75; i_1++) v[i + 1][i_1] = m.ty[i_1];

			}
			*/

			// Совсем без предобуславливателя.
			//for (i_1 = 0; i_1 < n75; i_1++) Z[i][i_1] = v[i][i_1];

			// Закоментировано без предобуславливания.
			//w = A * Z[i];
			MatrixCRSByVector(val75, col_ind75, row_ptr75, Z[i], w, n75); // результат занесён в  w

			for (k = 0; k <= i; k++) {
				H[k][i] = Scal(w, v[k], n75);

				for (integer j_1 = 0; j_1 < n75; j_1++)
				{
					w[j_1] -= H[k][i] * v[k][j_1];
				}
			}
			H[i + 1][i] = NormaV_for_gmres(w, n75);



			for (integer j_1 = 0; j_1 < n75; j_1++)
			{
				v[i + 1][j_1] = w[j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			}
			// Окончание ортогонализации Арнольди.
			// В v - хранится ортонормированный базис подпространства Крылова размерности m_restart.
			// H - Верхнетреугольная матрица Хессенберга - матрица коэффициентов ортогонализации.


			// 26.11.2017
			// Это проверенный и испытанный кусок кода.
			for (k = 0; k < i; k++)
				ApplyPlaneRotation(H[k][i], H[k + 1][i], cs[k], sn[k]);

			GeneratePlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);



			// Вручную устраняем случай полного совпадения невязок на двух соседних итерациях,
			// т.к. иначе это приводит к развалу решения.
			//if (fabs(s[i] - s[i + 1]) < 1.0e-37) s[i + 1] = 1.05*s[i];

			printf("%d %e \n", j, fabs(s[i + 1]) / normb);
			//printf("%d %e \n", j, beta*fabs(s[i + 1]));
			//getchar();

			resid = fabs(s[i + 1]) / normb;
			//resid = beta*fabs(s[i + 1]);

			if ((resid) < dterminatedTResudual) {
				printf("dosrochnji vjhod\n");
				bOkfgmres_amg1r5 = true;
				//getchar();				
				Update(dx, i, n75, H, s, Z);
				//tol = resid;
				//maxit = j;

				for (i_1 = 0; i_1<n75; i_1++) {
					u[i_1 + 1] = dx[i_1];
				}

				for (i_1 = 0; i_1 <= m_restart + 1; i_1++) delete[] v[i_1];
				delete[] v;
				for (i_1 = 0; i_1 <= m_restart + 1; i_1++) delete[] Z[i_1];
				delete[] Z;
				for (integer i_1 = 0; i_1 < m_restart + 2; i_1++) delete[] H[i_1];
				delete[] H;
				delete[] dx;
				delete[] buffer;
				delete[] r;
				delete[] w;
				delete[] s;
				delete[] cs;
				delete[] sn;
				++Zcopy;
				++vCopy;
				delete[] Zcopy;
				delete[] vCopy;

				// Освобождение оперативной памяти.
				if (val75 != NULL) {
					delete[] val75;
					val75 = NULL;
				}
				if (col_ind75 != NULL) {
					delete[] col_ind75;
					col_ind75 = NULL;
				}
				if (row_ptr75 != NULL) {
					delete[] row_ptr75;
					row_ptr75 = NULL;
				}

				goto LABEL_FGMRES_CONTINUE;

			}
		}



		// i-1 -> m_restart-1
		Update(dx, i - 1, n75, H, s, Z);//i-1 //ERROR

										//r = M.solve(b - A * x);
		MatrixCRSByVector(val75, col_ind75, row_ptr75, dx, r, n75); // Результат занесён в r
		for (i_1 = 0; i_1 < n75; i_1++) r[i_1] = f[i_1 + 1] - r[i_1];

		//beta = norm(r);
		beta = NormaV_for_gmres(r, n75);

		resid = beta / normb;
		//resid = beta;

		if ((resid) < dterminatedTResudual) {
			//tol = resid;
			//maxit = j;

			printf("end\n");
			//getchar();

			for (i_1 = 0; i_1<n75; i_1++) {
				u[i_1 + 1] = dx[i_1];

			}

			for (i_1 = 0; i_1 <= m_restart + 1; i_1++) delete[] v[i_1];
			delete[] v;
			for (i_1 = 0; i_1 <= m_restart + 1; i_1++) delete[] Z[i_1];
			delete[] Z;
			for (integer i_1 = 0; i_1 < m_restart + 2; i_1++) delete[] H[i_1];
			delete[] H;
			delete[] dx;
			delete[] buffer;
			delete[] r;
			delete[] w;
			delete[] s;
			delete[] cs;
			delete[] sn;
			++Zcopy;
			++vCopy;
			delete[] Zcopy;
			delete[] vCopy;

			// Освобождение оперативной памяти.
			if (val75 != NULL) {
				delete[] val75;
				val75 = NULL;
			}
			if (col_ind75 != NULL) {
				delete[] col_ind75;
				col_ind75 = NULL;
			}
			if (row_ptr75 != NULL) {
				delete[] row_ptr75;
				row_ptr75 = NULL;
			}

			goto LABEL_FGMRES_CONTINUE;


		}
	}

	//tol = resid;
	for (i_1 = 0; i_1<n75; i_1++) {
		u[i_1 + 1] = dx[i_1];

	}

	for (i_1 = 0; i_1 <= m_restart + 1; i_1++) delete[] v[i_1];
	delete[] v;
	for (i_1 = 0; i_1 <= m_restart + 1; i_1++) delete[] Z[i_1];
	delete[] Z;
	for (integer i_1 = 0; i_1 < m_restart + 2; i_1++) delete[] H[i_1];
	delete[] H;
	delete[] dx;
	delete[] buffer;
	delete[] r;
	delete[] w;
	delete[] s;
	delete[] cs;
	delete[] sn;
	++Zcopy;
	++vCopy;
	delete[] Zcopy;
	delete[] vCopy;

	// Освобождение оперативной памяти.
	if (val75 != NULL) {
		delete[] val75;
		val75 = NULL;
	}
	if (col_ind75 != NULL) {
		delete[] col_ind75;
		col_ind75 = NULL;
	}
	if (row_ptr75 != NULL) {
		delete[] row_ptr75;
		row_ptr75 = NULL;
	}
	goto LABEL_FGMRES_CONTINUE;


LABEL_FGMRES_CONTINUE:

	/*
	solve_(madapt, ncyc, nrd, nsolco, nru, &kout, ierr, &a[1], &u[1], &f[1], &
	ia[1], &ja[1], &ig[1], eps, imin, imax, iminw, imaxw, &ig[icgst],
	&ig[icgst + ndicg], nstcol, iarr, time, &ncyc0, &irow0, &levels,
	nda, ndja, ndu, ndf, &mda, &mdja, &mdu, &mdf, &iup, &ium, resi, &
	res0, &res);
	*/
	if (*ierr > 0) {
		return 0;
	}
L10:
	wrkcnt_(&kout, &ia[1], &ig[1], imin, imax, iminw, &levels, time, &ncyc0, &
		iup, &mda, &mdia, &mdja, &mdu, &mdf, &mdig, &res0, &res);
	return 0;
L60:
	//io___29.ciunit = ium;
	//s_wsfe(&io___29);
	//e_wsfe();
	printf("*** ERROR IN AMG1R5: NDIG TOO SMALL ***\n");

	*ierr = 6;
	return 0;
L70:
	//io___30.ciunit = ium;
	//s_wsfe(&io___30);
	//e_wsfe();

	printf("*** ERROR IN AMG1R5: ILLEGAL PARAMETER LEVELX ***\n");

	*ierr = 18;
	return 0;
} /* amg1r5_fgmres_version */

  

#endif

//CUSP 0.5.1 library
// my_cusp_alg.cu

#pragma once
//#define GPU_LIB_INCLUDE_MY_PROJECT 0

// amgx
//#include "amgx_c.h"

/*
#include <cusp/csr_matrix.h>
#include <cusp/krylov/bicgstab.h>

#if defined(__cplusplus)
extern "C" {
#endif

void bicgstab_(integer * device_I, integer * device_J, float * device_V, float * device_x, float * device_b, integer N, integer NNZ){

// *NOTE* raw pointers must be wrapped with thrust::device_ptr!
thrust::device_ptr<int> wrapped_device_I(device_I);
thrust::device_ptr<int> wrapped_device_J(device_J);
thrust::device_ptr<float> wrapped_device_V(device_V);
thrust::device_ptr<float> wrapped_device_x(device_x);
thrust::device_ptr<float> wrapped_device_b(device_b);

// use array1d_view to wrap the individual arrays
typedef typename cusp::array1d_view< thrust::device_ptr<int> > DeviceIndexArrayView;
typedef typename cusp::array1d_view< thrust::device_ptr<float> > DeviceValueArrayView;

DeviceIndexArrayView row_indices (wrapped_device_I, wrapped_device_I + (N+1));
DeviceIndexArrayView column_indices(wrapped_device_J, wrapped_device_J + NNZ);
DeviceValueArrayView values (wrapped_device_V, wrapped_device_V + NNZ);
DeviceValueArrayView x (wrapped_device_x, wrapped_device_x + N);
DeviceValueArrayView b (wrapped_device_b, wrapped_device_b + N);

// combine the three array1d_views into a csr_matrix_view
typedef cusp::csr_matrix_view<DeviceIndexArrayView,
DeviceIndexArrayView,
DeviceValueArrayView> DeviceView;

// construct a csr_matrix_view from the array1d_views
DeviceView A(N, N, NNZ, row_indices, column_indices, values);

// set stopping criteria:
// iteration_limit = 100
// relative_tolerance = 1e-5
cusp::verbose_monitor<float> monitor(b, 100, 1e-5);

// solve the linear system A * x = b with the Conjugate Gradient method
cusp::krylov::bicgstab(A, x, b, monitor);

}

#if defined(__cplusplus)
}
#endif

*/

#if GPU_LIB_INCLUDE_MY_PROJECT == 1 

// CUSP 0.5.1 NVIDIA Includes

//#include "cusp_library\cusp\csr_matrix.h"
//#include "cusp_library\cusp/krylov/bicgstab.h"
//#include  <cusp\csr_matrix.h>
//#include  <cusp/krylov/bicgstab.h>
//#include <cusp/hyb_matrix.h>
#include <cusp/coo_matrix.h>
#include <cusp/csr_matrix.h>
//#include <cusp/monitor.h>
//-->//#include <cusp/array1d.h>
//#include <cusp/io/matrix_market.h>
//#include <cusp/krylov/cg.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/precond/ainv.h>
#include <cusp/krylov/bicgstab.h>
//---->//#include <cusp/krylov/gmres.h>
//#include <cusp/precond/aggregation/smoothed_aggregation_options.h>
//#include <cusp/csr_matrix.h>
//#include <cusp/blas/blas.h>
//#include <cusp/linear_operator.h>
//#include <cusp/gallery/poisson.h>
//#include <cusp/relaxation/gauss_seidel.h>


#include <iostream>

template <typename Monitor>
void report_status(Monitor& monitor)
{
	if (monitor.converged())
	{
		std::cout << "  Solver converged to " << monitor.tolerance() << " tolerance";
		std::cout << " after " << monitor.iteration_count() << " iterations";
		std::cout << " (" << monitor.residual_norm() << " final residual)" << std::endl;
	}
	else
	{
		std::cout << "  Solver reached iteration limit " << monitor.iteration_limit() << " before converging";
		std::cout << " to " << monitor.tolerance() << " tolerance ";
		std::cout << " (" << monitor.residual_norm() << " final residual)" << std::endl;
	}
}

#endif


#if GPU_LIB_INCLUDE_MY_PROJECT==1

//bool bgl_first_start_nonlinear_cusp_amg = true;

//bool bstart7 = true;

#if 1
// Это вызов библиотечного решателя систем линейных алгебраических уравнений.
// Библиотека Cusp версии 0.5.1. Это библиотека с открытым исходным кодом распространяющаяся 
// по open Source лицензии Apache license 2.0. 
// На хосте (центральном процессоре) подключён метод bicgstab и алгебраический многосеточный метод
// на основе сглаженной аггрегации samg.
// Дата подсоединения 12 октября 2016 года.
void cusp_solver_amghost(equation3D* &sl, equation3D_bon* &slb,
	integer maxelm, integer maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit,
	doublereal alpharelax, integer iVar)
{

	// maxit,  iVar - не используется.


	if (dX0 == NULL) {
		dX0 = new doublereal[maxelm + maxbound];
		for (integer i = 0; i < maxelm + maxbound; i++) {
			dX0[i] = 0.0;
		}
	}

	// TODO получить val, col_ind, row_ptr
	integer nna = 0; // количество ненулевых элементов в матрице СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	// подсчёт числа ненулевых элементов в матрице.
	nna = 0;
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB2>-1) && (fabs(sl[i].ab2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE2>-1) && (fabs(sl[i].ae2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN2>-1) && (fabs(sl[i].an2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS2>-1) && (fabs(sl[i].as2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT2>-1) && (fabs(sl[i].at2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW2>-1) && (fabs(sl[i].aw2) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB3>-1) && (fabs(sl[i].ab3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE3>-1) && (fabs(sl[i].ae3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN3>-1) && (fabs(sl[i].an3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS3>-1) && (fabs(sl[i].as3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT3>-1) && (fabs(sl[i].at3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW3>-1) && (fabs(sl[i].aw3) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB4>-1) && (fabs(sl[i].ab4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE4>-1) && (fabs(sl[i].ae4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN4>-1) && (fabs(sl[i].an4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS4>-1) && (fabs(sl[i].as4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT4>-1) && (fabs(sl[i].at4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW4>-1) && (fabs(sl[i].aw4) > nonzeroEPS)) (nna)++;

	}
	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	}

	integer nnu = 0; // число неизвестных.
	nnu = maxelm + maxbound;

	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
	//typedef float    ScalarType;
	typedef cusp::device_memory MemorySpace;


	// allocate storage for (nnu,nnu) matrix with nna nonzeros
	cusp::coo_matrix<int, ScalarType, cusp::host_memory> Ah(nnu, nnu, nna);
	//cusp::coo_matrix<int, float, cusp::device_memory> Ah(nnu, nnu, nna);

	//printf("0\n");
	//getchar();

	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iP; Ah.values[nna] = sl[i].ap / alpharelax;
			(nna)++;
		}
		if ((sl[i].iB > -1) && (fabs(sl[i].ab) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB; Ah.values[nna] = -sl[i].ab;
			(nna)++;
		}
		if ((sl[i].iE > -1) && (fabs(sl[i].ae) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE; Ah.values[nna] = -sl[i].ae;
			(nna)++;
		}
		if ((sl[i].iN > -1) && (fabs(sl[i].an) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN; Ah.values[nna] = -sl[i].an;
			(nna)++;
		}
		if ((sl[i].iS > -1) && (fabs(sl[i].as) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS; Ah.values[nna] = -sl[i].as;
			(nna)++;
		}
		if ((sl[i].iT > -1) && (fabs(sl[i].at) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT; Ah.values[nna] = -sl[i].at;
			(nna)++;
		}
		if ((sl[i].iW > -1) && (fabs(sl[i].aw) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW; Ah.values[nna] = -sl[i].aw;
			(nna)++;
		}

		if ((sl[i].iB2 > -1) && (fabs(sl[i].ab2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB2; Ah.values[nna] = -sl[i].ab2;
			(nna)++;
		}
		if ((sl[i].iE2 > -1) && (fabs(sl[i].ae2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE2; Ah.values[nna] = -sl[i].ae2;
			(nna)++;
		}
		if ((sl[i].iN2 > -1) && (fabs(sl[i].an2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN2; Ah.values[nna] = -sl[i].an2;
			(nna)++;
		}
		if ((sl[i].iS2 > -1) && (fabs(sl[i].as2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS2; Ah.values[nna] = -sl[i].as2;
			(nna)++;
		}
		if ((sl[i].iT2 > -1) && (fabs(sl[i].at2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT2; Ah.values[nna] = -sl[i].at2;
			(nna)++;
		}
		if ((sl[i].iW2 > -1) && (fabs(sl[i].aw2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW2; Ah.values[nna] = -sl[i].aw2;
			(nna)++;
		}

		if ((sl[i].iB3 > -1) && (fabs(sl[i].ab3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB3; Ah.values[nna] = -sl[i].ab3;
			(nna)++;
		}
		if ((sl[i].iE3 > -1) && (fabs(sl[i].ae3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE3; Ah.values[nna] = -sl[i].ae3;
			(nna)++;
		}
		if ((sl[i].iN3 > -1) && (fabs(sl[i].an3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN3; Ah.values[nna] = -sl[i].an3;
			(nna)++;
		}
		if ((sl[i].iS3 > -1) && (fabs(sl[i].as3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS3; Ah.values[nna] = -sl[i].as3;
			(nna)++;
		}
		if ((sl[i].iT3 > -1) && (fabs(sl[i].at3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT3; Ah.values[nna] = -sl[i].at3;
			(nna)++;
		}
		if ((sl[i].iW3 > -1) && (fabs(sl[i].aw3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW3; Ah.values[nna] = -sl[i].aw3;
			(nna)++;
		}

		if ((sl[i].iB4 > -1) && (fabs(sl[i].ab4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB4; Ah.values[nna] = -sl[i].ab4;
			(nna)++;
		}
		if ((sl[i].iE4 > -1) && (fabs(sl[i].ae4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE4; Ah.values[nna] = -sl[i].ae4;
			(nna)++;
		}
		if ((sl[i].iN4 > -1) && (fabs(sl[i].an4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN4; Ah.values[nna] = -sl[i].an4;
			(nna)++;
		}
		if ((sl[i].iS4 > -1) && (fabs(sl[i].as4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS4; Ah.values[nna] = -sl[i].as4;
			(nna)++;
		}
		if ((sl[i].iT4 > -1) && (fabs(sl[i].at4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT4; Ah.values[nna] = -sl[i].at4;
			(nna)++;
		}
		if ((sl[i].iW4 > -1) && (fabs(sl[i].aw4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW4; Ah.values[nna] = -sl[i].aw4;
			(nna)++;
		}

	}

	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) {
			Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iW; Ah.values[nna] = slb[i].aw;
			(nna)++;
		}
		if ((slb[i].iI > -1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iI; Ah.values[nna] = -slb[i].ai;
			(nna)++;
		}
	}


	cusp::array1d<ScalarType, cusp::host_memory> xh(nnu);
	cusp::array1d<ScalarType, cusp::host_memory> yh(nnu);
	//cusp::array1d<float, cusp::device_memory> xh(nnu);
	//cusp::array1d<float, cusp::device_memory> yh(nnu);

	//printf("0i5\n");
	//getchar();

	// yh = Ah x xh;
	for (integer i = 0; i < maxelm + maxbound; i++) {
		// b == dV[i];
		yh[i] = dV[i];
		xh[i] = 0.0;
		//xh[i] = dX0[i];
	}

	//printf("1\n");
	//getchar();
	// convert host->device
	//cusp::array1d<float, cusp::device_memory> xd = xh;
	//cusp::array1d<float, cusp::device_memory> xd(xh);
	//printf("2\n");
	//getchar();
	//cusp::array1d<float, cusp::device_memory> yd = yh;
	//cusp::array1d<float, cusp::device_memory> yd(yh);
	//printf("3\n");
	//getchar();

	//cusp::coo_matrix<int, float, cusp::device_memory> Ad = Ah;
	//cusp::coo_matrix<int, float, cusp::device_memory> Ad(Ah);
	//printf("4\n");
	//getchar();

	// Управление решателем 25 октября 2016.
	doublereal tolerance = 1e-8;
	if (bSIMPLErun_now_for_temperature) {
		// Эти значения невязок для CFD задач были 
		// успешно опробованы на задаче теплового расчёта
		// радиатора водяного охлаждения 3л/мин (совместное решение cfd + temperature 
		// + приближение Обербека-Буссинеска.).
		switch (iVar) {
		case VX: tolerance = 1e-5;  break; //5e-5
		case VY: tolerance = 1e-5;  break; // 5e-5
		case VZ: tolerance = 1e-5;  break; // 5e-5
		case TEMP: tolerance = 1e-7;  break; // 1e-7
		case PAM: tolerance = 1e-6; break; // 1e-6
		}
	}

	integer imaxiter = 6000;
	if (((adiabatic_vs_heat_transfer_coeff > 0) || (breakRUMBAcalc_for_nonlinear_boundary_condition))) {
		// для нелинейных задач. 
		imaxiter = 1;
		tolerance = 0.1;
	}

	if (bonly_solid_calculation) {
		tolerance = 1e-6;
		imaxiter = 6000;
	}

	/*
	if ((bonly_solid_calculation)&&(breakRUMBAcalc_for_nonlinear_boundary_condition)) {
	// 1 V - cycle,
	// относительное значение невязки 0.1.
	imaxiter = 1;
	tolerance = 0.1;
	getchar();
	}
	*/
	/*
	if (bstart7 && (((adiabatic_vs_heat_transfer_coeff > 0) || (breakRUMBAcalc_for_nonlinear_boundary_condition))))
	{
		imaxiter = 1;
		tolerance = 0.1;
		//cusp::relaxation::gauss_seidel(Ah, xh, yh, monitor); // 15_01_2017

		cusp::array1d<ScalarType, cusp::host_memory> rh(nnu);
		cusp::csr_matrix<int, ScalarType, cusp::host_memory> Ah7(Ah);

		// Construct gauss_seidel relaxation class
		cusp::relaxation::gauss_seidel<ScalarType, cusp::host_memory> Mh7(Ah7);
		// Compute initial residual
		cusp::multiply(Ah7, xh, rh);
		cusp::blas::axpy(yh, rh, ScalarType(-1));
		// Construct monitor with stopping criteria of 100 iterations or 1e-4 residual error
		cusp::monitor<ScalarType> monitor(yh, 100, 1e-4, 0, true);
		// Iteratively solve system
		while (!monitor.finished(rh))
		{
			Mh7(Ah7, yh, xh);
			cusp::multiply(Ah7, xh, rh);
			cusp::blas::axpy(yh, rh, ScalarType(-1));
			++monitor;
		}

		//bstart7 = false;

		if (bonly_solid_calculation) {
			// report status
			monitor.print();
		}
	}
	else {
	*/
		// set stopping criteria
		// iteration_limit = 100
		// relative_tolerance = 1e-6
		//  absolute_tolerance = 0
		//  verbose            = true
		//cusp::monitor<float> monitor(yd, 2000, 1e-6, 0, true);

		//cusp::default_monitor<float> monitor(yd, 6000, 1e-12);
		//-->cusp::default_monitor<ScalarType> monitor(yh, 6000, 1e-8);
		//cusp::default_monitor<ScalarType> monitor(yh, imaxiter, tolerance);
		cusp::monitor<ScalarType> monitor(yh, imaxiter, tolerance);
		//cusp::default_monitor<float> monitor(yd, 2000, 1e-6);

		// setup preconditioner
		//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> Md(Ad);
		// setup preconditioner
		cusp::precond::aggregation::smoothed_aggregation<int, ScalarType, cusp::host_memory> Mh(Ah);
		//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> Md(Mh);
		// Диагональный предобуславливатель.
		//cusp::precond::diagonal<float, cusp::device_memory> Md(Ad);
		// AINV (NS Bridson).
		//cusp::precond::scaled_bridson_ainv<float, cusp::device_memory> Md(Ad);

		if (bonly_solid_calculation) {
			Mh.print();
		}

		// solve A * x = y to default tolerance with BiCGStab
		// with preconditioned BiCGStab
		//--->//cusp::krylov::bicgstab(Ad, xd, yd, monitor, Md);
		//cusp::krylov::bicgstab(Ad, xd, yd, monitor);
		cusp::krylov::bicgstab(Ah, xh, yh, monitor, Mh); // 15_01_2017
		//integer imy_restart_gmres = 20; // 20 и 50, 2000, 4000 не работают.
		//cusp::krylov::gmres(Ah, xh, yh, imy_restart_gmres, monitor, Mh); // не работает.

		if (bonly_solid_calculation) {
			// report status
			monitor.print();
		}

	//}

	//cusp::array1d<float, cusp::host_memory> xh_ret = xd;
	//cusp::array1d<float, cusp::host_memory> xh_ret(xd);
	//printf("5\n");
	//getchar();



	// Возвращение результата.
	for (integer i = 0; i < maxelm + maxbound; i++) {
		//dX0[i]=xh_ret[i];
		dX0[i] = xh[i];
	}

} // cusp_solver_amghost

#endif

void cudasafe(int error, char* message, char* file, int line) {
	if (error != cudaSuccess) {
		fprintf(stderr, "CUDA Error: %s : %i. In %s line %d\n", message, error, file, line);
		exit(-1);
	}
}


// Это вызов библиотечного решателя систем линейных алгебраических уравнений.
// Библиотека Cusp версии 0.5.1. Это библиотека с открытым исходным кодом распространяющаяся 
// по open Source лицензии Apache license 2.0. 
// На хосте (центральном процессоре) подключён метод bicgstab и алгебраический многосеточный метод
// на основе сглаженной аггрегации samg.
// Дата подсоединения 12 октября 2016 года.
void cusp_solver_GPU_SAMG(equation3D* &sl, equation3D_bon* &slb,
	integer maxelm, integer maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit,
	doublereal alpharelax, integer iVar)
{

	// maxit,  iVar - не используется.


	if (dX0 == NULL) {
		dX0 = new doublereal[maxelm + maxbound];
		for (integer i = 0; i < maxelm + maxbound; i++) {
			dX0[i] = 0.0;
		}
	}

	// TODO получить val, col_ind, row_ptr
	integer nna = 0; // количество ненулевых элементов в матрице СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	// подсчёт числа ненулевых элементов в матрице.
	nna = 0;
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB2>-1) && (fabs(sl[i].ab2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE2>-1) && (fabs(sl[i].ae2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN2>-1) && (fabs(sl[i].an2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS2>-1) && (fabs(sl[i].as2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT2>-1) && (fabs(sl[i].at2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW2>-1) && (fabs(sl[i].aw2) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB3>-1) && (fabs(sl[i].ab3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE3>-1) && (fabs(sl[i].ae3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN3>-1) && (fabs(sl[i].an3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS3>-1) && (fabs(sl[i].as3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT3>-1) && (fabs(sl[i].at3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW3>-1) && (fabs(sl[i].aw3) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB4>-1) && (fabs(sl[i].ab4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE4>-1) && (fabs(sl[i].ae4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN4>-1) && (fabs(sl[i].an4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS4>-1) && (fabs(sl[i].as4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT4>-1) && (fabs(sl[i].at4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW4>-1) && (fabs(sl[i].aw4) > nonzeroEPS)) (nna)++;

	}
	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	}

	integer nnu = 0; // число неизвестных.
	nnu = maxelm + maxbound;

	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
	//typedef float ScalarType; // 11.04.2019
	typedef cusp::device_memory MemorySpace;


	// allocate storage for (nnu,nnu) matrix with nna nonzeros
	cusp::coo_matrix<int, ScalarType, cusp::host_memory> Ah(nnu, nnu, nna);
	//cusp::coo_matrix<int, float, cusp::device_memory> Ah(nnu, nnu, nna);

	//printf("0\n");
	//getchar();

	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iP; Ah.values[nna] = sl[i].ap / alpharelax;
			(nna)++;
		}
		if ((sl[i].iB > -1) && (fabs(sl[i].ab) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB; Ah.values[nna] = -sl[i].ab;
			(nna)++;
		}
		if ((sl[i].iE > -1) && (fabs(sl[i].ae) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE; Ah.values[nna] = -sl[i].ae;
			(nna)++;
		}
		if ((sl[i].iN > -1) && (fabs(sl[i].an) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN; Ah.values[nna] = -sl[i].an;
			(nna)++;
		}
		if ((sl[i].iS > -1) && (fabs(sl[i].as) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS; Ah.values[nna] = -sl[i].as;
			(nna)++;
		}
		if ((sl[i].iT > -1) && (fabs(sl[i].at) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT; Ah.values[nna] = -sl[i].at;
			(nna)++;
		}
		if ((sl[i].iW > -1) && (fabs(sl[i].aw) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW; Ah.values[nna] = -sl[i].aw;
			(nna)++;
		}

		if ((sl[i].iB2 > -1) && (fabs(sl[i].ab2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB2; Ah.values[nna] = -sl[i].ab2;
			(nna)++;
		}
		if ((sl[i].iE2 > -1) && (fabs(sl[i].ae2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE2; Ah.values[nna] = -sl[i].ae2;
			(nna)++;
		}
		if ((sl[i].iN2 > -1) && (fabs(sl[i].an2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN2; Ah.values[nna] = -sl[i].an2;
			(nna)++;
		}
		if ((sl[i].iS2 > -1) && (fabs(sl[i].as2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS2; Ah.values[nna] = -sl[i].as2;
			(nna)++;
		}
		if ((sl[i].iT2 > -1) && (fabs(sl[i].at2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT2; Ah.values[nna] = -sl[i].at2;
			(nna)++;
		}
		if ((sl[i].iW2 > -1) && (fabs(sl[i].aw2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW2; Ah.values[nna] = -sl[i].aw2;
			(nna)++;
		}

		if ((sl[i].iB3 > -1) && (fabs(sl[i].ab3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB3; Ah.values[nna] = -sl[i].ab3;
			(nna)++;
		}
		if ((sl[i].iE3 > -1) && (fabs(sl[i].ae3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE3; Ah.values[nna] = -sl[i].ae3;
			(nna)++;
		}
		if ((sl[i].iN3 > -1) && (fabs(sl[i].an3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN3; Ah.values[nna] = -sl[i].an3;
			(nna)++;
		}
		if ((sl[i].iS3 > -1) && (fabs(sl[i].as3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS3; Ah.values[nna] = -sl[i].as3;
			(nna)++;
		}
		if ((sl[i].iT3 > -1) && (fabs(sl[i].at3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT3; Ah.values[nna] = -sl[i].at3;
			(nna)++;
		}
		if ((sl[i].iW3 > -1) && (fabs(sl[i].aw3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW3; Ah.values[nna] = -sl[i].aw3;
			(nna)++;
		}

		if ((sl[i].iB4 > -1) && (fabs(sl[i].ab4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB4; Ah.values[nna] = -sl[i].ab4;
			(nna)++;
		}
		if ((sl[i].iE4 > -1) && (fabs(sl[i].ae4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE4; Ah.values[nna] = -sl[i].ae4;
			(nna)++;
		}
		if ((sl[i].iN4 > -1) && (fabs(sl[i].an4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN4; Ah.values[nna] = -sl[i].an4;
			(nna)++;
		}
		if ((sl[i].iS4 > -1) && (fabs(sl[i].as4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS4; Ah.values[nna] = -sl[i].as4;
			(nna)++;
		}
		if ((sl[i].iT4 > -1) && (fabs(sl[i].at4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT4; Ah.values[nna] = -sl[i].at4;
			(nna)++;
		}
		if ((sl[i].iW4 > -1) && (fabs(sl[i].aw4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW4; Ah.values[nna] = -sl[i].aw4;
			(nna)++;
		}

	}

	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) {
			Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iW; Ah.values[nna] = slb[i].aw;
			(nna)++;
		}
		if ((slb[i].iI > -1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iI; Ah.values[nna] = -slb[i].ai;
			(nna)++;
		}
	}


	cusp::array1d<ScalarType, cusp::host_memory> xh(nnu);
	cusp::array1d<ScalarType, cusp::host_memory> yh(nnu);
	//cusp::array1d<float, cusp::device_memory> xh(nnu);
	//cusp::array1d<float, cusp::device_memory> yh(nnu);

	//printf("0i5\n");
	//getchar();

	// yh = Ah x xh;
	for (integer i = 0; i < maxelm + maxbound; i++) {
		// b == dV[i];
		yh[i] = dV[i];
		xh[i] = 0.0;
		//xh[i] = dX0[i];
	}

	//printf("1\n");
	//getchar();
	// convert host->device
	//cusp::array1d<float, cusp::device_memory> xd = xh;
	//cusp::array1d<float, cusp::device_memory> xd(xh);
	//printf("2\n");
	//getchar();
	//cusp::array1d<float, cusp::device_memory> yd = yh;
	//cusp::array1d<float, cusp::device_memory> yd(yh);
	//printf("3\n");
	//getchar();

	//cusp::coo_matrix<int, float, cusp::device_memory> Ad = Ah;
	//cusp::coo_matrix<int, float, cusp::device_memory> Ad(Ah);
	//printf("4\n");
	//getchar();

	// Управление решателем 25 октября 2016.
	doublereal tolerance = 1e-8;
	if (bSIMPLErun_now_for_temperature) {
		// Эти значения невязок для CFD задач были 
		// успешно опробованы на задаче теплового расчёта
		// радиатора водяного охлажения (совместное решение cfd + temperature 
		// + приближение Обербека-Буссинеска.).
		switch (iVar) {
		case VX: tolerance = 1e-5;  break; //5e-5
		case VY: tolerance = 1e-5;  break; // 5e-5
		case VZ: tolerance = 1e-5;  break; // 5e-5
		case TEMP: tolerance = 1e-7;  break; // 1e-7
		case PAM: tolerance = 1e-6; break; // 1e-6
		}
	}

	integer imaxiter = 6000;
	if (((adiabatic_vs_heat_transfer_coeff > 0) || (breakRUMBAcalc_for_nonlinear_boundary_condition))) {
		// для нелинейных задач. Граничное условие Стефана-Больцмана.
		// 348878024515312.94
		// 1947852997768.51
		imaxiter = 1;
		tolerance = 0.1;
	}

	if (bonly_solid_calculation) {
		tolerance = 1e-6;
		imaxiter = 6000;
	}

	/*
	if ((bonly_solid_calculation)&&(breakRUMBAcalc_for_nonlinear_boundary_condition)) {
	// 1 V - cycle,
	// относительное значение невязки 0.1.
	imaxiter = 1;
	tolerance = 0.1;
	getchar();
	}
	*/
	/*
	if (bstart7 && (((adiabatic_vs_heat_transfer_coeff > 0) || (breakRUMBAcalc_for_nonlinear_boundary_condition))))
	{
		imaxiter = 1;
		tolerance = 0.1;
		//cusp::relaxation::gauss_seidel(Ah, xh, yh, monitor); // 15_01_2017

		cusp::array1d<ScalarType, cusp::host_memory> rh(nnu);
		cusp::csr_matrix<int, ScalarType, cusp::host_memory> Ah7(Ah);

		// Construct gauss_seidel relaxation class
		cusp::relaxation::gauss_seidel<ScalarType, cusp::host_memory> Mh7(Ah7);
		// Compute initial residual
		cusp::multiply(Ah7, xh, rh);
		cusp::blas::axpy(yh, rh, ScalarType(-1));
		// Construct monitor with stopping criteria of 100 iterations or 1e-4 residual error
		cusp::monitor<ScalarType> monitor(yh, 100, 1e-4, 0, true);
		// Iteratively solve system
		while (!monitor.finished(rh))
		{
			Mh7(Ah7, yh, xh);
			cusp::multiply(Ah7, xh, rh);
			cusp::blas::axpy(yh, rh, ScalarType(-1));
			++monitor;
		}

		//bstart7 = false;

		if (bonly_solid_calculation) {
			// report status
			monitor.print();
		}

		// Возвращение результата.
		for (integer i = 0; i < maxelm + maxbound; i++) {
			//dX0[i] = xh_ret[i];
			dX0[i] = xh[i];
		}

	}
	else {
	*/

	
		int deviceCount;

		cudasafe(cudaGetDeviceCount(&deviceCount), "GetDeviceCount", __FILE__, __LINE__);

		printf("Number of CUDA devices %d.\n", deviceCount);
		cudaDeviceProp *deviceProp = NULL;
		if (deviceCount > 0) {
			deviceProp = new cudaDeviceProp[deviceCount];
		}

		for (int dev = 0; dev < deviceCount; dev++) {
			
			cudaDeviceProp deviceProp1;

			cudasafe(cudaGetDeviceProperties(&deviceProp1, dev), "Get Device Properties", __FILE__, __LINE__);
			deviceProp[dev] = deviceProp1;

			if (dev == 0) {
				if (deviceProp1.major == 9999 && deviceProp1.minor == 9999) {
					printf("No CUDA GPU has been detected\n");
					system("PAUSE");
					exit(1);
				}
				else if (deviceCount == 1) {
					printf("There is 1 device supporting CUDA\n");
				}
				else {
					printf("There are %d devices supporting CUDA\n", deviceCount);
				}
			}

			printf("For device #%d\n", dev);
			printf("Device name:                %s\n", deviceProp1.name);
			printf("Major revision number:      %d\n", deviceProp1.major);
			printf("Minor revision Number:      %d\n", deviceProp1.minor);
			printf("Total Global Memory:        %lld\n", deviceProp1.totalGlobalMem);
			printf("Total shared mem per block: %lld\n", deviceProp1.sharedMemPerBlock);
			printf("Total const mem size:       %lld\n", deviceProp1.totalConstMem);
			printf("Warp size:                  %lld\n", deviceProp1.warpSize);
			printf("Maximum block dimensions:   %lld x %lld x %lld\n", deviceProp1.maxThreadsDim[0], \
				deviceProp1.maxThreadsDim[1], \
				deviceProp1.maxThreadsDim[2]);

			printf("Maximum grid dimensions:    %lld x %lld x %lld\n", deviceProp1.maxGridSize[0], \
				deviceProp1.maxGridSize[1], \
				deviceProp1.maxGridSize[2]);
			printf("Clock Rate:                 %d\n", deviceProp1.clockRate);
			printf("Number of muliprocessors:   %d\n", deviceProp1.multiProcessorCount);

		}

		//Select compute - device which best matches criteria.
		//__host__ ​cudaError_t cudaChooseDevice ( int* device, const cudaDeviceProp* prop )
		cudaChooseDevice(&deviceCount, deviceProp);
		delete[] deviceProp;
		deviceProp = NULL;

		int dev;
		cudaGetDevice(&dev);
		printf("ID of current CUDA device: %d\n", dev);
	

		// set stopping criteria
		// iteration_limit = 100
		// relative_tolerance = 1e-6
		//  absolute_tolerance = 0
		//  verbose            = true
		//cusp::monitor<float> monitor(yd, 2000, 1e-6, 0, true);

		//cusp::default_monitor<float> monitor(yd, 6000, 1e-12);
		//-->cusp::default_monitor<ScalarType> monitor(yh, 6000, 1e-8);
		//cusp::default_monitor<ScalarType> monitor(yh, imaxiter, tolerance);
		cusp::monitor<ScalarType> monitor(yh, imaxiter, tolerance);
		//cusp::default_monitor<float> monitor(yd, 2000, 1e-6);


		cusp::array1d<ScalarType, cusp::device_memory> xd(xh);
		cusp::array1d<ScalarType, cusp::device_memory> yd(yh);
		cusp::coo_matrix<int, ScalarType, cusp::device_memory> Ad(Ah);

		// setup preconditioner
		cusp::precond::aggregation::smoothed_aggregation<int, ScalarType, cusp::device_memory> Md(Ad);
		// setup preconditioner
		//host
		//cusp::precond::aggregation::smoothed_aggregation<int, ScalarType, cusp::host_memory> Mh(Ah);
		//end host
		//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> Md(Ah);
		// Диагональный предобуславливатель.
		//cusp::precond::diagonal<float, cusp::device_memory> Md(Ad);
		// AINV (NS Bridson).
		//cusp::precond::scaled_bridson_ainv<float, cusp::device_memory> Md(Ad);

		if (bonly_solid_calculation) {
			Md.print();
		}

		// solve A * x = y to default tolerance with BiCGStab
		// with preconditioned BiCGStab
		cusp::krylov::bicgstab(Ad, xd, yd, monitor, Md);
		//cusp::krylov::bicgstab(Ad, xd, yd, monitor);
		//host
		//cusp::krylov::bicgstab(Ah, xh, yh, monitor, Mh); // 15_01_2017
		//end host

		//cusp::krylov::bicgstab(Ad, xd, yd, monitor, Md); // 21_10_2017

		cusp::array1d<ScalarType, cusp::host_memory> xh_ret(xd);

		//integer imy_restart_gmres = 20; // 20 и 50, 2000, 4000 не работают.
		//cusp::krylov::gmres(Ah, xh, yh, imy_restart_gmres, monitor, Mh); // не работает.

		if (bonly_solid_calculation) {
			// report status
			monitor.print();
		}

		// Возвращение результата.
		for (integer i = 0; i < maxelm + maxbound; i++) {
			dX0[i] = xh_ret[i];
			//dX0[i] = xh[i];
		}

	//}

	//cusp::array1d<float, cusp::host_memory> xh_ret = xd;
	//cusp::array1d<float, cusp::host_memory> xh_ret(xd);
	//printf("5\n");
	//getchar();

} // cusp_solver_amgdevice

// Это вызов библиотечного решателя систем линейных алгебраических уравнений.
// Библиотека Cusp версии 0.5.1. Это библиотека с открытым исходным кодом распространяющаяся 
// по open Source лицензии Apache license 2.0. 
// На видеокарте (графическом процессоре) подключён метод bicgstab и AINV (NS Brigson) метод
// в качестве предобуславливателя.
// Дата подсоединения 12 октября 2016 года.
void cusp_solver_GPU_AINV_Bridson(equation3D* &sl, equation3D_bon* &slb,
	integer maxelm, integer maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit,
	doublereal alpharelax, integer iVar)
{

	// maxit,  iVar - не используется.


	if (dX0 == NULL) {
		dX0 = new doublereal[maxelm + maxbound];
		for (integer i = 0; i < maxelm + maxbound; i++) {
			dX0[i] = 0.0;
		}
	}

	// TODO получить val, col_ind, row_ptr
	integer nna = 0; // количество ненулевых элементов в матрице СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	// подсчёт числа ненулевых элементов в матрице.
	nna = 0;
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB2>-1) && (fabs(sl[i].ab2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE2>-1) && (fabs(sl[i].ae2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN2>-1) && (fabs(sl[i].an2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS2>-1) && (fabs(sl[i].as2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT2>-1) && (fabs(sl[i].at2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW2>-1) && (fabs(sl[i].aw2) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB3>-1) && (fabs(sl[i].ab3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE3>-1) && (fabs(sl[i].ae3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN3>-1) && (fabs(sl[i].an3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS3>-1) && (fabs(sl[i].as3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT3>-1) && (fabs(sl[i].at3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW3>-1) && (fabs(sl[i].aw3) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB4>-1) && (fabs(sl[i].ab4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE4>-1) && (fabs(sl[i].ae4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN4>-1) && (fabs(sl[i].an4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS4>-1) && (fabs(sl[i].as4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT4>-1) && (fabs(sl[i].at4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW4>-1) && (fabs(sl[i].aw4) > nonzeroEPS)) (nna)++;

	}
	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	}

	integer nnu = 0; // число неизвестных.
	nnu = maxelm + maxbound;

	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
	//typedef float ScalarType;

	// allocate storage for (nnu,nnu) matrix with nna nonzeros
	cusp::coo_matrix<int, ScalarType, cusp::host_memory> Ah(nnu, nnu, nna);
	//cusp::coo_matrix<int, float, cusp::device_memory> Ah(nnu, nnu, nna);

	//printf("0\n");
	//getchar();

	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iP; Ah.values[nna] = sl[i].ap / alpharelax;
			(nna)++;
		}
		if ((sl[i].iB > -1) && (fabs(sl[i].ab) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB; Ah.values[nna] = -sl[i].ab;
			(nna)++;
		}
		if ((sl[i].iE > -1) && (fabs(sl[i].ae) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE; Ah.values[nna] = -sl[i].ae;
			(nna)++;
		}
		if ((sl[i].iN > -1) && (fabs(sl[i].an) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN; Ah.values[nna] = -sl[i].an;
			(nna)++;
		}
		if ((sl[i].iS > -1) && (fabs(sl[i].as) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS; Ah.values[nna] = -sl[i].as;
			(nna)++;
		}
		if ((sl[i].iT > -1) && (fabs(sl[i].at) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT; Ah.values[nna] = -sl[i].at;
			(nna)++;
		}
		if ((sl[i].iW > -1) && (fabs(sl[i].aw) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW; Ah.values[nna] = -sl[i].aw;
			(nna)++;
		}

		if ((sl[i].iB2 > -1) && (fabs(sl[i].ab2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB2; Ah.values[nna] = -sl[i].ab2;
			(nna)++;
		}
		if ((sl[i].iE2 > -1) && (fabs(sl[i].ae2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE2; Ah.values[nna] = -sl[i].ae2;
			(nna)++;
		}
		if ((sl[i].iN2 > -1) && (fabs(sl[i].an2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN2; Ah.values[nna] = -sl[i].an2;
			(nna)++;
		}
		if ((sl[i].iS2 > -1) && (fabs(sl[i].as2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS2; Ah.values[nna] = -sl[i].as2;
			(nna)++;
		}
		if ((sl[i].iT2 > -1) && (fabs(sl[i].at2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT2; Ah.values[nna] = -sl[i].at2;
			(nna)++;
		}
		if ((sl[i].iW2 > -1) && (fabs(sl[i].aw2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW2; Ah.values[nna] = -sl[i].aw2;
			(nna)++;
		}

		if ((sl[i].iB3 > -1) && (fabs(sl[i].ab3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB3; Ah.values[nna] = -sl[i].ab3;
			(nna)++;
		}
		if ((sl[i].iE3 > -1) && (fabs(sl[i].ae3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE3; Ah.values[nna] = -sl[i].ae3;
			(nna)++;
		}
		if ((sl[i].iN3 > -1) && (fabs(sl[i].an3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN3; Ah.values[nna] = -sl[i].an3;
			(nna)++;
		}
		if ((sl[i].iS3 > -1) && (fabs(sl[i].as3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS3; Ah.values[nna] = -sl[i].as3;
			(nna)++;
		}
		if ((sl[i].iT3 > -1) && (fabs(sl[i].at3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT3; Ah.values[nna] = -sl[i].at3;
			(nna)++;
		}
		if ((sl[i].iW3 > -1) && (fabs(sl[i].aw3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW3; Ah.values[nna] = -sl[i].aw3;
			(nna)++;
		}

		if ((sl[i].iB4 > -1) && (fabs(sl[i].ab4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB4; Ah.values[nna] = -sl[i].ab4;
			(nna)++;
		}
		if ((sl[i].iE4 > -1) && (fabs(sl[i].ae4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE4; Ah.values[nna] = -sl[i].ae4;
			(nna)++;
		}
		if ((sl[i].iN4 > -1) && (fabs(sl[i].an4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN4; Ah.values[nna] = -sl[i].an4;
			(nna)++;
		}
		if ((sl[i].iS4 > -1) && (fabs(sl[i].as4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS4; Ah.values[nna] = -sl[i].as4;
			(nna)++;
		}
		if ((sl[i].iT4 > -1) && (fabs(sl[i].at4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT4; Ah.values[nna] = -sl[i].at4;
			(nna)++;
		}
		if ((sl[i].iW4 > -1) && (fabs(sl[i].aw4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW4; Ah.values[nna] = -sl[i].aw4;
			(nna)++;
		}

	}

	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) {
			Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iW; Ah.values[nna] = slb[i].aw;
			(nna)++;
		}
		if ((slb[i].iI > -1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iI; Ah.values[nna] = -slb[i].ai;
			(nna)++;
		}
	}

	int deviceCount;

	cudasafe(cudaGetDeviceCount(&deviceCount), "GetDeviceCount", __FILE__, __LINE__);

	printf("Number of CUDA devices %d.\n", deviceCount);
	cudaDeviceProp *deviceProp = NULL;
	if (deviceCount > 0) {
		deviceProp = new cudaDeviceProp[deviceCount];
	}

	for (int dev = 0; dev < deviceCount; dev++) {

		cudaDeviceProp deviceProp1;

		cudasafe(cudaGetDeviceProperties(&deviceProp1, dev), "Get Device Properties", __FILE__, __LINE__);
		deviceProp[dev] = deviceProp1;

		if (dev == 0) {
			if (deviceProp1.major == 9999 && deviceProp1.minor == 9999) {
				printf("No CUDA GPU has been detected\n");
				system("PAUSE");
				exit(1);
			}
			else if (deviceCount == 1) {
				printf("There is 1 device supporting CUDA\n");
			}
			else {
				printf("There are %d devices supporting CUDA\n", deviceCount);
			}
		}

		printf("For device #%d\n", dev);
		printf("Device name:                %s\n", deviceProp1.name);
		printf("Major revision number:      %d\n", deviceProp1.major);
		printf("Minor revision Number:      %d\n", deviceProp1.minor);
		printf("Total Global Memory:        %lld\n", deviceProp1.totalGlobalMem);
		printf("Total shared mem per block: %lld\n", deviceProp1.sharedMemPerBlock);
		printf("Total const mem size:       %lld\n", deviceProp1.totalConstMem);
		printf("Warp size:                  %lld\n", deviceProp1.warpSize);
		printf("Maximum block dimensions:   %lld x %lld x %lld\n", deviceProp1.maxThreadsDim[0], \
			deviceProp1.maxThreadsDim[1], \
			deviceProp1.maxThreadsDim[2]);

		printf("Maximum grid dimensions:    %lld x %lld x %lld\n", deviceProp1.maxGridSize[0], \
			deviceProp1.maxGridSize[1], \
			deviceProp1.maxGridSize[2]);
		printf("Clock Rate:                 %d\n", deviceProp1.clockRate);
		printf("Number of muliprocessors:   %d\n", deviceProp1.multiProcessorCount);

	}

	//Select compute - device which best matches criteria.
	//__host__ ​cudaError_t cudaChooseDevice ( int* device, const cudaDeviceProp* prop )
	cudaChooseDevice(&deviceCount, deviceProp);
	delete[] deviceProp;
	deviceProp = NULL;

	int dev;
	cudaGetDevice(&dev);
	printf("ID of current CUDA device: %d\n", dev);


	cusp::array1d<ScalarType, cusp::host_memory> xh(nnu);
	cusp::array1d<ScalarType, cusp::host_memory> yh(nnu);
	//cusp::array1d<float, cusp::device_memory> xh(nnu);
	//cusp::array1d<float, cusp::device_memory> yh(nnu);

	//printf("0i5\n");
	//getchar();

	// yh = Ah x xh;
	for (integer i = 0; i < maxelm + maxbound; i++) {
		// b == dV[i];
		yh[i] = dV[i];
		xh[i] = 0.0;
		//xh[i] = dX0[i];
	}


	//printf("1\n");
	//getchar();
	// convert host->device
	//cusp::array1d<float, cusp::device_memory> xd = xh;
	cusp::array1d<ScalarType, cusp::device_memory> xd(xh);
	//printf("2\n");
	//getchar();
	//cusp::array1d<float, cusp::device_memory> yd = yh;
	cusp::array1d<ScalarType, cusp::device_memory> yd(yh);
	//printf("3\n");
	//getchar();

	//cusp::coo_matrix<int, float, cusp::device_memory> Ad = Ah;
	cusp::coo_matrix<int, ScalarType, cusp::device_memory> Ad(Ah);
	//printf("4\n");
	//getchar();


	// set stopping criteria
	// iteration_limit = 100
	// relative_tolerance = 1e-6
	//  absolute_tolerance = 0
	//  verbose            = true
	//cusp::monitor<float> monitor(yd, 2000, 1e-6, 0, true);

	//cusp::default_monitor<float> monitor(yd, 6000, 1e-12);
	// 1e-6 достаточно для стационарных задач.
	//cusp::default_monitor<ScalarType> monitor(yh, 6000, 1e-8);
	//cusp::monitor<ScalarType> monitor(yh, 6000, 1e-8);
	//cusp::default_monitor<float> monitor(yd, 2000, 1e-6); 
	cusp::monitor<ScalarType> monitor(yh, 1000, 1e-6);

	// setup preconditioner
	//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> Md(Ad);
	// setup preconditioner
	//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::host_memory> Mh(Ah);
	//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> Md(Mh);
	// Диагональный предобуславливатель.
	//cusp::precond::diagonal<float, cusp::device_memory> Md(Ad);
	// AINV (NS Bridson).
	cusp::precond::scaled_bridson_ainv<ScalarType, cusp::device_memory> Md(Ad);


	// solve A * x = y to default tolerance with BiCGStab
	// with preconditioned BiCGStab
	cusp::krylov::bicgstab(Ad, xd, yd, monitor, Md);
	//cusp::krylov::bicgstab(Ad, xd, yd, monitor);
	//cusp::krylov::bicgstab(Ah, xh, yh, monitor, Mh);



	//cusp::array1d<float, cusp::host_memory> xh_ret = xd;
	cusp::array1d<ScalarType, cusp::host_memory> xh_ret(xd);
	//printf("5\n");
	//getchar();

	// Возвращение результата.
	for (integer i = 0; i < maxelm + maxbound; i++) {
		dX0[i] = xh_ret[i];
		//dX0[i] = xh[i];
	}

} // cusp_solver

  // Это вызов библиотечного решателя систем линейных алгебраических уравнений.
  // Библиотека Cusp версии 0.5.1. Это библиотека с открытым исходным кодом распространяющаяся 
  // по open Source лицензии Apache license 2.0. 
  // На видеокарте (графическом процессоре) подключён метод bicgstab и AINV (NS Brigson) метод
  // в качестве предобуславливателя.
  // Дата подсоединения 12 октября 2016 года.
void cusp_solver_host(equation3D* &sl, equation3D_bon* &slb,
	integer maxelm, integer maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit,
	doublereal alpharelax, integer iVar)
{

	// maxit,  iVar - не используется.


	if (dX0 == NULL) {
		dX0 = new doublereal[maxelm + maxbound];
		for (integer i = 0; i < maxelm + maxbound; i++) {
			dX0[i] = 0.0;
		}
	}

	// TODO получить val, col_ind, row_ptr
	integer nna = 0; // количество ненулевых элементов в матрице СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

										 // подсчёт числа ненулевых элементов в матрице.
	nna = 0;
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB2>-1) && (fabs(sl[i].ab2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE2>-1) && (fabs(sl[i].ae2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN2>-1) && (fabs(sl[i].an2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS2>-1) && (fabs(sl[i].as2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT2>-1) && (fabs(sl[i].at2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW2>-1) && (fabs(sl[i].aw2) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB3>-1) && (fabs(sl[i].ab3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE3>-1) && (fabs(sl[i].ae3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN3>-1) && (fabs(sl[i].an3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS3>-1) && (fabs(sl[i].as3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT3>-1) && (fabs(sl[i].at3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW3>-1) && (fabs(sl[i].aw3) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB4>-1) && (fabs(sl[i].ab4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE4>-1) && (fabs(sl[i].ae4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN4>-1) && (fabs(sl[i].an4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS4>-1) && (fabs(sl[i].as4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT4>-1) && (fabs(sl[i].at4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW4>-1) && (fabs(sl[i].aw4) > nonzeroEPS)) (nna)++;

	}
	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	}

	integer nnu = 0; // число неизвестных.
	nnu = maxelm + maxbound;

	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
									   //typedef float ScalarType;

									   // allocate storage for (nnu,nnu) matrix with nna nonzeros
	cusp::coo_matrix<int, ScalarType, cusp::host_memory> Ah(nnu, nnu, nna);
	//cusp::coo_matrix<int, float, cusp::device_memory> Ah(nnu, nnu, nna);

	//printf("0\n");
	//getchar();

	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iP; Ah.values[nna] = sl[i].ap / alpharelax;
			(nna)++;
		}
		if ((sl[i].iB > -1) && (fabs(sl[i].ab) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB; Ah.values[nna] = -sl[i].ab;
			(nna)++;
		}
		if ((sl[i].iE > -1) && (fabs(sl[i].ae) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE; Ah.values[nna] = -sl[i].ae;
			(nna)++;
		}
		if ((sl[i].iN > -1) && (fabs(sl[i].an) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN; Ah.values[nna] = -sl[i].an;
			(nna)++;
		}
		if ((sl[i].iS > -1) && (fabs(sl[i].as) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS; Ah.values[nna] = -sl[i].as;
			(nna)++;
		}
		if ((sl[i].iT > -1) && (fabs(sl[i].at) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT; Ah.values[nna] = -sl[i].at;
			(nna)++;
		}
		if ((sl[i].iW > -1) && (fabs(sl[i].aw) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW; Ah.values[nna] = -sl[i].aw;
			(nna)++;
		}

		if ((sl[i].iB2 > -1) && (fabs(sl[i].ab2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB2; Ah.values[nna] = -sl[i].ab2;
			(nna)++;
		}
		if ((sl[i].iE2 > -1) && (fabs(sl[i].ae2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE2; Ah.values[nna] = -sl[i].ae2;
			(nna)++;
		}
		if ((sl[i].iN2 > -1) && (fabs(sl[i].an2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN2; Ah.values[nna] = -sl[i].an2;
			(nna)++;
		}
		if ((sl[i].iS2 > -1) && (fabs(sl[i].as2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS2; Ah.values[nna] = -sl[i].as2;
			(nna)++;
		}
		if ((sl[i].iT2 > -1) && (fabs(sl[i].at2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT2; Ah.values[nna] = -sl[i].at2;
			(nna)++;
		}
		if ((sl[i].iW2 > -1) && (fabs(sl[i].aw2) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW2; Ah.values[nna] = -sl[i].aw2;
			(nna)++;
		}

		if ((sl[i].iB3 > -1) && (fabs(sl[i].ab3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB3; Ah.values[nna] = -sl[i].ab3;
			(nna)++;
		}
		if ((sl[i].iE3 > -1) && (fabs(sl[i].ae3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE3; Ah.values[nna] = -sl[i].ae3;
			(nna)++;
		}
		if ((sl[i].iN3 > -1) && (fabs(sl[i].an3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN3; Ah.values[nna] = -sl[i].an3;
			(nna)++;
		}
		if ((sl[i].iS3 > -1) && (fabs(sl[i].as3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS3; Ah.values[nna] = -sl[i].as3;
			(nna)++;
		}
		if ((sl[i].iT3 > -1) && (fabs(sl[i].at3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT3; Ah.values[nna] = -sl[i].at3;
			(nna)++;
		}
		if ((sl[i].iW3 > -1) && (fabs(sl[i].aw3) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW3; Ah.values[nna] = -sl[i].aw3;
			(nna)++;
		}

		if ((sl[i].iB4 > -1) && (fabs(sl[i].ab4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB4; Ah.values[nna] = -sl[i].ab4;
			(nna)++;
		}
		if ((sl[i].iE4 > -1) && (fabs(sl[i].ae4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE4; Ah.values[nna] = -sl[i].ae4;
			(nna)++;
		}
		if ((sl[i].iN4 > -1) && (fabs(sl[i].an4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN4; Ah.values[nna] = -sl[i].an4;
			(nna)++;
		}
		if ((sl[i].iS4 > -1) && (fabs(sl[i].as4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS4; Ah.values[nna] = -sl[i].as4;
			(nna)++;
		}
		if ((sl[i].iT4 > -1) && (fabs(sl[i].at4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT4; Ah.values[nna] = -sl[i].at4;
			(nna)++;
		}
		if ((sl[i].iW4 > -1) && (fabs(sl[i].aw4) > nonzeroEPS)) {
			Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW4; Ah.values[nna] = -sl[i].aw4;
			(nna)++;
		}

	}

	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) {
			Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iW; Ah.values[nna] = slb[i].aw;
			(nna)++;
		}
		if ((slb[i].iI > -1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iI; Ah.values[nna] = -slb[i].ai;
			(nna)++;
		}
	}


	cusp::array1d<ScalarType, cusp::host_memory> xh(nnu);
	cusp::array1d<ScalarType, cusp::host_memory> yh(nnu);
	//cusp::array1d<float, cusp::device_memory> xh(nnu);
	//cusp::array1d<float, cusp::device_memory> yh(nnu);

	//printf("0i5\n");
	//getchar();

	// yh = Ah x xh;
	for (integer i = 0; i < maxelm + maxbound; i++) {
		// b == dV[i];
		yh[i] = dV[i];
		xh[i] = 0.0;
		//xh[i] = dX0[i];
	}

	/*
	//printf("1\n");
	//getchar();
	// convert host->device
	//cusp::array1d<float, cusp::device_memory> xd = xh;
	cusp::array1d<ScalarType, cusp::device_memory> xd(xh);
	//printf("2\n");
	//getchar();
	//cusp::array1d<float, cusp::device_memory> yd = yh;
	cusp::array1d<ScalarType, cusp::device_memory> yd(yh);
	//printf("3\n");
	//getchar();

	//cusp::coo_matrix<int, float, cusp::device_memory> Ad = Ah;
	cusp::coo_matrix<int, ScalarType, cusp::device_memory> Ad(Ah);
	//printf("4\n");
	//getchar();
	*/

	// set stopping criteria
	// iteration_limit = 100
	// relative_tolerance = 1e-6
	//  absolute_tolerance = 0
	//  verbose            = true
	//cusp::monitor<float> monitor(yd, 2000, 1e-6, 0, true);

	//cusp::default_monitor<float> monitor(yd, 6000, 1e-12);
	// 1e-6 достаточно для стационарных задач.
	//cusp::default_monitor<ScalarType> monitor(yh, 6000, 1e-8);
	//cusp::monitor<ScalarType> monitor(yh, 6000, 1e-8);
	//cusp::default_monitor<float> monitor(yd, 2000, 1e-6); 
	cusp::monitor<ScalarType> monitor(yh, 1000, 1e-6);

	// setup preconditioner
	//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> Md(Ad);
	// setup preconditioner
	//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::host_memory> Mh(Ah);
	//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> Md(Mh);
	// Диагональный предобуславливатель.
	//cusp::precond::diagonal<float, cusp::device_memory> Md(Ad);
	// AINV (NS Bridson).
	cusp::precond::scaled_bridson_ainv<ScalarType, cusp::host_memory> Mh(Ah);


	// solve A * x = y to default tolerance with BiCGStab
	// with preconditioned BiCGStab
	cusp::krylov::bicgstab(Ah, xh, yh, monitor, Mh);
	//cusp::krylov::bicgstab(Ad, xd, yd, monitor);
	//cusp::krylov::bicgstab(Ah, xh, yh, monitor, Mh);



	//cusp::array1d<float, cusp::host_memory> xh_ret = xd;
	//--->cusp::array1d<ScalarType, cusp::host_memory> xh_ret(xd);
	//printf("5\n");
	//getchar();

	// Возвращение результата.
	for (integer i = 0; i < maxelm + maxbound; i++) {
		dX0[i] = xh[i];
		//dX0[i] = xh[i];
	}

} // cusp_solver_host


bool bcusp_gl_first_start = true;
cusp::coo_matrix<int, doublereal, cusp::host_memory> Ah_gl;
cusp::array1d<doublereal, cusp::host_memory> xh_gl;
cusp::array1d<doublereal, cusp::host_memory> yh_gl;
cusp::array1d<doublereal, cusp::device_memory> xd_gl;
cusp::array1d<doublereal, cusp::device_memory> yd_gl;
cusp::coo_matrix<int, doublereal, cusp::device_memory> Ad_gl;

// Это вызов библиотечного решателя систем линейных алгебраических уравнений.
// Библиотека Cusp версии 0.5.1. Это библиотека с открытым исходным кодом распространяющаяся 
// по open Source лицензии Apache license 2.0. 
// На видеокарте (графическом процессоре) подключён метод bicgstab и AINV (NS Brigson) метод
// в качестве предобуславливателя.
// Дата подсоединения 12 октября 2016 года.
void cusp_solver_global_allocate(equation3D* &sl, equation3D_bon* &slb,
	integer maxelm, integer maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit,
	doublereal alpharelax, integer iVar)
{

	// maxit,  iVar - не используется.


	if (dX0 == NULL) {
		dX0 = new doublereal[maxelm + maxbound];
		for (integer i = 0; i < maxelm + maxbound; i++) {
			dX0[i] = 0.0;
		}
	}

	// TODO получить val, col_ind, row_ptr
	integer nna = 0; // количество ненулевых элементов в матрице СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	// подсчёт числа ненулевых элементов в матрице.
	nna = 0;
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB2>-1) && (fabs(sl[i].ab2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE2>-1) && (fabs(sl[i].ae2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN2>-1) && (fabs(sl[i].an2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS2>-1) && (fabs(sl[i].as2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT2>-1) && (fabs(sl[i].at2) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW2>-1) && (fabs(sl[i].aw2) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB3>-1) && (fabs(sl[i].ab3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE3>-1) && (fabs(sl[i].ae3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN3>-1) && (fabs(sl[i].an3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS3>-1) && (fabs(sl[i].as3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT3>-1) && (fabs(sl[i].at3) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW3>-1) && (fabs(sl[i].aw3) > nonzeroEPS)) (nna)++;

		if ((sl[i].iB4>-1) && (fabs(sl[i].ab4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iE4>-1) && (fabs(sl[i].ae4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iN4>-1) && (fabs(sl[i].an4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iS4>-1) && (fabs(sl[i].as4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iT4>-1) && (fabs(sl[i].at4) > nonzeroEPS)) (nna)++;
		if ((sl[i].iW4>-1) && (fabs(sl[i].aw4) > nonzeroEPS)) (nna)++;

	}
	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	}

	integer nnu = 0; // число неизвестных.
	nnu = maxelm + maxbound;

	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
	//typedef float    ScalarType;

	// allocate storage for (nnu,nnu) matrix with nna nonzeros

	/*
	cusp::coo_matrix<int, ScalarType, cusp::host_memory> Ah(nnu, nnu, nna);

	*/
	if (bcusp_gl_first_start) {
		Ah_gl.resize(nnu, nnu, nna);
	}
	//cusp::coo_matrix<int, float, cusp::device_memory> Ah(nnu, nnu, nna);

	//printf("0\n");
	//getchar();

	// initialize matrix entries on host
	nna = 0;
	//Ah_gl.row_indices[0] = 0; Ah_gl.column_indices[0] = 0; Ah_gl.values[0] = 10.0; // demo interface
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iP; Ah_gl.values[nna] = sl[i].ap / alpharelax;
			(nna)++;
		}
		if ((sl[i].iB > -1) && (fabs(sl[i].ab) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iB; Ah_gl.values[nna] = -sl[i].ab;
			(nna)++;
		}
		if ((sl[i].iE > -1) && (fabs(sl[i].ae) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iE; Ah_gl.values[nna] = -sl[i].ae;
			(nna)++;
		}
		if ((sl[i].iN > -1) && (fabs(sl[i].an) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iN; Ah_gl.values[nna] = -sl[i].an;
			(nna)++;
		}
		if ((sl[i].iS > -1) && (fabs(sl[i].as) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iS; Ah_gl.values[nna] = -sl[i].as;
			(nna)++;
		}
		if ((sl[i].iT > -1) && (fabs(sl[i].at) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iT; Ah_gl.values[nna] = -sl[i].at;
			(nna)++;
		}
		if ((sl[i].iW > -1) && (fabs(sl[i].aw) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iW; Ah_gl.values[nna] = -sl[i].aw;
			(nna)++;
		}

		if ((sl[i].iB2 > -1) && (fabs(sl[i].ab2) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iB2; Ah_gl.values[nna] = -sl[i].ab2;
			(nna)++;
		}
		if ((sl[i].iE2 > -1) && (fabs(sl[i].ae2) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iE2; Ah_gl.values[nna] = -sl[i].ae2;
			(nna)++;
		}
		if ((sl[i].iN2 > -1) && (fabs(sl[i].an2) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iN2; Ah_gl.values[nna] = -sl[i].an2;
			(nna)++;
		}
		if ((sl[i].iS2 > -1) && (fabs(sl[i].as2) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iS2; Ah_gl.values[nna] = -sl[i].as2;
			(nna)++;
		}
		if ((sl[i].iT2 > -1) && (fabs(sl[i].at2) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iT2; Ah_gl.values[nna] = -sl[i].at2;
			(nna)++;
		}
		if ((sl[i].iW2 > -1) && (fabs(sl[i].aw2) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iW2; Ah_gl.values[nna] = -sl[i].aw2;
			(nna)++;
		}

		if ((sl[i].iB3 > -1) && (fabs(sl[i].ab3) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iB3; Ah_gl.values[nna] = -sl[i].ab3;
			(nna)++;
		}
		if ((sl[i].iE3 > -1) && (fabs(sl[i].ae3) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iE3; Ah_gl.values[nna] = -sl[i].ae3;
			(nna)++;
		}
		if ((sl[i].iN3 > -1) && (fabs(sl[i].an3) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iN3; Ah_gl.values[nna] = -sl[i].an3;
			(nna)++;
		}
		if ((sl[i].iS3 > -1) && (fabs(sl[i].as3) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iS3; Ah_gl.values[nna] = -sl[i].as3;
			(nna)++;
		}
		if ((sl[i].iT3 > -1) && (fabs(sl[i].at3) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iT3; Ah_gl.values[nna] = -sl[i].at3;
			(nna)++;
		}
		if ((sl[i].iW3 > -1) && (fabs(sl[i].aw3) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iW3; Ah_gl.values[nna] = -sl[i].aw3;
			(nna)++;
		}

		if ((sl[i].iB4 > -1) && (fabs(sl[i].ab4) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iB4; Ah_gl.values[nna] = -sl[i].ab4;
			(nna)++;
		}
		if ((sl[i].iE4 > -1) && (fabs(sl[i].ae4) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iE4; Ah_gl.values[nna] = -sl[i].ae4;
			(nna)++;
		}
		if ((sl[i].iN4 > -1) && (fabs(sl[i].an4) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iN4; Ah_gl.values[nna] = -sl[i].an4;
			(nna)++;
		}
		if ((sl[i].iS4 > -1) && (fabs(sl[i].as4) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iS4; Ah_gl.values[nna] = -sl[i].as4;
			(nna)++;
		}
		if ((sl[i].iT4 > -1) && (fabs(sl[i].at4) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iT4; Ah_gl.values[nna] = -sl[i].at4;
			(nna)++;
		}
		if ((sl[i].iW4 > -1) && (fabs(sl[i].aw4) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = sl[i].iP; Ah_gl.column_indices[nna] = sl[i].iW4; Ah_gl.values[nna] = -sl[i].aw4;
			(nna)++;
		}

	}

	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = slb[i].iW; Ah_gl.column_indices[nna] = slb[i].iW; Ah_gl.values[nna] = slb[i].aw;
			(nna)++;
		}
		if ((slb[i].iI > -1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			Ah_gl.row_indices[nna] = slb[i].iW; Ah_gl.column_indices[nna] = slb[i].iI; Ah_gl.values[nna] = -slb[i].ai;
			(nna)++;
		}
	}

	/*
	cusp::array1d<ScalarType, cusp::host_memory> xh(nnu);
	cusp::array1d<ScalarType, cusp::host_memory> yh(nnu);
	*/
	if (bcusp_gl_first_start) {
		xh_gl.resize(nnu);
		yh_gl.resize(nnu);
	}
	//cusp::array1d<float, cusp::device_memory> xh(nnu);
	//cusp::array1d<float, cusp::device_memory> yh(nnu);

	//printf("0i5\n");
	//getchar();

	// yh = Ah x xh;
	for (integer i = 0; i < maxelm + maxbound; i++) {
		// b == dV[i];
		yh_gl[i] = dV[i];
		xh_gl[i] = 0.0;
		//xh[i] = dX0[i];
	}

	/*
	//printf("1\n");
	//getchar();
	// convert host->device
	//cusp::array1d<float, cusp::device_memory> xd = xh;
	cusp::array1d<ScalarType, cusp::device_memory> xd(xh);
	//printf("2\n");
	//getchar();
	//cusp::array1d<float, cusp::device_memory> yd = yh;
	cusp::array1d<ScalarType, cusp::device_memory> yd(yh);
	//printf("3\n");
	//getchar();

	//cusp::coo_matrix<int, float, cusp::device_memory> Ad = Ah;
	cusp::coo_matrix<int, ScalarType, cusp::device_memory> Ad(Ah);
	//printf("4\n");
	//getchar();
	*/
	if (bcusp_gl_first_start) {
		xd_gl.resize(nnu);
		yd_gl.resize(nnu);
		Ad_gl.resize(nnu, nnu, nna);
	}

	cusp::copy(xh_gl, xd_gl);
	cusp::copy(yh_gl, yd_gl);
	cusp::copy(Ah_gl, Ad_gl);

	// set stopping criteria
	// iteration_limit = 100
	// relative_tolerance = 1e-6
	//  absolute_tolerance = 0
	//  verbose            = true
	//cusp::monitor<float> monitor(yd, 2000, 1e-6, 0, true);

	//cusp::default_monitor<float> monitor(yd, 6000, 1e-12);
	// 1e-6 достаточно для стационарных задач.
	//cusp::default_monitor<ScalarType> monitor(yh_gl, 6000, 1e-8);
	cusp::monitor<ScalarType> monitor(yh_gl, 6000, 1e-8);
	//cusp::default_monitor<float> monitor(yd, 2000, 1e-6);

	// setup preconditioner
	//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> Md(Ad);
	// setup preconditioner
	//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::host_memory> Mh(Ah);
	//cusp::precond::aggregation::smoothed_aggregation<int, float, cusp::device_memory> Md(Mh);
	// Диагональный предобуславливатель.
	//cusp::precond::diagonal<float, cusp::device_memory> Md(Ad);
	// AINV (NS Bridson).
	cusp::precond::scaled_bridson_ainv<ScalarType, cusp::device_memory> Md(Ad_gl);


	// solve A * x = y to default tolerance with BiCGStab
	// with preconditioned BiCGStab
	cusp::krylov::bicgstab(Ad_gl, xd_gl, yd_gl, monitor, Md);
	//cusp::krylov::bicgstab(Ad, xd, yd, monitor);
	//cusp::krylov::bicgstab(Ah, xh, yh, monitor, Mh);



	//cusp::array1d<float, cusp::host_memory> xh_ret = xd;
	//cusp::array1d<ScalarType, cusp::host_memory> xh_ret(xd);
	cusp::copy(xd_gl, xh_gl);
	//printf("5\n");
	//getchar();

	// Возвращение результата.
	for (integer i = 0; i < maxelm + maxbound; i++) {
		//dX0[i] = xh_ret[i];
		dX0[i] = xh_gl[i];
	}

	bcusp_gl_first_start = false;

} // cusp_solver_global_allocate

#endif
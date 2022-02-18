
// my_vienna_alg.cu Алгоритмы решения СЛАУ из библиотеки ViennaCL 1.7.1.
// ВНИМАНИЕ!!! На настоящий момент версия библиотеки Vienna1.7.1 имеет статус  БЕТА.
// Последний выпуск 20 января 2016 года.
#pragma once


// ublas headers
/*
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
*/

// ViennaCL 1.7.1 Includes
//#include <vector>
//#include <cmath>

//#include "viennacl/forwards.h"

#include "viennacl/vector.hpp"
#include "viennacl/coordinate_matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/tools/matrix_generation.hpp"
/*
#ifndef _In_
#define _In_
#endif

#ifndef _Out_
#define _Out_
#endif

#ifndef PVOID
#define PVOID
#endif

#ifndef _In_opt_
#define _In_opt_
#endif

*/


//#include "viennacl/tools/tools.hpp"
//#include "viennacl/linalg/sparse_matrix_operations.hpp"
//#include "viennacl/linalg/row_scaling.hpp"
//#include "viennacl/ocl/kernel.hpp"
//#include "viennacl/ocl/platform.hpp"
//#include "viennacl/ocl/utils.hpp"
//#include "viennacl/linalg/spai.hpp"


//#include <map>

// Boost headers:
/*
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/io.hpp"
*/

// device info
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/device_specific/builtin_database/common.hpp"

/**
* Import the AMG functionality:
**/
#include "viennacl/linalg/amg.hpp"

/**
* Some more includes:
**/
#include <iostream>
#include <vector>
#include <ctime>
#include "vector-io.hpp"
#include "viennacl/tools/timer.hpp"

// Старейшая библиотека hypre содержащая Boomer amg. 18.01.2017.
//#include "hypre/utilities/_hypre_utilities.h"
//#include "hypre/krylov/krylov.h"
//#include "hypre/HYPRE.h"

//#include "krylov.h"
//#include "HYPRE.h"


//Chow-Pattel ilu
// Параллельное ilu разложение.
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/detail/ilu/common.hpp"
#include "viennacl/linalg/ilu_operations.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/backend/memory.hpp"





// ViennaCL open code:
/** <h2>Part 1: Worker routines</h2>
*
*  <h3>Run the Solver</h3>
*   Runs the provided solver specified in the `solver` object with the provided preconditioner `precond`
**/
template<typename MatrixType, typename VectorType, typename SolverTag, typename PrecondTag>
void run_solver(MatrixType const & matrix, VectorType const & rhs, VectorType & ref_result, SolverTag const & solver, PrecondTag const & precond)
{
	VectorType result(rhs);
	//VectorType my_result(rhs);
	VectorType residual(rhs);
	VectorType xo(ref_result);

	std::cout <<  "rhs, x0:" << viennacl::linalg::norm_2(residual) << "  " << viennacl::linalg::norm_2(xo)  << std::endl;
	xo= viennacl::linalg::prod(matrix, result);
	std::cout << "test matrix product: "<< viennacl::linalg::norm_2(xo) << std::endl;

	viennacl::tools::timer timer;
	timer.start();
	result = viennacl::linalg::solve(matrix, rhs, solver, precond);
	//my_result = result;
	viennacl::backend::finish();
	std::cout << "  > Solver time: " << timer.get() << std::endl;
	residual -= viennacl::linalg::prod(matrix, result);
	std::cout << "  > Relative residual: res rhs res/rhs" << viennacl::linalg::norm_2(residual) << "  " << viennacl::linalg::norm_2(rhs) << "  " << viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(rhs) << std::endl;
	std::cout << "  > Iterations: " << solver.iters() << std::endl;
	//result -= ref_result;
	ref_result = result; // Возвращение результата вычисления.
	std::cout << "  > Relative deviation from result: " << viennacl::linalg::norm_2(result) / viennacl::linalg::norm_2(ref_result) << std::endl;
}

/** <h3>Compare AMG preconditioner for uBLAS and ViennaCL types</h3>
*
*  The AMG implementations in ViennaCL can be used with uBLAS types as well as ViennaCL types.
*  This function compares the two in terms of execution time.
**/
template<typename ScalarType>
void run_amg(viennacl::linalg::bicgstab_tag & bicgstab_solver,
	viennacl::vector<ScalarType> & vcl_vec,
	viennacl::vector<ScalarType> & vcl_result,
	viennacl::compressed_matrix<ScalarType> & vcl_compressed_matrix,
	std::string info,
	viennacl::linalg::amg_tag & amg_tag)
{
	std::cout << "-- BiCGStab with AMG preconditioner, " << info << " --" << std::endl;

	viennacl::linalg::amg_precond<viennacl::compressed_matrix<ScalarType> > vcl_amg(vcl_compressed_matrix, amg_tag);
	std::cout << " * Setup phase (ViennaCL types)..." << std::endl;
	viennacl::tools::timer timer;
	timer.start();
	vcl_amg.setup();
	viennacl::backend::finish();
	std::cout << "  > Setup time: " << timer.get() << std::endl;

	std::cout << " * BiCGStab solver (ViennaCL types)..." << std::endl;
	run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, bicgstab_solver, vcl_amg);
}

/** <h3>Compare AMG preconditioner for uBLAS and ViennaCL types</h3>
*
*  The AMG implementations in ViennaCL can be used with uBLAS types as well as ViennaCL types.
*  This function compares the two in terms of execution time.
**/
template<typename ScalarType>
void run_amg(viennacl::linalg::gmres_tag & gmresm_solver,
	viennacl::vector<ScalarType> & vcl_vec,
	viennacl::vector<ScalarType> & vcl_result,
	viennacl::compressed_matrix<ScalarType> & vcl_compressed_matrix,
	std::string info,
	viennacl::linalg::amg_tag & amg_tag)
{
	std::cout << "-- FGMRES with AMG preconditioner, " << info << " --" << std::endl;

	viennacl::linalg::amg_precond<viennacl::compressed_matrix<ScalarType> > vcl_amg(vcl_compressed_matrix, amg_tag);
	std::cout << " * Setup phase (ViennaCL types)..." << std::endl;
	viennacl::tools::timer timer;
	timer.start();
	vcl_amg.setup();
	viennacl::backend::finish();
	std::cout << "  > Setup time: " << timer.get() << std::endl;

	std::cout << " * FGMRES solver (ViennaCL types)..." << std::endl;
	run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, gmresm_solver, vcl_amg);
}

/** <h3>Compare AMG preconditioner for uBLAS and ViennaCL types</h3>
*
*  The AMG implementations in ViennaCL can be used with uBLAS types as well as ViennaCL types.
*  This function compares the two in terms of execution time.
**/
template<typename ScalarType>
void run_amg(viennacl::linalg::cg_tag & cg_solver,
	viennacl::vector<ScalarType> & vcl_vec,
	viennacl::vector<ScalarType> & vcl_result,
	viennacl::compressed_matrix<ScalarType> & vcl_compressed_matrix,
	std::string info,
	viennacl::linalg::amg_tag & amg_tag)
{
	std::cout << "-- CG with AMG preconditioner, " << info << " --" << std::endl;

	viennacl::linalg::amg_precond<viennacl::compressed_matrix<ScalarType> > vcl_amg(vcl_compressed_matrix, amg_tag);
	std::cout << " * Setup phase (ViennaCL types)..." << std::endl;
	viennacl::tools::timer timer;
	timer.start();
	vcl_amg.setup();
	viennacl::backend::finish();
	std::cout << "  > Setup time: " << timer.get() << std::endl;
	
	

	std::cout << " * CG solver (ViennaCL types)..." << std::endl;
	run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, cg_solver, vcl_amg);
}



// Это вызов библиотечного решателя систем линейных алгебраических уравнений.
// Библиотека ViennaCL версии 1.7.1. Это библиотека с открытым исходным кодом распространяющаяся 
// по open Source лицензии MIT (X11) license. 
// Начало присоединения 15 октября 2016.
void viennacl_solver(equation3D* &sl, equation3D_bon* &slb,
	integer maxelm, integer maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit,
	doublereal alpharelax, integer iVar)
{


	// Для организации параллельных вычислений на центральном процессоре необходимо 
	// прописать #define VIENNACL_WITH_OPENMP 1 перед подключением заголовков VIENNACL.
	// Добавить опцию /openmp в параметры visual Studio.
	// использовать не ilu0 предобуславливание, а специальное распараллеленное chow_patel_ilu
	// предобуславливание. Стандарт ilu0 не подразумевает распараллеливания.

	// maxit,  iVar - не используется.




	if (dX0 == NULL) {
		dX0 = new doublereal[maxelm + maxbound];
		for (integer i = 0; i < maxelm + maxbound; ++i) {
			dX0[i] = 0.0;
		}
	}

	// TODO получить val, col_ind, row_ptr
	int nna = 0; // количество ненулевых элементов в матрице СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	// подсчёт числа ненулевых элементов в матрице.
	nna = 0;
	for (integer i = 0; i<maxelm; ++i) {
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
	for (integer i = 0; i<maxbound; ++i) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	}

	int nnu = 0; // число неизвестных.
	nnu = (int)(maxelm + maxbound);

	/**
	* Printeger some device info at the beginning. If there is more than one OpenCL device available, use the second device.
	**/
	std::cout << std::endl;
	std::cout << "----------------------------------------------" << std::endl;
	std::cout << "               Device Info" << std::endl;
	std::cout << "----------------------------------------------" << std::endl;

	/*
	typedef std::vector< viennacl::ocl::platform > platforms_type;
	platforms_type platforms = viennacl::ocl::get_platforms();

	bool is_first_element = true;
	for (platforms_type::iterator platform_iter = platforms.begin();
	platform_iter != platforms.end();
	++platform_iter)
	{
	typedef std::vector< viennacl::ocl::device> devices_type;
	devices_type devices = platform_iter->devices(CL_DEVICE_TYPE_ALL);


	}
	*/

#ifdef VIENNACL_WITH_OPENCL
	// Optional: Customize OpenCL backend
	viennacl::ocl::platform pf = viennacl::ocl::get_platforms()[0];
	std::vector<viennacl::ocl::device> const & devices = pf.devices();

	// Optional: Set first device to first context:
	viennacl::ocl::setup_context(0, devices[0]);

	// Optional: Set second device for second context (use the same device for the second context if only one device available):
	if (devices.size() > 1)
		viennacl::ocl::setup_context(1, devices[1]);
	else
		viennacl::ocl::setup_context(1, devices[0]);

	std::cout << viennacl::ocl::current_device().info() << std::endl;
	viennacl::context ctx(viennacl::ocl::get_context(0));
#else

#ifdef VIENNACL_WITH_CUDA

	// Optional: Customize OpenCL backend
	//viennacl::ocl::platform pf = viennacl::ocl::get_platforms()[0];
	//std::vector<viennacl::ocl::device> const & devices = pf.devices();

	std::vector< viennacl::ocl::device > devices =
		viennacl::ocl::platform().devices();

	//std::vector< viennacl::ocl::device > my_devices;
	//my_devices.push_back(devices[0]);
	//my_devices.push_back(devices[1]);

	// Optional: Set first device to first context:
	//viennacl::ocl::setup_context(1, devices);

	// Optional: Set second device for second context (use the same device for the second context if only one device available):
	/*
	if (devices.size() > 1)
	viennacl::ocl::setup_context(1, devices[1]);
	else
	viennacl::ocl::setup_context(1, devices[0]);
	*/

	std::cout << viennacl::ocl::current_device().info() << std::endl;


	//viennacl::context ctx(viennacl::ocl::get_context(0));

	viennacl::context ctx;
#else
	viennacl::context ctx;
	//std::size_t tid = 8;
	//std::size_t thread_id_ = tid;
	//viennacl::context ctx(viennacl::ocl::get_context(static_cast<long>(thread_id_)));

#endif

#endif



	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
	//typedef float    ScalarType;


	/**
	* Set up the matrices and vectors for the iterative solvers (cf. iterative.cpp)
	**/
	viennacl::vcl_size_t rows = nnu;
	viennacl::vcl_size_t cols = nnu;
	viennacl::vcl_size_t nonzeros = nna;


	// Прочитать матрицу:
	std::cout << "Reading matrix..." << std::endl;
	//std::vector< std::map<unsigned int, ScalarType> > vcl_compressed_matrix1;// read_in_matrix;
	//boost::numeric::ublas::matrix<ScalarType> vcl_compressed_matrix1(rows, cols);// read_in_matrix;

	//const void* row_jumper = new integer[nnu + 1];
	//const void* col_buffer = new integer[nonzeros];
	//const ScalarType* elements = new ScalarType[nonzeros];
	int* row_jumper = new int[nnu + 1];
	int* col_buffer = new int[nna];
	ScalarType* elements = new ScalarType[nna];

	row_jumper[nnu] = nna;
	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	for (int i = 0; i<(int)(maxelm); ++i) {
		row_jumper[i] = nna;

		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iP; Ah.values[nna] = sl[i].ap / alpharelax;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iP) = static_cast<ScalarType> (sl[i].ap / alpharelax);
			col_buffer[nna] = (int)(sl[i].iP);
			elements[nna] = static_cast<ScalarType> (sl[i].ap / alpharelax);
			(nna)++;
		}
		if ((sl[i].iB > -1) && (fabs(sl[i].ab) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB; Ah.values[nna] = -sl[i].ab;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB) = static_cast<ScalarType> (-sl[i].ab);
			col_buffer[nna] = (int)(sl[i].iB);
			elements[nna] = static_cast<ScalarType> (-sl[i].ab);
			(nna)++;
		}
		if ((sl[i].iE > -1) && (fabs(sl[i].ae) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE; Ah.values[nna] = -sl[i].ae;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE) = static_cast<ScalarType> (-sl[i].ae);
			col_buffer[nna] = (int)(sl[i].iE);
			elements[nna] = static_cast<ScalarType> (-sl[i].ae);
			(nna)++;
		}
		if ((sl[i].iN > -1) && (fabs(sl[i].an) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN; Ah.values[nna] = -sl[i].an;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN) = static_cast<ScalarType> (-sl[i].an);
			col_buffer[nna] = (int)(sl[i].iN);
			elements[nna] = static_cast<ScalarType> (-sl[i].an);
			(nna)++;
		}
		if ((sl[i].iS > -1) && (fabs(sl[i].as) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS; Ah.values[nna] = -sl[i].as;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS) = static_cast<ScalarType> (-sl[i].as);
			col_buffer[nna] = (int)(sl[i].iS);
			elements[nna] = static_cast<ScalarType> (-sl[i].as);
			(nna)++;
		}
		if ((sl[i].iT > -1) && (fabs(sl[i].at) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT; Ah.values[nna] = -sl[i].at;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT) = static_cast<ScalarType> (-sl[i].at);
			col_buffer[nna] = (int)(sl[i].iT);
			elements[nna] = static_cast<ScalarType> (-sl[i].at);
			(nna)++;
		}
		if ((sl[i].iW > -1) && (fabs(sl[i].aw) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW; Ah.values[nna] = -sl[i].aw;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW) = static_cast<ScalarType> (-sl[i].aw);
			col_buffer[nna] = (int)(sl[i].iW);
			elements[nna] = static_cast<ScalarType> (-sl[i].aw);
			(nna)++;
		}

		if ((sl[i].iB2 > -1) && (fabs(sl[i].ab2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB2; Ah.values[nna] = -sl[i].ab2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB2) = static_cast<ScalarType> (-sl[i].ab2);
			col_buffer[nna] = (int)(sl[i].iB2);
			elements[nna] = static_cast<ScalarType> (-sl[i].ab2);
			(nna)++;
		}
		if ((sl[i].iE2 > -1) && (fabs(sl[i].ae2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE2; Ah.values[nna] = -sl[i].ae2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE2) = static_cast<ScalarType> (-sl[i].ae2);
			col_buffer[nna] = (int)(sl[i].iE2);
			elements[nna] = static_cast<ScalarType> (-sl[i].ae2);
			(nna)++;
		}
		if ((sl[i].iN2 > -1) && (fabs(sl[i].an2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN2; Ah.values[nna] = -sl[i].an2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN2) = static_cast<ScalarType> (-sl[i].an2);
			col_buffer[nna] = (int)(sl[i].iN2);
			elements[nna] = static_cast<ScalarType> (-sl[i].an2);
			(nna)++;
		}
		if ((sl[i].iS2 > -1) && (fabs(sl[i].as2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS2; Ah.values[nna] = -sl[i].as2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS2) = static_cast<ScalarType>( -sl[i].as2);
			col_buffer[nna] = (int)(sl[i].iS2);
			elements[nna] = static_cast<ScalarType> (-sl[i].as2);
			(nna)++;
		}
		if ((sl[i].iT2 > -1) && (fabs(sl[i].at2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT2; Ah.values[nna] = -sl[i].at2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT2) = static_cast<ScalarType> (-sl[i].at2);
			col_buffer[nna] = (int)(sl[i].iT2);
			elements[nna] = static_cast<ScalarType> (-sl[i].at2);
			(nna)++;
		}
		if ((sl[i].iW2 > -1) && (fabs(sl[i].aw2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW2; Ah.values[nna] = -sl[i].aw2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW2) = static_cast<ScalarType> (-sl[i].aw2);
			col_buffer[nna] = (int)(sl[i].iW2);
			elements[nna] = static_cast<ScalarType> (-sl[i].aw2);
			(nna)++;
		}

		if ((sl[i].iB3 > -1) && (fabs(sl[i].ab3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB3; Ah.values[nna] = -sl[i].ab3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB3) = static_cast<ScalarType>(-sl[i].ab3);
			col_buffer[nna] = (int)(sl[i].iB3);
			elements[nna] = static_cast<ScalarType> (-sl[i].ab3);
			(nna)++;
		}
		if ((sl[i].iE3 > -1) && (fabs(sl[i].ae3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE3; Ah.values[nna] = -sl[i].ae3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE3) = static_cast<ScalarType>(-sl[i].ae3);
			col_buffer[nna] = (int)(sl[i].iE3);
			elements[nna] = static_cast<ScalarType> (-sl[i].ae3);
			(nna)++;
		}
		if ((sl[i].iN3 > -1) && (fabs(sl[i].an3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN3; Ah.values[nna] = -sl[i].an3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN3) = static_cast<ScalarType> (-sl[i].an3);
			col_buffer[nna] = (int)(sl[i].iN3);
			elements[nna] = static_cast<ScalarType> (-sl[i].an3);
			(nna)++;
		}
		if ((sl[i].iS3 > -1) && (fabs(sl[i].as3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS3; Ah.values[nna] = -sl[i].as3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS3) = static_cast<ScalarType>(-sl[i].as3);
			col_buffer[nna] = (int)(sl[i].iS3);
			elements[nna] = static_cast<ScalarType> (-sl[i].as3);
			(nna)++;
		}
		if ((sl[i].iT3 > -1) && (fabs(sl[i].at3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT3; Ah.values[nna] = -sl[i].at3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT3) = static_cast<ScalarType>(-sl[i].at3);
			col_buffer[nna] = (int)(sl[i].iT3);
			elements[nna] = static_cast<ScalarType> (-sl[i].at3);
			(nna)++;
		}
		if ((sl[i].iW3 > -1) && (fabs(sl[i].aw3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW3; Ah.values[nna] = -sl[i].aw3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW3) = static_cast<ScalarType>(-sl[i].aw3);
			col_buffer[nna] = (int)(sl[i].iW3);
			elements[nna] = static_cast<ScalarType> (-sl[i].aw3);
			(nna)++;
		}

		if ((sl[i].iB4 > -1) && (fabs(sl[i].ab4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB4; Ah.values[nna] = -sl[i].ab4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB4) = static_cast<ScalarType>(-sl[i].ab4);
			col_buffer[nna] = (int)(sl[i].iB4);
			elements[nna] = static_cast<ScalarType> (-sl[i].ab4);
			(nna)++;
		}
		if ((sl[i].iE4 > -1) && (fabs(sl[i].ae4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE4; Ah.values[nna] = -sl[i].ae4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE4) = static_cast<ScalarType>(-sl[i].ae4);
			col_buffer[nna] = (int)(sl[i].iE4);
			elements[nna] = static_cast<ScalarType> (-sl[i].ae4);
			(nna)++;
		}
		if ((sl[i].iN4 > -1) && (fabs(sl[i].an4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN4; Ah.values[nna] = -sl[i].an4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN4) = static_cast<ScalarType>(-sl[i].an4);
			col_buffer[nna] = (int)(sl[i].iN4);
			elements[nna] = static_cast<ScalarType> (-sl[i].an4);
			(nna)++;
		}
		if ((sl[i].iS4 > -1) && (fabs(sl[i].as4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS4; Ah.values[nna] = -sl[i].as4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS4) = static_cast<ScalarType>(-sl[i].as4);
			col_buffer[nna] = (int)(sl[i].iS4);
			elements[nna] = static_cast<ScalarType> (-sl[i].as4);
			(nna)++;
		}
		if ((sl[i].iT4 > -1) && (fabs(sl[i].at4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT4; Ah.values[nna] = -sl[i].at4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT4) = static_cast<ScalarType>(-sl[i].at4);
			col_buffer[nna] = (int)(sl[i].iT4);
			elements[nna] = static_cast<ScalarType> (-sl[i].at4);
			(nna)++;
		}
		if ((sl[i].iW4 > -1) && (fabs(sl[i].aw4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW4; Ah.values[nna] = -sl[i].aw4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW4) = static_cast<ScalarType>(-sl[i].aw4);
			col_buffer[nna] = (int)(sl[i].iW4);
			elements[nna] = static_cast<ScalarType> (-sl[i].aw4);
			(nna)++;
		}


	}

	for (int i = 0; i< (int)(maxbound); ++i) {
		row_jumper[maxelm + i] = nna;
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) {
			//Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iW; Ah.values[nna] = slb[i].aw;
			//vcl_compressed_matrix1(slb[i].iW, slb[i].iW) = static_cast<ScalarType> (slb[i].aw);
			col_buffer[nna] = (int)(slb[i].iW);
			elements[nna] = static_cast<ScalarType> (slb[i].aw);
			(nna)++;
		}
		if ((slb[i].iI > -1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			//Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iI; Ah.values[nna] = -slb[i].ai;
			//vcl_compressed_matrix1(slb[i].iW, slb[i].iI) = static_cast<ScalarType>(-slb[i].ai);
			col_buffer[nna] = (int)(slb[i].iI);
			elements[nna] = static_cast<ScalarType> (-slb[i].ai);
			(nna)++;
		}
	}

	//viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(rows, cols, nonzeros, ctx);
	viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(rows, cols, nonzeros, ctx);
	//viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(row_jumper, col_buffer, elements, rows, cols, nonzeros);
	vcl_compressed_matrix.set(row_jumper, col_buffer, elements, rows, cols, nonzeros);

	// free
	delete[] row_jumper;
	delete[] col_buffer;
	delete[] elements;

	//viennacl::copy(read_in_matrix, vcl_compressed_matrix);
	//viennacl::copy(vcl_compressed_matrix1, vcl_compressed_matrix);

	std::cout << "Reading matrix completed." << std::endl;

	// Вектора:
	//viennacl::vector<ScalarType> vcl_vec(vcl_compressed_matrix.size1(), ctx);
	//viennacl::vector<ScalarType> vcl_result(vcl_compressed_matrix.size1(), ctx);
	viennacl::vector<ScalarType> vcl_vec(rows, ctx);
	viennacl::vector<ScalarType> vcl_result(rows, ctx);

	std::vector<ScalarType> std_vec, std_result;

	// rhs and result vector:
	//std_vec.resize(vcl_compressed_matrix.size1());
	//std_result.resize(vcl_compressed_matrix.size1());
	std_vec.resize(rows);
	std_result.resize(rows);
	for (std::size_t i = 0; i < std_result.size(); ++i) {
		//std_result[i] = ScalarType(1);
		std_result[i] = dX0[i];
		std_vec[i] = dV[i];
	}

	// Copy to GPU
	viennacl::copy(std_vec, vcl_vec);
	viennacl::copy(std_result, vcl_result);

	//vcl_vec = viennacl::linalg::prod(vcl_compressed_matrix, vcl_result);





	/**
	* Instantiate a tag for the conjugate gradient solver, the AMG preconditioner tag, and create an AMG preconditioner object:
	**/
	viennacl::linalg::bicgstab_tag bicgstab_solver(1e-6, 2000); // 1e-8
	//viennacl::linalg::cg_tag cg_solver(1e-6, 10000);

	viennacl::context host_ctx(viennacl::MAIN_MEMORY);
	//viennacl::context target_ctx = viennacl::traits::context(vcl_compressed_matrix);
	viennacl::context target_ctx(viennacl::MAIN_MEMORY);

	/**
	* Run solver without preconditioner. This serves as a baseline for comparison.
	* Note that iterative solvers without preconditioner on GPUs can be very efficient because they map well to the massively parallel hardware.
	**/
	//std::cout << "-- BiCGStab solver (no preconditioner, warmup) --" << std::endl;
	//run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, bicgstab_solver, viennacl::linalg::no_precond());

	if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 5) {
		/**
		* Generate the setup for an AMG preconditioner of Ruge-Stueben type with only one pass and direct interpolation (ONEPASS+DIRECT)
		**/
		/*
		viennacl::linalg::amg_tag amg_tag_direct;
		amg_tag_direct.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_ONEPASS);
		amg_tag_direct.set_interpolation_method(viennacl::linalg::AMG_INTERPOLATION_METHOD_DIRECT);
		amg_tag_direct.set_strong_connection_threshold(0.25);
		amg_tag_direct.set_jacobi_weight(0.67);
		amg_tag_direct.set_presmooth_steps(1);
		amg_tag_direct.set_postsmooth_steps(1);
		amg_tag_direct.set_setup_context(host_ctx);    // run setup on host
		amg_tag_direct.set_target_context(target_ctx); // run solver cycles on device
		run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "ONEPASS COARSENING, DIRECT INTERPOLATION", amg_tag_direct);
		*/
		/**
		* Generate the setup for a smoothed aggregation AMG preconditioner
		**/
		viennacl::linalg::amg_tag amg_tag_sa_pmis;
		amg_tag_sa_pmis.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_MIS2_AGGREGATION);
		amg_tag_sa_pmis.set_interpolation_method(viennacl::linalg::AMG_INTERPOLATION_METHOD_SMOOTHED_AGGREGATION);
		// Минимальное количество ячеек сетки на самом грубом уровне должно быть меньше этой величины coarse_cutoff_(50);

		// cutoff очень большое, агрегация проходит плохо, сходимости нет. Тестировал на модели tgf01.
		amg_tag_sa_pmis.set_coarsening_cutoff(250);
		amg_tag_sa_pmis.set_strong_connection_threshold(0.25);
		std::cout << " threshold=" << amg_tag_sa_pmis.get_strong_connection_threshold() << std::endl;
		//system("PAUSE");
		//amg_tag_sa_pmis.set_strong_connection_threshold(0.4);//0.1; 0.25
		amg_tag_sa_pmis.set_jacobi_weight(0.6667);
		amg_tag_sa_pmis.set_presmooth_steps(2);
		amg_tag_sa_pmis.set_postsmooth_steps(2);
		//viennacl::linalg::gmres_tag gmresm_solver(1e-6, 100, my_amg_manager.m_restart); // 1e-6 1e-8
		run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), SA INTERPOLATION", amg_tag_sa_pmis);
		//run_amg(gmresm_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), SA INTERPOLATION", amg_tag_sa_pmis);
		//run_amg(cg_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), SA INTERPOLATION", amg_tag_sa_pmis);

		printf("number iteration %d finish residual %e\n", bicgstab_solver.iters(), bicgstab_solver.error());
	}
	else if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 9) {

		printf("Vienna CL 1.7.1 algorithms collection.\n");
		printf("ILU0: ilu0 preconditioner + bicgstab, ilu0 zero fill in\n");

		//ILU0
		viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<ScalarType> > vcl_ilu0(vcl_compressed_matrix, viennacl::linalg::ilu0_tag());
		//viennacl::linalg::block_ilu_precond< viennacl::compressed_matrix<ScalarType>, viennacl::linalg::ilu0_tag > vcl_ilu0(vcl_compressed_matrix, viennacl::linalg::ilu0_tag());
		//viennacl::linalg::jacobi_precond<viennacl::compressed_matrix<ScalarType>> vcl_jacobi(vcl_compressed_matrix, viennacl::linalg::jacobi_tag);
		// using viennacl objects on GPU
		//vcl_result = solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::bicgstab_tag(2000, 1e-6), vcl_ilu0);

		printf("ILU0 preconditioner found.\n");

		//bicgstab
		viennacl::linalg::bicgstab_tag mybicgstab_solver(1e-6, 10000); // 1e-8

		//vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::bicgstab_tag(1e-6, 2000), vcl_ilu0); //with preconditioner
		vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, mybicgstab_solver, vcl_ilu0); //with preconditioner
		//vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::bicgstab_tag(1e-6, 2000), vcl_ilu0); //with preconditioner

		//vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, mybicgstab_solver, vcl_jacobi); //with preconditioner
		//gmres_tag(1e-6, 2000, 30)
		//vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::gmres_tag(1e-6, 2000, 30), vcl_ilu0); //with preconditioner

		printf("number iteration %d finish residual %e\n", mybicgstab_solver.iters(), mybicgstab_solver.error());
	}
	else if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 10) {

		//printf("Vienna CL 1.7.1 algorithms collection.\n");
		//printf("ILUT: ilu threshold preconditioner + bicgstab, ilut(p,tau)\n");
		
		//viennacl::linalg::ilut_tag ilut_conf(10, 1e-5, false); // 10 entries, rel. tol. 1e-5
		// расходится.
		//viennacl::linalg::ilut_tag ilut_conf(20, 1e-4); // 20 entries (всего 40 для L и U совместно), rel. tol. 1e-4
		//--->viennacl::linalg::ilut_tag ilut_conf(200, 1e-8, false); // 40 entries (all 80), rel. tol. 1e-4
		//typedef viennacl::linalg::ilut_precond<viennacl::compressed_matrix<ScalarType> > vcl_ilut_t;

		//vcl_ilut_t vcl_ilut(vcl_compressed_matrix, ilut_conf);
		/*
		//---->viennacl::linalg::ilut_precond< viennacl::compressed_matrix<ScalarType> > vcl_ilut(vcl_compressed_matrix, viennacl::linalg::ilut_tag(200, 1e-8, false));
		viennacl::linalg::ilut_precond< viennacl::compressed_matrix<ScalarType> > vcl_ilut(vcl_compressed_matrix, viennacl::linalg::ilut_tag());


		// using viennacl objects on GPU
		vcl_result = solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::bicgstab_tag(1e-6, 2000), vcl_ilut);
		
		//vcl_result = solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::gmres_tag(1e-6, 2000, 30), vcl_ilut);
		*/
		
		
		// setup SPAI preconditioner, purely CPU-based
		//viennacl::linalg::chow_patel_icc_precond< viennacl::compressed_matrix<ScalarType> > spai_cpu(vcl_compressed_matrix, viennacl::linalg::spai_tag(1.0e-3, 3, 5e-2));
		//viennacl::linalg::spai_precond< viennacl::compressed_matrix<ScalarType> > spai_cpu(vcl_compressed_matrix, viennacl::linalg::spai_tag(1.0e-3, 3, 5e-2));


		// solve (e.g. using stab, Bi-conjugate gradient solver)
		//vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::bicgstab_tag(1e-6, 2000), spai_cpu);
		

		/*
		// Chow - Pattel ilu0 параллельное ilu0 разложение.
		// chow_patel_tag(vcl_size_t num_sweeps=3, vcl_size_t num_jacobi_iters=2)
		// Нет сходимости. Тестировался на 10мм tgf. Видимо это очень слабый предобуславливатель.
		viennacl::linalg::chow_patel_tag chow_patel_ilu_config;
		chow_patel_ilu_config.sweeps(3); //20 five nonlinear sweeps
		chow_patel_ilu_config.jacobi_iters(2); //19 Jacobi iterations per triangular 'solve' Rx=r

		// create and compute preconditioner:
		viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<ScalarType> > chow_patel_ilu(vcl_compressed_matrix, chow_patel_ilu_config);
		// solve (e.g. using stab, Bi-conjugate gradient solver)
		vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::bicgstab_tag(1e-6, 2000), chow_patel_ilu);
		*/


		printf("Vienna CL 1.7.1 algorithms collection.\n");
		printf("ILU0: ilu0 preconditioner + gmres(my_amg_manager.memory_size_Stress), ilu0 zero fill in\n");

		viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<ScalarType> > vcl_ilu0(vcl_compressed_matrix, viennacl::linalg::ilu0_tag());
		//gmres_tag(1e-6, 2000, my_amg_manager.m_restart)

		//my_amg_manager.memory_size_Stress
		viennacl::linalg::gmres_tag gmresm_solver(1e-6, 10000, my_amg_manager.m_restart); // 1e-6 1e-8
		//vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::gmres_tag(1e-6, 2000, my_amg_manager.memory_size_Stress), vcl_ilu0); //with preconditioner
		vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, gmresm_solver, vcl_ilu0); //with preconditioner


		printf("number iteration %d finish residual %e\n", gmresm_solver.iters(), gmresm_solver.error());

	}

	/**
	* Generate the setup for an aggregation-based AMG preconditioner with unsmoothed aggregation
	**/
	/*
	viennacl::linalg::amg_tag amg_tag_agg_pmis;
	amg_tag_agg_pmis.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_MIS2_AGGREGATION);
	amg_tag_agg_pmis.set_interpolation_method(viennacl::linalg::AMG_INTERPOLATION_METHOD_AGGREGATION);
	run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), AG INTERPOLATION", amg_tag_agg_pmis);
	*/

	/**
	* Generate the setup for a smoothed aggregation AMG preconditioner
	**/
	/*
	viennacl::linalg::amg_tag amg_tag_sa_pmis;
	amg_tag_sa_pmis.set_coarsening_method(viennacl::cuda::linalg::AMG_COARSENING_METHOD_MIS2_AGGREGATION);
	amg_tag_sa_pmis.set_interpolation_method(viennacl::cuda::linalg::AMG_INTERPOLATION_METHOD_SMOOTHED_AGGREGATION);
	run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), SA INTERPOLATION", amg_tag_sa_pmis);
	*/

	//std::cout << std::endl;
	//std::cout << " -------------- Benchmark runs -------------- " << std::endl;
	//std::cout << std::endl;

	//std::cout << "-- BiCGStab solver (no preconditioner) --" << std::endl;
	//run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, cg_solver, viennacl::linalg::no_precond());
	//run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "ONEPASS COARSENING, DIRECT INTERPOLATION", amg_tag_direct);
	//run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), AG INTERPOLATION", amg_tag_agg_pmis);
	//run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), SA INTERPOLATION", amg_tag_sa_pmis);

	/**
	*  That's it.
	**/
	std::cout << "!!!! CALCULATION COMPLETED SUCCESSFULLY !!!!" << std::endl;

	// Обратное копирование:
	viennacl::copy(vcl_result, std_result);

	for (std::size_t i = 0; i < rows; ++i) {
		// Возвращение результата вычисления.
		dX0[i] = std_result[i];
		//dX0[i] = vcl_result(i);
	}
}

// Это вызов библиотечного решателя систем линейных алгебраических уравнений.
// Библиотека ViennaCL версии 1.7.1. Это библиотека с открытым исходным кодом распространяющаяся 
// по open Source лицензии MIT (X11) license. 
// Начало присоединения 15 октября 2016.
// Здесь оставлен однопоточный работающий вариант.
void viennacl_solver_serial(equation3D* &sl, equation3D_bon* &slb,
	integer maxelm, integer maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit,
	doublereal alpharelax, integer iVar)
{

	// maxit,  iVar - не используется.




	if (dX0 == NULL) {
		dX0 = new doublereal[maxelm + maxbound];
		for (integer i = 0; i < maxelm + maxbound; ++i) {
			dX0[i] = 0.0;
		}
	}

	// TODO получить val, col_ind, row_ptr
	integer nna = 0; // количество ненулевых элементов в матрице СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	// подсчёт числа ненулевых элементов в матрице.
	nna = 0;
	for (integer i = 0; i<maxelm; ++i) {
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
	for (integer i = 0; i<maxbound; ++i) {
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	}

	integer nnu = 0; // число неизвестных.
	nnu = maxelm + maxbound;

	/**
	* Printeger some device info at the beginning. If there is more than one OpenCL device available, use the second device.
	**/
	std::cout << std::endl;
	std::cout << "----------------------------------------------" << std::endl;
	std::cout << "               Device Info" << std::endl;
	std::cout << "----------------------------------------------" << std::endl;

	/*
	typedef std::vector< viennacl::ocl::platform > platforms_type;
	platforms_type platforms = viennacl::ocl::get_platforms();

	bool is_first_element = true;
	for (platforms_type::iterator platform_iter = platforms.begin();
	platform_iter != platforms.end();
	++platform_iter)
	{
	typedef std::vector< viennacl::ocl::device> devices_type;
	devices_type devices = platform_iter->devices(CL_DEVICE_TYPE_ALL);
	

	}
	*/

	/*
#ifdef VIENNACL_WITH_OPENCL
	// Optional: Customize OpenCL backend
	viennacl::ocl::platform pf = viennacl::ocl::get_platforms()[0];
	std::vector<viennacl::ocl::device> const & devices = pf.devices();

	// Optional: Set first device to first context:
	viennacl::ocl::setup_context(0, devices[0]);

	// Optional: Set second device for second context (use the same device for the second context if only one device available):
	if (devices.size() > 1)
		viennacl::ocl::setup_context(1, devices[1]);
	else
		viennacl::ocl::setup_context(1, devices[0]);

	std::cout << viennacl::ocl::current_device().info() << std::endl;
	viennacl::context ctx(viennacl::ocl::get_context(0));
#else

#ifdef VIENNACL_WITH_CUDA

	// Optional: Customize OpenCL backend
	//viennacl::ocl::platform pf = viennacl::ocl::get_platforms()[0];
	//std::vector<viennacl::ocl::device> const & devices = pf.devices();

	std::vector< viennacl::ocl::device > devices =
		viennacl::ocl::platform().devices();

	//std::vector< viennacl::ocl::device > my_devices;
	//my_devices.push_back(devices[0]);
	//my_devices.push_back(devices[1]);

	// Optional: Set first device to first context:
	//viennacl::ocl::setup_context(1, devices);

	// Optional: Set second device for second context (use the same device for the second context if only one device available):
	
	//if (devices.size() > 1)
	//viennacl::ocl::setup_context(1, devices[1]);
	//else
	//viennacl::ocl::setup_context(1, devices[0]);
	

	std::cout << viennacl::ocl::current_device().info() << std::endl;


	//viennacl::context ctx(viennacl::ocl::get_context(0));

	viennacl::context ctx;
#else
	viennacl::context ctx;

#endif

#endif
	*/

	viennacl::context ctx;

	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
	//typedef float    ScalarType;


	/**
	* Set up the matrices and vectors for the iterative solvers (cf. iterative.cpp)
	**/
	viennacl::vcl_size_t rows = nnu;
	viennacl::vcl_size_t cols = nnu;
	viennacl::vcl_size_t nonzeros = nna;


	// Прочитать матрицу:
	std::cout << "Reading matrix..." << std::endl;
	//std::vector< std::map<unsigned int, ScalarType> > vcl_compressed_matrix1;// read_in_matrix;
	//boost::numeric::ublas::matrix<ScalarType> vcl_compressed_matrix1(rows, cols);// read_in_matrix;

	//const void* row_jumper = new integer[nnu + 1];
	//const void* col_buffer = new integer[nonzeros];
	//const ScalarType* elements = new ScalarType[nonzeros];
	integer* row_jumper = new integer[nnu + 1];
	integer* col_buffer = new integer[nna];
	ScalarType* elements = new ScalarType[nna];

	row_jumper[nnu] = nna;
	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	for (integer i = 0; i<maxelm; ++i) {
		row_jumper[i] = nna;

		// внутренность матрицы.
		if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iP; Ah.values[nna] = sl[i].ap / alpharelax;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iP) = static_cast<ScalarType> (sl[i].ap / alpharelax);
			col_buffer[nna] = sl[i].iP;
			elements[nna] = static_cast<ScalarType> (sl[i].ap / alpharelax);
			(nna)++;
		}
		if ((sl[i].iB > -1) && (fabs(sl[i].ab) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB; Ah.values[nna] = -sl[i].ab;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB) = static_cast<ScalarType> (-sl[i].ab);
			col_buffer[nna] = sl[i].iB;
			elements[nna] = static_cast<ScalarType> (-sl[i].ab);
			(nna)++;
		}
		if ((sl[i].iE > -1) && (fabs(sl[i].ae) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE; Ah.values[nna] = -sl[i].ae;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE) = static_cast<ScalarType> (-sl[i].ae);
			col_buffer[nna] = sl[i].iE;
			elements[nna] = static_cast<ScalarType> (-sl[i].ae);
			(nna)++;
		}
		if ((sl[i].iN > -1) && (fabs(sl[i].an) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN; Ah.values[nna] = -sl[i].an;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN) = static_cast<ScalarType> (-sl[i].an);
			col_buffer[nna] = sl[i].iN;
			elements[nna] = static_cast<ScalarType> (-sl[i].an);
			(nna)++;
		}
		if ((sl[i].iS > -1) && (fabs(sl[i].as) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS; Ah.values[nna] = -sl[i].as;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS) = static_cast<ScalarType> (-sl[i].as);
			col_buffer[nna] = sl[i].iS;
			elements[nna] = static_cast<ScalarType> (-sl[i].as);
			(nna)++;
		}
		if ((sl[i].iT > -1) && (fabs(sl[i].at) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT; Ah.values[nna] = -sl[i].at;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT) = static_cast<ScalarType> (-sl[i].at);
			col_buffer[nna] = sl[i].iT;
			elements[nna] = static_cast<ScalarType> (-sl[i].at);
			(nna)++;
		}
		if ((sl[i].iW > -1) && (fabs(sl[i].aw) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW; Ah.values[nna] = -sl[i].aw;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW) = static_cast<ScalarType> (-sl[i].aw);
			col_buffer[nna] = sl[i].iW;
			elements[nna] = static_cast<ScalarType> (-sl[i].aw);
			(nna)++;
		}

		if ((sl[i].iB2 > -1) && (fabs(sl[i].ab2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB2; Ah.values[nna] = -sl[i].ab2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB2) = static_cast<ScalarType> (-sl[i].ab2);
			col_buffer[nna] = sl[i].iB2;
			elements[nna] = static_cast<ScalarType> (-sl[i].ab2);
			(nna)++;
		}
		if ((sl[i].iE2 > -1) && (fabs(sl[i].ae2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE2; Ah.values[nna] = -sl[i].ae2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE2) = static_cast<ScalarType> (-sl[i].ae2);
			col_buffer[nna] = sl[i].iE2;
			elements[nna] = static_cast<ScalarType> (-sl[i].ae2);
			(nna)++;
		}
		if ((sl[i].iN2 > -1) && (fabs(sl[i].an2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN2; Ah.values[nna] = -sl[i].an2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN2) = static_cast<ScalarType> (-sl[i].an2);
			col_buffer[nna] = sl[i].iN2;
			elements[nna] = static_cast<ScalarType> (-sl[i].an2);
			(nna)++;
		}
		if ((sl[i].iS2 > -1) && (fabs(sl[i].as2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS2; Ah.values[nna] = -sl[i].as2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS2) = static_cast<ScalarType>( -sl[i].as2);
			col_buffer[nna] = sl[i].iS2;
			elements[nna] = static_cast<ScalarType> (-sl[i].as2);
			(nna)++;
		}
		if ((sl[i].iT2 > -1) && (fabs(sl[i].at2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT2; Ah.values[nna] = -sl[i].at2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT2) = static_cast<ScalarType> (-sl[i].at2);
			col_buffer[nna] = sl[i].iT2;
			elements[nna] = static_cast<ScalarType> (-sl[i].at2);
			(nna)++;
		}
		if ((sl[i].iW2 > -1) && (fabs(sl[i].aw2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW2; Ah.values[nna] = -sl[i].aw2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW2) = static_cast<ScalarType> (-sl[i].aw2);
			col_buffer[nna] = sl[i].iW2;
			elements[nna] = static_cast<ScalarType> (-sl[i].aw2);
			(nna)++;
		}

		if ((sl[i].iB3 > -1) && (fabs(sl[i].ab3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB3; Ah.values[nna] = -sl[i].ab3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB3) = static_cast<ScalarType>(-sl[i].ab3);
			col_buffer[nna] = sl[i].iB3;
			elements[nna] = static_cast<ScalarType> (-sl[i].ab3);
			(nna)++;
		}
		if ((sl[i].iE3 > -1) && (fabs(sl[i].ae3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE3; Ah.values[nna] = -sl[i].ae3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE3) = static_cast<ScalarType>(-sl[i].ae3);
			col_buffer[nna] = sl[i].iE3;
			elements[nna] = static_cast<ScalarType> (-sl[i].ae3);
			(nna)++;
		}
		if ((sl[i].iN3 > -1) && (fabs(sl[i].an3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN3; Ah.values[nna] = -sl[i].an3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN3) = static_cast<ScalarType> (-sl[i].an3);
			col_buffer[nna] = sl[i].iN3;
			elements[nna] = static_cast<ScalarType> (-sl[i].an3);
			(nna)++;
		}
		if ((sl[i].iS3 > -1) && (fabs(sl[i].as3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS3; Ah.values[nna] = -sl[i].as3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS3) = static_cast<ScalarType>(-sl[i].as3);
			col_buffer[nna] = sl[i].iS3;
			elements[nna] = static_cast<ScalarType> (-sl[i].as3);
			(nna)++;
		}
		if ((sl[i].iT3 > -1) && (fabs(sl[i].at3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT3; Ah.values[nna] = -sl[i].at3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT3) = static_cast<ScalarType>(-sl[i].at3);
			col_buffer[nna] = sl[i].iT3;
			elements[nna] = static_cast<ScalarType> (-sl[i].at3);
			(nna)++;
		}
		if ((sl[i].iW3 > -1) && (fabs(sl[i].aw3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW3; Ah.values[nna] = -sl[i].aw3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW3) = static_cast<ScalarType>(-sl[i].aw3);
			col_buffer[nna] = sl[i].iW3;
			elements[nna] = static_cast<ScalarType> (-sl[i].aw3);
			(nna)++;
		}

		if ((sl[i].iB4 > -1) && (fabs(sl[i].ab4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB4; Ah.values[nna] = -sl[i].ab4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB4) = static_cast<ScalarType>(-sl[i].ab4);
			col_buffer[nna] = sl[i].iB4;
			elements[nna] = static_cast<ScalarType> (-sl[i].ab4);
			(nna)++;
		}
		if ((sl[i].iE4 > -1) && (fabs(sl[i].ae4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE4; Ah.values[nna] = -sl[i].ae4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE4) = static_cast<ScalarType>(-sl[i].ae4);
			col_buffer[nna] = sl[i].iE4;
			elements[nna] = static_cast<ScalarType> (-sl[i].ae4);
			(nna)++;
		}
		if ((sl[i].iN4 > -1) && (fabs(sl[i].an4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN4; Ah.values[nna] = -sl[i].an4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN4) = static_cast<ScalarType>(-sl[i].an4);
			col_buffer[nna] = sl[i].iN4;
			elements[nna] = static_cast<ScalarType> (-sl[i].an4);
			(nna)++;
		}
		if ((sl[i].iS4 > -1) && (fabs(sl[i].as4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS4; Ah.values[nna] = -sl[i].as4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS4) = static_cast<ScalarType>(-sl[i].as4);
			col_buffer[nna] = sl[i].iS4;
			elements[nna] = static_cast<ScalarType> (-sl[i].as4);
			(nna)++;
		}
		if ((sl[i].iT4 > -1) && (fabs(sl[i].at4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT4; Ah.values[nna] = -sl[i].at4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT4) = static_cast<ScalarType>(-sl[i].at4);
			col_buffer[nna] = sl[i].iT4;
			elements[nna] = static_cast<ScalarType> (-sl[i].at4);
			(nna)++;
		}
		if ((sl[i].iW4 > -1) && (fabs(sl[i].aw4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW4; Ah.values[nna] = -sl[i].aw4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW4) = static_cast<ScalarType>(-sl[i].aw4);
			col_buffer[nna] = sl[i].iW4;
			elements[nna] = static_cast<ScalarType> (-sl[i].aw4);
			(nna)++;
		}


	}

	for (integer i = 0; i<maxbound; ++i) {
		row_jumper[maxelm + i] = nna;
		// граничные узлы.
		if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) {
			//Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iW; Ah.values[nna] = slb[i].aw;
			//vcl_compressed_matrix1(slb[i].iW, slb[i].iW) = static_cast<ScalarType> (slb[i].aw);
			col_buffer[nna] = slb[i].iW;
			elements[nna] = static_cast<ScalarType> (slb[i].aw);
			(nna)++;
		}
		if ((slb[i].iI > -1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			//Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iI; Ah.values[nna] = -slb[i].ai;
			//vcl_compressed_matrix1(slb[i].iW, slb[i].iI) = static_cast<ScalarType>(-slb[i].ai);
			col_buffer[nna] = slb[i].iI;
			elements[nna] = static_cast<ScalarType> (-slb[i].ai);
			(nna)++;
		}
	}

	//viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(rows, cols, nonzeros, ctx);
	viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(rows, cols, nonzeros, ctx);
	//viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(row_jumper, col_buffer, elements, rows, cols, nonzeros);
	vcl_compressed_matrix.set(row_jumper, col_buffer, elements, rows, cols, nonzeros);

	// free
	delete[] row_jumper;
	delete[] col_buffer;
	delete[] elements;

	//viennacl::copy(read_in_matrix, vcl_compressed_matrix);
	//viennacl::copy(vcl_compressed_matrix1, vcl_compressed_matrix);

	std::cout << "Reading matrix completed." << std::endl;

	// Вектора:
	//viennacl::vector<ScalarType> vcl_vec(vcl_compressed_matrix.size1(), ctx);
	//viennacl::vector<ScalarType> vcl_result(vcl_compressed_matrix.size1(), ctx);
	viennacl::vector<ScalarType> vcl_vec(rows, ctx);
	viennacl::vector<ScalarType> vcl_result(rows, ctx);

	std::vector<ScalarType> std_vec, std_result;

	// rhs and result vector:
	//std_vec.resize(vcl_compressed_matrix.size1());
	//std_result.resize(vcl_compressed_matrix.size1());
	std_vec.resize(rows);
	std_result.resize(rows);
	for (std::size_t i = 0; i < std_result.size(); ++i) {
		//std_result[i] = ScalarType(1);
		std_result[i] = dX0[i];
		std_vec[i] = dV[i];
	}

	// Copy to GPU
	viennacl::copy(std_vec, vcl_vec);
	viennacl::copy(std_result, vcl_result);

	//vcl_vec = viennacl::linalg::prod(vcl_compressed_matrix, vcl_result);





	/**
	* Instantiate a tag for the conjugate gradient solver, the AMG preconditioner tag, and create an AMG preconditioner object:
	**/
	

	

	/**
	* Run solver without preconditioner. This serves as a baseline for comparison.
	* Note that iterative solvers without preconditioner on GPUs can be very efficient because they map well to the massively parallel hardware.
	**/
	//std::cout << "-- BiCGStab solver (no preconditioner, warmup) --" << std::endl;
	//run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, bicgstab_solver, viennacl::linalg::no_precond());

	if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 5) {
		/**
		* Generate the setup for an AMG preconditioner of Ruge-Stueben type with only one pass and direct interpolation (ONEPASS+DIRECT)
		**/
		viennacl::context host_ctx(viennacl::MAIN_MEMORY);
		//viennacl::context target_ctx = viennacl::traits::context(vcl_compressed_matrix);
		viennacl::context target_ctx(viennacl::CUDA_MEMORY);

		viennacl::linalg::bicgstab_tag my_bicgstab_solver(1e-6, 10000); // 1e-8

		std::cout << "Start calculation" << std::endl;
		viennacl::linalg::amg_tag amg_tag_direct;
		amg_tag_direct.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_ONEPASS);
		amg_tag_direct.set_interpolation_method(viennacl::linalg::AMG_INTERPOLATION_METHOD_DIRECT);
		amg_tag_direct.set_strong_connection_threshold(0.25);
		amg_tag_direct.set_jacobi_weight(0.67);
		amg_tag_direct.set_presmooth_steps(1);
		amg_tag_direct.set_postsmooth_steps(1);
		amg_tag_direct.set_setup_context(host_ctx);    // run setup on host
		amg_tag_direct.set_target_context(target_ctx); // run solver cycles on device
		run_amg(my_bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "ONEPASS COARSENING, DIRECT INTERPOLATION", amg_tag_direct);

		std::cout << "No. of iters: " << my_bicgstab_solver.iters() << std::endl;
		std::cout << "Est. error: " << my_bicgstab_solver.error() << std::endl;
	}
	else if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 9) {

		std::cout << "Start calculation" << std::endl;

		viennacl::linalg::bicgstab_tag my_bicgstab_solver(1e-6, 2000); // 1e-8

		//ILU0
		viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<ScalarType> > vcl_ilu0(vcl_compressed_matrix, viennacl::linalg::ilu0_tag());
		// using viennacl objects on GPU
		//vcl_result = solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::bicgstab_tag(1e-6, 2000), vcl_ilu0);

		std::cout << "ILU0 decomposition found" << std::endl;

		vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec, my_bicgstab_solver, vcl_ilu0); //with preconditioner

																																	 // print number of iterations taken and estimated error:
		std::cout << "No. of iters: " << my_bicgstab_solver.iters() << std::endl;
		std::cout << "Est. error: " << my_bicgstab_solver.error() << std::endl;

	}
	else if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 10) {

		//viennacl::linalg::ilut_tag ilut_conf(10, 1e-5); // 10 entries, rel. tol. 1e-5
		//viennacl::linalg::ilut_tag ilut_conf(40, 1e-2); // 40 entries, rel. tol. 1e-2
		//typedef viennacl::linalg::ilut_precond<viennacl::compressed_matrix<ScalarType> > vcl_ilut_t;

		//vcl_ilut_t vcl_ilut(vcl_compressed_matrix, ilut_conf);

		//---->viennacl::linalg::ilut_precond< viennacl::compressed_matrix<ScalarType> > vcl_ilut(vcl_compressed_matrix, viennacl::linalg::ilut_tag());

		// using viennacl objects on GPU
		//---->vcl_result = solve(vcl_compressed_matrix, vcl_vec, viennacl::linalg::bicgstab_tag(1e-6, 2000), vcl_ilut);

		// compute ILU0 preconditioner:
		viennacl::linalg::ilu0_tag ilu0_config;
		viennacl::linalg::ilu0_precond< viennacl::compressed_matrix<ScalarType>  > vcl_ilut(vcl_compressed_matrix, ilu0_config);
		// solve (e.g. using conjugate gradient solver)
		vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_vec,
			viennacl::linalg::bicgstab_tag(1e-6, 2000), // solver here
			vcl_ilut);                        // preconditioner here

	}

	/**
	* Generate the setup for an aggregation-based AMG preconditioner with unsmoothed aggregation
	**/
	/*
	viennacl::linalg::amg_tag amg_tag_agg_pmis;
	amg_tag_agg_pmis.set_coarsening_method(viennacl::linalg::AMG_COARSENING_METHOD_MIS2_AGGREGATION);
	amg_tag_agg_pmis.set_interpolation_method(viennacl::linalg::AMG_INTERPOLATION_METHOD_AGGREGATION);
	run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), AG INTERPOLATION", amg_tag_agg_pmis);
	*/

	/**
	* Generate the setup for a smoothed aggregation AMG preconditioner
	**/
	/*
	viennacl::linalg::amg_tag amg_tag_sa_pmis;
	amg_tag_sa_pmis.set_coarsening_method(viennacl::cuda::linalg::AMG_COARSENING_METHOD_MIS2_AGGREGATION);
	amg_tag_sa_pmis.set_interpolation_method(viennacl::cuda::linalg::AMG_INTERPOLATION_METHOD_SMOOTHED_AGGREGATION);
	run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), SA INTERPOLATION", amg_tag_sa_pmis);
	*/

	//std::cout << std::endl;
	//std::cout << " -------------- Benchmark runs -------------- " << std::endl;
	//std::cout << std::endl;

	//std::cout << "-- BiCGStab solver (no preconditioner) --" << std::endl;
	//run_solver(vcl_compressed_matrix, vcl_vec, vcl_result, cg_solver, viennacl::linalg::no_precond());
	//run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "ONEPASS COARSENING, DIRECT INTERPOLATION", amg_tag_direct);
	//run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), AG INTERPOLATION", amg_tag_agg_pmis);
	//run_amg(bicgstab_solver, vcl_vec, vcl_result, vcl_compressed_matrix, "AG COARSENING (PMIS), SA INTERPOLATION", amg_tag_sa_pmis);

	/**
	*  That's it.
	**/
	
	std::cout << "!!!! CALCULATION COMPLETED SUCCESSFULLY !!!!" << std::endl;

	// Обратное копирование:
	viennacl::copy(vcl_result, std_result);

	for (std::size_t i = 0; i < rows; ++i) {
		// Возвращение результата вычисления.
		dX0[i] = std_result[i];
		//dX0[i] = vcl_result(i);
	}
}


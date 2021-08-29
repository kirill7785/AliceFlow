// my_amgcl_alg_builtin.cpp
// Алгебраический многосеточный метод Дениса Демидова.
// дата подключения библиотеки amgcl в проект AliceFlow 07-05-2019, 08-05-2019.
// Для присоединения библиотеки amgcl к проекту в папке
// C:\Program Files (x86)\Windows Kits\10\Include\10.0.15063.0\um
// Для VS2017 должны лежать исходные коды библиотек amgcl (папка amgcl)
// и библиотека boost (папка boost). Кроме того в папке с кодом данного проекта должен лежать
// файл amgcl.h из папки lib. А код файла amgcl.cpp должен быть добавлен в код проекта AliceFlow
// перед функцией передачи матрицы и вызова решателя. 
// Данная связка собралась не из под cuda а как обычное консольное 
// приложение на С++.
// 
// Достоинство свежесть и постоянная обновляемость библиотеки, 
// поддержка современными компиляторами.

#pragma once
#ifndef _DENIS_DEMIDOV_MY_AMGCL_ALG_CPU_ 
#define _DENIS_DEMIDOV_MY_AMGCL_ALG_CPU_ 1


#include <iostream>
#include <vector>


//#include "lib/amgcl.cpp"
//#include "lib/amgcl.h"
//#include "sample_problem.hpp"



#include <type_traits>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>



#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/solver/bicgstab.hpp>

#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/block_matrix.hpp>

#include <amgcl/io/mm.hpp>
#include <amgcl/profiler.hpp>

#include "amgcl.h"


#ifdef AMGCL_PROFILING
#include <amgcl/profiler.hpp>
namespace amgcl {
	profiler<> prof;
}
#endif

//---------------------------------------------------------------------------
typedef amgcl::backend::builtin<double>           Backend;
typedef amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper> AMG;
typedef amgcl::runtime::solver::wrapper<Backend>  ISolver;
typedef amgcl::make_solver<AMG, ISolver>          Solver;
typedef boost::property_tree::ptree               Params;

//---------------------------------------------------------------------------
amgclHandle STDCALL amgcl_params_create() {
	return static_cast<amgclHandle>(new Params());
}

//---------------------------------------------------------------------------
void STDCALL amgcl_params_seti(amgclHandle prm, const char *name, int value) {
	static_cast<Params*>(prm)->put(name, value);
}

//---------------------------------------------------------------------------
void STDCALL amgcl_params_setf(amgclHandle prm, const char *name, float value) {
	static_cast<Params*>(prm)->put(name, value);
}

//---------------------------------------------------------------------------
void STDCALL amgcl_params_sets(amgclHandle prm, const char *name, const char *value) {
	static_cast<Params*>(prm)->put(name, value);
}

//---------------------------------------------------------------------------
void STDCALL amgcl_params_read_json(amgclHandle prm, const char *fname) {
	read_json(fname, *static_cast<Params*>(prm));
}

//---------------------------------------------------------------------------
void STDCALL amgcl_params_destroy(amgclHandle prm) {
	delete static_cast<Params*>(prm);
}

//---------------------------------------------------------------------------
amgclHandle STDCALL amgcl_precond_create(
	int           n,
	const int    *ptr,
	const int    *col,
	const double *val,
	amgclHandle   prm
)
{
	auto A = std::make_tuple(n,
		boost::make_iterator_range(ptr, ptr + n + 1),
		boost::make_iterator_range(col, col + ptr[n]),
		boost::make_iterator_range(val, val + ptr[n])
	);

	if (prm)
		return static_cast<amgclHandle>(new AMG(A, *static_cast<Params*>(prm)));
	else
		return static_cast<amgclHandle>(new AMG(A));
}

//---------------------------------------------------------------------------
amgclHandle STDCALL amgcl_precond_create_f(
	int           n,
	const int    *ptr,
	const int    *col,
	const double *val,
	amgclHandle   prm
)
{
	auto ptr_c = boost::make_transform_iterator(ptr, [](int i) { return i - 1; });
	auto col_c = boost::make_transform_iterator(col, [](int i) { return i - 1; });

	auto A = std::make_tuple(n,
		boost::make_iterator_range(ptr_c, ptr_c + n + 1),
		boost::make_iterator_range(col_c, col_c + ptr[n]),
		boost::make_iterator_range(val, val + ptr[n])
	);

	if (prm)
		return static_cast<amgclHandle>(new AMG(A, *static_cast<Params*>(prm)));
	else
		return static_cast<amgclHandle>(new AMG(A));
}

//---------------------------------------------------------------------------
void STDCALL amgcl_precond_apply(amgclHandle handle, const double *rhs, double *x)
{
	AMG *amg = static_cast<AMG*>(handle);

	size_t n = amgcl::backend::rows(amg->system_matrix());

	boost::iterator_range<double*> x_range =
		boost::make_iterator_range(x, x + n);

	amg->apply(boost::make_iterator_range(rhs, rhs + n), x_range);
}

//---------------------------------------------------------------------------
void STDCALL amgcl_precond_report(amgclHandle handle) {
	std::cout << *static_cast<AMG*>(handle) << std::endl;
}

//---------------------------------------------------------------------------
void STDCALL amgcl_precond_destroy(amgclHandle handle) {
	delete static_cast<AMG*>(handle);
}

//---------------------------------------------------------------------------
amgclHandle STDCALL amgcl_solver_create(
	int           n,
	const int    *ptr,
	const int    *col,
	const double *val,
	amgclHandle   prm
)
{
	auto A = std::make_tuple(n,
		boost::make_iterator_range(ptr, ptr + n + 1),
		boost::make_iterator_range(col, col + ptr[n]),
		boost::make_iterator_range(val, val + ptr[n])
	);

	if (prm)
		return static_cast<amgclHandle>(new Solver(A, *static_cast<Params*>(prm)));
	else
		return static_cast<amgclHandle>(new Solver(A));
}

//---------------------------------------------------------------------------
amgclHandle STDCALL amgcl_solver_create_f(
	int           n,
	const int    *ptr,
	const int    *col,
	const double *val,
	amgclHandle   prm
)
{
	auto ptr_c = boost::make_transform_iterator(ptr, [](int i) { return i - 1; });
	auto col_c = boost::make_transform_iterator(col, [](int i) { return i - 1; });

	auto A = std::make_tuple(n,
		boost::make_iterator_range(ptr_c, ptr_c + n + 1),
		boost::make_iterator_range(col_c, col_c + ptr[n]),
		boost::make_iterator_range(val, val + ptr[n])
	);

	if (prm)
		return static_cast<amgclHandle>(new Solver(A, *static_cast<Params*>(prm)));
	else
		return static_cast<amgclHandle>(new Solver(A));
}

//---------------------------------------------------------------------------
void STDCALL amgcl_solver_report(amgclHandle handle) {
	std::cout << static_cast<Solver*>(handle)->precond() << std::endl;
}

//---------------------------------------------------------------------------
void STDCALL amgcl_solver_destroy(amgclHandle handle) {
	delete static_cast<Solver*>(handle);
}

//---------------------------------------------------------------------------
conv_info STDCALL amgcl_solver_solve(
	amgclHandle handle,
	const double *rhs,
	double *x
)
{
	Solver *slv = static_cast<Solver*>(handle);

	size_t n = slv->size();

	conv_info cnv;

	boost::iterator_range<double*> x_range = boost::make_iterator_range(x, x + n);

	std::tie(cnv.iterations, cnv.residual) = (*slv)(
		boost::make_iterator_range(rhs, rhs + n), x_range
		);

	return cnv;
}

//---------------------------------------------------------------------------
conv_info STDCALL amgcl_solver_solve_mtx(
	amgclHandle handle,
	int    const * A_ptr,
	int    const * A_col,
	double const * A_val,
	const double *rhs,
	double *x
)
{
	Solver *slv = static_cast<Solver*>(handle);

	size_t n = slv->size();

	conv_info cnv;

	boost::iterator_range<double*> x_range = boost::make_iterator_range(x, x + n);

	std::tie(cnv.iterations, cnv.residual) = (*slv)(
		std::make_tuple(
			n,
			boost::make_iterator_range(A_ptr, A_ptr + n),
			boost::make_iterator_range(A_col, A_col + A_ptr[n]),
			boost::make_iterator_range(A_val, A_val + A_ptr[n])
		),
		boost::make_iterator_range(rhs, rhs + n), x_range
		);

	return cnv;
}


void amgcl_networkT_solver_cpu(doublereal* &val, integer* &col_ind, integer* &row_ptr,
	integer n,	doublereal *dV, doublereal* &dX0,
	bool bprintmessage)
{
	// размерность квадратной матрицы.
	// n

#ifdef _OPENMP 
	// Узнаёт количество ядер в системе.
	// 15млн неизвестных время параллельного кода 21мин 39с.
	// Время однопоточного кода 27мин 48с.
	unsigned int nthreads = number_cores();
	omp_set_num_threads(nthreads); // установка числа потоков
#endif

	 // Разреженная матрица СЛАУ
	 // в CRS формате.
	 // val; col_ind; row_ptr;


	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
									   //typedef float    ScalarType;
	typedef int Myindextype;

	integer nnu = n;
	integer nna = row_ptr[n];
	if (bprintmessage) {
		std::cout << "nna=" << nna << " " << row_ptr[n - 1] << "\n";
	}
	//system("pause");

	/**
	* Set up the matrices and vectors for the iterative solvers (cf. iterative.cpp)
	**/
	Myindextype rows = static_cast<Myindextype>(nnu);
	Myindextype cols = static_cast<Myindextype>(nnu);
	Myindextype nonzeros = static_cast<Myindextype>(nna);

	// Прочитать матрицу:
	if (bprintmessage) {
		std::cout << "Reading matrix..." << std::endl;
	}

	//integer* row_jumper = new integer[nnu + 1];
	//integer* col_buffer = new integer[nna];
	//ScalarType* elements = new ScalarType[nna];

	std::vector<int> row_jumper(nnu + 1); // только тип int !!!
	std::vector<int> col_buffer(nna); // только тип int !!!
	std::vector<double> elements(nna);
	std::vector<double> rhs(nnu);

	for (integer i = 0; i < nnu; i++) {
		rhs[i] = dV[i];
	}

	row_jumper[nnu] = static_cast<Myindextype>(nna);
	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	integer iscan = 0;
	for (integer i = 0; i < n; i++) {

		row_jumper[i] = static_cast<Myindextype>(row_ptr[i]);


		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
			if (col_ind[j] == i) {
				elements[iscan] = static_cast<ScalarType>(val[j]);
				col_buffer[iscan] = static_cast<Myindextype>(col_ind[j]);
				iscan++;
				//debug
				//std::cout << "col_buffer[j]=" << col_ind[j] << std::endl;
				//system("pause");
			}
		}
		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
			if (col_ind[j] != i) {
				elements[iscan] = static_cast<ScalarType>(val[j]);
				col_buffer[iscan] = static_cast<Myindextype>(col_ind[j]);
				iscan++;
				// debug
				//std::cout << "col_buffer[j]=" << col_ind[j] << std::endl;
				//system("pause");
			}
		}

	}

	row_jumper[n] = static_cast<Myindextype>(row_ptr[n]);

	

	if (bprintmessage) {
		printf("Matrix load succsefull...\n");
	}

	amgclHandle prm = amgcl_params_create();

	// число ячеек на самом грубом уровне.
	amgcl_params_seti(prm, "precond.coarse_enough", 40);// Для графовой модели.

	switch (my_amg_manager.amgcl_selector) {
	case 0: // Ruge - Stueben (аналог amg1r5)
		amgcl_params_sets(prm, "precond.coarsening.type", "ruge_stuben");
		amgcl_params_setf(prm, "precond.coarsening.eps_strong", 0.9f);
		//{type = ruge_stuben, eps_strong = 0.9}
		break;
	case 1: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		break;
	default: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		break;
	}

	switch (my_amg_manager.amgcl_smoother) {
	case 0: // spai0
		amgcl_params_sets(prm, "precond.relax.type", "spai0");
		break;
	case 1: // ilu0
		amgcl_params_sets(prm, "precond.relax.type", "ilu0");
		break;
	case 2: // gauss_seidel
		amgcl_params_sets(prm, "precond.relax.type", "gauss_seidel");
		break;
	case 3: // damped_jacobi
		amgcl_params_sets(prm, "precond.relax.type", "damped_jacobi");
		amgcl_params_setf(prm, "precond.relax.damping", 0.8f);
		break;
	case 4: // spai1
		amgcl_params_sets(prm, "precond.relax.type", "spai1");
		break;
	case 5: // chebyshev
		amgcl_params_sets(prm, "precond.relax.type", "chebyshev");
		break;
	case 6: // ilu1
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 1);
		break;
	case 7: // ilu2
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 2);
		break;
	case 8: // ilu4
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 4);
		printf("precond.relax.type==ilu(k==4).\n");
		break;
	case 9: // ilu6
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 6);
		printf("precond.relax.type==ilu(k==6).\n");
		break;
	case 10: // ilu8
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 8);
		printf("precond.relax.type==ilu(k==8).\n");
		break;
	case 11: // ilu10
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 10);
		printf("precond.relax.type==ilu(k==10).\n");
		break;
	default:	amgcl_params_sets(prm, "precond.relax.type", "spai0");
		break;
	}
	//amgcl_params_sets(prm, "solver.type", "bicgstabl");
	//amgcl_params_seti(prm, "solver.L", 1);
	switch (my_amg_manager.amgcl_iterator) {
	case AMGCL_ITERATOR_ALG::BiCGStab: // BiCGStab
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		break;
	case AMGCL_ITERATOR_ALG::FGMRes: // FGMRes
		amgcl_params_sets(prm, "solver.type", "fgmres");
		break;
	default:
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		break;
	}

	if (bglobal_unsteady_temperature_determinant) {
		// Нестационарные задачи требуется считать до меньших значений невязки.
		// 14.05.2019
		// 1.0e-12 точность достаточна.
		if (my_amg_manager.amgcl_iterator == AMGCL_ITERATOR_ALG::FGMRes) {
			amgcl_params_setf(prm, "solver.tol", 1.0e-12f);// пробовал 1e-7 для fgmres полотно АФАР расходится.
		}
		else {
			amgcl_params_setf(prm, "solver.tol", 1.0e-12f);
		}
		// Нестационарные задачи требуется считать до меньших значений невязки.
		// Может локально не сойтись поэтому не надо делать очень большого числа итераций.
		// 14.05.2019
		amgcl_params_seti(prm, "solver.maxiter", 3000);
	}


	

	if (bprintmessage) {
		printf("Setup phase start...\n");
	}

	bool bprint_preconditioner = bprintmessage;

	//****
	if (bprint_preconditioner) {
		// Библиотека Дениса Демидова работает только с 32 битным 
		// типом int. Его хватает даже для числа неизвестных 80 млн
		// контрольных объемов.
		int n_loc = (int)(n);
		amgclHandle amg_precond = amgcl_precond_create(
			n_loc, row_jumper.data(), col_buffer.data(), elements.data(), prm
		);
		amgcl_precond_report(amg_precond);//печать samg предобуславливателя.
		amgcl_precond_destroy(amg_precond);
	}
	//*****

	// Библиотека Дениса Демидова работает только с 32 битным 
	// типом int. Его хватает даже для числа неизвестных 80 млн
	// контрольных объемов.
	int n_loc = (int)(n);
	amgclHandle solver = amgcl_solver_create(
		n_loc, row_jumper.data(), col_buffer.data(), elements.data(), prm
	);



	amgcl_params_destroy(prm);

	if (bprintmessage) {
		printf("Solution phase start...\n");
	}

	std::vector<double> x(n, 0);

	if (1) {
		// Инициализация
		for (integer i = 0; i < n; i++) {
			x[i] = dX0[i];
		}
	}

	conv_info cnv = amgcl_solver_solve(solver, rhs.data(), x.data());

	// Solve same problem again, but explicitly provide the matrix this time:
	//std::fill(x.begin(), x.end(), 0);
	//cnv = amgcl_solver_solve_mtx(
	//solver, row_jumper.data(), col_buffer.data(), elements.data(),
	//rhs.data(), x.data()
	//);

	if (bprintmessage) {
		std::cout << "Iterations: " << cnv.iterations << std::endl
			<< "Error:      " << cnv.residual << std::endl;
	}

	amgcl_solver_destroy(solver);


	for (integer i = 0; i < nnu; i++) {
		dX0[i] = x[i];
	}

#ifdef _OPENMP 
	omp_set_num_threads(1); // установка числа потоков
#endif

	if (bprintmessage) {
		getchar();
	}

} //amgcl_network_T_solver

// Так и не заработал сколько не пытался. 18,07,2021
void amgcl_secondT_solverStructural_block3x3_cpu(SIMPLESPARSE& sparseM, integer n,
	doublereal* dV, doublereal*& dX0,
	bool bprintmessage, WALL*& w, int& lw, bool*& bondary,
	bool bStructural_Mechanic)
{
	// размерность квадратной матрицы.
	// n

#ifdef _OPENMP 
	// Узнаёт количество ядер в системе.
	// 15млн неизвестных время параллельного кода 21мин 39с.
	// Время однопоточного кода 27мин 48с.
	unsigned int nthreads = number_cores();
	omp_set_num_threads(nthreads); // установка числа потоков
#endif

								   // Разреженная матрица СЛАУ
								   // в CRS формате.
	doublereal* val = nullptr;
	integer_mix_precision* col_ind = nullptr;
	integer_mix_precision* row_ptr = nullptr;




	simplesparsetoCRS(sparseM, val, col_ind, row_ptr, n); // преобразование матрицы из одного формата хранения в другой.
																//m.ballocCRScfd = true;
	simplesparsefree(sparseM, n);


	typedef double    ScalarType;  // feel free to change this to double if supported by your device
									   //typedef float    ScalarType;
	typedef int Myindextype;

	integer nnu = n;
	integer nna = row_ptr[n];
	printf("nna=%lld %d\n", nna, row_ptr[n - 1]);
	//system("pause");

	/**
	* Set up the matrices and vectors for the iterative solvers (cf. iterative.cpp)
	**/
	Myindextype rows = static_cast<Myindextype>(nnu);
	Myindextype cols = static_cast<Myindextype>(nnu);
	Myindextype nonzeros = static_cast<Myindextype>(nna);

	// Прочитать матрицу:
	std::cout << "Reading matrix..." << std::endl;

	//integer* row_jumper = new integer[nnu + 1];
	//integer* col_buffer = new integer[nna];
	//ScalarType* elements = new ScalarType[nna];

	std::vector<ptrdiff_t> row_jumper(nnu + 1); // только тип int !!!
	std::vector<ptrdiff_t> col_buffer(nna); // только тип int !!!
	
	std::vector<double> elements(nna);
	std::vector<double> rhs(nnu);

	for (integer i = 0; i < nnu; i++) {
		rhs[i] = dV[i];
	}

	std::vector<double> x(n, 0);

	if (1) {
		// Инициализация
		for (integer i = 0; i < n; i++) {
			x[i] = dX0[i];
		}
	}

	doublereal alpharelax = 1.0;

	if (!bStructural_Mechanic) {
		// Уравнение теплопераедачи.

		// Это не специальная нелинейная версия кода amg1r5 CAMG.
		for (int k = 0; k < lw; k++) {
			if ((w[k].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
				(w[k].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY)) {
				alpharelax = 0.99999; // Для того чтобы СЛАУ сходилась.
				// 0.9999 - недостаточное значение, температуры не те получаются.
			}
		}
		if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {
			alpharelax = 0.99999;
		}
		//if (adiabatic_vs_heat_transfer_coeff == ADIABATIC_WALL_BC) {
		//printf("ADIABATIC WALL BC"); getchar();
		//}

	}

	row_jumper[nnu] = static_cast<Myindextype>(nna);
	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	integer iscan = 0;
	for (integer i = 0; i < n; i++) {

		row_jumper[i] = static_cast<Myindextype>(row_ptr[i]);


		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
			if (col_ind[j] == i) {
				if ((!bStructural_Mechanic) && ((bondary != nullptr) && (!bondary[i])) && (row_ptr[i + 1] > row_ptr[i] + 1)) {
					// Уравнение теплопередачи.
					// Релаксация к предыдущему значению.
					elements[iscan] = static_cast<ScalarType>(val[j] / alpharelax);
					rhs[i] += (1.0 - alpharelax) * val[j] * x[i] / alpharelax;
				}
				else {
					elements[iscan] = static_cast<ScalarType>(val[j]); // Условие Дирихле.
				}
				col_buffer[iscan] = static_cast<Myindextype>(col_ind[j]);
				iscan++;
				//debug
				//std::cout << "col_buffer[j]=" << col_ind[j] << std::endl;
				//system("pause");
			}
		}
		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
			if (col_ind[j] != i) {
				elements[iscan] = static_cast<ScalarType>(val[j]);
				col_buffer[iscan] = static_cast<Myindextype>(col_ind[j]);
				iscan++;
				// debug
				//std::cout << "col_buffer[j]=" << col_ind[j] << std::endl;
				//system("pause");
			}
		}

	}

	row_jumper[n] = static_cast<Myindextype>(row_ptr[n]);

	delete[] val;
	val = nullptr;

	delete[] col_ind;
	col_ind = nullptr;

	delete[] row_ptr;
	row_ptr = nullptr;



	printf("Matrix load succsefull...\n");

	/*
	amgclHandle prm = amgcl_params_create();

	// число ячеек на самом грубом уровне.
	amgcl_params_seti(prm, "precond.coarse_enough", 1000);

	switch (my_amg_manager.amgcl_selector) {
	case 0: // Ruge - Stueben (аналог amg1r5)
		amgcl_params_sets(prm, "precond.coarsening.type", "ruge_stuben");
		amgcl_params_setf(prm, "precond.coarsening.eps_strong", 0.9f);
		//{type = ruge_stuben, eps_strong = 0.9}
		break;
	case 1: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		break;
	default: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		break;
	}

	switch (my_amg_manager.amgcl_smoother) {
	case 0: // spai0
		amgcl_params_sets(prm, "precond.relax.type", "spai0");
		break;
	case 1: // ilu0
		amgcl_params_sets(prm, "precond.relax.type", "ilu0");
		break;
	case 2: // gauss_seidel
		amgcl_params_sets(prm, "precond.relax.type", "gauss_seidel");
		break;
	case 3: // damped_jacobi
		amgcl_params_sets(prm, "precond.relax.type", "damped_jacobi");
		amgcl_params_setf(prm, "precond.relax.damping", 0.8f);
		break;
	case 4: // spai1
		amgcl_params_sets(prm, "precond.relax.type", "spai1");
		break;
	case 5: // chebyshev
		amgcl_params_sets(prm, "precond.relax.type", "chebyshev");
		break;
	case 6: // ilu1
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 1);
		break;
	case 7: // ilu2
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 2);
		break;
	case 8: // ilu4
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 4);
		printf("precond.relax.type==ilu(k==4).\n");
		break;
	case 9: // ilu6
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 6);
		printf("precond.relax.type==ilu(k==6).\n");
		break;
	case 10: // ilu8
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 8);
		printf("precond.relax.type==ilu(k==8).\n");
		break;
	case 11: // ilu10
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 10);
		printf("precond.relax.type==ilu(k==10).\n");
		break;
	default:	amgcl_params_sets(prm, "precond.relax.type", "spai0");
		break;
	}
	//amgcl_params_sets(prm, "solver.type", "bicgstabl");
	//amgcl_params_seti(prm, "solver.L", 1);
	switch (my_amg_manager.amgcl_iterator) {
	case AMGCL_ITERATOR_ALG::BiCGStab: // BiCGStab
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		break;
	case AMGCL_ITERATOR_ALG::FGMRes: // FGMRes
		amgcl_params_sets(prm, "solver.type", "fgmres");
		break;
	default:
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		break;
	}

	if (bglobal_unsteady_temperature_determinant) {
		// Нестационарные задачи требуется считать до меньших значений невязки.
		// 14.05.2019
		// 1.0e-12 точность достаточна.
		amgcl_params_setf(prm, "solver.tol", 1.0e-12f);
		// Нестационарные задачи требуется считать до меньших значений невязки.
		// Может локально не сойтись поэтому не надо делать очень большого числа итераций.
		// 14.05.2019
		amgcl_params_seti(prm, "solver.maxiter", 3000);
	}

	if (bStructural_Mechanic) {
		// Механика.
		// 1.0e-12 точность достаточна.
		amgcl_params_setf(prm, "solver.tol", 1.0e-14f);
		amgcl_params_seti(prm, "solver.maxiter", 3000);
	}

	printf("Setup phase start...\n");

	bool bprint_preconditioner = true;

	//****
	if (bprint_preconditioner) {
		// Библиотека Дениса Демидова работает только с 32 битным 
		// типом int. Его хватает даже для числа неизвестных 80 млн
		// контрольных объемов.
		int n_loc = (int)(n);
		amgclHandle amg_precond = amgcl_precond_create(
			n_loc, row_jumper.data(), col_buffer.data(), elements.data(), prm
		);
		amgcl_precond_report(amg_precond);//печать samg предобуславливателя.
		amgcl_precond_destroy(amg_precond);
	}
	//*****

	// Библиотека Дениса Демидова работает только с 32 битным 
	// типом int. Его хватает даже для числа неизвестных 80 млн
	// контрольных объемов.
	int n_loc = (int)(n);
	amgclHandle solver = amgcl_solver_create(
		n_loc, row_jumper.data(), col_buffer.data(), elements.data(), prm
	);



	amgcl_params_destroy(prm);

	printf("Solution phase start...\n");




	conv_info cnv = amgcl_solver_solve(solver, rhs.data(), x.data());

	// Solve same problem again, but explicitly provide the matrix this time:
	//std::fill(x.begin(), x.end(), 0);
	//cnv = amgcl_solver_solve_mtx(
	//solver, row_jumper.data(), col_buffer.data(), elements.data(),
	//rhs.data(), x.data()
	//);

	std::cout << "Iterations: " << cnv.iterations << std::endl
		<< "Error:      " << cnv.residual << std::endl;

	amgcl_solver_destroy(solver);

	*/

int n_loc = (int)(n);


/*
// Начинается с нуля всё верно. Проверено. 17,07,2021
for (ptrdiff_t i = 0; i < n_loc; ++i) {
	std::cout << row_jumper[i] << "  " << col_buffer[i] << std::endl;
	getchar();
}*/

// Scale the matrix so that it has the unit diagonal.
   // First, find the diagonal values:
std::vector<double> D(n_loc, 1.0);
for (ptrdiff_t i = 0; i < n_loc; ++i) {
	for (ptrdiff_t j = row_jumper[i], e = row_jumper[i + 1]; j < e; ++j) {
		if (col_buffer[j] == i) {
			if (elements[j] <= 0.0) {

				std::cout << elements[j] << " " << i << " " << n_loc << " rhs=" << rhs[j] << std::endl;
				getchar();
			}
			D[i] = 1 / sqrt(fabs(elements[j]));
			break;
		}
	}
}

// Then, apply the scaling in-place:
for (ptrdiff_t i = 0; i < n_loc; ++i) {
	for (ptrdiff_t j = row_jumper[i], e = row_jumper[i + 1]; j < e; ++j) {
		elements[j] *= D[i] * D[col_buffer[j]];
	}
	rhs[i] *= D[i];
}


	// We use the tuple of CRS arrays to represent the system matrix.
	   // Note that std::tie creates a tuple of references, so no data is actually
	   // copied here:
auto A = std::tie(n_loc, row_jumper, col_buffer, elements);

// Compose the solver type
typedef amgcl::static_matrix<double, 3, 3> dmat_type; // matrix value type in double precision
typedef amgcl::static_matrix<double, 3, 1> dvec_type; // the corresponding vector value type
typedef amgcl::static_matrix<float, 3, 3> smat_type; // matrix value type in single precision

typedef amgcl::backend::builtin<dmat_type> SBackend; // the solver backend
typedef amgcl::backend::builtin<smat_type> PBackend; // the preconditioner backend

typedef amgcl::make_solver<
	amgcl::amg<
	SBackend,
	amgcl::coarsening::smoothed_aggregation,
	amgcl::relaxation::spai0
	>,
	amgcl::solver::cg<SBackend>
> Solver;

// Solver parameters
Solver::params prm;
prm.solver.maxiter = 3000;

//amgcl_params_sets(prm, "precond.relax.type", "iluk");
//amgcl_params_seti(prm, "precond.relax.k", 4);

// The profiler:
amgcl::profiler<> prof("Serena");

// Initialize the solver with the system matrix.
// Use the block_matrix adapter to convert the matrix into
// the block format on the fly:
prof.tic("setup");
auto Ab = amgcl::adapter::block_matrix<dmat_type>(A);
Solver solve54(Ab, prm);
prof.toc("setup");

// Show the mini-report on the constructed solver:
std::cout << solve54 << std::endl;

// Solve the system with the zero initial approximation:
int iters;
double error;
//std::vector<double> x(rows, 0.0);

// Reinterpret both the RHS and the solution vectors as block-valued:
auto f_ptr = reinterpret_cast<dvec_type*>(rhs.data());
auto x_ptr = reinterpret_cast<dvec_type*>(x.data());
auto F = amgcl::make_iterator_range(f_ptr, f_ptr + rows / 3);
auto X = amgcl::make_iterator_range(x_ptr, x_ptr + rows / 3);
//auto F = amgcl::make_iterator_range(f_ptr, f_ptr + rows);
//auto X = amgcl::make_iterator_range(x_ptr, x_ptr + rows);

prof.tic("solve");
std::tie(iters, error) = solve54(Ab, F, X);
prof.toc("solve");

// Output the number of iterations, the relative error,
// and the profiling data:
std::cout << "Iters: " << iters << std::endl
<< "Error: " << error << std::endl
<< prof << std::endl;


	for (integer i = 0; i < nnu; i++) {
		dX0[i] = x[i]* D[i];
	}

#ifdef _OPENMP 
	omp_set_num_threads(1); // установка числа потоков
#endif


}


// Солвер настроен на оптимально быструю работу с механикой. 
// Настройки из интерфейса менять нельзя. Настроен для чистой механики без теплопередачи.
// Рецепт: bicgstab, ruge_stuben, gauss_seidel. 18.07.2021
// Рецепт рабочий и на задаче Фламана показывает рекордное время 26с и число итераций 72 шт.
void amgcl_secondT_solverStructural_cpu(SIMPLESPARSE& sparseM, integer n,
	doublereal* dV, doublereal*& dX0,
	bool bprintmessage, WALL*& w, int& lw, bool*& bondary,
	bool bStructural_Mechanic)
{
	// размерность квадратной матрицы.
	// n

#ifdef _OPENMP 
	// Узнаёт количество ядер в системе.
	// 15млн неизвестных время параллельного кода 21мин 39с.
	// Время однопоточного кода 27мин 48с.
	unsigned int nthreads = number_cores();
	omp_set_num_threads(nthreads); // установка числа потоков
#endif

								   // Разреженная матрица СЛАУ
								   // в CRS формате.
	doublereal* val = nullptr;
	integer_mix_precision* col_ind = nullptr;
	integer_mix_precision* row_ptr = nullptr;

#if !bStableVersion

	ell_to_CRS(val, col_ind, row_ptr, n); // преобразование матрицы из одного формата хранения в другой.


	for (int i = 0; i < n; ++i) {
		delete[] data_ell[i];
		delete[] coll_ell[i];
	}
	delete[] data_ell;
	delete[] coll_ell;

#else

	simplesparsetoCRS(sparseM, val, col_ind, row_ptr, n); // преобразование матрицы из одного формата хранения в другой.
																//m.ballocCRScfd = true;
	simplesparsefree(sparseM, n);

#endif

	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
									   //typedef float    ScalarType;
	typedef int Myindextype;

	integer nnu = n;
	integer nna = row_ptr[n];
	printf("nna=%lld %d\n", nna, row_ptr[n - 1]);
	//system("pause");

	/**
	* Set up the matrices and vectors for the iterative solvers (cf. iterative.cpp)
	**/
	Myindextype rows = static_cast<Myindextype>(nnu);
	Myindextype cols = static_cast<Myindextype>(nnu);
	Myindextype nonzeros = static_cast<Myindextype>(nna);

	// Прочитать матрицу:
	std::cout << "Reading matrix..." << std::endl;

	//integer* row_jumper = new integer[nnu + 1];
	//integer* col_buffer = new integer[nna];
	//ScalarType* elements = new ScalarType[nna];

	std::vector<int> row_jumper(nnu + 1); // только тип int !!!
	std::vector<int> col_buffer(nna); // только тип int !!!
	std::vector<double> elements(nna);
	std::vector<double> rhs(nnu);

	for (integer i = 0; i < nnu; i++) {
		rhs[i] = dV[i];
	}

	std::vector<double> x(n, 0);

	if (1) {
		// Инициализация
		for (integer i = 0; i < n; i++) {
			x[i] = dX0[i];
		}
	}

	doublereal alpharelax = 1.0;

	if (!bStructural_Mechanic) {
		// Уравнение теплопераедачи.

		// Это не специальная нелинейная версия кода amg1r5 CAMG.
		for (int k = 0; k < lw; k++) {
			if ((w[k].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
				(w[k].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY)) {
				alpharelax = 0.99999; // Для того чтобы СЛАУ сходилась.
				// 0.9999 - недостаточное значение, температуры не те получаются.
			}
		}
		if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {
			alpharelax = 0.99999;
		}
		//if (adiabatic_vs_heat_transfer_coeff == ADIABATIC_WALL_BC) {
		//printf("ADIABATIC WALL BC"); getchar();
		//}

	}

	row_jumper[nnu] = static_cast<Myindextype>(nna);
	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	integer iscan = 0;
	for (integer i = 0; i < n; i++) {

		row_jumper[i] = static_cast<Myindextype>(row_ptr[i]);


		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
			if (col_ind[j] == i) {
				if ((!bStructural_Mechanic) && ((bondary != nullptr) && (!bondary[i])) && (row_ptr[i + 1] > row_ptr[i] + 1)) {
					// Уравнение теплопередачи.
					// Релаксация к предыдущему значению.
					elements[iscan] = static_cast<ScalarType>(val[j] / alpharelax);
					rhs[i] += (1.0 - alpharelax) * val[j] * x[i] / alpharelax;
				}
				else {
					elements[iscan] = static_cast<ScalarType>(val[j]); // Условие Дирихле.
				}
				col_buffer[iscan] = static_cast<Myindextype>(col_ind[j]);
				iscan++;
				//debug
				//std::cout << "col_buffer[j]=" << col_ind[j] << std::endl;
				//system("pause");
			}
		}
		for (integer j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
			if (col_ind[j] != i) {
				elements[iscan] = static_cast<ScalarType>(val[j]);
				col_buffer[iscan] = static_cast<Myindextype>(col_ind[j]);
				iscan++;
				// debug
				//std::cout << "col_buffer[j]=" << col_ind[j] << std::endl;
				//system("pause");
			}
		}

	}

	row_jumper[n] = static_cast<Myindextype>(row_ptr[n]);

	delete[] val;
	val = nullptr;

	delete[] col_ind;
	col_ind = nullptr;

	delete[] row_ptr;
	row_ptr = nullptr;



	printf("Matrix load succsefull...\n");

	amgclHandle prm = amgcl_params_create();

	// число ячеек на самом грубом уровне.
	amgcl_params_seti(prm, "precond.coarse_enough", 100);

	/*
	switch (my_amg_manager.amgcl_selector) {
	case 0: // Ruge - Stueben (аналог amg1r5)
		amgcl_params_sets(prm, "precond.coarsening.type", "ruge_stuben");
		amgcl_params_setf(prm, "precond.coarsening.eps_strong", 0.9f);
		//{type = ruge_stuben, eps_strong = 0.9}
		break;
	case 1: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		break;
	default: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		break;
	}

	switch (my_amg_manager.amgcl_smoother) {
	case 0: // spai0
		amgcl_params_sets(prm, "precond.relax.type", "spai0");
		break;
	case 1: // ilu0
		amgcl_params_sets(prm, "precond.relax.type", "ilu0");
		break;
	case 2: // gauss_seidel
		amgcl_params_sets(prm, "precond.relax.type", "gauss_seidel");
		break;
	case 3: // damped_jacobi
		amgcl_params_sets(prm, "precond.relax.type", "damped_jacobi");
		amgcl_params_setf(prm, "precond.relax.damping", 0.8f);
		break;
	case 4: // spai1
		amgcl_params_sets(prm, "precond.relax.type", "spai1");
		break;
	case 5: // chebyshev
		amgcl_params_sets(prm, "precond.relax.type", "chebyshev");
		break;
	case 6: // ilu1
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 1);
		break;
	case 7: // ilu2
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 2);
		break;
	case 8: // ilu4
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 4);
		printf("precond.relax.type==ilu(k==4).\n");
		break;
	case 9: // ilu6
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 6);
		printf("precond.relax.type==ilu(k==6).\n");
		break;
	case 10: // ilu8
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 8);
		printf("precond.relax.type==ilu(k==8).\n");
		break;
	case 11: // ilu10
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 10);
		printf("precond.relax.type==ilu(k==10).\n");
		break;
	default:	amgcl_params_sets(prm, "precond.relax.type", "spai0");
		break;
	}
	//amgcl_params_sets(prm, "solver.type", "bicgstabl");
	//amgcl_params_seti(prm, "solver.L", 1);
	switch (my_amg_manager.amgcl_iterator) {
	case AMGCL_ITERATOR_ALG::BiCGStab: // BiCGStab
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		break;
	case AMGCL_ITERATOR_ALG::FGMRes: // FGMRes
		amgcl_params_sets(prm, "solver.type", "fgmres");
		break;
	default:
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		break;
	}

	if (bglobal_unsteady_temperature_determinant) {
		// Нестационарные задачи требуется считать до меньших значений невязки.
		// 14.05.2019
		// 1.0e-12 точность достаточна.
		amgcl_params_setf(prm, "solver.tol", 1.0e-12f);
		// Нестационарные задачи требуется считать до меньших значений невязки.
		// Может локально не сойтись поэтому не надо делать очень большого числа итераций.
		// 14.05.2019
		amgcl_params_seti(prm, "solver.maxiter", 3000);
	}
	*/

	if (bStructural_Mechanic)
	{
		// Механика.
		// 1.0e-8 точность достаточна.
		// число ячеек на самом грубом уровне.
		amgcl_params_seti(prm, "precond.coarse_enough", 100);

		amgcl_params_sets(prm, "precond.coarsening.type", "ruge_stuben");
		amgcl_params_setf(prm, "precond.coarsening.eps_strong", 0.9f);
		//{type = ruge_stuben, eps_strong = 0.9}

		amgcl_params_sets(prm, "precond.relax.type", "gauss_seidel");
		//amgcl_params_sets(prm, "solver.type", "bicgstab");
		amgcl_params_sets(prm, "solver.type", "cg");

		amgcl_params_setf(prm, "solver.tol", 1.0e-8f);
		amgcl_params_seti(prm, "solver.maxiter", 3000);
	}

	printf("Setup phase start...\n");

	bool bprint_preconditioner = true;

	int n_loc = (int)(n);

	/*
// Начинается с нуля всё верно. Проверено. 17,07,2021
for (ptrdiff_t i = 0; i < n_loc; ++i) {
	std::cout << row_jumper[i] << "  " << col_buffer[i] << std::endl;
	getchar();
}*/

	/*
	
	Вы решаете задачу

D^{-1} A y = {D^-1} f.

В примере Serena решается задача 

D^{-1/2} A {D^-1/2} x = {D^-1/2} f.

В вашем случае теряется симметричность матрицы (если она была симметричной). 
В моем случае симметричность сохраняется, но чтобы получить решение исходной системы y, решение масштабированной системы нужно домножить на D^{-1/2}:
y = D^{-1/2}x.

Об этом кстати написано в туториале, под листингом 8.
	*/

// Scale the matrix so that it has the unit diagonal.
   // First, find the diagonal values:
	std::vector<double> D(n_loc, 1.0);
	for (ptrdiff_t i = 0; i < n_loc; ++i) {
		for (ptrdiff_t j = row_jumper[i], e = row_jumper[i + 1]; j < e; ++j) {
			if (col_buffer[j] == i) {
				if (elements[j] <= 0.0) {

					std::cout << elements[j] << " " << i << " " << n_loc << " rhs=" << rhs[j] << std::endl;
					getchar();
				}
				D[i] = 1 / sqrt(fabs(elements[j]));
				//D[i] = 1 / elements[j];
				break;
			}
		}
	}

	// Then, apply the scaling in-place:
	for (ptrdiff_t i = 0; i < n_loc; ++i) {
		for (ptrdiff_t j = row_jumper[i], e = row_jumper[i + 1]; j < e; ++j) {
			elements[j] *= D[i] * D[col_buffer[j]];
		}
		rhs[i] *= D[i];
	}


	//****
	if (bprint_preconditioner) {
		// Библиотека Дениса Демидова работает только с 32 битным 
		// типом int. Его хватает даже для числа неизвестных 80 млн
		// контрольных объемов.

		amgclHandle amg_precond = amgcl_precond_create(
			n_loc, row_jumper.data(), col_buffer.data(), elements.data(), prm
		);
		amgcl_precond_report(amg_precond);//печать samg предобуславливателя.
		amgcl_precond_destroy(amg_precond);
	}
	//*****

	// Библиотека Дениса Демидова работает только с 32 битным 
	// типом int. Его хватает даже для числа неизвестных 80 млн
	// контрольных объемов.
	//int n_loc = (int)(n);
	amgclHandle solver = amgcl_solver_create(
		n_loc, row_jumper.data(), col_buffer.data(), elements.data(), prm
	);



	amgcl_params_destroy(prm);

	printf("Solution phase start...\n");




	conv_info cnv = amgcl_solver_solve(solver, rhs.data(), x.data());

	// Solve same problem again, but explicitly provide the matrix this time:
	//std::fill(x.begin(), x.end(), 0);
	//cnv = amgcl_solver_solve_mtx(
	//solver, row_jumper.data(), col_buffer.data(), elements.data(),
	//rhs.data(), x.data()
	//);

	std::cout << "Iterations: " << cnv.iterations << std::endl
		<< "Error:      " << cnv.residual << std::endl;

	amgcl_solver_destroy(solver);


	for (integer i = 0; i < nnu; i++) {
		dX0[i] = D[i]* x[i];
	}

#ifdef _OPENMP 
	omp_set_num_threads(1); // установка числа потоков
#endif


}


void amgcl_secondT_solver_cpu(SIMPLESPARSE &sparseM, integer n,
	doublereal *dV, doublereal* &dX0, 
	bool bprintmessage, WALL* &w, int &lw, bool* &bondary,
	bool bStructural_Mechanic)
{
	// размерность квадратной матрицы.
	// n

#ifdef _OPENMP 
	// Узнаёт количество ядер в системе.
	// 15млн неизвестных время параллельного кода 21мин 39с.
	// Время однопоточного кода 27мин 48с.
	unsigned int nthreads = number_cores();
	omp_set_num_threads(nthreads); // установка числа потоков
#endif

								   // Разреженная матрица СЛАУ
								   // в CRS формате.
	doublereal* val=nullptr;
	integer_mix_precision* col_ind = nullptr;
	integer_mix_precision* row_ptr = nullptr;


	simplesparsetoCRS(sparseM, val, col_ind, row_ptr, n); // преобразование матрицы из одного формата хранения в другой.
																//m.ballocCRScfd = true;
	simplesparsefree(sparseM, n);


	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
									   //typedef float    ScalarType;
	typedef int Myindextype;

	integer nnu = n;
	integer nna = row_ptr[n];
	printf("nna=%lld %d\n",nna, row_ptr[n-1]);
	//system("pause");

	/**
	* Set up the matrices and vectors for the iterative solvers (cf. iterative.cpp)
	**/
	Myindextype rows = static_cast<Myindextype>(nnu);
	Myindextype cols = static_cast<Myindextype>(nnu);
	Myindextype nonzeros = static_cast<Myindextype>(nna);

	// Прочитать матрицу:
	std::cout << "Reading matrix..." << std::endl;

	//integer* row_jumper = new integer[nnu + 1];
	//integer* col_buffer = new integer[nna];
	//ScalarType* elements = new ScalarType[nna];

	std::vector<int> row_jumper(nnu + 1); // только тип int !!!
	std::vector<int> col_buffer(nna); // только тип int !!!
	std::vector<double> elements(nna);
	std::vector<double> rhs(nnu);

	for (integer i = 0; i < nnu; i++) {
		rhs[i] = dV[i];
	}

	std::vector<double> x(n, 0);

	if (1) {
		// Инициализация
		for (integer i = 0; i < n; i++) {
			x[i] = dX0[i];
		}
	}

	doublereal alpharelax = 1.0;

	if (!bStructural_Mechanic) {
		// Уравнение теплопераедачи.

		// Это не специальная нелинейная версия кода amg1r5 CAMG.
		for (int k = 0; k < lw; k++) {
			if ((w[k].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
				(w[k].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY)) {
				alpharelax = 0.99999; // Для того чтобы СЛАУ сходилась.
				// 0.9999 - недостаточное значение, температуры не те получаются.
			}
		}
		if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {
			alpharelax = 0.99999;
		}
		//if (adiabatic_vs_heat_transfer_coeff == ADIABATIC_WALL_BC) {
		//printf("ADIABATIC WALL BC"); getchar();
		//}

	}

	row_jumper[nnu] = static_cast<Myindextype>(nna);
	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	integer iscan = 0;
	for (integer i = 0; i < n; i++) {
		
			row_jumper[i] = static_cast<Myindextype>(row_ptr[i]);
		
			
			for (integer j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
				if (col_ind[j] == i) {
					if ((!bStructural_Mechanic)&&((bondary!=nullptr)&&(!bondary[i]))&&(row_ptr[i + 1] > row_ptr[i] + 1)) {
						// Уравнение теплопередачи.
						// Релаксация к предыдущему значению.
						elements[iscan] = static_cast<ScalarType>(val[j] / alpharelax);
						rhs[i] += (1.0 - alpharelax)* val[j]*x[i] / alpharelax;
					}
					else {
						elements[iscan] = static_cast<ScalarType>(val[j]); // Условие Дирихле.
					}
					col_buffer[iscan] = static_cast<Myindextype>(col_ind[j]);
					iscan++;
					//debug
					//std::cout << "col_buffer[j]=" << col_ind[j] << std::endl;
					//system("pause");
				}
			}
			for (integer j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
				if (col_ind[j] != i) {
					elements[iscan] = static_cast<ScalarType>(val[j]);
					col_buffer[iscan] = static_cast<Myindextype>(col_ind[j]);
					iscan++;
					// debug
					//std::cout << "col_buffer[j]=" << col_ind[j] << std::endl;
					//system("pause");
				}
			}

	}

	row_jumper[n] = static_cast<Myindextype>(row_ptr[n]);

	delete[] val;
	val = nullptr;
	
	delete[] col_ind;
	col_ind = nullptr;
		
	delete[] row_ptr;
	row_ptr = nullptr;
	


	printf("Matrix load succsefull...\n");

	amgclHandle prm = amgcl_params_create();

	// число ячеек на самом грубом уровне.
	amgcl_params_seti(prm, "precond.coarse_enough", 1000);

	switch (my_amg_manager.amgcl_selector) {
	case 0: // Ruge - Stueben (аналог amg1r5)
		amgcl_params_sets(prm, "precond.coarsening.type", "ruge_stuben");
		amgcl_params_setf(prm, "precond.coarsening.eps_strong", 0.9f);
		//{type = ruge_stuben, eps_strong = 0.9}
		break;
	case 1: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		break;
	default: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		break;
	}

	switch (my_amg_manager.amgcl_smoother) {
	case 0: // spai0
		amgcl_params_sets(prm, "precond.relax.type", "spai0");
		break;
	case 1: // ilu0
		amgcl_params_sets(prm, "precond.relax.type", "ilu0");
		break;
	case 2: // gauss_seidel
		amgcl_params_sets(prm, "precond.relax.type", "gauss_seidel");
		break;
	case 3: // damped_jacobi
		amgcl_params_sets(prm, "precond.relax.type", "damped_jacobi");
		amgcl_params_setf(prm, "precond.relax.damping", 0.8f);
		break;
	case 4: // spai1
		amgcl_params_sets(prm, "precond.relax.type", "spai1");
		break;
	case 5: // chebyshev
		amgcl_params_sets(prm, "precond.relax.type", "chebyshev");
		break;
	case 6: // ilu1
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 1);
		break;
	case 7: // ilu2
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 2);
		break;
	case 8: // ilu4
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 4);
		printf("precond.relax.type==ilu(k==4).\n");
		break;
	case 9: // ilu6
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 6);
		printf("precond.relax.type==ilu(k==6).\n");
		break;
	case 10: // ilu8
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 8);
		printf("precond.relax.type==ilu(k==8).\n");
		break;
	case 11: // ilu10
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 10);
		printf("precond.relax.type==ilu(k==10).\n");
		break;
	default:	amgcl_params_sets(prm, "precond.relax.type", "spai0");
		break;
	}
	//amgcl_params_sets(prm, "solver.type", "bicgstabl");
	//amgcl_params_seti(prm, "solver.L", 1);
	switch (my_amg_manager.amgcl_iterator) {
	case AMGCL_ITERATOR_ALG::BiCGStab: // BiCGStab
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		break;
	case AMGCL_ITERATOR_ALG::FGMRes: // FGMRes
		amgcl_params_sets(prm, "solver.type", "fgmres");
		break;
	default:
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		break;
	}

	if (bglobal_unsteady_temperature_determinant) {
		// Нестационарные задачи требуется считать до меньших значений невязки.
		// 14.05.2019
		// 1.0e-12 точность достаточна.
		amgcl_params_setf(prm, "solver.tol", 1.0e-12f);
		// Нестационарные задачи требуется считать до меньших значений невязки.
		// Может локально не сойтись поэтому не надо делать очень большого числа итераций.
		// 14.05.2019
		amgcl_params_seti(prm, "solver.maxiter", 3000);
	}

	if (bStructural_Mechanic)
	{
		// Механика.
		// 1.0e-8 точность достаточна.
		// число ячеек на самом грубом уровне.
		amgcl_params_seti(prm, "precond.coarse_enough", 45);

		amgcl_params_setf(prm, "solver.tol", 1.0e-8f);
		amgcl_params_seti(prm, "solver.maxiter", 3000);
	}

	printf("Setup phase start...\n");

	bool bprint_preconditioner = true;

	int n_loc = (int)(n);

	/*
// Начинается с нуля всё верно. Проверено. 17,07,2021
for (ptrdiff_t i = 0; i < n_loc; ++i) {
	std::cout << row_jumper[i] << "  " << col_buffer[i] << std::endl;
	getchar();
}*/

// Scale the matrix so that it has the unit diagonal.
   // First, find the diagonal values:
	std::vector<double> D(n_loc, 1.0);
	for (ptrdiff_t i = 0; i < n_loc; ++i) {
		for (ptrdiff_t j = row_jumper[i], e = row_jumper[i + 1]; j < e; ++j) {
			if (col_buffer[j] == i) {
				if (elements[j] <= 0.0) {

					std::cout << elements[j] << " " << i << " " << n_loc << " rhs=" << rhs[j] << std::endl;
					getchar();
				}
				D[i] = 1 / sqrt(fabs(elements[j]));
				//D[i] = 1 / elements[j];
				break;
			}
		}
	}

	// Then, apply the scaling in-place:
	for (ptrdiff_t i = 0; i < n_loc; ++i) {
		for (ptrdiff_t j = row_jumper[i], e = row_jumper[i + 1]; j < e; ++j) {
			elements[j] *= D[i];
		}
		rhs[i] *= D[i];
	}


	//****
	if (bprint_preconditioner) {
		// Библиотека Дениса Демидова работает только с 32 битным 
	    // типом int. Его хватает даже для числа неизвестных 80 млн
	    // контрольных объемов.
		
		amgclHandle amg_precond = amgcl_precond_create(
			n_loc, row_jumper.data(), col_buffer.data(), elements.data(), prm
		);
		amgcl_precond_report(amg_precond);//печать samg предобуславливателя.
		amgcl_precond_destroy(amg_precond);
	}
	//*****

	// Библиотека Дениса Демидова работает только с 32 битным 
	// типом int. Его хватает даже для числа неизвестных 80 млн
	// контрольных объемов.
	//int n_loc = (int)(n);
	amgclHandle solver = amgcl_solver_create(
		n_loc, row_jumper.data(), col_buffer.data(), elements.data(), prm
	);



	amgcl_params_destroy(prm);

	printf("Solution phase start...\n");


	

	conv_info cnv = amgcl_solver_solve(solver, rhs.data(), x.data());

	// Solve same problem again, but explicitly provide the matrix this time:
	//std::fill(x.begin(), x.end(), 0);
	//cnv = amgcl_solver_solve_mtx(
	//solver, row_jumper.data(), col_buffer.data(), elements.data(),
	//rhs.data(), x.data()
	//);

	std::cout << "Iterations: " << cnv.iterations << std::endl
		<< "Error:      " << cnv.residual << std::endl;

	amgcl_solver_destroy(solver);


	for (integer i = 0; i < nnu; i++) {
		dX0[i] = x[i];
	}

#ifdef _OPENMP 
	omp_set_num_threads(1); // установка числа потоков
#endif


}




// Это вызов библиотечного решателя систем линейных алгебраических уравнений.
// Библиотека AMGCL Дениса Демидова. Это библиотека с открытым исходным кодом распространяющаяся 
// по open Source лицензии MIT (X11) license. 
// Начало присоединения 7 мая 2019.
void amgcl_solver_cpu(equation3D* &sl, equation3D_bon* &slb,
	integer maxelm, integer maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit,
	doublereal alpharelax, integer iVar,
	bool bprint_preconditioner, 
	doublereal dgx, doublereal dgy, doublereal dgz, 
	integer inumber_iteration_SIMPLE, WALL* &w, int &lw)
{

	bool bprint_log = false;

	if ((iVar == TEMP) && (maxelm > 3000000)) {
		// Если мы решаем чистую теплопередачу и задача довольно большая.
		bprint_log = true;
	}
	if ((bmyconvective7248) && (iVar == TEMP) && (bSIMPLErun_now_for_temperature==false)) {
		// Теплопередача с конвекцией не cfd расчёт.


		bprint_log = true;
	}

	bprint_log = true;

	// maxit - не используется.
	// bprint_preconditioner==true печать иерархии матриц на консоль.

#ifdef _OPENMP 
// Узнаёт количество ядер в системе.
// 15млн неизвестных время параллельного кода 21мин 39с.
// Время однопоточного кода 27мин 48с.
	unsigned int nthreads = number_cores();
	omp_set_num_threads(nthreads); // установка числа потоков
#endif

	if (deltP_init == nullptr) {
		deltP_init = new doublereal[maxelm+maxbound];
		for (integer i = 0; i < maxelm + maxbound; i++) {
			deltP_init[i] = 0.0;
		}
	}

	if (dX0 == nullptr) {
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
	for (integer i = 0; i < maxelm; i++) {
		// внутренность матрицы.
		if ((sl[i].iP > -1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;
		if (sl[i].ap <= nonzeroEPS) {
			printf("bad diagonal in string %lld ap=%e\n",i, sl[i].ap);
			system("PAUSE");
		}
		if (sl[i].iP == -1) {
			printf("bad diagonal in string %lld iP=%lld\n", i, sl[i].iP);
			system("PAUSE");
		}

		if (b_on_adaptive_local_refinement_mesh) {

			/*
			if (0 && (1.00005 * sl[i].ap < sl[i].ab + sl[i].at + sl[i].ae + sl[i].aw + sl[i].an + sl[i].as +
				sl[i].ab2 + sl[i].at2 + sl[i].ae2 + sl[i].aw2 + sl[i].an2 + sl[i].as2 +
				sl[i].ab3 + sl[i].at3 + sl[i].ae3 + sl[i].aw3 + sl[i].an3 + sl[i].as3 +
				sl[i].ab4 + sl[i].at4 + sl[i].ae4 + sl[i].aw4 + sl[i].an4 + sl[i].as4))
			{
				printf("bad diagonal preobladanie in string %lld ap=%e sum_anb=%e \n", i, sl[i].ap,
					sl[i].ab + sl[i].at + sl[i].ae + sl[i].aw + sl[i].an + sl[i].as +
					sl[i].ab2 + sl[i].at2 + sl[i].ae2 + sl[i].aw2 + sl[i].an2 + sl[i].as2 +
					sl[i].ab3 + sl[i].at3 + sl[i].ae3 + sl[i].aw3 + sl[i].an3 + sl[i].as3 +
					sl[i].ab4 + sl[i].at4 + sl[i].ae4 + sl[i].aw4 + sl[i].an4 + sl[i].as4);

				switch (iVar) {
				case TEMP: printf("TEMP iVar\n");
					break;
				case VELOCITY_X_COMPONENT: 
					printf("VELOCITY_X_COMPONENT iVar\n");
					break;
				case VELOCITY_Y_COMPONENT:
					printf("VELOCITY_Y_COMPONENT iVar\n");
					break;
				case VELOCITY_Z_COMPONENT:
					printf("VELOCITY_Z_COMPONENT iVar\n");
					break;
				case PAM :
					printf("PAM iVar\n");
					break;
				}
				system("PAUSE");
			}
			*/
			if ((sl[i].iB2 > -1) && (fabs(sl[i].ab2) > nonzeroEPS)) (nna)++;
			if (sl[i].ab2 < -nonzeroEPS) {
				printf("bad ab2 in string %lld %e\n", i, -sl[i].ab2);
			}
			if ((sl[i].iE2 > -1) && (fabs(sl[i].ae2) > nonzeroEPS)) (nna)++;
			if (sl[i].ae2 < -nonzeroEPS) {
				printf("bad ae2 in string %lld %e\n", i, -sl[i].ae2);
			}
			if ((sl[i].iN2 > -1) && (fabs(sl[i].an2) > nonzeroEPS)) (nna)++;
			if (sl[i].an2 < -nonzeroEPS) {
				printf("bad an2 in string %lld %e\n", i, -sl[i].an2);
			}
			if ((sl[i].iS2 > -1) && (fabs(sl[i].as2) > nonzeroEPS)) (nna)++;
			if (sl[i].as2 < -nonzeroEPS) {
				printf("bad as2 in string %lld %e\n", i, -sl[i].as2);
			}
			if ((sl[i].iT2 > -1) && (fabs(sl[i].at2) > nonzeroEPS)) (nna)++;
			if (sl[i].at2 < -nonzeroEPS) {
				printf("bad at2 in string %lld %e\n", i, -sl[i].at2);
			}
			if ((sl[i].iW2 > -1) && (fabs(sl[i].aw2) > nonzeroEPS)) (nna)++;
			if (sl[i].aw2 < -nonzeroEPS) {
				printf("bad aw2 in string %lld %e\n", i, -sl[i].aw2);
			}

			if ((sl[i].iB3 > -1) && (fabs(sl[i].ab3) > nonzeroEPS)) (nna)++;
			if (sl[i].ab3 < -nonzeroEPS) {
				printf("bad ab3 in string %lld %e\n", i, -sl[i].ab3);
			}
			if ((sl[i].iE3 > -1) && (fabs(sl[i].ae3) > nonzeroEPS)) (nna)++;
			if (sl[i].ae3 < -nonzeroEPS) {
				printf("bad ae3 in string %lld %e\n", i, -sl[i].ae3);
			}
			if ((sl[i].iN3 > -1) && (fabs(sl[i].an3) > nonzeroEPS)) (nna)++;
			if (sl[i].an3 < -nonzeroEPS) {
				printf("bad an3 in string %lld %e\n", i, -sl[i].an3);
			}
			if ((sl[i].iS3 > -1) && (fabs(sl[i].as3) > nonzeroEPS)) (nna)++;
			if (sl[i].as3 < -nonzeroEPS) {
				printf("bad as3 in string %lld %e\n", i, -sl[i].as3);
			}
			if ((sl[i].iT3 > -1) && (fabs(sl[i].at3) > nonzeroEPS)) (nna)++;
			if (sl[i].at3 < -nonzeroEPS) {
				printf("bad at3 in string %lld %e\n", i, -sl[i].at3);
			}
			if ((sl[i].iW3 > -1) && (fabs(sl[i].aw3) > nonzeroEPS)) (nna)++;
			if (sl[i].aw3 < -nonzeroEPS) {
				printf("bad aw3 in string %lld %e\n", i, -sl[i].aw3);
			}

			if ((sl[i].iB4 > -1) && (fabs(sl[i].ab4) > nonzeroEPS)) (nna)++;
			if (sl[i].ab4 < -nonzeroEPS) {
				printf("bad ab4 in string %lld %e\n", i, -sl[i].ab4);
			}
			if ((sl[i].iE4 > -1) && (fabs(sl[i].ae4) > nonzeroEPS)) (nna)++;
			if (sl[i].ae4 < -nonzeroEPS) {
				printf("bad ae4 in string %lld %e\n", i, -sl[i].ae4);
			}
			if ((sl[i].iN4 > -1) && (fabs(sl[i].an4) > nonzeroEPS)) (nna)++;
			if (sl[i].an4 < -nonzeroEPS) {
				printf("bad an4 in string %lld %e\n", i, -sl[i].an4);
			}
			if ((sl[i].iS4 > -1) && (fabs(sl[i].as4) > nonzeroEPS)) (nna)++;
			if (sl[i].as4 < -nonzeroEPS) {
				printf("bad as4 in string %lld %e\n", i, -sl[i].as4);
			}
			if ((sl[i].iT4 > -1) && (fabs(sl[i].at4) > nonzeroEPS)) (nna)++;
			if (sl[i].at4 < -nonzeroEPS) {
				printf("bad at4 in string %lld %e\n", i, -sl[i].at4);
			}
			if ((sl[i].iW4 > -1) && (fabs(sl[i].aw4) > nonzeroEPS)) (nna)++;
			if (sl[i].aw4 < -nonzeroEPS) {
				printf("bad aw4 in string %lld %e\n", i, -sl[i].aw4);
			}

		}
		else {
		/*
			if (0&&((iTEMPScheme == UDS)) && (1.00005 * sl[i].ap < sl[i].ab + sl[i].at + sl[i].ae + sl[i].aw + sl[i].an + sl[i].as))
			{
				printf("bad diagonal preobladanie in string %lld ap=%e sum_anb=%e \n", i, sl[i].ap,
					sl[i].ab + sl[i].at + sl[i].ae + sl[i].aw + sl[i].an + sl[i].as);

				switch (iVar) {
				case TEMP: printf("TEMP iVar\n");
					break;
				case VELOCITY_X_COMPONENT:
					printf("VELOCITY_X_COMPONENT iVar\n");
					break;
				case VELOCITY_Y_COMPONENT:
					printf("VELOCITY_Y_COMPONENT iVar\n");
					break;
				case VELOCITY_Z_COMPONENT:
					printf("VELOCITY_Z_COMPONENT iVar\n");
					break;
				case PAM:
					printf("PAM iVar\n");
					break;
				}
				system("PAUSE");
			}
			*/
		}

		if ((sl[i].iB > -1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
		if (sl[i].ab < -nonzeroEPS) {
			printf("bad ab in string %lld %e\n", i, -sl[i].ab);
		}
		if ((sl[i].iE > -1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
		if (sl[i].ae < -nonzeroEPS) {
			printf("bad ae in string %lld %e\n", i, -sl[i].ae);
		}
		if ((sl[i].iN > -1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
		if (sl[i].an < -nonzeroEPS) {
			printf("bad an in string %lld %e\n", i, -sl[i].an);
		}
		if ((sl[i].iS > -1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
		if (sl[i].as < -nonzeroEPS) {
			printf("bad as in string %lld %e\n", i, -sl[i].as);
		}
		if ((sl[i].iT > -1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
		if (sl[i].at < -nonzeroEPS) {
			printf("bad at in string %lld %e\n", i, -sl[i].at);
		}
		if ((sl[i].iW > -1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;
		if (sl[i].aw < -nonzeroEPS) {
			printf("bad aw in string %lld %e\n", i, -sl[i].aw);
		}

	}
	for (integer i = 0; i < maxbound; i++) {
		// граничные узлы.
		if ((slb[i].iW > -1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		if (slb[i].aw <= nonzeroEPS) {
			printf("bad diagonal in string %lld ap=%e maxelm=%lld\n", maxelm + i, slb[i].aw,maxelm);
			system("PAUSE");
		}
		if (slb[i].iW == -1) {
			printf("bad diagonal in string %lld iP=%lld maxelm=%lld\n", maxelm + i, slb[i].iW, maxelm);
			system("PAUSE");
		}
		if ((slb[i].iI > -1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
		if (slb[i].ai < -nonzeroEPS) {
			printf("bad ai in string %lld ai=%e maxelm=%lld\n", maxelm + i, -slb[i].ai, maxelm);
			system("PAUSE");
		}
	}

	integer nnu = 0; // число неизвестных.
	nnu = maxelm + maxbound;

	/**
	* Printeger some device info at the beginning. If there is more than one OpenCL device available, use the second device.
	**/
	// Мы не используем видеокарту, а используем один поток центрального процессора.
	//std::cout << std::endl;
	//std::cout << "----------------------------------------------" << std::endl;
	//std::cout << "               Device Info" << std::endl;
	//std::cout << "----------------------------------------------" << std::endl;




	typedef doublereal    ScalarType;  // feel free to change this to double if supported by your device
									   //typedef float    ScalarType;
	typedef int Myindextype;


									   /**
									   * Set up the matrices and vectors for the iterative solvers (cf. iterative.cpp)
									   **/
	Myindextype rows = static_cast<Myindextype>(nnu);
	Myindextype cols = static_cast<Myindextype>(nnu);
	Myindextype nonzeros = static_cast<Myindextype>(nna);


	// Прочитать матрицу:
	if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
		if (bprint_log) {
			std::cout << "Reading matrix..." << std::endl;
		}
	}
	
	//integer* row_jumper = new integer[nnu + 1];
	//integer* col_buffer = new integer[nna];
	//ScalarType* elements = new ScalarType[nna];

	std::vector<int> row_jumper(nnu + 1); // только тип int !!!
	std::vector<int> col_buffer(nna); // только тип int !!!
	std::vector<double> elements(nna);
	std::vector<double> rhs(nnu);

	for (integer i = 0; i < nnu; i++) {
		rhs[i] = dV[i];
	}

	doublereal alpharelax1 = alpharelax; // Запоминаем первоначальное значение.
	
	if (iVar == TEMP) {

		// Это не специальная нелинейная версия кода amg1r5 CAMG.
		for (int k = 0; k < lw; k++) {
			if ((w[k].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
				(w[k].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY)) {
				alpharelax = 0.99999; // Для того чтобы СЛАУ сходилась.
									  // 0.9999 - недостаточное значение, температуры не те получаются.
				//printf("incomming\n");
			}
		}
		if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {
			alpharelax = 0.99999; //free_debug_parametr1; // 0.99999;
			//printf("temperature alharelax = %e\n", free_debug_parametr1);
		}
		//if (adiabatic_vs_heat_transfer_coeff == ADIABATIC_WALL_BC) {
			//printf("ADIABATIC WALL BC"); getchar();
		//}

		if (bSIMPLErun_now_for_natural_convection) {
			//alpharelax = alpharelax1;
			//alpharelax = 0.99999;// 9995;// 9;
			//alpharelax = 1.0;
		}

	}
	
	//if (inumber_iteration_SIMPLE <= 5) 
	{
		if (iVar == PAM) {

			if (bSIMPLErun_now_for_natural_convection) {
				//alpharelax = 0.99999;// 999;// 99;
				//alpharelax = 0.3;
			}
		}
	}
	

	//alpharelax = 0.99999;

	if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
		if (bprint_log) {
			printf("alpharelax=%e\n", alpharelax);
			//system("PAUSE");
		}
	}

	

	row_jumper[nnu] = static_cast<Myindextype>(nna);
	// initialize matrix entries on host
	nna = 0;
	//Ah.row_indices[0] = 0; Ah.column_indices[0] = 0; Ah.values[0] = 10.0; // demo interface
	for (integer i = 0; i < maxelm; i++) {
		row_jumper[i] = static_cast<Myindextype>(nna);

		// внутренность матрицы.
		if ((sl[i].iP > -1) && (fabs(sl[i].ap) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iP; Ah.values[nna] = sl[i].ap / alpharelax;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iP) = static_cast<ScalarType> (sl[i].ap / alpharelax);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iP);
			if (/*(bSIMPLErun_now_for_natural_convection==false)&&*/(iVar == TEMP)) {
				
				if ((bmyconvective7248)||(bSIMPLErun_now_for_natural_convection)) {
					// Натуральная конвекция.
					//rhs[i] += (1 - alpharelax) * sl[i].ap * dX0[i] / alpharelax;
					elements[nna] = static_cast<ScalarType> (sl[i].ap / alpharelax);
					//elements[nna] = static_cast<ScalarType> (sl[i].ap);
				}
				else {
					rhs[i] += (1 - alpharelax) * sl[i].ap * dX0[i] / alpharelax;
					elements[nna] = static_cast<ScalarType> (sl[i].ap / alpharelax);
				}
			}
			else if (/*(inumber_iteration_SIMPLE <= 5)&&*/(iVar == PAM)) 
			{
				//if (iVar == PAM) {
					//rhs[i] += (1 - alpharelax) * sl[i].ap * /*dX0[i]*/deltP_init[i] / alpharelax;
					//elements[nna] = static_cast<ScalarType> (sl[i].ap / alpharelax);
				//}
				//else {
					//rhs[i] += (1 - alpharelax) * sl[i].ap * dX0[i]/*deltP_init[i]*/ / alpharelax;
					//elements[nna] = static_cast<ScalarType> (sl[i].ap / alpharelax);
					//elements[nna] = static_cast<ScalarType> (sl[i].ap );
				//}
				
				elements[nna] = static_cast<ScalarType> (sl[i].ap);
				//rhs[i] += (1 - alpharelax) * sl[i].ap * deltP_init[i] / alpharelax;
				//elements[nna] = static_cast<ScalarType> (sl[i].ap / alpharelax);
				
			}
			else if ((iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA)||(iVar== TURBULENT_KINETIK_ENERGY)||
				(iVar== GAMMA_LANGTRY_MENTER)||(iVar== RE_THETA_LANGTRY_MENTER)) {
				// Справляется только метод сглаженной аггрегации и нижняя релаксация безусловно нужна.
				// С нижней релаксацией сходимость монотоная без всплексов.
				elements[nna] = static_cast<ScalarType> (sl[i].ap / alpharelax);
			}
			else
			{
				elements[nna] = static_cast<ScalarType> (sl[i].ap/ alpharelax);
			}
			(nna)++;
		}
		if ((sl[i].iB > -1) && (fabs(sl[i].ab) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB; Ah.values[nna] = -sl[i].ab;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB) = static_cast<ScalarType> (-sl[i].ab);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iB);
			elements[nna] = static_cast<ScalarType> (-sl[i].ab);
			(nna)++;
		}
		if ((sl[i].iE > -1) && (fabs(sl[i].ae) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE; Ah.values[nna] = -sl[i].ae;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE) = static_cast<ScalarType> (-sl[i].ae);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iE);
			elements[nna] = static_cast<ScalarType> (-sl[i].ae);
			(nna)++;
		}
		if ((sl[i].iN > -1) && (fabs(sl[i].an) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN; Ah.values[nna] = -sl[i].an;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN) = static_cast<ScalarType> (-sl[i].an);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iN);
			elements[nna] = static_cast<ScalarType> (-sl[i].an);
			(nna)++;
		}
		if ((sl[i].iS > -1) && (fabs(sl[i].as) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS; Ah.values[nna] = -sl[i].as;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS) = static_cast<ScalarType> (-sl[i].as);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iS);
			elements[nna] = static_cast<ScalarType> (-sl[i].as);
			(nna)++;
		}
		if ((sl[i].iT > -1) && (fabs(sl[i].at) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT; Ah.values[nna] = -sl[i].at;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT) = static_cast<ScalarType> (-sl[i].at);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iT);
			elements[nna] = static_cast<ScalarType> (-sl[i].at);
			(nna)++;
		}
		if ((sl[i].iW > -1) && (fabs(sl[i].aw) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW; Ah.values[nna] = -sl[i].aw;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW) = static_cast<ScalarType> (-sl[i].aw);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iW);
			elements[nna] = static_cast<ScalarType> (-sl[i].aw);
			(nna)++;
		}

		if ((sl[i].iB2 > -1) && (fabs(sl[i].ab2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB2; Ah.values[nna] = -sl[i].ab2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB2) = static_cast<ScalarType> (-sl[i].ab2);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iB2);
			elements[nna] = static_cast<ScalarType> (-sl[i].ab2);
			(nna)++;
		}
		if ((sl[i].iE2 > -1) && (fabs(sl[i].ae2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE2; Ah.values[nna] = -sl[i].ae2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE2) = static_cast<ScalarType> (-sl[i].ae2);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iE2);
			elements[nna] = static_cast<ScalarType> (-sl[i].ae2);
			(nna)++;
		}
		if ((sl[i].iN2 > -1) && (fabs(sl[i].an2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN2; Ah.values[nna] = -sl[i].an2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN2) = static_cast<ScalarType> (-sl[i].an2);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iN2);
			elements[nna] = static_cast<ScalarType> (-sl[i].an2);
			(nna)++;
		}
		if ((sl[i].iS2 > -1) && (fabs(sl[i].as2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS2; Ah.values[nna] = -sl[i].as2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS2) = static_cast<ScalarType>( -sl[i].as2);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iS2);
			elements[nna] = static_cast<ScalarType> (-sl[i].as2);
			(nna)++;
		}
		if ((sl[i].iT2 > -1) && (fabs(sl[i].at2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT2; Ah.values[nna] = -sl[i].at2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT2) = static_cast<ScalarType> (-sl[i].at2);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iT2);
			elements[nna] = static_cast<ScalarType> (-sl[i].at2);
			(nna)++;
		}
		if ((sl[i].iW2 > -1) && (fabs(sl[i].aw2) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW2; Ah.values[nna] = -sl[i].aw2;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW2) = static_cast<ScalarType> (-sl[i].aw2);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iW2);
			elements[nna] = static_cast<ScalarType> (-sl[i].aw2);
			(nna)++;
		}

		if ((sl[i].iB3 > -1) && (fabs(sl[i].ab3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB3; Ah.values[nna] = -sl[i].ab3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB3) = static_cast<ScalarType>(-sl[i].ab3);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iB3);
			elements[nna] = static_cast<ScalarType> (-sl[i].ab3);
			(nna)++;
		}
		if ((sl[i].iE3 > -1) && (fabs(sl[i].ae3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE3; Ah.values[nna] = -sl[i].ae3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE3) = static_cast<ScalarType>(-sl[i].ae3);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iE3);
			elements[nna] = static_cast<ScalarType> (-sl[i].ae3);
			(nna)++;
		}
		if ((sl[i].iN3 > -1) && (fabs(sl[i].an3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN3; Ah.values[nna] = -sl[i].an3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN3) = static_cast<ScalarType> (-sl[i].an3);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iN3);
			elements[nna] = static_cast<ScalarType> (-sl[i].an3);
			(nna)++;
		}
		if ((sl[i].iS3 > -1) && (fabs(sl[i].as3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS3; Ah.values[nna] = -sl[i].as3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS3) = static_cast<ScalarType>(-sl[i].as3);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iS3);
			elements[nna] = static_cast<ScalarType> (-sl[i].as3);
			(nna)++;
		}
		if ((sl[i].iT3 > -1) && (fabs(sl[i].at3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT3; Ah.values[nna] = -sl[i].at3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT3) = static_cast<ScalarType>(-sl[i].at3);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iT3);
			elements[nna] = static_cast<ScalarType> (-sl[i].at3);
			(nna)++;
		}
		if ((sl[i].iW3 > -1) && (fabs(sl[i].aw3) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW3; Ah.values[nna] = -sl[i].aw3;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW3) = static_cast<ScalarType>(-sl[i].aw3);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iW3);
			elements[nna] = static_cast<ScalarType> (-sl[i].aw3);
			(nna)++;
		}

		if ((sl[i].iB4 > -1) && (fabs(sl[i].ab4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iB4; Ah.values[nna] = -sl[i].ab4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iB4) = static_cast<ScalarType>(-sl[i].ab4);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iB4);
			elements[nna] = static_cast<ScalarType> (-sl[i].ab4);
			(nna)++;
		}
		if ((sl[i].iE4 > -1) && (fabs(sl[i].ae4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iE4; Ah.values[nna] = -sl[i].ae4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iE4) = static_cast<ScalarType>(-sl[i].ae4);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iE4);
			elements[nna] = static_cast<ScalarType> (-sl[i].ae4);
			(nna)++;
		}
		if ((sl[i].iN4 > -1) && (fabs(sl[i].an4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iN4; Ah.values[nna] = -sl[i].an4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iN4) = static_cast<ScalarType>(-sl[i].an4);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iN4);
			elements[nna] = static_cast<ScalarType> (-sl[i].an4);
			(nna)++;
		}
		if ((sl[i].iS4 > -1) && (fabs(sl[i].as4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iS4; Ah.values[nna] = -sl[i].as4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iS4) = static_cast<ScalarType>(-sl[i].as4);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iS4);
			elements[nna] = static_cast<ScalarType> (-sl[i].as4);
			(nna)++;
		}
		if ((sl[i].iT4 > -1) && (fabs(sl[i].at4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iT4; Ah.values[nna] = -sl[i].at4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iT4) = static_cast<ScalarType>(-sl[i].at4);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iT4);
			elements[nna] = static_cast<ScalarType> (-sl[i].at4);
			(nna)++;
		}
		if ((sl[i].iW4 > -1) && (fabs(sl[i].aw4) > nonzeroEPS)) {
			//Ah.row_indices[nna] = sl[i].iP; Ah.column_indices[nna] = sl[i].iW4; Ah.values[nna] = -sl[i].aw4;
			//vcl_compressed_matrix1(sl[i].iP, sl[i].iW4) = static_cast<ScalarType>(-sl[i].aw4);
			col_buffer[nna] = static_cast<Myindextype> (sl[i].iW4);
			elements[nna] = static_cast<ScalarType> (-sl[i].aw4);
			(nna)++;
		}


	}

	for (integer i = 0; i < maxbound; i++) {
		row_jumper[maxelm + i] = static_cast<Myindextype>(nna);
		// граничные узлы.
		if ((slb[i].iW > -1) && (fabs(slb[i].aw) > nonzeroEPS)) {
			//Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iW; Ah.values[nna] = slb[i].aw;
			//vcl_compressed_matrix1(slb[i].iW, slb[i].iW) = static_cast<ScalarType> (slb[i].aw);
			col_buffer[nna] = static_cast<Myindextype> (slb[i].iW);
			elements[nna] = static_cast<ScalarType> (slb[i].aw);
			(nna)++;
		}
		if ((slb[i].iI > -1) && (fabs(slb[i].ai) > nonzeroEPS)) {
			//Ah.row_indices[nna] = slb[i].iW; Ah.column_indices[nna] = slb[i].iI; Ah.values[nna] = -slb[i].ai;
			//vcl_compressed_matrix1(slb[i].iW, slb[i].iI) = static_cast<ScalarType>(-slb[i].ai);
			col_buffer[nna] = static_cast<Myindextype> (slb[i].iI);
			elements[nna] = static_cast<ScalarType> (-slb[i].ai);
			(nna)++;
		}
	}

	if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
		if (bprint_log) {
			printf("Matrix load succsefull...\n");
		}
	}

	int n = static_cast<Myindextype>(nnu);
	

	amgclHandle prm = amgcl_params_create();

	// число ячеек на самом грубом уровне.
	amgcl_params_seti(prm, "precond.coarse_enough", 1000);
	
	switch (my_amg_manager.amgcl_selector) {
	case 0: // Ruge - Stueben (аналог amg1r5)
		amgcl_params_sets(prm, "precond.coarsening.type", "ruge_stuben");
		amgcl_params_setf(prm, "precond.coarsening.eps_strong", 0.9f);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.coarsening.type==ruge_stuben.  precond.coarsening.eps_strong=0.9f\n");
			}
		}
		//{type = ruge_stuben, eps_strong = 0.9}
		break;
	case 1: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.coarsening.type==smoothed_aggregation. precond.coarsening.aggr.eps_strong=1e-3f.\n");
			}
		}
		break;
	default: // smoothed aggregation
		amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
		amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.coarsening.type==smoothed_aggregation. precond.coarsening.aggr.eps_strong=1e-3f.\n");
			}
		}
		break;
	}
	
	switch (my_amg_manager.amgcl_smoother) {
	case 0: // spai0
		amgcl_params_sets(prm, "precond.relax.type", "spai0");
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==spai0. \n");
			}
		}
		break;
	case 1: // ilu0
		amgcl_params_sets(prm, "precond.relax.type", "ilu0");
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==ilu0. \n");
			}
		}
		break;
	case 2: // gauss_seidel
		amgcl_params_sets(prm, "precond.relax.type", "gauss_seidel");
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==gauss_seidel\n");
			}
		}
		break;
	case 3: // damped_jacobi
		amgcl_params_sets(prm, "precond.relax.type", "damped_jacobi");
		amgcl_params_setf(prm, "precond.relax.damping", 0.8f);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==damped_jacobi. precond.relax.damping=0.8f\n");
			}
		}
		break;
	case 4: // spai1
		amgcl_params_sets(prm, "precond.relax.type", "spai1");
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==spai1.\n");
			}
		}
		break;
	case 5: // chebyshev
		amgcl_params_sets(prm, "precond.relax.type", "chebyshev");
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==chebyshev.\n");
			}
		}
		break;
	case 6: // ilu1
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 1);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==ilu(k==1).\n");
			}
		}
		break;
	case 7: // ilu2
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 2);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==ilu(k==2).\n");
			}
		}
		break;
	case 8: // ilu4
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 4);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==ilu(k==4).\n");
			}
		}
		break;
	case 9: // ilu6
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 6);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==ilu(k==6).\n");
			}
		}
		break;
	case 10: // ilu8
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 8);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==ilu(k==8).\n");
			}
		}
		break;
	case 11: // ilu10
		amgcl_params_sets(prm, "precond.relax.type", "iluk");
		amgcl_params_seti(prm, "precond.relax.k", 10);
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==ilu(k==10).\n");
			}
		}
		break;
	default:	amgcl_params_sets(prm, "precond.relax.type", "spai0");
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("precond.relax.type==spai0. \n");
			}
		}
		break;
	}
	//amgcl_params_sets(prm, "solver.type", "bicgstabl");
	//amgcl_params_seti(prm, "solver.L", 1);
	switch (my_amg_manager.amgcl_iterator) {
	case AMGCL_ITERATOR_ALG::BiCGStab: // BiCGStab
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("solver.type==bicgstab.\n");
			}
		}
		break;
	case AMGCL_ITERATOR_ALG::FGMRes: // FGMRes
		amgcl_params_sets(prm, "solver.type", "fgmres");
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("solver.type==fgmres.\n");
			}
		}
		break;
	default:
		amgcl_params_sets(prm, "solver.type", "bicgstab");
		if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
			if (bprint_log) {
				printf("solver.type=bicgstab.\n");
			}
		}
		break;
	}
	

	if (bglobal_unsteady_temperature_determinant) {
		// Нестационарные задачи требуется считать до меньших значений невязки.
		// 14.05.2019
		// 1.0e-12 точность достаточна.
		if (my_amg_manager.amgcl_iterator == AMGCL_ITERATOR_ALG::FGMRes) {
			//fgmres
			// 29.12.2019 проверено что достаточно уровня невязки 1Е-7
			// при использовании fgmres для нестационарного моделирования на АЛИС сетке.
			// Выигрыш 3м 30с против 4м 8с.
			amgcl_params_setf(prm, "solver.tol", 1.0e-7f);
		}
		else {
			amgcl_params_setf(prm, "solver.tol", 1.0e-12f);
		}
	}

	if (iVar == PAM) {
		// Поправка давления.
		/*if (bSIMPLErun_now_for_natural_convection) {
			if (inumber_iteration_SIMPLE < 3) {
				amgcl_params_seti(prm, "solver.maxiter", 3);
			}
			else {
				amgcl_params_seti(prm, "solver.maxiter", 30);
			}
		}
		else
		*/{
			amgcl_params_seti(prm, "solver.maxiter", 1000);
			amgcl_params_setf(prm, "solver.tol", 1.0e-6f);
		}
		//if (inumber_iteration_SIMPLE < 20) {
		// Не получается подобны образом для amgcl добиться 
		// сходимости при решении системы уравнений Навье-Стокса
		// в задаче моделирования естественной конвекции с 
		// opening границами.
			//amgcl_params_seti(prm, "solver.maxiter", 2);
		//}
	}
	else if (iVar == TEMP) {
		// Температура.
		// Для задач большой размерности
		if (bglobal_unsteady_temperature_determinant) {
			// Нестационарные задачи требуется считать до меньших значений невязки.
			// Может локально не сойтись поэтому не надо делать очень большого числа итераций.
			// 14.05.2019
			amgcl_params_seti(prm, "solver.maxiter", 1300);
		}
		else {
			if ((fabs(dgx) > 1.0e-20) || (fabs(dgy) > 1.0e-20) || (fabs(dgz) > 1.0e-20)) {
				if (bSIMPLErun_now_for_temperature) {
					// Натуральная конвекция в cfd.
					//amgcl_params_setf(prm, "solver.tol", 1.0e-3f);
					if (bSIMPLErun_now_for_natural_convection) {
						amgcl_params_seti(prm, "solver.maxiter", 100);
						amgcl_params_setf(prm, "solver.tol", 1.0e-9f);
					}
					else {
						amgcl_params_seti(prm, "solver.maxiter", 1000);
						amgcl_params_setf(prm, "solver.tol", 1.0e-12f);
					}
					
				}
			}
			else {
				amgcl_params_seti(prm, "solver.maxiter", 3000);				
			}
		}
	}
	else {
		
		/*if ((bSIMPLErun_now_for_natural_convection) && ((iVar == VELOCITY_X_COMPONENT) ||
			(iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT))) {
			amgcl_params_setf(prm, "solver.tol", 1.0e-20f);
			if (inumber_iteration_SIMPLE < 10) {
				amgcl_params_seti(prm, "solver.maxiter", 1);
			}
			else {
				amgcl_params_seti(prm, "solver.maxiter", 100);
			}
		}
		else*/ {

			// K-Omega SST Turbulence Model
			if (iVar == TURBULENT_KINETIK_ENERGY) {
				if (my_amg_manager.amgcl_selector == 0) {
					// Руге Штубен.
					// Работает только с методом сглаженной агрегации.
					amgcl_params_setf(prm, "solver.tol", 1.0e-5f);
					amgcl_params_seti(prm, "solver.maxiter", 100);
					amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
					amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
					if (bprint_log) {
						printf("TURBULENCE k-omega redirect to precond.coarsening.type==smoothed_aggregation. precond.coarsening.aggr.eps_strong=1e-3f.\n");
					}

				}
				else {
					amgcl_params_setf(prm, "solver.tol", 1.0e-5f);
					amgcl_params_seti(prm, "solver.maxiter", 100);
				}
			}
			else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
				if (my_amg_manager.amgcl_selector == 0) {
					// Руге Штубен.
					// Работает только с методом сглаженной агрегации.
					amgcl_params_setf(prm, "solver.tol", 1.0e-5f);//1.0e-20f
					amgcl_params_seti(prm, "solver.maxiter", 100);
					amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
					amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
					if (bprint_log) {
						printf("TURBULENCE k-omega redirect to precond.coarsening.type==smoothed_aggregation. precond.coarsening.aggr.eps_strong=1e-3f.\n");
					}
				}
				else {
					amgcl_params_setf(prm, "solver.tol", 1.0e-5f);//1.0e-20f
					amgcl_params_seti(prm, "solver.maxiter", 100);
				}
			}
			else if (iVar == GAMMA_LANGTRY_MENTER) {
				if (my_amg_manager.amgcl_selector == 0) {
					// Руге Штубен.
					// Работает только с методом сглаженной агрегации.
					amgcl_params_setf(prm, "solver.tol", 1.0e-5f);
					amgcl_params_seti(prm, "solver.maxiter", 100);
					amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
					amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
					if (bprint_log) {
						printf("TURBULENCE k-omega redirect to precond.coarsening.type==smoothed_aggregation. precond.coarsening.aggr.eps_strong=1e-3f.\n");
					}

				}
				else {
					amgcl_params_setf(prm, "solver.tol", 1.0e-5f);
					amgcl_params_seti(prm, "solver.maxiter", 100);
				}
			}
			else if (iVar == RE_THETA_LANGTRY_MENTER) {
				if (my_amg_manager.amgcl_selector == 0) {
					// Руге Штубен.
					// Работает только с методом сглаженной агрегации.
					amgcl_params_setf(prm, "solver.tol", 1.0e-4f);
					amgcl_params_seti(prm, "solver.maxiter", 100);
					amgcl_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
					amgcl_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
					if (bprint_log) {
						printf("TURBULENCE k-omega redirect to precond.coarsening.type==smoothed_aggregation. precond.coarsening.aggr.eps_strong=1e-3f.\n");
					}

				}
				else {
					amgcl_params_setf(prm, "solver.tol", 1.0e-5f);
					amgcl_params_seti(prm, "solver.maxiter", 100);
				}
			}
			else {

				// Velocity or Spalart Allmares Turbulence model.
				if ((iVar == VELOCITY_X_COMPONENT) ||
					(iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) {
					// Velocity component
					amgcl_params_setf(prm, "solver.tol", 1.0e-6f);
				}
				amgcl_params_seti(prm, "solver.maxiter", 100);
			}
		}

	}

	// Итерируем тривиальное решение при нулевой правой части и нулевых граничных условиях.
	// 24.01.2021
	amgcl_params_seti(prm, "solver.ns_search", 1);

	if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
		if (bprint_log) {
			printf("Setup phase start...\n");


			//****
			if (bprint_preconditioner) {
				amgclHandle amg_precond = amgcl_precond_create(
					n, row_jumper.data(), col_buffer.data(), elements.data(), prm
				);
				amgcl_precond_report(amg_precond);//печать samg предобуславливателя.
				amgcl_precond_destroy(amg_precond);
			}
			//*****
		}
	}

	amgclHandle solver = amgcl_solver_create(
		n, row_jumper.data(), col_buffer.data(), elements.data(), prm
	);	

	
	
	amgcl_params_destroy(prm);

	if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
		if (bprint_log) {
			printf("Solution phase start...\n");
		}
	}

	
	std::vector<double> x(n, 0);

	if (iVar != PAM) {
		// Инициализация
		for (integer i = 0; i < maxelm + maxbound; i++) {
			x[i] = dX0[i];
		}
	}

	conv_info cnv = amgcl_solver_solve(solver, rhs.data(), x.data());

	// Solve same problem again, but explicitly provide the matrix this time:
	//std::fill(x.begin(), x.end(), 0);
	//cnv = amgcl_solver_solve_mtx(
		//solver, row_jumper.data(), col_buffer.data(), elements.data(),
		//rhs.data(), x.data()
	//);

	if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
		if (bprint_log) {
			std::cout << "Iterations: " << cnv.iterations << std::endl
				<< "Error:      " << cnv.residual << std::endl;
		}
	}

	if (cnv.residual < 0.1) {

		for (integer i = 0; i < nnu; i++) {
			dX0[i] = x[i];
		}

		if (iVar == PAM) {
			for (integer i = 0; i < maxelm + maxbound; i++) {
				deltP_init[i] = x[i];
			}
		}
	}


	amgcl_solver_destroy(solver);

	

#ifdef _OPENMP 
	omp_set_num_threads(1); // установка числа потоков
#endif

}

#endif
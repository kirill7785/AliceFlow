// classic_aglomerative_amg6_2018year.cpp
// Очищенный от лишнего кода алгебраический многосеточный метод РУМБА_v0_14.
// Код причёсан чтобы выглядеть компактнее и читабельнее.
// Данная версия classic_aglomerative_amg6() amg решателя отделилась от файла 
// my_agregat_amg.cpp основана версии amg алгоритма 
// classic_aglomerative_amg4(). Версия amg алгоритма
// classic_aglomerative_amg5() стала тупиковой и больше мной
// давно не поддерживалась. Я поддерживал какое то время 
// после версии classic_aglomerative_amg5() версию
// classic_aglomerative_amg4(). Потом от версии 
// classic_aglomerative_amg4() отделился файл
// classic_aglomerative_amg6_2018year.cpp в котором 
// развивается версия classic_aglomerative_amg6().
// Текущая развиваемая и поддерживаемая версия это
// classic_aglomerative_amg6() расположенная в 
// файле classic_aglomerative_amg6_2018year.cpp.
// Версии amg решателя classic_aglomerative_amg5() и 
// classic_aglomerative_amg6() а также существенно более
// ранние версии в файле my_agregat_amg.cpp не поддерживаются и
// их код оставлен для истории.

#pragma once
#ifndef CLASSIC_AGLOMERATIVE_AMG6_2018YEAR_CPP
#define CLASSIC_AGLOMERATIVE_AMG6_2018YEAR_CPP 1

// для функции setw().
//#include <conio.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <random>

#include "basic_interpolation_procedure_my_agregat_amg.cpp" // interpolation procedure for Ak2 type
#include "basic_functions_my_agregat_amg.cpp"

// Дерево ван Эмде Боаса.
#define VEB_FLAG 0 

#if VEB_FLAG
// Исправлено 21.03.2019
// Был конфликт имён члена класса min с функцией min(a,b) языка СИ.
// Аналогично для поля класса max. Поля классов min и max переименованы 
// в данные класса veb_min, veb_max. Теперь конфликт имён отсутствует.
// Дерево ван Эмде Боаса.
// Дерево Ван Эмде Боаса НАЧАЛО РЕАЛИЗАЦИИ 30.06.2018 - окончание 21.09.2018
// Все операции за log(log(U))
//#include "veb.h"
#include "veb.cpp"
#endif

// Строки с отрицательной диагональю запоминаем,
// домножаем на минус 1.0. При итерировании
// домножаем правую часть на минус 1.0.
typedef struct TBAD_STRING_PATCHING {
	integer ilevel;
	integer istring_number;
} BAD_STRING_PATCHING;

// Алгоритмы разреженного умножения матрицы на матрицу.
#include "spgemm_Matrix_by_Matrix_sparse_multiplication.cpp"



// row_ptr - один из трёх массивов составляющих
// хранение матрицы в Column Row Storage формате.
// два других col_ind - номера столбцов, и val -
// сами ненулевые коэффициенты матрицы.
// row_ptr - это номер начала новой строки в CRS формате
// представления матрицы в котором матрица СЛАУ хранится построчно.
// Функция calculate_row_ptr вычисляет массив начала каждой строки матрицы
// row_ind_SA, а также массив конца каждой строки row_ind_EA матрицы Amat
// хранящейся по строкам в Ak2 типе хранения.
void calculate_row_ptr(integer istart, integer iend,
	integer* &row_ind_SA, integer* & row_ind_EA,
	bool* &flag, const Ak2& Amat)
{

	integer i_size_75 = 0;
	// Это нельзя распараллелить.
	for (integer ii = istart; ii <= iend; ++ii) {
		// cnd - condition.
		bool cnd = (flag[Amat.i[ii]] == false);
		if (!cnd) continue;
		
		row_ind_SA[Amat.i[ii]] = ii;
		flag[Amat.i[ii]] = true;
		i_size_75++;
		
	}

	integer ilen = (i_size_75 + 1);

#pragma omp parallel for
	for (integer istr = 1; istr < ilen; ++istr)
	{
		integer kf = row_ind_SA[istr];
		while ((kf <= iend) && (Amat.i[kf] == istr)) {
			kf++;
		}
		kf--;
		row_ind_EA[istr] = kf;

	}

}

// row_ptr - один из трёх массивов составляющих
// хранение матрицы в Column Row Storage формате.
// два других col_ind - номера столбцов, и val -
// сами ненулевые коэффициенты матрицы.
// row_ptr - это номер начала новой строки в CRS формате
// представления матрицы в котором матрица СЛАУ хранится построчно.
// Функция calculate_row_ptr вычисляет массив начала каждой строки матрицы
// row_ind_SA, а также массив конца каждой строки row_ind_EA матрицы Amat
// хранящейся по строкам в Ak1 типе хранения.
void calculate_row_ptr(integer istart, integer iend,
	integer*& row_ind_SA, integer*& row_ind_EA,
	bool*& flag, Ak1 const *const Amat)
{

	integer i_size_75 = 0;
	// Это нельзя распараллелить.
	for (integer ii = istart; ii <= iend; ++ii) {
		bool cnd = (flag[Amat[ii].i] == false);
		if (!cnd) continue;
		
		row_ind_SA[Amat[ii].i] = ii;
		flag[Amat[ii].i] = true;
		i_size_75++;
		
	}

	integer ilen = (i_size_75 + 1);

#pragma omp parallel for
	for (integer istr = 1; istr < ilen; ++istr)
	{
		integer kf = row_ind_SA[istr];
		while ((kf <= iend) && (Amat[kf].i == istr)) {
			kf++;
		}
		kf--;
		row_ind_EA[istr] = kf;

	}

}// calculate_row_ptr


// row_ptr - один из трёх массивов составляющих
// хранение матрицы в Column Row Storage формате.
// два других col_ind - номера столбцов, и val -
// сами ненулевые коэффициенты матрицы.
// row_ptr - это номер начала новой строки в CRS формате
// представления матрицы в котором матрица СЛАУ хранится построчно.
// Функция calculate_row_ptr_j вычисляет массив начала каждой строки матрицы
// row_ind_SA, а также массив конца каждой строки row_ind_EA транспонированной матрицы Amat
// хранящейся по столбцам в результате операции транспонирования в Ak1 типе хранения.
// Т.е. здесь матрица хранится по столбцам как бы в CCS (Compress Column Storage) формате.
// Транспонирование матрицы осуществляется с помощью экономичного алгоритма сортировки
// за линейное время.
void calculate_row_ptr_j(integer istart, integer iend,
	integer*& row_ind_SA, integer*& row_ind_EA,
	bool*& flag, Ak1 const *const Amat)
{

	
	// Это нельзя распараллелить.
	for (integer ii = istart; ii <= iend; ++ii) {
		bool cnd = (flag[Amat[ii].j] == false);
		if (!cnd) continue;

		
		// сканируем построчно.
		integer istr = Amat[ii].j;
		integer ic = ii;

		integer kf = ic;

		while ((kf <= iend) && (Amat[kf].j == istr)) {
			kf++;
		}
		kf--;
		row_ind_SA[istr] = ic;
		row_ind_EA[istr] = kf;
		flag[Amat[ii].j] = true;
		ii = kf;
		
	}// for

}// calculate_row_ptr_j


// row_ptr - один из трёх массивов составляющих
// хранение матрицы в Column Row Storage формате.
// два других col_ind - номера столбцов, и val -
// сами ненулевые коэффициенты матрицы.
// row_ptr - это номер начала новой строки в CRS формате
// представления матрицы в котором матрица СЛАУ хранится построчно.
// Функция calculate_row_ptr вычисляет массив начала каждой строки матрицы
// row_ind_SA для матрицы Amat
// хранящейся по строкам в Ak2 типе хранения.
// Это быстрая версия функции.
void calculate_row_ptr(integer istart, integer iend,
	integer*& row_ind_SA, 
	bool*& flag, const Ak2& Amat)
{
	
	// Это нельзя распараллелить.
	for (integer ii = istart; ii <= iend; ++ii) {
		bool cnd = (flag[Amat.i[ii]] == false);
		if (!cnd) continue;
		
		row_ind_SA[Amat.i[ii]] = ii;
		flag[Amat.i[ii]] = true;
		
	}// for

	row_ind_SA[Amat.i[iend] + 1] = iend + 1;

}//calculate_row_ptr

// row_ptr - один из трёх массивов составляющих
// хранение матрицы в Column Row Storage формате.
// два других col_ind - номера столбцов, и val -
// сами ненулевые коэффициенты матрицы.
// row_ptr - это номер начала новой строки в CRS формате
// представления матрицы в котором матрица СЛАУ хранится построчно.
// Функция calculate_row_ptr вычисляет массив начала каждой строки матрицы
// row_ind_SA для матрицы Amat
// хранящейся по строкам в Ak2 типе хранения.
// Это быстрая версия функции.
void calculate_row_ptr_cuda(integer istart, integer iend,
	integer*& row_ind_SA,
	bool*& flag, const Ak2& Amat)
{

	// Это нельзя распараллелить.
	for (integer ii = istart; ii <= iend; ++ii) {
		bool cnd = (flag[Amat.i[ii]] == false);
		if (!cnd) continue;

		row_ind_SA[Amat.i[ii]-1] = ii;
		flag[Amat.i[ii]] = true;

	}// for

	row_ind_SA[Amat.i[iend]] = iend + 1;

}//calculate_row_ptr

// row_ptr - один из трёх массивов составляющих
// хранение матрицы в Column Row Storage формате.
// два других col_ind - номера столбцов, и val -
// сами ненулевые коэффициенты матрицы.
// row_ptr - это номер начала новой строки в CRS формате
// представления матрицы в котором матрица СЛАУ хранится построчно.
// Функция calculate_row_ptr вычисляет массив начала каждой строки матрицы
// row_ind_SA для матрицы Amat
// хранящейся по строкам в Ak1 типе хранения.
// Это быстрая версия функции.
void calculate_row_ptr(integer istart, integer iend,
	integer*& row_ind_SA, 
	bool*& flag, Ak1* &Amat)
{
		
	// Это нельзя распараллелить.
	for (integer ii = istart; ii <= iend; ++ii) {
		if (flag[Amat[ii].i] == false) {
			row_ind_SA[Amat[ii].i] = ii;
			flag[Amat[ii].i] = true;
		}
	}
	row_ind_SA[Amat[iend].i + 1] = iend + 1;

} // calculate_row_ptr

// для связи Solution и Setup phase
// Эти переменные общие для обоих фаз
// С помощью этой структуры данных 
// осуществляется передача данных из
// одной фазы amg в другую.
typedef struct Tamg_precond_param {
	integer maxlevel;
	integer ilevel;
	bool bprint_mesage_diagnostic;
	bool bonly_serial;
	integer bILU2smoother;
	integer nnz_P_memo_0;
	integer nnz_P_memo_all;
	doublereal dr_grid_complexity;
	bool debug_reshime;
	doublereal dapply_ilu_max_pattern_size;
	doublereal RealZERO;
	bool identiti;
	MY_AMG_SPLITTING_COARSENING_ALGORITHM  memo_icoarseningtype;
} amg_precond_param;

// Фаза решения в amg методе 14.04.2020.
// Вызывает V цикл amg решателя, устанавливает выбранную
// пользователем во внешнем интерфейсе процедуру сглаживания.
// В зависимости от настроек пользователя вызывает алгоритмы 
// BiCGStab или FGMRes(m)  с V-цикл amg в качестве предобуславливателя.
// V- цикл также может быть вызван самостоятельно.
template <typename doublerealT>
bool solution_phase(Ak2& Amat, // Матрица СЛАУ в CRS формате.
	Ak2& Amat_copy, // Не подвергнутая сортировке матрица СЛАУ в CRS формате.
	integer nsizeA, // количество ячеек выделенное извне для хранилища матриц А	
	integer nnz, // number of non zero elements исходной матрицы нулевого уровня.
	integer n, // dimension of vectors x and b.	
	doublereal*& x, // solution (решение) 
	doublereal*& b, // rthdsd (правая часть).
	real_mix_precision& ret74,
	integer iVar,// Идентификатор физической переменной для которой вызван алгебраический многосеточный метод. Например iVar==TEMP температура. 
	bool bmemory_savings,	
	BLOCK*& my_body, int& lb, // Блоки (типов HOLLOW, SOLID или FLUID) из которых состоит рассчитываемая модель. lb - число блоков.
	integer maxelm_out,
	// Параметры нужные только для solution phase.
	amg_precond_param& amg_pp,
	// исходный уровень (начальный уровень с самой подробной матрицей) имеет номер 0.
	integer* &n_a, // количество неизвестных на каждом из построенных уровней.
	integer*& nnz_a, // количество ненулевых коэффициентов в матрице на каждом из построенных уровней.
	integer*& nnz_aRP, // Количество ненулевых коэффициентов в матрице оператора интерполяции на каждом из построенных уровней.
	integer& ibsp_length,	
	BAD_STRING_PATCHING*& bsp,
	integer i_bsp_LIMIT,
	bool* &flag,	
	bool* &F_false_C_true, // C/F разбиение.
	Ak1* &P, // оператор интерполяции.
	int * &whot_is_block,
	INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14 imyinit
) {


	// Фаза решения: 30.01.2022
	// Здесь описаны некоторые соглашения касательно реализованной фазы решения.
	// Все переданные матрицы СЛАУ на всех уровнях включая нулевой уровень
	// преобразуются в начале функции solution_phase так чтобы на диагонали стоял обратный диагональному элемент,
	// причем позицией диагонали становится первый элемент в каждой строке матрицы СЛАУ.
	// При построении ilu разложения, ilu разложение знает что ему передано обратное 
	// значение в качестве диагонального, ilu возвращает диагональ^(-1) к диагонали в своей копии
	// и строит ilu разложение на основе обычной матрицы в которой первый элемент в каждой строке диагональ.
	// Для сглаживателей Зейделя, Роуч SOR, Damped Jacoby, Рунге-Кутты 3 и 5 порядков и Чебышев используются
	// матрицы где на диагонали стоит обратная диагонали величина в первой позиции каждой строки. 



	//std::mt19937 gen(time(0));
	std::mt19937 gen{ std::random_device()() };
	std::uniform_real_distribution<> urd(0, 1);

	 integer maxlevel = amg_pp.maxlevel;
	 integer ilevel = amg_pp.ilevel;
	 bool bprint_mesage_diagnostic = amg_pp.bprint_mesage_diagnostic;
	 bool bonly_serial = amg_pp.bonly_serial;
	 integer bILU2smoother = amg_pp.bILU2smoother;
	 integer nnz_P_memo_0 = amg_pp.nnz_P_memo_0;
	 integer nnz_P_memo_all = amg_pp.nnz_P_memo_all;
	 doublerealT dr_grid_complexity = static_cast<doublerealT>(amg_pp.dr_grid_complexity);
	 bool debug_reshime = amg_pp.debug_reshime;
	 doublerealT dapply_ilu_max_pattern_size = static_cast<doublerealT>(amg_pp.dapply_ilu_max_pattern_size);
	 doublerealT RealZERO = static_cast<doublerealT>(amg_pp.RealZERO);
	 bool identiti = amg_pp.identiti;
	 MY_AMG_SPLITTING_COARSENING_ALGORITHM  memo_icoarseningtype = amg_pp.memo_icoarseningtype;

	 //std::wcout << "my_amg_manager.bCFJacoby  " << my_amg_manager.bCFJacoby << " \n";
	 //std::wcout << "my_amg_manager.b_spai0  " << my_amg_manager.b_spai0 << " \n";
	 //std::wcout << "my_amg_manager.gold_const  " << my_amg_manager.gold_const << " \n";
	 //std::wcout << "my_amg_manager.b_gmres  " << my_amg_manager.b_gmres << " \n";
	 //std::wcout << "my_amg_manager.iRunge_Kutta_smoother  " << my_amg_manager.iRunge_Kutta_smoother << " \n";
	 //system("pause");


	// 4-5-6 30-31 dec 2016 Поддерживается не более 50 уровней вложенности
	// 5.06.2017 Поддерживается не более 100 уровней вложенности включительно. const integer maxlevel=101;
	// 16.02.2019.
	doublereal** diag = nullptr;
	diag = new doublereal * [maxlevel];
	if (diag == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for diag my_gregat_amg.cpp..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		exit(1);
	}
	if (ilevel == 0) {
		// Если уровней ноль то повидимому произошел сбой построения иерархии уровней,
		// но мы всё равно должны хоть как то отработать.
		std::cout << "ERROR: number of level=0;\n";
		diag[0] = (doublereal*)malloc((n_a[0] + 1) * sizeof(doublereal));
		handle_error<doublereal>(diag[0], "diag[", 0, "]", "classic_aglomerative_amg_6", (n_a[0] + 1));
	}
	else {
		for (integer i_id_level_local = 0; i_id_level_local < maxlevel; ++i_id_level_local) {
			diag[i_id_level_local] = nullptr; // инициализация.
			if (ilevel > i_id_level_local) {
				diag[i_id_level_local] = (doublereal*)malloc((n_a[i_id_level_local] + 1) * sizeof(doublereal));
				handle_error<doublereal>(diag[i_id_level_local], "diag[", i_id_level_local, "]", "classic_aglomerative_amg_6", (n_a[i_id_level_local] + 1));
			}
		}
	}

	// именно тип double, так требует функция seidelq.
	doublereal** diag_minus_one = nullptr;
	diag_minus_one = new doublereal * [maxlevel];
	if (diag_minus_one == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout <<  "Problem: not enough memory on your equipment for diag_minus_one my_gregat_amg.cpp..." << std::endl;
		std::cout <<  "Please any key to exit..." << std::endl;
		exit(1);
	}
	//std::cout << "ilevel=" <<  ilevel << std::endl;
	for (integer i_id_level_local = 0; i_id_level_local < maxlevel; ++i_id_level_local) {
		diag_minus_one[i_id_level_local] = nullptr; // инициализация.
		
		if (ilevel > i_id_level_local) {
			diag_minus_one[i_id_level_local] = (doublereal*)malloc((n_a[i_id_level_local] + 1) * sizeof(doublereal));
			handle_error<doublereal>(diag_minus_one[i_id_level_local], "diag_minus_one[", i_id_level_local, "]", "classic_aglomerative_amg_6", (n_a[i_id_level_local] + 1));
		}
		if ((ilevel!=0)&&(ilevel == i_id_level_local)) {
			diag_minus_one[i_id_level_local] = (doublereal*)malloc((n_a[ilevel-1] + 1) * sizeof(doublereal));
			handle_error<doublereal>(diag_minus_one[i_id_level_local], "diag_minus_one[", i_id_level_local, "]", "classic_aglomerative_amg_6", (n_a[i_id_level_local-1] + 1));

		}
	}

	for (integer i_id_level_local = 0; i_id_level_local < maxlevel; ++i_id_level_local) {
		if (ilevel > i_id_level_local) {
			for (integer i_96 = 1; i_96 <= n_a[i_id_level_local]; ++i_96) {
				diag_minus_one[i_id_level_local][i_96] = 1.0;
			}
		}
		if ((ilevel!=0)&&(ilevel == i_id_level_local)) {
			for (integer i_96 = 1; i_96 <= n_a[i_id_level_local-1]; ++i_96) {
				diag_minus_one[i_id_level_local][i_96] = 1.0;
			}
		}
	}
	// Помечаем сигналы для правой части где её необходимо домножить на минус 1,0.
	for (integer i_96 = 0; i_96 < ibsp_length; ++i_96) {
		diag_minus_one[bsp[i_96].ilevel][bsp[i_96].istring_number] = -1.0;
	}

	bnested_desection_global_amg = nullptr;
	bool** nested_desection = nullptr;
	nested_desection = new bool* [maxlevel];
	if (nested_desection == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout <<  "Problem: not enough memory on your equipment for nested_desection my_gregat_amg.cpp..." << std::endl;
		std::cout <<  "Please any key to exit..." << std::endl;
		exit(1);
	}
	for (integer i_id_level_local = 0; i_id_level_local < maxlevel; ++i_id_level_local) {
		nested_desection[i_id_level_local] = nullptr;
	}

	if (!bonly_serial) {
		// nested desection start
		bnested_desection_global_amg = (bool*)malloc((n_a[0] + 1) * sizeof(bool));
		handle_error<bool>(bnested_desection_global_amg, "bnested_desection_global_amg", "classic_aglomerative_amg_6", (n_a[0] + 1));


		nested_desection[0] = (bool*)malloc((n_a[0] + 1) * sizeof(bool));
		handle_error<bool>(nested_desection[0], "nested_desection[0]", "classic_aglomerative_amg_6", (n_a[0] + 1));

		//maxlevel==201
		for (integer i_17 = 1; i_17 <= maxlevel - 1; ++i_17) {
			if (ilevel > i_17) {
				nested_desection[i_17] = (bool*)malloc((n_a[i_17] + 1) * sizeof(bool));
				handle_error<bool>(nested_desection[i_17], "nested_desection[i_17]", "classic_aglomerative_amg_6", (n_a[i_17] + 1));
			}
		}

	}
	// nested_desection_end


	const integer isize_row_ptr = 4 * n_a[0] + 1;

	integer* row_ptr_start = nullptr;
	//row_ptr_start = new integer[4 * n_a[0] + 1];
	row_ptr_start = (integer*)malloc((isize_row_ptr) * sizeof(integer));
	handle_error<integer>(row_ptr_start, " row_ptr_start", "classic_aglomerative_amg_6", (isize_row_ptr));

	integer* row_ptr_end = nullptr;
	//row_ptr_end = new integer[4 * n_a[0] + 1];
	row_ptr_end = (integer*)malloc((isize_row_ptr) * sizeof(integer));
	handle_error<integer>(row_ptr_end, " row_ptr_end", "classic_aglomerative_amg_6", (isize_row_ptr));

	// ILU2
	LEVEL_ADDITIONAL_DATA* milu2 = nullptr;
	// инициализация.
	init_level_additional_data(milu2, ilevel);

	// ILU0
	LEVEL_ADDITIONAL_DATA0* milu0 = nullptr;
	// инициализация.
	init_level_additional_data(milu0, ilevel);

	// Освобождение общей памяти в ILU буфере.
	if (milu_gl_buffer.alu_copy != nullptr) delete[] milu_gl_buffer.alu_copy;
	if (milu_gl_buffer.jlu_copy != nullptr) delete[] milu_gl_buffer.jlu_copy;
	if (milu_gl_buffer.ju_copy != nullptr) delete[] milu_gl_buffer.ju_copy;
	milu_gl_buffer.alu_copy = nullptr;
	milu_gl_buffer.jlu_copy = nullptr;
	milu_gl_buffer.ju_copy = nullptr;

	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
#pragma omp parallel for
	for (integer i = 1; i <= n; ++i) {
		flag[i] = false;
	}
	for (integer ii = 1; ii <= nnz_a[0]; ++ii) {
		if (flag[Amat.i[ii]] == false) {
			integer istr = Amat.i[ii];
			integer ic = ii;
			integer icdiag = ii;
			if (istr >= isize_row_ptr) {
				std::cout << "need to increase isize_row_ptr " <<  istr << std::endl;

				//system("pause");
				system("PAUSE");
				exit(1);
			}
			row_ptr_start[istr] = ii;
			doublerealT ap = 0.0; // значение на диагонали.
								  //x[istr] = b[istr];
			while ((ic <= nnz_a[0]) && (Amat.i[ic] == istr)) {
				if (Amat.j[ic] != istr) {
					//x[istr] += -Amat.aij[ic]*x[Amat.j[ic]];
					// Все внедиагональные элементы должны быть строго отрицательны.
					// Если это не так то надо выдавать предупреждение о логической ошибке пользователю.
					if (Amat.aij[ic] >= 0.0) {
						//std::cout << "polochitelnji vnediagonalnj element " <<  Amat.aij[ic] << " in matrix level 0 in string " << istr << "..." << std::endl;

						// Вдруг это не страшно 26 октября 2016.
						// Ну да на задача с конвекцией встречается даже и на нулевом уровне вложенности.
						//system("PAUSE");
					}
				}
				else {
					// диагональный элемент строго положителен.
					ap = Amat.aij[ic];
					icdiag = ic;
				}
				ic++;
			}
			if (istr >= isize_row_ptr) {
				std::cout <<  "need to increase isize_row_ptr " << istr << std::endl;
				//system("pause");
				system("PAUSE");
				exit(1);
			}
			row_ptr_end[istr] = ic - 1;
			if (ap < RealZERO) {
				std::cout <<  "zero or negative diagonal elements in string " << istr <<"in basic matrix" << std::endl;
				system("PAUSE");
				exit(1);
			}
			else {
				//x[istr] /= ap;
			}

			flag[Amat.i[ii]] = true;
			swap(Amat, ii, icdiag);

			if (bmemory_savings) {
				// По исходному номеру получаем текущий,
				// но теперь два текущих поменялись.
				the_original_order_of_values[the_original_order_of_values_reverse[ii]] = icdiag;
				the_original_order_of_values[the_original_order_of_values_reverse[icdiag]] = ii;
			}

			diag[0][Amat.i[ii]] = ap; // для ускорения вычисления невязки.
			Amat.aij[ii] = static_cast<real_mix_precision>(1.0 / ap); // умножение быстрей деления.
		}
	}

	/* {
		doublereal s231 = 0.0;
		for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {
			doublereal s232 = b[i_1];
			for (integer i_2 = row_ptr_start[i_1]; i_2 <= row_ptr_end[i_1]; ++i_2) {
				if (Amat.j[i_2] == i_1) {
					s232 -= (1.0 / Amat.aij[i_2]) * x[Amat.j[i_2]];
				}
				else {
					s232 -= Amat.aij[i_2] * x[Amat.j[i_2]];
				}
			}
			s232 = s232 * s232;
			s231 += s232;
		}
		// отладка. Причина использование миксовой точности в задачах на излучение. Нужно использовать только double.
		std::cout <<"aprior 0 "<< sqrt(s231)<< "\n";
		std::cout << "x=" << norma(x, n_a[0]) <<" b=" << norma(b, n_a[0]) << "\n";
		system("pause");
	}*/

	if (bILU2smoother == 2) {
		// ILU2
		if (bprint_mesage_diagnostic) {
			std::cout << "apply ilu" << my_amg_manager.lfil <<" smoother for number 0 level\n";
		}
		equation3DtoCRSRUMBA1(milu2[0], true, Amat, 1, n_a[0], row_ptr_start, row_ptr_end, 0, 0, iVar);
	}

	

	// 14 сентября 2015 понедельник четвёртый уровень вложенности.
	// Уровни вложенности с первого по седьмой сразу. 12.07.2016.

	// Заголовок 29.10.2016.
	if (bprint_mesage_diagnostic) {
		std::cout <<  "1. (number positive connections)/Nnz [%], 2. max positive/ diagonal [average (positive connection/ diagonal)]" << std::endl;
	}

	for (integer ilevel_detector = 1; ilevel_detector <= maxlevel-1; ++ilevel_detector) {

		// Обработка матрицы действует до 99 уровня включительно, но
		// сбор статистики желательно сделать для всех уровней.
		const integer istop_level_scan = maxlevel - 2;

		if (ilevel > ilevel_detector) {

			integer istart_row_ptr0 = 0;
			for (integer ilev = 0; ilev < ilevel_detector; ++ilev) {
				istart_row_ptr0 += n_a[ilev];
			}

			doublerealT inum_vnediagonal_all = 0.0;
			doublerealT inum_only_positive_vnediagonal = 0.0;
			doublerealT memo_diagonal_element = 0.0;
			doublerealT max_positive_connections_element = -1.0;
			doublerealT ratio_positive_connections_by_diagonalelement = 0.0;
			doublerealT ratio_positive_connections_by_diagonalelement_avg = 0.0;
			bool b_ne_menee_2_positive_con_in_string = false;
			doublerealT inum_only_positive_vnediagonal_ne_menee2_in_string = 0.0;

#pragma omp parallel for
			for (integer i = 1; i <= n; ++i) {
				flag[i] = false;
			}
			integer ist = 1;
			for (integer ilev = 0; ilev < ilevel_detector; ++ilev) {
				ist += nnz_a[ilev];
			}
			integer iend = 0;
			for (integer ilev = 0; ilev <= ilevel_detector; ++ilev) {
				iend += nnz_a[ilev];
			}
			integer istPR = 1;
			for (integer ilev = 0; ilev < ilevel_detector; ++ilev) {
				istPR += nnz_aRP[ilev];
			}
			integer iendPR = 0;
			for (integer ilev = 0; ilev <= ilevel_detector; ++ilev) {
				iendPR += nnz_aRP[ilev];
			}

			const doublerealT theta7 = static_cast<doublerealT>(theta(ilevel)); // передаётся в функцию извне.

			double dn_num = 0.0;
			for (integer ii = ist; ii <= iend; ++ii) {
				if (flag[Amat.i[ii]] == false) {

					integer istr = Amat.i[ii];
					integer ic = ii;
					integer icdiag = ii;
					integer istart_row_ptr = istr + istart_row_ptr0;

					max_positive_connections_element = 0.0;
					dn_num += 1.0;


					if (istart_row_ptr >= isize_row_ptr) {
						std::cout <<  "need to increase isize_row_ptr " << istart_row_ptr << std::endl;
						//system("pause");
						system("PAUSE");
						exit(1);
					}
					if (ilevel_detector <= istop_level_scan) {
						row_ptr_start[istart_row_ptr] = ii;
					}
					doublerealT ap = 0.0;
					//doublerealT sum_4 = 0.0;




					
					b_ne_menee_2_positive_con_in_string = false;
					integer inum_pos_con_in_string = 0;
					doublerealT threshold7 = -1.0;
					integer ic7 = ic;
					while ((ic7 <= iend) && (Amat.i[ic7] == istr)) {
						if (Amat.j[ic7] != istr) {
							if (Amat.aij[ic7] >= 0.0) {
								inum_pos_con_in_string++;
								if (fabs(Amat.aij[ic7]) > threshold7) threshold7 = fabs(Amat.aij[ic7]);
							}
						}
						ic7++;
					}
					// мы обнаружили не менее двух положительных связей в данной строке.
					if (inum_pos_con_in_string >= 2) {
						inum_pos_con_in_string = 0;
						ic7 = ic;
						while ((ic7 <= iend) && (Amat.i[ic7] == istr)) {
							if (Amat.j[ic7] != istr) {
								if ((Amat.aij[ic7] >= 0.0) && (fabs(Amat.aij[ic7]) >= theta7 * threshold7)) {
									inum_pos_con_in_string++;
								}
							}
							ic7++;
						}

						if (inum_pos_con_in_string >= 2) {
							b_ne_menee_2_positive_con_in_string = true;
						}
					}


					//x[istr] = b[istr];
					while ((ic <= iend) && (Amat.i[ic] == istr)) {
						if (Amat.j[ic] != istr) {
							//x[istr] += -Amat.aij[ic]*x[Amat.j[ic]];
							inum_vnediagonal_all += 1.0;
							// Все внедиагональные элементы должны быть строго отрицательны.
							// Если это не так то надо выдавать предупреждение о логической ошибке пользователю.
							if (Amat.aij[ic] >= 0.0) {
								//std::cout << "polochitelnji vnediagonalnj element " << Amat.aij[ic] << " in matrix level " << ilevel_detector << " in string " << istr << "..."<< std::endl;

								//system("PAUSE");
								inum_only_positive_vnediagonal += 1.0;

								if (b_ne_menee_2_positive_con_in_string) {
									if (fabs(Amat.aij[ic7]) >= theta7 * threshold7) {
										inum_only_positive_vnediagonal_ne_menee2_in_string += 1.0;
									}
								}

								// Определение величины максимальной внедиагональной связи.
								if (Amat.aij[ic] > 0.0) {
									if (max_positive_connections_element < Amat.aij[ic]) {
										max_positive_connections_element = Amat.aij[ic];
									}
								}
							}
						}
						else {
							ap = Amat.aij[ic];
							memo_diagonal_element = ap;
							icdiag = ic;
						}
						ic++;
					}


					if (istart_row_ptr >= isize_row_ptr) {
						std::cout <<  "need to increase isize_row_ptr "<<  istart_row_ptr << std::endl;
						//system("pause");
						system("PAUSE");
						exit(1);
					}
					if (ilevel_detector <= istop_level_scan) {
						row_ptr_end[istart_row_ptr] = ic - 1;
					}
					if (fabs(ap) < RealZERO) {
						std::cout <<  "zero diagonal elements in string " << istr << " in level "<< ilevel <<" matrix";
						system("PAUSE");
						exit(1);
					}
					else {
						//x[istr] /= ap;
					}


					ratio_positive_connections_by_diagonalelement_avg += fabs(max_positive_connections_element / memo_diagonal_element);
					if (ratio_positive_connections_by_diagonalelement < fabs(max_positive_connections_element / memo_diagonal_element)) {
						ratio_positive_connections_by_diagonalelement = fabs(max_positive_connections_element / memo_diagonal_element);
					}
					flag[Amat.i[ii]] = true;
					if (ilevel_detector <= istop_level_scan) {
						swap(Amat,ii,icdiag);
						diag[ilevel_detector][Amat.i[ii]] = ap;// для ускорения вычисления невязки.						
						Amat.aij[ii] = static_cast<real_mix_precision>(1.0 / ap); // умножение быстрей деления.
					}


				}
			}

			integer iadd_now = 0;
			for (integer i54 = 1; i54 <= ilevel_detector; ++i54) {
				iadd_now += n_a[i54 - 1];
			}
			if (ilevel_detector <= istop_level_scan) {
				if (bILU2smoother == 2) {
					if (bprint_mesage_diagnostic) {
						std::cout <<  "apply ilu" << my_amg_manager.lfil <<" smoother for number "<< ilevel_detector <<" level\n";
					}

					equation3DtoCRSRUMBA1(milu2[ilevel_detector], true,
						Amat, 1, n_a[ilevel_detector], row_ptr_start, row_ptr_end, iadd_now, ilevel_detector, iVar);
				}
			}

			// statistic log:
			if (bprint_mesage_diagnostic) {
				//std::cout << "procent positive connections " << 100.0*inum_only_positive_vnediagonal / inum_vnediagonal_all << " " << std::endl;
				//std::cout << "the ratio of the maximum positive connections to the diagonal" << std::endl;
				//std::cout << "element in the row, in procent " << 100.0*ratio_positive_connections_by_diagonalelement << std::endl;
				//std::cout << "\n";

				
				std::cout << ilevel_detector << " ";
				//std::cout.precision(4);
				std::cout.width(12);
				std::cout << 1.00 * inum_only_positive_vnediagonal / inum_vnediagonal_all;
				std::cout << "  [ ";
				//std::cout.precision(2);
				std::cout.width(12);
				std::cout << 100.0 * inum_only_positive_vnediagonal / inum_vnediagonal_all;
				std::cout << " % ] ";
				std::cout.width(12);
				std::cout << ratio_positive_connections_by_diagonalelement;
				std::cout << "  [";
				//std::cout.precision(3);
				std::cout.width(12);
				std::cout  << ratio_positive_connections_by_diagonalelement_avg / dn_num;
				std::cout << " ]" << std::endl;
			}
		}


	}
	
	if (bILU2smoother > 0) {
		// Пауза только в случае применения ILU декомпозиции.
		//system("PAUSE");
		if (bILU2smoother == 2) {
			// Осторожно возможно код быстро устареет.
			// Выделение оперативной памяти под централизованное хранилище 
			// для ILU.
			memory_allocation_apostoriory_buffer_ilu(milu2, ilevel - 1);
			//memory_allocation_apostoriory_buffer_ilu(milu2, ilevel-1);// 4.01.2017
		}
	}


	if (!bonly_serial) {

		// Готовим nested desection
		// для двух потоков.
		// Самая подробная матрица 1.
		// nested_desection[1]		
		// maxlevel==201
		for (integer i_17 = 0; i_17 <= maxlevel - 1; ++i_17) {
			if (ilevel > i_17) {
				integer inasum = 0;
				for (integer i_18 = 0; i_18 < i_17; ++i_18) inasum += n_a[i_18];
				nested_desection_patch(Amat, n_a[i_17], nested_desection[i_17], row_ptr_start, row_ptr_end, inasum);
				if (bprint_mesage_diagnostic) {

					std::cout <<  "part" << i_17 <<std::endl; 
					std::cout <<  "nested desection is finish" << std::endl;
				}
			}
		}


	}


	if (Amat.i != nullptr) {
		// Освобождаем целую треть памяти для иерархии матриц, т.к. вместо 
		// обращения к индексу i у нас есть row_ptr 
		// row_ptr_start, row_ptr_end (ссылки на начало и конец каждой строки).
		free(Amat.i);
		Amat.i = nullptr;
	}
	//if (R != nullptr) {
		// Используется только оператор P, оператор R точно такой же что и P с точностью до сортировки.
		// Методы restriction и prolongation не чувствительны к сортировке.
		//free(R);
		//R = nullptr;
	//}

	if (flag != nullptr) {
		free(flag);
		flag = nullptr;
	}

	if (bprint_mesage_diagnostic) {
		std::cout <<  "cycling: V cycle." << std::endl;
		std::cout << "level=" << ilevel << std::endl;

		std::cout << "multigrid R.P.Fedorenko 1961." << std::endl;
		std::cout << "standart aglomerative algebraic multigrid method." << std::endl;
	}
	if (debug_reshime) system("pause");
	//exit(1);

	// 10 11 21 multigrid tutorial Вильм Бригг.
	// Высокорейнольдсовое обтекание квадрата в DavisTest,
	// решатель работал на x-компоненте скорости. Сетка сгущалась
	// к поверхности квадрата достаточно сильно. Это дало расходимость
	// amg v0.08 решателя на данной задаче с параметрами nu1=4, nu2=3.
	// Параметры nu1=8, nu2=7 обеспечили сходимость вычислительного процесса.

	// Возможно имеет смысл сделать управляемый выход из сглаживателя, допустим
	// если невязка опустилась ниже первоначальной в 0.1 раз то имеет смысл досрочно 
	// оборвать итерации сглаживателя.


	// nu1=0; nu2=2; nFinestSweeps=2 is recomended 
	// Masashi Imano. Optimization of parameter setting for GAMG
	// solver in simple solver. 
	// Aug 26th 2012. OpeanFOAM study. 
	// nu1=0 имеем расходимость на BSK Dmitrii.
	// nFinestSweeps=2 имеем расходимость на BSK Dmitrii.
	// BSK Dmitrii сходится при nu1=1, nu2=2, nFinestSweeps=3.

	integer nu1 = 1; // minimum value 1 // 4 // 8
	integer nu2 = 2; // minimum value 2 // 3 // 7

					 // на задаче Finned Heat Sync из первого туториала Icepak была обнаружена расходимость 
					 // для Y скорости и поправки давления. При этом обтекание куба отлично считалось на равномерной
					 // сетки с nu1=1, nu2=2 даже при весьма больших числах Рейнольдса.
					 // при nu1=10, nu2=10 скорости разрешаются хорошо и проблем с ними нет, но поправка давления по прежнему даёт сбой.
					 // при nu==20 сбой всё равно есть.
					 // не помогло.
					 //nu1 = 40;
					 //nu2 = 40;

	integer nFinestSweeps = 2;


	// с 26 октября 2016 мы передаём настройки из интерфейса AliceMesh_v0_39.
	// Т.к. есть трудносходящиеся задачи, то эти настройки должны помочь.
	nu1 = my_amg_manager.nu1;
	nu2 = my_amg_manager.nu2;
	nFinestSweeps = my_amg_manager.nFinnest;

	//if (iVar == PAM) {
	//nFinestSweeps = 300;
	//nu1 = 0;
	//nu2 = 20;
	//}
	// для Finner Heat Sink надо усилить сглаживания.
	// Это не помогает будет перенаправление на другой алгоритм.
	//if (iVar == PAM) {
	//nu1 = 3;
	//nu2 = 3;
	//nFinestSweeps = 6;
	//}
	const bool btheoryGuideANSYSFluent = false;
	if (iVar != PAM) {
		if (btheoryGuideANSYSFluent) {
			// Так написано в Theory Guide ANSYS Fluent.
			nu1 = 0;
			nu2 = 1;
			nFinestSweeps = 1;
		}
	}



	// Двойной вакуумный промежуток вызывает проблемы сходимости:
	//nu1 = 10;
	//nu2 = 20;

	// Смысл этих параметров в том что они экономят ресурсы процессора
	// в теории осуществляя досрочный выход из пред и пост сглаживаний.
	// Т.е. параметры nu1,nu2 задаются с запасом и алгоритм сам использует
	// сколько ему взять итераций для оптимальной работы (сходимости). 
	// Пользователь не ломает голову какие задавать параметры nu1, nu2 
	// а задаёт их верхние предельные значения. 
	doublerealT process_flow_beta = static_cast <doublerealT>(0.7);
	doublerealT process_flow_alpha = static_cast <doublerealT>(0.1);
	bool process_flow_logic = false;

	if (process_flow_logic) {
		nu1 = 40;
		nu2 = 40;
	}

	d_my_optimetric2_6_12_2019 = dr_grid_complexity;
	if (bprint_mesage_diagnostic) {
		std::cout << "operator complexity is " << dr_grid_complexity << std::endl;
		doublerealT gridCG = 1.0;
		for (int i93 = 1; i93 <= ilevel; ++i93) gridCG += 1.0 * n_a[i93] / (1.0 * n_a[0]);
		std::cout << "grid complexity is Cg = " <<  gridCG << std::endl;
		std::cout << "Prolongation operator complexity is |Psigma|/|P1|=" << static_cast<doublerealT>(1.0 * nnz_P_memo_all / nnz_P_memo_0) <<"  "<< static_cast<doublerealT>(1.0 * nnz_P_memo_all / n_a[0]) <<"*n" << std::endl;
		std::cout << "number pre-cycles nu1=" << nu1 <<", number post-cycles nu2=" << nu2 <<"\n";
	}

	//ilevel = 1; //debug
	doublerealT rho = 1.0;
	doublerealT dres = 1.0;
	integer iiter = 1;
	//const doublerealT tolerance = 1.0e-12;
	// 13 февраля 2016 калибруем точность солвера с целью ускорения получения результата.
	// Т.к. нам ненужна точность выше чем десятая доля градуса по температуре.
	// начальное значение невязки составляет примерно 7000.0.
	doublerealT tolerance = static_cast <doublerealT>(0.0001); // точность выхода по классическому определению L2 нормы.
	
	// 23 октября 2016
	if (bSIMPLErun_now_for_temperature) {
		// Решаем cfd задачи.
		tolerance = static_cast <doublerealT>(1.0e-8);
	}

	doublereal** residual_fine = nullptr;
	residual_fine = new doublereal * [maxlevel];
	if (residual_fine == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for residual_fine my_gregat_amg6.cpp..." << std::endl;
		std::cout <<  "Please any key to exit..." << std::endl;
		exit(1);
	}
	doublereal** residual_coarse = nullptr;
	residual_coarse = new doublereal * [maxlevel];
	if (residual_coarse == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout <<  "Problem: not enough memory on your equipment for residual_coarse my_gregat_amg6.cpp..."  << std::endl;
		std::cout <<  "Please any key to exit..." << std::endl;
		exit(1);
	}
	doublereal** error_approx_coarse = nullptr;
	error_approx_coarse = new doublereal * [maxlevel];
	if (error_approx_coarse == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout <<  "Problem: not enough memory on your equipment for error_approx_coarse my_gregat_amg6.cpp..."  << std::endl;
		std::cout <<  "Please any key to exit..."  << std::endl;
		exit(1);
	}
	doublereal** error_approx_fine = nullptr;
	error_approx_fine = new doublereal * [maxlevel];
	if (error_approx_fine == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout <<  "Problem: not enough memory on your equipment for error_approx_fine my_gregat_amg6.cpp..." << std::endl;
		std::cout <<  "Please any key to exit..." << std::endl;
		exit(1);
	}
	for (integer i_id_level_local = 0; i_id_level_local < maxlevel; ++i_id_level_local) {
		residual_fine[i_id_level_local] = nullptr;
		residual_coarse[i_id_level_local] = nullptr;
		error_approx_coarse[i_id_level_local] = nullptr;
		error_approx_fine[i_id_level_local] = nullptr;
	}


	// Устаревший код инициализации значением nullptr 4 декабря 2016. 

	// Закомментированный код безнадёжно устарел. В данный момент 
	//5.06.2017 поддерживается 100 уровней вложенности.

	// 25.04.2018 На этом месте удалён большой фрагмент устаревшего кода.

	// лучше выделять оперативную память небольшими блоками т.к.
	// оперативная память фрагментирована системными dll и
	// большого свободного блока может не найтись.




	// maxlevel==201
	for (integer i_17 = 1; i_17 <= maxlevel - 1; ++i_17)
	{
		// 05.06.2017
		integer i_17_prev = i_17 - 1;

		if (ilevel + 2 > i_17) {
			//std::cout << "i_17_prev=" << i_17_prev << std::endl;

			// residual
			residual_fine[i_17_prev] = new doublereal[n_a[i_17_prev] + 1];
			for (integer i_s1 = 0; i_s1 <= n_a[i_17_prev]; ++i_s1) {
				residual_fine[i_17_prev][i_s1] = 0.0;
			}
			//residual_fine[i_17_prev] = (doublereal*)malloc((n_a[i_17_prev] + 1) * sizeof(doublereal));
			//handle_error<doublereal>(residual_fine[i_17_prev], "residual_fine[", i_17_prev, "]", "classic_aglomerative_amg_6", (n_a[i_17_prev] + 1));

			if (ilevel + 1 > i_17) {
				residual_coarse[i_17_prev] = new doublereal[n_a[i_17] + 1];
				for (integer i_s1 = 0; i_s1 <= n_a[i_17]; ++i_s1) {
					residual_coarse[i_17_prev][i_s1] = 0.0;
				}
				//residual_coarse[i_17_prev] = (doublereal*)malloc((n_a[i_17] + 1) * sizeof(doublereal));
				//handle_error<doublereal>(residual_coarse[i_17_prev], "residual_coarse[", i_17_prev, "]", "classic_aglomerative_amg_6", (n_a[i_17] + 1));

				error_approx_coarse[i_17_prev] = new doublereal[n_a[i_17] + 1];
				for (integer i_s1 = 0; i_s1 <= n_a[i_17]; ++i_s1) {
					error_approx_coarse[i_17_prev][i_s1] = 0.0;
				}
				//error_approx_coarse[i_17_prev] = (doublereal*)malloc((n_a[i_17] + 1) * sizeof(doublereal));
				//handle_error<doublereal>(error_approx_coarse[i_17_prev], "error_approx_coarse[", i_17_prev, "]", "classic_aglomerative_amg_6", (n_a[i_17] + 1));
			}
			else {
				residual_coarse[i_17_prev] = nullptr;
				error_approx_coarse[i_17_prev] = nullptr;
			}


			error_approx_fine[i_17_prev] = new doublereal[n_a[i_17_prev] + 1];
			for (integer i_s1 = 0; i_s1 <= n_a[i_17_prev]; ++i_s1) {
				error_approx_fine[i_17_prev][i_s1] = 0.0;
			}
			//error_approx_fine[i_17_prev] = (doublereal*)malloc((n_a[i_17_prev] + 1) * sizeof(doublereal));
			//handle_error<doublereal>(error_approx_fine[i_17_prev], "error_approx_fine[", i_17_prev, "]", "classic_aglomerative_amg_6", (n_a[i_17_prev] + 1));
		}
	}




	// 2 jan 2016. 
	integer igam = 1; // 1-V; 2-W
	//INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14 imyinit = ZERO_INIT; // ZERO_INIT optimum

	doublerealT* x_copy = nullptr;
	x_copy = (doublerealT*)malloc((n_a[0] + 1) * sizeof(doublerealT));
	handle_error<doublerealT>(x_copy, "x_copy", "classic_aglomerative_amg_6", (n_a[0] + 1));

	// для ускорения счёта в вакуумном промежутке.
	doublerealT* x_old = nullptr;
	x_old = (doublerealT*)malloc((n_a[0] + 1) * sizeof(doublerealT));
	handle_error<doublerealT>(x_old, "x_old", "classic_aglomerative_amg_6", (n_a[0] + 1));

#pragma omp parallel for
	for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
		x_copy[i47] = static_cast<doublerealT>(x[i47]);
		x_old[i47] = static_cast<doublerealT>(x[i47]);
		//x_copy[i47] = 0.0; // 28.07.2016
	}

	doublerealT* x_best_search = nullptr;
	x_best_search = (doublerealT*)malloc((n_a[0] + 1) * sizeof(doublerealT));
	handle_error<doublerealT>(x_best_search, "x_best_search", "classic_aglomerative_amg_6", (n_a[0] + 1));

	doublerealT res_best_search = static_cast <doublerealT>(1e37);
#pragma omp parallel for
	for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
		x_best_search[i47] = static_cast<doublerealT>(x[i47]);
		//x_best_search[i47] = 0.0; // 28.07.2016
	}


	// Для поправки давления возникает задача когда на всех границах стоит условие Неймана,
	// это приводит к тому что метод работает бесконечно долго и не может сойтись, поэтому нужно 
	// заложить критерий останова по превышению количества допустимых итераций (не более 1000 итераций).
	// 1000 итераций это очень долго поэтому для поправки давления надо подобрать разумное количество 
	// итераций т.к. от этого существенным образом зависит быстродействие гидродинамического алгоритма.
	integer iter_limit = 0;
	integer istop_porog_reconst = 5000;// 50

	bool ret_value = false;
	doublerealT dres_previos = static_cast <doublerealT>(1.0e37);
	integer icount_bad_convergence_Vcycles = 0;
	integer i_count_stagnation = 0;
	doublerealT res0start = static_cast <doublerealT>(1.0e-36);
	bool bfirst_divergence = true;

	residualq2(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
	doublerealT dres_initial = static_cast<doublerealT>(norma(residual_fine[0], n_a[0]));
	if (((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) && (dres_initial > 20.0)) {
		// Это признак ошибки в сборке матрицы СЛАУ на компоненты скорости.
		std::cout <<  "may be problem convergence Speed Flow: very big dres0=" << dres_initial << std::endl;
		std::cout <<  "run residualq2 analysys." << std::endl;
		residualq2_analysys(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
	}
	if (((iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)) && (dres_initial > 20.0)) {
		// Это признак ошибки в сборке матрицы СЛАУ на турбулентные характеристики.
		std::cout <<  "may be problem convergence Turbulence equations: very big dres0=" << dres_initial << std::endl;
		std::cout <<  "run residualq2 analysys." << std::endl;
		residualq2_analysys(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
	}
	if ((iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) && (dres_initial > 6.0e6)) {
		// Это признак ошибки в сборке матрицы СЛАУ на турбулентные характеристики - удельную скорость диссипации.
		std::cout <<  "may be problem convergence Turbulence equations: very big dres0=" << dres_initial << std::endl;
		std::cout <<  "run residualq2 analysys." << std::endl;
		residualq2_analysys(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
	}
	if (((iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) && (dres_initial > 20.0)) {
		// 22.01.2021
		// Это признак ошибки в сборке матрицы СЛАУ на турбулентные характеристики.
		std::cout << "may be problem convergence Turbulence equations Langtry Menter model: very big dres0=" << dres_initial << std::endl;
		std::cout << "run residualq2 analysys." << std::endl;
		residualq2_analysys(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
	}
	/*
	// код заимствованный из amg5:
	integer iflag_cont = 1;
	if (iVar != PAM) {
	dres = fabs(dres_initial);

	if (iVar != TEMP) {
	if (dres < dterminatedTResudual) {
	// Вектор и так точно удовлетворяет решению, его не надо уточнять из решения СЛАУ.
	iflag_cont = 0;
	}
	}
	else {
	if (dres < 1.0e-4*dterminatedTResudual) {
	// Вектор и так точно удовлетворяет решению, его не надо уточнять из решения СЛАУ.
	iflag_cont = 0;
	}
	}
	}
	iflag_cont = 1;
	*/

	if (bprint_mesage_diagnostic) {
		// start residual.
		std::cout <<  0 << " " << dres_initial <<" rho=" << dres_initial / rho << std::endl;
	}

	// 25 10 2016
	integer iflag_cont = 1;
	if (iVar != PAM) {
		dres = fabs(dres_initial);
	}


	integer count_iter_for_film_coef = 0;
	// Если число расходимостей превысит оговорённую константу то произойдёт выход из алгоритма.
	integer i_signal_break_pam_opening = 0;
	// x хорошее значение.
	const integer i_limit_signal_pam_break_opening = 1000; // 8
	//doublerealT delta_old_iter = 1.0e10;



	//if (iVar == PAM) {// бред
	//for (integer iter = 0; iter < 2; ++iter) {
	//seidelq(Amat, 1, n_a[0], b, x, row_ptr_start, row_ptr_end, 0);
	//}
	//}
	integer icount_V_cycle = 0;

	//doublerealT dres_initial_ = 1e-6;


	doublerealT maxold = static_cast <doublerealT>(-1.0e30);
	for (integer i = 1; i <= n_a[0]; ++i) {
		if (x[i] > maxold) maxold = static_cast<doublerealT>(x[i]);
	}

	// с 26 октября 2016 мы передаём настройки из интерфейса AliceMesh_v0_39.
	// Т.к. есть трудносходящиеся задачи, то эти настройки должны помочь.
	// Отсекаем уровни которые выше порогового значения указанного пользователем.
	//if (ilevel > my_amg_manager.maximum_levels) {
	//ilevel = my_amg_manager.maximum_levels;
	//}
	ilevel -= my_amg_manager.maximum_delete_levels;



	doublereal* x_best_search2 = nullptr;
	x_best_search2 = new doublereal[n_a[0] + 1];
	if (x_best_search2 == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for x_best_search2 my_agregat_amg.cpp..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		exit(1);
	}
	doublereal* x_best_search_init = nullptr;
	x_best_search_init = new doublereal[n_a[0] + 1];
	if (x_best_search_init == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for x_best_search_init my_agregat_amg.cpp..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		exit(1);
	}
#pragma omp parallel for
	for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
		x_best_search_init[i47] = x[i47];
		x_best_search2[i47] = x[i47];
	}

	dbgmres_smoother = new DATA_BASE_GMRES[maxlevel];
	for (int i_52 = 0; i_52 < maxlevel; ++i_52) {
		dbgmres_smoother[i_52].val = nullptr;
		dbgmres_smoother[i_52].col_ind = nullptr;
		dbgmres_smoother[i_52].row_ptr = nullptr;
	}

	

	//std::cout << "my_amg_manager.istabilization=" << my_amg_manager.istabilization << std::endl;
	//system("pause");

	if ((my_amg_manager.istabilization == 0) || ((iVar == TEMP) && (my_amg_manager.istabilization == 3))) {


		// ((iVar==TEMP)&&(my_amg_manager.istabilization == 3)) - нелинейное граничное условие в уравнении теплопередачи.

		// Только алгебраический многосеточный метод.

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) tolerance *= static_cast <doublerealT>(1e-11);
		if ((iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
			(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
			(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER)||
			(iVar == RE_THETA_LANGTRY_MENTER)) tolerance *= static_cast <doublerealT>(1e-11);
		if (iVar == PAM) tolerance *= static_cast <doublerealT>(1e-14);
		if (iVar == TEMP) tolerance *= static_cast <doublerealT>(1e-6);
		if (iVar == TOTALDEFORMATIONVAR) tolerance = static_cast <doublerealT>(1.0e-17);
		doublereal minx_gl = 1.0e36;
		doublereal maxx_gl = -1.0e36;

		integer istop_speed_cycling = 10;


		//for (integer iprohod = 0; iprohod < 20; ++iprohod) {
		//while ((iflag_cont == 1) && ((dres>tolerance) || ((iVar != TEMP) && (icount_V_cycle<5)))) {
		///	while ((iflag_cont == 1) && ((dres>tolerance) )) {
		while (((iflag_cont == 1) && ((dres > tolerance))) ||
			((iVar == TEMP) && bSIMPLErun_now_for_temperature && (icount_V_cycle < 9))
			|| ((iVar == TOTALDEFORMATIONVAR) && (icount_V_cycle < 9))) {

			// Обеспечивает колосальное быстродействие без потери сходимости.

			if (bSIMPLErun_now_for_temperature) {
				// гидродинамика.

				//  Этот код непонятен, надо тестировать.
				if (icount_V_cycle > istop_speed_cycling) {



					if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
						(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
						(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) ||
						(iVar == GAMMA_LANGTRY_MENTER) ||
						(iVar == RE_THETA_LANGTRY_MENTER))
					{
						if (dres < 1.0e-3 * dres_initial) {
							break;
						}
						else {
							istop_speed_cycling += 2;
						}
					}
					if (iVar == PAM) {
						if (dres < 1.0e-1 * dres_initial) {
							//break;
						}
						else {
							istop_speed_cycling += 2;
						}
					}
				}
			}




			if (fabs(dres / rho) < 1.0) {
#pragma omp parallel for
				for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
					x_best_search2[i47] = x[i47];
				}
			}

			if (bPhysics_stop) {
				if (icount_V_cycle > 0) {
					doublerealT maxnew = static_cast <doublerealT>(-1.0e30);
					for (integer i = 1; i <= n_a[0]; ++i) {
						if (x[i] > maxnew) maxnew = static_cast<doublerealT>(x[i]);
					}
					if (iVar == TOTALDEFORMATIONVAR) {
						if ((fabs(dres) < 1.0e-2) && (fabs(maxnew - maxold) < 1.0e-9)) {
							std::cout <<  "break bPhysics_stop, dres<1e-2 && (fabs(maxnew - maxold) < 1.0e-9)" << std::endl;
							break;
						}
						else {
							maxold = maxnew;
						}
					}
					else {
						if ((fabs(dres) < 1.0e-2) && (fabs(maxnew - maxold) < 0.0005)) {
							std::cout <<  "break bPhysics_stop, dres<1e-2 && (fabs(maxnew - maxold) < 0.0005)" << std::endl;
							break;
						}
						else {
							maxold = maxnew;
						}
					}

				}
			}


			if (icount_V_cycle > 0) {
				// установить 0 в случае отката на предыдущую стабильную локально-линейную версию алгоритма.
				// главная причина установки значения 1 является сокращение числа проходов для устранения
				// нелинейности в системе с 26 до 4. При установке 1 в данном месте кода надо в модуле
				// mysolver_v0_03 установить fHORF=1.0; 
				if ((iVar == TEMP) && (my_amg_manager.istabilization == 3)) {
					if (bonly_solid_calculation) {
						if (bvacuumPrism) {
							// предполагается неизменный порядок следования позиций в x
							// и rthdsd.

							doublereal* x_temper = nullptr;
							//x_temper = new doublerealT[n_a[0] + 1];
							x_temper = (doublereal*)malloc((static_cast<integer>(n_a[0]) + 1) * sizeof(doublereal));
							if (x_temper == nullptr) {
								// недостаточно памяти на данном оборудовании.
								std::cout <<  "Problem: not enough memory on your equipment for x_temper my_agregat_amg.cpp..." << std::endl;
								std::cout <<  "Please any key to exit..." << std::endl;
								exit(1);
							}
#pragma omp parallel for
							for (integer i23 = 0; i23 < n_a[0]; ++i23) {
								if (x[i23 + 1] < -272.15) x[i23 + 1] = -272.15;
								//x_temper[i23] = x[i23 + 1];
								// 0.01 параметр нижней релаксации.
								// 0.25; 0.2; 0.01.
								// 0.005
								// etalon 0.01 (1250it; 891it; 2it; 7s 770ms)
								// experimental 0.02 (1250it; 28it; 2it; 5s 120ms)
								// experimental 0.04 (740it; 2it; 3s 360ms)
								// experimental 0.05 (618it; 2it; 3s 50ms)
								// experimental 0.06 (533it; 2it; 2s 900ms)
								// experimental 0.07 (469it; 2it; 2s 800ms)
								// experimental 0.08 (420it; 2it; 2s 570ms)
								// experimental 0.09 (381 it; 2it; 2s 390ms)
								// experimental 0.1 (1250it и не сходится, переборщил).
								if (fabs(x[i23 + 1] - x_old[i23 + 1]) > 5.0) {
									// Порог 5С оптимум.
									// 0.04 14s 440ms 358it
									// 0.09 13s 290ms 314it optimum
									x_temper[i23] = x_old[i23 + 1] + 0.084 * (x[i23 + 1] - x_old[i23 + 1]);
								}
								else if (fabs(x[i23 + 1] - x_old[i23 + 1]) > 1.0) {
									// Порог 1С оптимум.
									// 0.09 13s 290ms 314it
									// 0.095 17s 10ms 318; 102; 22;
									// 0.091 316it; 10; 10; 13s 350ms;
									// 0.092 optimum
									x_temper[i23] = x_old[i23 + 1] + 0.092 * (x[i23 + 1] - x_old[i23 + 1]);
								}
								else {
									// 0.09 0.25 перебор.
									// 0,09; 0.11; 314it; 9it; 9it; 2s 310ms optimum
									// 0.09; 0.12; 300it; 33it; 20it; 2s 530ms; 
									// 0.09; 0.1; 308it; 35it; 2s 400ms;
									// 0.11 optimum
									x_temper[i23] = x_old[i23 + 1] + 0.11 * (x[i23 + 1] - x_old[i23 + 1]);
								}
								if (x_temper[i23] < -272.15) x_temper[i23] = -272.15;
								x[i23 + 1] = x_temper[i23];
							}

							// На старте мы блокируем Стефана Больцмана дав сойтись лучистым потокам.
							// Вычисление осреднённых температур в К на границах вакуумных промежутков:
							for (integer i23 = 0; i23 < lb; ++i23) {
								update_avg_temperatures(x_temper, my_body[i23]);
							}
							// Вычисление плотностей радиационных тепловых потоков:
							for (integer i23 = 0; i23 < lb; ++i23) {
								calculation_density_radiation_heat_flux(my_body[i23]);
							}


							doublereal* rthdsd_loc123 = nullptr;
							//rthdsd_loc123 = new doublerealT[n_a[0] + 1];
							rthdsd_loc123 = (doublereal*)malloc((static_cast<integer>(n_a[0]) + 1) * sizeof(doublereal));
							if (rthdsd_loc123 == nullptr) {
								// недостаточно памяти на данном оборудовании.
								std::cout <<  "Problem: not enough memory on your equipment for rthdsd_loc123 my_agregat_amg.cpp..." << std::endl;
								std::cout <<  "Please any key to exit..." << std::endl;
								exit(1);
							}

							for (integer i23 = 0; i23 < n_a[0]; ++i23) {
								rthdsd_loc123[i23] = rthdsd_no_radiosity_patch[i23];
								if ((i23 >= iadd_qnbc_maxelm) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
									(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on)) {
									// 0.25 13s 60ms 314; 10; 13; old value
									// 0.26 13s 0ms; 310; 21; 13;
									// 0.27 318; 13; 3;
									// 0.28 325; 42; 13;
									// 0.35 13s 560ms 347; //new optimum
									// 0.09 1250; 1250; время запредельно большое.
									doublerealT alpha_relax142 = static_cast <doublerealT>(0.35);// d_my_optimetric1_6_12_2019;// 0.35;
									rthdsd_loc123[i23] = alpha_relax142 *
										(-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST *
											((273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) *
												(273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) -
												(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
												(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
												(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
												(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb))) +
										(1.0 - alpha_relax142) *
										(-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST *
											((273.15 + x_old[i23 + 1]) * (273.15 + x_old[i23 + 1]) *
												(273.15 + x_old[i23 + 1]) * (273.15 + x_old[i23 + 1]) -
												(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
												(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
												(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) *
												(273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb)));
									rthdsd_loc123[i23] *= qnbc[i23 - iadd_qnbc_maxelm].dS;
								}
							}

							radiosity_patch_for_vacuum_Prism_Object_(rthdsd_loc123, my_body, lb, maxelm_out, whot_is_block);
#pragma omp parallel for
							for (integer i23 = 0; i23 < n_a[0]; ++i23) {
								x_old[i23 + 1] = static_cast<doublerealT>(x_temper[i23]);
								//x_old[i23 + 1] = x[i23 + 1];
							}
#pragma omp parallel for
							for (integer i23 = 0; i23 < n_a[0]; ++i23) {
								b[i23 + 1] = rthdsd_loc123[i23];
							}

							if (rthdsd_loc123 != nullptr) {
								free(rthdsd_loc123);
							}
							rthdsd_loc123 = nullptr;

							if (x_temper != nullptr) {
								free(x_temper);
							}
							x_temper = nullptr;
						}
						else if (b_sign_on_nonlinear_bc) {
							//  25 декабря 2015. Ускорение сходимости при использовании 
							// нелинейных граничных условий.


							bool bNewtonRichman = false;
							bool bStefanBolcman = false;

							for (integer i23 = 0; i23 < n_a[0]; ++i23) {
								if (i23 >= iadd_qnbc_maxelm) {
									if ((qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
										(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on) &&
										(qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on == false)) {
										// Стефан-Больцман.
										bStefanBolcman = true;
									}
									if ((qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
										(qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on) &&
										(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on == false)) {

										// Ньютон-Рихман.
										bNewtonRichman = true;
									}
									if ((qnbc[i23 - iadd_qnbc_maxelm].bactive) &&
										(qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on) &&
										(qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on)) {
										// Условие смешанного типа.
										bStefanBolcman = true;
										bNewtonRichman = true;
									}
								}

							}


							doublerealT* x_temper = nullptr;
							//x_temper = new doublerealT[n_a[0] + 1];
							x_temper = (doublerealT*)malloc((static_cast<integer>(n_a[0]) + 1) * sizeof(doublerealT));
							if (x_temper == nullptr) {
								// недостаточно памяти на данном оборудовании.
								std::cout <<  "Problem: not enough memory on your equipment for x_temper my_agregat_amg.cpp..." << std::endl;
								std::cout <<  "Please any key to exit..."<< std::endl;
								exit(1);
							}
#pragma omp parallel for
							for (integer i23 = 0; i23 < n_a[0]; ++i23) {
								if (x[i23 + 1] < -272.15) x[i23 + 1] = -272.15;
								// x_temper[i23] = x[i23 + 1];
								// 0.01 параметр нижней релаксации.
								// 0.25
								// 0.2
								// 10 июня 2018 года заменил на коэффициент нижней релаксации равный 0.9.
								// 0.01
								// 0.005

								doublereal alpha_relax_nonlinear_boundary_condition = 0.01;
								if (bNewtonRichman && (bStefanBolcman == false)) {
									alpha_relax_nonlinear_boundary_condition = 0.9;
								}

								x_temper[i23] = static_cast<doublerealT>(x_old[i23 + 1] + alpha_relax_nonlinear_boundary_condition * (x[i23 + 1] - x_old[i23 + 1]));
								if (x_temper[i23] < -272.15) x_temper[i23] = static_cast <doublerealT>(-272.15);
								x[i23 + 1] = x_temper[i23];
							}

							doublerealT* rthdsd_loc123 = nullptr;
							//rthdsd_loc123 = new doublerealT[n_a[0] + 1];
							rthdsd_loc123 = (doublerealT*)malloc((static_cast<integer>(n_a[0]) + 1) * sizeof(doublerealT));
							if (rthdsd_loc123 == nullptr) {
								// недостаточно памяти на данном оборудовании.
								std::cout <<  "Problem: not enough memory on your equipment for rthdsd_loc123 my_agregat_amg.cpp..." << std::endl;
								std::cout <<  "Please any key to exit..." << std::endl;
								exit(1);
							}

							for (integer i23 = 0; i23 < n_a[0]; ++i23) {
								rthdsd_loc123[i23] = static_cast<doublerealT>(rthdsd_no_radiosity_patch[i23]);
								if ((i23 >= iadd_qnbc_maxelm) && (qnbc[i23 - iadd_qnbc_maxelm].bactive) && (qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on) && (qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on == false)) {
									// Стефан Больцман.
									doublerealT alpha_relax142 = 0.25;
									rthdsd_loc123[i23] = static_cast<doublerealT>(alpha_relax142 * (-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST * ((273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) - (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb))) +
										(1.0 - alpha_relax142) * (-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST * ((273.15 + x_old[i23 + 1]) * (273.15 + x_old[i23 + 1]) * (273.15 + x_old[i23 + 1]) * (273.15 + x_old[i23 + 1]) - (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb))));
									rthdsd_loc123[i23] *= static_cast<doublerealT>(qnbc[i23 - iadd_qnbc_maxelm].dS);
								}
								if ((i23 >= iadd_qnbc_maxelm) && (qnbc[i23 - iadd_qnbc_maxelm].bactive) && (qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on) && (qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on == false)) {
									// Ньютон-Рихман.
									//doublerealT alpha_relax142 = 0.25;
									rthdsd_loc123[i23] = static_cast<doublerealT>(-qnbc[i23 - iadd_qnbc_maxelm].film_coefficient * (x_temper[i23] - qnbc[i23 - iadd_qnbc_maxelm].Tamb));
									rthdsd_loc123[i23] *= static_cast<doublerealT>(qnbc[i23 - iadd_qnbc_maxelm].dS);
								}
								if ((i23 >= iadd_qnbc_maxelm) && (qnbc[i23 - iadd_qnbc_maxelm].bactive) && (qnbc[i23 - iadd_qnbc_maxelm].bNewtonRichman_q_on) && (qnbc[i23 - iadd_qnbc_maxelm].bStefanBolcman_q_on)) {
									// Условие смешанного типа.
									//doublerealT alpha_relax142 = 0.25;
									rthdsd_loc123[i23] = static_cast<doublerealT>((-qnbc[i23 - iadd_qnbc_maxelm].emissivity * STEFAN_BOLCMAN_CONST * ((273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) * (273.15 + x_temper[i23]) - (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb) * (273.15 + qnbc[i23 - iadd_qnbc_maxelm].Tamb))));
									rthdsd_loc123[i23] += static_cast<doublerealT>(-qnbc[i23 - iadd_qnbc_maxelm].film_coefficient * (x_temper[i23] - qnbc[i23 - iadd_qnbc_maxelm].Tamb));
									rthdsd_loc123[i23] *= static_cast<doublerealT>(qnbc[i23 - iadd_qnbc_maxelm].dS);
								}
							}
#pragma omp parallel for
							for (integer i23 = 0; i23 < n_a[0]; ++i23) {
								x_old[i23 + 1] = x_temper[i23];
								//x_old[i23 + 1] = x[i23 + 1];
							}
#pragma omp parallel for
							for (integer i23 = 0; i23 < n_a[0]; ++i23) {
								b[i23 + 1] = rthdsd_loc123[i23];
							}

							if (rthdsd_loc123 != nullptr) {
								free(rthdsd_loc123);
							}
							rthdsd_loc123 = nullptr;

							if (x_temper != nullptr) {
								free(x_temper);
							}
							x_temper = nullptr;

						}
					}
				}
			}

			//system("pause");
			icount_V_cycle++;
			//if (icount_V_cycle > 8) break;
			count_iter_for_film_coef++;
			// В случае задачи Ньютона - Рихмана, Стефана-Больцмана и миксового условия не итерируем до конца обрываем, 
			// т.к. нам требуется частая пересборка матрицы. 13 марта 2016.
			if (((adiabatic_vs_heat_transfer_coeff > DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC) ||
				(breakRUMBAcalc_for_nonlinear_boundary_condition)) &&
				(count_iter_for_film_coef > 1250)) break;

			// 1 dec 2016.
			//  Прерывание после 2 или 5 V циклов обязательно необходимо иначе не будет сходимости.
			if (bvacuumPrism) {
				// 5
				// 250
				if (icount_V_cycle > 2250) break;
			}


			if ((iter_limit == 5000) || ((iVar == PAM) && (fabs(dres) > 7.0e3))) {
#pragma omp parallel for
				for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
					x[i47] = x_copy[i47];
				}
				if (iVar == PAM) {
					std::cout <<  "pressure amendment divergence..." << std::endl;
				}
				std::cout <<  "amg divergence detected dres="<< dres <<"..."<< std::endl;
				std::cout <<  "number V-cycle=" << icount_V_cycle << " initial residual dres0=" << dres_initial << std::endl;

				std::cout <<  "Operator A complexity=" <<dr_grid_complexity<<"  Opertator P complexity="<< static_cast<doublerealT>(nnz_P_memo_all / n_a[0]) <<"..."<<std::endl;
				std::cout <<  "res_best_search=" << res_best_search << std::endl;
				//system("pause");
				// пауза убрана 22 12 2016
				//system("PAUSE");
				break;
			}

			if (iter_limit == 1) {
				// начальная невязка.
				res0start = fabs(dres);
			}

			// Невязка по температуре:
			// НЕТ сходимости для поля температур в гидродинамическом решателе и параметры не помогают.
			//if (iVar == TEMP) std::cout << "temp res=" << fabs(dres) << std::endl;

			if (fabs(dres) < res_best_search)
			{
				// Запоминаем лучшую попытку.
				res_best_search = fabs(dres);
#pragma omp parallel for
				for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
					x_best_search[i47] = static_cast<doublerealT>(x[i47]);
				}
			}
			/*
			if (iVar == PAM) {
			if (fabs(dres) < 1.0) {
			// Идея в том что нам нужна хоть какая-то поправка давления,
			// всё лучше чем тождественно нулевое распределение.
			// невязка при этом у нас менее 1.0 что гарантирует что мы не сильно улетели.
			for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
			x_best_search[i47] = x[i47];
			}
			}
			}
			*/

			// debug 7 июня 2016
			//if (iter_limit > 300) {
			//std::cout << "amg divergense detected...9 june 2016" << std::endl;
			//system("pause");
			//break;
			//}

			//if (dres < 1.0e-14) break;

			// 100
			if (iter_limit > 5000) { // Finned Heat Sink



				if (bfirst_divergence) {
					iter_limit = 3;
					nu1 += 2;
					nu2 += 2;
					nFinestSweeps += 2;
					bfirst_divergence = false;
				}
				else {
					if ((fabs(res_best_search / res0start) < 0.23) && (fabs(res_best_search) < 1.0e-3 * sqrt(static_cast<doublereal>(n_a[0])))) {
						// Если невязка меньше первоначальной на два порядка.
						// Обратное копирование и выход и алгоритма.
#pragma omp parallel for
						for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
							x[i47] = x_best_search[i47];
						}
						break;
					}
					else if ((fabs(res_best_search / res0start) <= 1.0) && (fabs(res_best_search) < 1.0e-4 * sqrt(static_cast<doublereal>(n_a[0])))) {
						// Если невязка меньше первоначальной на два порядка.
						// Обратное копирование и выход и алгоритма.
#pragma omp parallel for
						for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
							x[i47] = x_best_search[i47];
						}
						break;
					}
					// Обратное копирование и выход и алгоритма.
#pragma omp parallel for
					for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
						x[i47] = x_best_search[i47];
					}
					break;
					// закомментировал 11.01.2020
					// как недостижимый код.
					// Эта ветвь кода вообще никогда не вызовется.
					/*
					std::cout << "Fatal amg error: Strong divergence amg solver..."<<  fabs(res_best_search / res0start) << std::endl;
					std::cout << "res_best_search=" << fabs(res_best_search)<<", res0start=" <<  fabs(res0start) << std::endl;
					std::cout << "BiCGStab+ILU2 is start now..." << std::endl;
					std::cout << "please wait...";
					system("pause");
					break; // досрочный выход из while цикла.
					*/
				}
			}
			iter_limit++;

			if (fabs(dres) < fabs(dres_previos)) {
				// все нормально процесс сходится.
				icount_bad_convergence_Vcycles = 0;
			}
			else {
				icount_bad_convergence_Vcycles++;
			}

			//if (_finite(dres) == 0) {
			//if (fabs(dres) > 1.0e30)
			//{
			//std::cout << "\ndivergence AMG detected...solver will be restart... please wait... ... ...\n";
			//std::cout << "\a\a\a\a\a\a\a\a";
			//system("pause");
			//exit(1);
			//return true;
			//for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
			//	x[i47] = x_copy[i47];
			//}
			//if (iter_limit > 100) {
			//	ret_value = true;
			//	break;
			//}
			//else {
			// Увеличение количества сглаживающих итераций ни коим образом не 
			// исправляет факт расходимости. 
			//	nu1++;
			//	nu2++;
			//	nFinestSweeps++;
			// По видимому надо действовать очень тонкой настройкой параметра верхней релаксации omega optimal.
			// Настройка omega optimal должна быть самообучающейся (адаптированной к задаче).
			//}
			//}

			// 24 10 2016
			if (icount_bad_convergence_Vcycles > 40) break;

			if ((icount_bad_convergence_Vcycles >= istop_porog_reconst) || (fabs(dres) / sqrt(static_cast<doublereal>(n_a[0])) > 1.0e30)) {
				// детектировано 10 шагов расходимости подряд по-видимому метод расходится.
				// Также о расходимости говорит невязка большая 1.0e30.

				//if (fabs(dres) < 1.0e-3) break; // Будем считать сходимость достигнута успешно.
				if ((fabs(static_cast<doublereal>(res_best_search / res0start)) < 1.0e-1) && (fabs(dres) / sqrt(static_cast<doublereal>(n_a[0])) < 1.0e-3)) {
					// Если невязка меньше первоначальной на два порядка.
					// Обратное копирование и выход и алгоритма.
#pragma omp parallel for
					for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
						x[i47] = x_best_search[i47];
					}
					std::cout <<  "stagnaion break out" << std::endl;
					break;
				}
				i_count_stagnation++;

				std::cout << "\ndivergence AMG detected...solver will be restart... please wait... ... ..." << std::endl;
				//std::cout << "\a\a\a\a\a\a\a\a";
				//system("pause");
				//exit(1);
				//return true;
				if (i_count_stagnation < 20) {
#pragma omp parallel for
					for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
						//x[i47] = x_copy[i47];
						x[i47] = x_best_search[i47]; // лучшее найденное решение
					}
				}
				if (i_count_stagnation == 20 || i_count_stagnation == 21) gold_const = 0.2;
				if ((i_count_stagnation >= 20) && (i_count_stagnation < 30)) {
					// urd(gen) не работает в параллельном режиме.
//#pragma omp parallel for
					for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
						//x[i47] = x_copy[i47];
						// Можно еще единократно немного улучшить nu1 и nu2.
						doublerealT signumnow = 1.0;
						if (rand() % 2 == 0) signumnow = -1.0;
						//doublereal drand = static_cast<doublereal>(((double)(rand()))/((double)(RAND_MAX+1)));
						doublereal drand = urd(gen);
						x[i47] = signumnow * 1.0 * drand; // Случайное число в интервале от 0 до 1.
					}
				}
				if (i_count_stagnation == 30 || i_count_stagnation == 31) gold_const = 0.2;
				if (i_count_stagnation >= 30) {
#pragma omp parallel for
					for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
						x[i47] = 1.0;
					}
				}
				if (bproblem_amg_convergence1) {
					if (bproblem_amg_convergence2) {
						if (bproblem_amg_convergence3) {
							// выход к вызову BiCGStab+ILU2.
							ret_value = true;
							break;
						}
						else {
							// смена omega.
							bproblem_amg_convergence3 = true;
							icount_bad_convergence_Vcycles = 0;
							buffers3omega = dres / dres_previos;
							std::cout << "buffers1omega="<< buffers1omega <<", buffers2omega="<< buffers2omega <<", buffers3omega="<< buffers3omega <<std::endl; 
						}
					}
					else {
						// смена omega.
						bproblem_amg_convergence2 = true;
						icount_bad_convergence_Vcycles = 0;
						buffers2omega = dres / dres_previos;
						std::cout << "buffers1omega=" << buffers1omega <<", buffers2omega="<< buffers2omega <<std::endl; 
						//istop_porog_reconst += 50; // 10, 20, 30, 40
						// Увеличение количества сглаживающих итераций ничего не даёт.
						//nu1++;
						//nu2++;
						//nFinestSweeps++;
					}
				}
				else {

					bproblem_amg_convergence1 = true; // переход с SOR на стабильный Зейдель.
					icount_bad_convergence_Vcycles = 0;
					buffers1omega = dres / dres_previos;
				}
			}


			dres_previos = dres;


			doublerealT R0_0 = 0.0;
			doublerealT Rprev_0 = 0.0, Rnext_0 = 0.0;
			if (process_flow_logic) {

				//std::cout << "process_flow_logic \n";
				//system("pause");

				// calculate initial residual.
				//residualq<doublereal>(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
				R0_0 = static_cast<doublerealT>(norma(residual_fine[0], n_a[0]));
				Rprev_0 = R0_0;

				// smother
				integer iter = 0;
				for (iter = 0; iter < nu1; ++iter) {
					//quick seidel
					if (bonly_serial) {
						seidelq<doublereal>(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);
					}
					else {
						seidelq<doublereal>(Amat, 1, n_a[0], x, b, nested_desection[0], row_ptr_start, row_ptr_end, 0);
					}
					//residualq(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
					residualq2(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
					Rnext_0 = static_cast<doublerealT>(norma(residual_fine[0], n_a[0]));
					// this is process flow logic
					if (Rnext_0 > process_flow_beta * Rprev_0) {
						// Смысл модификации в том что мы экономим итерации на пресмутере.
						break; // досрочно опускаемся на следующий уровень если он есть конечно.
					}
					else {
						Rprev_0 = Rnext_0;
					}
				}
				if (iter == nu1) {
					std::cout << "level 0 limit presmother iteration is reached" << std::endl;
				}

			}
			else {
				


				// smoother
				for (integer iter = 0; iter < nu1; ++iter) {
					//seidel(Amat, 1, nnz_a[0], x, b, flag, n_a[0]);
					//quick seidel
					if (bonly_serial) {
						if (bILU2smoother == 2) {

						

							seidelq<doublereal>(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, iVar, 0);
							residualq2(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
#pragma omp parallel for
							for (integer i43 = 0; i43 < n_a[0]; ++i43) {
								milu2[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
							}
							lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
							for (integer i43 = 0; i43 < n_a[0]; ++i43) {
								x[i43 + 1] += milu2[0].zbuf2[i43 + 1];
							}
						}
						else {
							seidelq<doublereal>(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, iVar, 0);
						}
					}
					else {
						seidelq<doublereal>(Amat, 1, n_a[0], x, b, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 0);
					}
				}
			}

			//exporttecplot(x, n);

			move_down(nu1, nu2);

			if (!process_flow_logic) {
				// residual_r
				//doublerealT *residual_fine[0] = new doublerealT[n_a[0] + 1];
				//residual(Amat, 1, nnz_a[0], x, b, flag, n_a[0], residual_fine[0]);
				//residualq(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
			}
			dres = static_cast<doublerealT>(norma(residual_fine[0], n_a[0]));
			ret74 += static_cast<real_mix_precision>(fabs(dres));
			if (bprint_mesage_diagnostic) {


				doublereal minx = 1.0e30;
				doublereal maxx = -1.0e30;
				for (integer i_83 = 1; i_83 <= n_a[0]; ++i_83) {
					if (x[i_83] < minx) minx = x[i_83];
					if (x[i_83] > maxx) maxx = x[i_83];
				}

				if (((iVar == TEMP) && (my_amg_manager.istabilization == 3))) {
					//  Сходимость достигнута - досрочный выход из решения нелинейной задачи.
					if ((fabs(minx - minx_gl) < 2.0e-3) && (fabs(maxx - maxx_gl) < 2.0e-3)) {
						if (bprint_mesage_diagnostic) {
							std::cout << "Solution nonlinear problem converged succsefull. Ok..." << std::endl;
						}
						break;
					}
				}

				minx_gl = minx;
				maxx_gl = maxx;
				if (bprint_mesage_diagnostic) {
					if (iiter < 10) {
						std::cout << iiter << "   " << dres << " rho=" << dres / rho << " min=" << minx << " max=" << maxx << "\n";
					}
					else if (iiter < 100) {
						std::cout << iiter << "  " << dres << " rho=" << dres / rho << " min=" << minx << " max=" << maxx << "\n";
					}
					else {
						std::cout << iiter << " " << dres << " rho=" << dres / rho << " min=" << minx << " max=" << maxx << "\n";
					}
				}

				if (!((iVar == TEMP) && (my_amg_manager.istabilization == 3))) {
					if (fabs(1.0 - fabs(dres / rho)) < 1.0e-3) {
						std::cout <<  "stagnation in amg solver determinate ..." << std::endl;
						// Решение не сходится. Получаются сильно заниженные неверные значекния температуры.
						// Возвращаем divergence detected с рекомендацией усилить процесс сходимости 
						// внешним итерационным процессом Крыловского типа а также усилить сглаживатель.
						std::cout << "Recomended: Set up an external Krylov-type\n";
						std::cout << "iterative process BiCGStab or FGMRes.\n";
						return true;
						// 28_10_2016.
						// Осуществляем досрочный выход из итерирования, 
						// т.к. невязка перестала меняться.
						//break;
					}
				}
				if (icount_V_cycle == 1) {
					if (fabs(dres / rho) < 1.0) {
#pragma omp parallel for
						for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
							x_best_search[i47] = static_cast<doublerealT>(x[i47]);
						}
					}
				}
			}
			iiter++;
			// 28.07.2016

			if (fabs(dres) > 1.0e9) {

				std::cout <<  "amg solver divergence detected." << std::endl;
				system("pause");

#pragma omp parallel for
				for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
					///x[i47] = x_best_search[i47];
					//x_copy[i47] = x[i47]; // 4 ноября 2016.
					x[i47] = x_copy[i47];
				}
				residualq2_analysys(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);

				std::cout <<  "dres_initial="<< dres_initial <<" res_best_search="<< res_best_search <<" dres="<< dres <<" current="<< norma(residual_fine[0], n_a[0]) <<std::endl; 
				std::cout << "break. amg divergence detected. fabs(dres) > 1.0e7" << std::endl;
				//system("pause");
				if ((bILU2smoother == 2) || (bILU2smoother == 0)) {
					std::cout <<  "apply ilu" <<  my_amg_manager.lfil <<" smoother for number 0 level" << std::endl;
					equation3DtoCRSRUMBA1(milu2[0], true, Amat, 1, n_a[0], row_ptr_start, row_ptr_end, 0, 0, iVar);

				}
				if (bILU2smoother == 0) {
					// переключение.
					memory_allocation_apostoriory_buffer_ilu(milu2, ilevel - 1);
					bILU2smoother = 2;
				}

				// Это по умолчанию для поправки давления.
				doublerealT dresfinish_probably = static_cast<doublerealT>(0.1 * norma(residual_fine[0], n_a[0]));
				if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) ||
					(iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
					(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
					(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
					(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
					(iVar == GAMMA_LANGTRY_MENTER) ||
					(iVar == RE_THETA_LANGTRY_MENTER)) {
					// Это по умолчанию для компонент скорости внутри SIMPLE алгоритма.
					dresfinish_probably = static_cast<doublerealT>(1.0e-3 * norma(residual_fine[0], n_a[0]));
				}
				if (bSIMPLErun_now_for_temperature) {
					if (iVar == TEMP) {
						// Для поля температур при гидродинамическом расчёте.
						// В BiCGStab Internal 3 домножается на 1e-10.
						dresfinish_probably = static_cast<doublerealT>(1.0e-3 * norma(residual_fine[0], n_a[0]));
					}
				}
				if (iVar == TOTALDEFORMATIONVAR) {
					// Для механических деформаций
					dresfinish_probably = static_cast<doublerealT>(1.0e-3 * norma(residual_fine[0], n_a[0]));
				}
				if (bonly_solid_calculation) {
					dresfinish_probably = static_cast<doublerealT>(1.0e-5 * norma(residual_fine[0], n_a[0]));
				}
				integer i943 = 0;
				for (integer i_prob_detect_i = 0; i_prob_detect_i < 1000; ++i_prob_detect_i) {
					i943 = i_prob_detect_i;
					seidelq<doublereal>(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, iVar,  0);
					residualq2(Amat, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);

					doublereal minx = 1.0e30;
					doublereal maxx = -1.0e30;
					for (integer i_83 = 1; i_83 <= n_a[0]; ++i_83) {
						if (x[i_83] < minx) minx = x[i_83];
						if (x[i_83] > maxx) maxx = x[i_83];
					}

					std::cout <<  i_prob_detect_i << " residual=" << norma(residual_fine[0], n_a[0]) << " min=" << minx << " max=" << maxx << " \n";

					// Досрочный выход из цикла.
					if (norma(residual_fine[0], n_a[0]) < dresfinish_probably) {
						std::cout << "Ok!!! calculation local compleate... " << std::endl;
						break;
					}
#pragma omp parallel for
					for (integer i43 = 0; i43 < n_a[0]; ++i43) {
						milu2[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
					}
					lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
					for (integer i43 = 0; i43 < n_a[0]; ++i43) {
						x[i43 + 1] += milu2[0].zbuf2[i43 + 1];
					}
				}

				// Детектируем возможные проблемы со сходимостью:
				if (norma(residual_fine[0], n_a[0]) >= dresfinish_probably) {
					std::cout << "Fatal error !!! ilu2 divergence detected... "<< std::endl;
					std::cout << "residual curent=" << norma(residual_fine[0], n_a[0]) <<" target residual="<< dresfinish_probably <<std::endl;
					if (i943 < 997) {
						break;
					}
				}
#pragma omp parallel for
				for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
					x_best_search[i47] = static_cast<doublerealT>(x[i47]);
					x_copy[i47] = static_cast<doublerealT>(x[i47]); // 4 ноября 2016.
				}
				system("PAUSE");

				goto FULL_DIVERGENCE_DETECTED;
				//break;
			}

			if (iVar == PAM) {
				if ((fabs(dres / rho) > 0.99999) || (fabs(dres) > 1.0e7)) {
					// Выход из мультигрида если достигнуто 20 циклов расходимости.
					//delta_old_iter = fabs(dres);
					i_signal_break_pam_opening++;
					if (i_signal_break_pam_opening > i_limit_signal_pam_break_opening) {

						if (bprint_mesage_diagnostic) {
							std::cout << "iter = " << iiter << std::endl;
						}

						// Обратное копирование и выход и алгоритма.
#pragma omp parallel for
						for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
							x[i47] = x_best_search[i47];
						}
						break;
					}
				}
			}


			//rho=norma(residual_fine[0], n_a[0]);
			rho = dres;
			// start 08.01.2018
			V_cycle_solve<doublerealT>(Amat, x, b, process_flow_logic, row_ptr_start,
				row_ptr_end, residual_fine, diag, diag_minus_one, n_a, bonly_serial,
				process_flow_beta, F_false_C_true, nu1, nu2, bILU2smoother,
				ilevel, 1, imyinit, maxlevel, milu2, milu0, nested_desection,
				P, nnz_aRP, residual_coarse, igam, nnz_a,
				error_approx_coarse, dapply_ilu_max_pattern_size,
				process_flow_alpha,
				error_approx_fine, nFinestSweeps, ibsp_length, iVar);
			// end 08.01.2018

			//if (bfirst_start_nonlinear_process) {
			// Во избежании расходимости по начальному условию в двойном 
			// вакуумном промежутке.
			//bfirst_start_nonlinear_process = false;
			//break;
			//}
			if (iVar != PAM) {
				if (btheoryGuideANSYSFluent) break; // Делаем лишь один V  цикл.
			}
			//system("pause");
		}
	} // bBiCGStab_plus_RUMBA_camg if (my_amg_manager.istabilization == 1)
	else if (my_amg_manager.istabilization == 1) {
		// Рекомендуется использовать гибридную точность: двойную для BiCGStab и одинарную для предобуславливания с помощью V - цикла.
		// Алгебраический Многосеточный Метод как предобуславливатель
		// к алгоритму Крыловского типа Хенка ван дер Ворста BiCGStab
		// со стабилизацией.
		// Требует ещё одну память под матрицу А на самом подробном уровне.
		// 5.01.2017 Алгоритм BiCGStab изобретён в 1992 году.

		integer inumberVcyclelocbicgstab = 1;

		// нумерация векторов начинается с нуля.
		integer n75 = n_a[0]; // число неизвестных на подробном уровне.
		doublereal* val75 = nullptr;
		val75 = new doublereal[nnz_a[0]];
		integer* col_ind75 = nullptr;
		col_ind75 = new integer[nnz_a[0]];
		integer* row_ptr75 = nullptr;
		row_ptr75 = new integer[n_a[0] + 1];
		if ((val75 == nullptr) || (col_ind75 == nullptr) || (row_ptr75 == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for val, col_ind or row_ptr: bicgStab + camg..." <<std::endl;
			std::cout <<  "Please any key to exit..." <<std::endl;
			exit(1);
		}


		// Преобразовано к формату CRS.

#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON
		// Преобразование к формату CRS.

#pragma omp parallel for
		for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {

			for (integer i_2 = row_ptr_start[i_1]; i_2 <= row_ptr_end[i_1]; ++i_2) {
				//if (Amat.i[i_2] == Amat.j[i_2]) {
					//if (i_1 != Amat.i[i_2]) {
						//std::cout << "err i!=i"<< std::endl;
						//system("PAUSE");
					//}
				if (i_1 == Amat.j[i_2]) {
					val75[i_2 - 1] = diag[0][i_1];
					col_ind75[i_2 - 1] = i_1 - 1;
				}
				else {
					val75[i_2 - 1] = Amat.aij[i_2];
					col_ind75[i_2 - 1] = Amat.j[i_2] - 1;
				}
			}
			row_ptr75[i_1 - 1] = row_ptr_start[i_1] - 1;
		}

		row_ptr75[n_a[0]] = row_ptr_end[n_a[0]];
#else

#if MY_GPU_MATH_LIBRARY_CU_ON

		if (flag != nullptr) {
			free(flag);
			flag = nullptr;
		}
		flag = my_declaration_array<bool>(n_a[0] + 1, false, "flag");

#pragma omp parallel for
		for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {
			flag[i_1] = false;
		}
		calculate_row_ptr_cuda(0, nnz_a[0] - 1,
			row_ptr75,
			flag, Amat_copy);

		//#pragma omp parallel for
			//	for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {
				//	flag[i_1] = false;
				//}

		if (flag != nullptr) {
			free(flag);
			flag = nullptr;
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < nnz_a[0]; ++i_1) {
			val75[i_1] = Amat_copy.aij[i_1];
			col_ind75[i_1] = Amat_copy.j[i_1] - 1;
			//if (i_1 < n_a[0] + 1) {
				//row_ptr75[i_1]--;
			//}
		}

		if (Amat_copy.aij != nullptr) {
			delete[] Amat_copy.aij;
			Amat_copy.aij = nullptr;
		}
		if (Amat_copy.i != nullptr) {
			delete[] Amat_copy.i;
			Amat_copy.i = nullptr;
		}
		if (Amat_copy.j != nullptr) {
			delete[] Amat_copy.j;
			Amat_copy.j = nullptr;
		}

		

#else

		// Преобразование к формату CRS.

#pragma omp parallel for
		for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {

			for (integer i_2 = row_ptr_start[i_1]; i_2 <= row_ptr_end[i_1]; ++i_2) {
				//if (Amat.i[i_2] == Amat.j[i_2]) {
					//if (i_1 != Amat.i[i_2]) {
						//std::cout << "err i!=i"<< std::endl;
						//system("PAUSE");
					//}
				if (i_1 == Amat.j[i_2]) {
					val75[i_2 - 1] = diag[0][i_1];
					col_ind75[i_2 - 1] = i_1 - 1;
				}
				else {
					val75[i_2 - 1] = Amat.aij[i_2];
					col_ind75[i_2 - 1] = Amat.j[i_2] - 1;
				}
			}
			row_ptr75[i_1 - 1] = row_ptr_start[i_1] - 1;
		}

		row_ptr75[n_a[0]] = row_ptr_end[n_a[0]];


#endif
#endif

		


		// Вектора необходимые для работы BiCGStab.
		doublereal* ri75 = nullptr;
		doublereal* roc75 = nullptr;
		doublereal* s75 = nullptr;
		doublereal* t75 = nullptr;
		doublereal* vec75 = nullptr;
		doublereal* vi75 = nullptr;
		doublereal* pi75 = nullptr;
		doublereal* dx75 = nullptr;
		doublereal* dax75 = nullptr;
		doublereal* y75 = nullptr;
		doublereal* z75 = nullptr;
		// Первое предобуславливание:
		doublereal* y76 = nullptr;
		doublereal* pi76 = nullptr;
		y76 = new doublereal[n75 + 1];
		pi76 = new doublereal[n75 + 1];
		// Второе предобуславливание:
		doublereal* z76 = nullptr;
		doublereal* s76 = nullptr;
		z76 = new doublereal[n75 + 1];
		s76 = new doublereal[n75 + 1];

		ri75 = new doublereal[n75];
		roc75 = new doublereal[n75];
		s75 = new doublereal[n75];
		t75 = new doublereal[n75];
		vec75 = new doublereal[n75];
		vi75 = new doublereal[n75];
		pi75 = new doublereal[n75];
		dx75 = new doublereal[n75];
		dax75 = new doublereal[n75];
		y75 = new doublereal[n75];
		z75 = new doublereal[n75];
		if ((ri75 == nullptr) || (roc75 == nullptr) || (s75 == nullptr) || (t75 == nullptr) ||
			(vi75 == nullptr) || (pi75 == nullptr) || (dx75 == nullptr) || (dax75 == nullptr) ||
			(y75 == nullptr) || (z75 == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for: bicgStab + camg..."<< std::endl;
			std::cout << "Please any key to exit..." << std::endl;
			exit(1);
		}

		integer iflag75 = 1, icount75 = 0;
		doublereal delta075 = 1.0e30, deltai75 = 1.0e30;
		doublereal bet75 = 0.0, roi75 = 0.0;
		doublereal roim175 = 1.0, al75 = 1.0, wi75 = 1.0;

		doublereal epsilon75 = dterminatedTResudual;  // точность вычисления
		if (iVar == TEMP) {
			epsilon75 *= 1.0e-4; // 1.0e-4
		}
		if (iVar == TOTALDEFORMATIONVAR) {
			epsilon75 *= 1.0e-4; // 1.0e-4
								 //epsilon75 *= 1.0e-12;
		}
		

		// initialize
#pragma omp parallel for 
		for (integer i75 = 0; i75 < n75; ++i75) {
			s75[i75] = 0.0;
			t75[i75] = 0.0;
			vi75[i75] = 0.0;
			pi75[i75] = 0.0;
			// инициализатор массивов для предобуславливания
			y75[i75] = 0.0;
			z75[i75] = 0.0;
			// результат умножения матрицы на вектор.
			dax75[i75] = 0.0;
			// Начальное приближение.
			dx75[i75] = x[i75 + 1];
		}


#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON

		// Умножение матрицы на вектор. Нумерации векторов начинаются с нуля.
		MatrixCRSByVector(val75, col_ind75, row_ptr75, dx75, dax75, n75); // результат занесён в  dax75

#else
#if MY_GPU_MATH_LIBRARY_CU_ON

		//30.07.2021

		//format sparse matrix;  time
		// CRS 3min 17s 610ms
		// ELL 3min 21s 310ms
		// CPU CRS 3min 20s 10ms




		if (init_b_first_device) {

			int device;
			device = idevice_Tesla;

			// Choose which GPU to run on, change this on a multi-GPU system.
			cudaStatus = cudaSetDevice(device);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
				system("PAUSE");
				exit(1);
			}

			cudaDeviceProp prop;
			cudaGetDeviceProperties(&prop, device);
			std::cout <<  prop.name << std::endl;
			// GeForce 840M 384потока 1ГГц каждый март 2014 года. 28нм.

			init_b_first_device = false;
		}


#if EllFormat
		doublereal* data = nullptr;
		integer* indices = nullptr;
		integer string_size = -1;

		doublereal* data_bound = nullptr;
		integer* indices_bound = nullptr;
		integer string_size_bound = 2;

		for (integer i77 = 0; i77 < n75; ++i77) {
			if (row_ptr75[i77 + 1] - row_ptr75[i77] > string_size) {
				string_size = row_ptr75[i77 + 1] - row_ptr75[i77];
			}
			//if (i77>= maxelm) {
				//std::cout << "\nstart hvost=" << i77 << "  maxelm=" << maxelm << "start "<< row_ptr75[i77]  << "end=" << row_ptr75[i77 + 1] << std::endl;
				//system("pause");
			//}
		}

		data = new doublereal[string_size * maxelm_out];
		indices = new integer[string_size * maxelm_out];

		data_bound = new doublereal[string_size_bound * (n - maxelm_out)];
		indices_bound = new integer[string_size_bound * (n - maxelm_out)];


		// Перепаковка в Ell формат.
#pragma omp parallel for
		for (integer i77 = 0; i77 < string_size * maxelm_out; ++i77) {
			data[i77] = 0.0;
		}

#pragma omp parallel for
		for (integer i77 = 0; i77 < string_size_bound * (n - maxelm_out); ++i77) {
			data_bound[i77] = 0.0;
		}

#pragma omp parallel for
		for (integer i77 = 0; i77 < maxelm_out; ++i77) {

			const integer rowend = row_ptr75[i77 + 1];
			const integer rowbeg = row_ptr75[i77];

			integer jstart = 0;
			for (integer j77 = rowbeg; j77 < rowend; ++j77)
			{
				const integer element_offset = i77 + jstart * maxelm_out;
				data[element_offset] = val75[j77];
				indices[element_offset] = col_ind75[j77];

				++jstart;
			}
			for (integer j77 = jstart; j77 < string_size; ++j77)
			{
				const integer element_offset = i77 + j77 * maxelm_out;

				if (indices[element_offset - 1] < n75 - 2) {
					indices[element_offset] = indices[element_offset - 1] + 1;
					data[element_offset] = 0.0;
				}
				else {
					indices[element_offset] = indices[element_offset - 1];
					data[element_offset] = 0.0;
				}

			}
		}

#pragma omp parallel for
		for (integer i77 = maxelm_out; i77 < n; ++i77) {

			const integer rowend = row_ptr75[i77 + 1];
			const integer rowbeg = row_ptr75[i77];

			integer jstart = 0;
			for (integer j77 = rowbeg; j77 < rowend; ++j77)
			{
				const integer element_offset = i77 - maxelm_out + jstart * (n - maxelm_out);
				data_bound[element_offset] = val75[j77];
				indices_bound[element_offset] = col_ind75[j77];

				++jstart;
			}
			for (integer j77 = jstart; j77 < string_size_bound; ++j77)
			{
				const integer element_offset = i77 - maxelm_out + j77 * (n - maxelm_out);

				if (indices_bound[element_offset - 1] < n75 - 2) {
					indices_bound[element_offset] = indices_bound[element_offset - 1] + 1;
					data_bound[element_offset] = 0.0;
				}
				else {
					indices_bound[element_offset] = indices_bound[element_offset - 1];
					data_bound[element_offset] = 0.0;
				}

			}
		}


		// Освобождение оперативной памяти.
		if (val75 != nullptr) {
			delete[] val75;
			val75 = nullptr;
		}
		if (col_ind75 != nullptr) {
			delete[] col_ind75;
			col_ind75 = nullptr;
		}
		if (row_ptr75 != nullptr) {
			delete[] row_ptr75;
			row_ptr75 = nullptr;
		}


		doublereal* dev_data = nullptr;
		integer* dev_indices = nullptr;

		doublereal* dev_data_bound = nullptr;
		integer* dev_indices_bound = nullptr;

#else
		doublereal* dev_val = nullptr;
		integer* dev_col_ind = nullptr;
		integer* dev_row_ptr = nullptr;
#endif


		doublereal* dev_V = nullptr;
		doublereal* dev_tmp = nullptr;

#if EllFormat
		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_data, (maxelm_out * string_size) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_val failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_indices, (maxelm_out * string_size) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_col_ind failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_data_bound, ((n - maxelm_out) * string_size_bound) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_val failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_indices_bound, ((n - maxelm_out) * string_size_bound) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_col_ind failed!");
			system("PAUSE");
			exit(1);
		}

#else

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_val, (nnz) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_val failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_col_ind, (nnz) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_col_ind failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_row_ptr, (n75 + 1) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_row_ptr failed!");
			system("PAUSE");
			exit(1);
		}

#endif

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_V, (n75) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_V failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_tmp, (n75) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_tmp failed!");
			system("PAUSE");
			exit(1);
		}

#if EllFormat

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_data, data, maxelm_out * string_size * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_val, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_indices, indices, maxelm_out * string_size * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_col_ind, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		delete[] data;
		data = nullptr;
		delete[] indices;
		indices = nullptr;

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_data_bound, data_bound, (n - maxelm_out) * string_size_bound * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_data_bound, data_bound HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_indices_bound, indices_bound, (n - maxelm_out) * string_size_bound * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_indices_bound, indices_bound HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		delete[] data_bound;
		data_bound = nullptr;
		delete[] indices_bound;
		indices_bound = nullptr;

#else

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_val, val75, nnz * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_val, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_col_ind, col_ind75, nnz * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_col_ind, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}


		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_row_ptr, row_ptr75, (n75 + 1) * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_rho_ptr, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Освобождение оперативной памяти.
		if (val75 != nullptr) {
			delete[] val75;
			val75 = nullptr;
		}
		if (col_ind75 != nullptr) {
			delete[] col_ind75;
			col_ind75 = nullptr;
		}
		if (row_ptr75 != nullptr) {
			delete[] row_ptr75;
			row_ptr75 = nullptr;
		}

#endif

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_V, dx75, n75 * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_V, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		/*cudaStatus = cudaMemcpy(dev_tmp, dax75, n75 * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_tmp, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}*/

#if EllFormat

		MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data, dev_indices, string_size, dev_V, dev_tmp, maxelm_out, 0);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data_bound, dev_indices_bound, string_size_bound, dev_V, dev_tmp, (n - maxelm_out), maxelm_out);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

#else

		MatrixCRSByVectorKernel << <128, 128 >> > (dev_val, dev_col_ind, dev_row_ptr, dev_V, dev_tmp, n75);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "MatrixCRSByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}


#endif

		cudaStatus = cudaMemcpy(dax75, dev_tmp, n75 * sizeof(doublereal), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy tmp, dev_tmp DeviceToHost failed!");
			system("PAUSE");
			exit(1);
		}


#else

		// Умножение матрицы на вектор. Нумерации векторов начинаются с нуля.
		MatrixCRSByVector(val75, col_ind75, row_ptr75, dx75, dax75, n75); // результат занесён в  dax75
																		 
#endif
#endif

		 // Вычисление ri75 и roc75.

		// Вычисление ri75 и roc75.
#pragma omp parallel for
		for (integer i75 = 0; i75 < n75; ++i75) {
			ri75[i75] = b[i75 + 1] - dax75[i75];
			roc75[i75] = 1.0;
		}
		delta075 = NormaV(ri75, n75);


		// Если решение сразу хорошее то не считать:
		if (iVar == TEMP) {
			if (fabs(delta075) < 1.0e-4 * dterminatedTResudual) iflag75 = 0;
		}
		else {
			if (fabs(delta075) < dterminatedTResudual) iflag75 = 0;
		}
		integer iflag175 = 1;
		if (iVar == TURBULENT_KINETIK_ENERGY) {
			if (fabs(delta075) < 1e-21) iflag175 = 0;
		}
		else {
			if (fabs(delta075) < 1e-14) iflag175 = 0;
		}
		if ((iVar == TEMP) && (iflag75 == 0) && (iflag175 == 0)) {
			std::cout << "bicgStab+camg: iflag="<< iflag75 <<", iflag1="<< iflag175 <<", delta0="<< delta075 <<std::endl; 
			system("PAUSE");
		}
		if ((iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
			(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
			(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) ||
			(iVar == RE_THETA_LANGTRY_MENTER)) {
			std::cout << "Turbulence equations: bicgStab+camg: iflag=" << iflag75 << ", iflag1=" << iflag175 << ", delta0=" << delta075 <<std::endl; 
		}
		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) {
			std::cout << "Velocity equations: bicgStab+camg: iflag=" << iflag75 << ", iflag1=" << iflag175 << ", delta0=" << delta075 << std::endl; 
		}

		integer iN75 = 10;
		if (n75 < 30000) {
			// задача очень малой размерности !
			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) ||
				(iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
				(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
				(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
				(iVar == GAMMA_LANGTRY_MENTER) ||
				(iVar == RE_THETA_LANGTRY_MENTER)) {
				iN75 = 1; // обязательно нужна хотя бы одна итерация.
						  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
				if (1.0e-3 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-3 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
			}
			if (iVar == TEMP) {
				iN75 = 2;
				epsilon75 = fmin(0.1 * fabs(delta075), epsilon75);
				if (bSIMPLErun_now_for_temperature) {
					//std::cout << "epsilon=" << epsilon << std::endl;
					//system("PAUSE");
					// Экспериментальным образом обнаружена недоэтерированность по температуре для гидродинамического решателя.
					// поэтому точность было решено увеличить на 5 порядков.
					// 27.07.2016
					epsilon75 *= 1e-10;
					iN75 = 20;
					//epsilon75 *= 1e-16;
					//iN75 = 30;
				}
			}
			if (iVar == PAM) {
				iN75 = 12; // решение для поправки давления должно быть получено точно.

				if (1.0e-10 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-10 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
				//std::cout << epsilon75 << std::endl; system("pause");
			}
		}
		else if (n75 < 100000) {
			// Здесь я немного увеличил число итераций и 
			// скорректировал условие окончания чтобы считало 
			// поточнее, но это не повлияло.
			// Главный вопрос в том что невязка по температуре почему-то не меняется.
			// задача небольшой размерности.
			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) ||
				(iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
				(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
				(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
				(iVar == GAMMA_LANGTRY_MENTER) ||
				(iVar == RE_THETA_LANGTRY_MENTER)) {
				iN75 = 3; // обязательно нужна хотя бы одна итерация.
						  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
				if (1.0e-3 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-3 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
				// 27.07.2016
				iN75 = 12;
				epsilon75 *= 1e-2;
			}
			if (iVar == TEMP) {
				iN75 = 4;
				epsilon75 = fmin(0.1 * fabs(delta075), epsilon75);
				if (bSIMPLErun_now_for_temperature) {
					//std::cout << "epsilon75="<< epsilon75 << std::endl;
					//system("PAUSE");
					// Экспериментальным образом обнаружена недоэтерированность по температуре для гидродинамического решателя.
					// поэтому точность было решено увеличить на 5 порядков.
					// 27.07.2016
					epsilon75 *= 1e-10;
					iN75 = 20;
					//epsilon75 *= 1e-16;
					//iN75 = 30;
				}
			}
			if (iVar == PAM) {
				iN75 = 20; // решение для поправки давления должно быть получено точно.

				if (1.0e-10 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-10 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
				//std::cout << epsilon75; system("PAUSE");
				// 27.07.2016.
				//epsilon75 *= 1e-2;
				//iN75 = 20;
			}
		}
		else if (n75 < 300000) {
			// задача небольшой средней размерности.
			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) ||
				(iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
				(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
				(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
				(iVar == GAMMA_LANGTRY_MENTER) ||
				(iVar == RE_THETA_LANGTRY_MENTER)) {
				iN75 = 3; // обязательно нужна хотя бы одна итерация.
						  // Вообще говоря невязка для скоростей падает очень быстро поэтому всегда достаточно iN итераций для скорости.
						  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
				if (1.0e-3 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-3 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
			}
			if (iVar == TEMP) {
				iN75 = 4;
				epsilon75 = fmin(0.1 * fabs(delta075), epsilon75);
				if (bSIMPLErun_now_for_temperature) {
					//std::cout << "epsilon75=" << epsilon75 << std::endl;
					//system("PAUSE");
					// Экспериментальным образом обнаружена недоэтерированность по температуре для гидродинамического решателя.
					// поэтому точность было решено увеличить на 5 порядков.
					// 27.07.2016
					epsilon75 *= 1e-10;
					iN75 = 20;
					//epsilon75 *= 1e-16;
					//iN75 = 30;
				}
			}
			if (iVar == PAM) {
				iN75 = 30; // решение для поправки давления должно быть получено точно.
				if (1.0e-10 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-10 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
				//std::cout << epsilon75; system("PAUSE");
			}
		}
		else if (n75 < 1000000) {
			// задача истинно средней размерности.
			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) ||
				(iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
				(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
				(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
				(iVar == GAMMA_LANGTRY_MENTER) ||
				(iVar == RE_THETA_LANGTRY_MENTER)) {
				iN75 = 3; // обязательно нужна хотя бы одна итерация.
						  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
				if (1.0e-3 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-3 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
			}
			if (iVar == TEMP) {
				iN75 = 4;
				epsilon75 = 1e-5 * fmin(0.1 * fabs(delta075), epsilon75);
				if (bSIMPLErun_now_for_temperature) {
					//std::cout << "epsilon75=" << epsilon75 << std::endl;
					//system("PAUSE");
					// Экспериментальным образом обнаружена недоэтерированность по температуре для гидродинамического решателя.
					// поэтому точность было решено увеличить на 5 порядков.
					// 27.07.2016
					epsilon75 *= 1e-8;
					iN75 = 20;
					//epsilon75 *= 1e-16;
					//iN75 = 30;
				}
			}
			if (iVar == PAM) {
				iN75 = 16; // решение для поправки давления должно быть получено точно.
				if (1.0e-4 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-4 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
				//std::cout << epsilon75; system("PAUSE");
			}
		}
		else if (n75 < 3000000) {
			// задача достаточно большой размерности.
			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) ||
				(iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
				(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
				(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
				(iVar == GAMMA_LANGTRY_MENTER) ||
				(iVar == RE_THETA_LANGTRY_MENTER)) {
				iN75 = 6; // обязательно нужна хотя бы одна итерация.
						  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
				if (1.0e-3 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-3 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
			}
			if (iVar == TEMP) {
				iN75 = 8;
				//epsilon75 = 1e-5 * fmin(0.1 * fabs(delta075), epsilon75);

				// Рекомендация Федюшкина А.И. и И.М. Аболдуева.
				epsilon75 = 1.0e-10;
				if (bSIMPLErun_now_for_temperature) {
					//std::cout << "epsilon75="<< epsilon75 << std::endl;
					//system("PAUSE");
					// Экспериментальным образом обнаружена недоэтерированность по температуре для гидродинамического решателя.
					// поэтому точность было решено увеличить на 5 порядков.
					// 27.07.2016
					epsilon75 *= 1e-8;
					iN75 = 20;
					//epsilon75 *= 1e-16;
					//iN75 = 30;
				}
			}
			if (iVar == PAM) {
				iN75 = 23; // решение для поправки давления должно быть получено точно.
				if (1.0e-4 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-4 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
				//std::cout << epsilon75; system("PAUSE");
			}
		}
		else if (n75 >= 3000000) {
			// задача очень большой размерности.
			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) ||
				(iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
				(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
				(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
				(iVar == GAMMA_LANGTRY_MENTER) ||
				(iVar == RE_THETA_LANGTRY_MENTER)) {
				iN75 = 6; // обязательно нужна хотя бы одна итерация.
						  // если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
				if (1.0e-3 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-3 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
			}
			if (iVar == TEMP) {
				iN75 = 8;
				//epsilon75 = 1e-10 * fmin(0.1 * fabs(delta075), epsilon75);
				
				
				// Рекомендация Федюшкина А.И. и И.М. Аболдуева.
				epsilon75 = 1.0e-10;
			}
			if (iVar == PAM) {
				iN75 = 36; // решение для поправки давления должно быть получено точно.
				if (1.0e-4 * fabs(delta075) < epsilon75) {
					epsilon75 = 1.0e-4 * fabs(delta075);
				}
				if (iflag175 == 1) {
					iflag75 = 1;
				}
				//std::cout << epsilon75; system("PAUSE");
			}
		}

		integer maxit75 = 2000;
		if (iVar == TEMP) {
			maxit75 = 2000;
		}
		if (iVar == PAM) {
			maxit75 = 2000; // 2000
		}
		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) {
			maxit75 = 100; // 100
		}
		if ((iVar == NUSHA) || (iVar == TURBULENT_KINETIK_ENERGY) ||
			(iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) ||
			(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) ||
			(iVar == RE_THETA_LANGTRY_MENTER))
		{
			maxit75 = 100; // 100
		}
		if (iVar == TOTALDEFORMATIONVAR) {
			maxit75 = 800; // 2000
			if (1.0e-4 * fabs(delta075) < epsilon75) {
				epsilon75 = 1.0e-4 * fabs(delta075);
			}
			epsilon75 = 1.0e-12;
			iN75 = 8; // Количество обязательных итераций.
			if (iflag175 == 1) {
				iflag75 = 1;
			}

		}

		// Если число расходимостей превысит оговорённую константу то произойдёт выход из алгоритма.
		integer i_signal_break_pam_opening75 = 0;
		// x хорошее значение.
		const integer i_limit_signal_pam_break_opening75 = 8000;//20
		doublereal delta_old_iter75 = 1.0e10;

		integer count_iter_for_film_coef75 = 0;

		// Используется для досрочного прерывания вычислительного процесса
		// как в алгоритме FGMRES Юсефа Саада и Мартина Г. Шульца.
		doublereal norma_b = NormaV_for_gmres(b, n75);

		// Мы обязательно должны сделать несколько итераций. (не менее 10).
		// Если только решение не удовлетворяет уравнению тождественно.
		//if (iVar == TURBULENT_KINETIK_ENERGY) {
			//std::cout << "TURBULENT_KINETIK_ENERGY: iN75=="<< iN75 <<" iflag75==" << iflag75 << " iflag175==" << iflag175 << " maxit75=" << maxit75 << "\n delta075="<<delta075<<"  epsilon75=" << epsilon75<< "\n";
			//system("pause");
		//}
		//if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
			//std::cout << "TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: iN75=="<< iN75 <<" iflag75==" << iflag75 << " iflag175==" << iflag175 << " maxit75=" << maxit75 << "\n delta075="<<delta075<<"  epsilon75=" << epsilon75<< "\n";
			//system("pause");
		//}
		while (((icount75 < iN75) && (iflag175 != 0)) || (iflag75 != 0 && icount75 < maxit75)) {

			// 6.01.2017: Body BiCGStab + AMG. (BiCGStab_internal4).

			++icount75;

			++count_iter_for_film_coef75;
			// В случае задачи Ньютона - Рихмана, Стефана-Больцмана и миксового условия не итерируем до конца обрываем, 
			// т.к. нам требуется частая пересборка матрицы. 13 марта 2016.
			//if (((adiabatic_vs_heat_transfer_coeff > ADIABATIC_WALL_BC) || (breakRUMBAcalc_for_nonlinear_boundary_condition)) && (count_iter_for_film_coef75>5)) break;

			roi75 = Scal(roc75, ri75, n75);
			bet75 = (roi75 / roim175) * (al75 / wi75);


			//std::cout << "roi75=="<< roi75 << ", roim175=="<< roim175<< ", al75==" << al75 << ", wi75==" << wi75 << std::endl;
			//system("pause");

#pragma omp parallel for 
			for (integer i75 = 0; i75 < n75; ++i75) {
				doublereal pibuf75 = ri75[i75] + (pi75[i75] - vi75[i75] * wi75) * bet75;
				pi75[i75] = pibuf75;
			}

			// Первое предобуславливание.
			// Ky=pi
#pragma omp parallel for 
			for (integer i75 = 0; i75 < n75; ++i75) {
				y75[i75] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.
				y76[i75 + 1] = 0.0;
				pi76[i75 + 1] = pi75[i75];
			}

			// multigrid RUMBA preconditioner
			// Вставлено 6.01.2017 begin
			// одного V цикла недостаточно.
			// A*y76=pi76;
			V_cycle_solve<doublerealT>(Amat, y76, pi76, process_flow_logic, row_ptr_start,
				row_ptr_end, residual_fine, diag, diag_minus_one, n_a, bonly_serial,
				process_flow_beta, F_false_C_true, nu1, nu2, bILU2smoother,
				ilevel, inumberVcyclelocbicgstab, imyinit, maxlevel, milu2, milu0, nested_desection,
				P, nnz_aRP, residual_coarse, igam, nnz_a,
				error_approx_coarse, dapply_ilu_max_pattern_size,
				process_flow_alpha,
				error_approx_fine, nFinestSweeps, ibsp_length, iVar);
			// Вставлено 6.01.2017 end

			// Возвращение результата.
#pragma omp parallel for
			for (integer i75 = 0; i75 < n75; ++i75) {
				y75[i75] = y76[i75 + 1];
			}


#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON
			MatrixCRSByVector(val75, col_ind75, row_ptr75, y75, vi75, n75); // vi==A*y;

#else

#if MY_GPU_MATH_LIBRARY_CU_ON

			// Copy input vectors from host memory to GPU buffers.
			cudaStatus = cudaMemcpy(dev_V, y75, n75 * sizeof(doublereal), cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy dev_tmp, val HostToDevice failed!");
				system("PAUSE");
				exit(1);
			}

#if EllFormat

			MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data, dev_indices, string_size, dev_V, dev_tmp, maxelm_out, 0);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data_bound, dev_indices_bound, string_size_bound, dev_V, dev_tmp, (n - maxelm_out), maxelm_out);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

#else

			MatrixCRSByVectorKernel << <128, 128 >> > (dev_val, dev_col_ind, dev_row_ptr, dev_V, dev_tmp, n75);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "MatrixCRSByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

#endif

			cudaStatus = cudaMemcpy(vi75, dev_tmp, n75 * sizeof(doublereal), cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy tmp, dev_tmp DeviceToHost failed!");
				system("PAUSE");
				exit(1);
			}

#else

			MatrixCRSByVector(val75, col_ind75, row_ptr75, y75, vi75, n75); // vi==A*y;

#endif
#endif


			const doublereal scal_val1 = Scal(roc75, vi75, n75);

			if ((fabs(roi75) < 1e-30) && (fabs(scal_val1) < 1e-30)) {
				al75 = 1.0;
			}
			else if (fabs(roi75) < 1e-30) {
				al75 = 0.0;
			}
			else {
				al75 = roi75 / scal_val1;
			}


#pragma omp parallel for 
			for (integer i75 = 0; i75 < n75; ++i75) {
				s75[i75] = ri75[i75] - al75 * vi75[i75];
			}

			// Второе предобуславливание.
			// Kz=s

#pragma omp parallel for 
			for (integer i75 = 0; i75 < n75; ++i75) {
				z75[i75] = 0.0; // Если начинать не с нуля то не будет сходимости для PAM !.
			}

#pragma omp parallel for 
			for (integer i75 = 0; i75 < n75; ++i75) {
				vec75[i75] = s75[i75];
				z76[i75 + 1] = 0.0;
				s76[i75 + 1] = s75[i75];
			}

			// multigrid RUMBA preconditioner
			// Вставлено 6.01.2017 begin
			// одного V цикла недостаточно.
			// A*z76=s76;
			V_cycle_solve<doublerealT>(Amat, z76, s76, process_flow_logic, row_ptr_start,
				row_ptr_end, residual_fine, diag, diag_minus_one, n_a, bonly_serial,
				process_flow_beta, F_false_C_true, nu1, nu2, bILU2smoother,
				ilevel, inumberVcyclelocbicgstab, imyinit, maxlevel, milu2, milu0, nested_desection,
				P, nnz_aRP, residual_coarse, igam, nnz_a,
				error_approx_coarse, dapply_ilu_max_pattern_size,
				process_flow_alpha,
				error_approx_fine, nFinestSweeps, ibsp_length, iVar);
			// Вставлено 6.01.2017 end

#pragma omp parallel for 
			for (integer i75 = 0; i75 < n75; ++i75) {
				s75[i75] = vec75[i75];
				// Возвращаем результат.
				z75[i75] = z76[i75 + 1];
			}


#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON
			MatrixCRSByVector(val75, col_ind75, row_ptr75, z75, t75, n75); // t==A*z;

#else

#if MY_GPU_MATH_LIBRARY_CU_ON

			// Copy input vectors from host memory to GPU buffers.
			cudaStatus = cudaMemcpy(dev_V, z75, n75 * sizeof(doublereal), cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy dev_tmp, val HostToDevice failed!");
				system("PAUSE");
				exit(1);
			}

#if EllFormat

			MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data, dev_indices, string_size, dev_V, dev_tmp, maxelm_out, 0);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data_bound, dev_indices_bound, string_size_bound, dev_V, dev_tmp, (n - maxelm_out), maxelm_out);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

#else

			MatrixCRSByVectorKernel << <128, 128 >> > (dev_val, dev_col_ind, dev_row_ptr, dev_V, dev_tmp, n75);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "MatrixCRSByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

#endif

			cudaStatus = cudaMemcpy(t75, dev_tmp, n75 * sizeof(doublereal), cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy tmp, dev_tmp DeviceToHost failed!");
				system("PAUSE");
				exit(1);
			}

#else

			MatrixCRSByVector(val75, col_ind75, row_ptr75, z75, t75, n75); // t==A*z;

#endif
#endif

			

			wi75 = Scal(t75, s75, n75) / Scal(t75, t75, n75);
			// std::cout << "Scal(t75,s75,n75)==" << Scal(t75,s75,n75) << ", Scal(t75,t75,n75)=="<< Scal(t75,t75,n75) << std::endl;

#pragma omp parallel for 
			for (integer i75 = 0; i75 < n75; ++i75) {
				//dx75[i75]+=al75*pi75[i75]+wi75*s75[i75]; // так было без предобуславливателя
				dx75[i75] += al75 * y75[i75] + wi75 * z75[i75]; // так стало с предобуславливателем
				ri75[i75] = s75[i75] - wi75 * t75[i75];
			}
			deltai75 = NormaV(ri75, n75);

			//std::cout << "deltai75=" << deltai75 << std::endl; system("PAUSE");

			// печать невязки на консоль
			if (bprint_mesage_diagnostic) {
				if ((icount75 % 10) == 0) {
					std::cout << "iter  residual" << std::endl;
				}

			std::cout << icount75 << " " << deltai75 << std::endl;
			}
			else {

				// 28.07.2016.
				std::cout << icount75 << " " << deltai75 << std::endl;
			}

			//system("pause");
			if (deltai75 > delta_old_iter75) i_signal_break_pam_opening75++;
			delta_old_iter75 = deltai75;
			if (iVar == PAM) {
				if (i_signal_break_pam_opening75 > i_limit_signal_pam_break_opening75) {
					// досрочный выход из цикла.
					if (bprint_mesage_diagnostic) {
						std::cout << "icount PAM=" << icount75 << "\n";
					}

					break;
				}
			}

			if (deltai75 < epsilon75) iflag75 = 0; // конец вычисления
			else roim175 = roi75;

			if (iVar == TEMP) {
				//std::cout << "epsilon=" << epsilon75 <<" deltai="<< deltai75<< " icount="<< icount75 << std::endl;
				//system("pause");
			}

			//04.04.2019
			// Успешное условие окончания вычислительного процесса следуя алгоритму FGMRES Ю.Саада.
			if (0 && ((NormaV_for_gmres(ri75, n75) / norma_b) <= 0.1 * dterminatedTResudual)) {
				iflag75 = 0; // конец вычисления
				if (bprint_mesage_diagnostic) {
					std::cout << "dosrochnji vjhod" << std::endl;
				}
				icount_V_cycle = icount75; // количество итераций в BiCGStabP для лога.
				break;
			}

			icount_V_cycle = icount75; // количество итераций в BiCGStabP для лога.

			if (icount75 > 2600) break; // 15.02.2017

		}


		if (icount75 == maxit75) {
			// amg solver divergence detected 
			// Превышено допустимое число итераций.
			return true;
		}

		// Возвращение результата вычислений.
#pragma omp parallel for
		for (integer i75 = 0; i75 < n75; ++i75) {
			x[i75 + 1] = dx75[i75];
			x_best_search[i75 + 1] = static_cast<doublerealT>(dx75[i75]);
		}

		// Освобождение оперативной памяти.
		
#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON

#else
#if MY_GPU_MATH_LIBRARY_CU_ON

#if EllFormat

		cudaFree(dev_data);
		cudaFree(dev_indices);

		delete[] data;
		delete[] indices;

		cudaFree(dev_data_bound);
		cudaFree(dev_indices_bound);

		delete[] data_bound;
		delete[] indices_bound;

#else


		cudaFree(dev_row_ptr);
		cudaFree(dev_val);
		cudaFree(dev_col_ind);

#endif

		cudaFree(dev_V);
		cudaFree(dev_tmp);

#endif
#endif
		
		// Первое предобуславливание
		
			delete[] pi76;
			pi76 = nullptr;
			delete[] y76;
			y76 = nullptr;
			
		// Второе предобуславливание
		
			delete[] z76;
			z76 = nullptr;		
			delete[] s76;
			s76 = nullptr;
		
			delete[] ri75;
			ri75 = nullptr;		
			delete[] roc75;
			roc75 = nullptr;
		
			delete[] s75;
			s75 = nullptr;		
			delete[] t75;
			t75 = nullptr;
		
			delete[] vec75;
			vec75 = nullptr;		
			delete[] vi75;
			vi75 = nullptr;
		
		
			delete[] pi75;
			pi75 = nullptr;
			delete[] dx75;
			dx75 = nullptr;
			
		
			delete[] dax75;
			dax75 = nullptr;		
			delete[] y75;
			y75 = nullptr;
		
			delete[] z75;
			z75 = nullptr;
		

		// Освобождение оперативной памяти.
		
			delete[] val75;
			val75 = nullptr;		
			delete[] col_ind75;
			col_ind75 = nullptr;		
			delete[] row_ptr75;
			row_ptr75 = nullptr;
		

	}
	else if (my_amg_manager.istabilization == 2)
	{   //  09.01.2018
		// Рекомендуется использовать гибридную точность: двойную для FGMRES и одинарную для предобуславливания с помощью V - цикла.
		// FGMRes if (my_amg_manager.istabilization == 2)
		// Гибкий вариант обобщённого метода минимальных невязок в котором на каждой итерации
		// однократно применяется многосеточный предобуславливатель. Алгорим Юсефа Саада и Мартина Г. Шульца 1986 года.

		integer inumberVcyclelocbicgstab = 1;

		// нумерация векторов начинается с нуля.
		integer n75 = n_a[0]; // число неизвестных на подробном уровне.
		doublereal* val75 = nullptr;
		val75 = new doublereal[nnz_a[0]];
		integer* col_ind75 = nullptr;
		col_ind75 = new integer[nnz_a[0]];
		integer* row_ptr75 = nullptr;
		row_ptr75 = new integer[n_a[0] + 1];
		if ((val75 == nullptr) || (col_ind75 == nullptr) || (row_ptr75 == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for val, col_ind or row_ptr: FGMRes + camg..." << std::endl;
			std::cout << "Please any key to exit..." << std::endl;
			exit(1);
		}

		// инициализация матрицы.
		
		
		// Преобразовано к формату CRS.

#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON

		// Преобразование к формату CRS.

#pragma omp parallel for
		for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {

			for (integer i_2 = row_ptr_start[i_1]; i_2 <= row_ptr_end[i_1]; ++i_2) {
				//if (Amat.i[i_2] == Amat.j[i_2]) {
					//if (i_1 != Amat.i[i_2]) {
						//std::cout << "err i!=i"<< std::endl;
						//system("PAUSE");
					//}
				if (i_1 == Amat.j[i_2]) {
					val75[i_2 - 1] = diag[0][i_1];
					col_ind75[i_2 - 1] = i_1 - 1;
				}
				else {
					val75[i_2 - 1] = Amat.aij[i_2];
					col_ind75[i_2 - 1] = Amat.j[i_2] - 1;
				}
			}
			row_ptr75[i_1 - 1] = row_ptr_start[i_1] - 1;
		}

		row_ptr75[n_a[0]] = row_ptr_end[n_a[0]];


#else

#if MY_GPU_MATH_LIBRARY_CU_ON

		if (flag != nullptr) {
			free(flag);
			flag = nullptr;
		}
		flag = my_declaration_array<bool>(n_a[0]+1, false, "flag");

#pragma omp parallel for
		for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {
			flag[i_1] = false;
		}
		calculate_row_ptr_cuda(0, nnz_a[0]-1,
			row_ptr75,
			flag, Amat_copy);

//#pragma omp parallel for
	//	for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {
		//	flag[i_1] = false;
		//}

		if (flag != nullptr) {
			free(flag);
			flag = nullptr;
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < nnz_a[0]; ++i_1) {
			val75[i_1] = Amat_copy.aij[i_1];
			col_ind75[i_1] = Amat_copy.j[i_1] - 1;
			//if (i_1 < n_a[0] + 1) {
				//row_ptr75[i_1]--;
			//}
		}

		if (Amat_copy.aij != nullptr) {
			delete[] Amat_copy.aij;
			Amat_copy.aij = nullptr;
		}
		if (Amat_copy.i != nullptr) {
			delete[] Amat_copy.i;
			Amat_copy.i = nullptr;
		}
		if (Amat_copy.j != nullptr) {
			delete[] Amat_copy.j;
			Amat_copy.j = nullptr;
		}


#else

		// Преобразование к формату CRS.

#pragma omp parallel for
		for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {

			for (integer i_2 = row_ptr_start[i_1]; i_2 <= row_ptr_end[i_1]; ++i_2) {
				//if (Amat.i[i_2] == Amat.j[i_2]) {
					//if (i_1 != Amat.i[i_2]) {
						//std::cout << "err i!=i"<< std::endl;
						//system("PAUSE");
					//}
				if (i_1 == Amat.j[i_2]) {
					val75[i_2 - 1] = diag[0][i_1];
					col_ind75[i_2 - 1] = i_1 - 1;
				}
				else {
					val75[i_2 - 1] = Amat.aij[i_2];
					col_ind75[i_2 - 1] = Amat.j[i_2] - 1;
				}
				//std::cout << val75[i_2 - 1] << "  ";
			}
			//std::cout << "\n";

			//for (integer i_2 = row_ptr_start[i_1]; i_2 <= row_ptr_end[i_1]; ++i_2) {
				//std::cout << Amat.aij[i_2] << "  ";
			//}
			//std::cout << "\n";
			//system("pause");
			row_ptr75[i_1 - 1] = row_ptr_start[i_1] - 1;
		}

		row_ptr75[n_a[0]] = row_ptr_end[n_a[0]];

		/*
		// отладка. Причина использование миксовой точности в задачах на излучение. Нужно использовать только double.
		doublereal s231 = 0.0;
		for (integer i_1 = 0; i_1 < n_a[0]; ++i_1) {
			doublereal s232 = b[i_1+1];
			for (integer i_2 = row_ptr75[i_1]; i_2 <= row_ptr75[i_1+1]-1; ++i_2) {
				s232 -= val75[i_2] * x[col_ind75[i_2] + 1];
				
			}
			s232 = s232 * s232;
			s231 += s232;
		}
		std::cout << "aprior " <<  sqrt(s231) << std::endl;
		s231 = 0.0;
		for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {
			doublereal s232 = b[i_1];
			for (integer i_2 = row_ptr_start[i_1]; i_2 <= row_ptr_end[i_1]; ++i_2) {
				if (Amat.j[i_2] == i_1) {
					s232 -= (1.0/Amat.aij[i_2]) * x[Amat.j[i_2]];
				}
				else {
					s232 -= Amat.aij[i_2] * x[Amat.j[i_2]];
				}
			}
			s232 = s232 * s232;
			s231 += s232;
		}
		std::cout << "aprior 2 " <<  sqrt(s231)  << std::endl;
		std::cout << "x = " << norma(x, n_a[0]) << "  b = " << norma(b, n_a[0]) << std::endl;
		system("pause");
		*/

#endif
#endif

		//bool bnorelax = true; // Для уравнения теплопроводности не используется релаксация.
		integer m_restart = my_amg_manager.m_restart;

		doublereal resid;
		integer  j = 1, k;
		//Vector s(m + 1), cs(m + 1), sn(m + 1), w;
		doublereal* w = new doublereal[n75];
		doublereal* s = new doublereal[m_restart + 2];
		doublereal* cs = new doublereal[m_restart + 2];
		doublereal* sn = new doublereal[m_restart + 2];

		doublereal* dx = new doublereal[n75];
		doublereal* buffer = new doublereal[n75];
		doublereal* Zcopy = new doublereal[n75 + 1];
		doublereal* vCopy = new doublereal[n75 + 1];

		// A*x=b, x - решение, b - правая часть. 
		// Индексация в х и b начинается с единицы.

		// начальное приближение
		#pragma omp parallel for
		for (integer i = 0; i < n75; ++i) dx[i] = x[i + 1];


		//doublereal normb = norm(M.solve(b));
		doublereal normb = 0.0;
		// здесь реализованы все три нормы
		// вообще говоря они все эквивалентны



		normb = NormaV_for_gmres(&b[1], n75);
		//normb = NormaV(buffer, n75);

		//Vector r = &b[1] - A * x;
		doublereal* r = new doublereal[n75];


		bool bfreecuda = false;


#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON
		MatrixCRSByVector(val75, col_ind75, row_ptr75, dx, r, n75); // результат занесён в  r
#else

#if MY_GPU_MATH_LIBRARY_CU_ON

		//30.07.2021

		//format sparse matrix;  time
		// CRS 3min 17s 610ms
		// ELL 3min 21s 310ms
		// CPU CRS 3min 20s 10ms
		// HYB 2ELL: 3min 16s 10ms




		if (init_b_first_device) {

			int device;
			device = idevice_Tesla;

			// Choose which GPU to run on, change this on a multi-GPU system.
			cudaStatus = cudaSetDevice(device);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
				system("PAUSE");
				exit(1);
			}

			cudaDeviceProp prop;
			cudaGetDeviceProperties(&prop, device);
			std::cout <<  prop.name << std::endl;
			// GeForce 840M 384потока 1ГГц каждый март 2014 года. 28нм.

			init_b_first_device = false;
		}


#if EllFormat
		doublereal* data = nullptr;
		integer* indices = nullptr;
		integer string_size = -1;

		doublereal* data_bound = nullptr;
		integer* indices_bound = nullptr;
		integer string_size_bound = 2;

		for (integer i77 = 0; i77 < n75; ++i77) {
			if (row_ptr75[i77 + 1] - row_ptr75[i77] > string_size) {
				string_size = row_ptr75[i77 + 1] - row_ptr75[i77];
			}
			//if (i77>= maxelm) {
				//std::cout << "\nstart hvost=" << i77 << "  maxelm=" << maxelm << "start "<< row_ptr75[i77]  << "end=" << row_ptr75[i77 + 1] << std::endl;
				//system("pause");
			//}
		}

		data = new doublereal[string_size * maxelm_out];
		indices = new integer[string_size * maxelm_out];

		data_bound = new doublereal[string_size_bound * (n- maxelm_out)];
		indices_bound = new integer[string_size_bound * (n- maxelm_out)];


		// Перепаковка в Ell формат.
#pragma omp parallel for
		for (integer i77 = 0; i77 < string_size * maxelm_out; ++i77) {
			data[i77] = 0.0;
		}

#pragma omp parallel for
		for (integer i77 = 0; i77 < string_size_bound * (n - maxelm_out); ++i77) {
			data_bound[i77] = 0.0;
		}

#pragma omp parallel for
		for (integer i77 = 0; i77 < maxelm_out; ++i77) {

			const integer rowend = row_ptr75[i77 + 1];
			const integer rowbeg = row_ptr75[i77];

			integer jstart = 0;
			for (integer j77 = rowbeg; j77 < rowend; ++j77)
			{
				const integer element_offset = i77 + jstart * maxelm_out;
				data[element_offset] = val75[j77];
				indices[element_offset] = col_ind75[j77];

				++jstart;
			}
			for (integer j77 = jstart; j77 < string_size; ++j77)
			{
				const integer element_offset = i77 + j77 * maxelm_out;

				if (indices[element_offset - 1] < n75 - 2) {
					indices[element_offset] = indices[element_offset - 1] + 1;
					data[element_offset] = 0.0;
				}
				else {
					indices[element_offset] = indices[element_offset - 1];
					data[element_offset] = 0.0;
				}

			}
		}

#pragma omp parallel for
		for (integer i77 = maxelm_out; i77 < n; ++i77) {

			const integer rowend = row_ptr75[i77 + 1];
			const integer rowbeg = row_ptr75[i77];

			integer jstart = 0;
			for (integer j77 = rowbeg; j77 < rowend; ++j77)
			{
				const integer element_offset = i77 - maxelm_out + jstart * (n - maxelm_out);
				data_bound[element_offset] = val75[j77];
				indices_bound[element_offset] = col_ind75[j77];

				++jstart;
			}
			for (integer j77 = jstart; j77 < string_size_bound; ++j77)
			{
				const integer element_offset = i77 - maxelm_out + j77 * (n - maxelm_out);

				if (indices_bound[element_offset - 1] < n75 - 2) {
					indices_bound[element_offset] = indices_bound[element_offset - 1] + 1;
					data_bound[element_offset] = 0.0;
				}
				else {
					indices_bound[element_offset] = indices_bound[element_offset - 1];
					data_bound[element_offset] = 0.0;
				}

			}
		}


		// Освобождение оперативной памяти.
		if (val75 != nullptr) {
			delete[] val75;
			val75 = nullptr;
		}
		if (col_ind75 != nullptr) {
			delete[] col_ind75;
			col_ind75 = nullptr;
		}
		if (row_ptr75 != nullptr) {
			delete[] row_ptr75;
			row_ptr75 = nullptr;
		}


		doublereal* dev_data = nullptr;
		integer* dev_indices = nullptr;

		doublereal* dev_data_bound = nullptr;
		integer* dev_indices_bound = nullptr;

#else
		doublereal* dev_val = nullptr;
		integer* dev_col_ind = nullptr;
		integer* dev_row_ptr = nullptr;
#endif


		doublereal* dev_V = nullptr;
		doublereal* dev_tmp = nullptr;

#if EllFormat
		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_data, (maxelm_out * string_size) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_val failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_indices, (maxelm_out * string_size) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_col_ind failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_data_bound, ((n- maxelm_out) * string_size_bound) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_val failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_indices_bound, ((n- maxelm_out) * string_size_bound) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_col_ind failed!");
			system("PAUSE");
			exit(1);
		}

#else

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_val, (nnz) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_val failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_col_ind, (nnz) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_col_ind failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_row_ptr, (n75 + 1) * sizeof(integer));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_row_ptr failed!");
			system("PAUSE");
			exit(1);
		}

#endif

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_V, (n75) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_V failed!");
			system("PAUSE");
			exit(1);
		}

		// Allocate GPU buffers for three vectors (two input, one output)    .
		cudaStatus = cudaMalloc((void**)&dev_tmp, (n75) * sizeof(doublereal));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc dev_tmp failed!");
			system("PAUSE");
			exit(1);
		}

#if EllFormat

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_data, data, maxelm_out * string_size * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_val, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_indices, indices, maxelm_out * string_size * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_col_ind, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		delete[] data;
		data = nullptr;
		delete[] indices;
		indices = nullptr;

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_data_bound, data_bound, (n- maxelm_out) * string_size_bound * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_data_bound, data_bound HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_indices_bound, indices_bound, (n- maxelm_out) * string_size_bound * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_indices_bound, indices_bound HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		delete[] data_bound;
		data_bound = nullptr;
		delete[] indices_bound;
		indices_bound = nullptr;

#else

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_val, val75, nnz * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_val, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_col_ind, col_ind75, nnz * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_col_ind, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}


		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_row_ptr, row_ptr75, (n75 + 1) * sizeof(integer), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_rho_ptr, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Освобождение оперативной памяти.
		if (val75 != nullptr) {
			delete[] val75;
			val75 = nullptr;
		}
		if (col_ind75 != nullptr) {
			delete[] col_ind75;
			col_ind75 = nullptr;
		}
		if (row_ptr75 != nullptr) {
			delete[] row_ptr75;
			row_ptr75 = nullptr;
		}

#endif


		// Copy input vectors from host memory to GPU buffers.
		cudaStatus = cudaMemcpy(dev_V, dx, n75 * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_V, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}

		// Copy input vectors from host memory to GPU buffers.
		/*cudaStatus = cudaMemcpy(dev_tmp, dax75, n75 * sizeof(doublereal), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy dev_tmp, val HostToDevice failed!");
			system("PAUSE");
			exit(1);
		}*/

#if EllFormat

		MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data, dev_indices, string_size, dev_V, dev_tmp, maxelm_out, 0);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

		MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data_bound, dev_indices_bound, string_size_bound, dev_V, dev_tmp, n- maxelm_out, maxelm_out);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}

#else

		MatrixCRSByVectorKernel << <128, 128 >> > (dev_val, dev_col_ind, dev_row_ptr, dev_V, dev_tmp, n75);
		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "MatrixCRSByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			system("PAUSE");
			exit(1);
		}


#endif

		cudaStatus = cudaMemcpy(r, dev_tmp, n75 * sizeof(doublereal), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy tmp, dev_tmp DeviceToHost failed!");
			system("PAUSE");
			exit(1);
		}



#else

		MatrixCRSByVector(val75, col_ind75, row_ptr75, dx, r, n75); // результат занесён в  r


#endif
#endif

		

#pragma omp parallel for 
		for (integer i = 0; i < n75; ++i) r[i] = b[i + 1] - r[i];

		//  calculate residual precontidioning;


		//doublereal beta = norm(r);
		doublereal beta = 0.0;



		beta = NormaV_for_gmres(r, n75);

		if (fabs(normb) < 1.0e-30)
			normb = 1;

		doublereal norm_r = beta;


		integer maxit = 2000;

		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			// Нормальная система СЛАУ должна сходится не более чем за 2000итераций.
			maxit = 2000;
		}

		resid = norm_r / normb;

		if (0&&((resid) <= dterminatedTResudual)) {
			//tol = resid;
			//maxit = 0;
			//return 0;
			delete[] val75;
			delete[] col_ind75;
			delete[] row_ptr75;
			delete[] w;
			delete[] s;
			delete[] cs;
			delete[] sn;
			delete[] dx;
			delete[] buffer;
			delete[] Zcopy;
			delete[] vCopy;
			delete[] r;
			goto LABEL_FGMRES_CONTINUE;
		}

		//integer i_1 = 0; // счётчик цикла for

		doublereal** H = new doublereal * [m_restart + 2]; // Hessenberg
		for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) H[i_1] = new doublereal[m_restart + 2];

#pragma omp parallel for 
		for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1)
		{
			for (integer j_1 = 0; j_1 < m_restart + 2; ++j_1)
			{
				H[i_1][j_1] = 0.0;
			}
		}

		//Vector *v = new Vector[m_restart + 1];
		doublereal** v = new doublereal * [m_restart + 2];
		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) v[i_1] = new doublereal[n75];

#pragma omp parallel for
		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
			for (integer j_1 = 0; j_1 < n75; ++j_1)
			{
				v[i_1][j_1] = 0.0;
			}
		}

		doublereal** Z = new doublereal * [m_restart + 2];
		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) Z[i_1] = new doublereal[n75];

#pragma omp parallel for
		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
			for (integer j_1 = 0; j_1 < n75; ++j_1)
			{
				Z[i_1][j_1] = 0.0;
			}
		}

		j = 1; // номер первой итерации
			   //doublereal delta = 1.0e-3;// DOPOLNENIE

		//integer i_copy;

		//if ((bnonlinear) && (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY)) {
			//maxit = 5; // 5 не подходит, надо решать до сходимости.
		//}

		while (j <= maxit) {

			//v[0] = r * (1.0 / beta);    // ??? r / beta
#pragma omp parallel for firstprivate(beta)
			for (integer j_1 = 0; j_1 < n75; ++j_1)
			{
				v[0][j_1] = r[j_1] * (1.0 / beta);
			}

			//s = 0.0;
#pragma omp parallel for
			for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
				s[i_1] = 0.0;
			}
			s[0] = beta;
			//s[0] = 1.0;

#pragma omp parallel for
			for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1)
			{ // DOPOLNENIE
				for (integer j_1 = 0; j_1 < m_restart + 2; ++j_1)
				{
					H[i_1][j_1] = 0.0;
				}
			}

			integer i_global = 0;

			// Ортогонализация Арнольди.
			for (i_global = 0; i_global < m_restart && j <= maxit; ++i_global, ++j) {

				//i_copy = i_global;

				integer i = i_global;// i - локальная копия.

				// KZ[i]=v[i]

				// (LU)Z[i]=v[i];

				// multigrid Ruge and Stuben preconditioning [1986].
				// достаточно одного V цикла.
				// K*Z = v;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < n75; ++i_1) {
					Zcopy[i_1 + 1] = 0.0; // Только ноль единственно верно.
					//Zcopy[i_1 + 1] = urd(gen);
					//Zcopy[i_1 + 1] = 0.001 * (rand() % (1001));
					vCopy[i_1 + 1] = v[i][i_1];
				}

				
					// Предобуславливание с помощью V цикла многосеточного метода.
					// Нулевое начальное приближение
					for (integer i_numberV_cycle = 0; i_numberV_cycle < 1; ++i_numberV_cycle) {
						// достаточно одного V цикла.
						// A*Zcopy=vCopy;
						// В Zcopy и vCopy нумерация начинается с единицы.
						V_cycle_solve<doublerealT>(Amat, Zcopy, vCopy, process_flow_logic, row_ptr_start,
							row_ptr_end, residual_fine, diag, diag_minus_one, n_a, bonly_serial,
							process_flow_beta, F_false_C_true, nu1, nu2, bILU2smoother,
							ilevel, inumberVcyclelocbicgstab, imyinit, maxlevel, milu2, milu0, nested_desection,
							P, nnz_aRP, residual_coarse, igam, nnz_a,
							error_approx_coarse, dapply_ilu_max_pattern_size,
							process_flow_alpha,
							error_approx_fine, nFinestSweeps, ibsp_length, iVar);
						//system("pause");
					}
				

#pragma omp parallel for
				for (integer i_1 = 0; i_1 < n75; ++i_1) {
					Z[i][i_1] = Zcopy[i_1 + 1];
				}


				// Совсем без предобуславливателя.
				//for (i_1 = 0; i_1 < n75; ++i_1) Z[i][i_1] = v[i][i_1];

				// Закомментировано без предобуславливания.
				//w = A * Z[i];

#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON
				MatrixCRSByVector(val75, col_ind75, row_ptr75, Z[i], w, n75); // результат занесён в  w
#else

#if MY_GPU_MATH_LIBRARY_CU_ON

// Copy input vectors from host memory to GPU buffers.
				cudaStatus = cudaMemcpy(dev_V, Z[i], n75 * sizeof(doublereal), cudaMemcpyHostToDevice);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy dev_tmp, val HostToDevice failed!");
					system("PAUSE");
					exit(1);
				}

#if EllFormat

				MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data, dev_indices, string_size, dev_V, dev_tmp, maxelm_out, 0);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}

				MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data_bound, dev_indices_bound, string_size_bound, dev_V, dev_tmp, n- maxelm_out, maxelm_out);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}

#else

				MatrixCRSByVectorKernel << <128, 128 >> > (dev_val, dev_col_ind, dev_row_ptr, dev_V, dev_tmp, n75);
				// Check for any errors launching the kernel
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "MatrixCRSByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
					system("PAUSE");
					exit(1);
				}

#endif

				cudaStatus = cudaMemcpy(w, dev_tmp, n75 * sizeof(doublereal), cudaMemcpyDeviceToHost);
				if (cudaStatus != cudaSuccess) {
					fprintf(stderr, "cudaMemcpy tmp, dev_tmp DeviceToHost failed!");
					system("PAUSE");
					exit(1);
				}

#else

				MatrixCRSByVector(val75, col_ind75, row_ptr75, Z[i], w, n75); // результат занесён в  w

#endif
#endif
				

				for (k = 0; k <= i; ++k) {
					doublereal hscal= Scal(w, v[k], n75);
					H[k][i] = hscal;

					// если раскоментировать приложение падает
#pragma omp parallel for firstprivate(hscal, n75, k)
					for (integer j_1 = 0; j_1 < n75; ++j_1)
					{
						w[j_1] -= hscal * v[k][j_1];
					}
				}
				H[i + 1][i] = NormaV_for_gmres(w, n75);


				doublereal mult_for_w = (1.0 / H[i + 1][i]);

				// если раскоментировать приложение падает
#pragma omp parallel for firstprivate(n75, i, mult_for_w)
				for (integer j_1 = 0; j_1 < n75; ++j_1)
				{
					v[i + 1][j_1] = w[j_1] * mult_for_w; // ??? w / H(i+1, i)
				}
				// Окончание ортогонализации Арнольди.
				// В v - хранится ортонормированный базис подпространства Крылова размерности m_restart.
				// H - Верхнетреугольная матрица Хессенберга - матрица коэффициентов ортогонализации.


				// 26.11.2017
				// Это проверенный и испытанный кусок кода.
				for (k = 0; k < i; ++k) {

					ApplyPlaneRotation(H[k][i], H[k + 1][i], cs[k], sn[k]);
				}

				GeneratePlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
#pragma omp parallel sections
				{// распараллеливание на два потока.
#pragma omp section
					{
						ApplyPlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
					}
#pragma omp section
					{
						ApplyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);
					}
				}// omp parallel sections

				


				// Вручную устраняем случай полного совпадения невязок на двух соседних итерациях,
				// т.к. иначе это приводит к развалу решения.
				//if (fabs(s[i] - s[i + 1]) < 1.0e-37) s[i + 1] = 1.05*s[i];

				if (bprint_mesage_diagnostic) {
					std::cout << j << " " << fabs(s[i + 1]) / normb << " \n";
					//std::cout << j <<" "<< beta*fabs(s[i + 1]) << " \n";
					//system("pause");
				}

				resid = fabs(s[i + 1]) / normb;
				//resid = beta*fabs(s[i + 1]);

				bool b_dosrochnji_vjhod_Mechanical = false;
				if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
				{
					b_dosrochnji_vjhod_Mechanical = true;
				}

				if (((j + 2 == maxit) && (b_dosrochnji_vjhod_Mechanical)) || (((((i > 4)&&((iVar==VELOCITY_X_COMPONENT)|| (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)))||
					((iVar != VELOCITY_X_COMPONENT) && (iVar != VELOCITY_Y_COMPONENT) && (iVar != VELOCITY_Z_COMPONENT)))||
					((iVar==TEMP)&&(((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY)&&(i > 4))||
						(steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY))))
					&&((resid) < dterminatedTResudual))) {


					//std::cout << "j==" << j;
			//if (b_dosrochnji_vjhod_Mechanical) {
				//std::cout << "  b_dosrochnji_vjhod_Mechanical";
			//}

					if (bprint_mesage_diagnostic)
					{
						std::cout << "dosrochnji vjhod" << std::endl;
						std::cout << "resid=" << resid << " dterminatedTResudual=" << dterminatedTResudual;
						//system("pause");				
					}
					Update(dx, i, n75, H, s, Z);
					//tol = resid;
					//maxit = j;

#pragma omp parallel for
					for (integer i_1 = 0; i_1 < n75; ++i_1) {
						x[i_1 + 1] = dx[i_1];
						x_best_search[i_1 + 1] = static_cast<doublerealT>(dx[i_1]);
					}

					for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] v[i_1];
					delete[] v;
					for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] Z[i_1];
					delete[] Z;
					for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) delete[] H[i_1];
					delete[] H;
					delete[] dx;
					delete[] buffer;
					delete[] r;
					delete[] w;
					delete[] s;
					delete[] cs;
					delete[] sn;
					delete[] Zcopy;
					delete[] vCopy;

					// Освобождение оперативной памяти.
					if (val75 != nullptr) {
						delete[] val75;
						val75 = nullptr;
					}
					if (col_ind75 != nullptr) {
						delete[] col_ind75;
						col_ind75 = nullptr;
					}
					if (row_ptr75 != nullptr) {
						delete[] row_ptr75;
						row_ptr75 = nullptr;
					}

					if (!bfreecuda) {

#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON

#else
#if MY_GPU_MATH_LIBRARY_CU_ON

#if EllFormat

						cudaFree(dev_data);
						cudaFree(dev_indices);

						delete[] data;
						delete[] indices;

						cudaFree(dev_data_bound);
						cudaFree(dev_indices_bound);

						delete[] data_bound;
						delete[] indices_bound;

#else


						cudaFree(dev_row_ptr);
						cudaFree(dev_val);
						cudaFree(dev_col_ind);

#endif

						cudaFree(dev_V);
						cudaFree(dev_tmp);

#endif
#endif
						bfreecuda = true;
					}


					goto LABEL_FGMRES_CONTINUE;

				}
			}



			// i_global-1 -> m_restart-1
			Update(dx, i_global - 1, n75, H, s, Z);//i_global-1 

											//r = M.solve(b - A * x);

#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON

			MatrixCRSByVector(val75, col_ind75, row_ptr75, dx, r, n75); // Результат занесён в r

#else

#if MY_GPU_MATH_LIBRARY_CU_ON

// Copy input vectors from host memory to GPU buffers.
			cudaStatus = cudaMemcpy(dev_V, dx, n75 * sizeof(doublereal), cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy dev_tmp, val HostToDevice failed!");
				system("PAUSE");
				exit(1);
			}

#if EllFormat

			MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data, dev_indices, string_size, dev_V, dev_tmp, maxelm_out, 0);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

			MatrixEllPackItpackByVectorKernel << <128, 128 >> > (dev_data_bound, dev_indices_bound, string_size_bound, dev_V, dev_tmp, n- maxelm_out, maxelm_out);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, " MatrixEllPackItpackByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

#else

			MatrixCRSByVectorKernel << <128, 128 >> > (dev_val, dev_col_ind, dev_row_ptr, dev_V, dev_tmp, n75);
			// Check for any errors launching the kernel
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "MatrixCRSByVectorKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
				system("PAUSE");
				exit(1);
			}

#endif

			cudaStatus = cudaMemcpy(r, dev_tmp, n75 * sizeof(doublereal), cudaMemcpyDeviceToHost);
			if (cudaStatus != cudaSuccess) {
				fprintf(stderr, "cudaMemcpy tmp, dev_tmp DeviceToHost failed!");
				system("PAUSE");
				exit(1);
			}

#else

			MatrixCRSByVector(val75, col_ind75, row_ptr75, dx, r, n75); // Результат занесён в r

#endif
#endif

			
#pragma omp parallel for
			for (integer i_1 = 0; i_1 < n75; ++i_1) {
				r[i_1] = b[i_1 + 1] - r[i_1];
			}

			//beta = norm(r);
			beta = NormaV_for_gmres(r, n75);

			resid = beta / normb;
			//resid = beta;

			

			

			if ((resid) < dterminatedTResudual) {

				// В случае механической задачи сойтись полностью бывает проблематично, но нужно вернуть хоть какое-то решение
				// механической задачи, тем более что метод FGMRes более стабилен чем метод BiCGStab.

				//tol = resid;
				//maxit = j;

				std::cout << "end" << std::endl;
				//system("pause");

#pragma omp parallel for
				for (integer i_1 = 0; i_1 < n75; ++i_1) {
					x[i_1 + 1] = dx[i_1];
					x_best_search[i_1 + 1] = static_cast<doublerealT>(dx[i_1]);
				}

				for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] v[i_1];
				delete[] v;
				for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] Z[i_1];
				delete[] Z;
				for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) delete[] H[i_1];
				delete[] H;
				delete[] dx;
				delete[] buffer;
				delete[] r;
				delete[] w;
				delete[] s;
				delete[] cs;
				delete[] sn;
				delete[] Zcopy;
				delete[] vCopy;

				// Освобождение оперативной памяти.
				if (val75 != nullptr) {
					delete[] val75;
					val75 = nullptr;
				}
				if (col_ind75 != nullptr) {
					delete[] col_ind75;
					col_ind75 = nullptr;
				}
				if (row_ptr75 != nullptr) {
					delete[] row_ptr75;
					row_ptr75 = nullptr;
				}


				if (!bfreecuda) {

#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON

#else
#if MY_GPU_MATH_LIBRARY_CU_ON

#if EllFormat

					cudaFree(dev_data);
					cudaFree(dev_indices);

					delete[] data;
					delete[] indices;

					cudaFree(dev_data_bound);
					cudaFree(dev_indices_bound);

					delete[] data_bound;
					delete[] indices_bound;

#else


					cudaFree(dev_row_ptr);
					cudaFree(dev_val);
					cudaFree(dev_col_ind);

#endif

					cudaFree(dev_V);
					cudaFree(dev_tmp);

#endif
#endif
					bfreecuda = true;
				}

				goto LABEL_FGMRES_CONTINUE;


			}
		}

		if (j-1 == maxit) {
			// amg solver divergence detected 
			// Превышено допустимое число итераций.
			return true;
		}


		//tol = resid;
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < n75; ++i_1) {
			x[i_1 + 1] = dx[i_1];
			x_best_search[i_1 + 1] = static_cast<doublerealT>(dx[i_1]);
		}

		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] v[i_1];
		delete[] v;
		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] Z[i_1];
		delete[] Z;
		for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) delete[] H[i_1];
		delete[] H;
		delete[] dx;
		delete[] buffer;
		delete[] r;
		delete[] w;
		delete[] s;
		delete[] cs;
		delete[] sn;
		delete[] Zcopy;
		delete[] vCopy;

		// Освобождение оперативной памяти.
		if (val75 != nullptr) {
			delete[] val75;
			val75 = nullptr;
		}
		if (col_ind75 != nullptr) {
			delete[] col_ind75;
			col_ind75 = nullptr;
		}
		if (row_ptr75 != nullptr) {
			delete[] row_ptr75;
			row_ptr75 = nullptr;
		}

		if (!bfreecuda) {

#if MY_GPU_MATH_LIBRARY_CU_REALLOC_ON

#else
#if MY_GPU_MATH_LIBRARY_CU_ON

#if EllFormat

			cudaFree(dev_data);
			cudaFree(dev_indices);

			delete[] data;
			delete[] indices;

			cudaFree(dev_data_bound);
			cudaFree(dev_indices_bound);

			delete[] data_bound;
			delete[] indices_bound;

#else


			cudaFree(dev_row_ptr);
			cudaFree(dev_val);
			cudaFree(dev_col_ind);

#endif

			cudaFree(dev_V);
			cudaFree(dev_tmp);

#endif
#endif
			bfreecuda = true;
		}


		goto LABEL_FGMRES_CONTINUE;




	}
	else if (((iVar == TOTALDEFORMATIONVAR) || (iVar ==  PAM)) && (my_amg_manager.istabilization == 3)) {
	     // Метод сопряженных градиентов Хестенса и Штифеля 1952 года
		 // для симметричной положительно определенной механической системы.

		 integer inumberVcyclelocbicgstab = 1;

		 // 02.02.2022 - Работает в единственной матрице. Экономим оперативную память.

		 // нумерация векторов начинается с нуля.
		 integer n75 = n_a[0]; // число неизвестных на подробном уровне.
		 //doublereal* val75 = nullptr;
		 //val75 = new doublereal[nnz_a[0]];
		 //integer* col_ind75 = nullptr;
		 //col_ind75 = new integer[nnz_a[0]];
		 //integer* row_ptr75 = nullptr;
		 //row_ptr75 = new integer[n_a[0] + 1];
		 //if ((val75 == nullptr) || (col_ind75 == nullptr) || (row_ptr75 == nullptr)) {
			 // недостаточно памяти на данном оборудовании.
			 //std::cout << "Problem: not enough memory on your equipment for val, col_ind or row_ptr: bicgStab + camg..." << std::endl;
			 //std::cout << "Please any key to exit..." << std::endl;
			 //exit(1);
		 //}


		 // Преобразование к формату CRS.

//#pragma omp parallel for
		 //for (integer i_1 = 1; i_1 <= n_a[0]; ++i_1) {

			 //for (integer i_2 = row_ptr_start[i_1]; i_2 <= row_ptr_end[i_1]; ++i_2) {
				 //C//if (Amat.i[i_2] == Amat.j[i_2]) {
					 //ifC// (i_1 != Amat.i[i_2]) {
						 //C//std::cout << "err i!=i"<< std::endl;
						 //C//system("PAUSE");
					 //C//}
				 //if (i_1 == Amat.j[i_2]) {
					// val75[i_2 - 1] = diag[0][i_1];
					// col_ind75[i_2 - 1] = i_1 - 1;
				 //}
				 //else {
					 //val75[i_2 - 1] = Amat.aij[i_2];
					 //col_ind75[i_2 - 1] = Amat.j[i_2] - 1;
				 //}
			 //}
			 //row_ptr75[i_1 - 1] = row_ptr_start[i_1] - 1;
		 //}

		 //row_ptr75[n_a[0]] = row_ptr_end[n_a[0]];

		 // Вектора необходимые для работы предобусловленного метода сопряжённых градиентов.
		 doublereal* r75 = nullptr;
		 doublereal* p75 = nullptr;
		 doublereal* dx75 = nullptr;
		 doublereal* dax75 = nullptr;
		 doublereal* z76 = nullptr;
		 doublereal* s76 = nullptr;
		 doublereal* zold = nullptr;
		 doublereal* znew = nullptr;
		 //doublereal* rold75 = nullptr;

		 r75 = new doublereal[n75];
		 p75 = new doublereal[n75];
		 dx75 = new doublereal[n75];
		 dax75 = new doublereal[n75];
		 z76 = new doublereal[n75 + 1];
		 s76 = new doublereal[n75 + 1];
		 zold = new doublereal[n75];
		 znew = new doublereal[n75];
		 //rold75 = new doublereal[n75];

		 // initialize
#pragma omp parallel for 
		 for (integer i75 = 0; i75 < n75; ++i75) {
			 r75[i75] = 0.0;
			 p75[i75] = 0.0;

			 z76[i75] = 0.0;
			 s76[i75] = 0.0;
			
			 zold[i75] = 0.0;
			 znew[i75] = 0.0;

			 // результат умножения матрицы на вектор.
			 dax75[i75] = 0.0;
			 // Начальное приближение.
			 dx75[i75] = x[i75 + 1];
		 }
		 z76[n75] = 0.0;
		 s76[n75] = 0.0;

		 // Умножение матрицы на вектор. Нумерации векторов начинаются с нуля.
		 //MatrixCRSByVector(val75, col_ind75, row_ptr75, dx75, dax75, n75); // результат занесён в  dax75

		 // Вычисление ri75 и roc75.
//#pragma omp parallel for
	//	 for (integer i75 = 0; i75 < n75; ++i75) {
		//	 r75[i75] = b[i75 + 1] - dax75[i75];
		 //}

		 residual_for_Cheb(Amat, 1, n_a[0], b, x, r75, row_ptr_start, row_ptr_end, 0, 0);

		 // Предобуславливание.
		 
#pragma omp parallel for 
		 for (integer i75 = 0; i75 < n75; ++i75) {
			 
			 z76[i75 + 1] = 0.0;
			 s76[i75 + 1] = r75[i75];
			// r75[i75] = s76[i75 + 1];
		 }

		 int iamg_cg = 0;

		

		 if (iamg_cg==0) {

			 // multigrid RUMBA preconditioner
			 // Вставлено 6.01.2017 begin
			 // одного V цикла недостаточно.
			 // A*z76=s76;
			 V_cycle_solve<doublerealT>(Amat, z76, s76, process_flow_logic, row_ptr_start,
				 row_ptr_end, residual_fine, diag, diag_minus_one, n_a, bonly_serial,
				 process_flow_beta, F_false_C_true, nu1, nu2, bILU2smoother,
				 ilevel, inumberVcyclelocbicgstab, imyinit, maxlevel, milu2, milu0, nested_desection,
				 P, nnz_aRP, residual_coarse, igam, nnz_a,
				 error_approx_coarse, dapply_ilu_max_pattern_size,
				 process_flow_alpha,
				 error_approx_fine, nFinestSweeps, ibsp_length, iVar);
			 // Вставлено 6.01.2017 end

		 }
		 else {

			 seidel_for_cg(Amat, row_ptr_start, row_ptr_end, z76, s76,  n_a[0]);
		 }


#pragma omp parallel for 
		 for (integer i75 = 0; i75 < n75; ++i75) {			 
			 // Возвращаем результат.
			 zold[i75] = z76[i75 + 1];
			 p75[i75] = zold[i75];
		 }

		 doublereal alpha75, beta75;

		 doublereal epsilon75 = dterminatedTResudual;  // точность вычисления
		 if (iVar == TEMP) {
			 epsilon75 *= 1.0e-4; // 1.0e-4
		 }
		 if (iVar == TOTALDEFORMATIONVAR) {
			 //epsilon75 *= 1.0e-4; // 1.0e-4
			 epsilon75 *= 2.0e-2;
			  //epsilon75 *= 1.0e-12;
		 }

		 doublereal dnew = Scal(r75, zold, n75);

		 delete[] zold;

		 // Если число расходимостей превысит оговорённую константу то произойдёт выход из алгоритма.
		 integer i_signal_break_pam_opening75 = 0;
		 // x хорошее значение.
		 const integer i_limit_signal_pam_break_opening75 = 8000;//20
		 doublereal delta_old_iter75 = 1.0e10;

		 int k53 = 0;

		 // Метод cg расходится для поправки давления. Пока по непонятной причине. 03.02.2022.
		 do {

			 ++k53;

			

			 // Умножение матрицы на вектор. Нумерации векторов начинаются с нуля.
			 //MatrixCRSByVector(val75, col_ind75, row_ptr75, p75, dax75, n75); // результат занесён в  dax75
			 	
			 // нумерация  p75 и dax начинается с нуля.
			 Matrix_by_vector_for_Cheb(Amat, 1, n_a[0], p75, dax75, row_ptr_start, row_ptr_end, 0, 0);// результат занесён в  dax75

			 doublereal scal_val1 = Scal(p75, dax75, n75);

			 if ((fabs(dnew) < 1e-30) && (fabs(scal_val1) < 1e-30)) {
				 alpha75 = 1.0;
			 }
			 else if (fabs(dnew) < 1e-30) {
				 alpha75 = 0.0;
			 }
			 else {
				 alpha75 = dnew / scal_val1;
			 }

			 

#pragma omp parallel for
			 for (integer i75 = 0; i75 < n75; ++i75)
			 {
				 dx75[i75] += alpha75 * p75[i75];
				 r75[i75] -= alpha75 * dax75[i75];
			 }


			 doublereal deltai75 = NormaV_for_gmres(r75, n75);
			 if (deltai75 < epsilon75) break;

			 if (k53 > 2600) break; // 15.02.2017

			 // 02.02.2022
			// if (fabs(delta_old_iter75 - deltai75) < 1.0e-4) break;

			 if (deltai75 > delta_old_iter75) i_signal_break_pam_opening75++;
			
			/* if (iVar == PAM) {
				 if (i_signal_break_pam_opening75 > i_limit_signal_pam_break_opening75) {
					 // досрочный выход из цикла.
					 if (bprint_mesage_diagnostic) {
						 std::cout << "icount PAM=" << k53 << "\n";
					 }

					 break;
				 }
			 }*/

			 

			 // Если номер итерации больше третьей и процесс расходится и мы решаем поправку давления то стоп.
			 // 31.10.2019
			 //if ((iVar == PAM) && (k53 > 3) && (deltai75 > delta_old_iter75)) {
				 // Естественная конвекция на АЛИС сетке в случае кода заданы только температуры а не мощности.
				 //if (inumiterSIMPLE371 < 10)
				 //{
					// break;
				 //}
			 //}

			 delta_old_iter75 = deltai75;

			 if (iVar == PAM) {
				 if (i_signal_break_pam_opening75 > i_limit_signal_pam_break_opening75) {
					 // досрочный выход из цикла.
#if doubleintprecision == 1
					 printf("icount PAM=%d\n", k53);
#else
					 printf("icount PAM=%d\n", k53);
#endif
					 if (!b_on_adaptive_local_refinement_mesh) {
						 break;
					 }
					 //system("pause");
				 }
			 }

			 // Предобуславливание.
#pragma omp parallel for 
			 for (integer i75 = 0; i75 < n75; ++i75) 
			 {
				 z76[i75 + 1] = 0.0;
				 s76[i75 + 1] = r75[i75];
			 }

			 if (iamg_cg==0) {

				 // multigrid RUMBA preconditioner
				 // Вставлено 6.01.2017 begin
				 // одного V цикла недостаточно.
				 // A*z76=s76;
				 V_cycle_solve<doublerealT>(Amat, z76, s76, process_flow_logic, row_ptr_start,
					 row_ptr_end, residual_fine, diag, diag_minus_one, n_a, bonly_serial,
					 process_flow_beta, F_false_C_true, nu1, nu2, bILU2smoother,
					 ilevel, inumberVcyclelocbicgstab, imyinit, maxlevel, milu2, milu0, nested_desection,
					 P, nnz_aRP, residual_coarse, igam, nnz_a,
					 error_approx_coarse, dapply_ilu_max_pattern_size,
					 process_flow_alpha,
					 error_approx_fine, nFinestSweeps, ibsp_length, iVar);
				 // Вставлено 6.01.2017 end
			 }
			 else {
				 seidel_for_cg(Amat, row_ptr_start, row_ptr_end, z76, s76,  n_a[0]);
			 }

			 

#pragma omp parallel for 
			 for (integer i75 = 0; i75 < n75; ++i75) {

				 // Возвращаем результат.
				 znew[i75] = z76[i75 + 1];
			 }

			 doublereal dold = dnew;
			 dnew = Scal(znew, r75, n75);

			 //beta75= scr75 / scr75_old;
			 //beta75 = Scal(r75, znew,n75) / Scal(rold75, zold, n75);

			 if ((fabs(dnew) < 1e-30) && (fabs(dold) < 1e-30)) {
				 beta75 = 1.0;
			 }
			 else if (fabs(dnew) < 1e-30) {
				 beta75 = 0.0;
			 }
			 else {
				 beta75 = dnew / dold;
			 }

			  

#pragma omp parallel for
			 for (integer i75 = 0; i75 < n75; ++i75) {
				// p75[i75] = r75[i75] + beta75 * p75[i75];
				 p75[i75] = znew[i75] + beta75 * p75[i75];
			 }

			 std::cout << k53 << "  " << NormaV(r75, n75) << std::endl;

		 } while (NormaV_for_gmres(r75, n75) >= epsilon75);

		 // Возвращение результата вычислений.
#pragma omp parallel for
		 for (integer i75 = 0; i75 < n75; ++i75) 
		 {
			 x[i75 + 1] = dx75[i75];
			 x_best_search[i75 + 1] = static_cast<doublerealT>(dx75[i75]);
		 }

		 

		 // Освобождение оперативной памяти.
		 delete[] r75;
		 delete[] p75;
		 delete[] dx75;
		 delete[] dax75;
		 delete[] znew;
		 
		 delete[] z76;
		 delete[] s76;
		 //delete[] rold75;

		 // Освобождение оперативной памяти из под матрицы СЛАУ.

		 //delete[] val75;
		 //val75 = nullptr;
		 //delete[] col_ind75;
		 //col_ind75 = nullptr;
		 //delete[] row_ptr75;
		 //row_ptr75 = nullptr;

    }
	else {
	    std::cout << "error in classic_aglomeration_amg6_2018year.cpp in function solution phase undefined solver select." << std::endl;
		system("pause");
		exit(1);
	}



	for (int i_52 = 0; i_52 < maxlevel; ++i_52) {
		if (dbgmres_smoother[i_52].val != nullptr) {
			delete[] dbgmres_smoother[i_52].val;
			dbgmres_smoother[i_52].val = nullptr;
		}
		if (dbgmres_smoother[i_52].col_ind != nullptr) {
			delete[] dbgmres_smoother[i_52].col_ind;
			dbgmres_smoother[i_52].col_ind = nullptr;
		}
		if (dbgmres_smoother[i_52].row_ptr != nullptr) {
			delete[] dbgmres_smoother[i_52].row_ptr;
			dbgmres_smoother[i_52].row_ptr = nullptr;
		}
	}
	delete[] dbgmres_smoother;
	dbgmres_smoother = nullptr;

LABEL_FGMRES_CONTINUE:

	if (dbgmres_smoother != nullptr) {
		for (int i_52 = 0; i_52 < maxlevel; ++i_52) {
			if (dbgmres_smoother[i_52].val != nullptr) {
				delete[] dbgmres_smoother[i_52].val;
				dbgmres_smoother[i_52].val = nullptr;
			}
			if (dbgmres_smoother[i_52].col_ind != nullptr) {
				delete[] dbgmres_smoother[i_52].col_ind;
				dbgmres_smoother[i_52].col_ind = nullptr;
			}
			if (dbgmres_smoother[i_52].row_ptr != nullptr) {
				delete[] dbgmres_smoother[i_52].row_ptr;
				dbgmres_smoother[i_52].row_ptr = nullptr;
			}
		}
		delete[] dbgmres_smoother;
		dbgmres_smoother = nullptr;
	}

	if (debug_reshime) system("pause");





	// Внимание: именно эта строчка обеспечивает сходимость.
#pragma omp parallel for
	for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
		x[i47] = x_best_search[i47];
	}

	identiti = true;
	for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
		if (fabs(x[i47] - x_best_search_init[i47]) > 1.0e-5) {
			identiti = false;
		}
	}

	GLOBAL_identity_situation_in_RUMBA_for_NEtworkT_solver = identiti;

	if (identiti) {
		if ((iVar != TOTALDEFORMATIONVAR) && (iVar != TURBULENT_KINETIK_ENERGY)) {
			std::cout << "identity situation" << std::endl;
			// если техника x_best_search вообще не дала результатов.
#pragma omp parallel for
			for (integer i47 = 1; i47 <= n_a[0]; ++i47) {
				x[i47] = x_best_search2[i47];
			}
		}
	}

	// Метка к которой мы приходим если значение невязки превысило 1.0e7.
FULL_DIVERGENCE_DETECTED:

	if (bprint_mesage_diagnostic) {
		std::cout << "number of negative diagonals: ibsp_length=" << ibsp_length << std::endl;
	}

	// диагностическое сообщение какую переменную мы решаем.
	/*if (bprint_mesage_diagnostic) {
		switch (iVar) {
		case PAM: std::cout << "PAM" << std::endl;  break;
		case VELOCITY_X_COMPONENT:  std::cout << "VELOCITY_X_COMPONENT"<< std::endl; break;
		case VELOCITY_Y_COMPONENT:  std::cout << "VELOCITY_Y_COMPONENT"<< std::endl; break;
		case VELOCITY_Z_COMPONENT:  std::cout << "VELOCITY_Z_COMPONENT"<< std::endl; break;
		case NUSHA: std::cout << "NU"<< std::endl;  break;
		case TURBULENT_KINETIK_ENERGY: std::cout << "TURBULENT_KINETIK_ENERGY"<< std::endl;  break;
		case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: std::cout << "TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA"<< std::endl;  break;
		case TURBULENT_KINETIK_ENERGY_STD_K_EPS:  std::cout << " TURBULENT_KINETIK_ENERGY_STD_K_EPS"<< std::endl; break;
		case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:  std::cout << "TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS"<< std::endl; break;
		case GAMMA_LANGTRY_MENTER: std::cout << "TURBULENT_GAMMA_LANGTRY_MENTER" << std::endl; break;
		case RE_THETA_LANGTRY_MENTER: std::cout << "TURBULENT_RE_THETA_LANGTRY_MENTER" << std::endl; break;
		case TEMP:  std::cout << "TEMP"<< std::endl; break;
		case TOTALDEFORMATIONVAR: std::cout << "Stress system"<< std::endl; break;
		}
	}
	else {*/

	if (bprint_mesage_diagnostic) {
		//switch (iVar) {
		// Радиатор водяного охлаждения 3л/мин ilevel_VX_VY_VZ=10, ilevel_PAM=5 или 6.

		//case PAM: std::cout << "PAM "<< ilevel<<" "<< n_a[ilevel - 4] / n_a[ilevel - 3]<<" "<< n_a[ilevel - 3] / n_a[ilevel-2]<< " "<< n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl;  break;
		//case VELOCITY_X_COMPONENT:  std::cout << "VELOCITY_X_COMPONENT "<< ilevel<<" "<< n_a[ilevel - 4] / n_a[ilevel - 3] <<" "<< n_a[ilevel - 3] / n_a[ilevel - 2] << " " << n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//case VELOCITY_Y_COMPONENT:  std::cout << "VELOCITY_Y_COMPONENT "<< ilevel<<" " << n_a[ilevel - 4] / n_a[ilevel - 3]<<" "<< n_a[ilevel - 3] / n_a[ilevel - 2]<<" "<< n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//case VELOCITY_Z_COMPONENT:  std::cout << "VELOCITY_Z_COMPONENT "<< ilevel<< " "<< n_a[ilevel - 4] / n_a[ilevel - 3]<<" "<< n_a[ilevel - 3] / n_a[ilevel - 2]<<" "<< n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//case NUSHA:  std::cout << "NUSHA %lld %e %e %e %e\n"<< ilevel<< " "<<  n_a[ilevel - 4] / n_a[ilevel - 3]<<" "<< n_a[ilevel - 3] / n_a[ilevel - 2]<<" "<< n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//case TURBULENT_KINETIK_ENERGY:  std::cout << "TURBULENT_KINETIK_ENERGY "<< ilevel<<" "<< n_a[ilevel - 4] / n_a[ilevel - 3]<<" " << n_a[ilevel - 3] / n_a[ilevel - 2]<<" "<< n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:  std::cout << "TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA "<< ilevel<< " " << n_a[ilevel - 4] / n_a[ilevel - 3]<<" "<< n_a[ilevel - 3] / n_a[ilevel - 2]<<" "<< n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//case TURBULENT_KINETIK_ENERGY_STD_K_EPS:  std::cout << " TURBULENT_KINETIK_ENERGY_STD_K_EPS " << ilevel<< " " <<  n_a[ilevel - 4] / n_a[ilevel - 3]<<" " << n_a[ilevel - 3] / n_a[ilevel - 2]<< " " << n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:  std::cout << "TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS  " << ilevel<<" "<< n_a[ilevel - 4] / n_a[ilevel - 3]<<" "<< n_a[ilevel - 3] / n_a[ilevel - 2]<< " " << n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//case GAMMA_LANGTRY_MENTER: std::cout << "TURBULENT_GAMMA_LANGTRY_MENTER"  << ilevel<< " " <<  n_a[ilevel - 4] / n_a[ilevel - 3]<<" " << n_a[ilevel - 3] / n_a[ilevel - 2]<< " " << n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel] << std::endl; break;
		//case RE_THETA_LANGTRY_MENTER: std::cout << "TURBULENT_RE_THETA_LANGTRY_MENTER"  << ilevel<< " " <<  n_a[ilevel - 4] / n_a[ilevel - 3]<<" " << n_a[ilevel - 3] / n_a[ilevel - 2]<< " " << n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel] << std::endl; break;
		//case TEMP:  std::cout << "TEMP "<< ilevel<< " "<<  n_a[ilevel - 4] / n_a[ilevel - 3]<< " "<<  n_a[ilevel - 3] / n_a[ilevel - 2]<< " " << n_a[ilevel - 2] / n_a[ilevel - 1]<<" "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//case TOTALDEFORMATIONVAR:  std::cout << "Stress system "<< ilevel<< " "<<  n_a[ilevel - 4] / n_a[ilevel - 3]<< " "<< n_a[ilevel - 3] / n_a[ilevel - 2]<<" "<< n_a[ilevel - 2] / n_a[ilevel - 1]<< " "<< n_a[ilevel - 1] / n_a[ilevel]<<std::endl; break;
		//}

		switch (iVar) {
			// Радиатор водяного охлаждения 3л/мин ilevel_VX_VY_VZ=10, ilevel_PAM=5 или 6.

		case PAM: std::cout << "PAM level=" << ilevel << "  CopA=" << dr_grid_complexity << " CopP=" << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV=" << icount_V_cycle << " res0=" << dres_initial << " n_a[ilevel - 2]=" << n_a[ilevel - 2] << " n_a[ilevel - 1]=" << n_a[ilevel - 1] << " n_a[ilevel]=" << n_a[ilevel] << std::endl;  break;
		case VELOCITY_X_COMPONENT:  std::cout << "VELOCITY_X_COMPONENT level=" << ilevel << " CopA=" << dr_grid_complexity << " CopP=" << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV=" << icount_V_cycle << " res0=" << dres_initial << "  n_a[ilevel - 2]=" << n_a[ilevel - 2] << "  n_a[ilevel - 1]=" << n_a[ilevel - 1] << " n_a[ilevel]=" << n_a[ilevel] << std::endl; break;
		case VELOCITY_Y_COMPONENT:  std::cout << "VELOCITY_Y_COMPONENT level=" << ilevel << " CopA=" << dr_grid_complexity << " CopP=" << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV=" << icount_V_cycle << " res0=" << dres_initial << "  n_a[ilevel - 2]=" << n_a[ilevel - 2] << "  n_a[ilevel - 1]=" << n_a[ilevel - 1] << " n_a[ilevel]=" << n_a[ilevel] << std::endl; break;
		case VELOCITY_Z_COMPONENT:  std::cout << "VELOCITY_Z_COMPONENT level=" << ilevel << " CopA=" << dr_grid_complexity << " CopP=" << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV=" << icount_V_cycle << " res0=" << dres_initial << "  n_a[ilevel - 2]=" << n_a[ilevel - 2] << "  n_a[ilevel - 1]=" << n_a[ilevel - 1] << " n_a[ilevel]=" << n_a[ilevel] << std::endl; break;
		case NUSHA:  std::cout << "NUSHA level=" << ilevel << " CopA=" << dr_grid_complexity << " CopP=" << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV=" << icount_V_cycle << " res0=" << dres_initial << "  n_a[ilevel - 2]=" << n_a[ilevel - 2] << "  n_a[ilevel - 1]=" << n_a[ilevel - 1] << " n_a[ilevel]=" << n_a[ilevel] << std::endl; break;
		case TURBULENT_KINETIK_ENERGY:  std::cout << "TURBULENT_KINETIK_ENERGY level=" << ilevel << " CopA=" << dr_grid_complexity << " CopP=" << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV=" << icount_V_cycle << " res0=" << dres_initial << "  n_a[ilevel - 2]=" << n_a[ilevel - 2] << "  n_a[ilevel - 1]=" << n_a[ilevel - 1] << " n_a[ilevel]=" << n_a[ilevel] << std::endl; break;
		case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:  std::cout << "TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA level=" << ilevel << " CopA=" << dr_grid_complexity << " CopP=" << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV=" << icount_V_cycle << " res0=" << dres_initial << "  n_a[ilevel - 2]=" << n_a[ilevel - 2] << "  n_a[ilevel - 1]=" << n_a[ilevel - 1] << " n_a[ilevel]=" << n_a[ilevel] << std::endl; break;
		case TURBULENT_KINETIK_ENERGY_STD_K_EPS:  std::cout << " TURBULENT_KINETIK_ENERGY_STD_K_EPS level = " << ilevel << " CopA = " << dr_grid_complexity << " CopP = " << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV = " << icount_V_cycle << " res0 = " << dres_initial << "  n_a[ilevel - 2] =" << n_a[ilevel - 2] << "  n_a[ilevel - 1] = " << n_a[ilevel - 1] << " n_a[ilevel] = " << n_a[ilevel] << std::endl; break;
		case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:  std::cout << "TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS level = " << ilevel << " CopA = " << dr_grid_complexity << " CopP = " << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV = " << icount_V_cycle << " res0 =" << dres_initial << "  n_a[ilevel - 2] = " << n_a[ilevel - 2] << "  n_a[ilevel - 1] = " << n_a[ilevel - 1] << " n_a[ilevel] = " << n_a[ilevel] << std::endl;  break;
		case GAMMA_LANGTRY_MENTER: std::cout << "TURBULENT_GAMMA_LANGTRY_MENTER" << ilevel << " CopA = " << dr_grid_complexity << " CopP = " << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV = " << icount_V_cycle << " res0 =" << dres_initial << "  n_a[ilevel - 2] = " << n_a[ilevel - 2] << "  n_a[ilevel - 1] = " << n_a[ilevel - 1] << " n_a[ilevel] = " << n_a[ilevel] << std::endl; break;
		case RE_THETA_LANGTRY_MENTER: std::cout << "TURBULENT_RE_THETA_LANGTRY_MENTER" << ilevel << " CopA = " << dr_grid_complexity << " CopP = " << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV = " << icount_V_cycle << " res0 =" << dres_initial << "  n_a[ilevel - 2] = " << n_a[ilevel - 2] << "  n_a[ilevel - 1] = " << n_a[ilevel - 1] << " n_a[ilevel] = " << n_a[ilevel] << std::endl; break;
		case TEMP:  std::cout << "TEMP level=" << ilevel << " CopA=" << dr_grid_complexity << " CopP=" << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV=" << icount_V_cycle << " res0=" << dres_initial << "  n_a[ilevel - 2]=" << n_a[ilevel - 2] << "  n_a[ilevel - 1]=" << n_a[ilevel - 1] << " n_a[ilevel]=" << n_a[ilevel] << std::endl; break;
		case TOTALDEFORMATIONVAR:  std::cout << "Stress system  level=" << ilevel << " CopA=" << dr_grid_complexity << " CopP=" << static_cast<doublereal>(nnz_P_memo_all / n_a[0]) << " nV=" << icount_V_cycle << " res0=" << dres_initial << "  n_a[ilevel - 2]=" << n_a[ilevel - 2] << "  n_a[ilevel - 1]=" << n_a[ilevel - 1] << " n_a[ilevel]=" << n_a[ilevel] << std::endl; break;
		}

	}
	//}


	// free
	if (x_best_search2 != nullptr) {
		delete[] x_best_search2;
		x_best_search2 = nullptr;
	}
	if (x_best_search_init != nullptr) {
		delete[] x_best_search_init;
		x_best_search_init = nullptr;
	}


	// free	
	if (bnested_desection_global_amg != nullptr) {
		free(bnested_desection_global_amg);  // Глобальная память.
		bnested_desection_global_amg = nullptr;
	}
	if (diag_minus_one[ilevel + 1] != nullptr) {
		free(diag_minus_one[ilevel + 1]);
		diag_minus_one[ilevel + 1] = nullptr;
	}
	for (integer i_scan_levels = 0; i_scan_levels <= maxlevel - 1; ++i_scan_levels) {
		if (ilevel + 1 > i_scan_levels) {
			// free			
			if (diag[i_scan_levels] != nullptr) {
				free(diag[i_scan_levels]);
				diag[i_scan_levels] = nullptr;
			}
			if (diag_minus_one[i_scan_levels] != nullptr) {
				free(diag_minus_one[i_scan_levels]);
				diag_minus_one[i_scan_levels] = nullptr;
			}
			if (nested_desection[i_scan_levels] != nullptr) {
				free(nested_desection[i_scan_levels]);
				nested_desection[i_scan_levels] = nullptr;
			}
		}
		if (ilevel + 2 > i_scan_levels) {
				integer i_scan_levels_prev = i_scan_levels - 1;
				if (i_scan_levels_prev >= 0) {
					if (error_approx_fine[i_scan_levels_prev] != nullptr) {
						delete[] error_approx_fine[i_scan_levels_prev];
						//free(error_approx_fine[i_scan_levels_prev]);
						error_approx_fine[i_scan_levels_prev] = nullptr;
					}
					if (error_approx_coarse[i_scan_levels_prev] != nullptr) {
						delete[] error_approx_coarse[i_scan_levels_prev];
						//free(error_approx_coarse[i_scan_levels_prev]);
						error_approx_coarse[i_scan_levels_prev] = nullptr;
					}
					if (residual_coarse[i_scan_levels_prev] != nullptr) {
						delete[] residual_coarse[i_scan_levels_prev];
						//free(residual_coarse[i_scan_levels_prev]);
						residual_coarse[i_scan_levels_prev] = nullptr;
					}
					if (residual_fine[i_scan_levels_prev] != nullptr) {
						delete[] residual_fine[i_scan_levels_prev];
						//free(residual_fine[i_scan_levels]);
						residual_fine[i_scan_levels_prev] = nullptr;
					}
				}
				
			
		}
	}


	// метод огрубления.
	my_amg_manager.icoarseningtype = memo_icoarseningtype;


	

	if (diag != nullptr) {
		delete[] diag;
		diag = nullptr;
	}
	if (diag_minus_one != nullptr) {
		delete[] diag_minus_one;
		diag_minus_one = nullptr;
	}
	if (nested_desection != nullptr) {
		delete[] nested_desection;
		nested_desection = nullptr;
	}

	

	

	

	// освобождение оперативной памяти.
	free_level_additional_data(milu0, ilevel);
	free_level_additional_data(milu2, ilevel);

	// Освобождение общей памяти в ILU буфере.
	if (milu_gl_buffer.alu_copy != nullptr) delete[] milu_gl_buffer.alu_copy;
	if (milu_gl_buffer.jlu_copy != nullptr) delete[] milu_gl_buffer.jlu_copy;
	if (milu_gl_buffer.ju_copy != nullptr) delete[] milu_gl_buffer.ju_copy;
	milu_gl_buffer.alu_copy = nullptr;
	milu_gl_buffer.jlu_copy = nullptr;
	milu_gl_buffer.ju_copy = nullptr;


	//if (residual_fine[0] != nullptr) {
		//delete[] residual_fine[0];
		//residual_fine[0] = nullptr;
		//free(residual_fine[0]);
		//residual_fine[0] = nullptr;
	//}


	if (residual_fine != nullptr) {
		delete[] residual_fine;
		residual_fine = nullptr;
	}


	if (error_approx_fine != nullptr) {
		delete[] error_approx_fine;
		error_approx_fine = nullptr;
	}


	if (error_approx_coarse != nullptr) {
		delete[] error_approx_coarse;
		error_approx_coarse = nullptr;
	}


	if (residual_coarse != nullptr) {
		delete[] residual_coarse;
		residual_coarse = nullptr;
	}


	if (row_ptr_start != nullptr) {
		free(row_ptr_start);
		row_ptr_start = nullptr;
	}
	if (row_ptr_end != nullptr) {
		free(row_ptr_end);
		row_ptr_end = nullptr;
	}

	
	if (x_copy != nullptr) {
		free(x_copy);
		x_copy = nullptr;
	}
	if (x_old != nullptr) {
		free(x_old);
		x_old = nullptr;
	}
	if (x_best_search != nullptr) {
		free(x_best_search);
		x_best_search = nullptr;
	}

	// Для подстраховки:

	if (row_ptr_start != nullptr) {
		free(row_ptr_start);
		row_ptr_start = nullptr;
	}
	if (row_ptr_end != nullptr) {
		free(row_ptr_end);
		row_ptr_end = nullptr;
	}


	if (x_jacoby_buffer != nullptr) {
		delete[] x_jacoby_buffer;
		x_jacoby_buffer = nullptr;
	}


	free_hash_table_Gus_struct01();

	return ret_value;
}


// Разбиение множества узлов первоначальной сетки на С-узлы и 
// F - узлы. С -узлы составят грубую сетку следующего уровня вложенности.
// Значения же функции в F узлах должно быть восстановлено по значению функции
// в ближайших С узлах (см. задачу интерполяции). 
template <typename doublerealT>
void Ruge_and_Stuben_CF_decomposition(Ak2& Amat, bool*& this_is_F_node,
	bool*& this_is_C_node, integer ilevel,
	integer*& count_neighbour, integer*& n_a,
	integer& newCcount,
	doublerealT*& threshold_quick_only_negative,
	integer const *const row_startA,
	integer* &hash_StrongTranspose_collection1Eco,
	integer* & isize_hash_StrongTranspose_collection,
	bool bprint_mesage_diagnostic,
	bool*& flag, integer iadd, 
	bool bpositive_connections_CF_decomp,
	bool bStrongTransposeON,
	bool* &hash_table2,
	integer* &istack, integer const *const nnz_a,
	bool debug_reshime)
{

	

	integer ii_end1 = n_a[ilevel - 1];
	integer max_neighbour = 0;
	integer icandidate = 0;

	// Находим узел с наибольшим числом соседей и запоминаем его.
	// Это первый встретившийся узел с наибольшим числом соседей.
	// Это требуется для того чтобы стартовал алгоритм C/F разбиения.

	for (integer i7 = 1; i7 <= ii_end1; ++i7) {
		if (count_neighbour[i7] > max_neighbour) {
			max_neighbour = count_neighbour[i7];
			icandidate = row_startA[i7];
		}
	}




	// нужно выделить кто попал в coarse, а кто в этот раз попал в Fine Выделить всех кто соседствует
	// с новыми Fine увеличить им счётчик соседей.


	//const integer NULL_NEIGHBOUR = -1;
	// Построение C/F разбиения.
	//integer vacant = NULL_NEIGHBOUR;
	bool bcontinue_gl_1 = true;
	// Построение C/F разбиения.
	integer icountprohod = 0;

	// храним те узлы которые уже были пройдены при конструировании.
	// поначалу все узлы помечены как непосещённые.
	bool* bmarkervisit = nullptr;
	//if (bmarkervisit != nullptr) {
		//free(bmarkervisit);
		//bmarkervisit = nullptr;
	//}
	bmarkervisit = my_declaration_array<bool>(n_a[ilevel - 1], false, "bmarkervisit");


	// увеличение быстродействия достигается 
	// сокращением пределов сканирования
	// здесь хранится индекс начала сканирования flag.
	integer istartflag_scan = 1;
	
	// Задача 12mm hfet thermal resistance. 1.7млн неизвестных.
	// AVL_TREE_ID   3мин 29с 590мс      {5}
	// SPLAY_TREE_ID  3мин 16с 430мс {2}
	// BINARY_HEAP 3мин 4с 0мс {1 *самая быстрая.}
	// RANDOM_TREE_ID (Дерамида) 3мин 28с 90мс {4}
	// RED_BLACK_TREE_ID 3мин 27с 210мс {3}


	//const integer AVL_TREE_ID = 0;   // АВЛ дерево поиска. 12.12.2015.
	//const integer SPLAY_TREE_ID = 1; // Скошенное дерево поиска.
	//const integer BINARY_HEAP = 2; // Двоичная куча. 16.06.2017.
	//const integer RANDOM_TREE_ID = 3; // (Дерамида) Рандомизированное дерево поиска. 24.08.2017.
	//const integer RED_BLACK_TREE_ID = 4; // Красно-Чёрное дерево поиска. 22.06.2018.
	//const integer FIBONACCI_HEAP_ID = 5; // Фибоначчиева куча. 11.07.2018.
	//const integer VAN_EMDE_BOAS_TREE_ID = 6; // ван Эмде Боас дерево поиска. 30.06.2018
	//integer id_tree = BINARY_HEAP; // AVL_TREE_ID; // SPLAY_TREE_ID; // BINARY_HEAP; // RANDOM_TREE_ID; // RED_BLACK_TREE_ID;
	// 28.01.2018 На выбор пользователя.

	RS_COARSENING_KERNEL_DATA_STRUCTURE id_tree = my_amg_manager.iCFalgorithm_and_data_structure;

	integer n = n_a[0];

	// Выделяем память под двоичную кучу.
	// Деструктор вызывается автоматом при уходе из области видимости области определения.
	const integer isize_priority_queue01 = static_cast<integer>(0.4 * n); // 0.238
	integer ikonst1 = isize_priority_queue01, ikonst2 = n;
	if (id_tree != RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP) {
		ikonst1 = 0;
		ikonst2 = 0;
	}
	PQ<integer> binary_heap(ikonst1, ikonst2); // 500K для 2.1M

	// Фибоначчиева куча.
	FibonacciHeap<integer> fibo_heap;

	if (id_tree == RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP) {
		fibo_heap.WakeUp2(n + 1);// alloc memory hash table
	}

	node_AVL* root = 0;
	Tree_splay* root_splay = 0;
	size_splay_Tree = 0;
	TreapNode* random_tree_root = nullptr;
	RBtree RBroot; // Корень Красно-Чёрного дерева.

	root = 0;
	root_splay = 0;
	size_splay_Tree = 0;
	random_tree_root = nullptr;

	if (id_tree == RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP) {
		binary_heap.clear();
	}


	RBroot.Clear();
	if (id_tree == RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP) {
		fibo_heap.UpdateSize(n_a[ilevel - 1] + 1);
	}

#if VEB_FLAG
	int64_t res_vanEMDE_BOAS_Tree;
	int64_t universe = 4294967296; // 2 ^32=2^(2^5) (4294 млн) работает
	//int64_t universe = 67108864; // 2^26 не работает
	//int64_t universe = 134217728; // 2^27 не работает
	TvEB* vanEMDE_BOAS_Tree = nullptr;

	if (id_tree == RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE) {
		vanEMDE_BOAS_Tree = new TvEB(universe);
	}
#endif

	newCcount = 0;

	


	// Нехорошо постоянно выделять и уничтожать память в длинном цикле, 
	// более быстро выделить её один раз. См. выделение памяти под set.
	// 23.04.2017

	// Вынес на уровень выше 26,06,2021
	doublereal* max_vnediagonal33_arr = new doublereal[n_a[ilevel - 1] + 1];

#pragma omp parallel for
	for (integer i33 = 0; i33 <= n_a[ilevel - 1]; ++i33) {
		max_vnediagonal33_arr[i33] = -1.0e30;

		integer istart72 = row_startA[i33];
		integer istopmarker2 = row_startA[i33 + 1] - 1;

		// 22 _12_2016
					// Это лучший вариант: обеспечивает корректное построение иерархии
					// уровней на задаче passive module 6 в то время как все остальные 
					// отличные от этого способа давали сбой.

		for (integer is01 = istart72; (is01 <= istopmarker2); ++is01) {
			if (Amat.j[is01] != Amat.i[is01]) {
				if ((Amat.aij[is01] < 0.0) && (Amat.abs_aij[is01] > max_vnediagonal33_arr[i33])) {
					max_vnediagonal33_arr[i33] = Amat.abs_aij[is01];
				}
			}
		}
	}
	

	if (bprint_mesage_diagnostic) {
		std::cout << "   ***   CAMG SELECTOR "<< ilevel <<"  ***\n";
	}
	while (bcontinue_gl_1)
	{

		
		
		integer ii = icandidate;
		if (icandidate < 0) {
			// 15.09.2020
			icandidate = 0;
			bcontinue_gl_1 = false;
			break;
		}
		if (flag[Amat.i[ii]] == false) {

			integer ic = 0; // Обязательная инициализация.


			integer ic_end_F_SiTranspose = 0;
			integer set0 = Amat.i[ii];




			//A20.05.2017//this_is_C_node[set[0]] = true;
			//A20.05.2017//bmarkervisit[set[0]] = true;
			this_is_C_node[set0] = true;
			bmarkervisit[set0] = true;

			doublerealT max_vnediagonal = -1.0; // максимальное значение модуля вне диагонального элемента. 
												// добавляем диагональный элемент.
												// узел set[0]==Amat.i[is0].
												// Нахождение значения максимального внедиагонального элемента, с 
												// учётом того что даже узел Дирихле связан с одним внутренним узлом расчётной области.
												// 17 января 2016 правильное определение максимального внедиагонального элемента.
												// Обязательная перемотка в самое начало строки.
			integer ii_back = ii;
			while ((ii_back > iadd) && (Amat.i[ii_back] == set0)) ii_back--;
			ii_back++;

			doublerealT max_vnediagonal1 = static_cast <doublerealT>(-1.0e30);
			
			// Если делать по максимальному внедиагональному элементу то мы получим очень много элементов на грубых уровнях,
			// и чрезвычайно медленную сходимость.

			if (bpositive_connections_CF_decomp) {
				// 23_10_2016
				for (integer is0 = ii_back; (is0 <= row_startA[set0 + 1] - 1); ++is0) {
					if (Amat.j[is0] != set0) {
						if (Amat.abs_aij[is0] > max_vnediagonal1) {
							max_vnediagonal1 = Amat.abs_aij[is0]; //i,j
							// Большое количество элементов на грубых уровнях,
							// очень медленная сходимость.
							//if (Amat.j[is0] == set[0]) break; 
						}						
					}
				}
			}
			else {
				for (integer is0 = ii_back; (is0 <= row_startA[set0 + 1] - 1); ++is0) {
					if (Amat.j[is0] != set0) {
						if (Amat.aij[is0] < 0.0) {
							if (Amat.abs_aij[is0] > max_vnediagonal1) {
								max_vnediagonal1 = Amat.abs_aij[is0]; //i,j
								// Большое количество элементов на грубых уровнях,
								// очень медленная сходимость.
								//if (Amat.j[is0] == set[0]) break; 
							}							
						}
					}
				}
			}
			
			
			//max_vnediagonal = max_vnediagonal1;  // 1			
			// наиболее близка к оптимальной. -85%. но несомненно лучше max_vnediagonal = -1.0;
			// 19 января 2016 установлено что важны все связи, не нужно учитывать threshold
			// max_vnediagonal должно быть -1.0. Именно это значение обеспечивает наилучшую 
			// скорость агломерации и наилучшую скорость сходимости.
			max_vnediagonal = static_cast <doublerealT>(-1.0e30);  // все связи!!!											

			ic++;


			//  В set начиная с единицы и до <ic лежат кандидаты чтобы стать F.
			// 5.01.2017
			// 01.04.2017 Дополняемся F узлами из Si_Transpose связей.
			if ((my_amg_manager.ipatch_number == 7) && (bStrongTransposeON)) {

				integer imarker75_scan = 0;

				// обычный линейный список.
				/*formirate_F_SiTranspose_hash_table_Gus2_struct03(hash_StrongTranspose_collection1[Amat.i[ii]],
					isize_hash_StrongTranspose_collection[Amat.i[ii]], imarker75_scan, this_is_F_node, this_is_C_node);*/

				formirate_F_SiTranspose_hash_table_Gus2_struct04(hash_StrongTranspose_collection1Eco, Amat.i[ii], n_a[ilevel-1],
					isize_hash_StrongTranspose_collection[Amat.i[ii]], imarker75_scan, this_is_F_node, this_is_C_node);

				ic = imarker75_scan + 1;
			}

			ic_end_F_SiTranspose = ic; // С этой позиции заканчиваются F которые из Si_Transpose.

									   // если узел j ещё не был добавлен в агрегат.
			if (bpositive_connections_CF_decomp) {
				if (flag[Amat.j[ii]] == false) {
					if ((Amat.j[ii] != set0) && (Amat.abs_aij[ii] >= theta(ilevel) * max_vnediagonal)) {
						// 21.05.2017
						bool bfound_vacant = false;

						bfound_vacant = isfound_hash_table_Gus_struct01(Amat.j[ii]);
						if (!bfound_vacant) {
							insert_hash_table_Gus_struct01(Amat.j[ii]);
							ic++;
						}

					}
				}
			}
			else {
				if (flag[Amat.j[ii]] == false) {
					if ((Amat.j[ii] != set0) && (Amat.aij[ii] < 0.0) && (Amat.abs_aij[ii] >= theta(ilevel) * max_vnediagonal)) {
						// 21.05.2017
						bool bfound_vacant = false;

						bfound_vacant = isfound_hash_table_Gus_struct01(Amat.j[ii]);
						if (!bfound_vacant) {
							insert_hash_table_Gus_struct01(Amat.j[ii]);
							ic++;
						}

					}
				}
			}

			//std::cout<<"sboi start";

			// iscan = ii+1; // устаревший код
			integer iscan = ii_back + 1; // важная модификация 19 января 2016г.
			// 19 jan 2016.

			if (bpositive_connections_CF_decomp) {
				while ((iscan <= nnz_a[ilevel - 1] + iadd) && (Amat.i[iscan] == set0)) {
					// 14 февраля 2016 код иногда приводящий к сбою.
					//while (iscan <= row_startA[set0 + 1] - 1) { // код иногда приводящий к сбою по непонятной причине.
					// если узел j ещё не был добавлен в агрегат.
					if (flag[Amat.j[iscan]] == false) {
						if ((Amat.j[iscan] != set0) && (Amat.abs_aij[iscan] >= theta(ilevel) * max_vnediagonal)) {
							// 21.05.2017
							bool bfound_vacant = false;

							bfound_vacant = isfound_hash_table_Gus_struct01(Amat.j[iscan]);
							if (!bfound_vacant) {
								insert_hash_table_Gus_struct01(Amat.j[iscan]);
								ic++;
							}

							/*
							// Медленная версия с линейным поиском.
							vacant = Amat.j[iscan];
							for (integer js = 0; js < ic; ++js) {
							if (vacant == set[js]) {
							vacant = NULL_NEIGHBOUR;
							}
							}
							if (vacant != NULL_NEIGHBOUR) {
							set[ic] = vacant;

							ic++;

							}
							*/
						}
					}

					iscan++;

				} // while
			}
			else {
				while ((iscan <= nnz_a[ilevel - 1] + iadd) && (Amat.i[iscan] == set0)) {
					// 14 февраля 2016 код иногда приводящий к сбою.
					//while (iscan <= row_startA[set0 + 1] - 1) { // код иногда приводящий к сбою по непонятной причине.
					// если узел j ещё не был добавлен в агрегат.
					if (flag[Amat.j[iscan]] == false) {
						if ((Amat.j[iscan] != set0) && (Amat.aij[iscan] < 0.0) && (Amat.abs_aij[iscan] >= theta(ilevel) * max_vnediagonal)) {
							// 21.05.2017
							bool bfound_vacant = false;

							bfound_vacant = isfound_hash_table_Gus_struct01(Amat.j[iscan]);
							if (!bfound_vacant) {
								insert_hash_table_Gus_struct01(Amat.j[iscan]);
								ic++;
							}

							/*
							vacant = Amat.j[iscan];
							for (integer js = 0; js < ic; ++js) {
							if (vacant == set[js]) {
							vacant = NULL_NEIGHBOUR;
							}
							}
							if (vacant != NULL_NEIGHBOUR) {
							set[ic] = vacant;

							ic++;

							}
							*/
						}
					}

					iscan++;

				} // while
			}

			//std::cout<<"sboi end";
			// Это была учтена только связь i,j






		// В этом месте множество set успешно сформировано:
		// 1. Перепаковка из root_Gus_set в set.
		// 2. root_Gus_set больше не используется.
		// 3. Именно здесь надо выделить данные под set.
			integer* set = nullptr;
			set = new integer[ic + 2];
			//if (set == nullptr) {
				//std::cout<<"error!!! memory for set is nullptr. Problem allocate detected."<<std::endl;
				//std::cout<<"in function classic_aglomerative_amg6."<<std::endl;
				//system("pause");
				//exit(1);
			//}

			integer ic_986 = 1;
			set[0] = set0;



			formirate_hash_table_Gus_struct01__2__set(set, ic_986);

			clear_hash_table_Gus_struct01();


			for (integer isc = 1; isc < ic; ++isc) {
				this_is_F_node[set[isc]] = true; // это только новые F узлы.
				bmarkervisit[set[isc]] = true;
			}




			// Помечаем узлы как включённые в агрегат.
			for (integer js = 0; js < ic; ++js) {
				flag[set[js]] = true;
			}






			// Алгоритм (5 декабря 2015 revised) 
			// 1. Сканируем все F которые соседи данного С на данном проходе.
			// 2. Для каждого фиксированного F сканируем его "строчных" соседей.
			// 3. Если узел еще не был включён в агрегат то ищем всех соседей данного узла на предмет 
			// соседства с фиксированным набором F из пункта 1.




			TreapNode* nrt_temp = nullptr;
			TreapNode* save_root = nullptr;

			/*
			if (id_tree == FIBONACCI_HEAP_ID) {
			if (!fibo_heap.isEmpty()) {
			fibo_heap.removeMinimum();
			}
			}
			*/
			// 12 декабря 2015.
			// Надо удалить из АВЛ дерева C и F узлы.
			// Это удаление очищает АВЛ дерево и приводит его к
			// рабочему состоянию. Удаление несуществующих в дереве узлов
			// производится корректно. Удаление производится за логарифмическое
			// по основанию 2  время от количества элементов в дереве
			// сбалансированность дерева при этом сохраняется.
			for (integer js = 0; js < ic; ++js) {
				data_BalTree ddel;
				ddel.i = set[js];
				ddel.count_neighbour = count_neighbour[set[js]];
				// Уникальный ключ для дерева ван Эмде Боаса.
#if VEB_FLAG
				integer  veb_del_key = (count_neighbour[set[js]]) * (n_a[ilevel - 1] + 1) + (set[js]);
				if (id_tree == RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE) {
					if (veb_del_key > universe - 2) {
						std::cout<<"overflow veb-Van Emde Boas 2^2^5"<<std::endl;
						system("PAUSE");
					}
					if (veb_del_key < 1) {
						std::cout<<"overflow veb-Van Emde Boas < 1"<<std::endl;
						system("PAUSE");
					}
				}
#endif
				//ddel.ii = row_startA[ddel.i];
				switch (id_tree) {
				case RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE: root = remove_AVL(root, ddel);
					break;
				case RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE: root_splay = delete_splay_Tree(ddel, root_splay);
					break;
				case RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP:
					// Уникальным ключом удаления является set[js].
					binary_heap.remove(set[js]);
					break;
				case RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP:
					save_root = random_tree_root;
					nrt_temp = search(random_tree_root, ddel);
					random_tree_root = save_root;
					save_root = nullptr;
					if (nrt_temp != nullptr) {
						nrt_temp = nullptr;
						random_tree_root = deleteNode(random_tree_root, ddel);
					}
					break;
				case RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE:
					RBroot.Remove(ddel);
					break;
				case RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP:
					if (!fibo_heap.isEmpty()) {
						fibo_heap.deleteKey(ddel);
					}
					break;
				case RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE:
#if VEB_FLAG
					// Если элемент присутствует то мы его удалим.
					res_vanEMDE_BOAS_Tree = vEB_find(vanEMDE_BOAS_Tree, veb_del_key);
					if (!res_vanEMDE_BOAS_Tree) {

					}
					else {
						res_vanEMDE_BOAS_Tree = vEB_delete(vanEMDE_BOAS_Tree, veb_del_key);
						if (!res_vanEMDE_BOAS_Tree) {
							std::cout<<"cannot be deleted post factum delete ddel.count_neighbour=="<<ddel.count_neighbour<<", ddel.i=="<<ddel.i<<std::endl;
							system("PAUSE");
						}
					}
#endif
					break;
				default: root = remove_AVL(root, ddel);
					break;
				}
			}




			//std::cout<<"additional and modify new neighbour"<<std::endl;

			// 10 января 2016. Новая логика.
			// Устраним некоторые повторные модификации (это должно снизить нагрузку на АВЛ дерево).
			// Эта модификация даёт сокращение количества V циклов которые требуются до сходимости
			// Эта модификация наиболее близка к классической описанной в литературе чем все предыдущие.
			// На момент 13 января 2016 это лучший вариант по скорости вычислений.
			integer itop_stack2 = 0;


		


			// 10 января 2016. Старый вариант просто очищенный от устаревшего кода.
			for (integer js = 1; js < ic; ++js) {

				integer i_11 = set[js];
				integer ii_11 = row_startA[i_11];
				//integer iend5 = nnz_a[ilevel - 1] + iadd;
				//integer istart73 = ii_11;
				//while ((istart73 >= 1 + iadd) && (Amat.i[istart73] == Amat.i[ii_11])) istart73--;
				//istart73++;
				integer istart73 = row_startA[Amat.i[ii_11]];
				integer iend73 = row_startA[Amat.i[ii_11] + 1] - 1;
				//bool bvisitsos = false;
				for (integer is0 = istart73; (is0 <= iend73); ++is0) {
					//for (integer is0 = istart73; (is0 <= iend5) && (Amat.i[is0] == Amat.i[ii_11]); ++is0) {
					// В пересечении с U!!!
					if (flag[Amat.j[is0]] == false) {


						integer isc = Amat.j[is0];



						// Избавляемся от повторных инкрементаций.
						// В 2D на пятиточечном шаблоне повторные инкрементации составляют
						// около 33%.
						// Это даёт стандартный алгоритм сгрубления, описаный в статьях, но
						// на ряде тестовых задач при таком подходе агломерация проходила очень
						// плохо (переполнение по памяти, не хватало даже семикратного размера исходной матрицы).
						// Эта проблема проявилась на задачах:
						// CGHV1J006D, Потенциал тора, Электрический потенциал в FET, Module 2.
						// Плохая скорость агломерации получается главным образом из-за шестого способа интерполяции.
						// Проблема не в этом месте кода.
						if (hash_table2[isc] == false) {
							hash_table2[isc] = true;
							istack[itop_stack2] = isc;
							itop_stack2++;
							// закомментированный лучше.
							//}

							//21_12_2016
							integer ii_2 = row_startA[isc];


							integer ic2 = 0;
							//integer iend2loc = nnz_a[ilevel - 1] + iadd;
							//integer istart72 = ii_2;										
							integer istart72 = row_startA[Amat.i[ii_2]];
							integer istopmarker2 = row_startA[Amat.i[ii_2] + 1] - 1;

							// 22 _12_2016
							// Это лучший вариант: обеспечивает корректное построение иерархии
							// уровней на задаче passive module 6 в то время как все остальные 
							// отличные от этого способа давали сбой.

							// Вынес на уровень выше 26,06,2021
							/*doublerealT max_vnediagonal33 = -1.0e30;
							for (integer is01 = istart72; (is01 <= istopmarker2); ++is01) {
								if (Amat.j[is01] != Amat.i[is01]) {
									if ((Amat.aij[is01] < 0.0) && (Amat.abs_aij[is01] > max_vnediagonal33)) {
										max_vnediagonal33 = Amat.abs_aij[is01];
									}
								}
							}*/

							const doublereal max_vnediagonal33_limit = 0.2375 * max_vnediagonal33_arr[Amat.i[ii_2]];

#pragma omp parallel for reduction(+: ic2)
							for (integer is01 = istart72; (is01 <= istopmarker2); ++is01) {
								// 0.2375 импирически подобрана для passive module 6.
								if ((Amat.aij[is01] < 0.0) && (Amat.abs_aij[is01] > max_vnediagonal33_limit)) {
									if (Amat.j[is01] == set[js]) {
										if ((my_amg_manager.ipatch_number == 7) && (bStrongTransposeON)) {
											if (js < ic_end_F_SiTranspose) {
												// Увеличиваем счётчики только тех соседей F узлов которые
												// являются соседями F узлов которые были получены из Si_Transpose связей.
												// Именно так написано у Джона Руге и Клауса Штубена.
												ic2 = ic2 + 1;// ic2++;
											}
										}
										else {
											ic2 = ic2 + 1; //ic2++;
										}
									}
								}
								// уменьшить счетчик слабого (weakly) соседа ?
							}


							data_BalTree dsearch;
							dsearch.count_neighbour = count_neighbour[isc];
							//dsearch.ii = ii_2;
							dsearch.i = isc;
							// Увеличиваем на количество связей с новыми F узлами.
							count_neighbour[isc] += ic2;
							data_BalTree dadd;
							dadd.count_neighbour = count_neighbour[isc];
							//dadd.ii = ii_2;
							dadd.i = isc;

							// Уникальный ключ для дерева ван Эмде Боаса.
							integer  veb_dadd_key = (dadd.count_neighbour) * (n_a[ilevel - 1] + 1) + (dadd.i);
							integer  veb_dsearch_key = (dsearch.count_neighbour) * (n_a[ilevel - 1] + 1) + (dsearch.i);
							//integer  veb_dadd_key = (dadd.count_neighbour)*(n + 1) + (dadd.i);
							//integer  veb_dsearch_key = (dsearch.count_neighbour)*(n + 1) + (dsearch.i);

							if (id_tree == RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE) {

#if VEB_FLAG							
								if (veb_dadd_key > universe - 2) {
									std::cout<<"overflow veb-Van Emde Boas 2^2^5"<<std::endl;
									system("PAUSE");
								}
								if (veb_dsearch_key > universe - 2) {
									std::cout<<"overflow veb-Van Emde Boas 2^2^5"<<std::endl;
									system("PAUSE");
								}
							

							
								if (veb_dadd_key < 1) {
									std::cout << "overflow veb-Van Emde Boas <1 " << std::endl;
									system("PAUSE");
								}
								if (veb_dsearch_key < 1) {
									std::cout << "overflow veb-Van Emde Boas <1 " << std::endl;
									system("PAUSE");
								}

#endif
							}


							TreapNode* nrt_temp_1 = nullptr;
							TreapNode* save_root_1 = nullptr;

							// добавляем элемент в АВЛ дерево,
							// причём если элемент уже находился в дереве то он модифицируется.
							// 12 декабря 2015.
							// Добавление узла происходит за логарифмическое по основанию 2 время,
							// причём после добавления дерево остаётся сбалансированным.
							// Г.М. Адельсон-Вельский и Е.М. Ландис 1962.
							switch (id_tree)
							{
							case RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE: root = insert_and_modify(root, dadd, dsearch);
								break;
							case RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE: root_splay = insert_and_modify(root_splay, dadd, dsearch);
								break;
							case RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP:
								if (binary_heap.isfound(isc)) {
									// Найден
									// Удаляем существующий элемент и вставляем новый.
									binary_heap.remove(isc);
									// Осуществляем вставку нового элемента.
									binary_heap.insert(count_neighbour[isc], isc);
								}
								else {
									// отсутствует.
									// Осуществляем вставку нового элемента.
									binary_heap.insert(count_neighbour[isc], isc);
								}
								break;
							case RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP:
								nrt_temp_1 = nullptr;
								save_root_1 = random_tree_root;
								nrt_temp_1 = search(random_tree_root, dsearch);
								random_tree_root = save_root_1;
								save_root_1 = nullptr;
								if (nrt_temp_1 == nullptr) {
									// Элемент в дереве отсутствует.
									random_tree_root = insert(random_tree_root, dadd);
								}
								else {
									nrt_temp_1 = nullptr;
									// Удаление
									random_tree_root = deleteNode(random_tree_root, dsearch);
									// Вставка
									random_tree_root = insert(random_tree_root, dadd);
								}
								break;
							case RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE:
								RBroot.InsertAndModify(dadd, dsearch);
								break;
							case RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP:
								fibo_heap.insert_and_modify(-veb_dsearch_key, -veb_dadd_key);
								break;
							case RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE:
#if VEB_FLAG
								res_vanEMDE_BOAS_Tree = vEB_find(vanEMDE_BOAS_Tree, veb_dsearch_key);
								if (!res_vanEMDE_BOAS_Tree) {
									// не найден
									res_vanEMDE_BOAS_Tree = vEB_find(vanEMDE_BOAS_Tree, veb_dadd_key);
									if (!res_vanEMDE_BOAS_Tree) {
										// не найден
										res_vanEMDE_BOAS_Tree = vEB_insert(vanEMDE_BOAS_Tree, veb_dadd_key);
										if (!res_vanEMDE_BOAS_Tree) {
											std::cout<<"insert problem veb veb_dadd_key=="<< veb_dadd_key<<std::endl;
										}
									}
								}
								else {
									res_vanEMDE_BOAS_Tree = vEB_delete(vanEMDE_BOAS_Tree, veb_dsearch_key);
									if (!res_vanEMDE_BOAS_Tree) {
										std::cout<<"cannot be deleted post factum delete veb_dsearch_key=="<< veb_dsearch_key<<std::endl;
										system("PAUSE");
									}
									// найден, удален м вставлен == заменен.
									res_vanEMDE_BOAS_Tree = vEB_insert(vanEMDE_BOAS_Tree, veb_dadd_key);
									if (!res_vanEMDE_BOAS_Tree) {
										std::cout<<"insert problem veb veb_dadd_key=="<< veb_dadd_key<<std::endl;
									}
								}
#endif
								break;
							default: root = insert_and_modify(root, dadd, dsearch);
								break;
							}




						}


					}

				}
			}


			

			// Очистка (восстановление хеш-таблицы).
			// НИ в коем случае не параллелить по OPENMP в этом месте.!!!
			for (integer i_54 = 0; i_54 < itop_stack2; ++i_54) {
				hash_table2[istack[i_54]] = false;
			}
			itop_stack2 = 0; // стек снова готов к работе.





			if (set != nullptr) {
				delete[] set;
				set = nullptr;
			}

			newCcount++;
			// Один агрегат создан.

		} // узел не был ещё включён в агрегат.



		bcontinue_gl_1 = false;
		for (integer i_1 = istartflag_scan; i_1 <= n_a[ilevel - 1]; ++i_1) {
			if (flag[i_1] == false) {
				bcontinue_gl_1 = true;
				istartflag_scan = i_1; // сокращаем пределы сканирования.
				break; // досрочный выход из цикла for.
			}
		}

		// Вычисление узла с максимальным количеством соседей.
		icandidate = 0;


		// Данный код чрезвычайно компактен.
		// Надо найти максимальный элемент в АВЛ дереве.
		node_AVL* emax = 0;
		Tree_splay* emax_splay = 0;
		TreapNode* emax_random_tree = nullptr;
		TreapNode* save_root = nullptr;
		data_BalTree dbt_emax;

		integer ui_emax;

		switch (id_tree)
		{
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE: emax = findmax(root);
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE: emax_splay = findmax(root_splay);
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP:
			if (!binary_heap.empty()) {
				// Куча не пуста.
				icandidate = row_startA[binary_heap.readkeymaxelm()];
			}
			else {
				size_splay_Tree = 0;
				icandidate = 0;
				bcontinue_gl_1 = false;
			}
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP:
			save_root = random_tree_root;
			if (emax_random_tree != nullptr) {
				delete[] emax_random_tree;
				emax_random_tree = nullptr;
			}
			emax_random_tree = findmax_random_tree(random_tree_root);
			random_tree_root = save_root;
			save_root = nullptr;

			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE:
			dbt_emax = RBroot.GetMaxElm();
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP:
			if (fibo_heap.isEmpty()) {
				dbt_emax.i = -1;
			}
			else {
				ui_emax = -fibo_heap.getMinimum();
				dbt_emax.i = ((ui_emax) % (n_a[ilevel - 1] + 1));
				dbt_emax.count_neighbour = ((ui_emax) / (n_a[ilevel - 1] + 1));
			}
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE:
#if VEB_FLAG
			if (!((vanEMDE_BOAS_Tree == nullptr) || ((vanEMDE_BOAS_Tree->summary == nullptr) && (vanEMDE_BOAS_Tree->cluster == nullptr)))) {
				vEB_max(vanEMDE_BOAS_Tree, ui_emax);
				if (ui_emax <= 0) {
					// дерево ван Эмде Боаса пустое.
					dbt_emax.i = -1;
					dbt_emax.count_neighbour = -1;
				}
				else {
					dbt_emax.i = ((ui_emax) % (n_a[ilevel - 1] + 1));
					dbt_emax.count_neighbour = ((ui_emax) / (n_a[ilevel - 1] + 1));
				}

			}
			else {
				// дерево ван Эмде Боаса пустое.
				dbt_emax.i = -1;
				dbt_emax.count_neighbour = -1;
			}
#endif
			break;
		default: emax = findmax(root);
			break;
		}


		switch (id_tree) {
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE:
			// AVL tree
			if (emax != 0) {

				//icandidate = emax->key.ii; 23 jan 2016
				icandidate = row_startA[emax->key.i];
				emax = 0;
			}
			else {
				root = 0;
				icandidate = 0;
				bcontinue_gl_1 = false;

			}
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE:
			// SPLAY tree
			if (emax_splay != 0) {


				//icandidate = emax_splay->item.ii; 23 jan 2016
				icandidate = row_startA[emax_splay->item.i];
				emax_splay = 0;

			}
			else {
				RBroot.Clear();
				root_splay = 0;
				size_splay_Tree = 0;
				random_tree_root = nullptr;
				icandidate = 0;
				bcontinue_gl_1 = false;

			}
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP:
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP:
			// Random TREE
			if (emax_random_tree != nullptr) {
				icandidate = row_startA[emax_random_tree->key.i];
				delete emax_random_tree;
				emax_random_tree = nullptr;
			}
			else {
				RBroot.Clear();
				root_splay = 0;
				size_splay_Tree = 0;
				random_tree_root = nullptr;
				icandidate = 0;
				bcontinue_gl_1 = false;
			}
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE:
			if (RBroot.Find(dbt_emax)) {
				icandidate = row_startA[dbt_emax.i];
			}
			else {
				RBroot.Clear();
				root_splay = 0;
				size_splay_Tree = 0;
				random_tree_root = nullptr;
				icandidate = 0;
				bcontinue_gl_1 = false;
			}
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::FIBONACCI_HEAP:
			if (dbt_emax.i == -1)
			{
				// Дерево пусто.
				RBroot.Clear();
				root_splay = 0;
				size_splay_Tree = 0;
				random_tree_root = nullptr;
				icandidate = 0;
				bcontinue_gl_1 = false;
			}
			else {
				// искомый узел и дерево ван Эмде Боаса не пусто.
				if (!fibo_heap.isEmpty()) {
					fibo_heap.removeMinimum();
				}
				icandidate = row_startA[dbt_emax.i];
				//std::cout << "row_startA = " << icandidate<< " " << dbt_emax.i << std::endl;

			}
			break;
		case RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE:
			if (dbt_emax.i == -1)
			{
				// Дерево пусто.
				RBroot.Clear();
				root_splay = 0;
				size_splay_Tree = 0;
				random_tree_root = nullptr;
				icandidate = 0;
				bcontinue_gl_1 = false;
			}
			else {
				// искомый узел и дерево ван Эмде Боаса не пусто.
				icandidate = row_startA[dbt_emax.i];
				//std::cout << "row_startA = "<< icandidate << std::endl;
			}
			break;
		default:
			// AVL tree
			if (emax != 0) {

				//icandidate = emax->key.ii; 23 jan 2016
				icandidate = row_startA[emax->key.i];
				emax = 0;
			}
			else {
				root = 0;
				icandidate = 0;
				bcontinue_gl_1 = false;

			}
			break;
		}



		if (debug_reshime) system("pause");
		//system("pause");

		if (icandidate == 0) {
			bcontinue_gl_1 = false;
		}
		// 4 june 2016

		icountprohod++;

	} // Построение C/F разбиения. создано.

	delete[]  max_vnediagonal33_arr;

// Освобождение оперативной памяти из под АВЛ дерева.
			// 12 декабря 2015.
	switch (id_tree)
	{
	case RS_COARSENING_KERNEL_DATA_STRUCTURE::AVL_TREE: clear_AVL(root);
		root = 0;
		break;
	case RS_COARSENING_KERNEL_DATA_STRUCTURE::SPLAY_TREE:
		clear_SPLAY(root_splay);
		root_splay = 0;
		break;
	case RS_COARSENING_KERNEL_DATA_STRUCTURE::BINARY_HEAP:
		binary_heap.clear();
		break;
	case RS_COARSENING_KERNEL_DATA_STRUCTURE::TREAP:
		clear_random_tree(random_tree_root);
		random_tree_root = nullptr;
		break;
	case RS_COARSENING_KERNEL_DATA_STRUCTURE::RED_BLACK_TREE:
		RBroot.Clear();
		break;
	case RS_COARSENING_KERNEL_DATA_STRUCTURE::VAN_EMDE_BOAS_TREE:
#if VEB_FLAG
		if (!((vanEMDE_BOAS_Tree == nullptr) || ((vanEMDE_BOAS_Tree->summary == nullptr) && (vanEMDE_BOAS_Tree->cluster == nullptr)))) {
			vanEMDE_BOAS_Tree->~TvEB();
		}
#endif
		break;
	default: clear_AVL(root);
		root = 0;
		break;
	}

	//delete[] bmarkervisit;
	if (bmarkervisit != nullptr) {
		free(bmarkervisit);
		bmarkervisit = nullptr;
	}


} // Ruge_and_Stuben_CF_decomposition

// comparator function to make min heap 
struct greaters_Stuben {
	bool operator()(const std::pair<integer, integer>& a, const std::pair<integer, integer>& b) const {
		return a.second < b.second;
	}
};

#include <map>    //подключили библиотеку для работы с map
//#include <iomanip>

// Данная реализация корректна но непригодна из-за очень медленной скорости
// работы функции  std::max_element(m1.begin(), m1.end(), greaters_Stuben());
// По видимому метод std::max_element() работает очень неэффективно за линейное время.
// Данная  функция не рекомендуется к использованию из-за очень низкой производительности.
// Используйте для трёхмерных задач функцию PMIS_CF_decomposition() реализующую алгоритм 
// огрубления PMIS (Parallel Modified Independed Set) или функцию 
// Ruge_and_Stuben_CF_decomposition() для Руге-Штубена огрубления на основе быстродействующего 
// АТД (Абстрактного типа данных) Фибоначчиевой кучи.
// Разбиение множества узлов первоначальной сетки на С-узлы и 
// F - узлы. С -узлы составят грубую сетку следующего уровня вложенности.
// Значения же функции в F узлах должно быть восстановлено по значению функции
// в ближайших С узлах (см. задачу интерполяции). 
template <typename doublerealT>
void Ruge_and_Stuben_CF_decomposition_std(Ak2& Amat, bool*& this_is_F_node,
	bool*& this_is_C_node, integer ilevel,
	integer*& count_neighbour, integer*& n_a,
	integer& newCcount,
	doublerealT*& threshold_quick_only_negative,
	integer*& row_startA,
	integer*& hash_StrongTranspose_collection1Eco,
	integer* & isize_hash_StrongTranspose_collection,
	bool bprint_mesage_diagnostic,
	bool*& flag, integer iadd,
	bool bpositive_connections_CF_decomp,
	bool bStrongTransposeON,
	bool*& hash_table2,
	integer*& istack, integer*& nnz_a,
	bool debug_reshime)
{



	integer ii_end1 = n_a[ilevel - 1];
	integer max_neighbour = 0;
	integer icandidate = 0;

	// Находим узел с наибольшим числом соседей и запоминаем его.
	// Это первый встретившийся узел с наибольшим числом соседей.
	// Это требуется для того чтобы стартовал алгоритм C/F разбиения.

	for (integer i7 = 1; i7 <= ii_end1; ++i7) {
		if (count_neighbour[i7] > max_neighbour) {
			max_neighbour = count_neighbour[i7];
			icandidate = row_startA[i7];
		}
	}




	// нужно выделить кто попал в coarse, а кто в этот раз попал в Fine Выделить всех кто соседствует
	// с новыми Fine увеличить им счётчик соседей.


	const integer NULL_NEIGHBOUR = -1;
	// Построение C/F разбиения.
	integer vacant = NULL_NEIGHBOUR;
	bool bcontinue_gl_1 = true;
	// Построение C/F разбиения.
	integer icountprohod = 0;

	// храним те узлы которые уже были пройдены при конструировании.
	// поначалу все узлы помечены как непосещённые.
	bool* bmarkervisit = nullptr;
	if (bmarkervisit != nullptr) {
		free(bmarkervisit);
		bmarkervisit = nullptr;
	}
	bmarkervisit = my_declaration_array<bool>(n_a[ilevel - 1], false, "bmarkervisit");


	// увеличение быстродействия достигается 
	// сокращением пределов сканирования
	// здесь хранится индекс начала сканирования flag.
	integer istartflag_scan = 1;

	// Задача 12mm hfet thermal resistance. 1.7млн неизвестных.
	// AVL_TREE_ID   3мин 29с 590мс      {5}
	// SPLAY_TREE_ID  3мин 16с 430мс {2}
	// BINARY_HEAP 3мин 4с 0мс {1 *самая быстрая.}
	// RANDOM_TREE_ID (Дерамида) 3мин 28с 90мс {4}
	// RED_BLACK_TREE_ID 3мин 27с 210мс {3}

	integer n = n_a[0];

	
	std::map<integer, integer> m1;
	//std::vector<std::pair<integer, integer>> v1;
	
	newCcount = 0;

	// Нехорошо постоянно выделять и уничтожать память в длинном цикле, 
	// более быстро выделить её один раз. См. выделение памяти под set.
	// 23.04.2017


	if (bprint_mesage_diagnostic) {
		std::cout << "   ***   CAMG SELECTOR " << ilevel <<"  ***\n";
	}
	while (bcontinue_gl_1)
	{

		integer ic = 0;
		integer ic_end_F_SiTranspose = 0;

		integer ii = icandidate;
		if (flag[Amat.i[ii]] == false) {

			ic = 0; // Обязательная инициализация.


			ic_end_F_SiTranspose = 0;
			integer set0 = Amat.i[ii];




			//A20.05.2017//this_is_C_node[set[0]] = true;
			//A20.05.2017//bmarkervisit[set[0]] = true;
			this_is_C_node[set0] = true;
			bmarkervisit[set0] = true;

			doublerealT max_vnediagonal = -1.0; // максимальное значение модуля вне диагонального элемента. 
												// добавляем диагональный элемент.
												// узел set[0]==Amat.i[is0].
												// Нахождение значения максимального внедиагонального элемента, с 
												// учётом того что даже узел Дирихле связан с одним внутренним узлом расчётной области.
												// 17 января 2016 правильное определение максимального внедиагонального элемента.
												// Обязательная перемотка в самое начало строки.
			integer ii_back = ii;
			while ((ii_back > iadd) && (Amat.i[ii_back] == set0)) ii_back--;
			ii_back++;

			doublerealT max_vnediagonal1 = -1.0e30;
						

			// Если делать по максимальному внедиагональному элементу то мы получим очень много элементов на грубых уровнях,
			// и чрезвычайно медленную сходимость.

			if (bpositive_connections_CF_decomp) {
				// 23_10_2016
				for (integer is0 = ii_back; (is0 <= row_startA[set0 + 1] - 1); ++is0) {
					if (Amat.j[is0] != set0) {
						
						if (Amat.abs_aij[is0] > max_vnediagonal1) {
							max_vnediagonal1 = Amat.abs_aij[is0]; //i,j
							// Большое количество элементов на грубых уровнях,
							// очень медленная сходимость.
							//if (Amat.j[is0] == set[0]) break; 
						}						
					}
				}
			}
			else {
				for (integer is0 = ii_back; (is0 <= row_startA[set0 + 1] - 1); ++is0) {
					if (Amat.j[is0] != set0) {
						if (Amat.aij[is0] < 0.0) {
							
							if (Amat.abs_aij[is0] > max_vnediagonal1) {
								max_vnediagonal1 = Amat.abs_aij[is0]; //i,j
								// Большое количество элементов на грубых уровнях,
								// очень медленная сходимость.
								//if (Amat.j[is0] == set[0]) break; 
							}							
						}
					}
				}
			}
			
			
			//max_vnediagonal = max_vnediagonal1;  // 1			
			// наиболее близка к оптимальной. -85%. но несомненно лучше max_vnediagonal = -1.0;
			// 19 января 2016 установлено что важны все связи, не нужно учитывать threshold
			// max_vnediagonal должно быть -1.0. Именно это значение обеспечивает наилучшую 
			// скорость агломерации и наилучшую скорость сходимости.
			max_vnediagonal = -1.0e30;  // все связи!!!											

			ic++;


			//  В set начиная с единицы и до <ic лежат кандидаты чтобы стать F.
			// 5.01.2017
			// 01.04.2017 Дополняемся F узлами из Si_Transpose связей.
			if ((my_amg_manager.ipatch_number == 7) && (bStrongTransposeON)) {

				integer imarker75_scan = 0;

				// обычный линейный список.
				/*formirate_F_SiTranspose_hash_table_Gus2_struct03(hash_StrongTranspose_collection1[Amat.i[ii]],
					isize_hash_StrongTranspose_collection[Amat.i[ii]], imarker75_scan, this_is_F_node, this_is_C_node);*/

				// 10.06.2021
				// На основе таблицы значений в массиве.
				formirate_F_SiTranspose_hash_table_Gus2_struct04(hash_StrongTranspose_collection1Eco, Amat.i[ii], n_a[ilevel - 1],
					isize_hash_StrongTranspose_collection[Amat.i[ii]], imarker75_scan, this_is_F_node, this_is_C_node);

				ic = imarker75_scan + 1;
			}

			ic_end_F_SiTranspose = ic; // С этой позиции заканчиваются F которые из Si_Transpose.

									   // если узел j ещё не был добавлен в агрегат.
			if (bpositive_connections_CF_decomp) {
				if (flag[Amat.j[ii]] == false) {
					if ((Amat.j[ii] != set0) && (Amat.abs_aij[ii] >= theta(ilevel) * max_vnediagonal)) {
						// 21.05.2017
						bool bfound_vacant = false;

						bfound_vacant = isfound_hash_table_Gus_struct01(Amat.j[ii]);
						if (!bfound_vacant) {
							insert_hash_table_Gus_struct01(Amat.j[ii]);
							ic++;
						}

					}
				}
			}
			else {
				if (flag[Amat.j[ii]] == false) {
					if ((Amat.j[ii] != set0) && (Amat.aij[ii] < 0.0) && (Amat.abs_aij[ii] >= theta(ilevel) * max_vnediagonal)) {
						// 21.05.2017
						bool bfound_vacant = false;

						bfound_vacant = isfound_hash_table_Gus_struct01(Amat.j[ii]);
						if (!bfound_vacant) {
							insert_hash_table_Gus_struct01(Amat.j[ii]);
							ic++;
						}

					}
				}
			}

			//std::cout << "sboi start";

			// iscan = ii+1; // устаревший код
			integer iscan = ii_back + 1; // важная модификация 19 января 2016г.
			// 19 jan 2016.

			if (bpositive_connections_CF_decomp) {
				while ((iscan <= nnz_a[ilevel - 1] + iadd) && (Amat.i[iscan] == set0)) {
					// 14 февраля 2016 код иногда приводящий к сбою.
					//while (iscan <= row_startA[set0 + 1] - 1) { // код иногда приводящий к сбою по непонятной причине.
					// если узел j ещё не был добавлен в агрегат.
					if (flag[Amat.j[iscan]] == false) {
						if ((Amat.j[iscan] != set0) && (Amat.abs_aij[iscan] >= theta(ilevel) * max_vnediagonal)) {
							// 21.05.2017
							bool bfound_vacant = false;

							bfound_vacant = isfound_hash_table_Gus_struct01(Amat.j[iscan]);
							if (!bfound_vacant) {
								insert_hash_table_Gus_struct01(Amat.j[iscan]);
								ic++;
							}

							/*
							// Медленная версия с линейным поиском.
							vacant = Amat.j[iscan];
							for (integer js = 0; js < ic; ++js) {
							if (vacant == set[js]) {
							vacant = NULL_NEIGHBOUR;
							}
							}
							if (vacant != NULL_NEIGHBOUR) {
							set[ic] = vacant;

							ic++;

							}
							*/
						}
					}

					iscan++;

				} // while
			}
			else {
				while ((iscan <= nnz_a[ilevel - 1] + iadd) && (Amat.i[iscan] == set0)) {
					// 14 февраля 2016 код иногда приводящий к сбою.
					//while (iscan <= row_startA[set0 + 1] - 1) { // код иногда приводящий к сбою по непонятной причине.
					// если узел j ещё не был добавлен в агрегат.
					if (flag[Amat.j[iscan]] == false) {
						if ((Amat.j[iscan] != set0) && (Amat.aij[iscan] < 0.0) && (Amat.abs_aij[iscan] >= theta(ilevel) * max_vnediagonal)) {
							// 21.05.2017
							bool bfound_vacant = false;

							bfound_vacant = isfound_hash_table_Gus_struct01(Amat.j[iscan]);
							if (!bfound_vacant) {
								insert_hash_table_Gus_struct01(Amat.j[iscan]);
								ic++;
							}

							/*
							vacant = Amat.j[iscan];
							for (integer js = 0; js < ic; ++js) {
							if (vacant == set[js]) {
							vacant = NULL_NEIGHBOUR;
							}
							}
							if (vacant != NULL_NEIGHBOUR) {
							set[ic] = vacant;

							ic++;

							}
							*/
						}
					}

					iscan++;

				} // while
			}

			//std::cout << "sboi end";
			// Это была учтена только связь i,j






		// В этом месте множество set успешно сформировано:
		// 1. Перепаковка из root_Gus_set в set.
		// 2. root_Gus_set больше не используется.
		// 3. Именно здесь надо выделить данные под set.
			integer* set = nullptr;
			set = new integer[ic + 2];
			//if (set == nullptr) {
				//std::cout << "error!!! memory for set is nullptr. Problem allocate detected."<<std::endl;
				//std::cout << "in function classic_aglomerative_amg6."<<std::endl;
				//system("pause");
				//exit(1);
			//}

			integer ic_986 = 1;
			set[0] = set0;



			formirate_hash_table_Gus_struct01__2__set(set, ic_986);

			clear_hash_table_Gus_struct01();


			for (integer isc = 1; isc < ic; ++isc) {
				this_is_F_node[set[isc]] = true; // это только новые F узлы.
				bmarkervisit[set[isc]] = true;
			}

			// Помечаем узлы как включённые в агрегат.
			for (integer js = 0; js < ic; ++js) {
				flag[set[js]] = true;
			}

			// Алгоритм (5 декабря 2015 revised) 
			// 1. Сканируем все F которые соседи данного С на данном проходе.
			// 2. Для каждого фиксированного F сканируем его "строчных" соседей.
			// 3. Если узел еще не был включён в агрегат то ищем всех соседей данного узла на предмет 
			// соседства с фиксированным набором F из пункта 1.

			
			// 12 декабря 2015.
			// Надо удалить из АВЛ дерева C и F узлы.
			// Это удаление очищает АВЛ дерево и приводит его к
			// рабочему состоянию. Удаление несуществующих в дереве узлов
			// производится корректно. Удаление производится за логарифмическое
			// по основанию 2  время от количества элементов в дереве
			// сбалансированность дерева при этом сохраняется.
			for (integer js = 0; js < ic; ++js) {
				
				if (!m1.empty()) {
					m1.erase(set[js]);
					//v1.erase(std::pair<integer, integer>(set[js], count_neighbour[set[js]]));
					//v1.erase(set[js]);
				}		
			}




			//std::cout<<"additional and modify new neighbour"<<std::endl;

			// 10 января 2016. Новая логика.
			// Устраним некоторые повторные модификации (это должно снизить нагрузку на АВЛ дерево).
			// Эта модификация даёт сокращение количества V циклов которые требуются до сходимости
			// Эта модификация наиболее близка к классической описанной в литературе чем все предыдущие.
			// На момент 13 января 2016 это лучший вариант по скорости вычислений.
			integer itop_stack2 = 0;

			// 10 января 2016. Старый вариант просто очищенный от устаревшего кода.
			for (integer js = 1; js < ic; ++js) {

				integer i_11 = set[js];
				integer ii_11 = row_startA[i_11];
				//integer iend5 = nnz_a[ilevel - 1] + iadd;
				//integer istart73 = ii_11;
				//while ((istart73 >= 1 + iadd) && (Amat.i[istart73] == Amat.i[ii_11])) istart73--;
				//istart73++;
				integer istart73 = row_startA[Amat.i[ii_11]];
				integer iend73 = row_startA[Amat.i[ii_11] + 1] - 1;
				bool bvisitsos = false;
				for (integer is0 = istart73; (is0 <= iend73); ++is0) {
					//for (integer is0 = istart73; (is0 <= iend5) && (Amat.i[is0] == Amat.i[ii_11]); ++is0) {
					// В пересечении с U!!!
					if (flag[Amat.j[is0]] == false) {


						integer isc = Amat.j[is0];



						// Избавляемся от повторных инкрементаций.
						// В 2D на пятиточечном шаблоне повторные инкрементации составляют
						// около 33%.
						// Это даёт стандартный алгоритм сгрубления, описаный в статьях, но
						// на ряде тестовых задач при таком подходе агломерация проходила очень
						// плохо (переполнение по памяти, не хватало даже семикратного размера исходной матрицы).
						// Эта проблема проявилась на задачах:
						// CGHV1J006D, Потенциал тора, Электрический потенциал в FET, Module 2.
						// Плохая скорость агломерации получается главным образом из-за шестого способа интерполяции.
						// Проблема не в этом месте кода.
						if (hash_table2[isc] == false) {
							hash_table2[isc] = true;
							istack[itop_stack2] = isc;
							itop_stack2++;
							// закомментированный лучше.
							//}

							//21_12_2016
							integer ii_2 = row_startA[isc];


							integer ic2 = 0;
							integer iend2loc = nnz_a[ilevel - 1] + iadd;
							//integer istart72 = ii_2;										
							integer istart72 = row_startA[Amat.i[ii_2]];
							integer istopmarker2 = row_startA[Amat.i[ii_2] + 1] - 1;

							// 22 _12_2016
							// Это лучший вариант: обеспечивает корректное построение иерархии
							// уровней на задаче passive module 6 в то время как все остальные 
							// отличные от этого способа давали сбой.
							doublerealT max_vnediagonal33 = -1.0e30;
							for (integer is01 = istart72; (is01 <= istopmarker2); ++is01) {
								if (Amat.j[is01] != Amat.i[is01]) {
									if ((Amat.aij[is01] < 0.0) && (Amat.abs_aij[is01] > max_vnediagonal33)) {
										max_vnediagonal33 = Amat.abs_aij[is01];
									}
								}
							}
							for (integer is01 = istart72; (is01 <= istopmarker2); ++is01) {
								// 0.2375 импирически подобрана для passive module 6.
								if ((Amat.aij[is01] < 0.0) && (Amat.abs_aij[is01] > 0.2375 * max_vnediagonal33)) {
									if (Amat.j[is01] == set[js]) {
										if ((my_amg_manager.ipatch_number == 7) && (bStrongTransposeON)) {
											if (js < ic_end_F_SiTranspose) {
												// Увеличиваем счётчики только тех соседей F узлов которые
												// являются соседями F узлов которые были получены из Si_Transpose связей.
												// Именно так написано у Джона Руге и Клауса Штубена.
												ic2++;
											}
										}
										else {
											ic2++;
										}
									}
								}
								// уменьшить счетчик слабого (weakly) соседа ?
							}


							data_BalTree dsearch;
							dsearch.count_neighbour = count_neighbour[isc];
							//dsearch.ii = ii_2;
							dsearch.i = isc;
							// Увеличиваем на количество связей с новыми F узлами.
							count_neighbour[isc] += ic2;
							data_BalTree dadd;
							dadd.count_neighbour = count_neighbour[isc];
							//dadd.ii = ii_2;
							dadd.i = isc;

							std::map<integer, integer>::const_iterator it; // объявляем итератор
							//it = m1.begin(); // присваиваем ему начало map-а
							it = m1.find(isc);
							if (it == m1.end()) {
								// значения в map нету.
								m1.insert(std::pair<integer, integer>(isc, count_neighbour[isc]));
								//v1.push_back(std::pair<integer, integer>(isc, count_neighbour[isc]));
							}
							else {
								// значение в map присутствует.
								m1.erase(isc);
								//v1.erase(std::pair<integer, integer>(isc, dsearch.count_neighbour));
								
								m1.insert(std::pair<integer, integer>(isc, count_neighbour[isc]));
								//v1.push_back(std::pair<integer, integer>(isc, count_neighbour[isc]));
							}							

						}


					}

				}
			}

			// Очистка (восстановление хеш-таблицы).
			// НИ в коем случае не параллелить по OPENMP в этом месте.!!!
			for (integer i_54 = 0; i_54 < itop_stack2; ++i_54) {
				hash_table2[istack[i_54]] = false;
			}
			itop_stack2 = 0; // стек снова готов к работе.





			if (set != nullptr) {
				delete[] set;
				set = nullptr;
			}

			newCcount++;
			// Один агрегат создан.

		} // узел не был ещё включён в агрегат.



		bcontinue_gl_1 = false;
		for (integer i_1 = istartflag_scan; i_1 <= n_a[ilevel - 1]; ++i_1) {
			if (flag[i_1] == false) {
				bcontinue_gl_1 = true;
				istartflag_scan = i_1; // сокращаем пределы сканирования.
				break; // досрочный выход из цикла for.
			}
		}

		// Вычисление узла с максимальным количеством соседей.
		icandidate = 0;


		// Данный код чрезвычайно компактен.
		// Надо найти максимальный элемент в АВЛ дереве.
		

		if (!m1.empty()) {
			std::map<integer, integer>::const_iterator it; // объявляем итератор
			// Это очень медленный код.
			it = std::max_element(m1.begin(), m1.end(), greaters_Stuben());
			//std::cout << (*it).first << "= ";

			

			//std::make_heap(v1.begin(), v1.end(), greaters_Stuben());
			//std::pop_heap(v1.begin(), v1.end(), greaters_Stuben()); // удалить максимальный элемент из кучи.
			

			icandidate = row_startA[(*it).first];
			//icandidate = row_startA[(v1.back()).first];  // посмотреть максимальный элемент.
		}
		else {
			icandidate = 0;
			bcontinue_gl_1 = false;
		}
		


		if (debug_reshime) system("pause");
		//system("pause");

		if ((icandidate == 0) ) {
			bcontinue_gl_1 = false;
		}
		// 4 june 2016

		icountprohod++;

	} // Построение C/F разбиения. создано.


    // Освобождение оперативной памяти из под АВЛ дерева.
	// 12 декабря 2015.
	m1.clear();
	//v1.clear();

	//delete[] bmarkervisit;
	if (bmarkervisit != nullptr) {
		free(bmarkervisit);
		bmarkervisit = nullptr;
	}


} // Ruge_and_Stuben_CF_decomposition_std


// Разбиение множества узлов первоначальной сетки на С-узлы и 
// F - узлы. С -узлы составят грубую сетку следующего уровня вложенности.
// Значения же функции в F узлах должно быть восстановлено по значению функции
// в ближайших С узлах (см. задачу интерполяции). 
// Начало написания 20.09.2020. PMIS применённый к квадрату матрицы А.
// Распараллелил 08.10.2020
template <typename doublerealT>
void PMIS_CF_decomposition_applied_to_the_square_of_the_matrix(Ak2& Amat, bool*& this_is_F_node,
	bool*& this_is_C_node, integer ilevel,
	integer*& count_neighbour, integer*& n_a,
	integer& newCcount,
	doublerealT*& threshold_quick_only_negative,
	integer*& row_startA,
	integer*& hash_StrongTranspose_collection1Eco,
	integer* &isize_hash_StrongTranspose_collection)
{

	bool bprint_mesage_diagnostic = true;
	if (my_amg_manager.iprint_log == 0) {
		bprint_mesage_diagnostic = false;
	}

#ifdef _OPENMP
	int inum_core = number_cores();
	int i_my_num_core_parallelesation = omp_get_max_threads();
	omp_set_num_threads(inum_core); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
#endif

	integer ii_end1 = n_a[ilevel - 1];

	doublerealT* dcount_neighbour = new doublerealT[n_a[ilevel - 1] + 1];

	integer inumber_isolated_F_nodes = 0;

	//std::mt19937 gen(time(0));
	std::mt19937 gen{ std::random_device()() };
	std::uniform_real_distribution<> urd(0, 1);

	// urd(gen) не работает в многопотоке!!!
//#pragma omp parallel for reduction(+ : inumber_isolated_F_nodes)
	for (integer i7 = 1; i7 <= ii_end1; ++i7) {
		// число соседей плюс случайное число от нуля до единицы.
		if (count_neighbour[i7] == 0) {
			this_is_F_node[i7] = true;
			inumber_isolated_F_nodes++;
		}
		// drand случайное вещественное число от нуля до единицы.
		//-->doublerealT drand = static_cast<doublerealT>(((double)(1.0 * rand()) / ((double)(RAND_MAX + 1))));
		doublerealT drand = static_cast<doublerealT>(urd(gen));
		dcount_neighbour[i7] = count_neighbour[i7] + drand;
	}


	bool bflag_empty = false;
	bool bcontinue = true;
	while (bcontinue) {
		bcontinue = false;
		integer inumber_of_nodes_viewed = 0; // Число просмотренных узлов на новом проходе (сканировании).
		integer inumber_of_C_nodes = 0; // Число добавленных агрегатов.

#pragma omp parallel for reduction(+:inumber_of_nodes_viewed, inumber_of_C_nodes)
		for (integer i7 = 1; i7 <= ii_end1; ++i7) {
			bool cnd = ((this_is_C_node[i7] == false) && (this_is_F_node[i7] == false));
			if (!cnd) continue;

			
			inumber_of_nodes_viewed++;

			if (bflag_empty) {
				// Непонятные ошмётки, их небольшое количество делаем С узлами.
				inumber_of_C_nodes++;
				newCcount++;
				this_is_C_node[i7] = true;
				//std::cout<<"lambdai="<< count_neighbour[i7]<<std::endl;
				//system("PAUSE");
			}
			else {

				bcontinue = true;

				doublerealT id_diag = dcount_neighbour[i7];

				// Определение S(i)
				doublerealT max_vnediagonal = -1.0;
				//for (integer is0 = row_startA[i7]; (is0 <= row_startA[i7 + 1] - 1); ++is0) {
					//if (Amat.j[is0] != i7) {
						//if ((Amat.aij[is0]<0.0)&&(Amat.abs_aij[is0] > max_vnediagonal)) {
							//max_vnediagonal = Amat.abs_aij[is0];
						//}
					//}
				//}
				max_vnediagonal = threshold_quick_only_negative[i7];
				//  Забыли учесть transpose(S(i)).
				bool bunion_on = true;
				//bool bOk = false;
				integer is0_end = row_startA[i7 + 1] - 1;
				for (integer is0 = row_startA[i7]; (is0 <= is0_end); ++is0) {
					bool cnd1 = (Amat.j[is0] != i7);
					if (!cnd1) continue;
					bool cnd2 = ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > max_vnediagonal));
					if (!cnd2) continue;
					
						

					// Сосед еще не обрабатывался.
					if ((this_is_C_node[Amat.j[is0]] == false) && (this_is_F_node[Amat.j[is0]] == false)) {

						// Соседи соседей.
						integer i8 = Amat.j[is0];
						doublerealT max_vnediagonal2 = threshold_quick_only_negative[i8];
						integer is0_end2 = row_startA[i8 + 1] - 1;
						for (integer is02 = row_startA[i8]; (is02 <= is0_end2); ++is02) {
							if ((Amat.j[is02] != i7)&&(Amat.j[is02] != i8)) {
								if ((Amat.aij[is02] < 0.0) && (Amat.abs_aij[is02] > max_vnediagonal2)) {
									if (dcount_neighbour[Amat.j[is02]] > id_diag) {
										// Сосед еще не обрабатывался.
										if ((this_is_C_node[Amat.j[is02]] == false) && (this_is_F_node[Amat.j[is02]] == false)) {
											bunion_on = false;
										}
									}
								}
							}
						}


						if (dcount_neighbour[Amat.j[is0]] > id_diag) {									

							bunion_on = false;
						}
						//else {
							// Сосед еще не обрабатывался.
							//if ((this_is_C_node[Amat.j[is0]] == false) && (this_is_F_node[Amat.j[is0]] == false)) {
								//	bOk = true;
							//	}
						//}
					}						
					
				}
				// учитываем transpose(S(i)):
				//if (hash_StrongTranspose_collection1 != nullptr)
				{
					//if (hash_StrongTranspose_collection1[i7] != nullptr)
					{
						//Taccumulqtor_list* list_scan = hash_StrongTranspose_collection1[i7];
						integer i71 = 0;
						integer isize71 = isize_hash_StrongTranspose_collection[i7];
						while (i71<isize71) {

							//integer icandidate73 = list_scan->ikey;
							integer icandidate73 = hash_StrongTranspose_collection1Eco[i71 *n_a[ilevel-1]+i7];


							if (icandidate73 != i7) {
								
								// Сосед еще не обрабатывался.
								if ((this_is_C_node[icandidate73] == false) && (this_is_F_node[icandidate73] == false)) {

									// соседи соседей.
									//if (hash_StrongTranspose_collection1[icandidate73] != nullptr)
									{
										//Taccumulqtor_list* list_scan2 = hash_StrongTranspose_collection1[icandidate73];
										integer i72 = 0;
										integer isize72 = isize_hash_StrongTranspose_collection[icandidate73];

										while (i72<isize72) {
											//integer icandidate734 = list_scan2->ikey;

											integer icandidate734 = hash_StrongTranspose_collection1Eco[i72 * n_a[ilevel - 1] + icandidate73];

											if ((icandidate734 != i7)&&(icandidate734 != icandidate73)) {
												if (dcount_neighbour[icandidate734] > id_diag) {
													// Сосед соседей еще не обрабатывался.
													if ((this_is_C_node[icandidate734] == false) && (this_is_F_node[icandidate734] == false)) {
														bunion_on = false;
													}
												}
											}
											//list_scan2 = list_scan2->next;
											i72++;
										}
									}


									if (dcount_neighbour[icandidate73] > id_diag) {											

										bunion_on = false;
									}
									//else {
									//	if ((this_is_C_node[icandidate73] == false) && (this_is_F_node[icandidate73] == false)) {
									//		bOk = true;
									//	}
									//}
								}
							}
							//list_scan = list_scan->next;
							i71++;
						}
					}
				}
				// если bunion_on то на диагонали самая сильная связь и мы создаём агрегат.
				// Сюда в агрегат входят связи из S(i).
				/*
				if (bunion_on) {
					inumber_of_C_nodes++;
					this_is_C_node[i7] = true;
					// Здесь соседей брать только из transpose(S(i)).
					for (integer is0 = row_startA[i7]; (is0 <= row_startA[i7 + 1] - 1); ++is0) {
						// Сосед еще не обрабатывался.
						if ((this_is_C_node[Amat.j[is0]] == false) && (this_is_F_node[Amat.j[is0]] == false)) {
							if ((Amat.j[is0] != i7) && (Amat.abs_aij[is0] >= theta(ilevel) * max_vnediagonal)) {
								this_is_F_node[Amat.j[is0]] = true;
							}
						}
					}
					newCcount++;
					// Один агрегат создан.
				}
				*/
				// если bunion_on то на диагонали самая сильная связь и мы создаём агрегат.
				// Сюда в агрегат входят связи из transpose(S(i)).
				//if (true||bOk) 
				{
					if (bunion_on) {
						inumber_of_C_nodes++;
						this_is_C_node[i7] = true;

						// Увеличивает число итераций с 899 до 998.
						// Не уменьшая разреженность оператора.
						// Вообще всех соединенных с i7 делаем F узлами.
						//integer is0_end3 = row_startA[i7 + 1] - 1;
						//for (integer is02 = row_startA[i7]; (is02 <= is0_end3); ++is02) {
							//if (Amat.j[is02] != i7)  {
								//not//if ((Amat.aij[is02] < 0.0) && (Amat.abs_aij[is02] > max_vnediagonal2)) 
								//{

									// Сосед еще не обрабатывался.
									//if ((this_is_C_node[Amat.j[is02]] == false) && (this_is_F_node[Amat.j[is02]] == false)) {
										//this_is_F_node[Amat.j[is02]] = true;
									//}

								//}
							//}
						//}


						// Здесь соседей брать только из transpose(S(i)).
						//if (hash_StrongTranspose_collection1 != nullptr) 
						{
							//if (hash_StrongTranspose_collection1[i7] != nullptr)
							{
								//Taccumulqtor_list* list_scan = hash_StrongTranspose_collection1[i7];
								integer i71 = 0;
								integer isize71 = isize_hash_StrongTranspose_collection[i7];
								while (i71<isize71) {
									//integer icandidate72 = list_scan->ikey;
									integer icandidate72 = hash_StrongTranspose_collection1Eco[i71 * n_a[ilevel - 1] + i7];
									if (icandidate72 != i7) {

										// Сосед еще не обрабатывался.
										// закоментировать. Выяснено что без этого лучше сходится. //if ((this_is_C_node[icandidate72] == false) &&
										// закоментировать. Выяснено что без этого лучше сходится. //(this_is_F_node[icandidate72] == false)) 
										{

											// соседи соседей (квадрат матрицы).
											//if (hash_StrongTranspose_collection1[icandidate72] != nullptr)
											{
												//Taccumulqtor_list* list_scan2 = hash_StrongTranspose_collection1[icandidate72];
												integer i72 = 0;
												integer isize72 = isize_hash_StrongTranspose_collection[icandidate72];
												while (i72<isize72) {
													//integer icandidate723 = list_scan2->ikey;
													integer icandidate723 = hash_StrongTranspose_collection1Eco[i72 * n_a[ilevel - 1] + icandidate72];
													if ((icandidate723 != i7)&&(icandidate723 != icandidate72)) {
														// Сосед еще не обрабатывался.
														if ((this_is_C_node[icandidate723] == false) &&
															(this_is_F_node[icandidate723] == false)) {
															this_is_F_node[icandidate723] = true;
														}
													}
													//list_scan2 = list_scan2->next;
													i72++;
												}
											}

											integer i8 = icandidate72;
											integer is0_end2 = row_startA[i8 + 1] - 1;
											
											doublerealT max_vnediagonal2 = -1;

											if (0) {
												for (integer is02 = row_startA[i8]; (is02 <= is0_end2); ++is02) {
													if (Amat.j[is02] != i8) {
														if (Amat.abs_aij[is02] > max_vnediagonal2) {
															max_vnediagonal2 = Amat.abs_aij[is02];
														}
													}
												}
											}
											else {
												max_vnediagonal2 = threshold_quick_only_negative[i8];
											}

											max_vnediagonal2 *= static_cast <doublerealT>(0.23);
												
											for (integer is02 = row_startA[i8]; (is02 <= is0_end2); ++is02) {
												if ((Amat.j[is02] != i7) && (Amat.j[is02] != i8)) {
													if ((Amat.aij[is02] < 0.0) && (Amat.abs_aij[is02] > max_vnediagonal2)) 
													{
															
														// Сосед еще не обрабатывался.
														if ((this_is_C_node[Amat.j[is02]] == false) &&
															(this_is_F_node[Amat.j[is02]] == false)) {
															this_is_F_node[Amat.j[is02]] = true;
														}
															
													}
												}
											}


											
											this_is_F_node[icandidate72] = true;
										}
									}
									//list_scan = list_scan->next;
									i71++;
								}
							}
						}
						newCcount++;
						// Один агрегат создан.
					}
				}
				//else {
					//this_is_F_node[i7] = true;
					// приводит к большей операторной сложности и большему времени расчета
				//}
			}

			
		}

		if ((inumber_of_C_nodes == 0) && (inumber_of_nodes_viewed > 0)) {
			bflag_empty = true;
		}

		if (bprint_mesage_diagnostic) {
			std::cout << "view candidates=" << inumber_of_nodes_viewed << " additional aggregates=";
			std::cout << inumber_of_C_nodes << " inumber_isolated_F_nodes=" << inumber_isolated_F_nodes << "\n";
		}
		//system("pause");
	}

	delete[]  dcount_neighbour;

#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif

}


// Разбиение множества узлов первоначальной сетки на С-узлы и 
// F - узлы. С -узлы составят грубую сетку следующего уровня вложенности.
// Значения же функции в F узлах должно быть восстановлено по значению функции
// в ближайших С узлах (см. задачу интерполяции). 
// Начало написания 20.09.2020. PMIS применённый к квадрату матрицы А.
// Распараллелил 08.10.2020
template <typename doublerealT>
void PMIS_CF_decomposition_applied_to_the_square_of_the_matrix_aggressive(Ak2& Amat, bool*& this_is_F_node,
	bool*& this_is_C_node, integer ilevel,
	integer*& count_neighbour, integer*& n_a,
	integer& newCcount,
	doublerealT*& threshold_quick_only_negative,
	integer*& row_startA,
	integer*& hash_StrongTranspose_collection1Eco,
	integer*& isize_hash_StrongTranspose_collection,
	int imax_aggressive_pass)
{

	bool bprint_mesage_diagnostic = true;
	if (my_amg_manager.iprint_log == 0) {
		bprint_mesage_diagnostic = false;
	}

#ifdef _OPENMP
	int inum_core = number_cores();
	int i_my_num_core_parallelesation = omp_get_max_threads();
	omp_set_num_threads(inum_core); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
#endif

	integer ii_end1 = n_a[ilevel - 1];

	doublerealT* dcount_neighbour = new doublerealT[n_a[ilevel - 1] + 1];

	integer inumber_isolated_F_nodes = 0;

	//std::mt19937 gen(time(0));
	std::mt19937 gen{ std::random_device()() };
	std::uniform_real_distribution<> urd(0, 1);

	// urd(gen) не работает в многопотоке!!!
//#pragma omp parallel for reduction(+ : inumber_isolated_F_nodes)
	for (integer i7 = 1; i7 <= ii_end1; ++i7) {
		// число соседей плюс случайное число от нуля до единицы.
		if (count_neighbour[i7] == 0) {
			this_is_F_node[i7] = true;
			inumber_isolated_F_nodes++;
		}
		// drand случайное вещественное число от нуля до единицы.
		//--->doublerealT drand = static_cast<doublerealT>(((double)(1.0 * rand()) / ((double)(RAND_MAX + 1))));
		doublerealT drand = static_cast<doublerealT>(urd(gen));
		dcount_neighbour[i7] = count_neighbour[i7] + drand;
	}

	// Контролируем и уменьшаем количество пассов.
	int iglobal_number_scan = 0;

	bool bflag_empty = false;
	bool bcontinue = true;
	while (bcontinue) {

		iglobal_number_scan++;
		if (iglobal_number_scan > imax_aggressive_pass) break;

		bcontinue = false;
		integer inumber_of_nodes_viewed = 0; // Число просмотренных узлов на новом проходе (сканировании).
		integer inumber_of_C_nodes = 0; // Число добавленных агрегатов.

#pragma omp parallel for reduction(+:inumber_of_nodes_viewed, inumber_of_C_nodes)
		for (integer i7 = 1; i7 <= ii_end1; ++i7) {
			bool cnd = ((this_is_C_node[i7] == false) && (this_is_F_node[i7] == false));
			if (!cnd) continue;


			inumber_of_nodes_viewed++;

			if (bflag_empty) {
				// Непонятные ошмётки, их небольшое количество делаем С узлами.
				inumber_of_C_nodes++;
				newCcount++;
				this_is_C_node[i7] = true;
				//std::cout<<"lambdai="<< count_neighbour[i7]<<std::endl;
				//system("PAUSE");
			}
			else {

				bcontinue = true;

				doublerealT id_diag = dcount_neighbour[i7];

				// Определение S(i)
				doublerealT max_vnediagonal = -1.0;
				//for (integer is0 = row_startA[i7]; (is0 <= row_startA[i7 + 1] - 1); ++is0) {
					//if (Amat.j[is0] != i7) {
						//if ((Amat.aij[is0]<0.0)&&(Amat.abs_aij[is0] > max_vnediagonal)) {
							//max_vnediagonal = Amat.abs_aij[is0];
						//}
					//}
				//}
				max_vnediagonal = threshold_quick_only_negative[i7];
				//  Забыли учесть transpose(S(i)).
				bool bunion_on = true;
				//bool bOk = false;
				integer is0_end = row_startA[i7 + 1] - 1;
				for (integer is0 = row_startA[i7]; (is0 <= is0_end); ++is0) {
					bool cnd1 = (Amat.j[is0] != i7);
					if (!cnd1) continue;
					bool cnd2 = ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > max_vnediagonal));
					if (!cnd2) continue;



					// Сосед еще не обрабатывался.
					if ((this_is_C_node[Amat.j[is0]] == false) && (this_is_F_node[Amat.j[is0]] == false)) {

						// Соседи соседей.
						integer i8 = Amat.j[is0];
						doublerealT max_vnediagonal2 = threshold_quick_only_negative[i8];
						integer is0_end2 = row_startA[i8 + 1] - 1;
						for (integer is02 = row_startA[i8]; (is02 <= is0_end2); ++is02) {
							if ((Amat.j[is02] != i7) && (Amat.j[is02] != i8)) {
								if ((Amat.aij[is02] < 0.0) && (Amat.abs_aij[is02] > max_vnediagonal2)) {
									if (dcount_neighbour[Amat.j[is02]] > id_diag) {
										// Сосед еще не обрабатывался.
										if ((this_is_C_node[Amat.j[is02]] == false) && (this_is_F_node[Amat.j[is02]] == false)) {
											bunion_on = false;
										}
									}
								}
							}
						}


						if (dcount_neighbour[Amat.j[is0]] > id_diag) {

							bunion_on = false;
						}
						//else {
							// Сосед еще не обрабатывался.
							//if ((this_is_C_node[Amat.j[is0]] == false) && (this_is_F_node[Amat.j[is0]] == false)) {
								//	bOk = true;
							//	}
						//}
					}

				}
				// учитываем transpose(S(i)):
				//if (hash_StrongTranspose_collection1 != nullptr)
				{
					//if (hash_StrongTranspose_collection1[i7] != nullptr)
					{
						//Taccumulqtor_list* list_scan = hash_StrongTranspose_collection1[i7];
						integer i71 = 0;
						integer isize71 = isize_hash_StrongTranspose_collection[i7];
						while (i71 < isize71) {

							//integer icandidate73 = list_scan->ikey;
							integer icandidate73 = hash_StrongTranspose_collection1Eco[i71 * n_a[ilevel - 1] + i7];


							if (icandidate73 != i7) {

								// Сосед еще не обрабатывался.
								if ((this_is_C_node[icandidate73] == false) && (this_is_F_node[icandidate73] == false)) {

									// соседи соседей.
									//if (hash_StrongTranspose_collection1[icandidate73] != nullptr)
									{
										//Taccumulqtor_list* list_scan2 = hash_StrongTranspose_collection1[icandidate73];
										integer i72 = 0;
										integer isize72 = isize_hash_StrongTranspose_collection[icandidate73];

										while (i72 < isize72) {
											//integer icandidate734 = list_scan2->ikey;

											integer icandidate734 = hash_StrongTranspose_collection1Eco[i72 * n_a[ilevel - 1] + icandidate73];

											if ((icandidate734 != i7) && (icandidate734 != icandidate73)) {
												if (dcount_neighbour[icandidate734] > id_diag) {
													// Сосед соседей еще не обрабатывался.
													if ((this_is_C_node[icandidate734] == false) && (this_is_F_node[icandidate734] == false)) {
														bunion_on = false;
													}
												}
											}
											//list_scan2 = list_scan2->next;
											i72++;
										}
									}


									if (dcount_neighbour[icandidate73] > id_diag) {

										bunion_on = false;
									}
									//else {
									//	if ((this_is_C_node[icandidate73] == false) && (this_is_F_node[icandidate73] == false)) {
									//		bOk = true;
									//	}
									//}
								}
							}
							//list_scan = list_scan->next;
							i71++;
						}
					}
				}
				// если bunion_on то на диагонали самая сильная связь и мы создаём агрегат.
				// Сюда в агрегат входят связи из S(i).
				/*
				if (bunion_on) {
					inumber_of_C_nodes++;
					this_is_C_node[i7] = true;
					// Здесь соседей брать только из transpose(S(i)).
					for (integer is0 = row_startA[i7]; (is0 <= row_startA[i7 + 1] - 1); ++is0) {
						// Сосед еще не обрабатывался.
						if ((this_is_C_node[Amat.j[is0]] == false) && (this_is_F_node[Amat.j[is0]] == false)) {
							if ((Amat.j[is0] != i7) && (Amat.abs_aij[is0] >= theta(ilevel) * max_vnediagonal)) {
								this_is_F_node[Amat.j[is0]] = true;
							}
						}
					}
					newCcount++;
					// Один агрегат создан.
				}
				*/
				// если bunion_on то на диагонали самая сильная связь и мы создаём агрегат.
				// Сюда в агрегат входят связи из transpose(S(i)).
				//if (true||bOk) 
				{
					if (bunion_on) {
						inumber_of_C_nodes++;
						this_is_C_node[i7] = true;

						// Увеличивает число итераций с 899 до 998.
						// Не уменьшая разреженность оператора.
						// Вообще всех соединенных с i7 делаем F узлами.
						//integer is0_end3 = row_startA[i7 + 1] - 1;
						//for (integer is02 = row_startA[i7]; (is02 <= is0_end3); ++is02) {
							//if (Amat.j[is02] != i7)  {
								//not//if ((Amat.aij[is02] < 0.0) && (Amat.abs_aij[is02] > max_vnediagonal2)) 
								//{

									// Сосед еще не обрабатывался.
									//if ((this_is_C_node[Amat.j[is02]] == false) && (this_is_F_node[Amat.j[is02]] == false)) {
										//this_is_F_node[Amat.j[is02]] = true;
									//}

								//}
							//}
						//}


						// Здесь соседей брать только из transpose(S(i)).
						//if (hash_StrongTranspose_collection1 != nullptr) 
						{
							//if (hash_StrongTranspose_collection1[i7] != nullptr)
							{
								//Taccumulqtor_list* list_scan = hash_StrongTranspose_collection1[i7];
								integer i71 = 0;
								integer isize71 = isize_hash_StrongTranspose_collection[i7];
								while (i71 < isize71) {
									//integer icandidate72 = list_scan->ikey;
									integer icandidate72 = hash_StrongTranspose_collection1Eco[i71 * n_a[ilevel - 1] + i7];
									if (icandidate72 != i7) {

										// Сосед еще не обрабатывался.
										// закоментировать. Выяснено что без этого лучше сходится. //if ((this_is_C_node[icandidate72] == false) &&
										// закоментировать. Выяснено что без этого лучше сходится. //(this_is_F_node[icandidate72] == false)) 
										{

											// соседи соседей (квадрат матрицы).
											//if (hash_StrongTranspose_collection1[icandidate72] != nullptr)
											{
												//Taccumulqtor_list* list_scan2 = hash_StrongTranspose_collection1[icandidate72];
												integer i72 = 0;
												integer isize72 = isize_hash_StrongTranspose_collection[icandidate72];
												while (i72 < isize72) {
													//integer icandidate723 = list_scan2->ikey;
													integer icandidate723 = hash_StrongTranspose_collection1Eco[i72 * n_a[ilevel - 1] + icandidate72];
													if ((icandidate723 != i7) && (icandidate723 != icandidate72)) {
														// Сосед еще не обрабатывался.
														if ((this_is_C_node[icandidate723] == false) &&
															(this_is_F_node[icandidate723] == false)) {
															this_is_F_node[icandidate723] = true;
														}
													}
													//list_scan2 = list_scan2->next;
													i72++;
												}
											}

											integer i8 = icandidate72;
											integer is0_end2 = row_startA[i8 + 1] - 1;

											doublerealT max_vnediagonal2 = -1;

											if (0) {
												for (integer is02 = row_startA[i8]; (is02 <= is0_end2); ++is02) {
													if (Amat.j[is02] != i8) {
														if (Amat.abs_aij[is02] > max_vnediagonal2) {
															max_vnediagonal2 = Amat.abs_aij[is02];
														}
													}
												}
											}
											else {
												max_vnediagonal2 = threshold_quick_only_negative[i8];
											}

											max_vnediagonal2 *= static_cast <doublerealT>(0.23);

											for (integer is02 = row_startA[i8]; (is02 <= is0_end2); ++is02) {
												if ((Amat.j[is02] != i7) && (Amat.j[is02] != i8)) {
													if ((Amat.aij[is02] < 0.0) && (Amat.abs_aij[is02] > max_vnediagonal2))
													{

														// Сосед еще не обрабатывался.
														if ((this_is_C_node[Amat.j[is02]] == false) &&
															(this_is_F_node[Amat.j[is02]] == false)) {
															this_is_F_node[Amat.j[is02]] = true;
														}

													}
												}
											}



											this_is_F_node[icandidate72] = true;
										}
									}
									//list_scan = list_scan->next;
									i71++;
								}
							}
						}
						newCcount++;
						// Один агрегат создан.
					}
				}
				//else {
					//this_is_F_node[i7] = true;
					// приводит к большей операторной сложности и большему времени расчета
				//}
			}


		}

		if ((inumber_of_C_nodes == 0) && (inumber_of_nodes_viewed > 0)) {
			bflag_empty = true;
		}

		if (bprint_mesage_diagnostic) {
			std::cout << "view candidates=" << inumber_of_nodes_viewed << " additional aggregates=";
			std::cout << inumber_of_C_nodes << " inumber_isolated_F_nodes=" << inumber_isolated_F_nodes << "\n";
		}
		//system("pause");
	}


#pragma omp parallel for reduction(+ : inumber_isolated_F_nodes)
	for (integer i7 = 1; i7 <= ii_end1; ++i7) {
		if ((this_is_C_node[i7] == false) && (this_is_F_node[i7] == false)) {
			this_is_F_node[i7] = true;
		}
	}


	delete[]  dcount_neighbour;

#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif

}

// Разбиение множества узлов первоначальной сетки на С-узлы и 
// F - узлы. С -узлы составят грубую сетку следующего уровня вложенности.
// Значения же функции в F узлах должно быть восстановлено по значению функции
// в ближайших С узлах (см. задачу интерполяции). 
template <typename doublerealT>
void PMIS_CF_decomposition(Ak2& Amat, bool*& this_is_F_node,
	bool*& this_is_C_node, integer ilevel, 
	integer*& count_neighbour, integer*& n_a,
	integer &newCcount,
	doublerealT*& threshold_quick_only_negative,
	integer*& row_startA,
	integer* &hash_StrongTranspose_collection1Eco,
	integer* &isize_hash_StrongTranspose_collection)
{


	bool bprint_mesage_diagnostic = true;
	if (my_amg_manager.iprint_log == 0) {
		bprint_mesage_diagnostic = false;
	}

#ifdef _OPENMP
	int inum_core = number_cores();
	int i_my_num_core_parallelesation = omp_get_max_threads();
	omp_set_num_threads(inum_core); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
#endif

	integer ii_end1 = n_a[ilevel - 1];

	doublerealT* dcount_neighbour = new doublerealT[n_a[ilevel - 1] + 1];


	integer inumber_isolated_F_nodes = 0;

	//std::cout << "number thread=" <<  omp_get_max_threads() << std::endl;
	//system("PAUSE");

	//std::mt19937 gen(time(0));
	std::mt19937 gen{ std::random_device()() };
	std::uniform_real_distribution<> urd(0, 1);

	// urd(gen) не работает в многопотоке!!!
//#pragma omp parallel for reduction(+ : inumber_isolated_F_nodes)
	for (integer i7 = 1; i7 <= ii_end1; ++i7) {
		// число соседей плюс случайное число от нуля до единицы.
		if (count_neighbour[i7] == 0) {
			this_is_F_node[i7] = true;
			inumber_isolated_F_nodes++;
		}
		// drand случайное вещественное число от нуля до единицы.
		//-->doublerealT drand = static_cast<doublerealT>(((double)(1.0 * rand()) / ((double)(RAND_MAX + 1))));
		doublerealT drand = static_cast<doublerealT>(urd(gen));
		dcount_neighbour[i7] = count_neighbour[i7] + drand;
	}


	bool bflag_empty = false;
	bool bcontinue = true;
	while (bcontinue) {
		bcontinue = false;
		integer inumber_of_nodes_viewed = 0; // Число просмотренных узлов на новом проходе (сканировании).
		integer inumber_of_C_nodes = 0; // Число добавленных агрегатов.

#pragma omp parallel for reduction(+:inumber_of_nodes_viewed, inumber_of_C_nodes)
		for (integer i7 = 1; i7 <= ii_end1; ++i7) {
			if ((this_is_C_node[i7] == false) && (this_is_F_node[i7] == false)) {
				inumber_of_nodes_viewed++;

				if (bflag_empty) {
					// Непонятные ошмётки, их небольшое количество делаем С узлами.
					inumber_of_C_nodes++;
					/*newCcount++;*/
					this_is_C_node[i7] = true;
					//std::cout<<"lambdai="<< count_neighbour[i7]<<std::endl;
					//system("PAUSE");
				}
				else {

					bcontinue = true;

					doublerealT id_diag = dcount_neighbour[i7];

					// Определение S(i)
					doublerealT max_vnediagonal = -1.0;
					//for (integer is0 = row_startA[i7]; (is0 <= row_startA[i7 + 1] - 1); ++is0) {
						//if (Amat.j[is0] != i7) {
							//if ((Amat.aij[is0]<0.0)&&(Amat.abs_aij[is0] > max_vnediagonal)) {
								//max_vnediagonal = Amat.abs_aij[is0];
							//}
						//}
					//}
					max_vnediagonal = threshold_quick_only_negative[i7];
					//  Забыли учесть transpose(S(i)).
					bool bunion_on = true;
					//bool bOk = false;
					integer is0_end = row_startA[i7 + 1] - 1;
					for (integer is0 = row_startA[i7]; (is0 <= is0_end); ++is0) {
						if (Amat.j[is0] != i7) {
							if ((Amat.aij[is0] < 0.0) && (Amat.abs_aij[is0] > max_vnediagonal)) {
								if (dcount_neighbour[Amat.j[is0]] >= id_diag) {
									// Сосед еще не обрабатывался.
									if ((this_is_C_node[Amat.j[is0]] == false) && (this_is_F_node[Amat.j[is0]] == false)) {
										bunion_on = false;
									}
								}
								//else {
									// Сосед еще не обрабатывался.
									//if ((this_is_C_node[Amat.j[is0]] == false) && (this_is_F_node[Amat.j[is0]] == false)) {
									//	bOk = true;
								//	}
								//}
							}
						}
					}
					// учитываем transpose(S(i)):
					//if (hash_StrongTranspose_collection1 != nullptr)
					{
						//if (hash_StrongTranspose_collection1[i7] != nullptr)
						{
							//Taccumulqtor_list* list_scan = hash_StrongTranspose_collection1[i7];
							integer i71 = 0;
							integer isize71 = isize_hash_StrongTranspose_collection[i7];
							while (i71<isize71) {
								//integer icandidate73 = list_scan->ikey;
								integer icandidate73 = hash_StrongTranspose_collection1Eco[i71*n_a[ilevel-1]+i7];

								if (icandidate73 != i7) {
									if (dcount_neighbour[icandidate73] >= id_diag) {
										// Сосед еще не обрабатывался.
										if ((this_is_C_node[icandidate73] == false) && (this_is_F_node[icandidate73] == false)) {
											bunion_on = false;
										}
									}
									//else {
									//	if ((this_is_C_node[icandidate73] == false) && (this_is_F_node[icandidate73] == false)) {
									//		bOk = true;
									//	}
									//}
								}
								//list_scan = list_scan->next;
								i71++;
							}
						}
					}
					// если bunion_on то на диагонали самая сильная связь и мы создаём агрегат.
					// Сюда в агрегат входят связи из S(i).
					/*
					if (bunion_on) {
						inumber_of_C_nodes++;
						this_is_C_node[i7] = true;
						// Здесь соседей брать только из transpose(S(i)).
						for (integer is0 = row_startA[i7]; (is0 <= row_startA[i7 + 1] - 1); ++is0) {
							// Сосед еще не обрабатывался.
							if ((this_is_C_node[Amat.j[is0]] == false) && (this_is_F_node[Amat.j[is0]] == false)) {
								if ((Amat.j[is0] != i7) && (Amat.abs_aij[is0] >= theta(ilevel) * max_vnediagonal)) {
									this_is_F_node[Amat.j[is0]] = true;
								}
							}
						}
						//--->newCcount++;
						// Один агрегат создан.
					}
					*/
					// если bunion_on то на диагонали самая сильная связь и мы создаём агрегат.
					// Сюда в агрегат входят связи из transpose(S(i)).
					//if (true||bOk) 
					{
						if (bunion_on) {
							inumber_of_C_nodes++;
							this_is_C_node[i7] = true;
							// Здесь соседей брать только из transpose(S(i)).
							//if (hash_StrongTranspose_collection1 != nullptr)
							{
								//if (hash_StrongTranspose_collection1[i7] != nullptr) 
								{
									//Taccumulqtor_list* list_scan = hash_StrongTranspose_collection1[i7];
									integer i71 = 0;
									integer isize71 = isize_hash_StrongTranspose_collection[i7];
									while (i71<isize71) {
										//integer icandidate72 = list_scan->ikey;
										integer icandidate72 = hash_StrongTranspose_collection1Eco[i71 * n_a[ilevel - 1] + i7];

										if (icandidate72 != i7) {
											// Сосед еще не обрабатывался.
											if ((this_is_C_node[icandidate72] == false) &&
												(this_is_F_node[icandidate72] == false)) {
												this_is_F_node[icandidate72] = true;
											}
										}
										//list_scan = list_scan->next;
										i71++;
									}
								}
							}
							/*newCcount++;*/
							// Один агрегат создан.
						}
					}
					//else {
						//this_is_F_node[i7] = true;
						// приводит к большей операторной сложности и большему времени расчета
					//}
				}

			}
		}

		if ((inumber_of_C_nodes == 0) && (inumber_of_nodes_viewed > 0)) {
			bflag_empty = true;
		}

		if (bprint_mesage_diagnostic) {
			std::cout << "view candidates=" << inumber_of_nodes_viewed << " additional aggregates=" << inumber_of_C_nodes << " inumber_isolated_F_nodes=" << inumber_isolated_F_nodes << "\n";
		}
		//system("pause");
	}

	delete[]  dcount_neighbour;

#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif

}

// Обёртка для сишного realloc().
// В данной функции осуществляется проверка на NULL.
template <typename myARRT>
void my_realloc_memory(myARRT* &x, integer n) {
	if (x != NULL) {
		myARRT* x_temp = NULL;
		x_temp = (myARRT*)realloc(x, ((n) * sizeof(myARRT)));
		if (x_temp == NULL) {
			std::cout<<"application crash for bx array in procedure my_realloc_memory."<<std::endl;
			system("pause");
			exit(1);
		}
		else {
			x = x_temp;
		}
	}
}

// Визуализирует портрет матрицы с помощью графической библиотеки OpenGL.
int PortraitPrint(Ak2 &Amat, integer iadd, integer nnz, integer n) {

#ifndef NO_OPENGL_GLFW

	SCREEN_WIDTH = 500;
	SCREEN_HEIGHT = 500;

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "AliceFlow_v0_49", NULL, NULL);


	glfwSetKeyCallback(window, keyCallback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);


	glViewport(0.0f, 0.0f, screenWidth, screenHeight); // specifies the part of the window to which OpenGL will draw (in pixels), convert from
	// normalised to pixels
	glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame.
	// Here you typically set the zoom factor, aspect ratio and the near and far clipping planes
	glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrtho and
	// glRotate cumulate, besically puts us at (0, 0, 0)
	glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1); // essentially set coordinate system
	glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how objects are transformed (meaning translation, rotation
	// and scaling) in your world

	glLoadIdentity(); // same as above comment 

	GLfloat halfScreenWidth = SCREEN_WIDTH / 2;
	GLfloat halfScreenHeight = SCREEN_HEIGHT / 2;



	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{

		//glClearColor(0.7f, 1.0f, 0.7f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		/* Render OpenGL here */

		glPushMatrix();
		//glColor3d(1, 1, 1);
		glTranslatef(halfScreenWidth, halfScreenHeight, -500);
		glRotatef(rotationX, 1, 0, 0);
		glRotatef(rotationY, 0, 1, 0);
		glTranslatef(-halfScreenWidth, -halfScreenHeight, 500);

		//glPointSize(10);
		/*сообщаем что нужно рисовать точку или точки*/
		glBegin(GL_POINTS);
		/*передам коордианты что нарисовать точку в координатах X=0 и Y=0*/
		glColor3f(1.0, 1.0, 1.0);

		for (integer i58 = 1 + iadd; i58 <= nnz + iadd; ++i58) {
			//fprintf(fp_portrait, "%d %d\n", Amat.i[i58], Amat.j[i58]);

			//glVertex3f(0.5, 0.5, 0);
			glVertex2f(500.0 * Amat.i[i58] / (1.0 * n), 500.0 - 500.0 * Amat.j[i58] / (1.0 * n));
		}
		/*сообщаем что завершили рисовать точку или точки*/
		glEnd();

		glPopMatrix();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	glfwTerminate();

#endif

	return 0;
}



// 29.12.2021 ошибочное распараллеливание в строках 10730, 10826, 11033, 11125.
// Я его отключил-закоментировал. 10.01.2022.
// 08.06.2021 - xx.xx.xxxx Версия 6 на основе версии 4.
// 29.07.2018 - xx.xx.xxxx Версия 6 на основе версии 4.
// 08.06.2021 Всего используются три сортировки. Однократно исходной матрицы и
// двукратно в цикле одна сортировка для P и одна сортировка для R. 
// Оптимизирована скорость работы алгоритмов сортировки. Библиотечная сортировка 
// std::sort на основе quick sort существенно быстрее чем собственная реализация quick sort, 
// также std::sort быстрее чем пирамидальная и сортировка Тима,  поэтому
// имеет смысл переходить на неё. Осуществлён полный переход на std::sort во всех трех вышеописанных 
// вызовах. Сохранены пирамидальная и сортировка Тима. Сортировки std::sort, пирамидальная и Тима 
// работают в два потока с помощью omp parallel sections это ускоряет каждую из них в 1.6 раза.
// Достигнута экономия времени решения задачи не менее 7% на 2.7млн неизвестных в уравнении Пуассона.
// 24.04.2020 Выделил алгоритм селектора(C/F разбиения) в отдельную функцию.
// 14.04.2020 Разделил на две функции setup и solution phase.
// 13.04.2020 Убрал сортировку CountingSort из параллельного умножения
// разреженных матриц. Использование #pragma omp for schedule (static)
// группирует работу нитей последовательно. Данные разбиты на 8 частей (потоков),
// нулевая нить обрабатывает первую часть, первая вторую и т.д. Потом
// после последовательного объединения работы нитей итоговая матрица А
// окажется отсортированной по строкам. Это дало ускорение 12%.
// 12.04.2020 Умножение разреженных матриц вынесено в отдельные функции.
// 3.02.2019 Начало внедрения более гибкого типа данных Ak2, 
// что позволит сэкономить оперативную память.
// 25.04.2018 Версия №4 classic_aglomerative_amg4 это основная поддерживаемая версия.
// Пятая версия classic_aglomerative_amg5 давно не поддерживается (заморожена).
// июнь 2017 - добавлен Рунге-Кутта smoother, улучшена поддержка ilu0 разложения в алгоритме. 
// июнь 2017 - Поддерживается максимальное количество уровней вложенности 100 и менее. 
// зимние каникулы 2016-2017 года - добавлен bicgStab.
// Лето 2017 - алгебраический многосеточный метод теперь всё чаще и чаще используется как
// предобуславливатель к алгоритму Хенка ван дер Ворста BiCGStab.
// Эта связка показывает более стабильную и
// надежную работу чем просто отдельно amg.
// 4-6 ноября 2016. Добавлен ILU0. Полностью удалён устаревший код из Solution Phase.
// 9 августа 2016. Зейдель не справляется с большими спектральными радиусами матриц даже 
// в составе данного amg. Это же проявляется и на классическом amg1r5. 
// 9 августа решено уменьшить спектральный радиус в Зейделе 
// на каждом уровне вложенности с помощью ILU2 декомпозиции. Это подтверждает статья 
// Е.М.Андреева, Г.В.Муратова
// "Многосеточный метод решения сильно несимметричных систем" ЮГИНФО РГУ, Ростов-на-Дону,
// Россия. Там они
// показывают расходимость мультигрида на основе Зейделя для задач с существенным спектральным
// радиусом и рекомендуют заменить Зейделя на ТКМ2 метод
// (треугольный кососимметричный метод). В данной 
// программе у нас есть 
// успешный опыт использования ILU2 предобуславливателя из библиотеки SPARSKIT2 Ю.Саада 
// поэтому вместо ТКМ2 у нас 
// будет ILU2.
// 22 января 2016 текущий работоспособный вариант кода.
// Планы: 1. сделать версию amg3. 
// В ней: 2. заменить все проверки на невыделение оперативной памяти на универсальную 
// функцию единую для всего.
// это очевидно немного сократит программный код.
// В ней 3. код V цикла оформить в виде цикла for. Тогда же можно будет попробовать вставить
// direct метод для самого грубого уровня.
// Это же откроет возможности сделать из amg алгоритма предобуславливатель для BiCGStab.
// 15 января 2016 экономим память переходим на Ak1.
// 10 января 2016 двоичный поиск заменён на хеширование.
// 15 декабря 2015. Данная версия кода будет полностью очищена от устаревшего кода.
// 13 декабря 2015. Внедрено АВЛ дерево. При внедрении АВЛ дерева исправлена логическая ошибка в
// построении С/F разбиения, теперь C/F разбиение строится корректно.
// Исправлен и внедрён quicksort (qs,qsj)
// который в пять раз быстрее пирамидальной сортировки.
// Полный отказ от band_size!!!.
// Время работы алгоритма на 1.7М неизвестных в 3D составило ровно 1 минуту.
// 18 октября 2015. Полностью работоспособный мультигрид.
// Тестировалось на условиях Дирихле но должно работать на любых 
// краевых задачах. 18 октября 2015 датируется версия 0.04. Версия 0.04 на треть
// быстрее версии 0_03. Были ускорены как операции построения C/F разбиения, 
// так и нахождение оператора Галёркина. При нахождении С/F разбиения 
// учитывается уже построенная его часть и поэтому число сканирований на
// на поздних циклах сокращается охватывая только не построенную часть.
// При нахождении произведения Галёркина получена самая оптимальная по 
// быстродействию версия,
// Основанная на алгоритме слияния отсортированных списков.
// 4 октября 2015 правильное построение последовательности вложенных графов.
// 30 сентября 2015 продолжаем исправление метода. Делаем классический 
// алгебраический многосеточный метод на основе  C/F разбиения.
// 16 сентября 2015 года обнаружено что операции 
// сгрубления и интерполяции сделаны совершенно неверно,
// и если сгрубление еще в какой-то мере проецирует то интерполяция просто никакая.
// Операции сгрубления и интерполяции будут сделаны заново на основе статьи 
// К.Н. Волкова в новой версии солвера.
// 3 september 2015 Villa Borgese.
// Возвращает divergence detected.
// Использует внутри структуру данных my_amg_manager.
template <typename doublerealT>
void setup_phase_classic_aglomerative_amg6(Ak2& Amat,
	integer nsizeA, // количество ячеек выделенное извне для хранилища матриц А	
	integer nnz, // number of non zero elements
	integer n, // dimension of vectors x and b.	
	doublereal*& x, //solution (решение) 
	doublereal*& b, // rthdsd (правая часть).
	real_mix_precision& ret74,
	integer iVar,
	bool bmemory_savings,
	// Параметры нужные только для solution phase.
	amg_precond_param& amg_pp,
	integer*& n_a,  // количество неизвестных на каждом из построенных уровней.
	integer*& nnz_a, // количество ненулевых коэффициентов в матрице на каждом из построенных уровней.
	integer*& nnz_aRP, // Количество ненулевых коэффициентов в матрице оператора интерполяции на каждом из построенных уровней.
	integer& ibsp_length,
	BAD_STRING_PATCHING*& bsp,
	integer i_bsp_LIMIT,
	bool*& flag,
	bool*& F_false_C_true, // С/F разбиение построенное данной функцией на каждом уровне и перадаваемое в функцию решения solution_phase.
	Ak1*& P, // оператор интерполяции построенный данной функцией на каждом уровне и передаваемый ыв функцию solution_phase для решения.
	bool bQuick_sort_for_reorder
) {

	bool bprint_mesage_diagnostic = true;
	if (my_amg_manager.iprint_log == 0) {
		bprint_mesage_diagnostic = false;
	}

	// дополнительное упорядочивание по столбцам в строке в 
	// надежде что кеш будет лучше использоваться.
	// Его можно безболезненно выключить.
	//const bool bQuick_sort_for_reorder = false;

	// Parallel Modifined Independed Set
	bool bPMIS = false;
	// PMIS применённый к квадрату матрицы А. начало 20.09.2020.
	// Всё равно на этапе интерполяции для сходимости надо добавлять новые С узлы,
	// и тот выигрыш разреженности который достигнут при C/F разбиении 
	// примененного к квадрату матрицы А теряется на этапе интерполяции при
	// добавлении новых С узлов. Если на этапе интерполяции сохранять разреженность
	// и не добавлять новых С узлов то сходимость очень быстро сильно ухудщается и 
	// время расчёта неприемлемо возрастает. 21.09.2020.
	bool bPMIS_applied_to_the_square_of_the_matrix = false;// true;

	//doublereal theta = my_amg_manager.theta;
	//doublereal theta83 = my_amg_manager.theta;
	doublerealT magic82 = static_cast<doublerealT>(my_amg_manager.magic);
	//doublerealT magic83 = static_cast<doublerealT>(my_amg_manager.magic);


	const bool b_REALLOC = false;

	integer* C1 = (integer*)malloc((n + 1) * sizeof(integer));
	char c19[3] = "C1";
	char c29[14] = "Counting_Sort";
	handle_error<integer>(C1, c19, c29, (n + 1));

	integer* C2 = (integer*)malloc((n + 1) * sizeof(integer));
	char c39[3] = "C2";
	char c49[14] = "Counting_Sort";
	handle_error<integer>(C2, c39, c49, (n + 1));

	Ak1* Bm1 = nullptr;
	Ak1* Bm2 = nullptr;

	if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
		(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
		(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
		(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
	{
		Bm1 = (Ak1*)malloc((nnz + 2) * sizeof(Ak1));
		char c59[4] = "Bm1";
		char c69[14] = "Counting_Sort";
		handle_error<Ak1>(Bm1, c59, c69, (nnz + 2));

		Bm2 = (Ak1*)malloc((nnz + 2) * sizeof(Ak1));
		char c79[4] = "Bm2";
		char c89[14] = "Counting_Sort";
		handle_error<Ak1>(Bm2, c79, c89, (nnz + 2));

	}



	// Не более 100 ошибочных строк
	// с отрицательной диагональю.
	//const int i_bsp_LIMIT = 100;
	if (bsp != nullptr) {
		delete[] bsp;
		bsp = nullptr;
	}
	bsp = new BAD_STRING_PATCHING[i_bsp_LIMIT];
	ibsp_length = 0;

	//integer &nsizePR, // Память под P в количествах n.
	//Ak1* &R, // restriction
	//Ak1* &P, // prolongation  

	//Ak1* R = nullptr; 
	// Оператор проекции равен оператору интерполяции с точностью до транспонирования.
	//Ak1* P = nullptr;

	// Было значение 35 до 04.09.2020.
	integer nsizePR = 12; // 04,09,2020 - Память под P в количествах n.
	//nsizePR = 135;
	if (iVar == TEMP) {
		if (bonly_solid_calculation) {
			// 31 октября 2016.
			// По моим замерам с надёжным запасом в 30% для всех твёрдотельных 
			// задач должно хватить значения 12*n.
			//nsizePR = 135; // было 12. 01,02,2020 стало 35*n.

			bool bmyconvective = false;
			if (starting_speed_Vx * starting_speed_Vx + starting_speed_Vy * starting_speed_Vy + starting_speed_Vz * starting_speed_Vz > 1.0e-30) {
				if (f[0].maxelm > 0) {
					bmyconvective = true;
				}
			}
			else {
				// Загрузка распределения начальной скорости.

				FILE* fp_inicialization_data = NULL;
#ifdef MINGW_COMPILLER
				int err_inicialization_data = 0;
				fp_inicialization_data = fopen64("load.txt", "r");
				if (fp_inicialization_data == NULL) err_inicialization_data = 1;
#else
				errno_t err_inicialization_data = 0;
				err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif

				if (err_inicialization_data == 0) {
					// Открытие удачно и файл присутствует.
					if (f[0].maxelm > 0) {
						bmyconvective = true;
					}
					fclose(fp_inicialization_data);
				}
			}

			if (!bmyconvective) {
				if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
				{
					// PMIS, HMIS
					if (my_amg_manager.number_interpolation_procedure == 0) {
						nsizePR = 5;// 3;
					}
					else
					if (my_amg_manager.number_interpolation_procedure == 2) {
						nsizePR = 8;// 8;
					}
					else {
						nsizePR = 5; // 12.09.2020
					}
				}
				else {
					nsizePR = 9; // 04.09.2020
				}
			}
			else {
				if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
				{
					if (my_amg_manager.number_interpolation_procedure == 2) {
						nsizePR = 8;
					}
					else {
						nsizePR = 7;// 30.11.2021 5; // 12.09.2020
					}
				}
				else {
					nsizePR = 12; // 04.09.2020 Даже 9.1*n уже с 30% запасом.
				}
			}
		}
	}
	if (iVar == TOTALDEFORMATIONVAR) {
		// Механическая задача.
		nsizePR = 22; // 04.09.2020
	}



	//R = new Ak1[static_cast<integer>(35 * n) + 1]; // 3*nnz 2.4 // 35
	//--->R = (Ak1*)malloc(((nsizePR * n) + 1) * sizeof(Ak1));
	//char c7[2] = "R";
	//handle_error<Ak1>(R, c7, c1, ((nsizePR * n) + 1));

	//P = new Ak1[static_cast<integer>(35 * n) + 1]; // 3*nnz 2.4 // 35
	P = (Ak1*)malloc(((nsizePR * n) + 1) * sizeof(Ak1));
	char c1[22] = "my_agr_amg_loc_memory";
	char c8[2] = "P";
	handle_error<Ak1>(P, c8, c1, ((nsizePR * n) + 1));
	if (bprint_mesage_diagnostic) {
		std::cout << "Prolongation alloc succseful..." << std::endl;
	}

	//if (R != nullptr) {
		// Используется только оператор P, оператор R точно такой же что и P с точностью до сортировки.
		// Методы restriction и prolongation не чувствительны к сортировке.
		//free(R);
		//R = nullptr;
	//}

	//  используем хеш-таблицу.
	construct_hash_table_Gus_struct01(n);


	// 23.12.2016 ускорение счёта нелинейных задач:
	// лучистые потоки обновляются после каждого V цикла,
	// для этого внутрь передаётся 
	// b и lb.


	integer* row_ind_PE = nullptr;
	integer* row_ind_PS = nullptr;
	integer* row_ind_AS = nullptr;
	integer* row_ind_AE = nullptr;
	integer* index_visit = nullptr;
	doublerealT* vector_sum = nullptr;
	bool* hash_table = nullptr;

	integer* row_ind_ER = nullptr;
	integer* row_ind_SR = nullptr;
	integer iend_marker_position;
	doublerealT* ap_coarse = nullptr;
	integer icounter = 1;
	integer icount1;
	integer numberofcoarcenodes;
	integer* C_numerate = nullptr;
	bool bweSholdbeContinue = true;
	integer the_number_of_neighbors_that_are_not_C_nodes = 0;
	integer number_of_F_nodes_with_one_single_strong_C_neighbor = 0;
	integer number_of_F_nodes_with_one_single_strong_C_neighborF = 0;

	integer iadditionalCstatistic = 0;

	integer newCcount = 0;

	integer* row_startA = nullptr;
	integer* count_neighbour = nullptr;

	bool identiti = true;


	// Вершина технологии решения плохообусловленных разреженных СЛАУ: BiCGStab + camg(РУМБА).
	// 1. многосеточные технологии.
	// 2. предобуславливание.
	// 3. стабилизация.
	// Если my_amg_manager.istabilization == 1 то мы используем метод бисопряженных градиентов со стабилизацией с предобуславливанием 
	// классическим алгебраическим многосеточным методом РУМБА.
	// Начало реализации 5.01.2017.(more robust).
	// Если my_amg_manager.istabilization == 0 - То просто используется 
	// многосеточный решатель без какого либо метода Крыловского подпространства.
	// Если my_amg_manager.istabilization == 2 - То используется fgmres - 
	// алгоритм Юсефа Саада и Мартина Г. Шульца (гибкий вариант обобщённого метода минимальных невязок) в котором 
	// на каждой итерации алгоритма fgmres делается одно многосеточное предобуславливание (один V цикл). 
	//bool bBiCGStab_plus_RUMBA_camg = true;
	//if (my_amg_manager.istabilization == 0) {
	// Просто многосеточный метод без какого-либо Крыловского подпространства.
	// none
	//bBiCGStab_plus_RUMBA_camg = false;
	//}


	bfirst_jacoby_start = true;

	bool from_re_operation_protection0 = true;
	integer ifrom_re_operation_protection = 0;

	// Универсальные сглаживающие процедуры. 4 ноября 2016.
	// ILU2 smoother
	// 0 - ILU не используется. используется Gauss-Seidel.
	// 1 - ILUk(k==lfil) используется.
	// 2 - 
	integer bILU2smoother = 0;
	if (my_amg_manager.ilu2_smoother == 1) {
		// Включаем ILUk(k==lfil) сглаживатель. 
		// он ест больше памяти но более быстро сходится.
		// Есть надежда что он справится с гораздо более плохообусловленными матрицами.

		//---->bILU2smoother = 1; // ILU0

						   // По - видимому алгоритм 
						   // ilu0_(maxelm_plus_maxbound, milu0.val, milu0.col_ind, milu0.row_ptr, milu0.alu, milu0.jlu, milu0.ju, milu0.iw, ierr);
						   // является дефектным. Я не получил с ним сходимости как ни пытался. Зато алгоритм iluk с lfil=0 проявил себя наилучшим 
						   // образом и я его рекомендую к использованию. Это реализовано в ветке кода my_amg_manager.ilu2_smoother == 2.
						   // Причём iluk с lfil=0 работает на всех уровнях и прекрасно себя проявляет.

		// Включаем ILU2 сглаживатель. 
		// он ест больше памяти но более быстро сходится.

		// Его рекомендуется применять только для исходной матрицы - уровень ноль.
		// Если его применять на более глубоких уровнях то сходимость лишь замедляется.

		 // ILU2 ест слишком много оперативной памяти и я его заменил на ILU0 сглаживатель на каждом уровне: iluk с lfil=0.
		// Возможно я ещё вернусь к ilu2 хотя-бы на нулевом уровне, т.к. там он особенно хорош.

						   // Перенаправление.
		bILU2smoother = 2; // ILU0 // ILU2
	}

	if (my_amg_manager.ilu2_smoother == 2) {
		// Рунге Кутта сглаживатель третьего порядка.
		my_amg_manager.iRunge_Kutta_smoother = 3;
	}
	if (my_amg_manager.ilu2_smoother == 3) {
		// Рунге Кутта сглаживатель пятого порядка.
		my_amg_manager.iRunge_Kutta_smoother = 5;
	}

	if (my_amg_manager.ilu2_smoother == 6) {
		// FGMres smoother
		// gmres smmother. K-smoother.
		my_amg_manager.b_gmres = true;
		my_amg_manager.b_ChebyshevSmoother = false;
	}
	if (my_amg_manager.ilu2_smoother == 7) {
		// spai0 smoother.
		my_amg_manager.b_spai0 = true;
	}
	if (my_amg_manager.ilu2_smoother == 8) {
		// BiCGStab smoother
		// gmres smmother. K-smoother.
		my_amg_manager.b_gmres = true;
		my_amg_manager.b_ChebyshevSmoother = true;
	}

	//std::cout << "smoother type = " << my_amg_manager.ilu2_smoother << std::endl;
	//getchar();

	//bILU2smoother = 0; // only seidel sor smoother.
	const doublerealT dapply_ilu_max_pattern_size = static_cast <doublerealT>(9.2);

	// Параметры отвечающие за автоматическую настройку SOR.
	// По трём точкам мы построим параболу и на её основе 
	// спрогнозируем улучшенный параметр релаксации omega_optimal.
	// Парабола представляется намного лучшей чем простая линейная экстраполяция.
	bproblem_amg_convergence1 = false;
	bproblem_amg_convergence2 = false;
	bproblem_amg_convergence3 = false;
	gold_const = 0.2;




	bool bpositive_connections_CF_decomp = true;
	MY_AMG_SPLITTING_COARSENING_ALGORITHM  memo_icoarseningtype = my_amg_manager.icoarseningtype;

	if (my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_NEG_CONNECTION) {
		// only negative connections 
		// Внедиагональные положительные связи игнорируются при создании C/F разбиения.
		bpositive_connections_CF_decomp = false;
		my_amg_manager.icoarseningtype = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION;
	}
	if (my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_NEG_CONNECTION) {
		// only negative connections 
		// Внедиагональные положительные связи игнорируются при создании C/F разбиения.
		bpositive_connections_CF_decomp = false;
		my_amg_manager.icoarseningtype = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION;
	}
	if (my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_NEG_CONNECTION) {
		// only negative connections 
		// Внедиагональные положительные связи игнорируются при создании C/F разбиения.
		bpositive_connections_CF_decomp = false;
		my_amg_manager.icoarseningtype = MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ST_ALL_CONNECTION;
	}
	if (my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_NEG_CONNECTION) {
		// only negative connections 
		// Внедиагональные положительные связи игнорируются при создании C/F разбиения.
		bpositive_connections_CF_decomp = false;
		my_amg_manager.icoarseningtype = MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION;
	}


	if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
		(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
		(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS)) {
		// Parallel Modified Independed Set.
		bPMIS = true;
	}
	// 19.01.2016 Для построения C/F разбиения и интерполяции используется разная логика
	// в области игнорирования и не игнорирования positive connections.
	// Требует обсуждения следующий вопрос: 
	// 1. При построении процедуры интерполяции важны все связи как позитив так и негатив.
	// 2. При построении C/F декомпозиции важны только негатив связи. 
	// Это гипотеза требующая подтверждения.
	// Разделение между bpositive_connections_CF_decomp используемом при построении
	// C/F декомпозиции и bpositive_connections
	// Произошло 19.01.2017.


	bool bpositive_connections = true;
	// 23 октября 2016
	if (bSIMPLErun_now_for_temperature) {
		// Решаем cfd задачи.
		//bprint_mesage_diagnostic = false;

		// Гипотеза в том, что positive connections 
		// ускоряющие задачи теплопередачи в твёрдом теле приводят 
		// к расходимости в гидродинамических задачах:
		// гипотеза неверна, с убранными positive connections сходимость только хуже.
		//bpositive_connections = false;
	}


	// Для вычисления grid complexity оператора интерполяции:
	integer nnz_P_memo_0 = 0;
	integer nnz_P_memo_all = 0;

	// Надо заменить все new на malloc.
	//theta = 0.25;

	bool bonly_serial = true;

	// контроль числа сильных связей между переменными.
	// doublerealT theta = 0.25;  // 0.25 for 2D and 0.5 for 3D 


	//const integer MY_SORT_ALGORITHM::QUICK_SORT = 1; // Быстрая сортировка Ч. Хоара.
	// Использовать ли quicksort qs and qsj.
	// Сортировка с подсчётом быстрее quickSort.
	// Использовать ли сортировку подсчётом которая 
	//потребляет килотонну памяти (Короче для машин у которых море оперативки).
	//const integer MY_SORT_ALGORITHM:: COUNTING_SORT = 0; // Сортировка с подсчётом лучший выбор.
	// Сортировка с подсчётом подходит потому что ключи целочисленны и 
	// лежат в заданном интервале непрерывно.
	//const integer MY_SORT_ALGORITHM::HEAP_SORT = 2; // пирамидальная сортировка.
	// количество рекурсивных вызовов ограничено, поэтому QuickSort не подходит.
	// В компиляторе надо увеличить размер стека до 4Мб.
	//bmemory_savings =false при MY_SORT_ALGORITHM::QUICK_SORT и MY_SORT_ALGORITHM::HEAP_SORT;
	MY_SORT_ALGORITHM imy_sort_algorithm = my_amg_manager.imySortAlgorithm;// MY_SORT_ALGORITHM:: COUNTING_SORT;

	const doublerealT RealZERO = static_cast <doublerealT>(1.0e-30);// 1.0e-10;
	//const doublerealT divisionZERO = static_cast <doublerealT>(1.0e-30);
	const doublerealT RealMAXIMUM = static_cast <doublerealT>(1.0e30);
	// включение/отключение отладочного режима.
	bool debug_reshime = false;


	const integer maxlevel = 201; // (101 до 16.02.2019) (51 до 5.06.2017) 30
	integer ilevel = 1;
	if (n_a != nullptr) {
		delete[] n_a;
		n_a = nullptr;
	}
	if (nnz_a != nullptr) {
		delete[] nnz_a;
		nnz_a = nullptr;
	}
	n_a = new integer[maxlevel];
	nnz_a = new integer[maxlevel];
	nnz_a[0] = nnz;
	n_a[0] = n;

	// number_cores() 01.08.2021
	const int iKnumber_thread = number_cores(); // 8;

	//char *str_Function_name = "classic_aglomerative_amg_6";


	Ak1** AccumulqtorA_m = nullptr;
	doublerealT** vector_sum_m = nullptr;
	integer** index_visit_m = nullptr;
	bool** hash_table_m = nullptr;
	integer* index_size_m = nullptr;
	integer* istartAnew_m = nullptr;

	// Были обнаружены задачи теплопередачи с очень сложной пространственной формой для 
		// которых значение 5.55*А мало. Необходимое значение с запасом 11.0*А.
		// Реальное значение не превышало 7.82*А.
		//  Эти задачи - кольцевой радиатор для 1500Вт ВУМ и вторая задача - разномасштабная
		// геометрия в пять порядков: от 17.5см до 0.5мкм. 
		// было значение 5.55; 11; 30.
		//doublereal AccumulqtorA_m_SIZE8 = static_cast<doublereal>(1.3*my_amg_manager.memory_size); // было значение 2.0.
	doublereal AccumulqtorA_m_SIZE8 = 8.0; // было значение 2.0.


	if ((number_cores() == 1)) {

		// Однопоточный код. Никакой дополнительной памяти выделять ненадо. 
		// Это сильно экономит время.

	}
	else {

		//#ifdef _NONAME_STUB29_10_2017
#ifdef _OPENMP 
	// Данные используемые для частичного формирователя суммы.
	// iKnumber_thread=8 - Это число потоков.

		
		AccumulqtorA_m = new Ak1 * [iKnumber_thread];
		
		vector_sum_m = new doublerealT * [iKnumber_thread];
		
		index_visit_m = new integer * [iKnumber_thread];
		hash_table_m = new bool* [iKnumber_thread];
		
		istartAnew_m = new integer[iKnumber_thread];

		index_size_m = new integer[iKnumber_thread];
		
		//if (bPMIS) 
		if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
			(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
			(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
			(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
		{
			// PMIS or HMIS
			if (b_on_adaptive_local_refinement_mesh) {
				if (my_amg_manager.number_interpolation_procedure == 0) {

					AccumulqtorA_m_SIZE8 = 1.5; // Дальнобойная интерполяция.

				}
				if (my_amg_manager.number_interpolation_procedure == 2) {

					AccumulqtorA_m_SIZE8 = 2.2; // Дальнобойная интерполяция.
				}
				else {

					AccumulqtorA_m_SIZE8 = 2.2; // 2.0;
				}
			}
			else {
				if (my_amg_manager.number_interpolation_procedure == 0) {

					AccumulqtorA_m_SIZE8 = 1.5; // Дальнобойная интерполяция.

				}
				if (my_amg_manager.number_interpolation_procedure == 2) {

					AccumulqtorA_m_SIZE8 = 2.2; // Дальнобойная интерполяция.

				}
				else {
					AccumulqtorA_m_SIZE8 = 2.2;// 2.2 static;
				}
			}
		}
		//if (n > 20000000) {// Более 20М
			// Экономим память
			//AccumulqtorA_m_SIZE8 = 11.0;
		//}
		for (int i_9 = 0; i_9 < iKnumber_thread; ++i_9) {

			AccumulqtorA_m[i_9] = new Ak1[static_cast<integer>((1.0 / (1.0 * iKnumber_thread)) * AccumulqtorA_m_SIZE8 * nnz + 1)];
			if (bprint_mesage_diagnostic) {
				std::cout << "AccumulqtorA_m[" << i_9 << "] alloc succseful..." << std::endl;
			}
			vector_sum_m[i_9] = new doublerealT[n + 1];
			//vector_sum_m[i_9] = (doublerealT*)malloc((n + 1) * sizeof(doublerealT));
			//handle_error<doublerealT>(vector_sum_m[i_9], "vector_sum_m[i_9]", "classic_aglomerative_amg_6", (n + 1));

			index_visit_m[i_9] = new integer[n + 1];
			//index_visit_m[i_9] = (integer*)malloc((n + 1) * sizeof(integer));
			//handle_error<integer>(index_visit_m[i_9], "index_visit_m[i_9]", "classic_aglomerative_amg_6", (n + 1));

			hash_table_m[i_9] = new bool[(n + 2)];
			//hash_table_m[i_9] = (bool*)malloc((10 * n + 1) * sizeof(bool));
			//handle_error<bool>(hash_table_m[i_9], "hash_table_m[i_9]", "classic_aglomerative_amg_6", (10 * n + 1));

			/*integer i_91_end = n + 2;
	#pragma omp parallel for
			for (integer i_91 = 0; i_91 < i_91_end; ++i_91) {
				hash_table_m[i_9][i_91] = false;// inicialization
			}*/
			//#pragma omp parallel for
				//	for (integer i_91 = 0; i_91 < n + 1; ++i_91) {
						//vector_sum_m[i_9][i_91] = 0.0;
					//}
			index_size_m[i_9] = 0;
			istartAnew_m[i_9] = 0;
	}
#endif

}
	// quick - вычисляется один раз, а потом везде только используется.
	// threshold - пороговое значение контролирующее число сильных связей между переменными.
	// Определяется по модулям внедиагональных коэффициентов.
	// размер от 0 до n включительно.
	doublerealT* threshold_quick_all = my_declaration_array<doublerealT>(n, -1.0, "threshold_quick_all");

	// threshold - пороговое значение контролирующее число сильных связей между переменными.
	// Определяется только по отрицательным внедиагональным коэффициентам.
	// размер от 0 до n включительно.
	doublerealT* threshold_quick_only_negative = my_declaration_array<doublerealT>(n, -1.0, "threshold_quick_only_negative");
	bool btreshold_on_new_vetv = true; // false откат изменений назад на старую стабильную ветвь кода.


	// flag n+1	
	// размер от 0 до n включительно.
	if (flag != nullptr) {
		free(flag);
		flag = nullptr;
	}
	flag = my_declaration_array<bool>(n, false, "flag");
	

	// hash_table2  n+1
	// Для построения C/F декомпозиции нам тоже потребуется хеш-таблица
	// и стек для очистки хеш-таблицы.
	// Инициализация требуется и выполнена. 
	// размер от 0 до n включительно.
	bool* hash_table2 = my_declaration_array<bool>(n, false, "hash_table2");

	// istack n+1
	// И теперь стек для очистки хеш-таблицы.
	// размер от 0 до n включительно.
	integer* istack = my_declaration_array<integer>(n, -1, "istack");
	


	integer iadd = 0;
	integer nnzR = 1;
	integer iaddR = 0;
	if (nnz_aRP != nullptr) {
		delete[] nnz_aRP;
		nnz_aRP = nullptr;
	}
	nnz_aRP=new integer[maxlevel];
	bool bcontinue_global = true;
	
	bool* this_is_C_node = my_declaration_array<bool>(n, false, "this_is_C_node");
	bool* this_is_F_node = my_declaration_array<bool>(n, false, "this_is_F_node");

	// инициализация требуется обязательно.
	if (F_false_C_true != nullptr) {
		free(F_false_C_true);
		F_false_C_true = nullptr;
	}
	F_false_C_true = my_declaration_array<bool>(4 * n, false, "F_false_C_true");

	

	bool bStrongTransposeON = true; // Как в литературе используем Strong Transpose.
	if (my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::CLASSICAL_ALL_CONNECTION) {
		bStrongTransposeON = false;
	}
	
	node_AVLST** hash_StrongTranspose_collection = nullptr;
	//Taccumulqtor_list** hash_StrongTranspose_collection1 = nullptr;
	//integer isize_memory_alloc_hash_StrongTranspose_collection1 = -1;
	integer *isize_hash_StrongTranspose_collection = nullptr;

	integer isize_hash_StrongTranspose_collection1Eco = 20;
	if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY))
	{
		isize_hash_StrongTranspose_collection1Eco = 160;
	}
	else {
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			isize_hash_StrongTranspose_collection1Eco = 160; //40
		}
		else {
			if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
				(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) || 
				(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
				(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
			{
				if (bmyconvective7248) {
					// Конвекция.
					if (my_amg_manager.number_interpolation_procedure == 2) {
						// Long range interpolation distance 2.
						if (iswitchMeshGenerator != CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER) {
							isize_hash_StrongTranspose_collection1Eco = 520; //20
						}
						else {

							isize_hash_StrongTranspose_collection1Eco = 80; //20
						}
					}
					else {
						if (iswitchMeshGenerator != CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER) {
							isize_hash_StrongTranspose_collection1Eco = 320; //20
						}
						else {
							isize_hash_StrongTranspose_collection1Eco = 60; //20
						}
					}
				}
				else {
					if (b_on_adaptive_local_refinement_mesh) {
						isize_hash_StrongTranspose_collection1Eco = 101; //81
					}
					else {
						if (my_amg_manager.number_interpolation_procedure == 0) {
							isize_hash_StrongTranspose_collection1Eco = 10; //17
						}
						else {
							isize_hash_StrongTranspose_collection1Eco = 60; //17
						}
					}
				}
			}
		}
	}
	// Будем выделять память ровно один раз и хранить все данные  одномерном массиве:
	integer* hash_StrongTranspose_collection1Eco = new integer[isize_hash_StrongTranspose_collection1Eco *n];
	isize_hash_StrongTranspose_collection = new integer[n + 1];

#ifdef _OPENMP

	if ((number_cores() == 1)) {

		// Однопоточный код. Никакой дополнительной памяти выделять ненадо. 
		// Это сильно экономит время.

		vector_sum = my_declaration_array<doublerealT>(n, 0.0, "vector_sum");
		index_visit = my_declaration_array<integer>(n, 0, "index_visit");
		hash_table = my_declaration_array<bool>(n, false, "hash_table");

	}

#else
	vector_sum = my_declaration_array<doublerealT>(n, 0.0, "vector_sum");
	index_visit = my_declaration_array<integer>(n, 0, "index_visit");
	hash_table = my_declaration_array<bool>(n, false, "hash_table");
#endif

	// Минимальное количество неизвестных на грубой сетке
	integer MIN_COARSE_MESH_NODES = 50;

	if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) || (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T)) {
		if (n < 50) {
			MIN_COARSE_MESH_NODES = 5; //  27.06.2021
		}
	}

	while ((ilevel < maxlevel - 1) && (n_a[ilevel - 1] > MIN_COARSE_MESH_NODES) && (bcontinue_global)) {
		

		// защита от повторного срабатывания на добавление в интерполяции.
		from_re_operation_protection0 = true;
		ifrom_re_operation_protection = 0;		

		if (ilevel > 1) {

			if (my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS) {
				// HMIS
				// Hybrid Modifyed Independed Set.
				// RS на первом уровне, PMIS на остальных уровнях.
				bPMIS = true;
			}
			
			//threshold_quick_all
			my_realloc_memory<doublerealT>(threshold_quick_all, ((n_a[ilevel - 1]) + 1));
						
			//threshold_quick_only_negative
			my_realloc_memory<doublerealT>(threshold_quick_only_negative, ((n_a[ilevel - 1]) + 1));

			// istack
			my_realloc_memory<integer>(istack, ((n_a[ilevel - 1]) + 1));
			
			//hash_table2
			my_realloc_memory<bool>(hash_table2, ((n_a[ilevel - 1]) + 1));
			
			// this_is_C_node
			my_realloc_memory<bool>(this_is_C_node, ((n_a[ilevel - 1]) + 1));

			
			// this_is_F_node
			my_realloc_memory<bool>(this_is_F_node, ((n_a[ilevel - 1]) + 1));			
			
			// Выделяем оперативную память под хеш-таблицы экономно.
			construct_hash_table_Gus_struct01(n_a[ilevel - 1]+2);

		}

		// Константы размеров памяти 3.3 и 1.2 могут быть оспорены и изменены в последствии.
		// для этого требуется тестирование на большом числе рабочих задач.
		doublerealT dsize_memory_for_Amat = static_cast <doublerealT>(3.9); // на задачах с конвекцией тоже надо 3.9.

		if ((bmyconvective7248) && (my_amg_manager.number_interpolation_procedure == 2))
		{
			// Long range interpolation distance 2.
			if (iswitchMeshGenerator != CONFORMAL_MESH_GENERATOR_SELECTOR::COARSEMESHGEN_MESHER) {
				dsize_memory_for_Amat = static_cast <doublerealT>(5.9);
			}
		}


		if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION) ||
			((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION))) {
			// RS2 Активно. Джон Руге и Клаус Штубен второй проход.
			dsize_memory_for_Amat = static_cast <doublerealT>(3.9);
		}
		if (b_on_adaptive_local_refinement_mesh) {
			if ((1 && steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_STEADY)) {
				if (iVar == PAM) {
					// cfd для поправки давления на АЛИС сетке.
					//dsize_memory_for_Amat = 10;
					dsize_memory_for_Amat = static_cast <doublerealT>(4.9);
				}
			}
			else {
				dsize_memory_for_Amat = static_cast <doublerealT>(4.9);
			}
		}
		if (b_REALLOC) {
			// Уменьшение памяти отводимой под хранение матрицы А.
			// Матрица должна занимать в памяти не более чем под неё нужно и не мегабайтом больше.
			my_realloc_memory<real_mix_precision>(Amat.aij, ((iadd + 2 + static_cast<integer>(dsize_memory_for_Amat * nnz_a[ilevel - 1]))));
			
			my_realloc_memory<real_mix_precision>(Amat.abs_aij, ((iadd + 2 + static_cast<integer>(dsize_memory_for_Amat * nnz_a[ilevel - 1]))));
			
			my_realloc_memory<integer_mix_precision>(Amat.i, ((iadd + 2 + static_cast<integer>(dsize_memory_for_Amat * nnz_a[ilevel - 1]))));
			
			my_realloc_memory<integer_mix_precision>(Amat.j, ((iadd + 2 + static_cast<integer>(dsize_memory_for_Amat * nnz_a[ilevel - 1]))));
						
			if (bprint_mesage_diagnostic) {
				std::cout<<" 1 of 3 compleated.  OK!! ierarhion matrix Amat realloc successfully..."<<std::endl;
				std::cout<<"Prolongation ierarhion..."<<std::endl;
			}
			my_realloc_memory<Ak1>(P, ((nnz_P_memo_all + static_cast<integer>(1.2 * nnz_a[ilevel - 1])) + 2));
			
			if (bprint_mesage_diagnostic) {
				std::cout<<"2 of 3 compleated. OK!! ierarhion matrix Prolongation realloc successfully..."<<std::endl;
			}
		}
		

		if (ilevel > 1) {
			doublerealT procent_n = static_cast<doublerealT>((100.0*(abs(n_a[ilevel - 1] - n_a[ilevel - 2]))) / (1.0*n_a[ilevel - 2]));
			//doublerealT procent_nnz = (100.0 * (abs(nnz_a[ilevel - 1] - nnz_a[ilevel - 2]))) / (1.0 * nnz_a[ilevel - 2]);
			//if ((procent_n<2.0)&&(procent_nnz < 2.0)) break;
			if (procent_n < 2.0) break;
		}

		if (((ilevel > 1) && (nnz_a[ilevel - 1] > nnz_a[ilevel - 2]))) {
			//break;
		}

		// 19.04.2018
		print_control_volume_statistics(n_a, nnz_a, ilevel, bprint_mesage_diagnostic, debug_reshime);

		

#pragma omp parallel for
		for (integer ii = 1; ii <= n_a[ilevel - 1]; ++ii) {
			this_is_C_node[ii] = this_is_F_node[ii] = false;
		}

		const integer isize2 = (nnz_a[ilevel - 1] / 2);
		//integer ilab1 = nnz_a[ilevel - 1] + iadd + 1;

		// Сортировка нужна лишь на первом уровне, т.к.
		// результат алгоритма перемножения по Ф. Густавсону 1978 уже 
		// даёт на выходе отсортированную по строкам матрицу.
		if (ilevel == 1) {
			// сортировка исходной  А  по i.
			//heapsort(Amat, key=i*n_a[ilevel - 1] + j, 1, nnz_a[ilevel - 1]);

			int inth = number_cores();
			if (inth > 1) {
				omp_set_num_threads(inth);
			}
			//omp_set_num_threads(8);
			Ak1* Amattmp = nullptr;
		

			// 7 января 2016. Обязательно нужна эта сортировка.
			switch (imy_sort_algorithm) {
			case MY_SORT_ALGORITHM:: COUNTING_SORT:
				if (bmemory_savings == false) {
					// выделение памяти.
					Amattmp = new Ak1[nnz_a[ilevel - 1]];
					// прямое копирование.

#pragma omp parallel for
					for (integer i86 = 1 + iadd; i86 <= nnz_a[ilevel - 1] + iadd; ++i86) {
						integer j87 = i86 - 1 - iadd;
						Amattmp[j87].i = Amat.i[i86];
						Amattmp[j87].j = Amat.j[i86];
						Amattmp[j87].aij = Amat.aij[i86];
					}
					// сортировка

#pragma omp parallel for
					for (integer i = 0; i <= n_a[ilevel - 1]; ++i) {
						C1[i] = 0; // инициализация.
						//C2[i] = 0; // инициализация.
					}
					if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
						(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) || 
						(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
						(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
					{
						Counting_Sort_bmemo_false(Amattmp, 0, nnz_a[ilevel - 1] - 1, n_a[ilevel - 1], indx_comparei, C1, Bm1);// numberofcoarcenodes <-> n_a[ilevel - 1]
					}
					else {
						Counting_Sort_bmemo_false(Amattmp, 0, nnz_a[ilevel - 1] - 1, n_a[ilevel - 1], indx_comparei, C1);// numberofcoarcenodes <-> n_a[ilevel - 1]
					}

					// обратное копирование

#pragma omp parallel for
					for (integer i86 = 1 + iadd; i86 <= nnz_a[ilevel - 1] + iadd; ++i86) {
						integer j86 = i86 - 1 - iadd;
						Amat.i[i86] = Amattmp[j86].i;
						Amat.j[i86] = Amattmp[j86].j;
						Amat.aij[i86] = Amattmp[j86].aij;
					}
					// освобождение оперативной памяти.
					delete[] Amattmp;
				}
				else {

					Counting_Sort(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd, bmemory_savings, n_a[ilevel - 1]);	//подходит именно n_a[ilevel - 1]			
				}
				break;
			case MY_SORT_ALGORITHM::HEAP_SORT:

				if (inth == 1) {

					HeapSort(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					//LeftistHeapSort(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				}
				else {
#pragma omp parallel sections
					{
#pragma omp section
						{
							HeapSort(Amat, 1 + iadd, isize2 + iadd);
						}
#pragma omp section
						{
							HeapSort(Amat, isize2 + 1 + iadd, nnz_a[ilevel - 1] + iadd);
						}
					}
					mergeTim_amg(Amat, 1 + iadd, isize2 + iadd, nnz_a[ilevel - 1] + iadd);
				}
				break;
			case MY_SORT_ALGORITHM::QUICK_SORT:
				qs_abbys_heigh = 0;
				// quicksort
				if (0) {
					if (inth == 1) {

						qs(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd);

					}
					else {

#pragma omp parallel sections
						{
#pragma omp section
							{
								qs(Amat, 1 + iadd, isize2 + iadd);
							}
#pragma omp section
							{
								qs(Amat, isize2 + 1 + iadd, nnz_a[ilevel - 1] + iadd);
							}
						}
						mergeTim_amg(Amat, 1 + iadd, isize2 + iadd, nnz_a[ilevel - 1] + iadd);
					}
				}
				else {
					// выделение памяти.
					Amattmp = new Ak1[nnz_a[ilevel - 1]];
					// прямое копирование.
					
#pragma omp parallel for
					for (integer i86 = 1 + iadd; i86 <= nnz_a[ilevel - 1] + iadd; ++i86) {
						integer j87 = i86 - 1 - iadd;
						Amattmp[j87].i = Amat.i[i86];
						Amattmp[j87].j = Amat.j[i86];
						Amattmp[j87].aij = Amat.aij[i86];
					}
					// сортировка
					// Библиотечный алгоритм. O(n*log2(n)).
				    // Не использует лишней памяти.
					// Библиотечный std::sort намного быстрее собственной реализации QuickSort qs.
					// Распараллеленная версия библиотечного std::sort на два потока в 1.6 раза 
					// быстрее чем его последовательная версия (Проверено на размерности 10 млн значений в сортируемом массиве).
					// 08.06.2021
					if (inth == 1) {
						std::sort(Amattmp, Amattmp + nnz_a[ilevel - 1], compareAk1R);
					}
					else {
#pragma omp parallel sections
						{
#pragma omp section
							{
								std::sort(Amattmp, Amattmp + isize2 , compareAk1R);
							}
#pragma omp section
							{
								std::sort(Amattmp + isize2 + 1, Amattmp  + nnz_a[ilevel - 1], compareAk1R);
							}
						}
						mergeTim_amg(Amattmp, 0, isize2, nnz_a[ilevel - 1]-1, indx_comparei);
					}
					// обратное копирование
					
#pragma omp parallel for
					for (integer i86 = 1 + iadd; i86 <= nnz_a[ilevel - 1] + iadd; ++i86) {
						integer j86 = i86 - 1 - iadd;
						Amat.i[i86] = Amattmp[j86].i;
						Amat.j[i86] = Amattmp[j86].j;
						Amat.aij[i86] = Amattmp[j86].aij;
					}
					// освобождение оперативной памяти.
					delete[] Amattmp;
				}
				// Библиотечный алгоритм. O(n*log2(n)).
				// Не использует лишней памяти.
				//std::sort(Amat + (1 + iadd) * sizeof(Ak1), Amat + (nnz_a[ilevel - 1] + iadd + 1) * sizeof(Ak1), compAi);			

				

				//QuickSort(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				break;
			case MY_SORT_ALGORITHM::TIM_SORT:
				// Сортировка Тима Петерсома 2002г.
				if (inth == 1) {
					timSort_amg(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				}
				else {
#pragma omp parallel sections
					{
#pragma omp section
						{
							timSort_amg(Amat, 1 + iadd, isize2 + iadd);
						}
#pragma omp section
						{
							timSort_amg(Amat, isize2 + 1 + iadd, nnz_a[ilevel - 1] + iadd);
						}
					}
					mergeTim_amg(Amat, 1 + iadd, isize2 + iadd, nnz_a[ilevel - 1] + iadd);
				}
				break;
			default:
				Counting_Sort(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd, bmemory_savings, n_a[ilevel - 1]);//подходит именно n_a[ilevel - 1]
				break;
			}


			//omp_set_num_threads(1);

		} // ilevel == 1


		

		if (my_amg_manager.bMatrixPortrait) {
			// Печать портрета матрицы.

			PortraitPrint(Amat, iadd, nnz_a[ilevel - 1], n_a[ilevel - 1]);

			FILE* fp_portrait = NULL;
			
#ifdef MINGW_COMPILLER
			int err_portrait = 0;
			fp_portrait = fopen64("matrix_load.txt", "w");
			if (fp_portrait != NULL) {				
#if doubleintprecision == 1
				fprintf(fp_portrait, "%lld %lld\n", n_a[ilevel - 1], nnz_a[ilevel - 1]);
				for (integer i58 = 1 + iadd; i58 <= nnz_a[ilevel - 1] + iadd; ++i58) {
					fprintf(fp_portrait, "%d %d\n", Amat.i[i58], Amat.j[i58]);
				}
#else
				fprintf(fp_portrait, "%d %d\n", n_a[ilevel - 1], nnz_a[ilevel - 1]);
				for (integer i58 = 1 + iadd; i58 <= nnz_a[ilevel - 1] + iadd; ++i58) {
					fprintf(fp_portrait, "%d %d\n", Amat.i[i58], Amat.j[i58]);
				}
#endif
		    }
			else {
				err_portrait = 1;
			}
#else
			errno_t err_portrait = 0;
			err_portrait = fopen_s(&fp_portrait, "matrix_load.txt", "w");
#if doubleintprecision == 1
			fprintf_s(fp_portrait, "%lld %lld\n", n_a[ilevel - 1], nnz_a[ilevel - 1]);
			for (integer i58 = 1 + iadd; i58 <= nnz_a[ilevel - 1] + iadd; ++i58) {
				fprintf_s(fp_portrait, "%d %d\n", Amat.i[i58], Amat.j[i58]);
			}
#else
			fprintf_s(fp_portrait, "%d %d\n", n_a[ilevel - 1], nnz_a[ilevel - 1]);
			for (integer i58 = 1 + iadd; i58 <= nnz_a[ilevel - 1] + iadd; ++i58) {
				fprintf_s(fp_portrait, "%d %d\n", Amat.i[i58], Amat.j[i58]);
			}
#endif
#endif
			if (fp_portrait != NULL) {
				fclose(fp_portrait);
			}
			std::cout << "matrix portrait in level export\n";
			system("PAUSE");
		}
		
		integer ii_end1 = n_a[ilevel - 1];
#pragma omp parallel for
		for (integer ii = 1; ii <= ii_end1; ++ii) {
			flag[ii] = false;
		}

		
		// позиция начала каждой строки в матрице.
		if (row_startA != nullptr) {
			free(row_startA);
			row_startA = nullptr;
		}
		row_startA = my_declaration_array<integer>((n_a[ilevel - 1] + 1), 0, "row_startA");
		

		ii_end1 = nnz_a[ilevel - 1] + iadd;
		calculate_row_ptr(1 + iadd, ii_end1, row_startA, flag, Amat);


		ii_end1 = n_a[ilevel - 1];
#pragma omp parallel for
		for (integer ii = 1; ii <= ii_end1; ++ii) {
			flag[ii] = false;
		}

		if (bQuick_sort_for_reorder) {
			for (integer ii = 1; ii <= ii_end1; ++ii) {
				// 14.04.2020
				// В каждой строке i индексы столбца j отсортированы по возрастанию.		
				qs_abbys_heigh = 0;
				// quicksort
				qsj(Amat, row_startA[ii], row_startA[ii + 1] - 1);
			}
		}

		integer is_1 = row_startA[1];
		integer is_e = row_startA[n_a[ilevel - 1] + 1] - 1;
		// Заранее один раз вычисляем модуль элемента.
		if (is_e + 1 < ((iadd + 2 + static_cast<integer>(dsize_memory_for_Amat * nnz_a[ilevel - 1])))) {

			if (is_e + 1 < my_amg_manager.memory_size * nnz_a[0]) {

#pragma omp parallel for
				for (integer iscan = is_1; iscan <= is_e; ++iscan) {
					Amat.abs_aij[iscan] = fabs(Amat.aij[iscan]);
				}
			}
			else {
				std::cout << "overflow!!! You mast increase my_amg_manager.memory_size constant.\n";
				std::cout << "my_amg_manager.memory_size = " << my_amg_manager.memory_size << std::endl;
				system("PAUSE");
				exit(1);
			}
		}
		else {
			std::cout << "overflow!!! You mast increase dsize_memory_for_Amat constant.\n";
			system("PAUSE");
			exit(1);
		}

		// вычисляем для каждого узла число его соседей.
		// инициализация обязательна.
		if (count_neighbour != nullptr) {
			free(count_neighbour);
			count_neighbour = nullptr;
		}
		count_neighbour = my_declaration_array<integer>(n_a[ilevel - 1], 0, "count_neighbour"); // 0 - нет соседей.		

		if (bStrongTransposeON) 
		{
			// Освобождение ОЗУ.
			// Используем.
			//std::cout<<"usage strong transpose"<<std::endl;
			//system("PAUSE"); // debug
			
			// Эта ветвь активна лес АВЛ деревьев ненужен.

			// Обычный накопитель - линейный список с быстрой вставкой.
			/*
			if (hash_StrongTranspose_collection1 != nullptr) {

				//for (integer i_1 = 0; i_1 <= n_a[ilevel - 2]; ++i_1)
				//isize_memory_alloc_hash_StrongTranspose_collection1
				// Не надо операции с памятью распараллеливать т.к. они при
				// освобождении ресурса делается запрос к ядру ОС.
				// WARNING NO #pragma
				for (integer i_1 = 0; i_1 <= isize_memory_alloc_hash_StrongTranspose_collection1; ++i_1)
				{
					clear_list(hash_StrongTranspose_collection1[i_1]);
				}
				delete[] hash_StrongTranspose_collection1;
				hash_StrongTranspose_collection1 = nullptr;
			}
			if (isize_hash_StrongTranspose_collection != nullptr) {
				delete[] isize_hash_StrongTranspose_collection;
				isize_hash_StrongTranspose_collection = nullptr;
			}
			*/
			ii_end1 = n_a[ilevel - 1];

			// Выделяем память под лес линейных однонаправленных списков.
			/*hash_StrongTranspose_collection1 = new Taccumulqtor_list * [n_a[ilevel - 1] + 1];
#pragma omp parallel for
			for (integer i_1 = 0; i_1 <= ii_end1; ++i_1) {
				hash_StrongTranspose_collection1[i_1] = nullptr;
			}*/

			//isize_memory_alloc_hash_StrongTranspose_collection1 = ii_end1;
			/*if (isize_hash_StrongTranspose_collection != nullptr) {
				delete[] isize_hash_StrongTranspose_collection;
				isize_hash_StrongTranspose_collection = nullptr;
			}
			isize_hash_StrongTranspose_collection = new integer[n_a[ilevel - 1] + 1];
			*/

#pragma omp parallel for
			for (integer i_73 = 0; i_73 <= ii_end1; ++i_73) {
				isize_hash_StrongTranspose_collection[i_73] = 0;// Нет элементов.
			}
			
		}

		// При таком коде узел Дирихле тоже имеет соседа, сосед это 
		// внутренний узел который связан с этим узлом Дирихле.
		// Соседей вычисляем на самой первой матрице А (самой левой).
		ii_end1 = nnz_a[ilevel - 1] + iadd;

		if (bpositive_connections_CF_decomp) {

		
			integer ilen = (ii_end1 + 1);

#pragma omp parallel for
			for (integer ii = 1 + iadd; ii < ilen; ++ii) {
				bool cond = (flag[Amat.i[ii]] == false);
				if (!cond) continue;

				// Новейшая ветвь кода: 11.06.2017.
				// Введение новой ветви вызвано желанием ускорить код избегая повторных массовых вычислений threshold.
				// Ни в коем случае не ставить 0 в if.
				// Это новая единственно верная ветка. Её убирание приводит к неработоспособности всего приложения.
				threshold_quick_all[Amat.i[ii]] = -1.0;
				threshold_quick_only_negative[Amat.i[ii]] = -1.0;
				integer is0_end = row_startA[Amat.i[ii] + 1] - 1;


				for (integer is0 = ii; (is0 <= is0_end); ++is0) {
					bool cond_diagonal = (Amat.j[is0] == Amat.i[ii]);
					if (cond_diagonal) continue;


					if (Amat.abs_aij[is0] > threshold_quick_all[Amat.i[ii]]) {
						// Определяем максимальный внедиагональный элемент.
						threshold_quick_all[Amat.i[ii]] = Amat.abs_aij[is0];
					}

				}

				flag[Amat.i[ii]] = true;

			}

			//for (integer ii = 1 + iadd; ii <= ii_end1; ++ii) {
				//flag[Amat.i[ii]] = false;
			//}

			ilen = (n_a[ilevel - 1] + 1);

#pragma omp parallel for
			for (integer ii = 0; ii < ilen; ++ii) {
				flag[ii] = false;
		    }

			for (integer ii = 1 + iadd; ii <= ii_end1; ++ii) {
				bool cond = (flag[Amat.i[ii]] == false);
				if (!cond) continue;

				

				integer ic = -1;

				integer is0_end = row_startA[Amat.i[ii] + 1] - 1;

				//doublerealT theta_threshold3 = theta(ilevel)*threshold;
				doublerealT theta_threshold3 = static_cast<doublerealT>(theta(ilevel) * threshold_quick_all[Amat.i[ii]]);
				for (integer is0 = ii; (is0 <= is0_end); ++is0) {

					bool cond_diagonal = (Amat.j[is0] == Amat.i[ii]);
					if (cond_diagonal) continue;
					bool cond_strong_connection = (Amat.abs_aij[is0] > theta_threshold3);
					if (!cond_strong_connection) continue;


					// Учитываем только сильно связанных соседей.
					ic++; //i,j

					if (bStrongTransposeON) {
						// O(1) вставка в начало линейного списка.
						//insert_list(hash_StrongTranspose_collection1[Amat.j[is0]], Amat.i[ii]);
						hash_StrongTranspose_collection1Eco[isize_hash_StrongTranspose_collection[Amat.j[is0]]*n_a[ilevel-1]+ Amat.j[is0]] = Amat.i[ii];
						isize_hash_StrongTranspose_collection[Amat.j[is0]]++;
					}

				}

				

				count_neighbour[Amat.i[ii]] = ic;
				// 01.03.2019
				// Это делается ниже в 895 строке.
				//count_neighbour[Amat.i[ii]] = isize_hash_StrongTranspose_collection[Amat.i[ii]];
				// 22_12_2016
				if (ic == 0) {
					// Большой вопрос уместно ли так делать 8.апреля 2017 ???

					// До начала работы алгоритма все условия Дирихле становятся F узлами.
					this_is_C_node[Amat.i[ii]] = false;
					this_is_F_node[Amat.i[ii]] = true;

				}
				flag[Amat.i[ii]] = true;


			}

		}
		else {

			
			integer ilen = (ii_end1+1);

#pragma omp parallel for
			for (integer ii = 1 + iadd; ii < ilen; ++ii) {
				bool cond = (flag[Amat.i[ii]] == false);
				if (!cond) continue;

				// Новейшая ветвь кода: 11.06.2017.
				// Введение новой ветви вызвано желанием ускорить код избегая повторных массовых вычислений threshold.
				// Ни в коем случае не ставить 0 в if.
				// Это новая единственно верная ветка. Её убирание приводит к неработоспособности всего приложения.
				threshold_quick_all[Amat.i[ii]] = -1.0;
				threshold_quick_only_negative[Amat.i[ii]] = -1.0;
				integer is0_end = row_startA[Amat.i[ii] + 1] - 1;


				for (integer is0 = ii; (is0 <= is0_end); ++is0) {
					bool cond_diagonal = (Amat.j[is0] == Amat.i[ii]);
					if (cond_diagonal) continue;
					bool cond_negative_coefficient = (Amat.aij[is0] < 0.0);
					if (!cond_negative_coefficient) continue;


					if (Amat.abs_aij[is0] > threshold_quick_only_negative[Amat.i[ii]]) {
						// Определяем максимальный внедиагональный элемент.
						threshold_quick_only_negative[Amat.i[ii]] = Amat.abs_aij[is0];
					}

				}


				flag[Amat.i[ii]] = true;
			}

			//for (integer ii = 1 + iadd; ii <= ii_end1; ++ii) {
				//flag[Amat.i[ii]] = false;
			//}

			ilen = (n_a[ilevel - 1] + 1);

#pragma omp parallel for
			for (integer ii = 0; ii < ilen; ++ii) {
				flag[ii] = false;
			}

			for (integer ii = 1 + iadd; ii <= ii_end1; ++ii) {
				bool cond = (flag[Amat.i[ii]] == false);
				if (!cond) continue;

				integer ic = -1;

				
				integer is0_end = row_startA[Amat.i[ii] + 1] - 1;

				for (integer is0 = ii; (is0 <= is0_end); ++is0) {

					bool cond_diagonal = (Amat.j[is0] == Amat.i[ii]);
					if (cond_diagonal) continue;
					bool cond_negative_coefficient = (Amat.aij[is0] < 0.0);
					if (!cond_negative_coefficient) continue;
					doublereal barrier = theta(ilevel) * threshold_quick_only_negative[Amat.i[ii]];
					bool cond_strong_connection = (Amat.abs_aij[is0] > barrier);
					if (!cond_strong_connection) continue;


					// Учитываем только сильно связанных соседей.
					ic++; //i,j

					if (bStrongTransposeON) {
						// O(1) вставка в начало линейного списка.
						//insert_list(hash_StrongTranspose_collection1[Amat.j[is0]], Amat.i[ii]);
						hash_StrongTranspose_collection1Eco[isize_hash_StrongTranspose_collection[Amat.j[is0]] * n_a[ilevel - 1] + Amat.j[is0]] = Amat.i[ii];
						isize_hash_StrongTranspose_collection[Amat.j[is0]]++;
					}
				}
				




				count_neighbour[Amat.i[ii]] = ic;
				// 01.03.2019
				// Это делается ниже в 895 строке.
				//count_neighbour[Amat.i[ii]] = isize_hash_StrongTranspose_collection[Amat.i[ii]];
				// 22_12_2016
				if (ic == 0) {
					// Большой вопрос уместно ли так делать 8.апреля 2017 ???

					// До начала работы алгоритма все условия Дирихле становятся F узлами.
					this_is_C_node[Amat.i[ii]] = false;
					this_is_F_node[Amat.i[ii]] = true;

				}
				flag[Amat.i[ii]] = true;
			}
		}

		//integer imaximum_neighbour = 0;
		//for (integer i_73 = 0; i_73 <= n_a[ilevel - 1]; ++i_73) {
			//if (isize_hash_StrongTranspose_collection[i_73] > imaximum_neighbour) {
				//imaximum_neighbour = isize_hash_StrongTranspose_collection[i_73];
			//}
		//}

		integer imaximum_neighbour = get_max_array_elm_integer(isize_hash_StrongTranspose_collection, n_a[ilevel - 1]+1);

		//std::cout << imaximum_neighbour << "  " << imaximum_neighbour1 << std::endl;
		//system("pause");

		if (bprint_mesage_diagnostic) {
			std::cout << "\nilevel=" << ilevel << " imaximum_neighbour = " << imaximum_neighbour  << "\n";
		}

		if (bStrongTransposeON) {
			// 5.01.2017. StrongTranspose.
			// Счётчик lambda инициализирован согласно литературным описаниям через Strong Transpose.
			ii_end1 = n_a[ilevel - 1];

#pragma omp parallel for
			for (integer i_1 = 1; i_1 <= ii_end1; ++i_1) {

				// 20.05.2017 Добавлен быстрый доступ по ключу для количества элементов в дереве.
				//count_neighbour[i_1] = getnumber_AVL_node_global(hash_StrongTranspose_collection[i_1]);
				count_neighbour[i_1] = isize_hash_StrongTranspose_collection[i_1];
				if (count_neighbour[i_1] == 0) {
					// 14.04.2017 Важнейшая положительная модификация 
					// сокращающая количество итераций:
					// # задача; число ит. до; число ит. после;
					// 1. passiv_module6; 179; 97;
					// 2. CGHV 12mm HFET; 18, 8, 6, 3, 2; 17, 8, 6, 3, 2;
					// 3. PIONER; 77; 73;

					// До начала работы алгоритма все условия Дирихле становятся F узлами.
					this_is_C_node[i_1] = false;
					this_is_F_node[i_1] = true;

				}			
			}
		}


		
		
		ii_end1 = n_a[ilevel - 1];

#pragma omp parallel for
		for (integer ii = 1; ii <= ii_end1; ++ii) {
			flag[ii] = false; // init flag
		}

		

		if (bPMIS) {

			if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) &&
				(ilevel >= 0)) {

				// PMIS применённый к квадрату матрицы.

				PMIS_CF_decomposition_applied_to_the_square_of_the_matrix(Amat, this_is_F_node,
					this_is_C_node, ilevel,
					count_neighbour, n_a,
					newCcount,
					threshold_quick_only_negative,
					row_startA,
					hash_StrongTranspose_collection1Eco,
					isize_hash_StrongTranspose_collection);

				/*

				const int imax_aggressive_pass = 2;

				PMIS_CF_decomposition_applied_to_the_square_of_the_matrix_aggressive(Amat, this_is_F_node,
					this_is_C_node, ilevel,
					count_neighbour, n_a,
					newCcount,
					threshold_quick_only_negative,
					row_startA,
					hash_StrongTranspose_collection1Eco,
					isize_hash_StrongTranspose_collection, imax_aggressive_pass);
					*/
			}
			else
			if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) && 
				(ilevel >= 4)) {
				     // PMIS применённый к квадрату матрицы (Это что то типо аналога aggressive coarsening).
					 // Помоему имеет смысл применять aggresive coarsening только на "плотных" уровнях с большим шаблоном.
					 // 14.09.2021 Я решил применять начиная с четвёртого уровня.

					 // Применено к квадрату матрицы А.

					// Разбиение множества узлов первоначальной сетки на С-узлы и 
					// F - узлы. С -узлы составят грубую сетку следующего уровня вложенности.
					// Значения же функции в F узлах должно быть восстановлено по значению функции
					// в ближайших С узлах (см. задачу интерполяции).
					// начало написания 20.09.2020.
				if (1) {
					//my_amg_manager.truncation_interpolation_Temperature = 0.9;

					PMIS_CF_decomposition_applied_to_the_square_of_the_matrix(Amat, this_is_F_node,
						this_is_C_node, ilevel,
						count_neighbour, n_a,
						newCcount,
						threshold_quick_only_negative,
						row_startA,
						hash_StrongTranspose_collection1Eco,
						isize_hash_StrongTranspose_collection);
				}
				else {
					const int imax_aggressive_pass = 2;

					// Так тоже работает, но требует дополнительного изучения.

					PMIS_CF_decomposition_applied_to_the_square_of_the_matrix_aggressive(Amat, this_is_F_node,
						this_is_C_node, ilevel,
						count_neighbour, n_a,
						newCcount,
						threshold_quick_only_negative,
						row_startA,
						hash_StrongTranspose_collection1Eco,
						isize_hash_StrongTranspose_collection, imax_aggressive_pass);
				}

				}
			else {

				if (bPMIS_applied_to_the_square_of_the_matrix) {
					// Применено к квадрату матрицы А.

					// Разбиение множества узлов первоначальной сетки на С-узлы и 
					// F - узлы. С -узлы составят грубую сетку следующего уровня вложенности.
					// Значения же функции в F узлах должно быть восстановлено по значению функции
					// в ближайших С узлах (см. задачу интерполяции).
					// начало написания 20.09.2020.
					PMIS_CF_decomposition_applied_to_the_square_of_the_matrix(Amat, this_is_F_node,
						this_is_C_node, ilevel,
						count_neighbour, n_a,
						newCcount,
						threshold_quick_only_negative,
						row_startA,
						hash_StrongTranspose_collection1Eco,
						isize_hash_StrongTranspose_collection);

				}
				else {
					// Применено к матрице А.

					// Разбиение множества узлов первоначальной сетки на С-узлы и 
					// F - узлы. С -узлы составят грубую сетку следующего уровня вложенности.
					// Значения же функции в F узлах должно быть восстановлено по значению функции
					// в ближайших С узлах (см. задачу интерполяции). 
					PMIS_CF_decomposition(Amat, this_is_F_node,
						this_is_C_node, ilevel,
						count_neighbour, n_a,
						newCcount,
						threshold_quick_only_negative,
						row_startA,
						hash_StrongTranspose_collection1Eco,
						isize_hash_StrongTranspose_collection);
				}
			}
		}
		else {
			
			integer max_neighbour = 0;
			integer icandidate = 0;

			// Находим узел с наибольшим числом соседей и запоминаем его.
			// Это первый встретившийся узел с наибольшим числом соседей.
			// Это требуется для того чтобы стартовал алгоритм C/F разбиения.

			

#ifdef _OPENMP

			omp_set_num_threads(number_cores()); // установка числа потоков

			// Многопоточный код.
			max_neighbour = get_max_array_elm_integer1(count_neighbour, ii_end1);

			bool found63 = false;
#pragma omp parallel for
			for (integer i7 = 1; i7 <= ii_end1; ++i7) {
				if ((!found63)&&(count_neighbour[i7] == max_neighbour)) {

#pragma omp critical
					{
						found63 = true;
						icandidate = row_startA[i7];
					}

				}
			}

#else

			// Однопоточный код.
			for (integer i7 = 1; i7 <= ii_end1; ++i7) {
				if (count_neighbour[i7] > max_neighbour) {
					max_neighbour = count_neighbour[i7];
					icandidate = row_startA[i7];
				}
			}

#endif

			// 4 июля 2016.
	        // это случай когда следующий уровень вложенности просто не из чего строить и это 
	        // становится понятно только здесь.
			if ((icandidate == 0) && (max_neighbour == 0)) {

				//system("pause");
				// уровень построить нельзя поэтому досрочный выход из цикла.
				break;
			}

			// Разбиение множества узлов первоначальной сетки на С-узлы и 
            // F - узлы. С -узлы составят грубую сетку следующего уровня вложенности.
			// Значения же функции в F узлах должно быть восстановлено по значению функции
			// в ближайших С узлах (см. задачу интерполяции). 
			Ruge_and_Stuben_CF_decomposition<doublerealT>(Amat, this_is_F_node,
				this_is_C_node, ilevel,
				count_neighbour, n_a,
				newCcount,
				threshold_quick_only_negative,
				row_startA,
				hash_StrongTranspose_collection1Eco,
				isize_hash_StrongTranspose_collection,
				bprint_mesage_diagnostic,
				flag, iadd, bpositive_connections_CF_decomp,
				bStrongTransposeON, hash_table2,
				istack, nnz_a, debug_reshime);
		}

		  //delete[] count_neighbour;
		if (count_neighbour != nullptr) {
			free(count_neighbour);
			count_neighbour = nullptr;
		}

		

		if (bprint_mesage_diagnostic) {
			if (n_a[ilevel - 1] == 0) {
				std::cout<<"n_a is zero"<<std::endl;
				system("pause");
			}
			std::cout<<"additional C="<<static_cast<doublerealT>(100.0*newCcount / n_a[ilevel - 1])<<"%"<<std::endl;
			//system("pause");
		}

		
			
		if (!bPMIS) {

			// Добавление новых С узлов для удовлетворения правил RS интерполяции.

			// В методе стандартной интерполяции присутствует шаг уменьшения разреженности,
			// для того чтобы правильно аппроксимировать все F переменные C переменными надо
			// увеличить количество С переменных.
			the_number_of_neighbors_that_are_not_C_nodes = 0;
			number_of_F_nodes_with_one_single_strong_C_neighbor = 0;
			number_of_F_nodes_with_one_single_strong_C_neighborF = 0;

			iadditionalCstatistic = 0;

			bweSholdbeContinue = true;
			while (bweSholdbeContinue) {
				bweSholdbeContinue = false;

				bool* bvacant_candidates = nullptr;
				bvacant_candidates = new bool[n_a[ilevel - 1] + 1];
				if (bvacant_candidates == nullptr) {
					std::cout<<"error memory alloc in bvacant_candidates"<<std::endl;
					system("pause");
					exit(1);
				}
#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					bvacant_candidates[i_1] = false;
				}

				// Построение пролонгации для узлов которые составляют F-nodes.
				// Каждый F-nodes окружён C-nodes.
#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					bool cond1 = (this_is_F_node[i_1]);
					if (!cond1) continue;

					
					// Найти соседей данного F-node которые C-node.
					integer icsos = 0;
					// старая версия до 10 января 2016.
					//integer i_2 = BinarySearchAi(Amat, i_1, 1 + iadd, nnz_a[ilevel - 1] + iadd);
					// Быстрый вариант без поиска, просто индексирование на основе "хеш таблицы".
					// 10 января 2016. на основе хеширования.
					integer i_2 = row_startA[i_1];

					bool bvisit = false;
					//for (integer is0 = i_2; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[i_2]); ++is0) {
					integer iend_merker_position = row_startA[Amat.i[i_2] + 1] - 1;
					for (integer is0 = i_2; (is0 <= iend_merker_position); ++is0) {
						bool cond_diagonal = (Amat.j[is0] == Amat.i[i_2]);
						if (cond_diagonal) continue;

						
						bvisit = true;
						if (this_is_C_node[Amat.j[is0]]) {
							icsos++;
						}
						else {
							//the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
						}
						
					}
					//if (icsos == 1) number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным С соседом.
					// Если bvisit то внедиагональные элементы есть но они все F-nodes. Иначе там обособленное условие Дирихле.
					if ((icsos == 0) && (bvisit)) {

						// А если он F узел Дирихле без соседей, то сумма тоже может быть нулевой и это вызовет деление на ноль.
						// Узлы Дирихле могли быть без соседей на начальных уровнях, они располагались в конце списка и были
						// поглощены агломератами внутренних узлов и всё было в порядке.
						// Чтобы преодолеть это затруднение нужен алгоритм с обратной связью.

						// Нет С соседей, этот узел станет С узлом.
						bvacant_candidates[i_1] = true;
					}
					
				}


#ifdef _OPENMP
				// Параллельное исполнение не более чем в 40 потоков
				const unsigned int MAX_THREAD_LOC = 40;
				integer newCcount_arr[MAX_THREAD_LOC];
				integer the_number_of_neighbors_that_are_not_C_nodes_arr[MAX_THREAD_LOC];
				integer number_of_F_nodes_with_one_single_strong_C_neighbor_arr[MAX_THREAD_LOC];
				bool bweSholdbeContinue_arr[MAX_THREAD_LOC];

#pragma omp parallel for
				for (integer i_1 = 0; i_1 < MAX_THREAD_LOC; ++i_1) {
					newCcount_arr[i_1] = 0;
					the_number_of_neighbors_that_are_not_C_nodes_arr[i_1] = 0;
					number_of_F_nodes_with_one_single_strong_C_neighbor_arr[i_1] = 0;
					bweSholdbeContinue_arr[i_1] = false;
				}

				// Построение пролонгации для узлов которые составляют F-nodes.
				// Каждый F-nodes окружён C-nodes.
#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1)
				{
					if (this_is_F_node[i_1]) {

#ifdef _OPENMP 
						int tid = omp_get_thread_num();
#else
						int tid = 0;
#endif

						// Найти соседей данного F-node которые C-node.
						integer icsos = 0;
						// старая версия до 10 января 2016.
						//integer i_2 = BinarySearchAi(Amat, i_1, 1 + iadd, nnz_a[ilevel - 1] + iadd);
						// Быстрый вариант без поиска, просто индексирование на основе "хеш таблицы".
						// 10 января 2016. на основе хеширования.
						integer i_2 = row_startA[i_1];

						bool bvisit = false;
						//for (integer is0 = i_2; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[i_2]); ++is0) {
						integer iend_merker_position = row_startA[Amat.i[i_2] + 1] - 1;
						for (integer is0 = i_2; (is0 <= iend_merker_position); ++is0) {
							if (Amat.j[is0] != Amat.i[i_2]) {
								bvisit = true;
								if (this_is_C_node[Amat.j[is0]]) {
									icsos++;
								}
								else {
									//the_number_of_neighbors_that_are_not_C_nodes++; // подсчитываем проблемы интерполяции 
									the_number_of_neighbors_that_are_not_C_nodes_arr[tid]++;
								}
							}
						}
						if (icsos == 1) {
							//	number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным С соседом.
							number_of_F_nodes_with_one_single_strong_C_neighbor_arr[tid]++;
						}
						// Если bvisit то внедиагональные элементы есть но они все F-nodes. Иначе там обособленное условие Дирихле.
						if ((icsos == 0) && (bvisit)) {

							// А если он F узел Дирихле без соседей, то сумма тоже может быть нулевой и это вызовет деление на ноль.
							// Узлы Дирихле могли быть без соседей на начальных уровнях, они располагались в конце списка и были
							// поглощены агломератами внутренних узлов и всё было в порядке.
							// Чтобы преодолеть это затруднение нужен алгоритм с обратной связью.							

							// Нет С соседей, этот узел станет С узлом.
							this_is_F_node[i_1] = false;
							this_is_C_node[i_1] = true;
							// F-node стал C-node!!! Идея стандартной интерполяции 
							// приводит к уменьшению разреженности оператора Галёркина.
							bweSholdbeContinue_arr[tid] = true;
							newCcount_arr[tid]++;
						}

						// 1 января 2015 Один сосед это недостаточно.
						// Поэтому в случае одного соседа делаем такой узел С узлом.
						if ((false) && (icsos == 1)) {
							// bvisit и так true т.к. icsos==1.
							this_is_F_node[i_1] = false;
							this_is_C_node[i_1] = true;
							//bweSholdbeContinue = true;
							bweSholdbeContinue_arr[tid] = true;
						}

					}

				}

				for (integer i_1 = 0; i_1 < MAX_THREAD_LOC; ++i_1) {
					newCcount += newCcount_arr[i_1];
					the_number_of_neighbors_that_are_not_C_nodes += the_number_of_neighbors_that_are_not_C_nodes_arr[i_1];
					number_of_F_nodes_with_one_single_strong_C_neighbor += number_of_F_nodes_with_one_single_strong_C_neighbor_arr[i_1];
					if (bweSholdbeContinue_arr[i_1]) {
						bweSholdbeContinue = true;
					}
				}
#else

				// Параллельное исполнение не более чем в MAX_THREAD_LOC потоков
				integer newCcount_arr = 0;
				integer the_number_of_neighbors_that_are_not_C_nodes_arr = 0;
				integer number_of_F_nodes_with_one_single_strong_C_neighbor_arr = 0;
				bool bweSholdbeContinue_arr = false;


				// Построение пролонгации для узлов которые составляют F-nodes.
				// Каждый F-nodes окружён C-nodes.
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1)
				{
					if (this_is_F_node[i_1]) {


						// Найти соседей данного F-node которые C-node.
						integer icsos = 0;
						// старая версия до 10 января 2016.
						//integer i_2 = BinarySearchAi(Amat, i_1, 1 + iadd, nnz_a[ilevel - 1] + iadd);
						// Быстрый вариант без поиска, просто индексирование на основе "хеш таблицы".
						// 10 января 2016. на основе хеширования.
						integer i_2 = row_startA[i_1];

						bool bvisit = false;
						//for (integer is0 = i_2; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[i_2]); ++is0) {
						integer iend_merker_position = row_startA[Amat.i[i_2] + 1] - 1;
						for (integer is0 = i_2; (is0 <= iend_merker_position); ++is0) {
							if (Amat.j[is0] != Amat.i[i_2]) {
								bvisit = true;
								if (this_is_C_node[Amat.j[is0]]) {
									icsos++;
								}
								else {
									//++the_number_of_neighbors_that_are_not_C_nodes; // подсчитываем проблемы интерполяции 
									++the_number_of_neighbors_that_are_not_C_nodes_arr;
								}
							}
						}
						if (icsos == 1) {
							//	number_of_F_nodes_with_one_single_strong_C_neighbor++; // количество F узлов с одним единственным сильным С соседом.
							++number_of_F_nodes_with_one_single_strong_C_neighbor_arr;
						}
						// Если bvisit то внедиагональные элементы есть но они все F-nodes. Иначе там обособленное условие Дирихле.
						if ((icsos == 0) && (bvisit)) {

							// А если он F узел Дирихле без соседей, то сумма тоже может быть нулевой и это вызовет деление на ноль.
							// Узлы Дирихле могли быть без соседей на начальных уровнях, они располагались в конце списка и были
							// поглощены агломератами внутренних узлов и всё было в порядке.
							// Чтобы преодолеть это затруднение нужен алгоритм с обратной связью.							

							// Нет С соседей, этот узел станет С узлом.
							this_is_F_node[i_1] = false;
							this_is_C_node[i_1] = true;
							// F node стал C_node!!! Идея стандартной интерполяции 
							// приводит к уменьшению разреженности оператора Галёркина.
							bweSholdbeContinue_arr = true;
							++newCcount_arr;
						}

						// 1 января 2015 Один сосед это недостаточно.
						// Поэтому в случае одного соседа делаем такой узел С узлом.
						if ((false) && (icsos == 1)) {
							// bvisit и так true т.к. icsos==1.
							this_is_F_node[i_1] = false;
							this_is_C_node[i_1] = true;
							//bweSholdbeContinue = true;
							bweSholdbeContinue_arr = true;
						}

					}

				}

				newCcount += newCcount_arr;
				the_number_of_neighbors_that_are_not_C_nodes += the_number_of_neighbors_that_are_not_C_nodes_arr;
				number_of_F_nodes_with_one_single_strong_C_neighbor += number_of_F_nodes_with_one_single_strong_C_neighbor_arr;
				if (bweSholdbeContinue_arr) {
					bweSholdbeContinue = true;
				}

#endif



				if (bprint_mesage_diagnostic) {

					std::cout << "newCcount=" << newCcount <<", n_a=" << n_a[ilevel - 1] <<" "<< 100.0 * newCcount / n_a[ilevel - 1] << "\n";

				}
				
				delete[] bvacant_candidates;
				

				if (bprint_mesage_diagnostic) {
					if (bweSholdbeContinue) {
						std::cout<<" prohod succseful: bweSholdbeContinue==true\n";
					}
					else {
						std::cout<<" prohod empty: bweSholdbeContinue=false\n";
					}
				}

			}

			// ***

			// 01.01.2017 Алгоритм улучшения качества C/F разбиения. Проход 2. RS2
			// Цикл по всем F переменным, полученным после первого прохода.
			// Пусть Fi текущая F переменная и у неё множество соседей не пусто.
			// Сканируем строку элементов где Fi есть диагональный элемент.
			// Amat. Определяем порог - threshold для каждой строки.
			// В. Заносим всех сильных С соседей в специальный линейный список.
			// C. Если мы встретили сильного F соседа  (Fj), так что Fi и Fj сильно связаны,
			// то ищем всех сильных С соседей узла Fj и формируем из них линейный список.
			// С помощью алгоритма слияния за линейное время сравниваем два предварительно отсортированных линейных
			// списка на предмет общих С узлов.
			// D. Если общий С узел есть то ничего не меняем.
			// E. Если общего сильного С узла не обнаружено то один из узлов Fi или Fj становится С узлом.
			// Среди Fi и Fj тот становится С узлом у которого больше сильных F соседей. Если С узлом стал Fj 
			// то линейный список С соседей узла Fi обновляется. Если С узлом стал узел Fi то мы заканчиваем обработку Fi 
			// возвращая всех помеченных Fj снова в F тип.
			//  30.12.2016
			// 11.06.2017 Здесь для сортировки используется библиотечный std::sort на массиве.
			if (1) {
				if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION) ||
					((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION))) { // RS2 Проход 2.

					if (bprint_mesage_diagnostic) {
						std::cout << "   ***   CAMG SELECTOR RS2 " << ilevel << "  ***\n"; 
					}


					for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) if (this_is_F_node[i_1]) {
						// i_1 это F переменная Fi.
						//Amat.Определяем порог - threshold для каждой строки.
						doublerealT thresholdRS = -1.0;
						integer i_2 = row_startA[i_1];

						// Очистка хеш-таблицы.
						clear_hash_table_Gus_struct01();
						// занесение данных из линейного списка в хеш-таблицу для дерева с корнем в Amat.i[i_2].
						//integer imarker75_scan = 0;
						//formirate_F_SiTranspose_hash_table_Gus_struct02(hash_StrongTranspose_collection1[Amat.i[i_2]], imarker75_scan);

						integer iend_merker_position = row_startA[Amat.i[i_2] + 1] - 1;
						if (!btreshold_on_new_vetv) {
							for (integer is0 = i_2; (is0 <= iend_merker_position); ++is0) {
								if (Amat.j[is0] != Amat.i[i_2]) {
									if (Amat.aij[is0] < 0.0) {
										if (Amat.abs_aij[is0] > thresholdRS) thresholdRS = Amat.abs_aij[is0];
									}
								}
							}
						}
						else {
							// Новейшая ветвь кода: 11.06.2017.
							thresholdRS = threshold_quick_only_negative[Amat.i[i_2]];
						}
						if (thresholdRS > 0.0) {
							// Множество соседей не пусто а порог равен thresholdRS.
							hashlist_i* ivacant_F2C = nullptr;
							//  В. Заносим всех сильных С соседей в специальный линейный список.
							hashlist_i* ibuffer_strongC = nullptr;
							integer ibuffer_strongC_marker = -1;
							integer inumber_strongF_count_Fi = 0;
							hashlist_i* ibuffer_strongF = nullptr;
							integer ibuffer_strongF_marker = -1;
							for (integer is0 = i_2; (is0 <= iend_merker_position); ++is0) {
								bool cond_diagonal = (Amat.j[is0] == Amat.i[i_2]);
								if (cond_diagonal) continue;
								bool cond_negative = (Amat.aij[is0] < 0.0);
								if (!cond_negative) continue;
								bool cond_strong_connection = (fabs(Amat.aij[is0]) > theta(ilevel) * thresholdRS);
								if (!cond_strong_connection) continue;

								
								
								
								if (this_is_C_node[Amat.j[is0]]) {
									++ibuffer_strongC_marker;
									insertion_list_i(ibuffer_strongC, Amat.j[is0]);
									insert_hash_table_Gus_struct01(Amat.j[is0]);// 11.08.2018
								}
								if (this_is_F_node[Amat.j[is0]]) {

									//if (1) 19.01.2017
									if (0) {// if (0) 11.08.2018
										// Добавок 19.01.2017
										
										// Внимание hash_StrongTranspose_collection должна быть инициализирована сначала
										// в общем настроена для использования, а этого по видимому не сделано т.к. используется
										//hash_StrongTranspose_collection1. Кстати операции isfound для hash_StrongTranspose_collection1
										// нет т.к. она очень медленная т.к. он просто линейный список.

										if (hash_StrongTranspose_collection != nullptr) {
											data_BalTreeST dat_key;
											dat_key.i = Amat.j[is0];
											if (isfound(hash_StrongTranspose_collection[Amat.i[is0]], dat_key)) {
												// конец добавка 19.01.2017
												
												// Сильный Fj сосед найден.
												// Элементы Fi и Fj сильно связаны.
												++inumber_strongF_count_Fi;
												++ibuffer_strongF_marker;
												insertion_list_i(ibuffer_strongF, Amat.j[is0]);
											}
										}
										else {
											// Сильный Fj сосед найден.
											// Элементы Fi и Fj сильно связаны.
											++inumber_strongF_count_Fi;
											++ibuffer_strongF_marker;
											insertion_list_i(ibuffer_strongF, Amat.j[is0]);
										}

									}
									else {
										// Сильный Fj сосед найден.
										// Элементы Fi и Fj сильно связаны.
										++inumber_strongF_count_Fi;
										++ibuffer_strongF_marker;

										insertion_list_i(ibuffer_strongF, Amat.j[is0]);
									}
								}																	
								
							}
							// Очистка хеш-таблицы.
							clear_hash_table_Gus_struct01();
							// Сортировка буффера ibuffer_strongC по возрастанию.
							// рекомендуется использовать iusage_old_version = 0
							// при котором активируется использование быстродействующей хеш-таблицы.
							// Достигается ускорение полного цикла решения задачи при включённом RS2 coarsening
							// на 7.5% по сравнению с двоичным поиском на массиве. 
							// Полностью отпадает необходимость в использовании алгоритма сортировки.
							// 11.06.2017.
							//integer iusage_old_version = 0; // 1 старая рабочая версия. // 0 новая версия на основе хеш таблицы.
							
							
							// Вместо сортировки и двоичного поиска используем хеш-таблицу.
							clear_hash_table_Gus_struct01();
							hashlist_i* ibuffer_strongC_scan = ibuffer_strongC;
							while (ibuffer_strongC_scan != nullptr) {
								insert_hash_table_Gus_struct01(ibuffer_strongC_scan->item);
								ibuffer_strongC_scan = ibuffer_strongC_scan->next;
							}
							ibuffer_strongC_scan = nullptr;
							

							
							// Все сильные F-соседи занесены в буфер ibuffer_strongF. 
							hashlist_i* ibuffer_strongF_current = ibuffer_strongF;
							for (integer i_3 = 0; i_3 <= ibuffer_strongF_marker; ++i_3) {
								if (ibuffer_strongF_current != nullptr) {
									// Сканируем всех сильных F соседей последовательно.
									//1. Определяем threshold для Fj.
									doublerealT thresholdRS1 = -1.0;
									integer i_4 = row_startA[ibuffer_strongF_current->item];
									integer iend_merker_position1 = row_startA[Amat.i[i_4] + 1] - 1;
									if (!btreshold_on_new_vetv) {
										for (integer is01 = i_4; (is01 <= iend_merker_position1); ++is01) {
											if (Amat.j[is01] != Amat.i[i_4]) {
												if (Amat.aij[is01] < 0.0) {
													if (fabs(Amat.aij[is01]) > thresholdRS1) thresholdRS1 = fabs(Amat.aij[is01]);
												}
											}
										}
									}
									else {
										// Новейшая ветвь кода: 11.06.2017.
										thresholdRS1 = threshold_quick_only_negative[Amat.i[i_4]];
									}
									integer inumber_strongF_count_Fj = 0;
									// искомый порог thresholdRS1.

									hashlist_i* ibuffer_strongCFj = nullptr;
									integer ibuffer_strongCFj_marker = -1;
									for (integer is01 = i_4; (is01 <= iend_merker_position1); ++is01) {
										if (Amat.j[is01] != Amat.i[i_4]) {
											if (Amat.aij[is01] < 0.0) {
												if (fabs(Amat.aij[is01]) > theta(ilevel)* thresholdRS1) {
													if (this_is_C_node[Amat.j[is01]]) {
														++ibuffer_strongCFj_marker;
														insertion_list_i(ibuffer_strongCFj, Amat.j[is01]);
													}
													if (this_is_F_node[Amat.j[is01]]) {
														++inumber_strongF_count_Fj;
													}
												}
											}
										}
									}
									// В ibuffer_strongCFj список сильных С соседей.

									// Есть ли общие С узлы за линейное время.
									// Создаём на основе списка ibuffer_strongC
									// целочисленный массив.
									// Сортируем его. Делаем  ibuffer_strongCFj_marker 
									// двоичных поисков в этом отсортированном массиве 
									// до тех пор пока не встретится успешный поиск.
									bool bfound_32 = false;
									hashlist_i* ibuffer_strongCFj_scan = ibuffer_strongCFj;
									
									// Версия на основе хеш-таблицы.
									while ((bfound_32 == false) && (ibuffer_strongCFj_scan != nullptr)) {
										// Совпадение найдено мы ничего не делаем.
										bfound_32 = isfound_hash_table_Gus_struct01(ibuffer_strongCFj_scan->item);
										ibuffer_strongCFj_scan = ibuffer_strongCFj_scan->next;
									}
									
									ibuffer_strongCFj_scan = nullptr;
									if (bfound_32 == false) {
										// Один из них станет С узлом.
										if ((ibuffer_strongF_current->item > i_1) && (inumber_strongF_count_Fj >= inumber_strongF_count_Fi)) {
											// Если Fj находится в ещё непросмотренной части списка F узлов и
											// у него по сравнению с F узлом Fi больше сильных F связей.								

											// Fj становится С.
											insertion_list_i(ivacant_F2C, ibuffer_strongF_current->item);
											this_is_C_node[ibuffer_strongF_current->item] = true;
											this_is_F_node[ibuffer_strongF_current->item] = false;
											ibuffer_strongC_marker++;
											inumber_strongF_count_Fi--;
											insertion_list_i(ibuffer_strongC, ibuffer_strongF_current->item);
											
											hashlist_i* ibuffer_strongC_scan_2 = ibuffer_strongC;

											// Очищаем хеш-таблицу и заполняем её по новой.
											clear_hash_table_Gus_struct01();
											while (ibuffer_strongC_scan_2 != nullptr) {
												insert_hash_table_Gus_struct01(ibuffer_strongC_scan_2->item);
												ibuffer_strongC_scan_2 = ibuffer_strongC_scan_2->next;
											}

											ibuffer_strongC_scan_2 = nullptr;
											
										}
										else {
											// Fi становится С.
											this_is_C_node[i_1] = true;
											this_is_F_node[i_1] = false;
											// Возвращаем все Fj с С на F.
											hashlist_i* ivacant_F2C_marker = ivacant_F2C;
											while (ivacant_F2C_marker != nullptr) {
												this_is_F_node[ivacant_F2C_marker->item] = true;
												this_is_C_node[ivacant_F2C_marker->item] = false;
												ivacant_F2C_marker = ivacant_F2C_marker->next;
											}
											ivacant_F2C_marker = nullptr;
											if (ivacant_F2C != nullptr) {
												clear_hash_list_i(ivacant_F2C);
												ivacant_F2C = nullptr;
											}
											ivacant_F2C = nullptr;
											// Очищаем ОЗУ.
											if (ibuffer_strongCFj != nullptr) {
												clear_hash_list_i(ibuffer_strongCFj);
												ibuffer_strongCFj = nullptr;
											}
											ibuffer_strongCFj = nullptr;
											// Досрочно прерываем текущее сканирование 
											// списка сильных F узлов.
											break;
										}

									}


									// Очищаем ОЗУ.
									if (ibuffer_strongCFj != nullptr) {
										clear_hash_list_i(ibuffer_strongCFj);
										ibuffer_strongCFj = nullptr;
									}
									ibuffer_strongCFj = nullptr;
									// Переход к следующему кандидату.
									ibuffer_strongF_current = ibuffer_strongF_current->next;
								}
							}

							// Освобождение ОЗУ.
							if (ibuffer_strongC != nullptr) {
								clear_hash_list_i(ibuffer_strongC);
								ibuffer_strongC = nullptr;
							}
							ibuffer_strongC = nullptr;
							if (ibuffer_strongF != nullptr) {
								clear_hash_list_i(ibuffer_strongF);
								ibuffer_strongF = nullptr;
							}
							ibuffer_strongF = nullptr;
							if (ivacant_F2C != nullptr) {
								clear_hash_list_i(ivacant_F2C);
								ivacant_F2C = nullptr;
							}
							ivacant_F2C = nullptr;

						}
					}
				} // Алгоритм улучшения качества C/F разбиения. Проход 2.
			}
			else {
				// мертвая экспериментальная ветвь. Можно удалить.

				// 01.01.2017 Алгоритм улучшения качества C/F разбиения. Проход 2. 
				// Цикл по всем F переменным, полученным после первого прохода.
				// Пусть Fi текущая F переменная и у неё множество соседей не пусто.
				// Сканируем строку элементов где Fi есть диагональный элемент.
				// Amat. Определяем порог - threshold для каждой строки.
				// В. Заносим всех сильных С соседей в специальный линейный список.
				// C. Если мы встретили сильного F соседа  (Fj), так что Fi и Fj сильно связаны,
				// то ищем всех сильных С соседей узла Fj и формируем из них линейный список.
				// С помощью алгоритма слияния за линейное время сравниваем два предварительно отсортированных линейных
				// списка на предмет общих С узлов.
				// D. Если общий С узел есть то ничего не меняем.
				// E. Если общего сильного С узла не обнаружено то один из узлов Fi или Fj становится С узлом.
				// Среди Fi и Fj тот становится С узлом у которого больше сильных F соседей. Если С узлом стал Fj 
				// то линейный список С соседей узла Fi обновляется. Если С узлом стал узел Fi то мы заканчиваем обработку Fi 
				// возвращая всех помеченных Fj снова в F тип.
				//  30.12.2016
				if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ALL_CONNECTION) ||
					((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::RS2_ST_ALL_CONNECTION))) { // RS2 Проход 2.
					for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) if (this_is_F_node[i_1]) {
						// i_1 это F переменная Fi.
						//Amat.Определяем порог - threshold для каждой строки.
						doublerealT thresholdRS = -1.0;
						integer i_2 = row_startA[i_1];

						// Очистка хеш-таблицы.
						clear_hash_table_Gus_struct01();
						// занесение данных из линейного списка в хеш-таблицу для дерева с корнем в Amat.i[i_2].
						//!!! устарело 27.09.2020 integer imarker75_scan = 0;
						//!!! устарело 27.09.2020 formirate_F_SiTranspose_hash_table_Gus_struct02(hash_StrongTranspose_collection1[Amat.i[i_2]], imarker75_scan);


						integer iend_merker_position = row_startA[Amat.i[i_2] + 1] - 1;
						if (!btreshold_on_new_vetv) {
							for (integer is0 = i_2; (is0 <= iend_merker_position); ++is0) {
								if (Amat.j[is0] != Amat.i[i_2]) {
									if (Amat.aij[is0] < 0.0) {
										if (fabs(Amat.aij[is0]) > thresholdRS) thresholdRS = fabs(Amat.aij[is0]);
									}
								}
							}
						}
						else {
							// Новейшая ветвь кода: 11.06.2017.
							thresholdRS = threshold_quick_only_negative[Amat.i[i_2]];
						}
						if (thresholdRS > 0.0) {
							// Множество соседей не пусто а порог равен thresholdRS.
							hashlist_i* ivacant_F2C = nullptr;
							//  В. Заносим всех сильных С соседей в специальный линейный список.
							hashlist_i* ibuffer_strongC = nullptr;
							integer ibuffer_strongC_marker = -1;
							integer inumber_strongF_count_Fi = 0;
							hashlist_i* ibuffer_strongF = nullptr;
							integer ibuffer_strongF_marker = -1;
							for (integer is0 = i_2; (is0 <= iend_merker_position); ++is0) {
								if (Amat.j[is0] != Amat.i[i_2]) {
									if (Amat.aij[is0] < 0.0) {
										if (fabs(Amat.aij[is0]) > theta(ilevel)* thresholdRS) {
											if (this_is_C_node[Amat.j[is0]]) {
												ibuffer_strongC_marker++;
												insertion_list_i(ibuffer_strongC, Amat.j[is0]);
												insert_hash_table_Gus_struct01(Amat.j[is0]);// 11.08.2018
											}
											if (this_is_F_node[Amat.j[is0]]) {

												//if (1) 19.01.2017
												if (1) {// if (0) 11.08.2018
														// Добавок 19.01.2017

													//if (hash_StrongTranspose_collection1 != nullptr)
													{
														//data_BalTreeST dat_key;
														//dat_key.i = Amat.j[is0];
														bool bfound71 = false;
														for (integer i72 = 0; i72 < isize_hash_StrongTranspose_collection[Amat.i[i_2]]; ++i72) {
															if (Amat.j[is0] == hash_StrongTranspose_collection1Eco[i72*n_a[ilevel-1]+ Amat.i[i_2]]) {
																bfound71 = true;
																break;
															}
														}
														//if (isfound(hash_StrongTranspose_collection1[Amat.i[i_2]], Amat.j[is0])) 
														if (bfound71)
														{
															// конец добавка 19.01.2017
															// Сюда почему-то вообще не заходит код исполнения??? 

															// Сильный Fj сосед найден.
															// Элементы Fi и Fj сильно связаны.
															inumber_strongF_count_Fi++;
															ibuffer_strongF_marker++;
															insertion_list_i(ibuffer_strongF, Amat.j[is0]);
														}
													}
													/*else {
														// Сильный Fj сосед найден.
														// Элементы Fi и Fj сильно связаны.
														inumber_strongF_count_Fi++;
														ibuffer_strongF_marker++;
														insertion_list_i(ibuffer_strongF, Amat.j[is0]);
													}*/

												}
												else {
													// Сильный Fj сосед найден.
													// Элементы Fi и Fj сильно связаны.
													inumber_strongF_count_Fi++;
													ibuffer_strongF_marker++;

													insertion_list_i(ibuffer_strongF, Amat.j[is0]);
												}
											}
										}
									}
								}
							}
							
							// Сортировка буфера ibuffer_strongC по возрастанию.
							// рекомендуется использовать iusage_old_version = 0
							// при котором активируется использование быстродействующей хеш-таблицы.
							// Достигается ускорение полного цикла решения задачи при включённом RS2 coarsening
							// на 7.5% по сравнению с двоичным поиском на массиве. 
							// Полностью отпадает необходимость в использовании алгоритма сортировки.
							// 11.06.2017.
							//integer iusage_old_version = 0; // 1 старая рабочая версия. // 0 новая версия на основе хеш-таблицы.
							

							// Вместо сортировки и двоичного поиска используем хеш-таблицу.
							// Очистка хеш-таблицы.
							clear_hash_table_Gus_struct01();
							hashlist_i* ibuffer_strongC_scan = ibuffer_strongC;
							while (ibuffer_strongC_scan != nullptr) {
								insert_hash_table_Gus_struct01(ibuffer_strongC_scan->item);
								ibuffer_strongC_scan = ibuffer_strongC_scan->next;
							}
							ibuffer_strongC_scan = nullptr;


							// Сортировка целочисленного массива при индексации с нуля!!!

							// Все сильные F соседи занесены в буфер ibuffer_strongF. 
							hashlist_i* ibuffer_strongF_current = ibuffer_strongF;
							for (integer i_3 = 0; i_3 <= ibuffer_strongF_marker; ++i_3) {
								if (ibuffer_strongF_current != nullptr) {
									// Сканируем всех сильных F соседей последовательно.
									//1. Определяем threshold для Fj.
									doublerealT thresholdRS1 = -1.0;
									integer i_4 = row_startA[ibuffer_strongF_current->item];
									integer iend_merker_position1 = row_startA[Amat.i[i_4] + 1] - 1;
									if (!btreshold_on_new_vetv) {
										for (integer is01 = i_4; (is01 <= iend_merker_position1); ++is01) {
											if (Amat.j[is01] != Amat.i[i_4]) {
												if (Amat.aij[is01] < 0.0) {
													if (fabs(Amat.aij[is01]) > thresholdRS1) thresholdRS1 = fabs(Amat.aij[is01]);
												}
											}
										}
									}
									else {
										// Новейшая ветвь кода: 11.06.2017.
										thresholdRS1 = threshold_quick_only_negative[Amat.i[i_4]];
									}
									integer inumber_strongF_count_Fj = 0;
									// искомый порог thresholdRS1.

									hashlist_i* ibuffer_strongCFj = nullptr;
									integer ibuffer_strongCFj_marker = -1;
									for (integer is01 = i_4; (is01 <= iend_merker_position1); ++is01) {
										if (Amat.j[is01] != Amat.i[i_4]) {
											if (Amat.aij[is01] < 0.0) {
												if (fabs(Amat.aij[is01]) > theta(ilevel)* thresholdRS1) {
													if (this_is_C_node[Amat.j[is01]]) {
														ibuffer_strongCFj_marker++;
														insertion_list_i(ibuffer_strongCFj, Amat.j[is01]);
													}
													if (this_is_F_node[Amat.j[is01]]) {
														inumber_strongF_count_Fj++;
													}
												}
											}
										}
									}

									// В ibuffer_strongCFj список сильных С соседей.

									// Есть ли общие С узлы за линейное время.
									// Создаём на основе списка ibuffer_strongC
									// целочисленный массив.
									// Сортируем его. Делаем  ibuffer_strongCFj_marker 
									// двоичных поисков в этом отсортированном массиве 
									// до тех пор пока не встретится успешный поиск.
									bool bfound_32 = false;
									hashlist_i* ibuffer_strongCFj_scan = ibuffer_strongCFj;
									
									// Версия на основе хеш-таблицы.
									while ((bfound_32 == false) && (ibuffer_strongCFj_scan != nullptr)) {
										// Совпадение найдено мы ничего не делаем.
										bfound_32 = isfound_hash_table_Gus_struct01(ibuffer_strongCFj_scan->item);
										ibuffer_strongCFj_scan = ibuffer_strongCFj_scan->next;
									}
									
									ibuffer_strongCFj_scan = nullptr;

									if (bfound_32 == false) {
										// Один из них станет С узлом.
										if ((ibuffer_strongF_current->item > i_1) && (inumber_strongF_count_Fj >= inumber_strongF_count_Fi)) {
											// Если Fj находится в ещё не просмотренной части списка F узлов и
											// у него по сравнению с F узлом Fi больше сильных F связей.								

											// Fj становится С.
											insertion_list_i(ivacant_F2C, ibuffer_strongF_current->item);
											this_is_C_node[ibuffer_strongF_current->item] = true;
											this_is_F_node[ibuffer_strongF_current->item] = false;
											ibuffer_strongC_marker++;
											inumber_strongF_count_Fi--;
											insertion_list_i(ibuffer_strongC, ibuffer_strongF_current->item);											

											hashlist_i* ibuffer_strongC_scan_1 = ibuffer_strongC;
											
											// Очищаем хеш-таблицу и заполняем её по новой.
											clear_hash_table_Gus_struct01();
											while (ibuffer_strongC_scan_1 != nullptr) {
												insert_hash_table_Gus_struct01(ibuffer_strongC_scan_1->item);
												ibuffer_strongC_scan_1 = ibuffer_strongC_scan_1->next;
											}
											
											ibuffer_strongC_scan_1 = nullptr;
											// Сортировка целочисленного массива при индексации с нуля!!!

										}
										else {
											// Fi становится С.
											this_is_C_node[i_1] = true;
											this_is_F_node[i_1] = false;
											// Возвращаем все Fj с С на F.
											hashlist_i* ivacant_F2C_marker = ivacant_F2C;
											while (ivacant_F2C_marker != nullptr) {
												this_is_F_node[ivacant_F2C_marker->item] = true;
												this_is_C_node[ivacant_F2C_marker->item] = false;
												ivacant_F2C_marker = ivacant_F2C_marker->next;
											}
											ivacant_F2C_marker = nullptr;
											if (ivacant_F2C != nullptr) {
												clear_hash_list_i(ivacant_F2C);
												ivacant_F2C = nullptr;
											}
											ivacant_F2C = nullptr;
											// Очищаем ОЗУ.
											if (ibuffer_strongCFj != nullptr) {
												clear_hash_list_i(ibuffer_strongCFj);
												ibuffer_strongCFj = nullptr;
											}
											ibuffer_strongCFj = nullptr;
											// Досрочно прерываем текущее сканирование 
											// списка сильных F узлов.
											break;
										}

									}


									// Очищаем ОЗУ.
									if (ibuffer_strongCFj != nullptr) {
										clear_hash_list_i(ibuffer_strongCFj);
										ibuffer_strongCFj = nullptr;
									}
									ibuffer_strongCFj = nullptr;
									// Переход к следующему кандидату.
									ibuffer_strongF_current = ibuffer_strongF_current->next;
								}
							}

							// Освобождение ОЗУ.
							if (ibuffer_strongC != nullptr) {
								clear_hash_list_i(ibuffer_strongC);
								ibuffer_strongC = nullptr;
							}
							ibuffer_strongC = nullptr;
							if (ibuffer_strongF != nullptr) {
								clear_hash_list_i(ibuffer_strongF);
								ibuffer_strongF = nullptr;
							}
							ibuffer_strongF = nullptr;
							if (ivacant_F2C != nullptr) {
								clear_hash_list_i(ivacant_F2C);
								ivacant_F2C = nullptr;
							}
							ivacant_F2C = nullptr;

						}
					}
				} // Алгоритм улучшения качества C/F разбиения. Проход 2.
			}


		}
		// Освобождаем оперативную память как только можем. 
		// Как только C/F разбиение построено данные hash_StrongTranspose_collection1 и аналог уже не используются.
		if (bStrongTransposeON) {
			// Освобождение ОЗУ.

			
			// Обычный линейный список.
			/*
			if (hash_StrongTranspose_collection1 != nullptr) {
				//for (integer i_1 = 0; i_1 <= n_a[ilevel - 2]; ++i_1)
				//isize_memory_alloc_hash_StrongTranspose_collection1
				for (integer i_1 = 0; i_1 <= isize_memory_alloc_hash_StrongTranspose_collection1; ++i_1)
				{
					if (hash_StrongTranspose_collection1[i_1] != NULL) {
						if (hash_StrongTranspose_collection1[i_1] != nullptr) {
							clear_list(hash_StrongTranspose_collection1[i_1]);
						}
					}
				}
				delete[] hash_StrongTranspose_collection1;
				hash_StrongTranspose_collection1 = nullptr;
			}*/

			/*
			if (isize_hash_StrongTranspose_collection != nullptr) {
				delete[] isize_hash_StrongTranspose_collection;
				isize_hash_StrongTranspose_collection = nullptr;
			}*/
		}

		// Нужно корректно обработать узлы Дирихле,
		// Если F узел окажется узлом Дирихле без соседей то его надо сделать С узлом,
		// Но узнать такой узел можно лишь в процессе выполнения алгоритма дальше по ходу исполнения.
		// Поэтому может потребоваться вернуться и начать заново (обратная связь).

		if (C_numerate != nullptr) {
			free(C_numerate);
			C_numerate = nullptr;
		}
		C_numerate = my_declaration_array<integer>(n_a[ilevel - 1], 0, "C_numerate");

		icounter = 1;
		ap_coarse = nullptr;

		// Мысль в том чтобы избавится от перезапуска и сделать всё за один проход, но это 
		// сделает код менее понятным и главное ухудшит силу интерполяции.
		//bool no_FeedBack = true;
		

		if (bprint_mesage_diagnostic) {
			std::cout << "   ***   CAMG PROLONGATOR "<< ilevel <<"  ***\n";
		}

		if (!bPMIS) {

			bool bsuffix_work = true;

			bweSholdbeContinue = true;
			while (bweSholdbeContinue) {
				bweSholdbeContinue = false;

				integer n_coarce = 0;
#pragma omp parallel for reduction(+:n_coarce)
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					if (this_is_C_node[i_1]) {
						++n_coarce;
					}
				}
				
				// debug
				// проверка качества C/F разбиения.
				//doublerealT* exp1 = new doublerealT[n_a[ilevel - 1] + 1];
				//for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) exp1[i_1] = 0.0;
				//for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) if (this_is_C_node[i_1]) exp1[i_1] = 2.0;
				//for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) if (this_is_F_node[i_1]) exp1[i_1] = 1.0;
				//exporttecplot(exp1,n);
				//delete[] exp1;

				//std::cout<<"export ready";

				//system("pause");



				// C/F разбиение построено, самое время построить оператор интерполяции.
				// потом найти оператор проекции, как транспонированный оператор интерполяции.
				// Всё завершает построение матрицы нового сеточного уровня и можно запускать новый уровень.

				// Построение оператора интерполяции: 
				// coarse 2 fine.
				//P*coarse==fine


				// Занумеруем (упорядочим) узлы грубой сетки.
#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) C_numerate[i_1] = 0;
				icounter = 1;


				// 29.12.2021 ошибочное распараллеливание в строках 10730, 10826, 11033, 11125.
/*
#ifdef _OPENMP

				const int iSIZE_THD = number_cores();

				omp_set_num_threads(iSIZE_THD); // установка числа потоков

				integer* icounter_tid = new integer[iSIZE_THD];
				integer* icounter_tidSize = new integer[iSIZE_THD];

#pragma omp parallel for
				for (int i_1 = 0; i_1 < iSIZE_THD; ++i_1) {
					icounter_tid[i_1] = 1;
					icounter_tidSize[i_1] = 0;
				}

				integer* threadid = new integer[n_a[ilevel - 1]+1];

#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					threadid[i_1] = -1;
				}

#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					if (this_is_C_node[i_1]) {

						int tid = omp_get_thread_num();

						threadid[i_1] = tid;

						C_numerate[i_1] = icounter_tid[tid];

						++icounter_tid[tid];
					}
				}

				 

				// Очень короткий последовательный цикл.
				for (int i_tid = 1; i_tid < iSIZE_THD; ++i_tid)
				{
					if (i_tid == 1) {
						icounter_tidSize[i_tid-1]= icounter_tid[i_tid - 1] - 1;
					}
					else {
						icounter_tidSize[i_tid - 1] = icounter_tidSize[i_tid - 2] + icounter_tid[i_tid - 1] - 1;
					}
				}

					

#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					if (this_is_C_node[i_1]) {
						if (threadid[i_1] > 0) {
							C_numerate[i_1] += icounter_tidSize[threadid[i_1]-1];
						}
					}
				}
				

				icounter = icounter_tidSize[iSIZE_THD - 2] + icounter_tid[iSIZE_THD - 1];


				delete[] threadid;
				delete[] icounter_tid;
				delete[] icounter_tidSize;


				//for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) if (this_is_C_node[i_1]) {

					//std::cout<<"C ind= "<< i_1; system("PAUSE");

					//C_numerate[i_1] = icounter;
					//++icounter;
				//}

				//std::cout << icounter87 << " " << icounter << "\n";
				//system("pause");
#else
*/
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) if (this_is_C_node[i_1]) {

					//std::cout<<"C ind= "<< i_1; system("PAUSE");

					C_numerate[i_1] = icounter;
					++icounter;
				}

//#endif				

				// C_numerate - перенумерация на множестве Coarse узлов.
				// Построение пролонгации для узлов которые составляют грубую сетку.
				
// 29.12.2021 ошибочное распараллеливание в строках 10730, 10826, 11033, 11125.
/*#ifdef _OPENMP

				// Многопоточная версия 29.12.2021.

				icount1 = 0;

#pragma omp parallel for reduction(+: icount1)
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) 
				{
					if (this_is_C_node[i_1]) {

						integer icount62 = C_numerate[i_1] + iaddR;

						P[icount62].aij = 1.0;
						P[icount62].i = static_cast<integer_mix_precision>(C_numerate[i_1]); // coarse number
						P[icount62].j = static_cast<integer_mix_precision>(i_1); // fine number.
						++icount1;
						if (icount62 +1  >= nsizePR * n) {

#pragma omp critical
							{
								std::cout << "memory error!!!" << std::endl;
								std::cout << "not enough memory for the interpolation operator." << std::endl;
								system("PAUSE");
								//exit(1);
								deallocate_prolongation(nsizePR, n, P);
							}
						}
					}
				}

				icount1 += 1 + iaddR; 

#else*/

				// Однопоточная версия 29.12.2021.

				icount1 = 1 + iaddR; // nnz_R

				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1)
				{
					if (this_is_C_node[i_1]) {

						

						P[icount1].aij = 1.0;
						P[icount1].i = static_cast<integer_mix_precision>(C_numerate[i_1]); // coarse number
						P[icount1].j = static_cast<integer_mix_precision>(i_1); // fine number.
						++icount1;
						if (icount1 >= nsizePR * n) {
#
								std::cout << "memory error!!!" << std::endl;
								std::cout << "not enough memory for the interpolation operator." << std::endl;
								system("PAUSE");
								//exit(1);
								deallocate_prolongation(nsizePR, n, P);
						}
					}
				}

//#endif

				// значение icount1 нужно далее.НЕ трогать !!!.
				numberofcoarcenodes = icount1 - 1 - iaddR;



				// Для модификации R  надо transpose(P)/ap.
				if (bprint_mesage_diagnostic) {

					std::cout << "number of coarce nodes=" << numberofcoarcenodes << "\n"; 

					if (debug_reshime) system("pause");
				}


				if (ap_coarse != nullptr) {
					free(ap_coarse);
					ap_coarse = nullptr;
				}
				ap_coarse = my_declaration_array<doublerealT>(numberofcoarcenodes, 0.0, "ap_coarse");



				// Для каждого С узла запоминаем в ap_coarse[C_numerate[i8]] 
				// модуль диагонального элемента.
#pragma omp parallel for
				for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8) {
					if (this_is_C_node[i8]) {
						// Старая версия до 10 января 2016. Время O(log2(nnz))
						//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
						// 10 января 2016 новая версия на основе хеширования. Время O(1).
						integer ii1 = row_startA[i8];
						// бинарный поиск должен гарантирует нахождение самого левого представителя.
						//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); ++is0) {
						integer iend_marker_position_1 = row_startA[Amat.i[ii1] + 1] - 1;
						for (integer is0 = ii1; (is0 <= iend_marker_position_1); ++is0) {
							//std::cout<<"i="<< Amat.i[is0] <<" j=" << Amat.j[is0] <<" Amat.aij[is0]="<< Amat.aij[is0]<<" "<< std::endl;
							
							if (Amat.j[is0] == Amat.i[ii1]) {

								if (fabs(Amat.aij[is0]) > RealMAXIMUM) {
									std::cout<<"overflow error #1: fabs(Amat.aij["<< is0 <<"]) > RealMAXIMUM \n";
									system("pause");
								}
								ap_coarse[C_numerate[i8]] = fabs(Amat.aij[is0]);
								//std::cout<<"find = "<< fabs(Amat.aij[is0]) << std::endl;
							}
						}
					}
					//std::cout<<std::endl;
					//system("pause");
				}

				//std::cout<<"incoming="<< my_amg_manager.number_interpolation_procedure<< std::endl;

			//system("pause");


				

				// верно 2 октября.	
				doublerealT magic82f = static_cast<doublerealT>(magic82);

				my_interpolation_procedure_universal<doublerealT>(the_number_of_neighbors_that_are_not_C_nodes,
						number_of_F_nodes_with_one_single_strong_C_neighbor,
						n_a, this_is_F_node, row_startA,
						nnz_a, bpositive_connections, Amat,
						bweSholdbeContinue, this_is_C_node, iadditionalCstatistic,
						RealZERO, icount1, P, nsizePR, ilevel,
						iadd, n, C_numerate,
						number_of_F_nodes_with_one_single_strong_C_neighborF,
						btreshold_on_new_vetv, ifrom_re_operation_protection,
						from_re_operation_protection0, magic82f, threshold_quick_all,
						threshold_quick_only_negative, bsuffix_work);
								

				if (bprint_mesage_diagnostic) {
					std::cout<<"Additional C nodes in interpolation procedure. Statistics:"<< std::endl;
					std::cout << "number firstable C nodes=" << n_coarce <<", number secondary C nodes="<< iadditionalCstatistic <<"\n";
				}

				if (bweSholdbeContinue) {

					// На каждом уровне делается всего два прохода.
					// На первом проходе добавляются все необходимые недостающие С
					// узлы. Экспериментально установлено (опытным путем) что 
					// второй проход не добавляет новых С узлов и поэтому проверку
					// на добавление можно не делать. Отказ от проверки экономит флопы процессора и 
					// увеличивает быстродействие для Якоби интерполяции, для amg1r5 интерполяции 
					// ускорения не обнаружено.
					// bsuffix_work = false; - отказ от проверки.
					// bsuffix_work = true; - проверку всё равно делать.
					bsuffix_work = false; // recomended false for Jacoby interpolation.
					//delete[] ap_coarse;
					if (ap_coarse != nullptr) {
						free(ap_coarse);
						ap_coarse = nullptr;
					}
					if (bprint_mesage_diagnostic) {
						std::cout<<"Feedback restart..."<< std::endl;
					}
				}

				if (bprint_mesage_diagnostic) {
					// отношение добавленных узлов к количеству С узлов на предыдущем уровне.
					//std::cout << "addition C nodes "<< static_cast<doublerealT>(100.0*iadditionalCstatistic / n_a[ilevel - 1])<<"%\n";
					// отношение количества добавленных С узлов к первоначальному количеству С узлов на данном уровне.
					std::cout << "addition C nodes = "<< static_cast<doublerealT>(100.0 * iadditionalCstatistic / n_coarce) <<"% firstable C nodes,  level population "<<static_cast<doublerealT>(100.0 * (n_coarce + iadditionalCstatistic) / n_a[ilevel - 1]) <<"%\n";
				}
				iadditionalCstatistic = 0;
				//system("pause");


			}
		}
		else {
		    //интерполяция PMIS в том числе и дальнобойный вариант distance=3.
			
			bool bsuffix_work = true;

			bweSholdbeContinue = true;
			while (bweSholdbeContinue) {
				bweSholdbeContinue = false;

				integer n_coarce = 0;
#pragma omp parallel for reduction(+:n_coarce)
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					if (this_is_C_node[i_1]) {
						n_coarce++;
					}
				}
				
				// C/F разбиение построено, самое время построить оператор интерполяции.
				// потом найти оператор проекции, как транспонированный оператор интерполяции.
				// Всё завершает построение матрицы нового сеточного уровня и можно запускать новый уровень.

				// Построение оператора интерполяции: 
				// coarse -> fine (2 == ->).
				// P*coarse==fine


				// Занумеруем (упорядочим) узлы грубой сетки.
#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) C_numerate[i_1] = 0;
				icounter = 1;

				// 29.12.2021 ошибочное распараллеливание в строках 10730, 10826, 11033, 11125.
/*#ifdef _OPENMP

				const int iSIZE_THD = number_cores();

				omp_set_num_threads(iSIZE_THD); // установка числа потоков

				integer* icounter_tid = new integer[iSIZE_THD];
				integer* icounter_tidSize = new integer[iSIZE_THD];

#pragma omp parallel for
				for (int i_1 = 0; i_1 < iSIZE_THD; ++i_1) {
					icounter_tid[i_1] = 1;
					icounter_tidSize[i_1] = 0;
				}

				integer* threadid = new integer[n_a[ilevel - 1] + 1];

#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					threadid[i_1] = -1;
				}

#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					if (this_is_C_node[i_1]) {

						int tid = omp_get_thread_num();

						threadid[i_1] = tid;

						C_numerate[i_1] = icounter_tid[tid];

						++icounter_tid[tid];
					}
				}



				// Очень короткий последовательный цикл.
				for (int i_tid = 1; i_tid < iSIZE_THD; ++i_tid)
				{
					if (i_tid == 1) {
						icounter_tidSize[i_tid - 1] = icounter_tid[i_tid - 1] - 1;
					}
					else {
						icounter_tidSize[i_tid - 1] = icounter_tidSize[i_tid - 2] + icounter_tid[i_tid - 1] - 1;
					}
				}



#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					if (this_is_C_node[i_1]) {
						if (threadid[i_1] > 0) {
							C_numerate[i_1] += icounter_tidSize[threadid[i_1] - 1];
						}
					}
				}


				icounter = icounter_tidSize[iSIZE_THD - 2] + icounter_tid[iSIZE_THD - 1];


				delete[] threadid;
				delete[] icounter_tid;
				delete[] icounter_tidSize;


				//for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) if (this_is_C_node[i_1]) {

					//std::cout<<"C ind= "<< i_1; system("PAUSE");

					//C_numerate[i_1] = icounter;
					//++icounter;
				//}

				//std::cout << icounter87 << " " << icounter << "\n";
				//system("pause");
#else*/

				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) if (this_is_C_node[i_1]) {

					//std::cout<<"C ind= "<< i_1; system("PAUSE");

					C_numerate[i_1] = icounter;
					++icounter;
				}

//#endif

// 29.12.2021 ошибочное распараллеливание в строках 10730, 10826, 11033, 11125.				
/*#ifdef _OPENMP
				// Многопоточная версия 29.12.2021.

				icount1 = 0;

#pragma omp parallel for reduction(+: icount1)
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1)
				{
					if (this_is_C_node[i_1]) {

						integer icount62 = C_numerate[i_1] + iaddR;

						P[icount62].aij = 1.0;
						P[icount62].i = static_cast<integer_mix_precision>(C_numerate[i_1]); // coarse number
						P[icount62].j = static_cast<integer_mix_precision>(i_1); // fine number.
						++icount1;
						if (icount62 + 1 >= nsizePR * n) {

#pragma omp critical
							{
								std::cout << "memory error!!!" << std::endl;
								std::cout << "not enough memory for the interpolation operator." << std::endl;
								system("PAUSE");
								//exit(1);
								deallocate_prolongation(nsizePR, n, P);
							}
						}
					}
				}

				icount1 += 1 + iaddR;
#else*/

				// Однопоточная версия 29.12.2021.

				// C_numerate - перенумерация на множестве Coarse узлов.
				// Построение пролонгации для узлов которые составляют грубую сетку.
				icount1 = 1 + iaddR; // nnz_R
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) if (this_is_C_node[i_1]) {
					P[icount1].aij = 1.0;
					P[icount1].i = static_cast<integer_mix_precision>(C_numerate[i_1]); // coarse number
					P[icount1].j = static_cast<integer_mix_precision>(i_1); // fine number.
					++icount1;
					if (icount1 >= nsizePR * n) {
						std::cout << "memory error!!!"<< std::endl;
						std::cout << "not enough memory for the interpolation operator."<< std::endl;
						//system("PAUSE");
						//exit(1);
						deallocate_prolongation(nsizePR, n, P);
					}
				}

//#endif


				// значение icount1 нужно далее.НЕ трогать !!!.
				numberofcoarcenodes = icount1 - 1 - iaddR;



				// Для модификации R  надо transpose(P)/ap.
				if (bprint_mesage_diagnostic) {

					std::cout << "number of coarce nodes=" << numberofcoarcenodes << "\n";

					if (debug_reshime) system("pause");
				}


				if (ap_coarse != nullptr) {
					free(ap_coarse);
					ap_coarse = nullptr;
				}
				ap_coarse = my_declaration_array<doublerealT>(numberofcoarcenodes, 0.0, "ap_coarse");



				// Для каждого С узла запоминаем в ap_coarse[C_numerate[i8]] 
				// модуль диагонального элемента.
#pragma omp parallel for
				for (integer i8 = 1; i8 <= n_a[ilevel - 1]; ++i8) {
					if (this_is_C_node[i8]) {
						// Старая версия до 10 января 2016. Время O(log2(nnz))
						//integer ii1 = BinarySearchAi(Amat, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
						// 10 января 2016 новая версия на основе хеширования. Время O(1).
						integer ii1 = row_startA[i8];
						// бинарный поиск должен гарантирует нахождение самого левого представителя.
						//for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (Amat.i[is0] == Amat.i[ii1]); ++is0) {
						integer iend_marker_position_1 = row_startA[Amat.i[ii1] + 1] - 1;
						for (integer is0 = ii1; (is0 <= iend_marker_position_1); ++is0) {
							//std::cout << "i="<<Amat.i[is0]<<" j="<<Amat.j[is0]<<" Amat.aij[is0]="<< Amat.aij[is0]<<" "<<std::endl;

							if (Amat.j[is0] == Amat.i[ii1]) {

								if (Amat.abs_aij[is0] > RealMAXIMUM) {
									std::cout << "overflow error: fabs(Amat.aij["<< is0 <<"]) > RealMAXIMUM !" << std::endl;
									system("pause");
								}
								ap_coarse[C_numerate[i8]] = Amat.abs_aij[is0];
								//std::cout << "find = "<< fabs(Amat.aij[is0])<<std::endl;
							}
						}
					}
					//std::cout << std::endl;
					//system("PAUSE");
				}
				
				//std::cout << "incoming="<< my_amg_manager.number_interpolation_procedure <<std::endl;
				
				//system("PAUSE");

				// Было уже рассчитано и запомнено выше по коду.
				//integer is_1 = row_startA[1];
				//integer is_e = row_startA[n_a[ilevel - 1] + 1] - 1;
				// Заранее один раз вычисляем модуль элемента.
				//for (integer iscan = is_1; iscan <= is_e; ++iscan) {
					//Amat.abs_aij[iscan] = fabs(Amat.aij[iscan]);
				//}

				// Интерполяция PMIS
				// 1.04.2017; 28.04.2017;
				// Главная идея в том чтобы разделить интерполяцию по знакам, отдельно положительные коэффициенты и отдельно положительные,
				// в итоге учитывается и то и то.
					

				if (0&&(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) &&
					(ilevel == 3)) {

					// 18.09.2021
					// Для высокой скорости сходимости при интерполяции обязательно должны добавляться новые С узлы = baggressive_interpol = false;
					// Если необходимое количество С узлов не будет добавлено, то сходимость утрачивается. Поэтому ни в коем случае не задавать baggressive_interpol = true;
					// Регулировать operator complexity можно параметрами truncation of interpolation, F-2-F, threshold.
					// long range interpolation значительно более сильная и сокращает число итераций с 15 до 9. Увеличивая при этом Operator complexity до 3.42 вместо 2.47. 
					// Smoothed aggregation дает operator complexity 2.1. По видимому Operator  complexity 2.45 для PMIS наилучшая = наименьшая при которой сходимость отличная.

					bool baggressive_interpol = false; // Никакого разрежения, иначе сходимость деградирует сильнейшим образом.

					// Дальнобойная процедура интерполяции на
						// основе интерполяционной процедуры №3.
						// Улучшенный базовый вариант.
						// Присоединение 20.09.2020.
						// Интерполяция с расстоянием 3 (дальнобойная).
					doublerealT theta_loc = static_cast<doublerealT>(theta(ilevel));
					doublerealT theta83_loc = static_cast<doublerealT>(theta(ilevel));//theta83
					my_long_range3_interpolation_procedure_number3A_parallelx8(
						the_number_of_neighbors_that_are_not_C_nodes,
						number_of_F_nodes_with_one_single_strong_C_neighbor,
						n_a, this_is_F_node, row_startA,
						nnz_a, bpositive_connections, Amat,
						bweSholdbeContinue, this_is_C_node, iadditionalCstatistic,
						RealZERO, icount1, P, nsizePR, ilevel,
						iadd, theta_loc, n, C_numerate,
						number_of_F_nodes_with_one_single_strong_C_neighborF,
						theta83_loc, btreshold_on_new_vetv, ifrom_re_operation_protection,
						from_re_operation_protection0, magic82, threshold_quick_all,
						threshold_quick_only_negative, bsuffix_work, baggressive_interpol);

					bsuffix_work = false; // запускаем только один раз суффиксную работу.

				}
				else {

					// 18.09.2021
					// Для высокой скорости сходимости при интерполяции обязательно должны добавляться новые С узлы = baggressive_interpol = false;
					// Если необходимое количество С узлов не будет добавлено, то сходимость утрачивается. Поэтому ни в коем случае не задавать baggressive_interpol = true;
					// Регулировать operator complexity можно параметрами truncation of interpolation, F-2-F, threshold.
					// long range interpolation значительно более сильная и сокращает число итераций с 15 до 9. Увеличивая при этом Operator complexity до 3.42 вместо 2.47. 
					// Smoothed aggregation дает operator complexity 2.1. По видимому Operator  complexity 2.45 для PMIS наилучшая = наименьшая при которой сходимость отличная.


					if (my_amg_manager.number_interpolation_procedure == 4) {

						bool baggressive_interpol = false; // Никакого разрежения, иначе сходимость деградирует сильнейшим образом.

						// Дальнобойная процедура интерполяции на
						// основе змееевидной интерполляции.
						// Улучшенный базовый вариант.
						// Присоединение 20.10.2021.
						// Интерполяция с расстоянием 6 (дальнобойная).
						// Присутствует разумно небольшое добавление узлов без которого сходимость деградирует.
						// Расстояния 4 явно недостаточно по силе сходимости V циклов в случае PMIS2_full агломаратора.
						// Расстояние 5 напротив отлично работает имея такую же операторную сложность.
						// Расстояние 6 Future 20.10.2021.
						doublerealT theta_loc = static_cast<doublerealT>(theta(ilevel));
						doublerealT theta83_loc = static_cast<doublerealT>(theta(ilevel));//theta83
						my_long_rang_interpolation_procedure_FFFFFFC_number3A_parallelx8(
							the_number_of_neighbors_that_are_not_C_nodes,
							number_of_F_nodes_with_one_single_strong_C_neighbor,
							n_a, this_is_F_node, row_startA,
							nnz_a, bpositive_connections, Amat,
							bweSholdbeContinue, this_is_C_node, iadditionalCstatistic,
							RealZERO, icount1, P, nsizePR, ilevel,
							iadd, theta_loc, n, C_numerate,
							number_of_F_nodes_with_one_single_strong_C_neighborF,
							theta83_loc, btreshold_on_new_vetv, ifrom_re_operation_protection,
							from_re_operation_protection0, magic82, threshold_quick_all,
							threshold_quick_only_negative, bsuffix_work, baggressive_interpol);

						bsuffix_work = false; // запускаем только один раз суффиксную работу.
					}
					else
					if (my_amg_manager.number_interpolation_procedure == 3) {

						bool baggressive_interpol = false; // Никакого разрежения, иначе сходимость деградирует сильнейшим образом.

						// Дальнобойная процедура интерполяции на
						// основе змееевидной интерполляции.
						// Улучшенный базовый вариант.
						// Присоединение 3.10.2021.
						// Интерполяция с расстоянием 4 (дальнобойная).
						// Присутствует разумно небольшое добавление узлов без которого сходимость деградирует.
						// Расстояния 4 явно недостаточно по силе сходимости V циклов в случае PMIS2_full агломаратора.
						// Расстояние 5 напротив отлично работает имея такую же операторную сложность.
						doublerealT theta_loc = static_cast<doublerealT>(theta(ilevel));
						doublerealT theta83_loc = static_cast<doublerealT>(theta(ilevel));//theta83
						my_long_rang_interpolation_procedure_FFFC4_number3A_parallelx8(
							the_number_of_neighbors_that_are_not_C_nodes,
							number_of_F_nodes_with_one_single_strong_C_neighbor,
							n_a, this_is_F_node, row_startA,
							nnz_a, bpositive_connections, Amat,
							bweSholdbeContinue, this_is_C_node, iadditionalCstatistic,
							RealZERO, icount1, P, nsizePR, ilevel,
							iadd, theta_loc, n, C_numerate,
							number_of_F_nodes_with_one_single_strong_C_neighborF,
							theta83_loc, btreshold_on_new_vetv, ifrom_re_operation_protection,
							from_re_operation_protection0, magic82, threshold_quick_all,
							threshold_quick_only_negative, bsuffix_work, baggressive_interpol);

						bsuffix_work = false; // запускаем только один раз суффиксную работу.
					}
					else
					if (my_amg_manager.number_interpolation_procedure == 0) {

						bool baggressive_interpol = false; // Никакого разрежения, иначе сходимость деградирует сильнейшим образом.

						// Дальнобойная процедура интерполяции на
						// основе змееевидной интерполляции.
						// Улучшенный базовый вариант.
						// Присоединение 25-26.09.2021.
						// Интерполяция с расстоянием 5 (дальнобойная).
						// Присутствует разумно небольшое добавление узлов без которого сходимость деградирует.
						doublerealT theta_loc = static_cast<doublerealT>(theta(ilevel));
						doublerealT theta83_loc = static_cast<doublerealT>(theta(ilevel));//theta83
						my_long_rang_interpolation_procedure_FFFC_number3A_parallelx8(
							the_number_of_neighbors_that_are_not_C_nodes,
							number_of_F_nodes_with_one_single_strong_C_neighbor,
							n_a, this_is_F_node, row_startA,
							nnz_a, bpositive_connections, Amat,
							bweSholdbeContinue, this_is_C_node, iadditionalCstatistic,
							RealZERO, icount1, P, nsizePR, ilevel,
							iadd, theta_loc, n, C_numerate,
							number_of_F_nodes_with_one_single_strong_C_neighborF,
							theta83_loc, btreshold_on_new_vetv, ifrom_re_operation_protection,
							from_re_operation_protection0, magic82, threshold_quick_all,
							threshold_quick_only_negative, bsuffix_work, baggressive_interpol);

						bsuffix_work = false; // запускаем только один раз суффиксную работу.
					}
					else
					if (my_amg_manager.number_interpolation_procedure == 2) {

						bool baggressive_interpol = false; // Никакого разрежения, иначе сходимость деградирует сильнейшим образом.

						// Дальнобойная процедура интерполяции на
						// основе интерполяционной процедуры №3.
						// Улучшенный базовый вариант.
						// Присоединение 20.09.2020.
						// Интерполяция с расстоянием 3 (дальнобойная).
						doublerealT theta_loc = static_cast<doublerealT>(theta(ilevel));
						doublerealT theta83_loc = static_cast<doublerealT>(theta(ilevel));//theta83
						my_long_range3_interpolation_procedure_number3A_parallelx8(
							the_number_of_neighbors_that_are_not_C_nodes,
							number_of_F_nodes_with_one_single_strong_C_neighbor,
							n_a, this_is_F_node, row_startA,
							nnz_a, bpositive_connections, Amat,
							bweSholdbeContinue, this_is_C_node, iadditionalCstatistic,
							RealZERO, icount1, P, nsizePR, ilevel,
							iadd, theta_loc, n, C_numerate,
							number_of_F_nodes_with_one_single_strong_C_neighborF,
							theta83_loc, btreshold_on_new_vetv, ifrom_re_operation_protection,
							from_re_operation_protection0, magic82, threshold_quick_all,
							threshold_quick_only_negative, bsuffix_work, baggressive_interpol);

						bsuffix_work = false; // запускаем только один раз суффиксную работу.
					}
					else {

						// 18.09.2021
					    // Для высокой скорости сходимости при интерполяции обязательно должны добавляться новые С узлы = baggressive_interpol = false;
					    // Если необходимое количество С узлов не будет добавлено, то сходимость утрачивается. Поэтому ни в коем случае не задавать baggressive_interpol = true;
					    // Регулировать operator complexity можно параметрами truncation of interpolation, F-2-F, threshold.
					    // long range interpolation значительно более сильная и сокращает число итераций с 15 до 9. Увеличивая при этом Operator complexity до 3.42 вместо 2.47. 
					    // Smoothed aggregation дает operator complexity 2.1. По видимому Operator  complexity 2.45 для PMIS наилучшая = наименьшая при которой сходимость отличная.

						bool baggressive_interpol = false; // Никакого разрежения, иначе сходимость деградирует сильнейшим образом.

						// Интерполяционная процедура №3.
						// Улучшенный базовый вариант.
						doublerealT theta_loc = static_cast<doublerealT>(theta(ilevel));
						doublerealT theta83_loc = static_cast<doublerealT>(theta(ilevel));//theta83
						// 31,01,2020. Интерполяция при которой не добавляется ни одного нового С узла, приводит к тому
						// что у многих F нет strong C соседа, хотя бы одного. Это приводит к сильнейшей деградации 
						// сходимости. Обязательно нужно добавить C узлы чтобы у каждого F был хотя бы один Strong C узел.
						my_interpolation_procedure_number3A_PMIS_parallel8(the_number_of_neighbors_that_are_not_C_nodes,
							number_of_F_nodes_with_one_single_strong_C_neighbor,
							n_a, this_is_F_node, row_startA,
							nnz_a, bpositive_connections, Amat,
							bweSholdbeContinue, this_is_C_node, iadditionalCstatistic,
							RealZERO, icount1, P, nsizePR, ilevel,
							iadd, theta_loc, n, C_numerate,
							number_of_F_nodes_with_one_single_strong_C_neighborF,
							theta83_loc, btreshold_on_new_vetv, ifrom_re_operation_protection,
							from_re_operation_protection0, magic82, threshold_quick_all,
							threshold_quick_only_negative, baggressive_interpol);
						/*
						// Здесь происходит добавление необходимых С узлов и сходимость многократно улучшается.
						my_interpolation_procedure_number3A(the_number_of_neighbors_that_are_not_C_nodes,
							number_of_F_nodes_with_one_single_strong_C_neighbor,
							n_a, this_is_F_node, row_startA,
							nnz_a, bpositive_connections, Amat,
							bweSholdbeContinue, this_is_C_node, iadditionalCstatistic,
							RealZERO, icount1, P, nsizePR, ilevel,
							iadd, theta_loc, n, C_numerate,
							number_of_F_nodes_with_one_single_strong_C_neighborF,
							theta83_loc, btreshold_on_new_vetv, ifrom_re_operation_protection,
							from_re_operation_protection0, magic82, threshold_quick_all,
							threshold_quick_only_negative);
							*/
					}
				}

				if (bprint_mesage_diagnostic) {
					std::cout << "Additional C nodes in interpolation procedure. Statistics:\n";
					std::cout << "number firstable C nodes=" << n_coarce << ", number secondary C nodes="<< iadditionalCstatistic <<"\n";
				}

				if (bweSholdbeContinue) {
					//delete[] ap_coarse;
					if (ap_coarse != nullptr) {
						free(ap_coarse);
						ap_coarse = nullptr;
					}
					if (bprint_mesage_diagnostic) {
						std::cout << "Feedback restart...\n";
					}
					bsuffix_work = false; // recomended false for Jacoby interpolation.
				}

				if (bprint_mesage_diagnostic) {
					// отношение добавленных узлов к количеству С узлов на предыдущем уровне.
					//std::cout << "addition C nodes "<< static_cast<doublerealT>(100.0*iadditionalCstatistic / n_a[ilevel - 1]) <<"%\n";
					// отношение количества добавленных С узлов к первоначальному количеству С узлов на данном уровне.
					std::cout << "addition C nodes = " << static_cast<doublerealT>(100.0 * iadditionalCstatistic / n_coarce) <<"% firstable C nodes,  level population " << static_cast<doublerealT>(100.0 * (n_coarce + iadditionalCstatistic) / n_a[ilevel - 1])<<"%\n";
				}
				iadditionalCstatistic = 0;
				//system("pause");

			}

        }
		nnzR = icount1 - iaddR;
		if (bprint_mesage_diagnostic) {
			std::cout << "Prolongation operator complexity = "<< static_cast<doublerealT>(1.0*icount1 / n) << "*n\n";
		}
		//system("pause");

		// нужно определить nnzR количество ненулевых элементов в матрице R и P.

		// оператор restriction построен и он упорядочен по i.
		// число ненулевых элементов nnzR-1.
		// P=Copy(R);
		iend_marker_position = iaddR + nnzR - 1;

		// truncation of interpolation.
		// 30.04.2017.
		if (my_amg_manager.itruncation_interpolation == 1) {

			/*
			// Однопоточный вариант работает и без сортировки,
			// что говорит о том что оператор интерполяции уже предварительно был отсортирован по j.
			switch (imy_sort_algorithm ) {
			case MY_SORT_ALGORITHM:: COUNTING_SORT:
			//Counting_Sortj(P, 1 + iaddR, iaddR + nnzR - 1, false);
			qsj(P, 1 + iaddR, iaddR + nnzR - 1);
			//HeapSort(P, 1 + iaddR, iaddR + nnzR - 1, comparej);
			break;
			case MY_SORT_ALGORITHM::QUICK_SORT:
			qsj(P, 1 + iaddR, iaddR + nnzR - 1);
			// Библиотечный алгоритм. O(nlog(n)).
			// Не использует лишней памяти.
			//std::sort(P + (1 + iaddR) * sizeof(Ak1), P + (iaddR + nnzR - 1+1) * sizeof(Ak1), compAi);
			break;
			case MY_SORT_ALGORITHM::HEAP_SORT:
			HeapSort(P, 1 + iaddR, iaddR + nnzR - 1, comparej);
			break;
			case MY_SORT_ALGORITHM::TIM_SORT:
				// Сортировка Тима Петерсона.
				timSort_amg_j(P, 1 + iaddR, iaddR + nnzR - 1);
				break;
			default:
			//Counting_Sortj(P, 1 + iaddR, iaddR + nnzR - 1, false);
			qsj(P, 1 + iaddR, iaddR + nnzR - 1);
			//HeapSort(P, 1 + iaddR, iaddR + nnzR - 1, comparej);
			break;
			}
			*/


			// Большое число связей увеличивает сложность оператора Галёркина.
			// Наличие слабых связей в процедуре интерполяции, приводит к замедлению 
			// сходимости или расходимости.
			// Алгоритм устранения слабых связей:
			//doublerealT const alpha_truncation = 0.2;
			doublerealT alpha_truncation = static_cast<doublerealT>(my_amg_manager.truncation_interpolation);
			// Рассмотрим каждую строку оператора интерполяции.
			// Найдем сумму элементов данной строки каждого знака.
			// Найдём максимальный по модулю элемент каждого знака.
			// Удалим все элементы в операторе интерполяции каждого знака 
			// которые меньше максимального по модулю того-же знака * на alpha_truncation.
			// Проведём перемасштабирование чтобы сумма осталась неизменной.
			// Сделаем это в памяти P. 17.02.2019
#pragma omp parallel for
			for (integer i_1 = 1; i_1 <= n_a[ilevel-1]; ++i_1) {
				flag[i_1] = false; // init flag.
			}
			integer icounter_truncation = iend_marker_position + 1;//1 + iaddR;

			if (1) {
				// Многопоточная версия.

				bool bstrong_trunc = false;
				if (1&&((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS)))
				{
					// PMIS, HMIS
					bstrong_trunc = true; // не более 4 записей в каждой строке оператора интерполяции.
				}

				integer i_size_75 = 0;
				// Это нельзя распараллелить.
				for (integer ii = 1 + iaddR; ii <= iend_marker_position; ++ii) {
					if (flag[P[ii].j] == false) {
						//row_ind_SRloc[P[ii].j] = ii;
						flag[P[ii].j] = true;
						//i_size_75++;
						if (P[ii].j > i_size_75) i_size_75 = P[ii].j;
					}
				}

#pragma omp parallel for
				for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
					flag[i_1] = false; // init flag.
				}

				// inicialization обязательна.
				integer* row_ind_SRloc = my_declaration_array<integer>(i_size_75, -1, "row_ind_SRloc");

				//std::cout << "number of coarcenodes=" << numberofcoarcenodes << " i_size_75=" << i_size_75 <<"\n";

				//system("pause");

				// Это нельзя распараллелить.
				for (integer ii = 1 + iaddR; ii <= iend_marker_position; ++ii) {
					if (flag[P[ii].j] == false) {
						row_ind_SRloc[P[ii].j] = ii;
						flag[P[ii].j] = true;
					}
				}

				//for (integer ii = 1 + iaddR; ii <= iend_marker_position; ++ii) {
				// Это нельзя распараллелить. 06.07.2019
				// оператор интерполяции P заполняется строго последовательно.
				//#pragma omp parallel for
				for (integer i_75 = 1; i_75 <= i_size_75; ++i_75) {
					if (row_ind_SRloc[i_75] != -1) {
						integer ii = row_ind_SRloc[i_75];

						//if (flag[P[ii].j] == false) {
						//flag[P[ii].j] = true;
						integer istr_65 = P[ii].j;
						integer ii_65 = ii;
						doublerealT dsum_plus = 0.0;
						doublerealT dsum_minus = 0.0;
						doublerealT dmax_plus = -1.0;
						doublerealT dmax_minus = -1.0;
						while ((ii_65 <= iend_marker_position) && (P[ii_65].j == istr_65)) {
							if (P[ii_65].aij > 0) {
								dsum_plus += P[ii_65].aij;
								if (P[ii_65].aij > dmax_plus) dmax_plus = P[ii_65].aij;
							}
							if (P[ii_65].aij < 0) {
								dsum_minus += fabs(P[ii_65].aij);
								if (fabs(P[ii_65].aij) > dmax_minus) dmax_minus = fabs(P[ii_65].aij);
							}
							++ii_65;
						}
						ii_65 = ii;
						doublerealT dsum_plus_new = 0.0;
						doublerealT dsum_minus_new = 0.0;
						
						const integer IMAX_STRONG_TRUNC = 4;
						integer icount_strong_trunc = 0;

						while ((ii_65 <= iend_marker_position) && (P[ii_65].j == istr_65)) {
							if ((P[ii_65].aij > 0) && (fabs(P[ii_65].aij) > alpha_truncation*dmax_plus)) {
								if (bstrong_trunc) {
									if (icount_strong_trunc < IMAX_STRONG_TRUNC) {
										dsum_plus_new += fabs(P[ii_65].aij);
									}
									icount_strong_trunc++;
								}
								else {
									dsum_plus_new += fabs(P[ii_65].aij);
								}
							}
							if ((P[ii_65].aij < 0) && (fabs(P[ii_65].aij) > alpha_truncation*dmax_minus)) {
								if (bstrong_trunc) {
									if (icount_strong_trunc < IMAX_STRONG_TRUNC) {
										dsum_minus_new += fabs(P[ii_65].aij);
									}
									icount_strong_trunc++;
								}
								else {
									dsum_minus_new += fabs(P[ii_65].aij);
								}
							}
							++ii_65;
						}
						// заполнение перемасштабированными.

						icount_strong_trunc = 0;
						ii_65 = ii;
						while ((ii_65 <= iend_marker_position) && (P[ii_65].j == istr_65)) {
							if ((P[ii_65].aij > 0) && (fabs(P[ii_65].aij) > alpha_truncation*dmax_plus)) {
								if (bstrong_trunc) {
									if (icount_strong_trunc < IMAX_STRONG_TRUNC) {
										P[icounter_truncation] = P[ii_65];
										P[icounter_truncation].aij = static_cast<real_mix_precision>(fabs(dsum_plus / dsum_plus_new)*P[ii_65].aij);
										icounter_truncation++;
									}
									icount_strong_trunc++;
								}
								else {
									P[icounter_truncation] = P[ii_65];
									P[icounter_truncation].aij = static_cast<real_mix_precision>(fabs(dsum_plus / dsum_plus_new)*P[ii_65].aij);
									icounter_truncation++;
								}
							}
							if ((P[ii_65].aij < 0) && (fabs(P[ii_65].aij) > alpha_truncation*dmax_minus)) {
								if (bstrong_trunc) {
									if (icount_strong_trunc < IMAX_STRONG_TRUNC) {
										P[icounter_truncation] = P[ii_65];
										P[icounter_truncation].aij = static_cast<real_mix_precision>(fabs(dsum_minus / dsum_minus_new)*P[ii_65].aij);
										icounter_truncation++;
									}
									icount_strong_trunc++;
								}
								else {
									P[icounter_truncation] = P[ii_65];
									P[icounter_truncation].aij = static_cast<real_mix_precision>(fabs(dsum_minus / dsum_minus_new)*P[ii_65].aij);
									icounter_truncation++;
								}
							}
							++ii_65;
						}

					}
				}

				free(row_ind_SRloc);
				row_ind_SRloc = nullptr;
			}
			else {

				// Однопоточная версия.

				for (integer ii = 1 + iaddR; ii <= iend_marker_position; ++ii) {
					if (flag[P[ii].j] == false) {
						flag[P[ii].j] = true;
						integer istr_65 = P[ii].j;
						integer ii_65 = ii;
						doublerealT dsum_plus = 0.0;
						doublerealT dsum_minus = 0.0;
						doublerealT dmax_plus = -1.0;
						doublerealT dmax_minus = -1.0;
						while ((ii_65 <= iend_marker_position) && (P[ii_65].j == istr_65)) {
							if (P[ii_65].aij > 0) {
								dsum_plus += P[ii_65].aij;
								if (P[ii_65].aij > dmax_plus) dmax_plus = P[ii_65].aij;
							}
							if (P[ii_65].aij < 0) {
								dsum_minus += fabs(P[ii_65].aij);
								if (fabs(P[ii_65].aij) > dmax_minus) dmax_minus = fabs(P[ii_65].aij);
							}
							++ii_65;
						}
						ii_65 = ii;
						doublerealT dsum_plus_new = 0.0;
						doublerealT dsum_minus_new = 0.0;
						
						while ((ii_65 <= iend_marker_position) && (P[ii_65].j == istr_65)) {
							if ((P[ii_65].aij > 0) && (fabs(P[ii_65].aij) > alpha_truncation*dmax_plus)) {
								dsum_plus_new += fabs(P[ii_65].aij);
							}
							if ((P[ii_65].aij < 0) && (fabs(P[ii_65].aij) > alpha_truncation*dmax_minus)) {
								dsum_minus_new += fabs(P[ii_65].aij);
							}
							++ii_65;
						}
						// заполнение перемасштабированными.
						ii_65 = ii;
						while ((ii_65 <= iend_marker_position) && (P[ii_65].j == istr_65)) {
							if ((P[ii_65].aij > 0) && (fabs(P[ii_65].aij) > alpha_truncation*dmax_plus)) {
								P[icounter_truncation] = P[ii_65];
								P[icounter_truncation].aij = static_cast<real_mix_precision>(fabs(dsum_plus / dsum_plus_new)*P[ii_65].aij);
								icounter_truncation++;
							}
							if ((P[ii_65].aij < 0) && (fabs(P[ii_65].aij) > alpha_truncation*dmax_minus)) {
								P[icounter_truncation] = P[ii_65];
								P[icounter_truncation].aij = static_cast<real_mix_precision>(fabs(dsum_minus / dsum_minus_new)*P[ii_65].aij);
								icounter_truncation++;
							}
							++ii_65;
						}

					}
				}
			}

			// Ужатие (обратное копирование).
			integer ist_in_P = 1 + iaddR;
			// Мы дописывали новые коэффициенты в конец матрицы интерполяции P.
			// Это нельзя распараллеливать так явно.
//#pragma omp parallel for
			for (integer ii = iend_marker_position+1; ii <= icounter_truncation-1; ++ii) {
				P[ist_in_P++] = P[ii];
			}
			iend_marker_position = ist_in_P - 1;

			//iend_marker_position = iaddR + nnzR - 1;
			nnzR = iend_marker_position - iaddR + 1;
			//nnzR = icount1 - iaddR;
			icount1 = nnzR + iaddR;

#pragma omp parallel for
			for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; ++i_1) {
				flag[i_1] = false; // init flag.
			}
		}

		

		// Этот оператор нужен для вычисления grid complexity для оператора 
		// интерполяции и проекции. Данная информация важна для оптимизации количества выделяемой памяти.
		if (ilevel - 1 == 0) {
			nnz_P_memo_0 = iend_marker_position - (iaddR + 1) + 1;
		}
		nnz_P_memo_all = iend_marker_position;
		

		// где то надо разделить на ap, т.к. 
		// R=P/ap. ????  
		// НЕТ делить НЕ НАДО!!! т.к. в теории R=transpose(P).

		// Сортировка оператора интерполяции P по строкам 17.02.2018
		// heapsort(P,key==i,iaddR+1,iaddR+nnzR - 1);

		integer isize20 = nnzR / 2;
		int inth0 = number_cores();
		if (inth0 > 1) {
			omp_set_num_threads(inth0);
		}
		//omp_set_num_threads(8);

		switch (imy_sort_algorithm) {
		case MY_SORT_ALGORITHM:: COUNTING_SORT:
			//Counting_Sort(P, 1 + iaddR, iaddR + nnzR - 1, false, numberofcoarcenodes, indx_comparei);// numberofcoarcenodes <-> n_a[ilevel - 1]
			if (inth0 == 1) {
#pragma omp parallel for
				for (integer i = 0; i <= numberofcoarcenodes; ++i) {
					C1[i] = 0; // инициализация.
					//C2[i] = 0; // инициализация.
				}
				if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) || 
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
				{
					Counting_Sort_bmemo_false(P, 1 + iaddR, iaddR + nnzR - 1, numberofcoarcenodes, indx_comparei, C1, Bm1);// numberofcoarcenodes <-> n_a[ilevel - 1]
				}
				else {
					Counting_Sort_bmemo_false(P, 1 + iaddR, iaddR + nnzR - 1, numberofcoarcenodes, indx_comparei, C1);// numberofcoarcenodes <-> n_a[ilevel - 1]
				}
			}
			else {
#pragma omp parallel for
				for (int i = 0; i <= numberofcoarcenodes; ++i) {
					C1[i] = 0; // инициализация.
					C2[i] = 0; // инициализация.
				}

#pragma omp parallel sections
				{
#pragma omp section
					{
						if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
						{
							Counting_Sort_bmemo_false(P, 1 + iaddR, isize20 + iaddR, numberofcoarcenodes, indx_comparei, C1,Bm1);
						}
						else {
							Counting_Sort_bmemo_false(P, 1 + iaddR, isize20 + iaddR, numberofcoarcenodes, indx_comparei, C1);
						}
					}
#pragma omp section
					{
						if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) || 
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
						{
							Counting_Sort_bmemo_false(P, isize20 + 1 + iaddR, iaddR + nnzR - 1, numberofcoarcenodes, indx_comparei, C2,Bm2);
						}
						else {
							Counting_Sort_bmemo_false(P, isize20 + 1 + iaddR, iaddR + nnzR - 1, numberofcoarcenodes, indx_comparei, C2);
						}
					}
				}
				mergeTim_amg(P, 1 + iaddR, isize20 + iaddR, iaddR + nnzR - 1, indx_comparei);
			}
			break;
		case MY_SORT_ALGORITHM::QUICK_SORT:
			//qs_abbys_heigh = 0;
			//qs(P, 1 + iaddR, iaddR + nnzR - 1);
			// Библиотечный алгоритм. O(nlog(n)).
			// Не использует лишней памяти.
			if (inth0 == 1) {
				std::sort(P + 1 + iaddR, P + iaddR + nnzR - 1 + 1, compareAk1R);
			}
		else {
#pragma omp parallel sections
			{
#pragma omp section
				{
					std::sort(P + 1 + iaddR, P + isize20 + iaddR, compareAk1R);
				}
#pragma omp section
				{
					std::sort(P + isize20 + 1 + iaddR, P + iaddR + nnzR - 1 + 1, compareAk1R);
				}
			}
			mergeTim_amg(P, 1 + iaddR, isize20 + iaddR, iaddR + nnzR - 1, indx_comparei);

		}

			/*
#pragma omp parallel sections
		{
#pragma omp section
			{
				qs(Amat, 1 + iadd, isize2 + iadd);
			}
#pragma omp section
			{
				qs(Amat, isize2 + 1 + iadd, nnz_a[ilevel - 1] + iadd);
			}
		}
		mergeTim_amg(Amat, 1 + iadd, isize2 + iadd, nnz_a[ilevel - 1] + iadd);
		*/

			break;
		case MY_SORT_ALGORITHM::HEAP_SORT:
			if (inth0 == 1) {
				//HeapSort(P, 1 + iaddR, iaddR + nnzR - 1,comparei);
				//LeftistHeapSort(P, 1 + iaddR, iaddR + nnzR - 1);
				mySTDHeapSort(P, 1 + iaddR, iaddR + nnzR - 1, indx_comparei);
			}
			else {

#pragma omp parallel sections
				{
#pragma omp section
					{
						mySTDHeapSort(P, 1 + iaddR, isize20 + iaddR, indx_comparei);
					}
#pragma omp section
					{
						mySTDHeapSort(P,  isize20 + 1 + iaddR,  iaddR + nnzR - 1, indx_comparei);
					}
				}
				mergeTim_amg(P, 1 + iaddR, isize20 + iaddR, iaddR + nnzR - 1, indx_comparei);
			}
			break;
		case MY_SORT_ALGORITHM::TIM_SORT:
			// Сортировка Тима Петерсона 2002г.
			//timSort_amg(P, 1 + iaddR, iaddR + nnzR - 1, indx_comparei);
			if (inth0 == 1) {
				gfx::timsort(P + 1 + iaddR, P + iaddR + nnzR - 1 + 1, compareAk1R);
			}
			else {
				std::cout << "number threads=" <<  omp_get_max_threads() << std::endl;

				

#pragma omp parallel sections
				{
#pragma omp section
					{
						gfx::timsort(P + 1 + iaddR, P + isize20 + iaddR, compareAk1R);
					}
#pragma omp section
					{
						gfx::timsort(P + isize20 + 1 + iaddR, P + iaddR + nnzR - 1 + 1, compareAk1R);
					}
				}
				mergeTim_amg(P, 1 + iaddR, isize20 + iaddR, iaddR + nnzR - 1, indx_comparei);
			}
			break;
		default:
			Counting_Sort(P, 1 + iaddR, iaddR + nnzR - 1, false, numberofcoarcenodes, indx_comparei);// numberofcoarcenodes <-> n_a[ilevel - 1]
			break;
		}
		
		//omp_set_num_threads(1);


		if (bprint_mesage_diagnostic) {

			std::cout << "first level size n="<< n <<"; number of coarcenodes="<< numberofcoarcenodes <<", procent = "<<100.0*numberofcoarcenodes/n<<"%\n";
		}

		// Проверка Restriction нет ли пропусков строк при интерполяции: 
		// Роль R играет P сортированное по строкам (транспонированное).
		if (1) {// Проверка обязательна. Она страхует от ошибок в новых версиях
			// алгебраического многосеточного алгоритма.
#pragma omp parallel for
			for (integer i_1 = 1; i_1 <= numberofcoarcenodes; ++i_1) {
				flag[i_1] = false; // init flag.
			}
			for (integer i_1 = 1 + iaddR; i_1 <= iaddR + nnzR - 1; ++i_1) {
				if (flag[P[i_1].i] == false)
				{
					doublerealT dsum27 = 0.0;
					for (integer i_2 = i_1; (i_2 <= iaddR + nnzR - 1) && (P[i_2].i == P[i_1].i); ++i_2) {
						dsum27 += fabs(P[i_2].aij);
					}
					if (dsum27 < 1.0e-37) {
						std::cout << "fatal error!!! zero string R["<< P[i_1].i <<"][j]="<< dsum27 <<std::endl;

						system("PAUSE");
					}
					flag[P[i_1].i] = true;
				}
			}
			for (integer i_1 = 1; i_1 <= numberofcoarcenodes; ++i_1) {
				if (flag[i_1] == false) {
					//06.07.2019
					// пропуск строки номер i_1
					std::cout << "fatal error!!! string number " << i_1 << " propushena\n";

					system("PAUSE");
				}
			}
#pragma omp parallel for
			for (integer i_1 = 1; i_1 <= numberofcoarcenodes; ++i_1) {
				flag[i_1] = false; // init flag.
			}
		}

		
		if (ap_coarse != nullptr) {
			free(ap_coarse);
			ap_coarse = nullptr;
		}
	
		// Освобождение оперативной памяти из под хеш-таблицы.02.02.2019
		free_hash_table_Gus_struct01();

		if (bprint_mesage_diagnostic) {
			std::cout << "Prolongation ierarhion..." << std::endl;
		}
		if (b_REALLOC) {
			my_realloc_memory<Ak1>(P, (static_cast<integer>(nnz_P_memo_all)+2));
			
			if (bprint_mesage_diagnostic) {
				std::cout << "2 of 3 compleated. OK!! ierarhion matrix Prolongation realloc successfully..." << std::endl;
			}
		}

		// MARKER GUSTAVSON

		// Нахождение матрицы грубосеточного уровня:
		// Acorse=R*Afine*P;
		// часть 1: R*Afine.
		//         xxxxxx
		//         xxxxxx
		//  xxxxxx xxxxxx xxxxxx
		//  xxxxxx xxxxxx xxxxxx
		//         xxxxxx
		//         xxxxxx
		//    R       Amat     [RA]


		if (bprint_mesage_diagnostic) {
			std::cout << "   ***   CAMG GALERKIN MULTIPLICATOR " << ilevel <<"  ***\n";

			std::cout << "nnz left operand=" << nnzR << ", nnz right operand="<< nnz_a[ilevel - 1] <<"\n";

		}

		integer istartAnew;


		// Фред Густавсон IBM 1978.
		// 23 октября 2015 года.

		// часть 1: R*Afine.
		//         xxxxxx
		//         xxxxxx
		//  xxxxxx xxxxxx xxxxxx
		//  xxxxxx xxxxxx xxxxxx
		//         xxxxxx
		//         xxxxxx
		//    R       Amat     [RA]
		// Сортировка А по строкам.
		/*
		// 7 января 2016. Сортировка матрицы А была выполнена единожды при начале обработке данного уровня.
		switch (imy_sort_algorithm) {
		case MY_SORT_ALGORITHM:: COUNTING_SORT:
		Counting_Sort(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd);
		break;
		case MY_SORT_ALGORITHM::QUICK_SORT :
		qs(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd);
		// Библиотечный алгоритм. O(nlog(n)).
		// Не использует лишней памяти.
		//std::sort(Amat + (1 + iadd)*sizeof(Ak1), Amat + (nnz_a[ilevel - 1] + iadd)*sizeof(Ak1), compAi);
		break;
		case MY_SORT_ALGORITHM::HEAP_SORT:
		HeapSort(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd, comparei);
		break;
		case MY_SORT_ALGORITHM::TIM_SORT:
			// Сортировка Тима Петерсона.
			timSort_amg(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd);
			break;
		default:
		Counting_Sort(Amat, 1 + iadd, nnz_a[ilevel - 1] + iadd);
		break;
		}
		*/
		// Преобразование к формату CRS.

		if (row_ind_SR != nullptr) {
			free(row_ind_SR);
			row_ind_SR = nullptr;
		}
		row_ind_SR = my_declaration_array<integer>(numberofcoarcenodes, -1, "row_ind_SR");
		if (row_ind_ER != nullptr) {
			free(row_ind_ER);
			row_ind_ER = nullptr;
		}
		row_ind_ER = my_declaration_array<integer>(numberofcoarcenodes, -2, "row_ind_ER");

		
		
#pragma omp parallel for
		for (integer i = 1; i <= icounter - 1; ++i) {
			flag[i] = false;
		}

		integer istart1 = 1 + iaddR;
		integer iend1 = nnzR - 1 + iaddR;
		// Роль R играет транспонированный (сортированный по строкам) оператор интерполяции.
		calculate_row_ptr(istart1, iend1, row_ind_SR, row_ind_ER, flag, P);
		
			if (bQuick_sort_for_reorder) {
				for (integer ii = 1; ii <= numberofcoarcenodes; ++ii) {
					// 14.04.2020
					// В каждой строке i индексы столбца j отсортированы по возрастанию.		
					//qs_abbys_heigh = 0;
					// quicksort
					//qsj(P, row_ind_SR[ii], row_ind_ER[ii]);
					std::sort(P + row_ind_SR[ii], P + row_ind_ER[ii] + 1, compareAk1P);
				}
			}

#pragma omp parallel for
		for (integer i = 1; i <= n_a[ilevel - 1]; ++i) {
			flag[i] = false;
		}			

		istartAnew = nnz_a[ilevel - 1] + 1 + iadd;
		//integer istartAnew_mem = istartAnew;

		

		

		if ((number_cores() == 1)) {

#pragma omp parallel for
			for (integer i = 1; i <= n_a[ilevel - 1]; ++i) {
				vector_sum[i] = 0.0;
				hash_table[i] = false;
				index_visit[i] = 0;
			}


			// Умножение R*A. Результат пишется в А начиная с позиции istartAnew.
			// P - это R, Amat - это А.
			my_sparse_matrix_by_matrix_multiplication_RA<doublerealT>(Amat, P, istartAnew,
				row_ind_SR, row_ind_ER, row_startA, numberofcoarcenodes,
				// Вспомогательные данные.
				hash_table, index_visit, vector_sum, nsizeA);
		}
		else {

			//#ifdef _NONAME_STUB29_10_2017
#ifdef _OPENMP

	//	for (int i_9 = 0; i_9 < iKnumber_thread; ++i_9) {
//#pragma omp parallel for
			//for (integer i_91 = 0; i_91 < n + 1; ++i_91) {
				//vector_sum_m[i_9][i_91] = 0.0;
			//}
		//}

		// Параллельное умножение R*A. Результат пишется в А начиная с позиции istartAnew.
		// P - это R, Amat - это А.
			my_parallel8_sparse_matrix_by_matrix_multiplication_RA<doublerealT>(Amat,
				P, istartAnew, istartAnew_m,
				row_ind_SR, row_ind_ER, row_startA,
				numberofcoarcenodes, iKnumber_thread,
				hash_table_m, index_visit_m,
				vector_sum_m, index_size_m,
				nsizeA, n_a[ilevel - 1], nnz, ilevel,
				bprint_mesage_diagnostic, n_a, AccumulqtorA_m_SIZE8,
				AccumulqtorA_m);

#else

		// Данные используемые для частичного формирователя суммы.
		/*if (vector_sum != nullptr) {
			free(vector_sum);
			vector_sum = nullptr;
		}
		vector_sum = my_declaration_array<doublerealT>(n_a[ilevel - 1], 0.0, "vector_sum");*/


		/*
		// Храним индексы ненулевых элементов в отсортированном порядке.
		if (index_visit != nullptr) {
			free(index_visit);
			index_visit = nullptr;
		}
		index_visit = my_declaration_array<integer>(n_a[ilevel - 1], 0, "index_visit");


		// hash_table nnz+1
		// Огромного размера хеш-таблица.
		// Огромный размер поэтому инициализация делается лишь единожды.
		// размер от 0 до nnz включительно.
		bool* hash_table = my_declaration_array<bool>(n_a[ilevel - 1], false, "hash_table");
		*/

#pragma omp parallel for
			for (integer i = 1; i <= n_a[ilevel - 1]; ++i) {
				vector_sum[i] = 0.0;
				hash_table[i] = false;
				index_visit[i] = 0;
			}


			// Умножение R*A. Результат пишется в А начиная с позиции istartAnew.
			// P - это R, Amat - это А.
			my_sparse_matrix_by_matrix_multiplication_RA<doublerealT>(Amat, P, istartAnew,
				row_ind_SR, row_ind_ER, row_startA, numberofcoarcenodes,
				// Вспомогательные данные.
				hash_table, index_visit, vector_sum, nsizeA);

#endif

		}

		
		/*if (index_visit != nullptr) {
			free(index_visit);
			index_visit = nullptr;
		}*/
		if (row_ind_SR!=nullptr) {
			free(row_ind_SR);
			row_ind_SR = nullptr;
		}
		if (row_ind_ER != nullptr) {
			free(row_ind_ER);
			row_ind_ER = nullptr;
		}
		
		/*if (vector_sum != nullptr) {
			free(vector_sum);
			vector_sum = nullptr;
		}*/



		// Часть 2. [R*Afine]*P=Abuf*P.
		// Сортировка [R*А] по i.
		//heapsort(Amat, key=i*n_coarce + j, 1, nnz_a[ilevel - 1]);

		// В результате работы алгоритма разреженного матричного умножения по Ф. Густавсону,
		// мы и так имеем отсортированный по строкам результат, поэтому дополнительная 
		// сортировка не требуется. Это проверено 11 января 2016.
		// 11 января 2016.

		


		// Prolongation должна быть упорядочена по j.
		// Начальная позиция элементов матрицы грубосеточного уровня.
		integer istartAnew2 = istartAnew;



		// Быстрее этого кода на основе идеи слияния списков уже не будет.
		// 17 октября 2015. Нужно двигаться в сторону Писсанецки.
		if (bprint_mesage_diagnostic) {

			std::cout << "nnz left operand=" << istartAnew - (nnz_a[ilevel - 1] + 1 + iadd) << ", nnz right operand="<< nnzR <<"\n";
		}
		
		// Фред Густавсон IBM 1978
		// В ядре кода Густавсона нету ни одного ветвления,
		// а мы знаем что в результате профайлинга предыдущих версий кода:
		// (наивный, слияние, Писсанецки) львиная доля вычислительной работы уходила
		// на сравнения (ветвления) в отношении примерно 30 к 1. 30 сравнений на одно суммирование.
		// 23 октября 2015 года.
		// 6 января 2016 года Добавлено АВЛ дерево.


		// Рабочая версия алгоритма Фреда Густавсона.
		// IBM 1978 Sparse Matrix multiplication.

		integer isize21 = nnzR / 2;

		int inth = number_cores();
		if (inth > 1) {
				omp_set_num_threads(inth);
		}
		//omp_set_num_threads(8);

		// Сортировка обязательно требуется.
		// Преобразование обоих матриц в формат CRS.
		// Сортировка матрицы интерполяции по столбцам.
		switch (imy_sort_algorithm) {
		case MY_SORT_ALGORITHM::COUNTING_SORT:
			//Counting_Sort(P, 1 + iaddR, iaddR + nnzR - 1, false, n_a[ilevel - 1], indx_comparej);//подходит именно n_a[ilevel - 1]
			if (inth == 1) {
#pragma omp parallel for
				for (integer i = 0; i <= n_a[ilevel - 1]; ++i) {
					C1[i] = 0; // инициализация.
					//C2[i] = 0; // инициализация.
				}
				if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
					(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
				{
					Counting_Sort_bmemo_false(P, 1 + iaddR, iaddR + nnzR - 1, n_a[ilevel - 1], indx_comparej, C1, Bm1);//подходит именно n_a[ilevel - 1]
				}
				else {
					Counting_Sort_bmemo_false(P, 1 + iaddR, iaddR + nnzR - 1, n_a[ilevel - 1], indx_comparej, C1);//подходит именно n_a[ilevel - 1]
				}
			}
			else {

#pragma omp parallel for
				for (integer i = 0; i <= n_a[ilevel - 1]; ++i) {
					C1[i] = 0; // инициализация.
					C2[i] = 0; // инициализация.
				}

#pragma omp parallel sections
				{
#pragma omp section
					{
						if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
						{
							Counting_Sort_bmemo_false(P, 1 + iaddR, isize21 + iaddR, n_a[ilevel - 1], indx_comparej, C1,Bm1);
						}
						else {
							Counting_Sort_bmemo_false(P, 1 + iaddR, isize21 + iaddR, n_a[ilevel - 1], indx_comparej, C1);
						}
					}
#pragma omp section
					{
						if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
							(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
						{
							Counting_Sort_bmemo_false(P, isize21 + 1 + iaddR, iaddR + nnzR - 1, n_a[ilevel - 1], indx_comparej, C2,Bm2);
						}
						else {
							Counting_Sort_bmemo_false(P, isize21 + 1 + iaddR, iaddR + nnzR - 1, n_a[ilevel - 1], indx_comparej, C2);
						}
					}
				}
				mergeTim_amg(P, 1 + iaddR, isize21 + iaddR, iaddR + nnzR - 1, indx_comparej);
			}
			break;
		case MY_SORT_ALGORITHM::QUICK_SORT:
			//qs_abbys_heigh = 0;
			//qsj(P, 1 + iaddR, iaddR + nnzR - 1);
			// Библиотечный алгоритм. O(nlog(n)).
			// Не использует лишней памяти.
			if (inth == 1) {
				std::sort(P + 1 + iaddR, P + iaddR + nnzR - 1 + 1, compareAk1P);
			}
			else {
#pragma omp parallel sections
				{
#pragma omp section
					{
						std::sort(P + 1 + iaddR, P+ isize21 + iaddR, compareAk1P);
					}
#pragma omp section
					{
						std::sort(P + isize21 + 1 + iaddR, P + iaddR + nnzR - 1 + 1, compareAk1P);
					}
				}
				mergeTim_amg(P, 1 + iaddR, isize21 + iaddR, iaddR + nnzR - 1, indx_comparej);

			}
			break;
		case MY_SORT_ALGORITHM::HEAP_SORT:
			//HeapSort(P, 1 + iaddR, iaddR + nnzR - 1, comparej);
			//LeftistHeapSort_j(P, 1 + iaddR, iaddR + nnzR - 1);
			if (inth == 1) {
				mySTDHeapSort(P, 1 + iaddR, iaddR + nnzR - 1, indx_comparej);
			}
			else {
#pragma omp parallel sections
				{
#pragma omp section
					{
						mySTDHeapSort(P, 1 + iaddR, isize21 + iaddR, indx_comparej);
					}
#pragma omp section
					{
						mySTDHeapSort(P, isize21 + 1 + iaddR, iaddR + nnzR - 1, indx_comparej);
					}
				}
				mergeTim_amg(P, 1 + iaddR, isize21 + iaddR, iaddR + nnzR - 1, indx_comparej);
			}
			break;
		case MY_SORT_ALGORITHM::TIM_SORT:
			// Сортировка Тима Петерсома 2002.
			//timSort_amg(P, 1 + iaddR, iaddR + nnzR - 1, indx_comparej);
			 
			if (inth == 1) {

				gfx::timsort(P + 1 + iaddR, P + iaddR + nnzR - 1 + 1, compareAk1P);
			}
			else {
				std::cout << "number threads=" <<  omp_get_max_threads() << std::endl;

#pragma omp parallel sections
				{
#pragma omp section
					{
						gfx::timsort(P + 1 + iaddR, P + isize21 + iaddR, compareAk1P);
					}
#pragma omp section
					{
						gfx::timsort(P + isize21 + 1 + iaddR, P + iaddR + nnzR - 1 + 1, compareAk1P);
					}
				}
				mergeTim_amg(P, 1 + iaddR, isize21 + iaddR, iaddR + nnzR - 1, indx_comparej);
			}
			break;
		default:
			Counting_Sort(P, 1 + iaddR, iaddR + nnzR - 1, n_a[ilevel - 1], indx_comparej);//подходит именно n_a[ilevel - 1]
			break;
		}

		if (row_ind_AS != nullptr) {
			free(row_ind_AS);
			row_ind_AS = nullptr;
		}
		row_ind_AS = my_declaration_array<integer>(numberofcoarcenodes, -1, "row_ind_AS");
		if (row_ind_AE != nullptr) {
			free(row_ind_AE);
			row_ind_AE = nullptr;
		}
		row_ind_AE = my_declaration_array<integer>(numberofcoarcenodes, -2, "row_ind_AE");

		
#pragma omp parallel for
		for (integer i = 1; i <= n; ++i) {
			flag[i] = false;
		}
		integer istart2 = nnz_a[ilevel - 1] + 1 + iadd;
		integer iend2 = istartAnew - 1;
		calculate_row_ptr(istart2, iend2, row_ind_AS, row_ind_AE, flag, Amat);

		if (bQuick_sort_for_reorder) {
			for (integer ii = 1; ii <= numberofcoarcenodes; ++ii) {
				// 14.04.2020
				// В каждой строке i индексы столбца j отсортированы по возрастанию.		
				qs_abbys_heigh = 0;
				// quicksort
				qsj(Amat, row_ind_AS[ii], row_ind_AE[ii]);
			}
		}

		// Инициализация чрезвычайно важна, т.к. 
		// обязательно присутствуют пустые строки которые
		// надо корректно обрабатывать.
		if (row_ind_PS != nullptr) {
			free(row_ind_PS);
			row_ind_PS = nullptr;
		}
		row_ind_PS = my_declaration_array<integer>(n_a[ilevel - 1], -1, "row_ind_PS");
		if (row_ind_PE != nullptr) {
			free(row_ind_PE);
			row_ind_PE = nullptr;
		}
		row_ind_PE = my_declaration_array<integer>(n_a[ilevel - 1], -2, "row_ind_PE");

		

		integer istart4 = 1 + iaddR;
		integer iend4 = nnzR - 1 + iaddR;
#pragma omp parallel for
		for (integer i = 1; i <= n_a[ilevel - 1]; ++i) {//n
			flag[i] = false;
		}
		calculate_row_ptr_j(istart4, iend4, row_ind_PS, row_ind_PE, flag, P);

		if (bQuick_sort_for_reorder) {
			for (integer ii = 1; ii <= n_a[ilevel - 1]; ++ii) {
				// 14.04.2020
				// В каждой строке j индексы столбца i отсортированы по возрастанию.		
				//qs_abbys_heigh = 0;
				// quicksort
				//qs(P, row_ind_PS[ii], row_ind_PE[ii]);
				std::sort(P + row_ind_PS[ii], P + row_ind_PE[ii] + 1, compareAk1R);
			}
		}

		//omp_set_num_threads(1);





		if ((number_cores() == 1)) {

			index_visit[0] = 0;

#pragma omp parallel for
			for (integer i = 1; i <= numberofcoarcenodes; ++i) {
				vector_sum[i] = 0.0;
				hash_table[i] = false;
				index_visit[i] = 0;
			}

			// Умножение разреженной матрицы A на разреженную матрицу P.
			// Результат записывается в матрицу А начиная с позиции istartAnew2.
			my_sparse_matrix_by_matrix_multiplication_AP<doublerealT>(Amat,
				P, istartAnew2,
				row_ind_AS, row_ind_AE, row_ind_PS, row_ind_PE, numberofcoarcenodes,
				// Вспомогательные данные.
				hash_table, index_visit, vector_sum, nsizeA,
				// для корректировки ошибочных строк
				// с отрицательной диагональю
				ibsp_length, i_bsp_LIMIT, ilevel, bsp);

		}
		else {

			//#ifdef _NONAME_STUB29_10_2017
#ifdef _OPENMP

		//for (int i_9 = 0; i_9 < iKnumber_thread; ++i_9) {
//#pragma omp parallel for
			//for (integer i_91 = 0; i_91 < n_a[ilevel - 1] + 1; ++i_91) {
				//vector_sum_m[i_9][i_91] = 0.0;
			//}
		//}

			my_parallel8_sparse_matrix_by_matrix_multiplication_AP<doublerealT>(Amat,
				P, istartAnew, istartAnew_m,
				row_ind_AS, row_ind_AE, row_ind_PS, row_ind_PE,
				numberofcoarcenodes,
				iKnumber_thread,
				hash_table_m, index_visit_m,
				vector_sum_m, index_size_m,
				nsizeA, numberofcoarcenodes/*n*/, nnz, ilevel,
				bprint_mesage_diagnostic,
				n_a, AccumulqtorA_m_SIZE8,
				AccumulqtorA_m, istartAnew2,
				// для корректировки ошибочных строк
				// с отрицательной диагональю
				ibsp_length, i_bsp_LIMIT, bsp, iVar,
				bcontinue_global);

			if (!bcontinue_global) {
				std::cout << "error in my_parallel8_sparse_matrix_by_matrix_multiplication_AP\n";
				system("pause");
			}

#else

		// Накопитель результата.
		//vector_sum = new doublerealT[numberofcoarcenodes + 1];
		/*
		// Память уже была освобождена выше по тексту программы.
		if (vector_sum != nullptr) {
			free(vector_sum);
			vector_sum = nullptr;
		}
		*//*
		// Данные используемые для частичного формирователя суммы.
		if (vector_sum != nullptr) {
			free(vector_sum);
			vector_sum = nullptr;
		}
		vector_sum = my_declaration_array<doublerealT>(numberofcoarcenodes, 0.0, "vector_sum");
		*/

		/*
		//integer size_v = sizeof(doublerealT)*(1 + numberofcoarcenodes);
		// Храним индексы ненулевых элементов в отсортированном порядке.
		if (index_visit != nullptr) {
			free(index_visit);
			index_visit = nullptr;
		}
		index_visit = my_declaration_array<integer>(n_a[ilevel - 1], 0, "index_visit");

		if (index_visit != nullptr) {
			index_visit[0] = 0;
		}*/

			index_visit[0] = 0;

#pragma omp parallel for
			for (integer i = 1; i <= numberofcoarcenodes; ++i) {
				vector_sum[i] = 0.0;
				hash_table[i] = false;
				index_visit[i] = 0;
			}

			// Умножение разреженной матрицы A на разреженную матрицу P.
			// Результат записывается в матрицу А начиная с позиции istartAnew2.
			my_sparse_matrix_by_matrix_multiplication_AP<doublerealT>(Amat,
				P, istartAnew2,
				row_ind_AS, row_ind_AE, row_ind_PS, row_ind_PE, numberofcoarcenodes,
				// Вспомогательные данные.
				hash_table, index_visit, vector_sum, nsizeA,
				// для корректировки ошибочных строк
				// с отрицательной диагональю
				ibsp_length, i_bsp_LIMIT, ilevel, bsp);

#endif

			}
		/*if (hash_table != nullptr) {
			free(hash_table);
			hash_table = nullptr;
		}
		
		if (vector_sum != nullptr) {
			free(vector_sum);
			vector_sum = nullptr;
		}
		if (index_visit != nullptr) {
			free(index_visit);
			index_visit = nullptr;
		}*/
		
		if (row_ind_AS != nullptr) {
			free(row_ind_AS);
			row_ind_AS = nullptr;
		}
		if (row_ind_AE != nullptr) {
			free(row_ind_AE);
			row_ind_AE = nullptr;
		}
		if (row_ind_PS != nullptr) {
			free(row_ind_PS);
			row_ind_PS = nullptr;
		}
		if (row_ind_PE != nullptr) {
			free(row_ind_PE);
			row_ind_PE = nullptr;
		}

		
		// Копируем матрицу А следующего уровня влево вплотную к матрице первоначального уровня.
		//integer icounter3 = 1;
		 integer nsize = istartAnew2 - (istartAnew);
		//doublereal mH = pow(icounter - 1, 0.33333);
		//doublereal alphaH = 2.0*(mH - 1)*(mH - 1) / ((2.0*mH - 1.0)*(2.0*mH - 1.0));
		//std::cout << "alphaH=" << alphaH << "\n";
		integer i_1start1 = nnz_a[ilevel - 1] + 1 + iadd;

#pragma omp parallel for
		for (integer i_2 = 1; i_2 <= nsize; ++i_2)
		{
			integer i_1 = i_1start1 - 1 + i_2;
			integer i_right_position = istartAnew - 1 + i_2;
			//if (my_amg_manager.baglomeration_with_consistency_scaling > 0) {
				// agglomeration with consistency scaling 24.05.2019
				//Amat.aij[i_1] = alphaH*Amat.aij[i_right_position];
			//}
			//else {
				Amat.aij[i_1] = Amat.aij[i_right_position];
			//}
			Amat.i[i_1] = Amat.i[i_right_position];
			Amat.j[i_1] = Amat.j[i_right_position];
		}

		if (bprint_mesage_diagnostic) {
			std::cout << "Prolongation is construct.\n";			 

			// Общее количество узлов F у которых нет соседних С узлов.
			std::cout << "diagnostic: the number of neighbors that are not Coarse (C) nodes " << the_number_of_neighbors_that_are_not_C_nodes << "\n";
			// Количество F узлов у которых только один интерполяционный С сосед.
			std::cout << "diagnostic: the number of Fine (F) nodes with one single strong Coarse (C) neighbor=" << number_of_F_nodes_with_one_single_strong_C_neighbor << " \n";
			std::cout << "diagnostic: the number of Fine (F) nodes with one single strong Coarse (C) neighbor\n";
			std::cout << "and to the same not having strong Fine(F) neighbors " << number_of_F_nodes_with_one_single_strong_C_neighborF << "\n";
			//system("pause");

		}
		if (debug_reshime) system("pause");


		//delete[] C_numerate;
		if (C_numerate != nullptr) {
			free(C_numerate);
			C_numerate = nullptr;
		}

		// Использование упорядочивания типа F/C ускоряет сходимость вычислительного процесса,
		// сокращая число V циклов требуемых для достижения сходимости.
		integer iaddFCcolor = 0;
#pragma omp parallel for reduction(+ : iaddFCcolor)
		for (integer i_71 = 0; i_71 < ilevel - 1; ++i_71) {
			iaddFCcolor += n_a[i_71];
		}
		const integer istop_i_1 = n_a[ilevel - 1];
#pragma omp parallel for 
		for (integer i_1 = 1; i_1 <= istop_i_1; ++i_1)
		{
			if (this_is_C_node[i_1]) {
				F_false_C_true[iaddFCcolor + i_1] = true;
			}
		}

		nnz_aRP[ilevel - 1] = nnzR - 1;
		iaddR += nnzR - 1;
		n_a[ilevel] = icounter - 1;
		nnz_a[ilevel] = nsize;
		iadd += nnz_a[ilevel - 1];



		if (bcontinue_global) {
			// если bad string не встречалось.
			ilevel++;	
		}
		

		if (bStrongTransposeON) {
			// Освобождение ОЗУ.

			// Обычный линейный список.
			/*
			if (hash_StrongTranspose_collection1 != nullptr) {
				//for (integer i_1 = 0; i_1 <= n_a[ilevel - 2]; ++i_1)
				//isize_memory_alloc_hash_StrongTranspose_collection1
				for (integer i_1 = 0; i_1 <= isize_memory_alloc_hash_StrongTranspose_collection1; ++i_1)
				{
					clear_list(hash_StrongTranspose_collection1[i_1]);
				}
				delete[] hash_StrongTranspose_collection1;
				hash_StrongTranspose_collection1 = nullptr;
			}*/

			/*if (isize_hash_StrongTranspose_collection != nullptr) {
				delete[] isize_hash_StrongTranspose_collection;
				isize_hash_StrongTranspose_collection = nullptr;
			}*/
		}

	
		if (count_neighbour != nullptr) {
			free(count_neighbour);
			count_neighbour = nullptr;
		}
		
		if (row_startA != nullptr) {
			free(row_startA);
			row_startA = nullptr;
		}


		// построение иерархии уровней досрочно прекращено.
	//BAD_STRING_MARKER:

		
		if (bprint_mesage_diagnostic) {
			std::cout << "one level construct OK." << std::endl;
		}
		if (debug_reshime) system("pause");	

		//проверка конец

	} // иерархия сеток построена.

	
	delete[] hash_StrongTranspose_collection1Eco; // 10.06.2021
	if (isize_hash_StrongTranspose_collection != nullptr) {
		delete[] isize_hash_StrongTranspose_collection;
		isize_hash_StrongTranspose_collection = nullptr;
	}

	if (vector_sum != nullptr) {
		free(vector_sum);
		vector_sum = nullptr;
	}

	if (index_visit != nullptr) {
		free(index_visit);
		index_visit = nullptr;
	}

	if (hash_table != nullptr) {
		free(hash_table);
		hash_table = nullptr;
	}

	if (bprint_mesage_diagnostic) {
		std::cout << "   ***   CAMG ITERATOR   ***" << std::endl;
	}


	 // Освобождение памяти используемой на этапе построения иерархии матриц.
	 // Освобождение оперативной памяти.
	if (threshold_quick_all != nullptr) {
		free(threshold_quick_all);
		threshold_quick_all = nullptr;
	}

	if (threshold_quick_only_negative != nullptr) {
		free(threshold_quick_only_negative);
		threshold_quick_only_negative = nullptr;
	}	

	if (istack != nullptr) {
		free(istack);
		istack = nullptr;
	}

	if (hash_table2 != nullptr) {
		free(hash_table2);
		hash_table2 = nullptr;
	}

	if (this_is_C_node != nullptr) {
		free(this_is_C_node);
		this_is_C_node = nullptr;
	}
	if (this_is_F_node != nullptr) {
		free(this_is_F_node);
		this_is_F_node = nullptr;
	}

	ilevel--; // 4.01.2017
	if (n_a[ilevel] < 5) {
		// Чтобы не было последних уровней где меньше 5 узлов сетки.
		ilevel--;
	}

	//omp_set_num_threads(8);

	// Эта память в два вектора длины n_a[0] была нужна для алгоритма сортировки CountingSort.
	free(C1);
	free(C2);

	if ((my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2_full) ||
		(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS2) ||
		(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::PMIS) ||
		(my_amg_manager.icoarseningtype == MY_AMG_SPLITTING_COARSENING_ALGORITHM::HMIS))
	{
		free(Bm1);
		free(Bm2);
	}

	// Вычисляем и запоминаем grid complexity
	// Операторная сложность.
	doublerealT dr_grid_complexity = static_cast<doublerealT>((((double)(1.0*iadd)) / ((double)(1.0*nnz_a[0]))));
	if (bprint_mesage_diagnostic) {

		std::cout << "operator complexity is Co = " << dr_grid_complexity << std::endl;
		doublerealT gridCG = 1.0;
		for (int i93 = 1; i93 <= ilevel; ++i93) gridCG += 1.0 * n_a[i93] / (1.0*n_a[0]);
		std::cout << "grid complexity is Cg = " << gridCG << std::endl;
		std::cout << "Prolongation operator complexity is |Psigma|/|P1|=" << static_cast<doublerealT>(1.0 * nnz_P_memo_all / nnz_P_memo_0) << " " << static_cast<doublerealT>(1.0 * nnz_P_memo_all / n_a[0]) << "  * n\n"; 
		doublerealT sizegb = static_cast<doublerealT>(16 * iadd / 1.0e9);
		std::cout << "memory usage is "<< sizegb <<" Gb.reserved "<< 16 * nsizeA / 1.0e9 << " Gb.ratio is equal = " << sizegb / (16 * nsizeA / 1.0e9) << "\n"; 
	
		
	    // 31.224s [50.986] 2D m=81 debug x64 acumulqtor
	    // 13.792 [18.156] 2D m=81 realese x64 acumulqtor
	    // 8.028s 2D m=81 debug x64 rozetka
	    // 3.827 2D m=81 realese x64 rozetka
	
		std::cout << "number of levels=" << ilevel << std::endl;
		std::cout << "levels    unknowns			nonzeros     sample_pattern"<<std::endl;
		// <= ilevel 4.01.2017
		for (integer i_1 = 0; i_1 <= ilevel; ++i_1) {
			if (i_1 == 0) {
				std::cout << std::setw(2) << std::right << i_1 << " \t " << std::setw(9) << std::right << n_a[i_1]<<  "            " << std::setw(11) << std::right << nnz_a[i_1] <<  "         "<< "\t" << std::setfill(' ') << " " << std::setw(4) << std::right << static_cast<integer>(nnz_a[i_1] / n_a[i_1]) << "\n";
			}
			else {
				std::cout << std::setw(2) << std::right << i_1 << " \t " << std::setw(9) << std::right << n_a[i_1] << " " << std::setfill(' ') << std::setw(7) << std::right <<  (100.0 * n_a[i_1] / n_a[i_1 - 1]) << "%   " << std::setw(11) << std::right << nnz_a[i_1] << " " << std::setfill(' ') << std::setw(7) << std::right << (100.0 * nnz_a[i_1] / nnz_a[i_1 - 1]) << "% \t " << std::setw(4) << std::right << static_cast<integer>(nnz_a[i_1] / n_a[i_1]) << "\n";
			}
		}
		std::cout << "Graph(Mesh) ierarhion is construct sucsseful..." << std::endl;
	}

	if (debug_reshime) system("pause");
	//system("pause");
	//exit(1);

	if (bprint_mesage_diagnostic) {
		std::cout << "memory optimization 13 november 2016."<< std::endl;
		std::cout << "ierarhion matrix Amat..."<< std::endl;
	}
		// Уменьшение памяти отводимой под хранение матрицы А.
		// Матрица должна занимать в памяти не более чем под неё нужно и не мегабайтом больше.
		my_realloc_memory<real_mix_precision>(Amat.aij, ((iadd + 2)));
		
		
		if (Amat.abs_aij != nullptr) {
			//delete[] Amat;
			free(Amat.abs_aij);
			Amat.abs_aij = nullptr;
		}

		my_realloc_memory<integer_mix_precision>(Amat.i, ((iadd + 2)));
				
		my_realloc_memory<integer_mix_precision>(Amat.j, ((iadd + 2)));
		
		
		if (bprint_mesage_diagnostic) {
			std::cout << " 1 of 3 compleated.  OK!! ierarhion matrix Amat realloc successfully..." << std::endl;
			std::cout << "Prolongation ierarhion..." << std::endl;
		}
		my_realloc_memory<Ak1>(P, (static_cast<integer>(nnz_P_memo_all)+2));
				
		if (bprint_mesage_diagnostic) {
			std::cout << "2 of 3 compleated. OK!! ierarhion matrix Prolongation realloc successfully..." << std::endl;
		}

		if ((number_cores() == 1)) {

			// Однопоточный код. Никакой дополнительной памяти выделять(освобождать) ненадо. 
			// Это сильно экономит время.

		}
		else {

		//#ifdef	_NONAME_STUB29_10_2017
#ifdef _OPENMP
	// Освобождение озу ГУСТАВСОН умножение разреженных матриц.
	// Единожды!!!
		for (int i_9 = 0; i_9 < iKnumber_thread; ++i_9) {
			//free(vector_sum_m[i_9]);
			//free(index_visit_m[i_9]);
			//free(hash_table_m[i_9]);

			delete[] vector_sum_m[i_9];
			delete[] index_visit_m[i_9];
			delete[] hash_table_m[i_9];
		}
		delete[] vector_sum_m;
		delete[] index_visit_m;
		delete[] hash_table_m;
		delete[] index_size_m;
		vector_sum_m = nullptr;
		index_visit_m = nullptr;
		hash_table_m = nullptr;
		index_size_m = nullptr;

		for (int i_9 = 0; i_9 < iKnumber_thread; ++i_9) {
			delete[] AccumulqtorA_m[i_9];
			AccumulqtorA_m[i_9] = nullptr;
		}
		delete[] AccumulqtorA_m;
		AccumulqtorA_m = nullptr;
		delete[] istartAnew_m;
		istartAnew_m = nullptr;
#endif

		}

		amg_pp.maxlevel = maxlevel;
		amg_pp.ilevel = ilevel;
		amg_pp.bprint_mesage_diagnostic = bprint_mesage_diagnostic;
		amg_pp.bonly_serial = bonly_serial;
		amg_pp.bILU2smoother = bILU2smoother;
		amg_pp.nnz_P_memo_0 = nnz_P_memo_0;
		amg_pp.nnz_P_memo_all = nnz_P_memo_all;
		amg_pp.dr_grid_complexity = dr_grid_complexity;
		amg_pp.debug_reshime = debug_reshime;
		amg_pp.dapply_ilu_max_pattern_size = dapply_ilu_max_pattern_size;
		amg_pp.RealZERO = RealZERO;
		amg_pp.identiti = identiti;
		amg_pp.memo_icoarseningtype = memo_icoarseningtype;
		

} // setup_phase_classic_aglomerative_amg6

// Версия classic_aglomerative_amg6 на основе версии
// classic_aglomerative_amg4. 
// Setup и Solution фазы отделены друг от друга.
// Теперь в classic_aglomerative_amg6 содержится только вызов 
// setup и solution фаз. См. отдельные функции для подробностей.
// В версии classic_aglomerative_amg4 разделение Setup и Solution фаз
// произведено не было.
template <typename doublerealT>
bool classic_aglomerative_amg6(Ak2& Amat,
	integer nsizeA, // количество ячеек выделенное извне для хранилища матриц А	
	integer nnz, // number of non zero elements
	integer n, // dimension of vectors x and b.	
	doublereal*& x, //solution (решение) 
	doublereal*& b, // rthdsd (правая часть).
	real_mix_precision& ret74,
	integer iVar,
	bool bmemory_savings,
	BLOCK*& my_body, int& lb, integer maxelm_out,
	int * &whot_is_block
) {

	// Замер времени.
	unsigned int calculation_start_time_setup = 0; // начало счёта мс.
	unsigned int calculation_end_time_setup = 0; // окончание счёта мс.
	unsigned int calculation_seach_time_setup = 0; // время выполнения участка кода в мс.

	// Замер времени.
	unsigned int calculation_start_time_solve = 0; // начало счёта мс.
	unsigned int calculation_end_time_solve = 0; // окончание счёта мс.
	unsigned int calculation_seach_time_solve = 0; // время выполнения участка кода в мс.

	calculation_start_time_setup = clock(); // момент начала счёта.

	amg_precond_param amg_pp;

	//prolongation_gl = new doublereal[n* number_cores()+1];
	//bflag_equal= new bool[n + 1];
	//for (int i_62 = 0; i_62 <= n; ++i_62) {
		//bflag_equal[i_62] = false;
	//}
	//pool_prolongation = new integer[n + 1];
	//pool_prolongation_acum = new integer[n + 1];
	//pool_prolongation_acum2 = new integer[n + 1];
	
	// Не более 100 ошибочных строк
	// с отрицательной диагональю.
	const int i_bsp_LIMIT = 100;
	BAD_STRING_PATCHING* bsp = nullptr;
	integer ibsp_length = 0;

	// Оператор проекции равен оператору интерполяции с точностью до транспонирования.
	Ak1* P = nullptr; // оператор интерполяции.

	Ak2 Amat_copy;

#if MY_GPU_MATH_LIBRARY_CU_ON

	Amat_copy.aij = new real_mix_precision[nnz];
	Amat_copy.i = new integer_mix_precision[nnz];
	Amat_copy.j = new integer_mix_precision[nnz];

	for (integer i = 0; i < nnz; ++i) {
		Amat_copy.aij[i] = Amat.aij[i+1];
		Amat_copy.i[i] = Amat.i[i+1];
		Amat_copy.j[i] = Amat.j[i+1];
	}

#else

	Amat_copy.aij = nullptr;
	Amat_copy.i = nullptr;
	Amat_copy.j = nullptr;

#endif

	integer* n_a = nullptr; // число неизвестных на каждом из уровней.
	integer* nnz_a = nullptr; // число ненулевых коэффициентов в матрице А
	bool* flag = nullptr; 

	integer* nnz_aRP = nullptr; // число ненулевых элементов в операторе интерполяции P на каждом из уровней.
	bool* F_false_C_true = nullptr; // C/F - разбиение.

	// дополнительное упорядочивание по столбцам в строке в 
	// надежде что кеш будет лучше использоваться.
	// Его можно безболезненно выключить.
	const bool bQuick_sort_for_reorder = false;

	// отладка. Причина использование миксовой точности в задачах на излучение. Нужно использовать только double.
	//std::cout << "apriorishe x=" << norma(x, n) <<" b=" << norma(b, n) <<"\n";
	//system("pause");

	// Фаза подготовки вызывается один раз для фиксированной матрицы СЛАУ.
	setup_phase_classic_aglomerative_amg6<doublerealT>(Amat,
		nsizeA, // количество ячеек выделенное извне для хранилища матриц А	
		nnz, // number of non zero elements
		n, // dimension of vectors x and b.	
		x, //solution (решение) 
		b, // rthdsd (правая часть).
		ret74,
		iVar,
		bmemory_savings,
		// Параметры нужные только для solution phase.
		amg_pp,
		n_a, nnz_a, nnz_aRP,
		ibsp_length,
		bsp,
		i_bsp_LIMIT,
		flag,
		F_false_C_true, P,
		bQuick_sort_for_reorder);

	calculation_end_time_setup = clock();
	calculation_start_time_solve = clock(); // момент начала счёта.

	INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14 imyinit = INIT_SELECTOR_CASE_CAMG_RUMBAv_0_14::ZERO_INIT;

	// Для ускорения паралельной версии пролонгации 1.08.2021.
	//initQuickProlongation(n);

	level_Chebyshev_info = new LEVEL_CHEBYSHEV_INFO[amg_pp.maxlevel];
	for (int i = 0; i < amg_pp.maxlevel; ++i) {
		level_Chebyshev_info[i].bGershgorin_calc_Temp = false;
		level_Chebyshev_info[i].bstart_lo_Temp = false;
		level_Chebyshev_info[i].bGershgorin_calc_Vx = false;
		level_Chebyshev_info[i].bstart_lo_Vx = false;
		level_Chebyshev_info[i].bGershgorin_calc_Vy = false;
		level_Chebyshev_info[i].bstart_lo_Vy = false;
		level_Chebyshev_info[i].bGershgorin_calc_Vz = false;
		level_Chebyshev_info[i].bstart_lo_Vz = false;

		level_Chebyshev_info[i].bGershgorin_calc_Pressure = false;
		level_Chebyshev_info[i].bstart_lo_Pressure = false;

		level_Chebyshev_info[i].bGershgorin_calc_Stress = false;
		level_Chebyshev_info[i].bstart_lo_Stress = false;

		level_Chebyshev_info[i].bGershgorin_calc_Speed = false;
		level_Chebyshev_info[i].bstart_lo_Speed = false;

		level_Chebyshev_info[i].bGershgorin_calc_nu = false;
		level_Chebyshev_info[i].bstart_lo_nu = false;

		level_Chebyshev_info[i].bGershgorin_calc_k = false;
		level_Chebyshev_info[i].bstart_lo_k = false;

		level_Chebyshev_info[i].bGershgorin_calc_omega = false;
		level_Chebyshev_info[i].bstart_lo_omega = false;

		level_Chebyshev_info[i].bGershgorin_calc_gamma = false;
		level_Chebyshev_info[i].bstart_lo_gamma = false;

		level_Chebyshev_info[i].bGershgorin_calc_ReTheta = false;
		level_Chebyshev_info[i].bstart_lo_ReTheta = false;

		level_Chebyshev_info[i].bGershgorin_calc_epsilon = false;
		level_Chebyshev_info[i].bstart_lo_epsilon = false;

		level_Chebyshev_info[i].bGershgorin_calc_k_for_ke = false;
		level_Chebyshev_info[i].bstart_lo_k_for_ke = false;
	}


	// Фаза решения может вызываться многократно,
	// требуется обновлять правую часть СЛАУ.
	// Матрица СЛАУ не должна меняться.
	
	bool ret_value = false;
	ret_value = solution_phase<doublerealT>(Amat, Amat_copy,
		nsizeA, // количество ячеек выделенное извне для хранилища матриц А	
		nnz, // number of non zero elements
		n, // dimension of vectors x and b.	
		x, // solution (решение) 
		b, // rthdsd (правая часть).
		ret74,
		iVar,
		bmemory_savings,
		my_body, lb, maxelm_out,
		// данные для solution phase
		amg_pp,
		n_a, nnz_a, nnz_aRP,
		ibsp_length,
		bsp,
		i_bsp_LIMIT,
		flag,
		F_false_C_true, P,
		whot_is_block,
		imyinit);

	// Для ускорения паралельной версии пролонгации 1.08.2021.
	//freeQuickProlongation();

#if MY_GPU_MATH_LIBRARY_CU_ON

	if (Amat_copy.aij != nullptr) {
		delete[] Amat_copy.aij;
		Amat_copy.aij = nullptr;
	}
	if (Amat_copy.i != nullptr) {
		delete[] Amat_copy.i;
		Amat_copy.i = nullptr;
	}
	if (Amat_copy.j != nullptr) {
		delete[] Amat_copy.j;
		Amat_copy.j = nullptr;
	}

#endif
	
	if (P != nullptr) {
		//delete[] P;
		free(P);
		P = nullptr;
	}

	if (F_false_C_true != nullptr) {
		free(F_false_C_true);
		F_false_C_true = nullptr;
	}

	delete[] bsp;

	delete[] nnz_aRP;
	delete[] n_a;
	delete[] nnz_a;

	if (flag != nullptr) {
		free(flag);
		flag = nullptr;
	}

	//delete[] prolongation_gl;
	//prolongation_gl = nullptr;
	//delete[] bflag_equal;
	//bflag_equal = nullptr;
	//delete[] pool_prolongation;
	//pool_prolongation = nullptr;
	//delete[] pool_prolongation_acum;
	//pool_prolongation_acum = nullptr;
	//delete[] pool_prolongation_acum2;
	//pool_prolongation_acum2 = nullptr;

	calculation_end_time_solve = clock();

	
	calculation_seach_time_setup = calculation_end_time_setup - calculation_start_time_setup;
	unsigned int im_setup = 0, is_setup = 0, ims_setup = 0;
	im_setup = (unsigned int)(calculation_seach_time_setup / 60000); // минуты
	is_setup = (unsigned int)((calculation_seach_time_setup - 60000 * im_setup) / 1000); // секунды
	ims_setup = (unsigned int)((calculation_seach_time_setup - 60000 * im_setup - 1000 * is_setup) / 10); // миллисекунды делённые на 10

	calculation_seach_time_solve = calculation_end_time_solve - calculation_start_time_solve;
	unsigned int im_solve = 0, is_solve = 0, ims_solve = 0;
	im_solve = (unsigned int)(calculation_seach_time_solve / 60000); // минуты
	is_solve = (unsigned int)((calculation_seach_time_solve - 60000 * im_solve) / 1000); // секунды
	ims_solve = (unsigned int)((calculation_seach_time_solve - 60000 * im_solve - 1000 * is_solve) / 10); // миллисекунды делённые на 10

	if (my_amg_manager.iprint_log) {
		std::cout << "setup time algebraic multigrid RUMBA.v.0.14:  "<< im_setup <<" minute "<< is_setup <<" second "<< 10 * ims_setup <<" millisecond\n";
		std::cout << "solution time algebraic multigrid RUMBA.v.0.14:  "<< im_solve <<" minute "<< is_solve <<" second "<< 10 * ims_solve <<" millisecond\n";
	}

	return ret_value;
} // classic_aglomerative_amg6


// Передает пользовательские настройки решающему алгоритму РУМБА 0.14
void set_RUMBA_Classic_AMG_Setting(integer iVar)
{
	// настройка параметров РУМБА 0.14 решателя.
	switch (iVar) {
	case TEMP:
		my_amg_manager.theta = static_cast<real_mix_precision>(my_amg_manager.theta_Temperature);
		my_amg_manager.maximum_delete_levels = my_amg_manager.maximum_delete_levels_Temperature;
		my_amg_manager.nFinnest = my_amg_manager.nFinnest_Temperature;
		my_amg_manager.nu1 = my_amg_manager.nu1_Temperature;
		my_amg_manager.nu2 = my_amg_manager.nu2_Temperature;
		my_amg_manager.memory_size = my_amg_manager.memory_size_Temperature;
		//my_amg_manager.memory_size = 5.0; // вместо 9.0 даёт понижение используемой ОЗУ на 15.6%.
		my_amg_manager.ilu2_smoother = my_amg_manager.ilu2_smoother_Temperature;
		my_amg_manager.icoarseningtype = my_amg_manager.icoarseningTemp; // standart vs RS 2.
		my_amg_manager.istabilization = my_amg_manager.istabilizationTemp; // Stabilization: 0 - none, 1 - bicgstab + amg (РУМБА), 2 - FGMRes + amg (РУМБА).
		my_amg_manager.magic = my_amg_manager.F_to_F_Temperature; // magic
		my_amg_manager.number_interpolation_procedure = my_amg_manager.number_interpolation_procedure_Temperature;
		my_amg_manager.iCFalgorithm_and_data_structure = my_amg_manager.iCFalgorithm_and_data_structure_Temperature;
		my_amg_manager.iprint_log = my_amg_manager.iprint_log_Temperature;
		my_amg_manager.itruncation_interpolation = my_amg_manager.itruncation_interpolation_Temperature;
		my_amg_manager.truncation_interpolation = my_amg_manager.truncation_interpolation_Temperature;
		my_amg_manager.gold_const = my_amg_manager.gold_const_Temperature;
		my_amg_manager.b_gmres = my_amg_manager.b_gmresTemp;
		my_amg_manager.b_ChebyshevSmoother = my_amg_manager.b_ChebyshevSmootherTemp;
		my_amg_manager.bMatrixPortrait = my_amg_manager.bTemperatureMatrixPortrait;
		my_amg_manager.bcf_reorder = my_amg_manager.bcf_reorder_Temperature;
		my_amg_manager.bthreshold_auto = my_amg_manager.bthreshold_Temperature_auto;
		break;
	case PAM:
		my_amg_manager.theta = static_cast<real_mix_precision>(my_amg_manager.theta_Pressure);
		my_amg_manager.maximum_delete_levels = my_amg_manager.maximum_delete_levels_Pressure;
		my_amg_manager.nFinnest = my_amg_manager.nFinnest_Pressure;
		my_amg_manager.nu1 = my_amg_manager.nu1_Pressure;
		my_amg_manager.nu2 = my_amg_manager.nu2_Pressure;
		my_amg_manager.memory_size = my_amg_manager.memory_size_Pressure;
		my_amg_manager.ilu2_smoother = my_amg_manager.ilu2_smoother_Pressure;
		my_amg_manager.icoarseningtype = my_amg_manager.icoarseningPressure; // standart vs RS 2.
		my_amg_manager.istabilization = my_amg_manager.istabilizationPressure; // Stabilization: 0 - none, 1 - bicgstab + amg (РУМБА), 2 - FGMRes + amg (РУМБА).
		my_amg_manager.magic = my_amg_manager.F_to_F_Pressure; // magic
		my_amg_manager.number_interpolation_procedure = my_amg_manager.number_interpolation_procedure_Pressure;
		my_amg_manager.iCFalgorithm_and_data_structure = my_amg_manager.iCFalgorithm_and_data_structure_Pressure;
		my_amg_manager.iprint_log = my_amg_manager.iprint_log_Pressure;
		my_amg_manager.itruncation_interpolation = my_amg_manager.itruncation_interpolation_Pressure;
		my_amg_manager.truncation_interpolation = my_amg_manager.truncation_interpolation_Pressure;
		my_amg_manager.gold_const = my_amg_manager.gold_const_Pressure;
		my_amg_manager.b_gmres = my_amg_manager.b_gmresPressure;
		my_amg_manager.b_ChebyshevSmoother = my_amg_manager.b_ChebyshevSmootherPressure;
		my_amg_manager.bMatrixPortrait = my_amg_manager.bPressureMatrixPortrait;
		my_amg_manager.bcf_reorder = my_amg_manager.bcf_reorder_Pressure;
		my_amg_manager.bthreshold_auto = my_amg_manager.bthreshold_Pressure_auto;
		break;
		// 10.10.2019 Для турбулентных характеристик настройка решателя такая же как и для компонент скорости.
	case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT:
	case NUSHA: case TURBULENT_KINETIK_ENERGY: case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
	case TURBULENT_KINETIK_ENERGY_STD_K_EPS: case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
	case GAMMA_LANGTRY_MENTER: case RE_THETA_LANGTRY_MENTER:
		my_amg_manager.theta = static_cast<real_mix_precision>(my_amg_manager.theta_Speed);
		my_amg_manager.maximum_delete_levels = my_amg_manager.maximum_delete_levels_Speed;
		my_amg_manager.nFinnest = my_amg_manager.nFinnest_Speed;
		my_amg_manager.nu1 = my_amg_manager.nu1_Speed;
		my_amg_manager.nu2 = my_amg_manager.nu2_Speed;
		my_amg_manager.memory_size = my_amg_manager.memory_size_Speed;
		my_amg_manager.ilu2_smoother = my_amg_manager.ilu2_smoother_Speed;
		my_amg_manager.icoarseningtype = my_amg_manager.icoarseningSpeed; // standart vs RS 2.
		my_amg_manager.istabilization = my_amg_manager.istabilizationSpeed; // Stabilization: 0 - none, 1 - bicgstab + amg (РУМБА), 2 - FGMRes + amg (РУМБА).
		my_amg_manager.magic = my_amg_manager.F_to_F_Speed; // magic
		my_amg_manager.number_interpolation_procedure = my_amg_manager.number_interpolation_procedure_Speed;
		my_amg_manager.iCFalgorithm_and_data_structure = my_amg_manager.iCFalgorithm_and_data_structure_Speed;
		my_amg_manager.iprint_log = my_amg_manager.iprint_log_Speed;
		my_amg_manager.itruncation_interpolation = my_amg_manager.itruncation_interpolation_Speed;
		my_amg_manager.truncation_interpolation = my_amg_manager.truncation_interpolation_Speed;
		my_amg_manager.gold_const = my_amg_manager.gold_const_Speed;
		my_amg_manager.b_gmres = my_amg_manager.b_gmresSpeed;
		my_amg_manager.b_ChebyshevSmoother = my_amg_manager.b_ChebyshevSmootherSpeed;
		my_amg_manager.bMatrixPortrait = my_amg_manager.bSpeedMatrixPortrait;
		my_amg_manager.bcf_reorder = my_amg_manager.bcf_reorder_Speed;
		my_amg_manager.bthreshold_auto = my_amg_manager.bthreshold_Speed_auto;
		break;
	default :
		std::cout << "set_RUMBA_Classic_AMG_Setting unknown iVar\n";
		system("pause");
		break;
	}
}// set_RUMBA_Classic_AMG_Setting


bool classic_aglomerative_amg6_for_NetworkT(doublereal*& val, integer*& col_ind, integer*& row_ptr,
	integer n, integer nnz, doublereal*& rthdsd, doublereal* & potent, BLOCK* &b, int &lb	
) {

	set_RUMBA_Classic_AMG_Setting(TEMP);

	// 0 - Убираем подробную печать лога на консоль.
	my_amg_manager.iprint_log = 0;


	// РУМБА v.0.14

							// Свой собственный многосеточный метод, работает с alphaA строго меньше 0.9.
							// работает с alphaA=0.1, но очень медленно. 

	integer nsizeA = 13 * nnz;//13 слишком много, можно ужимать.
	Ak2 Amat;
	Amat.i = (integer_mix_precision*)malloc(nsizeA * sizeof(integer_mix_precision));
	Amat.j = (integer_mix_precision*)malloc(nsizeA * sizeof(integer_mix_precision));
	Amat.aij = (real_mix_precision*)malloc(nsizeA * sizeof(real_mix_precision)); ;
	Amat.abs_aij = (real_mix_precision*)malloc(nsizeA * sizeof(real_mix_precision));

	doublereal* x_copy = nullptr;
	x_copy = (doublereal*)malloc((n + 1) * sizeof(doublereal));

	doublereal* rthdsd_amg = nullptr;
	rthdsd_amg = (doublereal*)malloc((n + 1) * sizeof(doublereal));

	if ((x_copy != NULL) && (rthdsd_amg != NULL)) {

		for (integer i47 = 1; i47 <= n; ++i47) {
			x_copy[i47] = potent[i47 - 1];
			//x_copy[i47] = 0.001*(rand() % 100); // Случайное число от нуля до 0.1
			rthdsd_amg[i47] = rthdsd[i47 - 1];
		}
	}
	else {
		std::cout << "x_copy or rthdsd_amg == NULL memory alloc error in function classic_aglomerative_amg6_for_NetworkT.\n";
		system("pause");
		exit(1);
	}



	for (integer i_78 = 0; i_78 < n; ++i_78) {
		// Первая всегда диагональ.
		for (integer i_79 = row_ptr[i_78]; i_79 < row_ptr[i_78 + 1]; ++i_79) {
			if (col_ind[i_79] == i_78) {
				integer diagonal = row_ptr[i_78] + 1;
				Amat.aij[diagonal] = static_cast<real_mix_precision>(val[i_79]);
				if (val[i_79] <= 0.0) {
					std::cout << "diag aij=" << val[i_79] << " i=" << i_79 << std::endl;
				}
				Amat.abs_aij[diagonal] = static_cast<real_mix_precision>(fabs(val[i_79]));
				Amat.j[diagonal] = static_cast<integer_mix_precision>(col_ind[i_79] + 1); //  индексация столбца начинается с единицы
			}
		}

		integer post = row_ptr[i_78] + 1 + 1;
		for (integer i_79 = row_ptr[i_78]; i_79 < row_ptr[i_78 + 1]; ++i_79) {
			if (col_ind[i_79] != i_78) {
				Amat.aij[post] = static_cast<real_mix_precision>(val[i_79]);
				if (val[i_79] >= 0.0) {
					std::cout << "aij=" << val[i_79] << " i=" << i_79 << std::endl;
				}
				Amat.abs_aij[post] = static_cast<real_mix_precision>(fabs(val[i_79]));
				Amat.j[post] = static_cast<integer_mix_precision>(col_ind[i_79] + 1); //  индексация столбца начинается с единицы
				post++;
			}
		}

		for (integer i_79 = row_ptr[i_78]; i_79 < row_ptr[i_78 + 1]; ++i_79) {
			Amat.i[i_79 + 1] = static_cast<integer_mix_precision>(i_78 + 1); // номер строки который начинается с единиуцы.
		}
	}

	// Инициализация нулём.
	for (integer i_78 = nnz + 1; i_78 < nsizeA - 1; ++i_78) {

		Amat.i[i_78] = 0;
		Amat.j[i_78] = 0;
		Amat.aij[i_78] = 0.0;
		Amat.abs_aij[i_78] = 0.0;

	}


	for (integer i_78 = 1; i_78 <= nnz; ++i_78) {
		//std::cout << "i=" << Amat.i[i_78] << " j=" << Amat.j[i_78] << " aij=" << Amat.aij[i_78] << " abs_aij=" << Amat.abs_aij[i_78] << std::endl;
		// Проверка на корректность.
		if ((Amat.i[i_78] < 1) || (Amat.i[i_78] > n) || (Amat.j[i_78] < 1) || (Amat.j[i_78] > n)) {
			system("PAUSE");
		}
	}


	//system("PAUSE");


	real_mix_precision ret74 = 0.0;
	bool bmemory_savings = false;
	int* whot_is_block_out = nullptr;

	// Везде нумерация начинается с единицы
	classic_aglomerative_amg6<real_mix_precision>(Amat, nsizeA, nnz, n, x_copy, rthdsd_amg,
		ret74, TEMP, bmemory_savings, b, lb, -1, whot_is_block_out);


	for (integer i47 = 1; i47 <= n; ++i47) {
		potent[i47 - 1] = x_copy[i47];
		rthdsd[i47 - 1] = rthdsd_amg[i47];
	}

	free(Amat.i);
	free(Amat.j);
	free(Amat.aij);
	free(Amat.abs_aij);
	free(x_copy);
	free(rthdsd_amg);

	//system("pause");

	

	return true;
}

// Специальная нелинейная версия amg1r5 алгоритма.
#include "amg1r5_nonlinear.cpp"

#endif /*CLASSIC_AGLOMERATIVE_AMG6_2018YEAR_CPP*/
// Файл my_linalg.cpp
// самостоятельная реализация некоторых функций линейной алгебры.


#pragma once
#ifndef MY_LINALG_CPP
#define MY_LINALG_CPP 1

#include <stdio.h> // для функции getchar
#include <stdlib.h> // Для функции exit, atoi, atof
#include <math.h> // математические функции sqrt, fabs
#include "my_gpu_math_library.cpp"

#include "my_linalg.h" // самописные функции линейной алгебры
#include "ilut.cpp" // библиотека Юсефа Саада транслированная с помощью f2c.exe;
#include <ctime> // для замера времени выполнения.
#include <iostream> // для _finite




#include "my_cusp_alg.cpp" // Cusp 0.5.1
//#include "my_paralution.cpp" 
// закомментировать #include "my_vienna_alg.cpp"  если она не используется.
// Задать GPU_LIB_INCLUDE_MY_PROJECT_vienna = 0; Если viennacl 1.7.1 lib не используется.
const integer GPU_LIB_INCLUDE_MY_PROJECT_vienna = 0;
//#define AMGCL_INCLUDE_IN_MY_PROJECT 1
//#define VIENNA_CL_INCLUDE_IN_MY_PRJ 1

#ifdef VIENNA_CL_INCLUDE_IN_MY_PRJ
//#include "my_vienna_alg.cpp" // ViennaCL 1.7.1
#endif
#ifdef AMGCL_INCLUDE_IN_MY_PROJECT

doublereal* deltP_init = nullptr;

//#include "my_amgcl_alg_builtin.cpp" // Библиотека Дениса Демидова AMGCL.
#include "my_amgcl_alg.cpp" // Библиотека Дениса Демидова AMGCL. На видеокарте. Ограничена памятью видеокарты.
#include "my_amgcl_alg_builtin.cpp" // Библиотека Дениса Демидова AMGCL. На центральном процессоре. 
//#include "my_amgcl_alg_openMP.cpp"  // Библиотека Дениса Демидова AMGCL OpenMP.
#endif
// реализация алгебраического многосеточного метода 1985 года.
#include "amg1r5.cpp"
#include "my_agregat_amg.cpp"




/*
// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
void swap(NONZEROELEM* &v, integer i, integer j)
{
		NONZEROELEM temp;

		// change v[i] <-> v[j]
		temp = v[i];
		v[i] = v[j];
		v[j] = temp;
} // swap

// Вот алгоритм PivotList
integer PivotList(NONZEROELEM* &list, integer first, integer last) {
	// list обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue = list[first].key;
	integer PivotPointeger = first;

	for (integer index=(first+1); index<=last; index++) {
		if (list[index].key<PivotValue) {
			PivotPoint++;
			swap(list, PivotPoint, index);
		}
	}

	swap(list, first, PivotPoint);

	return PivotPoint;
} // PivotList


// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
void QuickSort(NONZEROELEM* &list, integer first, integer last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	integer pivot;

	if (first < last) {
		pivot = PivotList(list, first, last);
		QuickSort(list, first, pivot-1);
		QuickSort(list, pivot+1, last);
	}
} // QuickSort
*/
// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
void swap(integer*& v, integer i, integer j)
{
	integer temp;

	// change v[i] <-> v[j]
	temp = v[i];
	v[i] = v[j];
	v[j] = temp;
} // swap

// Вот алгоритм PivotList
integer PivotList(integer*& list, integer first, integer last) {
	// list обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue = list[first];
	integer PivotPoint = first;

	for (integer index = (first + 1); index <= last; index++) {
		if (list[index] < PivotValue) {
			PivotPoint++;
			swap(list, PivotPoint, index);
		}
	}

	swap(list, first, PivotPoint);

	return PivotPoint;
} // PivotList


// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
void QuickSort(integer*& list, integer first, integer last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	integer pivot;

	if (first < last) {
		pivot = PivotList(list, first, last);
		QuickSort(list, first, pivot - 1);
		QuickSort(list, pivot + 1, last);
	}
} // QuickSort

// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами v[i] и v[j]
template <typename MY_IND_TYPE>
void swapCSIR(MY_IND_TYPE*& v, doublereal*& dr, integer i, integer j)
{
	MY_IND_TYPE tempi;
	doublereal tempr;

	// change v[i] <-> v[j]
	tempi = v[i];
	v[i] = v[j];
	v[j] = tempi;
	// change dr[i] <-> dr[j]
	tempr = dr[i];
	dr[i] = dr[j];
	dr[j] = tempr;

} // swap

// Вот алгоритм PivotList
template <typename MY_IND_TYPE>
integer PivotListCSIR(MY_IND_TYPE*& jptr, doublereal*& altr, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue = jptr[first];
	integer PivotPoint = first;

	for (integer index = (first + 1); index <= last; index++) {
		if (jptr[index] < PivotValue) {
			PivotPoint++;
			swapCSIR(jptr, altr, PivotPoint, index);
		}
	}

	swapCSIR(jptr, altr, first, PivotPoint);

	return PivotPoint;
} // PivotListCSIR


// Быстрая сортировка Хоара.
// Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
// стр. 106.
template <typename MY_IND_TYPE>
void QuickSortCSIR(MY_IND_TYPE*& jptr, doublereal*& altr, integer first, integer last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	if (0) {
		// BubbleSort
		integer numberOfPairs = last - first + 1;
		bool swappedElements = true;
		while (swappedElements) {
			numberOfPairs--;
			swappedElements = false;
			for (integer i = first; i <= first + numberOfPairs - 1; ++i) {
				if (jptr[i] > jptr[i + 1]) {
					swapCSIR(jptr, altr, i, i + 1);
					swappedElements = true;
				}
			}
		}
	}
	else
	{
		integer pivot;

		if (first < last) {
			pivot = PivotListCSIR(jptr, altr, first, last);
			QuickSortCSIR(jptr, altr, first, pivot - 1);
			QuickSortCSIR(jptr, altr, pivot + 1, last);
		}
	}
} // QuickSortCSIR

/* Реализация на динамическом массиве.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CRS. Всего nodes - уравнений.
void simplesparsetoCRS(SIMPLESPARSE &M, doublereal* &val, integer* &col_ind, integer* &row_ptr, integer nodes) {
	if (M.n!=0) {
		val = new doublereal[M.n];
		col_ind = new integer[M.n];
		row_ptr = new integer[nodes+1];

		integer k; // счётчик
		// инициализация
		for (k=0; k<(M.n); ++k) {
		   val[k]=0.0;
		   col_ind[k]=0;
		}
		for (k=0; k<=nodes; ++k) {
			row_ptr[k]=M.n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}

		// Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		QuickSort(M.a, 0, M.n-1);

		// заполнение разреженной матрицы
		for (k=0; k<M.n; ++k) {
			val[k]=M.a[k].aij;
			col_ind[k]=M.a[k].j;
			row_ptr[M.a[k].i]=myi_min(k,row_ptr[M.a[k].i]);
		}
	}
} // simplesparsetoCRS
*/

// Пузырьковая сортировка одной строки. 10.04.2021
template <typename MY_IND_TYPE>
void BubbleSortCSIR(MY_IND_TYPE*& col_ind, doublereal*& val, integer start, integer end) {

	// BubbleSort
	integer numberOfPairs = end - start + 1;
	bool swappedElements = true;
	while (swappedElements) {
		numberOfPairs--;
		swappedElements = false;
		//for (integer i = start; i <= start + numberOfPairs - 1; ++i) 
		for (integer i = start; (i < end)/*&&(i <= start + numberOfPairs - 1)*/; ++i)
		{
			if (col_ind[i] > col_ind[i + 1]) {
				swapCSIR(col_ind, val, i, i + 1);
				swappedElements = true;
			}
		}
	}
}// BubbleSortCSIR


#if !bStableVersion

// 21.12.2021
// Реализация на связном списке
// Обнуляет внедиагональный коэффициент aij, Соответствующим образом корректирует диагональный элемент.
void set0elmsimplesparse_Stress_ell(integer i, integer j) {
	//const doublereal MY_ZERO_TOLERANCE = -1.0; // 1.0e-300; -1.0 Отключена проверка. Так нужно т.к. при условиях Дирихле для симметричности портрета нужны нули.

	doublereal vne_diagonal = 0.0;

	if (0) {

		// Могла использоваться в Solve_Thermal на ранних этапах. Не ставить 1.0 это ухудшает точность нахождения поля температур.
		// 22.12.2021 оставить выключенным.

		//NONZEROELEM* p;
		//p = M.root[i];// Корневой элемент строки i
		// линейный поиск элемента с ключом key
		//while ((p != nullptr) && (p->key != j)) p = p->next;

		// Чтобы сделать быструю обработку нужно для каждого i хранить АВЛ дерево ключей, чтобы быстро искать,
		// а чтобы быстро добавлять хранить размер - число добавленных индексов.
		// Данный линейный поиск занимает 6% времени работы алгоритма.
		int iscan = 0;
		while ((iscan < coll_size[i]) && (iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1) && (coll_ell[i][iscan] != j)) ++iscan;
		if (iscan == coll_size[i]) {
			std::cout << "data and coll ell limit index overflow in function 1 set0elmsimplesparse_Stress_ell\n";
			std::cout << "i==" << i << "size=" << coll_size[i] << std::endl;
			system("PAUSE");
		}
		if ((iscan < coll_size[i]) && (iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] == j)) {
			// элемент найден
			{
				// Удалить строку.
				vne_diagonal = fabs(data_ell[i][iscan]);
				data_ell[i][iscan] = 0.0;

				iscan = 0;
				while ((iscan < coll_size[i]) && (iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1) && (coll_ell[i][iscan] != i)) ++iscan;

				if (iscan == coll_size[i]) {
					std::cout << "data and coll ell limit index overflow in function 2 set0elmsimplesparse_Stress_ell\n";
					std::cout << "i==" << i << "size=" << coll_size[i] << std::endl;
					system("PAUSE");
				}

				if ((iscan < coll_size[i]) && (iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] == i)) {
					data_ell[i][iscan] -= vne_diagonal;
				}
			}

		}
	}
	
} // set0lmsimplesparse_ell


// Реализация на связном списке
 // Добавляет ненулевой элемент в
 // простейшую разряженную матрицу M
 // Проверки на равенство добавляемого элемента нулю нет, поэтому
 // может добавить и нулевой элемент.
void addelmsimplesparse_Stress_ell(doublereal aij, int i, int j, bool bset, bool bsetD) {
	const doublereal MY_ZERO_TOLERANCE = -1.0; // 1.0e-300; -1.0 Отключена проверка. Так нужно т.к. при условиях Дирихле для симметричности портрета нужны нули.

	//NONZEROELEM* p;
	//p = M.root[i];// Корневой элемент строки i
	// линейный поиск элемента с ключом key
	//while ((p != nullptr) && (p->key != j)) p = p->next;

	// Чтобы сделать быструю обработку нужно для каждого i хранить АВЛ дерево ключей, чтобы быстро искать,
	// а чтобы быстро добавлять хранить размер - число добавленных индексов.
	// Данный линейный поиск занимает 6% времени работы алгоритма.
	int iscan = 0;
	while ((iscan<coll_size[i])&&(iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1) && (coll_ell[i][iscan] != j)) ++iscan;
	if (iscan == coll_size[i]) {
		std::cout << "data and coll ell limit index overflow in function addelmsimplesparse_Stress_ell\n";
		std::cout << "i=="<< i << "size=" << coll_size[i] << std::endl;
		system("PAUSE");
	}
	if ((iscan < coll_size[i]) && (iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] == j)) {
		// элемент найден
		if (bsetD) {
			// Удалить строку.
			data_ell[i][iscan] = aij;
			if (iscan > 0) {
				std::cout << "Fatal error in addelmsimplesparse_Stress_ell\n";
				system("pause");
			}
		}
		else {
			if (fabs(aij) > MY_ZERO_TOLERANCE) {
				if (bset) data_ell[i][iscan] = aij; // установка
				else {
					data_ell[i][iscan] += aij; // добавление
				}
			}
		}
	}
	else
	{
		// если такого элемента нет в списке
		// то добавление элемента в начало списка.
		if (fabs(aij) > MY_ZERO_TOLERANCE) {
			if ((iscan < MAX_STRING_LENGTH_ELL)) {

				coll_ell[i][iscan] = j;
				data_ell[i][iscan] = aij;

#pragma omp atomic
				nnz_ell++;

			}
			else {
				std::cout << "ERROR!!! UVELICHTE MAX_STRING_LENGTH_ELL\n";
				system("pause");
				exit(1);
			}
		}
	}
} // addelmsimplesparse_ell

// Удаляет положительные внедиагональные элементы.
void patch_Ell_zero(int isize)
{
	for (int j = 0; j < isize; ++j) {
		int i = j;

		doublereal dsum = 0.0;

		int iscan = 0;
		int ics = 0;
		while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1)) {


			if (j != coll_ell[i][iscan]) {
				if (data_ell[i][iscan] > 0.0) {

					data_ell[i][iscan] = 0.0;

					/*
					for (int iscan3 = iscan; iscan3 < MAX_STRING_LENGTH_ELL - 1; ++iscan3) {
						if (coll_ell[i][iscan3] == -1) break;
						data_ell[i][iscan3] = data_ell[i][iscan3 + 1];
						coll_ell[i][iscan3] = coll_ell[i][iscan3 + 1];
					}
					data_ell[i][MAX_STRING_LENGTH_ELL - 1] = 0.0;
					coll_ell[i][MAX_STRING_LENGTH_ELL - 1] = -1;
					iscan--;
					*/
				}
				else {
					dsum += fabs(data_ell[i][iscan]);
					ics++;
				}
			}

			++iscan;
		}

		if (ics > 0) {
			data_ell[i][0] = dsum; // коррекция диагонали для соблюдения баланса.
		}

		/*iscan = 0;
		while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1)) {
			std::cout << "a[" << i << "," << coll_ell[i][iscan] << "]=" << data_ell[i][iscan] << "  ";
			iscan++;
		}
		std::cout << std::endl;
		system("pause");*/

	}
	//system("pause");
}


// Диагональ должна стоять первой в строке.
void patch_Ell(integer isize)
{
	for (integer j = 0; j < isize; ++j) {
		integer i = j;

		integer iscan = 0;
		while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1) && (coll_ell[i][iscan] != j)) ++iscan;
		if ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] == j)) {
			int itmp = coll_ell[i][0];
			coll_ell[i][0] = coll_ell[i][iscan];
			coll_ell[i][iscan] = itmp;

			doublereal dtmp = data_ell[i][0];
			data_ell[i][0] = data_ell[i][iscan];
			data_ell[i][iscan] = dtmp;
		}
		else {
			std::cout << "error patch_Ell: diagonal otsutstvuet. j==" << j << "  isize=" << isize << "\n";
			system("pause");
		}
		/*iscan = 0;
		while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[i][iscan] != -1)) {
			std::cout << "a[" << i << "," << coll_ell[i][iscan] << "]=" << data_ell[i][iscan] << "  ";
			iscan++;
		}
		std::cout << std::endl;
		system("pause");*/

	}
	//patch_Ell_zero(isize);
	//system("pause");
}



// При обнаружении сбоя включить stable версия данной функции.
// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CRS. Всего nodes - уравнений.
template <typename MY_IND_TYPE>
void ell_to_CRS(doublereal*& val, MY_IND_TYPE*& col_ind, MY_IND_TYPE*& row_ptr, integer nodes) {
	//bool flag = true;
	//for (MY_IND_TYPE k = 0; k < nodes; ++k) if (M.root[k] == nullptr) {
		//flag = false; break;
	//}

	//if (flag) {
	if ((val == nullptr) && (col_ind == nullptr) && (row_ptr == nullptr)) {
		val = new doublereal[nnz_ell];
		col_ind = new MY_IND_TYPE[nnz_ell];
		row_ptr = new MY_IND_TYPE[nodes + 1];
	}
	else {
		if (val != nullptr) {
			std::cout << "error allocate ell_to_CRS. for val...\n";
			system("PAUSE");
		}
		if (col_ind != nullptr) {
			std::cout << "error allocate ell_to_CRS. for col_ind...\n";
			system("PAUSE");
		}
		if (row_ptr != nullptr) {
			std::cout << "error allocate ell_to_CRS. for row_ptr...\n";
			system("PAUSE");
		}
	}

	/*bool* bcheck = new bool[M.n];
	for (MY_IND_TYPE i_1 = 0; i_1 < M.n; ++i_1) {
		bcheck[i_1] = false;
	}
	NONZEROELEM* p_1 = nullptr;
	for (MY_IND_TYPE k = 0; k < nodes; ++k) {
		p_1 = M.root[k];
		while (p_1 != nullptr) {
			if (bcheck[p_1->key]) {
				printf("ERROR MATRIX CHECK duplicate ja index string=%lld col_ind=%lld\n", k, p_1->key);
				system("pause");
			}
			bcheck[p_1->key] = true;

			p_1 = p_1->next;
		}
		p_1 = M.root[k];
		// Сброс.
		while (p_1 != nullptr) {
			bcheck[p_1->key] = false;
			p_1 = p_1->next;
		}
	}
	delete[] bcheck;*/

	// инициализация
#pragma omp parallel for
	for (int k = 0; k < (nnz_ell); ++k) {
		val[k] = 0.0;
		col_ind[k] = 0;
	}

#pragma omp parallel for
	for (integer k = 0; k <= nodes; ++k) {
		row_ptr[k] = static_cast<MY_IND_TYPE>(nnz_ell); // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	}

	// Быстрая Сортировка Хоара.
	// упорядочивание по строкам
	//QuickSort(...); не требуется,
	// т.к. сама структура хранения 
	// подразумевает упорядочивание по строкам.

	/*
	// заполнение разреженной матрицы
	for (MY_IND_TYPE k=0; k<M.n; ++k) {
		val[k]=M.a[k].aij;
		col_ind[k]=M.a[k].j;
		row_ptr[M.a[k].i]=myi_min(k,row_ptr[M.a[k].i]);
	}
	*/


	MY_IND_TYPE ik = 0; // счётчик ненулевых элементов СЛАУ
	//NONZEROELEM* p;

	MY_IND_TYPE nnz = 0;
	for (MY_IND_TYPE k = 0; k < nodes; ++k) {
		int iscan = 0;
		//p = M.root[k];
		while ((iscan < MAX_STRING_LENGTH_ELL)&&(iscan<coll_size[k]) && (coll_ell[k][iscan] != -1)) {
			//if (ik < M.n) {
				// Защита от записи нулевого коэффициента.
			//if (fabs(data_ell[k][iscan]) > 1.0e-300) При граничных условиях Дирихле для симметричности портрета нужно некоторое небольшое количество нулей.
			{
				nnz++;// 20.03.2021
				val[ik] = data_ell[k][iscan];
				/*if (p->key == k) {

					printf("diag=%e string %lld \n", val[ik], k);
					system("pause");
				}*/
				// p->key начинается с нуля
				col_ind[ik] = coll_ell[k][iscan];
				/*if (p->key < 0) {
					printf("%lld\n", row_ptr[k]);
					system("pause");
				}*/
				row_ptr[k] = std::min(ik, row_ptr[k]);
				ik++;
			}
			/*else {
				if (col_ind[ik] == k) {

					printf("diag=%e string %lld \n", val[ik], k);
					system("pause");
				}
			}*/
			//}
			/*else {
				printf("error: ik>=M.n in simplesparsetoCRS\n");
				printf("see module my_linalg.c\n");
				system("pause");
				exit(1);
			}*/
			//p = p->next;
			iscan++;
		}
	
		if (iscan == coll_size[k]) {
			std::cout << "data and coll ell limit index overflow in ell_to_CRS function\n";
			system("PAUSE");
		}
	}

	// Исправлена ошибка. Истинное количество ненулевых элементов меньше из за
	// того что берутся лишь ненулевые элементы выше некоторого порога.
	if (yes_print_amg) {
		if (my_amg_manager.iprint_log != 0) {
			std::cout << "reserve: M.n =" << nnz_ell << ", natural count : nnz = " << nnz << "\n";
		}
	}
	//system("pause"); // 20.03.2021
	row_ptr[nodes] = nnz;

	// в каждой строке элементы отсортированы по номерам столбцов:
	for (MY_IND_TYPE k = 0; k < nodes; ++k) {


		if (0) {
			// BubbleSort
			BubbleSortCSIR(col_ind, val, row_ptr[k], row_ptr[k + 1] - 1);
		}
		else {
			QuickSortCSIR(col_ind, val, row_ptr[k], row_ptr[k + 1] - 1);
		}

	}

	//std::cout << "Sortirovka proidena...\n";
	//system("pause");

	//}
} // ell_to_CRS

// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CSIR. Всего nodes - уравнений.
// Это работает только для SPD матриц.
// Симметричный положительно определённый случай,
// хранится только нижний треугольник.
void  ell_to_CSIR(doublereal*& adiag, doublereal*& altr, integer*& jptr, integer*& iptr, integer nodes) {
	bool flag = true;
	integer k; // счётчик
	for (k = 0; k < nodes; ++k) if (coll_ell[k][0] == -1) {

		std::cout << "error ell_to_CSIR:  k=" << k << "  nodes=" << nodes << std::endl;
		system("pause");

		std::cout << coll_ell[k][0] << " " << coll_ell[k][1] << " " << coll_ell[k][2] << " " << coll_ell[k][3] << " " << coll_ell[k][4] << " " << coll_ell[k][5] << " \n";

		flag = false;
		//break;
	}

	if (flag) {
		// поддиагональные элементы в altr хранятся построчно
		integer nz = static_cast<integer>((nnz_ell - nodes) / 2); // число ненулевых элементов
		adiag = new doublereal[nodes]; // диагональные элементы
		altr = new doublereal[nz]; // поддиагональные элементы
		jptr = new integer[nz]; // номера столцов для нижнего треугольника
		iptr = new integer[nodes + 1]; // указатели на следующую строку


		// инициализация
		for (k = 0; k < nodes; ++k) adiag[k] = 0.0;
		for (k = 0; k < (nz); ++k) {
			altr[k] = 0.0;
			jptr[k] = 0;
		}
		for (k = 0; k <= nodes; ++k) {
			iptr[k] = nz; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}

		// Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (k=0; k<M.n; ++k) {
			val[k]=M.a[k].aij;
			col_ind[k]=M.a[k].j;
			row_ptr[M.a[k].i]= myi_min(k,row_ptr[M.a[k].i]);
		}
		*/
		/*
		integer ik=0; // счётчик ненулевых элементов СЛАУ
		NONZEROELEM* p;
		for (k=0; k<nodes; ++k) {
			p=M.root[k];
			while (p!=nullptr) {
				val[ik]=p->aij;
				col_ind[ik]=p->key;
				row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
				p=p->next;
			}
		}
		*/

		integer ik = 0, imin = 1, k1; // счётчик ненулевых поддиагональных элементов СЛАУ
		bool bvisit;
		//NONZEROELEM* p;
		for (k = 0; k < nodes; ++k) {
			bvisit = false;
			//p = M.root[k];
			int iscan = 0;
			while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[k][iscan] != -1)) {
				if (coll_ell[k][iscan] == k) {
					adiag[k] = data_ell[k][iscan];
					bvisit = false;
				}
				else if (coll_ell[k][iscan] < k) {
					if (ik < (nz)) {
						altr[ik] = data_ell[k][iscan]; // ненулевое значение
						jptr[ik] = coll_ell[k][iscan]; // номер столбца
					}
					else {
						printf("non simmetric matrix ICCG. ell_to_CSIR\n");
						std::cout << "ik=" << ik << " nz=" << nz << " coll_ell[" << k << "][iscan]=" << coll_ell[k][iscan] << "n==" << nodes << std::endl;
						//system("pause");
						system("pause");
					}
					bvisit = true;
				}
				imin = myi_min(ik, iptr[k]);
#if doubleintprecision == 1
				//printf("imin=%lld\n",imin);
#else
				//printf("imin=%d\n",imin);
#endif

				iptr[k] = imin;
				if (imin == 0) for (k1 = 0; k1 < k; k1++) iptr[k1] = 0;
				if (bvisit) {
					ik++;
					bvisit = false;
				}
				//p = p->next;
				++iscan;
			}
			if (iscan >= MAX_STRING_LENGTH_ELL) {
				std::cout << "ERROR!!! UVELICHTE MAX_STRING_LENGTH_ELL\n";
				system("pause");
				exit(1);
			}
		}


		for (k = 0; k < nodes; ++k) QuickSortCSIR(jptr, altr, iptr[k], iptr[k + 1] - 1);

	}
} // ell_to_CSIR


// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CSIR_ITL. Всего nodes - уравнений.
// Это работает только для SPD матриц.
// Симметричный положительно определённый случай,
// хранится только верхний треугольник.
// Память выделяется внутри метода.
void ell_to_CSIR_ITLSPD(doublereal*& val, integer*& indx, integer*& pntr, integer nodes) {
	bool flag = true;
	integer k; // счётчик
	for (k = 0; k < nodes; ++k) if (coll_ell[k][0] == -1) {
		flag = false;

		std::cout << "error ell_to_CSIR_ITLSPD:  k=" << k << "  nodes=" << nodes << std::endl;
		system("pause");
		//break;
	}

	if (flag) {

		//printM_and_CSIR(M, nodes); // debug

		// поддиагональные элементы в altr хранятся построчно
		integer nz = static_cast<integer>((nnz_ell - nodes) / 2 + nodes); // число ненулевых элементов
		val = new doublereal[nz]; // диагональные элементы и наддиагональные элементы
		indx = new integer[nz]; // номера столцов для нижнего треугольника
		pntr = new integer[nodes + 1]; // указатели на следующую строку


		// инициализация
		for (k = 0; k < (nz); ++k) {
			val[k] = 0.0;
			indx[k] = 0;
		}
		for (k = 0; k <= nodes; ++k) {
			pntr[k] = nz; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}



		integer ik = 0; // счётчик ненулевых поддиагональных элементов СЛАУ
		//NONZEROELEM* p;
		for (k = 0; k < nodes; ++k) {

			//p = M.root[k];
			int iscan = 0;
			while ((iscan < MAX_STRING_LENGTH_ELL) && (coll_ell[k][iscan] != -1)) {

				// k - номер диагонального элемента
				if ((coll_ell[k][iscan] >= k) && (coll_ell[k][iscan] < nodes)) {
					if (ik < (nz)) {
						val[ik] = data_ell[k][iscan]; // ненулевое значение
						indx[ik] = coll_ell[k][iscan]; // номер столбца	
					}
					else {
						std::cout << "\nik=" << ik << " nz=" << nz << " nodes=" << nodes << " k=" << k << std::endl;
						printf(" Error non simmetric matrix ICCG. ell_to_CSIR_ITLSPD\n");
						//system("pause");
						system("pause");
					}
					pntr[k] = myi_min(ik, pntr[k]);

					ik++;
				}

				//p = p->next;
				iscan++;
			}

			if (iscan >= MAX_STRING_LENGTH_ELL) {
				std::cout << "ERROR!!! UVELICHTE MAX_STRING_LENGTH_ELL\n";
				system("pause");
				exit(1);
			}

		}

		for (k = 0; k < nodes; ++k) QuickSortCSIR(indx, val, pntr[k], pntr[k + 1] - 1);


		/*
		FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, "matr.txt", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
	#if doubleintprecision == 1
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");
		// debug
		for (k=0; k<=nodes; ++k) {
			fprintf(fp,"%lld ",pntr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<pntr[nodes]; ++k) {
			fprintf(fp, "%e %lld\n",val[k],indx[k]);
		}
		fprintf(fp, "nz==%lld\n", nz);
	#else
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");
		// debug
		for (k=0; k<=nodes; ++k) {
			fprintf(fp,"%d ",pntr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<pntr[nodes]; ++k) {
			fprintf(fp, "%e %d\n",val[k],indx[k]);
		}
		fprintf(fp, "nz==%d\n", nz);
	#endif


		fclose(fp);
	}
	printf("ready");
	system("pause");
	*/
	}
} // ell_to_CSIR_ITLSPD

#endif


doublereal rterminate_residual_ICCG_Oh2(FLOW floc) {
	// точность дискретизации O(h!2)
	doublereal* resterm = new doublereal[floc.maxelm + floc.maxbound];
	for (integer i = 0; i < floc.maxelm + floc.maxbound; ++i) {
		resterm[i] = 0.0; // инициализация.
	}

	for (integer iP = 0; iP < floc.maxelm; iP++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
		volume3D(iP, floc.nvtx, floc.pa, dx, dy, dz);
		doublereal dl = fmin(dx, fmin(dy, dz));
		resterm[iP] = 0.1*dl*dl; // O(h!2)
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
		iE = floc.neighbors_for_the_internal_node[E_SIDE][0][iP]; 
		iN = floc.neighbors_for_the_internal_node[N_SIDE][0][iP];
		iT = floc.neighbors_for_the_internal_node[T_SIDE][0][iP]; 
		iW = floc.neighbors_for_the_internal_node[W_SIDE][0][iP];
		iS = floc.neighbors_for_the_internal_node[S_SIDE][0][iP]; 
		iB = floc.neighbors_for_the_internal_node[B_SIDE][0][iP];
		// Если с одной из сторон стоит граница расчётной области
		// то соответствующая переменная равна true
		bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

		if (iE >= floc.maxelm) bE = true;
		if (iN >= floc.maxelm) bN = true;
		if (iT >= floc.maxelm) bT = true;
		if (iW >= floc.maxelm) bW = true;
		if (iS >= floc.maxelm) bS = true;
		if (iB >= floc.maxelm) bB = true;

		if ((bE) || (bW)) {
			dl = 0.5*dx;
			if (bE) resterm[iE] = 0.1*dl*dl; // O(h!2)
			if (bW) resterm[iW] = 0.1*dl*dl; // O(h!2)
		}
		if ((bN) || (bS)) {
			dl = 0.5*dy;
			if (bN) resterm[iN] = 0.1*dl*dl; // O(h!2)
			if (bS) resterm[iS] = 0.1*dl*dl; // O(h!2)
		}
		if ((bT) || (bB)) {
			dl = 0.5*dz;
			if (bT) resterm[iT] = 0.1*dl*dl; // O(h!2)
			if (bB) resterm[iB] = 0.1*dl*dl; // O(h!2)
		}
	}
	doublereal ret = Scal(resterm, resterm, floc.maxelm + floc.maxbound);
	delete[] resterm;
	resterm = nullptr;
	return ret;
} // rterminate_residual_ICCG_Oh2

doublereal rterminate_residual_LR1sk_Oh3(FLOW floc) {
	// точность дискретизации O(h!2)
	doublereal* resterm = new doublereal[floc.maxelm + floc.maxbound];
	for (integer i = 0; i < floc.maxelm + floc.maxbound; ++i) {
		resterm[i] = 0.0; // инициализация.
	}

	for (integer iP = 0; iP < floc.maxelm; iP++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
		volume3D(iP, floc.nvtx, floc.pa, dx, dy, dz);
		doublereal dl = fmin(dx, fmin(dy, dz));
		resterm[iP] = 0.1*dl*dl*dl; // O(h!3)
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
		iE = floc.neighbors_for_the_internal_node[E_SIDE][0][iP];
		iN = floc.neighbors_for_the_internal_node[N_SIDE][0][iP]; 
		iT = floc.neighbors_for_the_internal_node[T_SIDE][0][iP]; 
		iW = floc.neighbors_for_the_internal_node[W_SIDE][0][iP];
		iS = floc.neighbors_for_the_internal_node[S_SIDE][0][iP]; 
		iB = floc.neighbors_for_the_internal_node[B_SIDE][0][iP];
		// Если с одной из сторон стоит граница расчётной области
		// то соответствующая переменная равна true
		bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

		if (iE >= floc.maxelm) bE = true;
		if (iN >= floc.maxelm) bN = true;
		if (iT >= floc.maxelm) bT = true;
		if (iW >= floc.maxelm) bW = true;
		if (iS >= floc.maxelm) bS = true;
		if (iB >= floc.maxelm) bB = true;

		if ((bE) || (bW)) {
			dl = 0.5*dx;
			if (bE) resterm[iE] = 0.1*dl*dl*dl; // O(h!3)
			if (bW) resterm[iW] = 0.1*dl*dl*dl; // O(h!3)
		}
		if ((bN) || (bS)) {
			dl = 0.5*dy;
			if (bN) resterm[iN] = 0.1*dl*dl*dl; // O(h!3)
			if (bS) resterm[iS] = 0.1*dl*dl*dl; // O(h!3)
		}
		if ((bT) || (bB)) {
			dl = 0.5*dz;
			if (bT) resterm[iT] = 0.1*dl*dl*dl; // O(h!3)
			if (bB) resterm[iB] = 0.1*dl*dl*dl; // O(h!3)
		}
	}
	doublereal ret;
	//ret=Scal(resterm,resterm,floc.maxelm+floc.maxbound);
	ret = NormaV(resterm, floc.maxelm + floc.maxbound);
	// Освобождение оперативной памяти.
	if (resterm != nullptr) {
		delete[] resterm;
		resterm = nullptr;
	}
	return ret;
} // rterminate_residual_LR1sk_Oh3

doublereal rterminate_residual_LR1sk_temp_Oh3(TEMPER t) {
	// точность дискретизации O(h!2)
	// точность дискретизации теплопередачи в твёрдом теле O(h).
	doublereal* resterm = new doublereal[t.maxelm + t.maxbound];
	for (integer i = 0; i < t.maxelm + t.maxbound; ++i) {
		resterm[i] = 0.0; // инициализация.
	}

	for (integer iP = 0; iP < t.maxelm; iP++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
		volume3D(iP, t.nvtx, t.pa, dx, dy, dz);
		doublereal dl = fmin(dx, fmin(dy, dz));
		//resterm[iP]=0.1*dl*dl*dl; // O(h!3)
		resterm[iP] = dl; // O(h)
		integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
		iE = t.neighbors_for_the_internal_node[E_SIDE][0][iP]; 
		iN = t.neighbors_for_the_internal_node[N_SIDE][0][iP]; 
		iT = t.neighbors_for_the_internal_node[T_SIDE][0][iP]; 
		iW = t.neighbors_for_the_internal_node[W_SIDE][0][iP]; 
		iS = t.neighbors_for_the_internal_node[S_SIDE][0][iP];
		iB = t.neighbors_for_the_internal_node[B_SIDE][0][iP];
		// Если с одной из сторон стоит граница расчётной области
		// то соответствующая переменная равна true
		bool bE = false, bN = false, bT = false, bW = false, bS = false, bB = false;

		if (iE >= t.maxelm) bE = true;
		if (iN >= t.maxelm) bN = true;
		if (iT >= t.maxelm) bT = true;
		if (iW >= t.maxelm) bW = true;
		if (iS >= t.maxelm) bS = true;
		if (iB >= t.maxelm) bB = true;

		if ((bE) || (bW)) {
			dl = 0.5*dx;
			//if (bE) resterm[iE]=0.1*dl*dl*dl; // O(h!3)
			//if (bW) resterm[iW]=0.1*dl*dl*dl; // O(h!3)
			if (bE) resterm[iE] = dl; // O(h)
			if (bW) resterm[iW] = dl; // O(h)

		}
		if ((bN) || (bS)) {
			dl = 0.5*dy;
			//if (bN) resterm[iN]=0.1*dl*dl*dl; // O(h!3)
			//if (bS) resterm[iS]=0.1*dl*dl*dl; // O(h!3)
			if (bN) resterm[iN] = dl; // O(h)
			if (bS) resterm[iS] = dl; // O(h)

		}
		if ((bT) || (bB)) {
			dl = 0.5*dz;
			//if (bT) resterm[iT]=0.1*dl*dl*dl; // O(h!3)
			//if (bB) resterm[iB]=0.1*dl*dl*dl; // O(h!3)
			if (bT) resterm[iT] = dl; // O(h)
			if (bB) resterm[iB] = dl; // O(h)
		}
	}
	doublereal ret;
	//ret=Scal(resterm,resterm,f.maxelm+f.maxbound);
	ret = NormaV(resterm, t.maxelm + t.maxbound);
	// Освобождение оперативной памяти.
	if (resterm != nullptr) {
		delete[] resterm;
		resterm = nullptr;
	}
	return ret;
} // rterminate_residual_LR1sk_temp_Oh3

// модель вещественной арифметики определяется в главном модуле программы 
// AliceFlow_v0_27.cpp в самом начале программы.
//#define doublereal double // модель веществекнного числа

void isfinite_vec(integer n, doublereal* xtest, const char* sname) {
	for (integer i=0; i<n; ++i) {
		if (xtest[i]!= xtest[i]) {
#if doubleintprecision == 1
			printf(" problem infinity in vector %s in position %lld. size vector=%lld\n", sname, i, n);
#else
			printf(" problem infinity in vector %s in position %d. size vector=%d\n", sname, i, n);
#endif
			
		//	system("pause");
			system("pause");
		}
	}
}


//const doublereal dterminatedTResudual = 1e-40; // для МСГ Congruate Gradients

// данная параллельная версия метода гаусса не даёт отличий от серийной версии
// Это было проверено массовым тестированием (см. код eqsolve_simple_gauss).
// Данная версия правильная.
// Существует еще одна версия параллельной программы. В книге Анализ алгоритмов вводный курс
// есть параллельный алгоритмом SIMD. Он пригоден для очень большого числа потока и возможно
// более хорош чем данный код. 
void my_version_gauss(doublereal **A, integer nodes, doublereal *b, doublereal *x, bool bparallel) {

	// если bparallel==true то параллельный режим работы иначе серийный.
	// можно использовать как то так и другой код они оба правильные.

  
   const doublereal epsilon = 1e-100;
 

//#ifdef _OPENMP
   //omp_set_num_threads(inumcore); // установка числа потоков
//#endif

   // приведение к треугольному виду:
   for(integer k=0; k<nodes; ++k){
	   doublereal akk=A[k][k];
       if(fabs(akk)<epsilon){
		  // решение не может быть получено, т.к.
		  // на диагонали находится ноль.
	      printf("\nSolution is not exist! Gauss divizion by zero...\n");
	    //  system("pause");
		  system("pause");
	      exit(0);
	   }

	   if (bparallel) {
          #pragma omp parallel for shared(k, nodes, A, b)  firstprivate(akk)
          for(integer i=k+1; i<nodes; ++i) {
	      
	          doublereal M = A[i][k] / akk;
	          for(integer j=k; j<nodes; ++j){
	             A[i][j] -= M * A[k][j];
	          }
	          b[i] -= M*b[k];
          }
	   }
	   else {
          for(integer i=k+1; i<nodes; ++i) {
	      
	          doublereal M = A[i][k] / akk;
	          for(integer j=k; j<nodes; ++j){
	             A[i][j] -= M * A[k][j];
	          }
	          b[i] -= M*b[k];
          }
	   }
   }
   // процесс обратного исключения
   x[nodes-1]=b[nodes-1]/A[nodes-1][nodes-1];
   for(integer i=nodes-2; i>=0; --i){
       doublereal sum = 0.0;
	   if (bparallel) {
          #pragma omp parallel for shared(A,x,i,nodes)  reduction (+: sum)
          for(integer j = i+1; j<nodes; ++j){
	          sum+= A[i][j]*x[j];
          }
	   }
	   else {
          for(integer j = i+1; j<nodes; ++j){
	          sum+= A[i][j]*x[j];
          }
	   }
       x[i] = (b[i] - sum) / A[i][i];
   }
} // my_version_gauss

// Для тестирования метода исключения Гаусса.
// для целей тестирования создаёт клон матрицы СЛАУ А, правой части b, и вектора x.
void alloc_duplicate(doublereal ** &A, integer nodes, doublereal * &b,
	doublereal * &x,doublereal ** &A_copy, doublereal * &b_copy, doublereal * &x_copy) {
	A_copy=new doublereal*[nodes];
	for (integer i=0; i<nodes; ++i) {
		A_copy[i]=new doublereal[nodes];
	}
	b_copy=new doublereal[nodes];
	x_copy=new doublereal[nodes];
	for (integer i=0; i<nodes; ++i) {
		for (integer j=0; j<nodes; ++j) {
			A_copy[i][j]=A[i][j];
		}
		b_copy[i]=b[i];
		x_copy[i]=x[i];
	}
} // alloc duplicate.

// Для тестирования, освобождение памяти при тестировании метода Гаусса.
void free_duplicate(doublereal ** &A, integer nodes, doublereal * &b, doublereal * &x) {
	for (integer i=0; i<nodes; ++i) {
		delete A[i];
	}
	delete A;
	delete x;
	delete b;
}

/*  Решает систему уравнений для квадратной
 *  несимметричной матрицы коэффициентов A
 *        A*x=b;
 *  где A размеров nodes*nodes. Нумерация 
 *  элементов везде начинается с нуля.
 *  Процедура представляет собой метод Гаусса
 *  без выбора главного элемента и без учёта
 *  разреженности матрицы.
 *  A и b не сохраняются. 
*/
void eqsolve_simple_gauss(doublereal **A, integer nodes, doublereal *b, doublereal *x) {

   /*
   // Тестирование и верификация параллельного метода Гаусса.
   doublereal **A_copy;
   doublereal* b_copy;
   doublereal* x_copy;
   alloc_duplicate(A, nodes, b, x, A_copy, b_copy, x_copy);
   my_version_gauss(A_copy, nodes, b_copy, x_copy, false);

   my_version_gauss(A, nodes, b, x, true);
   doublereal *err=new doublereal[nodes];
   for (integer i=0; i<nodes; ++i) {
	   err[i]=fabs(x[i]-x_copy[i]);
   }
   doublereal er=NormaV(err,nodes);
   if (fabs(er)>1.0e-30) {
      printf("Error is: %e\n",er);
	  system("pause");
   }
   
   free_duplicate(A_copy, nodes, b_copy, x_copy);
   delete err;
   */
   bool bparallel=true;
   my_version_gauss(A, nodes, b, x, bparallel);

  
} // eqsolve_simple_gauss

/*  Решает систему уравнений для квадратной
 *  симметричной положительно определённой
 *  (с диагональным преобладанием) матрицы
 *  коэффициентов А:
 *        A*x=b;
 *  где A размеров nodes*nodes. Матрица А
 *  предполагается не разреженной. Нумерация 
 *  элементов везде начинается с нуля.
 *  Процедура представляет собой разложение Холесского:
 *        A=L*transpose(L),
 *  после которого выполняются прямое исключение и 
 *  обратная подстановка. A и b не сохраняются. 
*/
void eqsolv_simple_holesskii(doublereal **A, integer nodes, doublereal *b, doublereal *x) {
	// Разложение Холесского: замена A верхним и нижним 
	// треугольными множителями.
#if doubleprecision == 1 
		A[0][0] = sqrt(A[0][0]);
		A[1][0] /= A[0][0];
		A[0][1] = A[1][0];
		A[1][1] = sqrt(A[1][1] - A[1][0] * A[1][0]);
	
#else 
		A[0][0] = sqrtf(A[0][0]);
		A[1][0] /= A[0][0];
		A[0][1] = A[1][0];
		A[1][1] = sqrtf(A[1][1] - A[1][0] * A[1][0]);
#endif

//#ifdef _OPENMP
	//omp_set_num_threads(inumcore); // установка числа потоков
//#endif


	
	for (integer irow=2; irow<nodes; irow++) {
		integer irow1=irow-1;
		A[irow][0]/=A[0][0];
        A[0][irow]=A[irow][0];
        #pragma omp parallel for shared(irow1,A) 
		for (integer icol=1; icol<=irow1; icol++) {
			integer icol1=icol-1;
			doublereal  sum=0.0;
            for (integer k=0; k<=icol1; ++k) sum+=A[irow][k]*A[icol][k];
			A[irow][icol]=(A[irow][icol]-sum)/A[icol][icol];
			A[icol][irow]=A[irow][icol];
		}
		doublereal sum=0.0;
		#pragma omp parallel for shared(A,irow,irow1)  reduction (+: sum)
		for (integer k=0; k<=irow1; ++k) sum+=A[irow][k]*A[irow][k];
#if doubleprecision == 1 
			A[irow][irow] = sqrt(A[irow][irow] - sum);
		
#else 
			A[irow][irow] = sqrtf(A[irow][irow] - sum);
#endif
	}
    
	// Прямое исключение. Происходит разрушение правой части
	b[0]/=A[0][0];

	for (integer irow=1; irow<nodes; ++irow) {
		integer irow1=irow-1;
		doublereal sum=0.0;
		#pragma omp parallel for shared(A,b,irow,irow1) reduction (+: sum)
		for (integer icol=0; icol<=irow1; ++icol) sum+=A[irow][icol]*b[icol];
        b[irow]=(b[irow]-sum)/A[irow][irow];
	}

	// Обратная подстановка используется верхний треугольный множитель
	x[nodes-1]=b[nodes-1]/A[nodes-1][nodes-1];
	for (integer k=1; k<=nodes; ++k) {
		integer irow=nodes+1-k-1;
		integer irow1=irow+1;
		doublereal sum=0.0;
        #pragma omp parallel for shared(A,x,irow,irow1,nodes) reduction (+: sum)
		for (integer icol=irow1; icol<nodes; ++icol) sum+=A[irow][icol]*x[icol];
		x[irow]=(b[irow]-sum)/A[irow][irow];
	}

} // eqsolv_simple_holesskii

/* Находит обратную матрицу для 
*  квадратной матрицы A nodes*nodes с 
*  ненулевыми элементами на главной диагонали.
*  Решение производится путём метода исключения
*  Гаусса, а именно решая nodes СЛАУ. 
*          A*inv=e
*  Приведение  к треугольному виду делается
*  только один раз.
* Если flag==true, то матрица уже приведена к верхнетреугольному виду.
*/
void inverse_matrix_simple(doublereal** &A, integer nodes, bool flag) {

    const doublereal epsilon = 1e-100;

	doublereal **e=nullptr; // единичная матрица правых частей.
	doublereal **inv=nullptr; // будущая обратная матрица

	doublereal **acopy=nullptr; // копия матрицы А
	
	acopy = new doublereal* [nodes];
    for (integer i1=0; i1<nodes; i1++) acopy[i1]=new doublereal[nodes]; 
	

	if (acopy != nullptr) {

		integer i1 = 0, j1 = 0, k1 = 0;
		e = new doublereal*[nodes];
		for (i1 = 0; i1 < nodes; i1++) e[i1] = new doublereal[nodes];
		inv = new doublereal*[nodes];
		for (i1 = 0; i1 < nodes; i1++) inv[i1] = new doublereal[nodes];

		// инициализация
		for (i1 = 0; i1 < nodes; i1++) for (j1 = 0; j1 < nodes; j1++) {
			inv[i1][j1] = 0.0; // обратная матрица
			e[i1][j1] = 0.0; // правые части
			acopy[i1][j1] = A[i1][j1];
		}
		for (i1 = 0; i1 < nodes; i1++) e[i1][i1] = 1.0;




		if (!flag) { // если матрица ещё не приведена к верхнетреугольному виду
			doublereal M;
			// приведение к верхне треугольному виду:
			for (k1 = 0; k1 < nodes; k1++) {
				for (i1 = k1 + 1; i1 < nodes; i1++) {
					// Если на диагонали ноль:
					if (fabs(A[k1][k1]) < epsilon) {
						// решение не может быть получено, т.к.
						// на диагонали находится ноль.
						printf("\n inverse matrix simple ERROR !!! may be diagonal value is zero...\n");
						printf("\nSolution is not exist.\n");
						for (integer irow = 0; irow < nodes; irow++) {
							for (integer icol = 0; icol < nodes; icol++) {
								printf("%1.4e ", acopy[irow][icol]);
								//printf("%1.4e ", A[irow][icol]);
								
							}
							printf("\n");
						}
						//system("pause");
						system("pause");
						exit(0);
					}
					M = A[i1][k1] / A[k1][k1];
					for (j1 = k1; j1 < nodes; j1++) {
						A[i1][j1] -= M * A[k1][j1];
					}
					// преобразование правых частей:
					for (j1 = 0; j1 < nodes; j1++) e[i1][j1] -= M*e[k1][j1];
				}
			}
		}
		doublereal *sum = nullptr;
		sum = new doublereal[nodes];

		// процесс обратного исключения
		for (i1 = nodes - 1; i1 >= 0; i1--) {
			// инициализация
			for (k1 = 0; k1 < nodes; k1++) sum[k1] = 0.0;

			for (j1 = i1 + 1; j1 < nodes; j1++) {
				for (k1 = 0; k1 < nodes; k1++) {
					sum[k1] += A[i1][j1] * inv[j1][k1];
				}
			}
			for (k1 = 0; k1 < nodes; k1++) {
				inv[i1][k1] = (e[i1][k1] - sum[k1]) / A[i1][i1];
			}
		}
		if (e != nullptr) {
			for (i1 = nodes - 1; i1 >= 0; i1--) {
				if (e[i1] != nullptr) {
					delete[] e[i1];
				}
			}
			delete[] e;
			e = nullptr;
		}

		

		for (k1 = 0; k1 < nodes; k1++) {
			for (i1 = 0; i1 < nodes; i1++) {
				A[k1][i1] = inv[k1][i1];
			}
		}
		if (inv != nullptr) {
			for (i1 = nodes - 1; i1 >= 0; i1--) {
				if (inv[i1] != nullptr) {
					delete[] inv[i1];
				}
			}
			delete[] inv;
			inv = nullptr;
		}
		
		if (acopy != nullptr) {
			for (i1 = nodes - 1; i1 >= 0; i1--) {
				if (acopy[i1] != nullptr) {
					delete[] acopy[i1];
				}
			}
			delete[] acopy;
			acopy = nullptr;
		}
		
		if (sum != nullptr) {
			delete[] sum;
		}
	}
	else {
		printf("no allocate memory for acopy\n");
		printf("in function inverse_matrix_simple\n");
		system("pause");
		exit(1);
	}
} // inverse_matrix_simple
 
/* Находит произведение двух квадратных
* матриц A и B размерами nodes*nodes 
*             C=A*B. 
* Результат  записывается в матрицу B.
*/
void multiply_matrix_simple(doublereal **A1, doublereal **A2, integer nodes) {
	integer i1=0, j1=0, k1=0; // счётчики цикла for
	
	doublereal **c=nullptr;
	c = new doublereal* [nodes];
    for (i1=0; i1<nodes; i1++) c[i1]=new doublereal[nodes];

	for (i1=0; i1<nodes; i1++) for (j1=0; j1<nodes; j1++) c[i1][j1]=0.0; // инициализация

	// умножение C=A1*A2:
    for (i1=0; i1 < nodes; i1++)
        for (k1=0; k1 < nodes; k1++)
            for (j1=0; j1 < nodes; j1++)
                c[i1][k1]+=(A1[i1][j1])*(A2[j1][k1]);

	// копирование результата в A2:
    for (i1=0; i1<nodes; i1++) for (j1=0; j1<nodes; j1++) A2[i1][j1]=c[i1][j1];

	if (c != nullptr) {
		for (i1 = 0; i1 < nodes; i1++) 
			if (c[i1] != nullptr) {
				delete[] c[i1];
			}
		delete[] c;
	}
} // multiply_matrix_simple


// Следующие несколько функций (шесть, но теперь
// транспонирование не используется, а умножение заменено более быстрым). 
// используются как вспомогательные для решения
// полной проблемы собственных значений:

/* 1. Умножение квадратных матриц размера n*n:
*                t=m*p.
* Нумерация начинается с нуля.
* По окончании работы результат хранится в матрице t.
*/
void multi_m(doublereal **m, doublereal **p, doublereal **t, integer n) {
    for (integer i = 0; i < n; ++i)
       for (integer j = 0; j < n; ++j) {
           doublereal s = 0;
           for (integer l = 0; l < n; l++)
               s += m[i][l]*p[l][j];
           t[i][j] = s;
    }
} // multi_m 

/* 2. Транспонирование квадратной матрицы m
*  размером n*n. По окончании работы в матрице 
*  m хранится результат транспонирования.
*/
void tr_m(doublereal **m, integer n) {
    for (integer i = 1; i < n; ++i)
        for (integer j = 0; j < i; ++j) {
            doublereal buf = m[i][j];
            m[i][j] = m[j][i];
            m[j][i] = buf;
        }
} // tr_m

/* 3. Возвращает максимальный внедиагональный
* элемент для симметричной матрицы A размером 
* n*n. Позиция максимального элемента A[f][g].
* Это медленная реализация, т.к. она не использует
* информацию о предыдущих поисках максимального
* элемента в матрице А.
*/
doublereal max_el(doublereal **A, integer n, integer& f, integer& g) {
   doublereal max = A[0][1];
   f=0; g=1; // стартовое значение
   for (integer j = 1; j < n; ++j)
      for (integer i = 0; i < j; ++i) {
        if (A[i][j] > max) {
            max = A[i][j];
            f = i; g = j;
        }
    }
    return max;
 } // max_el

/* 4. Копирует вторую матрицу в первую: A=B.
* Матрицы квадратные размером n*n
*/
void matr_copy(doublereal **A1, doublereal **A2, integer n) {
   for (integer i = 0; i < n; ++i)
      for (integer j = 0; j < n; ++j)
		  A1[i][j]=A2[i][j];
}

/* 5. Быстрое умножение двух квадратных матриц специального 
* вида размера n*n (левое умножение):
*                A=hiс*A.
* Здесь hic -  несимметричная транспонированная матрица вращения:
* hic[f][f] = cosfi;
* hic[g][g] = cosfi;
* hic[f][g] = +sinfi;
* hic[g][f] = -sinfi;
* Здесь f и g позиции ненулевых элементов.
* Нумерация начинается с нуля.
* По окончании работы результат хранится в исходной матрице A.
* Теперь матрица hiс передаётся только как свои четыре особых элемента
* что позволяет существенно экономить память и быстродействие.
*/
void multi_m_left(doublereal **A, doublereal **rab, integer n, integer f, integer g, doublereal cosfi, doublereal sinfi) {
	/* Устаревший неэффективный но понятный код:
    for (integer i = 0; i < n; ++i)
       for (integer j = 0; j < n; ++j) {
		   if ((i!=f) && (i!=g)) {
			   t[i][j]=A[i][j];
		   }
		   else if (i==f) {
			   //t[i][j]=hic[f][f]*A[f][j]+hic[f][g]*A[g][j];
               t[i][j]=cosfi*A[f][j]+sinfi*A[g][j];
		   }
		   else if (i==g) {
			   //t[i][j]=hic[g][f]*A[f][j]+hic[g][g]*A[g][j];
			   t[i][j]=-sinfi*A[f][j]+cosfi*A[g][j];
		   }
    }
	*/
    
	// Теперь результат умножения возвращается прямо в матрице А
	// В качестве рабочего используется массив rab размерности
	// 2*n. Трудоёмкость операции всего 4*n умножений.
	for (integer j = 0; j < n; ++j) {
	   rab[0][j]=cosfi*A[f][j]+sinfi*A[g][j];
	   rab[1][j]=-sinfi*A[f][j]+cosfi*A[g][j];
	}
    for (integer j = 0; j < n; ++j) {
	   A[f][j]=rab[0][j];
	   A[g][j]=rab[1][j];
	}

} // multi_m_left 

/* 6. Быстрое умножение двух квадратных матриц специального 
* вида размера n*n (правое умножение):
*                A=A*hi.
* Здесь hi - несимметричная матрица вращения:
* hi[f][f] = cosfi;
* hi[g][g] = cosfi;
* hi[f][g] = -sinfi;
* hi[g][f] = +sinfi;
* Здесь f и g позиции ненулевых элементов.
* Нумерация начинается с нуля.
* По окончании работы результат хранится в исходной матрице A.
* Теперь матрица hi передаётся только как свои четыре особых элемента
* что позволяет существенно экономить память и быстродействие.
*/
void multi_m_right(doublereal **A, doublereal **rab, integer n, integer f, integer g, doublereal cosfi, doublereal sinfi) {
	/* Неэффективное умножение
    for (integer i = 0; i < n; ++i)
       for (integer j = 0; j < n; ++j) {
		   if ((j!=f) && (j!=g)) {
			   t[i][j]=A[i][j];
		   }
		   else if (j==f) {
			   //t[i][j]=A[i][f]*hi[f][f]+A[i][g]*hi[g][f];
               t[i][j]=A[i][f]*cosfi+A[i][g]*sinfi;
		   }
		   else if (j==g) {
			   //t[i][j]=A[i][f]*hi[f][g]+A[i][g]*hi[g][g];
			   t[i][j]=-A[i][f]*sinfi+A[i][g]*cosfi;
		   }
    }
	*/

	// Теперь результат умножения возвращается прямо в матрице А
	// В качестве рабочего используется массив rab размерности
	// 2*n. Трудоёмкость операции всего 4*n умножений.
	for (integer i = 0; i < n; ++i) {
	   rab[0][i]=A[i][f]*cosfi+A[i][g]*sinfi; // f
	   rab[1][i]=-A[i][f]*sinfi+A[i][g]*cosfi; // g
	}
    for (integer i = 0; i < n; ++i) {
		A[i][f]=rab[0][i];
		A[i][g]=rab[1][i];
	}

} // multi_m_right 


/* Оригинальный алгоритм Якоби от 1846 года. 
* Решает полную проблему собственных значений в форме
*              A-lambda_scal*E=0
*  методом вращений. См. Ревизников лекции книга.
*  Симметричная положительно определённая матрица A 
*  размером nodes*nodes. Матрица A в результате работы  
*  портится ( на диагонали у её испорченного варианта будут СЗ).
*  В квадратной матрице U по столбцам хранятся 
*  собственные векторы. В векторе lambda находится список
*  собственных значений.
*  Процесс нахождения векторов и СЗ является итерационным,
*  Его точность характеризуется значением epsilon.
*  На каждой итерации делается 12*nodes матричных умножений.
*  Дополнительная память равна 2*nodes.
*  EIGEM - метод Якоби.
*/
void jacobi_matrix_simple(doublereal **A, doublereal **U, doublereal *lambda, integer nodes, doublereal epsilon) {

	// значение этой постоянной нужно подобрать импирически.
    const doublereal eps=1e-10; // точность  с которой элементы проверяются на равенство,
	
	integer i,j; // счётчики цикла for
    integer im , jm; // позиция максимального элемента
    integer p = 1; // номер итерации
	doublereal maxij; // максимальный элемент
	doublereal fi; // значение угла
	doublereal cosfi, sinfi; // значение косинуса и синуса угла fi

	/* При быстром умножении матриц этот код не требуется
	// матрица  вращения
	doublereal **hi=new doublereal*[nodes];
    for (i = 0; i < nodes; ++i) hi[i]=new doublereal[nodes];

	// инициализация всей матрицы выполняется один раз:
    for (i = 0; i < nodes; ++i)
         for (j = 0; j < nodes; ++j) {
            if (i == j)
                hi[i][j] = 1.0;
            else hi[i][j] = 0.0;
    }

	// вспомогательная матрица вращения (копия)
	doublereal **hic=new doublereal*[nodes];
    for (i = 0; i < nodes; ++i) hic[i]=new doublereal[nodes];

	// инициализация вспомогательной матрицы один раз:
	matr_copy(hic,hi,nodes); // тоже единичная матрица 
	*/

	// вспомогательная матрица для умножения
    // Устаревший неэффективный по памяти вариант.
    //doublereal **b=new doublereal*[nodes];
    //for (i = 0; i < nodes; ++i) b[i]=new doublereal[nodes];

	// Рабочий массив.
	doublereal **rab=new doublereal*[2];
    for (i = 0; i < 2; ++i) rab[i]=new doublereal[nodes];

    maxij = max_el(A, nodes,im,jm);
    
	// каждую итерацию делается 12*nodes умножений.
    while (fabs(maxij) > epsilon) {
       
       
	   // Вычисление угла:
	   if (fabs(A[im][im]-A[jm][jm])<eps) {
		   // особый случай значения равны
		   fi=M_PI/4.0;
	   }
	   else fi= atan(2*maxij/(A[im][im]-A[jm][jm]))/2;
       
       // Нахождение тригонометрических функций
	   // от угла fi:
       cosfi = cos(fi);
	   sinfi = sin(fi);
 
	   /* При быстром умножении этот закомментированный код не используется.
	   // матрица вращения не является симметричной:
	   // инициализация матрицы вращения
       hi[im][im] = cosfi;
       hi[jm][jm] = cosfi;
       hi[im][jm] = -sinfi;
       hi[jm][im] = sinfi;
	   // транспонированный вариант: 
       hic[im][im] = cosfi;
       hic[jm][jm] = cosfi;
       hic[im][jm] = +sinfi; // транспонирование.
       hic[jm][im] = -sinfi;
	   */
 
       //  инициализация матрицы СВ которые хранятся по столбцам
	   if (p==1) {
		   //matr_copy(U,hi,nodes);
           for (i = 0; i < nodes; ++i)
               for (j = 0; j < nodes; ++j) {
                   if (i == j)
                      U[i][j] = 1.0;
                   else U[i][j] = 0.0;
               }
		   U[im][im] = cosfi;
           U[jm][jm] = cosfi;
           U[im][jm] = -sinfi;
           U[jm][im] = sinfi;

	   } else {
            //multi_m(U,hi,b, nodes);
            multi_m_right(U, rab, nodes, im, jm, cosfi, sinfi); // Быстрое умножение
			//matr_copy(U,b,nodes); // теперь вместо промежуточной матрицы b используется 
			// экономичная матрица rab. Эффективность операции 4xnodes умножений.
	   }
      
       //multi_m(hic, A, b, nodes); // b=transpose(H)*A
	   multi_m_left(A, rab, nodes, im, jm, cosfi, sinfi); // Быстрое умножение: 4xnodes операций умножения
       //multi_m(b, hi, A, nodes); // A=b*H.
	   multi_m_right(A, rab, nodes, im, jm, cosfi, sinfi); // Быстрое умножение: 4xnodes операций умножения
 
	   /* При быстром умножении этот закомментированный код не используется.
	   // восстановление матриц вращения:
       hi[im][im] = 1.0;
       hi[jm][jm] = 1.0;
       hi[im][jm] = 0.0;
       hi[jm][im] = 0.0;
	   // восстановление копии матрицы вращения:
       hic[im][im] = 1.0;
       hic[jm][jm] = 1.0;
       hic[im][jm] = 0.0;
       hic[jm][im] = 0.0;
	   */

	   maxij = max_el(A, nodes,im,jm); // определение максимума ресурсоёмкая операция.
       p++; // переход к следующей итерации

    } // while

    for (i = 0; i < nodes; ++i) lambda[i]=A[i][i]; //  СЗ

} // jacobi_matrix_simple

/* Следующие несколько функций применяются для GSEP
*  в целях упорядочивания по возрастанию набора 
*  собственных значений.
*/

// Пузырьковая сортировка.
void BubbleSortGSEP1(doublereal *a, integer *mask, integer n) {
   integer i=0, j=0, k=0;
   doublereal x;

   for (i=1; i<n; ++i) {
	   for (j=n-1; j>=i; j--) {
		   if (a[j-1] > a[j]) {
			   // swap
			   x=a[j-1];
			   a[j-1]=a[j];
			   a[j]=x;
               k=mask[j-1];
			   mask[j-1]=mask[j];
			   mask[j]=k;
		   }
	   }
   }
} // BubbleSortGSEP1


/* Первая обобщённая симметричная проблема собственных значений
*   GSEP1:  A*x-lambda_scal*B*x=0;
*   Разложение Холесского: B=L*transpose(L);
*   L - нижняя треугольная, transpose(L) - верхняя треугольная.
*/
void GSEP1(doublereal **A1, doublereal **A2, doublereal **U, doublereal *lambda, integer *mask, integer nodes, doublereal epsilon) {

	// Разложение Холесского: замена B верхним и нижним 
	// треугольными множителями.
#if doubleprecision == 1
		A2[0][0] = sqrt(A2[0][0]);
		A2[1][0] /= A2[0][0];
		A2[0][1] = A2[1][0];
		A2[1][1] = sqrt(A2[1][1] - A2[1][0] * A2[1][0]);
	
#else 
		A2[0][0] = sqrtf(A2[0][0]);
		A2[1][0] /= A2[0][0];
		A2[0][1] = A2[1][0];
		A2[1][1] = sqrtf(A2[1][1] - A2[1][0] * A2[1][0]);
#endif

	integer irow,irow1;
	integer icol, icol1;
	doublereal sum=0.0;
	integer k=0;
	for (irow=2; irow<nodes; irow++) {
		irow1=irow-1;
		A2[irow][0]/=A2[0][0];
        A2[0][irow]=A2[irow][0];
		for (icol=1; icol<=irow1; icol++) {
			icol1=icol-1;
            sum=0.0;
            for (k=0; k<=icol1; ++k) sum+=A2[irow][k]*A2[icol][k];
			A2[irow][icol]=(A2[irow][icol]-sum)/A2[icol][icol];
			A2[icol][irow]=A2[irow][icol];
		}
		sum=0.0;
		for (k=0; k<=irow1; ++k) sum+=A2[irow][k]*A2[irow][k];
#if doubleprecision == 1 
			A2[irow][irow] = sqrt(A2[irow][irow] - sum);
		
#else 
			A2[irow][irow] = sqrtf(A2[irow][irow] - sum);
#endif
	}

	printf("L*LT 10...\n");
   // TODO: дальше и до конца функции идёт чрезвычайно
	// неэффективный кусок кода:

    integer i=0, j=0;

	for (i=0; i<nodes; ++i) mask[i]=i;

	// нижняя треугольная матрица
	doublereal **L = nullptr;
	L=new doublereal*[nodes];
    for (i = 0; i < nodes; ++i) L[i]=new doublereal[nodes];

	/* Этот закомментированный кусок кода относится к медленной реализации:
	// Если использовать медленную реализацию то это надо раскоментировать.
	// инициализация всей матрицы выполняется один раз:
    for (i = 0; i < nodes; ++i)
         for (j = 0; j < nodes; ++j) {
            if (j > i)
                L[i][j] = 0.0;
            else L[i][j] = B[i][j];
    }
	*/

    /*
    // верхняя треугольная матрица
    doublereal **LT=new doublereal*[nodes];
    for (i = 0; i < nodes; ++i) LT[i]=new doublereal[nodes];

	// инициализация всей матрицы выполняется один раз:
    for (i = 0; i < nodes; ++i)
         for (j = 0; j < nodes; ++j) {
            if (j < i)
                LT[i][j] = 0.0;
            else LT[i][j] = B[i][j];
    }
	*/

    // вспомогательная матрица для умножения
	doublereal **b = nullptr;
	b=new doublereal*[nodes];
    for (i = 0; i < nodes; ++i) b[i]=new doublereal[nodes];

	// Ac копия матрицы А
	doublereal **Ac = nullptr;
	Ac=new doublereal*[nodes];
    for (i = 0; i < nodes; ++i) Ac[i]=new doublereal[nodes];
	matr_copy(Ac,A1,nodes); // сохранение А TODO временно потом удалить

	// Медленная реализация
    //inverse_matrix_simple(L,nodes); // нахождение L^(-1)
	//multi_m(L,A,b,nodes); // b=(L^(-1))*A;
	//matr_copy(A,b,nodes); // A=(L^(-1))*A;

	// Более быстрая реализация
    // A=(L^(-1))*A;
	for (i=0; i < nodes; ++i) {
		A1[0][i]/=A2[0][0];

	    for (irow=1; irow<nodes; irow++) {
		    irow1=irow-1;
		    sum=0.0;
		    for (icol=0; icol<=irow1; icol++) sum+=A2[irow][icol]*A1[icol][i];
            A1[irow][i]=(A1[irow][i]-sum)/A2[irow][irow];
	    }
	}

    printf("(L^(-1))*A 20...\n");
	//matr_copy(L,LT,nodes); // L=transpose(L); т.к. матрица L больше не нужна:
	// дальше везде под именем L используется transpose(L).
    
    // L=LT: L=transpose(L);
	// Теперь L верхняя треугольная матрица
    for (i = 0; i < nodes; ++i)
         for (j = 0; j < nodes; ++j) {
            if (j < i)
                L[i][j] = 0.0;
            else L[i][j] = A2[i][j];
    }

    // Эта матрица уже приведена к верхнетреугольному виду,
    // поэтому второй раз её к такому виду приводить не понадобиться поэтому true.
    inverse_matrix_simple(L,nodes,true); // нахождение (transpose(L))^(-1)
	 
	multi_m(A1,L,b,nodes); // b=(L^(-1))*A*(transpose(L))^(-1).
    matr_copy(A1,b,nodes); // A=(L^(-1))*A*(transpose(L))^(-1).

	printf("C 30...\n");

	jacobi_matrix_simple(A1,U,lambda,nodes,epsilon); // нахождение СВ и СЗ с заданной точностью.

    printf("C 90...\n");

	BubbleSortGSEP1(lambda,mask,nodes); // упорядочивание собственных значений.
	multi_m(L,U,b,nodes); // b=((transpose(L))^(-1))*U
    matr_copy(U,b,nodes); // собственные вектора.

	/* проверка найденных собственных значений.
    multi_m(Ac,U,b,nodes); // b=A1*U
    matr_copy(L,U,nodes); // L=U
	tr_m(L,nodes);
	multi_m(L,b,Ac,nodes); // Ac=transpose(U)*A1*U

	doublereal *test=new doublereal[nodes];
    for (integer i=0; i<nodes; ++i) test[i]=Ac[i][i];
    BubbleSortGSEP1(test,mask,nodes); 
    for (integer i=0; i<8; ++i) printf("%.2f ",test[i]/M_PI/M_PI); // собственные значения
	printf("\n");
	*/
	if (L != nullptr) {
		delete[] L;
	}
	if (b != nullptr) {
		delete[] b;
	}
	//delete LT;
} // GSEP1

/* метод Гаусса для ленточной матрицы A размером
*              nodes x 2*icolx+1, где
*   2*icolx+1 - ширина ленты. Под тем что матрица
*  A ленточная понимается то что ненулевые элементы
*  матрицы содержатся только внутри ленты.
*  b - вектор правой части СЛАУ, x - вектор решение.
*  Нумерация элементов начинается с нуля.
*  Для положительно определённых возможно несимметричных
*  матриц А, которые задаются своей лентой.
*  Гаусс Карл Фридрих 1777-1855.
*  В результате работы матрица А портится.
*/
void eqsolve_lenta_gauss(doublereal **A, integer nodes, integer icolx, doublereal *b, doublereal *x) {

	const doublereal eps=1e-300; // для сравнения с нулём
	doublereal dCik, dSum=0.0;
	integer max;

	integer *move=new integer[nodes]; // массив сдвигов.
	integer i=0, j=0, k=0; // счётчики цикла for
	for (i=0; i<nodes; ++i) move[i]=icolx-i; // инициализация массива сдвигов

	for (i=0; i<nodes; ++i) x[i]=0.0; // инициализация

	// прямой ход метода Гаусса
	// приведение к верхнему треугольному виду:

	// по всем столбцам слева направо
	for (k=0; k<nodes; ++k) {
        max=myi_min(k+icolx,nodes-1);
		// цикл по всем строкам ниже строки с номером k
		for (i=k+1; i<=max; ++i) {
			// применяется только в том случае
			// если элемент ненулевой
			// это должно несколько ускорить счёт.
			if ((i<nodes)&&(fabs(A[i][k+move[i]]) > eps)) {
               
                if(fabs(A[k][k+move[k]])<eps){
			          // решение не может быть получено, т.к.
			          // на диагонали находится ноль.
	                  printf("\nSolution is not exist! divizion by zero...\n");
	                //  system("pause");
					  system("pause");
		              exit(0);
	            }

                // обработка фиксированной строки с номером i
				dCik=A[i][k+move[i]]/A[k][k+move[k]];
				// преаобразование матрицы к верхнетреугольному виду:
				for (j=k; j<=max; ++j) A[i][j+move[i]] -= dCik*A[k][j+move[k]];
				b[i]-= dCik*b[k]; // преобразование правой части
			}
		}
	}

    // Теперь когда матрица приведена к верхнетреугольному виду
	// можно совершить обратный ход метода Гаусса:
	for (k=nodes-1; k>=0; k--) {
        dSum=0.0; // обнуление сумматора
		max= myi_min(k+icolx,nodes-1);
		for (i=k+1; i<=max; ++i) dSum+= A[k][i+move[k]]*x[i];
		x[k]=(b[k]-dSum)/A[k][k+move[k]];
	}

	delete[] move;
	move = nullptr;

}  // eqsolve_lenta_gauss

// Метод (Якоби) Гаусса-Зейделя
// для решения СЛАУ с матрицей А n*n
// возможно несимметричной, но с диагональным 
// преобладанием. Матрица А предполагается
// полностью заполненной (неразреженной).
// b - правая часть, x - уточняемое решение, 
// eps - точность определения решения.
// omega - импирически подобранный параметр релаксации.
void Seidel(doublereal **A, doublereal *b, doublereal *x, integer n, doublereal eps, doublereal omega) {
	integer i,j;
	doublereal s1, s2, s, v, m;
	bool bdiag=true;

	// Исследуем сходимость
	for (i=0; i<n; ++i) {
		s=0.0;
		for (j=0; j<n; ++j) {
			if (j!=i) s+=fabs(A[i][j]);
		}
		if (s>=fabs(A[i][i])) {
			bdiag=false;
		}
	}
	if (!bdiag) {
		printf("net diagonalnogo preobladaniq...");
		//system("pause");
		system("pause");
	}

	do {
		m=0.0;
		for (i=0; i<n; ++i) {
			// Вычисляем суммы
			s1=s2=0.0; 
			for (j=0; j<=i-1; ++j) s1+=A[i][j]*x[j];
			for (j=i+1; j<n; ++j) s2+=A[i][j]*x[j];
			// Вычисляем новое приближение и погрешность
			v=x[i];
			x[i]=omega*(b[i]-s1-s2)/A[i][i]+(1-omega)*x[i];

			if (fabs(v-x[i])>m) m=fabs(v-x[i]);
		}

	} while (m > eps);

} // Seidel

// возвращает максимальное из двух 
// вещественных чисел.
/*
doublereal fmax(doublereal fA, doublereal fB) {
	doublereal r=fB;
	if (fA > fB) r=fA;
	return r;
} // fmax 
*/

// применяется для уравнения поправки давления 
// в случае когда на всей границе стоят условия Неймана.
void SOR(equation* &sl, doublereal* &x, integer n) {
	doublereal rURF=1.855; // параметр верхней релаксации
	// пороговое значение невязки
	doublereal eps = 1e-3;
	doublereal ptilda;
	doublereal sE,sW,sN, sS;
	integer i=0,j=0, kend=3000; // счётчик цикла for
	doublereal dmax=1.0;
	while ((dmax>eps) && (j<kend)) {
		dmax=0.0;
	    for (i=0; i<n; ++i) {
            if (sl[i].iE>-1) sE=sl[i].ae*x[sl[i].iE]; else sE=0.0;
		    if (sl[i].iW>-1) sW=sl[i].aw*x[sl[i].iW]; else sW=0.0;
		    if (sl[i].iN>-1) sN=sl[i].an*x[sl[i].iN]; else sN=0.0;
		    if (sl[i].iS>-1) sS=sl[i].as*x[sl[i].iS]; else sS=0.0;
		    ptilda=(sE+sW+sN+sS+sl[i].b)/sl[i].ap;
		    //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			dmax+=fabs(sl[i].ap*(ptilda-x[sl[i].iP]));
		    x[sl[i].iP]=x[sl[i].iP]+rURF*(ptilda-x[sl[i].iP]);
	    }
		dmax/=n;
		printf("%e \n",dmax);
		j++;
	}

} // SOR

// Нормализация невязок, так чтобы максимально приблизится
// к выводу программы ANSYS Icepak 17.2.
typedef struct  TResidualNormalization {
	// 5.05.2017
	 integer iM;
	doublereal resVX0;
	doublereal resVY0;
	doublereal resVZ0;
	integer icVX, icVY, icVZ;
	// Спаларт Аллмарес [1992].
	doublereal resNUSHA0;
	integer icNUSHA;
	// SST k-omega Menter [1993].
	doublereal reskMenter0;
	integer ickMenter;
	doublereal resomegaMenter0;
	integer icomegaMenter;
	// Двухслойная модель на основе стандартной k-epsilon модели [2001].
	doublereal reskStandart_k_epsilon0;
	integer ickStandart_k_epsilon;
	doublereal resepsilonStandart_k_epsilon0;
	integer icepsilonStandart_k_epsilon;
	// Ментор Лантгрии
	doublereal resGamma0;
	integer icGamma0;
	doublereal resReTheta0;
	integer icReTheta0;


	TResidualNormalization() {
		// 5.05.2017
		iM = 3;
		resVX0 = 1.0;
		resVY0 = 1.0;
		resVZ0 = 1.0;
		icVX = 0; icVY = 0; icVZ = 0;
		// Спаларт Аллмарес [1992].
		resNUSHA0 = 1.0;
		icNUSHA = 0;
		// SST k-omega Menter [1993].
		 reskMenter0 = 1.0;
		 ickMenter = 0;
		 resomegaMenter0 = 1.0;
		icomegaMenter = 0;
		// Двухслойная модель на основе стандартной k-epsilon модели [2001].
		reskStandart_k_epsilon0 = 1.0;
		ickStandart_k_epsilon = 0;
		resepsilonStandart_k_epsilon0 = 1.0;
		icepsilonStandart_k_epsilon = 0;
		// Ментор Лантгрии [2009]
		resGamma0 = 1.0;
		icGamma0 = 0;
		resReTheta0 = 1.0;
		icReTheta0 = 0;

	}

} ResidualNormalization;

ResidualNormalization fluent_resformat;

// здесь будет вычислена невязка по формуле которая используется
// в комерческой программе ANSYS fluent. Имея невязки вычисляемые по 
// формуле Fluent можно будет произвести с ним сравнение и настроить быстродействие
// программы AliceFlow_v0_07.
// Данная невязка используется для всех полевых величин кроме поправки давления.
// Информация о формуле по которой вычисляются невязки взята из icepak user Guide. Threory chapter.
// 26.09.2016 теперь на АЛИС сетке.
doublereal fluent_residual_for_x(equation3D* &sl, equation3D_bon* &slb, doublereal* &x, integer maxelm, integer maxbound, integer iVar) {
	doublereal r=0.0;

	doublereal fsum1=0.0, fsum2=0.0;

	// внутренние контрольные объёмы.
	doublereal fsum1_loc1 = 0.0;
	doublereal fsum2_loc1 = 0.0;

#pragma omp parallel for reduction (+ : fsum1_loc1, fsum2_loc1)
	for (integer i=0; i<maxelm; ++i) {
		// числитель
		doublereal sE=0.0,sW=0.0,sN=0.0,sS=0.0,sT=0.0,sB=0.0;
		if (sl[i].iE>-1) sE=sl[i].ae*x[sl[i].iE]; else sE=0.0;
		if (sl[i].iW>-1) sW=sl[i].aw*x[sl[i].iW]; else sW=0.0;
		if (sl[i].iN>-1) sN=sl[i].an*x[sl[i].iN]; else sN=0.0;
		if (sl[i].iS>-1) sS=sl[i].as*x[sl[i].iS]; else sS=0.0;
        if (sl[i].iT>-1) sT=sl[i].at*x[sl[i].iT]; else sT=0.0;
		if (sl[i].iB>-1) sB=sl[i].ab*x[sl[i].iB]; else sB=0.0;

		doublereal fbuf1 = sE + sW + sN + sS + sT + sB;
		fbuf1 += sl[i].b;
		//fbuf1 += rthdsd[sl[i].iP];

		if (b_on_adaptive_local_refinement_mesh) {
			doublereal sE2 = 0.0, sW2 = 0.0, sN2 = 0.0, sS2 = 0.0, sT2 = 0.0, sB2 = 0.0;
			if (sl[i].bE2) {
				if (sl[i].iE2 > -1) sE2 = sl[i].ae2 * x[sl[i].iE2]; else sE2 = 0.0;
			}
			if (sl[i].bW2) {
				if (sl[i].iW2 > -1) sW2 = sl[i].aw2 * x[sl[i].iW2]; else sW2 = 0.0;
			}
			if (sl[i].bN2) {
				if (sl[i].iN2 > -1) sN2 = sl[i].an2 * x[sl[i].iN2]; else sN2 = 0.0;
			}
			if (sl[i].bS2) {
				if (sl[i].iS2 > -1) sS2 = sl[i].as2 * x[sl[i].iS2]; else sS2 = 0.0;
			}
			if (sl[i].bT2) {
				if (sl[i].iT2 > -1) sT2 = sl[i].at2 * x[sl[i].iT2]; else sT2 = 0.0;
			}
			if (sl[i].bB2) {
				if (sl[i].iB2 > -1) sB2 = sl[i].ab2 * x[sl[i].iB2]; else sB2 = 0.0;
			}
			doublereal sE3 = 0.0, sW3 = 0.0, sN3 = 0.0, sS3 = 0.0, sT3 = 0.0, sB3 = 0.0;
			if (sl[i].bE3) {
				if (sl[i].iE3 > -1) sE3 = sl[i].ae3 * x[sl[i].iE3]; else sE3 = 0.0;
			}
			if (sl[i].bW3) {
				if (sl[i].iW3 > -1) sW3 = sl[i].aw3 * x[sl[i].iW3]; else sW3 = 0.0;
			}
			if (sl[i].bN3) {
				if (sl[i].iN3 > -1) sN3 = sl[i].an3 * x[sl[i].iN3]; else sN3 = 0.0;
			}
			if (sl[i].bS3) {
				if (sl[i].iS3 > -1) sS3 = sl[i].as3 * x[sl[i].iS3]; else sS3 = 0.0;
			}
			if (sl[i].bT3) {
				if (sl[i].iT3 > -1) sT3 = sl[i].at3 * x[sl[i].iT3]; else sT3 = 0.0;
			}
			if (sl[i].bB3) {
				if (sl[i].iB3 > -1) sB3 = sl[i].ab3 * x[sl[i].iB3]; else sB3 = 0.0;
			}
			doublereal sE4 = 0.0, sW4 = 0.0, sN4 = 0.0, sS4 = 0.0, sT4 = 0.0, sB4 = 0.0;
			if (sl[i].bE4) {
				if (sl[i].iE4 > -1) sE4 = sl[i].ae4 * x[sl[i].iE4]; else sE4 = 0.0;
			}
			if (sl[i].bW4) {
				if (sl[i].iW4 > -1) sW4 = sl[i].aw4 * x[sl[i].iW4]; else sW4 = 0.0;
			}
			if (sl[i].bN4) {
				if (sl[i].iN4 > -1) sN4 = sl[i].an4 * x[sl[i].iN4]; else sN4 = 0.0;
			}
			if (sl[i].bS4) {
				if (sl[i].iS4 > -1) sS4 = sl[i].as4 * x[sl[i].iS4]; else sS4 = 0.0;
			}
			if (sl[i].bT4) {
				if (sl[i].iT4 > -1) sT4 = sl[i].at4 * x[sl[i].iT4]; else sT4 = 0.0;
			}
			if (sl[i].bB4) {
				if (sl[i].iB4 > -1) sB4 = sl[i].ab4 * x[sl[i].iB4]; else sB4 = 0.0;
			}


			
			// Учёт расширенного шаблона на АЛИС сетке.
			fbuf1 += sE2 + sW2 + sN2 + sS2 + sT2 + sB2;
			fbuf1 += sE3 + sW3 + sN3 + sS3 + sT3 + sB3;
			fbuf1 += sE4 + sW4 + sN4 + sS4 + sT4 + sB4;
		}
		fsum1_loc1 +=fabs(fbuf1-sl[i].ap*x[sl[i].iP]);
		fsum2_loc1 +=fabs(sl[i].ap*x[sl[i].iP]); // знаменатель.
	}

	
	fsum1 += fsum1_loc1;
	fsum2 += fsum2_loc1;

	// граничные контрольные объёмы.
	doublereal fsum1_loc = 0.0;
	doublereal fsum2_loc = 0.0;

#pragma omp parallel for reduction (+ : fsum1_loc, fsum2_loc)
	for (integer i=0; i<maxbound; ++i) {
		// числитель
		doublereal sI;
		if (slb[i].iI>-1) sI=slb[i].ai*x[slb[i].iI]; else sI=0.0;
		fsum1_loc+=fabs(sI+slb[i].b-slb[i].aw*x[slb[i].iW]);
		//fsum1 += fabs(sI + rthdsd[slb[i].iW] - slb[i].aw*x[slb[i].iW]);
		// знаменатель
		fsum2_loc+=fabs(slb[i].aw*x[slb[i].iW]);
	}

	fsum1 += fsum1_loc;
	fsum2 += fsum2_loc;
	
	switch (iVar) {
	  case VELOCITY_X_COMPONENT: fluent_resformat.icVX++;
		  if (fluent_resformat.icVX == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.resVX0 = fsum1 / fsum2;
			  }
		  }
		       break;
	  case VELOCITY_Y_COMPONENT: fluent_resformat.icVY++;
		  if (fluent_resformat.icVY == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.resVY0 = fsum1 / fsum2;
			  }
		  }
		  break;
	  case VELOCITY_Z_COMPONENT: fluent_resformat.icVZ++; 
		  if (fluent_resformat.icVZ == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.resVZ0 = fsum1 / fsum2;
			  }
		  }
		  break;
	  case NUSHA: fluent_resformat.icNUSHA++;
		  if (fluent_resformat.icNUSHA == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.resNUSHA0 = fsum1 / fsum2;
			  }
		  }
		  break;
	  case GAMMA_LANGTRY_MENTER: fluent_resformat.icGamma0++;
		  if (fluent_resformat.icGamma0 == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.resGamma0 = fsum1 / fsum2;
			  }
		  }
		  break;
	  case RE_THETA_LANGTRY_MENTER: fluent_resformat.icReTheta0++;
		  if (fluent_resformat.icReTheta0 == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.resReTheta0 = fsum1 / fsum2;
			  }
		  }
		  break;
	  case TURBULENT_KINETIK_ENERGY: fluent_resformat.ickMenter++;
		  if (fluent_resformat.ickMenter == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.reskMenter0 = fsum1 / fsum2;
			  }
		  }
		  break;
	  case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: fluent_resformat.icomegaMenter++;
		  if (fluent_resformat.icomegaMenter == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.resomegaMenter0 = fsum1 / fsum2;
			  }
		  }
		  break;
	  case TURBULENT_KINETIK_ENERGY_STD_K_EPS: fluent_resformat.ickStandart_k_epsilon++;
		  if (fluent_resformat.ickStandart_k_epsilon  == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.reskStandart_k_epsilon0 = fsum1 / fsum2;
			  }
		  }
		  break;
	  case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: fluent_resformat.icepsilonStandart_k_epsilon++;
		  if (fluent_resformat.icepsilonStandart_k_epsilon == fluent_resformat.iM) {
			  if (fsum2 > 1.0e-41) {
				  fluent_resformat.resepsilonStandart_k_epsilon0 = fsum1 / fsum2;
			  }
		  }
		  break;
	}

	if (fsum2<1.0e-41) {
		r=0.0;
	} 
	else {
		r=fsum1/fsum2;
		// Нормализация.
		switch (iVar) {
		  case VELOCITY_X_COMPONENT: r = r / fluent_resformat.resVX0; break;
		  case VELOCITY_Y_COMPONENT: r = r / fluent_resformat.resVY0; break;
		  case VELOCITY_Z_COMPONENT: r = r / fluent_resformat.resVZ0; break;
		  case NUSHA: r = r / fluent_resformat.resNUSHA0; break;
		  case GAMMA_LANGTRY_MENTER: r = r / fluent_resformat.resGamma0; break;
		  case RE_THETA_LANGTRY_MENTER: r = r / fluent_resformat.resReTheta0; break;
		  case TURBULENT_KINETIK_ENERGY: r = r / fluent_resformat.reskMenter0; break;
		  case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: r = r / fluent_resformat.resomegaMenter0; break;
		  case TURBULENT_KINETIK_ENERGY_STD_K_EPS: r = r / fluent_resformat.reskStandart_k_epsilon0; break;
		  case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: r = r / fluent_resformat.resepsilonStandart_k_epsilon0; break;
		}
		
	}
	return r;
} // fluent_residual_for_x

  // здесь будет вычислена невязка по формуле которая используется
  // в комерческой программе ANSYS fluent. Имея невязки вычисляемые по 
  // формуле Fluent можно будет произвести с ним сравнение и настроить быстродействие
  // программы AliceFlow_v0_07.
  // Данная невязка используется для всех полевых величин кроме поправки давления.
  // Информация о формуле по которой вычисляются невязки взята из icepak user Guide. Threory chapter.
  // 26.09.2016 теперь на АЛИС сетке.
doublereal fluent_residual_for_x_new(equation3D* &sl, equation3D_bon* &slb, doublereal* &x, integer maxelm, integer maxbound, doublereal* &rthdsd, doublereal alpha_relax) {
	doublereal r = 0.0;

	doublereal fsum1 = 0.0, fsum2 = 0.0;

	// внутренние контрольные объёмы.
	for (integer i = 0; i<maxelm; ++i) {
		// числитель
		doublereal sE = 0.0, sW = 0.0, sN = 0.0, sS = 0.0, sT = 0.0, sB = 0.0;
		if (sl[i].iE>-1) sE = sl[i].ae*x[sl[i].iE]; else sE = 0.0;
		if (sl[i].iW>-1) sW = sl[i].aw*x[sl[i].iW]; else sW = 0.0;
		if (sl[i].iN>-1) sN = sl[i].an*x[sl[i].iN]; else sN = 0.0;
		if (sl[i].iS>-1) sS = sl[i].as*x[sl[i].iS]; else sS = 0.0;
		if (sl[i].iT>-1) sT = sl[i].at*x[sl[i].iT]; else sT = 0.0;
		if (sl[i].iB>-1) sB = sl[i].ab*x[sl[i].iB]; else sB = 0.0;

		doublereal fbuf1 = sE + sW + sN + sS + sT + sB;
		//fbuf1 += sl[i].b;
		fbuf1 += rthdsd[sl[i].iP];

		if (b_on_adaptive_local_refinement_mesh) {

			doublereal sE2 = 0.0, sW2 = 0.0, sN2 = 0.0, sS2 = 0.0, sT2 = 0.0, sB2 = 0.0;
			if (sl[i].bE2) {
				if (sl[i].iE2 > -1) sE2 = sl[i].ae2 * x[sl[i].iE2]; else sE2 = 0.0;
			}
			if (sl[i].bW2) {
				if (sl[i].iW2 > -1) sW2 = sl[i].aw2 * x[sl[i].iW2]; else sW2 = 0.0;
			}
			if (sl[i].bN2) {
				if (sl[i].iN2 > -1) sN2 = sl[i].an2 * x[sl[i].iN2]; else sN2 = 0.0;
			}
			if (sl[i].bS2) {
				if (sl[i].iS2 > -1) sS2 = sl[i].as2 * x[sl[i].iS2]; else sS2 = 0.0;
			}
			if (sl[i].bT2) {
				if (sl[i].iT2 > -1) sT2 = sl[i].at2 * x[sl[i].iT2]; else sT2 = 0.0;
			}
			if (sl[i].bB2) {
				if (sl[i].iB2 > -1) sB2 = sl[i].ab2 * x[sl[i].iB2]; else sB2 = 0.0;
			}
			doublereal sE3 = 0.0, sW3 = 0.0, sN3 = 0.0, sS3 = 0.0, sT3 = 0.0, sB3 = 0.0;
			if (sl[i].bE3) {
				if (sl[i].iE3 > -1) sE3 = sl[i].ae3 * x[sl[i].iE3]; else sE3 = 0.0;
			}
			if (sl[i].bW3) {
				if (sl[i].iW3 > -1) sW3 = sl[i].aw3 * x[sl[i].iW3]; else sW3 = 0.0;
			}
			if (sl[i].bN3) {
				if (sl[i].iN3 > -1) sN3 = sl[i].an3 * x[sl[i].iN3]; else sN3 = 0.0;
			}
			if (sl[i].bS3) {
				if (sl[i].iS3 > -1) sS3 = sl[i].as3 * x[sl[i].iS3]; else sS3 = 0.0;
			}
			if (sl[i].bT3) {
				if (sl[i].iT3 > -1) sT3 = sl[i].at3 * x[sl[i].iT3]; else sT3 = 0.0;
			}
			if (sl[i].bB3) {
				if (sl[i].iB3 > -1) sB3 = sl[i].ab3 * x[sl[i].iB3]; else sB3 = 0.0;
			}
			doublereal sE4 = 0.0, sW4 = 0.0, sN4 = 0.0, sS4 = 0.0, sT4 = 0.0, sB4 = 0.0;
			if (sl[i].bE4) {
				if (sl[i].iE4 > -1) sE4 = sl[i].ae4 * x[sl[i].iE4]; else sE4 = 0.0;
			}
			if (sl[i].bW4) {
				if (sl[i].iW4 > -1) sW4 = sl[i].aw4 * x[sl[i].iW4]; else sW4 = 0.0;
			}
			if (sl[i].bN4) {
				if (sl[i].iN4 > -1) sN4 = sl[i].an4 * x[sl[i].iN4]; else sN4 = 0.0;
			}
			if (sl[i].bS4) {
				if (sl[i].iS4 > -1) sS4 = sl[i].as4 * x[sl[i].iS4]; else sS4 = 0.0;
			}
			if (sl[i].bT4) {
				if (sl[i].iT4 > -1) sT4 = sl[i].at4 * x[sl[i].iT4]; else sT4 = 0.0;
			}
			if (sl[i].bB4) {
				if (sl[i].iB4 > -1) sB4 = sl[i].ab4 * x[sl[i].iB4]; else sB4 = 0.0;
			}

			// Учёт расширенного шаблона на АЛИС сетке.
			fbuf1 += sE2 + sW2 + sN2 + sS2 + sT2 + sB2;
			fbuf1 += sE3 + sW3 + sN3 + sS3 + sT3 + sB3;
			fbuf1 += sE4 + sW4 + sN4 + sS4 + sT4 + sB4;
		}
		
		
		fsum1 += fabs(fbuf1 - sl[i].ap*x[sl[i].iP] / alpha_relax);
		fsum2 += fabs(sl[i].ap*x[sl[i].iP] / alpha_relax); // знаменатель.
	}


	// граничные контрольные объёмы.
	for (integer i = 0; i<maxbound; ++i) {
		// числитель
		doublereal sI;
		if (slb[i].iI>-1) sI = slb[i].ai*x[slb[i].iI]; else sI = 0.0;
		//fsum1+=fabs(sI+slb[i].b-slb[i].aw*x[slb[i].iW]);
		fsum1 += fabs(sI + rthdsd[slb[i].iW] - slb[i].aw*x[slb[i].iW]);
		// знаменатель
		fsum2 += fabs(slb[i].aw*x[slb[i].iW]);
	}

	if (fsum2<1.0e-41) {
		r = 0.0;
	}
	else {
		r = fsum1 / fsum2;
	}
	return r;
} // fluent_residual_for_x_new

// Возвращает норму источникового члена (несбалансированные источники массы)
// в уравнении для поправки давления в форме соответствующей программе ANSYS fluent.
// По этому критерию можно судить о сходимости всей системы уравнений Навье-Стокса.
// Информация о невязках получена из руководства к программе icepak.
doublereal no_balance_mass_flux_fluent(doublereal* &b, doublereal operating_value_b, integer n) {
	// b - вектор несбалансированных источников массы.
	// n - размерность этого вектора.

	doublereal r=0.0;

	for (integer i=0; i<n; ++i) {
		r+=fabs(b[i]);
	}

	if (fabs(operating_value_b)<1.0e-30) {
        r=r/1.0; // чтобы избежать деления на вещественный ноль.
	}
	else {
	   r=r/operating_value_b;
	}
	return r;
} // no_balance_mass_flux_fluent

// применяется для уравнения поправки давления 
// в случае когда на всей границе стоят условия Неймана.
// Скорость сходимости очень и очень медленная,
// поэтому этот метод используется НИКОГДА.
// Основываясь на идеях алгоритма Федоренко этот метод можно применять как сглаживатель.
void SOR3D(equation3D* &sl, equation3D_bon* &slb, doublereal* &x, doublereal* x_cor, integer maxelm, integer maxbound, integer iVar, doublereal alpha) {
	// К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
	// скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
	//printf("SOR3D incomming...\n"); // debug.
	//system("pause");

	// Параметр релаксации лежит в интервале от 0.0 до 2.0.
	// Верхняя релаксация способна существенно ускорить вычислительный процесс. 
	// Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
	// rURF=1.5 можно использовать в уравнении на поправку давления.
	// В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
	// к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
	// в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
	// можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
	// коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость)
	//, поэтому наверное лучше не вводить коээффициент 
	// верхней релаксации равный 1.5 в уравнениях на скорость. 

    doublereal rURF=1.0; // параметр верхней релаксации
	switch (iVar) {
		case PAM: rURF=1.0;//1.855; //1.855; 
			      break;
		case VELOCITY_X_COMPONENT: rURF=1.0; //1.0;  // уравнения по скорости нелинейны и им нужна нижняя релаксация. Иначе возможна расходимость.
			      break; // эта нижняя релаксация определяется при формировании матрицы СЛАУ.
		case VELOCITY_Y_COMPONENT: rURF=1.0; //1.0;
			      break;
		case VELOCITY_Z_COMPONENT: rURF=1.0; //1.0;
			      break;
		default: rURF=1.0; break; // в остальных случаях.
	}
	// пороговое значение невязки
	doublereal eps = 1e-40;
	doublereal ptilda;
	doublereal sE,sW,sN,sS,sT,sB,sI;
	integer i=0,j=0, kend=10;//100; // Для целей пост сглаживания должно хватить 40 итераций.
	/*
	if (iVar==PAM) {
		kend=10000;
	}*/
	doublereal dmax=1.0;
	while ((dmax>eps) && (j<kend)) {  
		dmax=0.0;
        //#pragma omp parallel for private(i,ptilda,sE,sW,sN,sS,sT,sB) shared(maxelm,x,rURF,sl) schedule (guided)
	    for (i=0; i<maxelm; ++i) {
            if (sl[i].iE>-1) sE=sl[i].ae*x[sl[i].iE]; else sE=0.0;
		    if (sl[i].iW>-1) sW=sl[i].aw*x[sl[i].iW]; else sW=0.0;
		    if (sl[i].iN>-1) sN=sl[i].an*x[sl[i].iN]; else sN=0.0;
		    if (sl[i].iS>-1) sS=sl[i].as*x[sl[i].iS]; else sS=0.0;
            if (sl[i].iT>-1) sT=sl[i].at*x[sl[i].iT]; else sT=0.0;
		    if (sl[i].iB>-1) sB=sl[i].ab*x[sl[i].iB]; else sB=0.0;

			switch (iVar) {
			case PAM:  if (fabs(sl[i].ap)<1.0e-30) {
#if doubleintprecision == 1
				printf("division by zero in i=%lld internal node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
				printf("division by zero in i=%d internal node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						 
						  printf("Matrix construct PAM error...\n");
						  printf("Please, press any key to halt calculation...\n");
						 // system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=(sE+sW+sN+sS+sT+sB+sl[i].b)/sl[i].ap;
		                  //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			              //dmax=fmax(dmax,fabs(sl[i].ap*(ptilda-x[sl[i].iP]))); // Чебышева
				          dmax+=fabs((ptilda-x[sl[i].iP])); // сумма модулей
		                  x[sl[i].iP]=x[sl[i].iP]+rURF*(ptilda-x[sl[i].iP]);
					  }
					   break;
			case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: // нижняя релаксация на компоненты скорости.
				      if (fabs(sl[i].ap)<1.0e-30) {
#if doubleintprecision == 1
						  printf("division by zero in i=%lld internal node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						  printf("division by zero in i=%d internal node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						  
						  printf("Matrix construct Velocity error...\n");
						  printf("Please, press any key to halt calculation...\n");
						 // system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=alpha*(sE+sW+sN+sS+sT+sB+sl[i].b+(1.0-alpha)*sl[i].ap*x_cor[sl[i].iP]/alpha)/sl[i].ap;
		                  //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			              //dmax=fmax(dmax,fabs(sl[i].ap*(ptilda-x[sl[i].iP])/alpha)); // Чебышева
					      dmax+=fabs((ptilda-x[sl[i].iP])); // сумма модулей
		                  x[sl[i].iP]=x[sl[i].iP]+rURF*(ptilda-x[sl[i].iP]);
					   }
				       break;
			default: // по умолчанию без нижней реласкации.
				       if (fabs(sl[i].ap)<1.0e-30) {
#if doubleintprecision == 1
						   printf("division by zero in i=%lld boundary node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						   printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						 
						  printf("Matrix construct error...\n");
						  printf("Please, press any key to halt calculation...\n");
						 // system("pause");
						  system("pause");
						  exit(0);
					   }
					   else {
				          ptilda=(sE+sW+sN+sS+sT+sB+sl[i].b)/sl[i].ap;
		                  //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			              //dmax=fmax(dmax,fabs(sl[i].ap*(ptilda-x[sl[i].iP]))); // Чебышева
					      dmax+=fabs((ptilda-x[sl[i].iP])); // сумма модулей
		                  x[sl[i].iP]=x[sl[i].iP]+rURF*(ptilda-x[sl[i].iP]);
					   }
					   break;
				break;
			}
			
			/*
			if (0&&j==0) {// debug
				printf("ae=%e, aw=%e, an=%e, as=%e, at=%e, ab=%e, ap=%e, b=%e\n",sl[i].ae,sl[i].aw,sl[i].an,sl[i].as,sl[i].at,sl[i].ab,sl[i].ap,sl[i].b);
				system("pause");
			}
			*/
	    }
		//#pragma omp parallel for private(i,ptilda,sI) shared(maxbound,x,rURF,slb) schedule (guided)
		for (i=0; i<maxbound; ++i) {
			if (slb[i].iI>-1) sI=slb[i].ai*x[slb[i].iI]; else sI=0.0;

			switch (iVar) {
			case PAM:
				      if (fabs(slb[i].aw)<1.0e-30) {
#if doubleintprecision == 1
						  printf("division by zero in i=%lld boundary node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						  printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						  
						  printf("Matrix construct PAM error...\n");
						  printf("Please, press any key to halt calculation...\n");
						 // system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=(sI+slb[i].b)/slb[i].aw;
			              //dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[slb[i].iW]))); // Чебышева
					      dmax+=fabs((ptilda-x[slb[i].iW])); // сумма модулей.
			              //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
                          if (slb[i].iI==-1) x[slb[i].iW]=(ptilda);
			              else x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
					  }
					  break;
			case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT:
				      if (fabs(slb[i].aw)<1.0e-30) {
#if doubleintprecision == 1
						  printf("division by zero in i=%lld boundary node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						  printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						  
						  printf("Matrix construct Velocity error...\n");
						  printf("Please, press any key to halt calculation...\n");
						//  system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=(sI+slb[i].b)/slb[i].aw;
			              //dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[slb[i].iW]))); // Чебышева
					      dmax+=fabs((ptilda-x[slb[i].iW])); // сумма модулей
			              //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
                          if (slb[i].iI==-1) x[slb[i].iW]=(ptilda);
			              else x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
				          /*
				          ptilda=alpha*(sI+slb[i].b+(1.0-alpha)*slb[i].aw*x_cor[slb[i].iW]/alpha)/slb[i].aw;
			              dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[slb[i].iW])/alpha));
			              //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
                          if (slb[i].iI==-1) x[slb[i].iW]=(ptilda);
			              else x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]); 
					      */
					  }
				      break;
			default: // по умолчанию без нижней релаксации.
				      if (fabs(slb[i].aw)<1.0e-30) {
#if doubleintprecision == 1
						  printf("division by zero in i=%lld boundary node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						  printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						 
						  printf("Matrix construct error...\n");
						  printf("Please, press any key to halt calculation...\n");
                           // system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=(sI+slb[i].b)/slb[i].aw;
			              //dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[slb[i].iW]))); // Чебышева
					      dmax+=fabs((ptilda-x[slb[i].iW])); // сумма модулей.
			              //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
                          if (slb[i].iI==-1) x[slb[i].iW]=(ptilda);
			              else x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
					  }
				break;
			}
		}
		
		/*
		if (iVar==PAM)  {
		  //dmax/=maxelm;
		  if (j%1==0) {
		  #if doubleintprecision == 1
				 printf("%d %e \n", j+1, dmax);
		  #else
				 printf("%d %e \n", j+1, dmax);
		  #endif
			  
		  }
		}
		*/
		
		j++;
	}
	
	/*
	if (iVar==PAM)  {
	   printf("calc complete...\n");
       system("pause");
	}
	*/
	//printf("4000 %e \n", dmax);
	//system("pause");

} // SOR3D

// применяется для уравнения поправки давления 
// в случае когда на всей границе стоят условия Неймана.
// Скорость сходимости очень и очень медленная,
// поэтому этот метод используется НИКОГДА.
// Основываясь на идеях алгоритма Федоренко этот метод можно применять как сглаживатель.
void SOR3Dnow(equation3D* &sl, equation3D_bon* &slb, doublereal*& b,
	doublereal* &x, integer maxelm, integer maxbound, integer iVar,
	bool bprintMessage) {
	// К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
	// скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
	//printf("SOR3D incomming...\n"); // debug.
	//system("pause");

	// Параметр релаксации лежит в интервале от 0.0 до 2.0.
	// Верхняя релаксация способна существенно ускорить вычислительный процесс. 
	// Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
	// rURF=1.5 можно использовать в уравнении на поправку давления.
	// В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
	// к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
	// в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
	// можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
	// коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
	// поэтому наверное лучше не вводить коээффициент 
	// верхней релаксации равный 1.5 в уравнениях на скорость. 

    doublereal rURF=1.0; // параметр верхней релаксации
	switch (iVar) {
		case PAM: rURF=1.0;//1.855; //1.855; 
			      break;
		case VELOCITY_X_COMPONENT: rURF=1.0; //1.0;  // уравнения по скорости нелинейны и им нужна нижняя релаксация. Иначе возможна расходимость.
			      break; // эта нижняя релаксация определяется при формировании матрицы СЛАУ.
		case VELOCITY_Y_COMPONENT: rURF=1.0; //1.0;
			      break;
		case VELOCITY_Z_COMPONENT: rURF=1.0; //1.0;
			      break;
		default: rURF=1.0; break; // в остальных случаях.
	}
	// пороговое значение невязки
	doublereal eps = 1e-40;
	doublereal ptilda;
	doublereal sE,sW,sN,sS,sT,sB,sI;
	integer i=0,j=0, kend=10;//100; // Для целей пост сглаживания должно хватить 40 итераций.
	
	if (iVar==PAM) {
		kend=100;
	}
	doublereal dmax=1.0;
	while ((dmax>eps) && (j<kend)) {  
		dmax=0.0;
        //#pragma omp parallel for private(i,ptilda,sE,sW,sN,sS,sT,sB) shared(maxelm,x,rURF,sl) schedule (guided)
	    for (i=0; i<maxelm; ++i) {
            if (sl[i].iE>-1) sE=sl[i].ae*x[sl[i].iE]; else sE=0.0;
		    if (sl[i].iW>-1) sW=sl[i].aw*x[sl[i].iW]; else sW=0.0;
		    if (sl[i].iN>-1) sN=sl[i].an*x[sl[i].iN]; else sN=0.0;
		    if (sl[i].iS>-1) sS=sl[i].as*x[sl[i].iS]; else sS=0.0;
            if (sl[i].iT>-1) sT=sl[i].at*x[sl[i].iT]; else sT=0.0;
		    if (sl[i].iB>-1) sB=sl[i].ab*x[sl[i].iB]; else sB=0.0;

			switch (iVar) {
			case PAM:  if (fabs(sl[i].ap)<1.0e-30) {
#if doubleintprecision == 1
				printf("division by zero in i=%lld internal node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
				printf("division by zero in i=%d internal node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						  
						  printf("Matrix construct PAM error...\n");
						  printf("Please, press any key to halt calculation...\n");
						 // system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=(sE+sW+sN+sS+sT+sB+b[sl[i].iP])/sl[i].ap;
		                  //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			              //dmax=fmax(dmax,fabs(sl[i].ap*(ptilda-x[sl[i].iP]))); // Чебышева
				          dmax+=fabs((ptilda-x[sl[i].iP])); // сумма модулей
		                  x[sl[i].iP]=x[sl[i].iP]+rURF*(ptilda-x[sl[i].iP]);
					  }
					   break;
			case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT: // нижняя релаксация на компоненты скорости.
				      if (fabs(sl[i].ap)<1.0e-30) {
#if doubleintprecision == 1
						  printf("division by zero in i=%lld internal node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						  printf("division by zero in i=%d internal node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						  
						  printf("Matrix construct Velocity error...\n");
						  printf("Please, press any key to halt calculation...\n");
						 // system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=(sE+sW+sN+sS+sT+sB+ b[sl[i].iP])/sl[i].ap;
		                  //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			              //dmax=fmax(dmax,fabs(sl[i].ap*(ptilda-x[sl[i].iP])/alpha)); // Чебышева
					      dmax+=fabs((ptilda-x[sl[i].iP])); // сумма модулей
		                  x[sl[i].iP]=x[sl[i].iP]+rURF*(ptilda-x[sl[i].iP]);
					   }
				       break;
			default: // по умолчанию без нижней реласкации.
				       if (fabs(sl[i].ap)<1.0e-30) {
#if doubleintprecision == 1
						   printf("division by zero in i=%lld boundary node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						   printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						 
						  printf("Matrix construct error...\n");
						  printf("Please, press any key to halt calculation...\n");
						 // system("pause");
						  system("pause");
						  exit(0);
					   }
					   else {
				          ptilda=(sE+sW+sN+sS+sT+sB+ b[sl[i].iP])/sl[i].ap;
		                  //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			              //dmax=fmax(dmax,fabs(sl[i].ap*(ptilda-x[sl[i].iP]))); // Чебышева
					      dmax+=fabs((ptilda-x[sl[i].iP])); // сумма модулей
		                  x[sl[i].iP]=x[sl[i].iP]+rURF*(ptilda-x[sl[i].iP]);
					   }
					   break;
				break;
			}
			
			/*
			if (0&&j==0) {// debug
				printf("ae=%e, aw=%e, an=%e, as=%e, at=%e, ab=%e, ap=%e, b=%e\n",sl[i].ae,sl[i].aw,sl[i].an,sl[i].as,sl[i].at,sl[i].ab,sl[i].ap,sl[i].b);
				system("pause");
			}
			*/
	    }
		//#pragma omp parallel for private(i,ptilda,sI) shared(maxbound,x,rURF,slb) schedule (guided)
		for (i=0; i<maxbound; ++i) {
			if (slb[i].iI>-1) sI=slb[i].ai*x[slb[i].iI]; else sI=0.0;

			switch (iVar) {
			case PAM:
				      if (fabs(slb[i].aw)<1.0e-30) {
#if doubleintprecision == 1
						  printf("division by zero in i=%lld boundary node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						  printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						  
						  printf("Matrix construct PAM error...\n");
						  printf("Please, press any key to halt calculation...\n");
						 // system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=(sI+b[slb[i].iW])/slb[i].aw;
			              //dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[slb[i].iW]))); // Чебышева
					      dmax+=fabs((ptilda-x[slb[i].iW])); // сумма модулей.
			              //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
                          if (slb[i].iI==-1) x[slb[i].iW]=(ptilda);
			              else x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
					  }
					  break;
			case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT:
				      if (fabs(slb[i].aw)<1.0e-30) {
#if doubleintprecision == 1
						  printf("division by zero in i=%lld boundary node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						  printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						 
						  printf("Matrix construct Velocity error...\n");
						  printf("Please, press any key to halt calculation...\n");
						 // system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=(sI+ b[slb[i].iW])/slb[i].aw;
			              //dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[slb[i].iW]))); // Чебышева
					      dmax+=fabs((ptilda-x[slb[i].iW])); // сумма модулей
			              //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
                          if (slb[i].iI==-1) x[slb[i].iW]=(ptilda);
			              else x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
				          /*
				          ptilda=alpha*(sI+slb[i].b+(1.0-alpha)*slb[i].aw*x_cor[slb[i].iW]/alpha)/slb[i].aw;
			              dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[slb[i].iW])/alpha));
			              //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
                          if (slb[i].iI==-1) x[slb[i].iW]=(ptilda);
			              else x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]); 
					      */
					  }
				      break;
			default: // по умолчанию без нижней релаксации.
				      if (fabs(slb[i].aw)<1.0e-30) {
#if doubleintprecision == 1
						  printf("division by zero in i=%lld boundary node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
						  printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
						  
						  printf("Matrix construct error...\n");
						  printf("Please, press any key to halt calculation...\n");
						//  system("pause");
						  system("pause");
						  exit(0);
					  }
					  else {
				          ptilda=(sI+ b[slb[i].iW])/slb[i].aw;
			              //dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[slb[i].iW]))); // Чебышева
					      dmax+=fabs((ptilda-x[slb[i].iW])); // сумма модулей.
			              //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
                          if (slb[i].iI==-1) x[slb[i].iW]=(ptilda);
			              else x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
					  }
				break;
			}
		}
		
		if (bprintMessage) {
			if (iVar == PAM) {
				//dmax/=maxelm;
				//if (j % 1 == 0) // Печатается каждый раз. Если нужно печатать реже нужно окорректировать и раскоментировать условие.
				{
#if doubleintprecision == 1
					printf("%lld %e \n", j + 1, dmax);
#else
					printf("%d %e \n", j + 1, dmax);
#endif

				}
			}
		}
		
		j++;
	}
	
	if (bprintMessage) {
		if (iVar == PAM) {
			printf("calc complete...\n");
			system("pause");
		}
	}
	//printf("4000 %e \n", dmax);
	//system("pause");

} // SOR3Dnow



// применяется для компонент скорости
// Достаточно относительной точности 0.1.
// Полностью параллельный код.
// Метод Якоби без использования параметра нижней релаксации.
// Судя по статьям программного комплекса ЛОГОС для компонент 
// скорости может использоваться метод Якоби т.к. число обусловленности матриц на компоненты скорости примерно 1.0e+1;
// Для сравнения число обусловленности для давления 1e+7 - 1e+8.
// 
// Работает только для Адаптивной локально измельчённой расчётной сетки.
//
void specialSpeed_Alice(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
	doublereal*& x, integer maxelm, integer maxbound, integer iVar,
	bool bprintMessage, doublereal alpharelax) 
{
	// К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
	// скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
	//printf("SOR3D incomming...\n"); // debug.
	//system("pause");

	// Параметр релаксации лежит в интервале от 0.0 до 2.0.
	// Верхняя релаксация способна существенно ускорить вычислительный процесс. 
	// Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
	// rURF=1.5 можно использовать в уравнении на поправку давления.
	// В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
	// к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
	// в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
	// можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
	// коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
	// поэтому наверное лучше не вводить коээффициент 
	// верхней релаксации равный 1.5 в уравнениях на скорость. 

#ifdef _OPENMP 
	omp_set_num_threads(number_cores()); // установка числа потоков
#endif

	// пороговое значение невязки
	doublereal eps = 1e-40;
	doublereal sRF = 1.0;// alpharelax;


	integer j = 0, kend = 10000000;//100; // Для целей пост сглаживания должно хватить 40 итераций.

	const integer size0 = maxelm + maxbound;

	doublereal* x_new = new doublereal[size0];
#pragma omp parallel for
	for (int i = 0; i < size0; ++i)
	{
		x_new[i] = 0.0;
	}

	doublereal* data_b = new doublereal[maxbound * 2];
	int* col_ind_b = new int[maxbound * 2];

#pragma omp parallel for 
	for (int i = 0; i < maxbound; ++i) {

		const int is = i * 2;
		data_b[is] = slb[i].ai;
		data_b[is + 1] = slb[i].aw;

		if (slb[i].iI > -1) {
			col_ind_b[is] = slb[i].iI;
		}
		else {
			data_b[is] = 0.0;
			col_ind_b[is] = slb[i].iW;
		}
		col_ind_b[is + 1] = slb[i].iW;
	}


	doublereal dmax = 1.0;
	while ((dmax > eps) && (j < kend)) {
		dmax = 0.0;

		doublereal dmax_loc = 0.0;

#pragma omp parallel for  reduction(+: dmax_loc)
		for (int i = 0; i < maxelm; ++i)
		{
			doublereal sE, sW, sN, sS, sT, sB;

			if (sl[i].iE > -1) sE = sl[i].ae * x[sl[i].iE]; else sE = 0.0;
			if (sl[i].iW > -1) sW = sl[i].aw * x[sl[i].iW]; else sW = 0.0;
			if (sl[i].iN > -1) sN = sl[i].an * x[sl[i].iN]; else sN = 0.0;
			if (sl[i].iS > -1) sS = sl[i].as * x[sl[i].iS]; else sS = 0.0;
			if (sl[i].iT > -1) sT = sl[i].at * x[sl[i].iT]; else sT = 0.0;
			if (sl[i].iB > -1) sB = sl[i].ab * x[sl[i].iB]; else sB = 0.0;

			doublereal sE2, sW2, sN2, sS2, sT2, sB2;
			
			if (sl[i].iE2 > -1) sE2 = sl[i].ae2 * x[sl[i].iE2]; else sE2 = 0.0;
			if (sl[i].iW2 > -1) sW2 = sl[i].aw2 * x[sl[i].iW2]; else sW2 = 0.0;
			if (sl[i].iN2 > -1) sN2 = sl[i].an2 * x[sl[i].iN2]; else sN2 = 0.0;
			if (sl[i].iS2 > -1) sS2 = sl[i].as2 * x[sl[i].iS2]; else sS2 = 0.0;
			if (sl[i].iT2 > -1) sT2 = sl[i].at2 * x[sl[i].iT2]; else sT2 = 0.0;
			if (sl[i].iB2 > -1) sB2 = sl[i].ab2 * x[sl[i].iB2]; else sB2 = 0.0;

			doublereal sE3, sW3, sN3, sS3, sT3, sB3;

			if (sl[i].iE3 > -1) sE3 = sl[i].ae3 * x[sl[i].iE3]; else sE3 = 0.0;
			if (sl[i].iW3 > -1) sW3 = sl[i].aw3 * x[sl[i].iW3]; else sW3 = 0.0;
			if (sl[i].iN3 > -1) sN3 = sl[i].an3 * x[sl[i].iN3]; else sN3 = 0.0;
			if (sl[i].iS3 > -1) sS3 = sl[i].as3 * x[sl[i].iS3]; else sS3 = 0.0;
			if (sl[i].iT3 > -1) sT3 = sl[i].at3 * x[sl[i].iT3]; else sT3 = 0.0;
			if (sl[i].iB3 > -1) sB3 = sl[i].ab3 * x[sl[i].iB3]; else sB3 = 0.0;

			doublereal sE4, sW4, sN4, sS4, sT4, sB4;

			if (sl[i].iE4 > -1) sE4 = sl[i].ae4 * x[sl[i].iE4]; else sE4 = 0.0;
			if (sl[i].iW4 > -1) sW4 = sl[i].aw4 * x[sl[i].iW4]; else sW4 = 0.0;
			if (sl[i].iN4 > -1) sN4 = sl[i].an4 * x[sl[i].iN4]; else sN4 = 0.0;
			if (sl[i].iS4 > -1) sS4 = sl[i].as4 * x[sl[i].iS4]; else sS4 = 0.0;
			if (sl[i].iT4 > -1) sT4 = sl[i].at4 * x[sl[i].iT4]; else sT4 = 0.0;
			if (sl[i].iB4 > -1) sB4 = sl[i].ab4 * x[sl[i].iB4]; else sB4 = 0.0;

			const int it6 = sl[i].iP;

			const long double ptildaQ = alpharelax*(sE + sW + sN +	sS + sT + sB +
				sE2 + sW2 + sN2 + sS2 + sT2 + sB2 +
				sE3 + sW3 + sN3 + sS3 + sT3 + sB3 +
				sE4 + sW4 + sN4 + sS4 + sT4 + sB4 +
				b[it6]) / sl[i].ap;

			const doublereal ptilda = static_cast<doublereal>(ptildaQ);
			dmax_loc += fabs((ptilda - x[it6])); // сумма модулей
			x_new[it6] = x[it6] + sRF * (ptilda - x[it6]);


			/*
			if (0&&j==0) {// debug
				printf("ae=%e, aw=%e, an=%e, as=%e, at=%e, ab=%e, ap=%e, b=%e\n",sl[i].ae,sl[i].aw,sl[i].an,sl[i].as,sl[i].at,sl[i].ab,sl[i].ap,sl[i].b);
				system("pause");
			}
			*/
		}

		dmax += dmax_loc;
		dmax_loc = 0.0;

#pragma omp parallel for reduction(+: dmax_loc)
		for (int i = 0; i < maxbound; ++i) {

			const int is = i * 2;
			const int it1 = col_ind_b[is];//iI
			const int it2 = col_ind_b[is + 1];//iW 

			const doublereal ptilda = (data_b[is] * x[it1] + b[it2]) / data_b[is + 1];
			dmax_loc += fabs((ptilda - x[slb[i].iW])); // сумма модулей
			if (it1 == it2) x_new[it2] = (ptilda);
			else x_new[it2] = x[it2] + sRF * (ptilda - x[it2]);

		}

		dmax += dmax_loc;

		if (bprintMessage)
		{
			if (iVar == PAM)
			{
				//dmax/=maxelm;
				//if (j % 1 == 0) // Печатаем каждый раз. Если нужно печатать реже надо раскоментировать и подрправить условие.
				{
#if doubleintprecision == 1
					printf("%lld %e \n", j + 1, dmax);
#else
					printf("%d %e \n", j + 1, dmax);
#endif

				}
			}
		}

		if (j == 0) eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.

#pragma omp parallel for
		for (integer i = 0; i < size0; ++i)
		{
			x[i] = x_new[i];
		}


		j++;
	}

	delete[] data_b;
	delete[] col_ind_b;

	delete[] x_new;

	if (bprintMessage) {
		if (iVar == PAM) {
			printf("calc complete...\n");
			system("pause");
		}
	}
	//printf("4000 %e \n", dmax);
	//system("pause");

} // specialSpeed_Alice

// применяется для компонент скорости
// Достаточно относительной точности 0.1.
// Полностью параллельный код.
// Метод Якоби без использования параметра нижней релаксации.
// Судя по статьям программного комплекса ЛОГОС для компонент 
// скорости может использоваться метод Якоби т.к. число обусловленности матриц на компоненты скорости примерно 1.0e+1;
// Для сравнения число обусловленности для давления 1e+7 - 1e+8.
// 
// Работает только для структурированной расчётной сетки.
//
void specialSpeed(equation3D*& sl, equation3D_bon*& slb, doublereal*& b,
	doublereal*& x, integer maxelm, integer maxbound, integer iVar,
	bool bprintMessage, doublereal alpharelax) {
	// К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
	// скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
	//printf("SOR3D incomming...\n"); // debug.
	//system("pause");

	// Параметр релаксации лежит в интервале от 0.0 до 2.0.
	// Верхняя релаксация способна существенно ускорить вычислительный процесс. 
	// Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
	// rURF=1.5 можно использовать в уравнении на поправку давления.
	// В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
	// к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
	// в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
	// можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
	// коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость),
	// поэтому наверное лучше не вводить коээффициент 
	// верхней релаксации равный 1.5 в уравнениях на скорость. 
	
#ifdef _OPENMP 
	omp_set_num_threads(number_cores()); // установка числа потоков
#endif

	// пороговое значение невязки
	doublereal eps = 1e-40;
	doublereal sRF = 1.0;// alpharelax;
	

	integer j = 0, kend = 10000000;//100; // Для целей пост сглаживания должно хватить 40 итераций.

	const integer size0 = maxelm + maxbound;

	doublereal* x_new = new doublereal[size0];
#pragma omp parallel for
	for (int i = 0; i < size0; ++i)
	{
		x_new[i] = 0.0;
	}

	// Улучшает попадание в кеш ускорение на 12.5%.
	doublereal* data = new doublereal[maxelm * 7];
	int* col_ind = new int[maxelm * 7];

#pragma omp parallel for 
	for (int i = 0; i < maxelm; ++i) {

		const int is = i * 7;
		data[is] = sl[i].ae;
		data[is + 1] = sl[i].aw;
		data[is + 2] = sl[i].an;
		data[is + 3] = sl[i].as;
		data[is + 4] = sl[i].at;
		data[is + 5] = sl[i].ab;
		data[is + 6] = sl[i].ap / alpharelax;

		col_ind[is] = sl[i].iE;
		col_ind[is + 1] = sl[i].iW;
		col_ind[is + 2] = sl[i].iN;
		col_ind[is + 3] = sl[i].iS;
		col_ind[is + 4] = sl[i].iT;
		col_ind[is + 5] = sl[i].iB;
		col_ind[is + 6] = sl[i].iP;
	}


	doublereal* data_b = new doublereal[maxbound * 2];
	int* col_ind_b = new int[maxbound * 2];

#pragma omp parallel for 
	for (int i = 0; i < maxbound; ++i) {

		const int is = i * 2;
		data_b[is] = slb[i].ai;
		data_b[is + 1] = slb[i].aw;

		if (slb[i].iI > -1) {
			col_ind_b[is] = slb[i].iI;
		}
		else {
			data_b[is] = 0.0;
			col_ind_b[is] = slb[i].iW;
		}
		col_ind_b[is + 1] = slb[i].iW;
	}


	doublereal dmax_init = 1.0;

	doublereal dmax = 1.0;
	while ((dmax > eps) && (j < kend)) {
		dmax = 0.0;
		
		doublereal dmax_loc = 0.0;

#pragma omp parallel for  reduction(+: dmax_loc)
		for (int i = 0; i < maxelm; ++i)
		{			
			
				const int is = i * 7;
				const int is1 = is + 1;
				const int is2 = is + 2;
				const int is3 = is + 3;
				const int is4 = is + 4;
				const int is5 = is + 5;
				const int is6 = is + 6;
				const int it6 = col_ind[is6];

				const long double ptildaQ = (data[is] * x[col_ind[is]] +
					data[is1] * x[col_ind[is1]] +
					data[is2] * x[col_ind[is2]] +
					data[is3] * x[col_ind[is3]] +
					data[is4] * x[col_ind[is4]] +
					data[is5] * x[col_ind[is5]] +
					b[it6]) / data[is6];

				const doublereal ptilda = static_cast<doublereal>(ptildaQ);
				dmax_loc += fabs((ptilda - x[it6])); // сумма модулей
				x_new[it6] = x[it6] + sRF * (ptilda - x[it6]);

             
			/*
			if (0&&j==0) {// debug
				printf("ae=%e, aw=%e, an=%e, as=%e, at=%e, ab=%e, ap=%e, b=%e\n",sl[i].ae,sl[i].aw,sl[i].an,sl[i].as,sl[i].at,sl[i].ab,sl[i].ap,sl[i].b);
				system("pause");
			}
			*/
		}

		dmax += dmax_loc;
		dmax_loc = 0.0;

#pragma omp parallel for reduction(+: dmax_loc)
		for (int i = 0; i < maxbound; ++i) {
			
			const int is = i * 2;
			const int it1 = col_ind_b[is];//iI
			const int it2 = col_ind_b[is + 1];//iW 

			const doublereal ptilda = (data_b[is] * x[it1] + b[it2]) / data_b[is+1];
			dmax_loc += fabs((ptilda - x[it2])); // сумма модулей
			if (it1 == it2) x_new[it2] = (ptilda);
			else x_new[it2] = x[it2] + sRF * (ptilda - x[it2]);
							
		}

		dmax += dmax_loc;

		if (bprintMessage)
		{
			if (iVar == PAM)
			{
				//dmax/=maxelm;
				//if (j % 1 == 0) // Печатаем каждый раз. Если нужно печатать реже надо раскоментировать и подрправить условие.
				{
#if doubleintprecision == 1
					printf("%lld %e \n", j + 1, dmax);
#else
					printf("%d %e \n", j + 1, dmax);
#endif

				}
			}
		}

		if (j == 0) {
			dmax_init = dmax;
			eps = 0.1 * dmax; // Важнейшее условие выхода из итерационного процесса.

		}

#pragma omp parallel for
		for (integer i = 0; i < size0; ++i)
		{
			x[i] = x_new[i];
		}


		j++;
	}


	std::cout << " " << j << " " << dmax_init << " " << dmax <<  " \n";

	delete[] data;
	delete[] col_ind;

	delete[] data_b;
	delete[] col_ind_b;

	delete[] x_new;

	if (bprintMessage) {
		if (iVar == PAM) {
			printf("calc complete...\n");
			system("pause");
		}
	}
	//printf("4000 %e \n", dmax);
	//system("pause");

} // specialSpeed


// Этот алгоритм используется как составляющая алгоритма Ван дер Ворса.
// применяется для уравнения поправки давления 
// в случае когда на всей границе стоят условия Неймана.
// Скорость сходимости очень и очень медленная,
// поэтому этот метод используется НИКОГДА.
// Основываясь на идеях алгоритма Федоренко этот метод можно применять как сглаживатель.
// Версия для переупорядочивания по методу nested desection.
void PAMGSPnd(equation3D* &sl, equation3D_bon* &slb, doublereal* &x, doublereal* &rthdsd, int maxelm, int maxbound, int* &ifrontregulationgl) {
	// К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
	// скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
	//printf("SOR3D incomming...\n"); // debug.
	//system("pause");

	// Параметр релаксации лежит в интервале от 0.0 до 2.0.
	// Верхняя релаксация способна существенно ускорить вычислительный процесс. 
	// Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
	// rURF=1.5 можно использовать в уравнении на поправку давления.
	// В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
	// к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
	// в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
	// можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
	// коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость)
	//, поэтому наверное лучше не вводить коээффициент 
	// верхней релаксации равный 1.5 в уравнениях на скорость. 

    doublereal rURF=1.0; // параметр верхней релаксации
	
	// пороговое значение невязки
	doublereal eps = 1e-40;
	
	
	integer j=0, kend=500;//100; // Для целей пост сглаживания должно хватить 40 итераций.

	//doublereal sigma1=0.0, sigma2=0.0, sigma;
	/*
	if (iVar==PAM) {
		kend=10000;
	}*/
	doublereal dmax=1.0, dmaxl;
	while ((dmax>eps) && (j<kend)) {  
		dmax=0.0;
		dmaxl=0.0;
		//sigma=0.0;
		// Нельзя использовать такое примитивное распараллеливание оно приводит к расходимости,
		// т.к. разрешающие способности параллельного метода Зейделя стали хуже.
        //#pragma omp parallel for shared(maxelm,x,rURF,sl,rthdsd,j) reduction(+: dmaxl/*, sigma*/)  schedule (static)
	    for (integer i=0; i<maxelm; ++i) {
			doublereal sE,sW,sN,sS,sT,sB;
			doublereal ptilda;
			// по iP получить i.
            if (sl[i].iE>-1) sE=sl[i].ae*x[ifrontregulationgl[sl[i].iE]]; else sE=0.0;
		    if (sl[i].iW>-1) sW=sl[i].aw*x[ifrontregulationgl[sl[i].iW]]; else sW=0.0;
		    if (sl[i].iN>-1) sN=sl[i].an*x[ifrontregulationgl[sl[i].iN]]; else sN=0.0;
		    if (sl[i].iS>-1) sS=sl[i].as*x[ifrontregulationgl[sl[i].iS]]; else sS=0.0;
            if (sl[i].iT>-1) sT=sl[i].at*x[ifrontregulationgl[sl[i].iT]]; else sT=0.0;
		    if (sl[i].iB>-1) sB=sl[i].ab*x[ifrontregulationgl[sl[i].iB]]; else sB=0.0;

			
			if (fabs(sl[i].ap)<1.0e-30) {
#if doubleintprecision == 1
				printf("division by zero in i=%lld internal node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#else
				printf("division by zero in i=%d internal node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
			   
				printf("Matrix construct PAM error...\n");
				printf("Please, press any key to halt calculation...\n");
			//	system("pause");
				system("pause");
				exit(0);
			}
			else {
			    ptilda=(sE+sW+sN+sS+sT+sB+rthdsd[ifrontregulationgl[i]])/sl[i].ap;
		        //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			    //dmax=fmax(dmax,fabs(sl[i].ap*(ptilda-x[sl[i].iP]))); // Чебышева
			    dmaxl+=fabs((ptilda-x[ifrontregulationgl[sl[i].iP]])); // сумма модулей
			    // if ((j==1)||(j==2)) sigma+=(ptilda-x[sl[i].iP]);
				doublereal xbuf=x[ifrontregulationgl[sl[i].iP]]+rURF*(ptilda-x[ifrontregulationgl[sl[i].iP]]);
		        x[ifrontregulationgl[sl[i].iP]]=xbuf;
			}
					  
			
			
			
			/*
			if (0&&j==0) {// debug
				printf("ae=%e, aw=%e, an=%e, as=%e, at=%e, ab=%e, ap=%e, b=%e\n",sl[i].ae,sl[i].aw,sl[i].an,sl[i].as,sl[i].at,sl[i].ab,sl[i].ap,sl[i].b);
				system("pause");
			}
			*/
	    }

		//if (j==1) sigma1+=sigma;
		//if (j==2) sigma2+=sigma;
		dmax+=dmaxl;

		//sigma=0.0;
		dmaxl=0.0;

		//#pragma omp parallel for shared(maxbound,maxelm,x,rURF,slb, rthdsd, j) reduction(+: dmaxl/*, sigma*/) schedule (static)
		for (integer i=0; i<maxbound; ++i) {

			doublereal sI;
			if (slb[i].iI>-1) sI=slb[i].ai*x[ifrontregulationgl[slb[i].iI]]; else sI=0.0;

			
			if (fabs(slb[i].aw)<1.0e-30) {
#if doubleintprecision == 1
				printf("division by zero in i=%lld boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#else
				printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
			   
			   printf("Matrix construct PAM error...\n");
			   printf("Please, press any key to halt calculation...\n");
			  // system("pause");
			   system("pause");
			   exit(0);
			}
			else {
			   doublereal ptilda;
			   ptilda=(sI+rthdsd[ifrontregulationgl[i+maxelm]])/slb[i].aw;
			   //dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[ifrontregulationgl[slb[i].iW]]))); // Чебышева
			   dmaxl+=fabs((ptilda-x[ifrontregulationgl[slb[i].iW]])); // сумма модулей.
			   //x[ifrontregulationgl[slb[i].iW]]=x[ifrontregulationgl[slb[i].iW]]+rURF*(ptilda-x[ifrontregulationgl[slb[i].iW]]);
               if (slb[i].iI==-1) {
				   // if ((j==1)||(j==2)) sigma+=(ptilda-x[slb[i].iW]);
				   x[ifrontregulationgl[slb[i].iW]]=(ptilda);
			  }
			  else {
			       //if ((j==1)||(j==2)) sigma+=(ptilda-x[slb[i].iW]);
				   //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
				   x[ifrontregulationgl[slb[i].iW]]=ptilda;
			  }
			}
		}


		//if (j==1) sigma1+=sigma;
		//if (j==2) sigma2+=sigma;	

		dmax+=dmaxl;

#if doubleintprecision == 1
		printf("%lld %e\n", j, dmax);
#else
		printf("%d %e\n", j, dmax);
#endif
		
		//system("pause");
		
		/*
		if (iVar==PAM)  {
		  //dmax/=maxelm;
		  if (j%1==0) {
		  #if doubleintprecision == 1
				 printf("%lld %e \n", j+1, dmax);
		  #else
				 printf("%d %e \n", j+1, dmax);
		  #endif
			  
		  }
		}
		*/
		
		j++;

		 /*// при использовании формулы Ильина возникают арифметические переполнения.
		if (j==3) {
			if (fabs(sigma1)>1.0e-20) {
			   rURF=2.0/(1.0+sqrt(fmax(1.0-sqrt(sigma2/sigma1),0.0)));
			   if (fabs(rURF-2.0)<1.0e-20) { rURF=1.0; }
			}
		}
		*/
		
		
	}
	
	
	
	
	//printf("4000 %e \n", dmax);
	//system("pause");

} // PAMGSPnd


// Этот алгоритм используется как составляющая алгоритма Ван дер Ворса.
// применяется для уравнения поправки давления 
// в случае когда на всей границе стоят условия Неймана.
// Скорость сходимости очень и очень медленная,
// поэтому этот метод используется НИКОГДА.
// Основываясь на идеях алгоритма Федоренко этот метод можно применять как сглаживатель.
void PAMGSP(equation3D* &sl, equation3D_bon* &slb, doublereal* &x, doublereal* &rthdsd, integer maxelm, integer maxbound) {
	// К величине xcor - будет осуществляться нижняя релаксация, т.е. x_cor - это
	// скорректированная компонента скорости удовлетворяющая уравнению неразрывности.
	//printf("SOR3D incomming...\n"); // debug.
	//system("pause");

	// Параметр релаксации лежит в интервале от 0.0 до 2.0.
	// Верхняя релаксация способна существенно ускорить вычислительный процесс. 
	// Патанкар рекомендует брать коэффициент верхней релаксации равный 1.5;
	// rURF=1.5 можно использовать в уравнении на поправку давления.
	// В уравнениях на скорость мы имеем двоякую релаксацию. С одной стороны это нижняя релаксация
	// к скорректированной скорости с коэффициентом alpha. Эта нижняя релаксация учитывается в матрице СЛАУ,
	// в результате чего мы получаем новую матрицу СЛАУ.  К этой новой преоблразованной матрице СЛАУ казалось бы в целях ускорения сходимости
	// можно применить верхнюю релаксацию, пусть она будет опять с коэффициентом rURF=1.5. Вычисления показывают что при введении в уравнения на скорость
	// коэффициента верхней релаксации 1.5 точность решения уравнений на скорость за 4000 итераций падает (или вообще имеем расходимость)
	//, поэтому наверное лучше не вводить коээффициент 
	// верхней релаксации равный 1.5 в уравнениях на скорость. 

    doublereal rURF=1.0; // параметр верхней релаксации
	
	// пороговое значение невязки
	doublereal eps = 1e-40;
	
	
	integer j=0, kend=500;//100; // Для целей пост сглаживания должно хватить 40 итераций.

	//doublereal sigma1=0.0, sigma2=0.0, sigma;
	/*
	if (iVar==PAM) {
		kend=10000;
	}*/
	doublereal dmax=1.0, dmaxl;
	while ((dmax>eps) && (j<kend)) {  
		dmax=0.0;
		dmaxl=0.0;
		//sigma=0.0;
		// Нельзя использовать такое примитивное распараллеливание оно приводит к расходимости,
		// т.к. разрешающие способности параллельного метода Зейделя стали хуже.
        //#pragma omp parallel for shared(maxelm,x,rURF,sl,rthdsd,j) reduction(+: dmaxl/*, sigma*/)  schedule (static)
	    for (integer i=0; i<maxelm; ++i) {
			doublereal sE,sW,sN,sS,sT,sB;
			doublereal ptilda;
            if (sl[i].iE>-1) sE=sl[i].ae*x[sl[i].iE]; else sE=0.0;
		    if (sl[i].iW>-1) sW=sl[i].aw*x[sl[i].iW]; else sW=0.0;
		    if (sl[i].iN>-1) sN=sl[i].an*x[sl[i].iN]; else sN=0.0;
		    if (sl[i].iS>-1) sS=sl[i].as*x[sl[i].iS]; else sS=0.0;
            if (sl[i].iT>-1) sT=sl[i].at*x[sl[i].iT]; else sT=0.0;
		    if (sl[i].iB>-1) sB=sl[i].ab*x[sl[i].iB]; else sB=0.0;

			
			if (fabs(sl[i].ap)<1.0e-30) {
#if doubleintprecision == 1
				printf("division by zero in i=%lld internal node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
				printf("division by zero in i=%d internal node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
			    
				printf("Matrix construct PAM error...\n");
				printf("Please, press any key to halt calculation...\n");
				//system("pause");
				system("pause");
				exit(0);
			}
			else {
			    ptilda=(sE+sW+sN+sS+sT+sB+rthdsd[i])/sl[i].ap;
		        //dmax=fmax(dmax,sl[i].ap*(ptilda-x[sl[i].iP]));
			    //dmax=fmax(dmax,fabs(sl[i].ap*(ptilda-x[sl[i].iP]))); // Чебышева
			    dmaxl+=fabs((ptilda-x[sl[i].iP])); // сумма модулей
			    // if ((j==1)||(j==2)) sigma+=(ptilda-x[sl[i].iP]);
				doublereal xbuf=x[sl[i].iP]+rURF*(ptilda-x[sl[i].iP]);
		        x[sl[i].iP]=xbuf;
			}
					  
			
			
			
			/*
			if (0&&j==0) {// debug
				printf("ae=%e, aw=%e, an=%e, as=%e, at=%e, ab=%e, ap=%e, b=%e\n",sl[i].ae,sl[i].aw,sl[i].an,sl[i].as,sl[i].at,sl[i].ab,sl[i].ap,sl[i].b);
				system("pause");
			}
			*/
	    }

		//if (j==1) sigma1+=sigma;
		//if (j==2) sigma2+=sigma;
		dmax+=dmaxl;

		//sigma=0.0;
		dmaxl=0.0;

		//#pragma omp parallel for shared(maxbound,maxelm,x,rURF,slb, rthdsd, j) reduction(+: dmaxl/*, sigma*/) schedule (static)
		for (integer i=0; i<maxbound; ++i) {

			doublereal sI;
			if (slb[i].iI>-1) sI=slb[i].ai*x[slb[i].iI]; else sI=0.0;

			
			if (fabs(slb[i].aw)<1.0e-30) {
#if doubleintprecision == 1
				printf("division by zero in i=%lld boundary node. maxelm=%lld, maxbound=%lld\n", i, maxelm, maxbound);
#else
				printf("division by zero in i=%d boundary node. maxelm=%d, maxbound=%d\n", i, maxelm, maxbound);
#endif
			  
			   printf("Matrix construct PAM error...\n");
			   printf("Please, press any key to halt calculation...\n");
			//   system("pause");
			   system("pause");
			   exit(0);
			}
			else {
			   doublereal ptilda;
			   ptilda=(sI+rthdsd[i+maxelm])/slb[i].aw;
			   //dmax=fmax(dmax,fabs(slb[i].aw*(ptilda-x[slb[i].iW]))); // Чебышева
			   dmaxl+=fabs((ptilda-x[slb[i].iW])); // сумма модулей.
			   //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
               if (slb[i].iI==-1) {
				   // if ((j==1)||(j==2)) sigma+=(ptilda-x[slb[i].iW]);
				   x[slb[i].iW]=(ptilda);
			  }
			  else {
			       //if ((j==1)||(j==2)) sigma+=(ptilda-x[slb[i].iW]);
				   //x[slb[i].iW]=x[slb[i].iW]+rURF*(ptilda-x[slb[i].iW]);
				   x[slb[i].iW]=ptilda;
			  }
			}
		}


		//if (j==1) sigma1+=sigma;
		//if (j==2) sigma2+=sigma;	

		dmax+=dmaxl;

#if doubleintprecision == 1
		printf("%lld %e\n", j, dmax);
#else
		printf("%d %e\n", j, dmax);
#endif
		
		//system("pause");
		
		/*
		if (iVar==PAM)  {
		  //dmax/=maxelm;
		  if (j%1==0) {
		  #if doubleintprecision == 1
				 printf("%lld %e \n", j+1, dmax);
		  #else
				 printf("%d %e \n", j+1, dmax);
		  #endif
			  
		  }
		}
		*/
		
		j++;

		 /*// при использовании формулы Ильина возникают арифметические переполнения.
		if (j==3) {
			if (fabs(sigma1)>1.0e-20) {
			   rURF=2.0/(1.0+sqrt(fmax(1.0-sqrt(sigma2/sigma1),0.0)));
			   if (fabs(rURF-2.0)<1.0e-20) { rURF=1.0; }
			}
		}
		*/
		
		
	}
	
	
	
	
	//printf("4000 %e \n", dmax);
	//system("pause");

} // PAMGSP

// Для сошедшейся задачи подходит уровень среднеквадратических невязок 1.0e-4
// Источник информации мануал по CFX на русском.
// Если норма Чебышёва то невязка сошедшаяся равна 1.0e-3.
// Источник опять же мануал по CFX.

/*
// Евклидова норма вектора
// данный код уже реализован в файле my_LR.c
doublereal NormaV(double *V, integer n) {
    doublereal r=0.0;
	doublereal dsize=static_cast<doublereal>(1.0*n);

    #pragma omp parallel for shared(V,dsize) schedule (guided) reduction (+:s)
    for (integer i=0; i<n; ++i) {
        r+=V[i]*V[i]/dsize;
	}

    return r;
}
*/

// Норма Чебышева для корректного сравнения с другими методами
doublereal NormaChebyshev(doublereal *V, integer n){
	doublereal norma=-1.0;
	integer i=0;
	for (i=0; i<n; ++i) if (fabs(V[i])>norma) norma=fabs(V[i]);
	return norma;
} // NormaChebyshev 

// применяется для уравнения поправки давления и в уравнениях на скорость.
// в случае когда на всей границе стоят условия Неймана.
// Здесь приводится реализация BT-процесса с управлением
// так как он описан в книге Фадеев-Фадеева. Данный процесс
// основывается на итерационном процессе SOR3D.
// данный метод также можно назвать методом Абрамова.
void BTrules(equation3D* &sl, equation3D_bon* &slb, doublereal* &x, doublereal* x_cor, integer maxelm, integer maxbound, integer iVar, doublereal alpha) {
    doublereal rURF=1.0; // параметр верхней релаксации
	switch (iVar) {
		case PAM: rURF=1.0; //1.855; 
			      break;
		case VELOCITY_X_COMPONENT: rURF=1.0;  // уравнения по скорости нелинейны и им нужна нижняя релаксация. Иначе возможна расходимость.
			      break; // эта нижняя релаксация определяется при формировании матрицы СЛАУ.
		case VELOCITY_Y_COMPONENT: rURF=1.0;
			      break;
		case VELOCITY_Z_COMPONENT: rURF=1.0;
			      break;
		default: rURF=1.0; break; // в остальных случаях.
	}
	// пороговое значение невязки
	doublereal eps = 1.0e-36;
	doublereal ptilda=0.0;
	doublereal sE=0.0,sW=0.0,sN=0.0,sS=0.0,sT=0.0,sB=0.0,sI=0.0;
	integer i=0,j=0, kend=4000;//100; // Для целей пост сглаживания должно хватить 40 итераций.
	/*
	if (iVar==PAM) {
		kend=10000;
	}
	*/

	doublereal** xarg=nullptr; // решение на двух последовательных итерациях.
	doublereal** resarg=nullptr; // невязка на двух последовательных итерациях.
	doublereal** sarg=nullptr; // вспомогательный вектор на двух соседних итерациях.

	xarg=new doublereal*[3];
	resarg=new doublereal*[3];
	sarg=new doublereal*[2];

	for (i=0; i<3; ++i) {
		xarg[i]=new doublereal[maxelm+maxbound];
		resarg[i]=new doublereal[maxelm+maxbound];
		if (i<2) {
			sarg[i]=new doublereal[maxelm+maxbound];
		}
	}

	for (i=0; i<maxelm+maxbound; ++i) {
		xarg[0][i]=x[i];
	}
	doublereal dmax=1.0;

	for (i=0; i<maxelm; ++i) {
        if (sl[i].iE>-1) sE=sl[i].ae*xarg[0][sl[i].iE]; else sE=0.0;
	    if (sl[i].iW>-1) sW=sl[i].aw*xarg[0][sl[i].iW]; else sW=0.0;
	    if (sl[i].iN>-1) sN=sl[i].an*xarg[0][sl[i].iN]; else sN=0.0;
	    if (sl[i].iS>-1) sS=sl[i].as*xarg[0][sl[i].iS]; else sS=0.0;
        if (sl[i].iT>-1) sT=sl[i].at*xarg[0][sl[i].iT]; else sT=0.0;
	    if (sl[i].iB>-1) sB=sl[i].ab*xarg[0][sl[i].iB]; else sB=0.0;
	    

		switch (iVar) {
		  case PAM: ptilda=(sE+sW+sN+sS+sT+sB+sl[i].b)/sl[i].ap; 
			       break;
		  case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT:
			      ptilda=alpha*(sE+sW+sN+sS+sT+sB+sl[i].b+(1.0-alpha)*sl[i].ap*x_cor[sl[i].iP]/alpha)/sl[i].ap;
			  break;
		  default:
			   ptilda=(sE+sW+sN+sS+sT+sB+sl[i].b)/sl[i].ap; 
			  break;
		}

	    resarg[0][sl[i].iP]=(ptilda-xarg[0][sl[i].iP]);// невязка
	    xarg[1][sl[i].iP]=xarg[0][sl[i].iP]+rURF*(ptilda-xarg[0][sl[i].iP]); // новое приближение.
	}
	for (i=0; i<maxbound; ++i) {
		   // граничные условия не стремятся к скорректированной скорости.
			if (slb[i].iI>-1) sI=slb[i].ai*xarg[0][slb[i].iI]; else sI=0.0;
			ptilda=(sI+slb[i].b)/slb[i].aw;
			
            if (slb[i].iI==-1) {
				resarg[0][slb[i].iW]=(ptilda)-xarg[0][slb[i].iW];
				xarg[1][slb[i].iW]=(ptilda);
			}
			else {
				resarg[0][slb[i].iW]=(ptilda-xarg[0][slb[i].iW]);
				xarg[1][slb[i].iW]=xarg[0][slb[i].iW]+rURF*(ptilda-xarg[0][slb[i].iW]);
			}
		}

	for (i=0; i<maxelm; ++i) {
        if (sl[i].iE>-1) sE=sl[i].ae*xarg[1][sl[i].iE]; else sE=0.0;
	    if (sl[i].iW>-1) sW=sl[i].aw*xarg[1][sl[i].iW]; else sW=0.0;
	    if (sl[i].iN>-1) sN=sl[i].an*xarg[1][sl[i].iN]; else sN=0.0;
	    if (sl[i].iS>-1) sS=sl[i].as*xarg[1][sl[i].iS]; else sS=0.0;
        if (sl[i].iT>-1) sT=sl[i].at*xarg[1][sl[i].iT]; else sT=0.0;
	    if (sl[i].iB>-1) sB=sl[i].ab*xarg[1][sl[i].iB]; else sB=0.0;
	    

		switch (iVar) {
		case PAM: ptilda=(sE+sW+sN+sS+sT+sB+sl[i].b)/sl[i].ap;
			    break;
		case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT:
			 ptilda=alpha*(sE+sW+sN+sS+sT+sB+sl[i].b+(1.0-alpha)*sl[i].ap*x_cor[sl[i].iP]/alpha)/sl[i].ap;   
			break;
		default:
			ptilda=(sE+sW+sN+sS+sT+sB+sl[i].b)/sl[i].ap;
			break;
		}

	    resarg[1][sl[i].iP]=(ptilda-xarg[1][sl[i].iP]);// невязка
	}
	// граничные условия заданы наперёд и не релаксируют к скорректированной скорости.
	for (i=0; i<maxbound; ++i) {
		if (slb[i].iI>-1) sI=slb[i].ai*xarg[1][slb[i].iI]; else sI=0.0;
		ptilda=(sI+slb[i].b)/slb[i].aw;
			
        if (slb[i].iI==-1) resarg[1][slb[i].iW]=(ptilda)-xarg[1][slb[i].iW];
	    else resarg[1][slb[i].iW]=(ptilda-xarg[1][slb[i].iW]);
    }
		


	dmax=1.0;
	while ((dmax>eps) && (j<kend)) {
		dmax=0.0;
        //#pragma omp parallel for private(i,ptilda,sE,sW,sN,sS,sT,sB) shared(maxelm,x,rURF,sl) schedule (guided)
	    for (i=0; i<maxelm; ++i) {
            if (sl[i].iE>-1) sE=sl[i].ae*resarg[1][sl[i].iE]; else sE=0.0;
		    if (sl[i].iW>-1) sW=sl[i].aw*resarg[1][sl[i].iW]; else sW=0.0;
		    if (sl[i].iN>-1) sN=sl[i].an*resarg[1][sl[i].iN]; else sN=0.0;
		    if (sl[i].iS>-1) sS=sl[i].as*resarg[1][sl[i].iS]; else sS=0.0;
            if (sl[i].iT>-1) sT=sl[i].at*resarg[1][sl[i].iT]; else sT=0.0;
		    if (sl[i].iB>-1) sB=sl[i].ab*resarg[1][sl[i].iB]; else sB=0.0;

			switch (iVar) {
		       case PAM: ptilda=(sE+sW+sN+sS+sT+sB)/sl[i].ap;
			             break;
		       case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT:
			             ptilda=alpha*(sE+sW+sN+sS+sT+sB)/sl[i].ap;   
			             break;
		       default:
			             ptilda=(sE+sW+sN+sS+sT+sB)/sl[i].ap;
			             break;
		    }
		    sarg[0][sl[i].iP]=ptilda;
			sarg[1][sl[i].iP]=2*ptilda-resarg[0][sl[i].iP];			
	    }
		//#pragma omp parallel for private(i,ptilda,sI) shared(maxbound,x,rURF,slb) schedule (guided)
		for (i=0; i<maxbound; ++i) {
			if (slb[i].iI>-1) sI=slb[i].ai*resarg[1][slb[i].iI]; else sI=0.0;
			ptilda=(sI)/slb[i].aw;
			sarg[0][slb[i].iW]=ptilda;
			sarg[1][slb[i].iW]=2*ptilda-resarg[0][slb[i].iW];
		}

		if (fabs(NormaChebyshev(sarg[0],maxelm+maxbound))<=fabs(NormaChebyshev(sarg[1],maxelm+maxbound))) {
			for (i=0; i<maxelm+maxbound; ++i) {
			     xarg[2][i]=xarg[1][i]+resarg[1][i];
				 resarg[2][i]=sarg[0][i];
			}
		}
		else {
			for (i=0; i<maxelm+maxbound; ++i) {
				xarg[2][i]=2.0*(xarg[1][i]+resarg[1][i])-xarg[0][i];
				resarg[2][i]=sarg[1][i];
			}
		}
		
		// Сдвиг влево для обновления вычислений.
		for (i=0; i<maxelm+maxbound; ++i) {
			xarg[0][i]=xarg[1][i];
			xarg[1][i]=xarg[2][i];
			resarg[0][i]=resarg[1][i];
			resarg[1][i]=resarg[2][i];
		}

		// Раз в пять итераций чтобы избежать накопления ошибок рекомендуется
		// вычислять невязку точно, а не по рекурентным формулам.
		if (j%5==0) {
			for (i=0; i<maxelm; ++i) {
                if (sl[i].iE>-1) sE=sl[i].ae*xarg[1][sl[i].iE]; else sE=0.0;
	            if (sl[i].iW>-1) sW=sl[i].aw*xarg[1][sl[i].iW]; else sW=0.0;
	            if (sl[i].iN>-1) sN=sl[i].an*xarg[1][sl[i].iN]; else sN=0.0;
	            if (sl[i].iS>-1) sS=sl[i].as*xarg[1][sl[i].iS]; else sS=0.0;
                if (sl[i].iT>-1) sT=sl[i].at*xarg[1][sl[i].iT]; else sT=0.0;
	            if (sl[i].iB>-1) sB=sl[i].ab*xarg[1][sl[i].iB]; else sB=0.0;

				switch (iVar) {
		           case PAM: ptilda=(sE+sW+sN+sS+sT+sB+sl[i].b)/sl[i].ap; 
			                  break;
		           case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT: case VELOCITY_Z_COMPONENT:
			                  ptilda=alpha*(sE+sW+sN+sS+sT+sB+sl[i].b+(1.0-alpha)*sl[i].ap*x_cor[sl[i].iP]/alpha)/sl[i].ap;
			                  break;
		           default:
			                  ptilda=(sE+sW+sN+sS+sT+sB+sl[i].b)/sl[i].ap; 
			                 break;
		        }

	            resarg[1][sl[i].iP]=(ptilda-xarg[1][sl[i].iP]);// невязка
	        }
	        for (i=0; i<maxbound; ++i) {
			    if (slb[i].iI>-1) sI=slb[i].ai*xarg[1][slb[i].iI]; else sI=0.0;
			    ptilda=(sI+slb[i].b)/slb[i].aw;
			
                if (slb[i].iI==-1) resarg[1][slb[i].iW]=(ptilda)-xarg[1][slb[i].iW];
			    else resarg[1][slb[i].iW]=(ptilda-xarg[1][slb[i].iW]);
		    }
		}

		dmax=NormaChebyshev(resarg[1],maxelm+maxbound);

		/*
		if (iVar==PAM)  {
		  //dmax/=maxelm;
		  if (j%200==0) {
		  #if doubleintprecision == 1
				printf("%lld %e \n", j+1, dmax);
		  #else
				printf("%d %e \n", j+1, dmax);
		  #endif
			  
		  }
		}
		*/
		j++;
	}
	
	/*
	if (iVar==PAM)  {
	   printf("calc complete...\n");
       system("pause");
	}
	*/

	for (i=0; i<maxelm+maxbound; ++i) {
		x[i]=xarg[1][i]; // финишное копирование.
	}

	// освобождение оперативной памяти.
	if (xarg != nullptr) {
		for (i = 0; i < 3; ++i) {
			if (xarg[i] != nullptr) {
				delete xarg[i];
				xarg[i] = nullptr;
			}
		}
		delete[] xarg;
		xarg = nullptr;
	}
	if (resarg != nullptr) {
		for (i = 0; i < 3; ++i) {
			if (resarg[i] != nullptr) {
				delete resarg[i];
				resarg[i] = nullptr;
			}
		}
		delete[] resarg;
		resarg = nullptr;
	}
	if (sarg != nullptr) {
		for (i = 0; i < 2; ++i) {
	        if (sarg[i]!=nullptr) {
		        delete[] sarg[i];
				sarg[i] = nullptr;
	        }
		}
		delete[] sarg;
		sarg = nullptr;
	}


} // BTrules

/* Метод Сопряжённых градиентов
*  без учёта разреженности матрицы СЛАУ.
*/

// умножение матрицы на вектор
doublereal* MatrixByVector(doublereal** H,doublereal* V,integer n){
	doublereal* tmp=new doublereal[n];
	doublereal sum=0.0;
	for (integer i=0;i<n;++i){
		for (integer j=0;j<n;++j)
			sum+=V[j]*H[i][j];
		tmp[i]=sum;
		sum=0.0;}
	return tmp;
} // MatrixByVector



// Евклидова норма вектора
// Отладочный вариант.
doublereal NormaVdebug(doublereal *V, integer n){
	doublereal norma;
	doublereal s=0.0;
	//#pragma omp parallel for shared(V) schedule (guided) reduction (+:s)
	for (integer i=0;i<n;++i) {
		s+=V[i]*V[i];
		if (!(V[i]==V[i])) {
#if doubleintprecision == 1
			printf("bitji vector i=%lld\n", i);
#else
			printf("bitji vector i=%d\n", i);
#endif
			
			printf("%e ",V[i]);
			//system("pause");
			system("pause");
		}
		//if (i%200==0) system("pause");
	}
	printf("%e\n",s);
	if (!(s==s)) {
		// Это NaN (Not a Number)
		// Предположительно NaN возникает при перемножении двух малых чисел
		// таких что получается настолько малое число что оно не умещается в модели вещественных чисел.
		// Это объяснение несостоятельно.
		norma=0.0;
	}
	else {
#if doubleprecision == 1 
			norma = sqrt(s);
		
#else 
			norma = sqrtf(s);
#endif
	}
	printf("%e\n",norma);
	return norma;
} // NormaV



 // Скалярное произведение двух векторов
// 14.01.2018.
// Перенесено в функцию my_aggregat_amg.cu т.к. она первая использует данный функционал.
//doublereal Scal(doublereal *v1, doublereal *v2, integer n);

//----------метод сопряженных градиентов---------------
/* Входные параметры:
*  A - неразреженная матрица СЛАУ,
*  dV - вектор правой части, 
*  x - начальное приближение к решению или nullptr.
*  n - размерность СЛАУ An*n.
*  Матрица A полагается положительно определённой и 
*  симметричной (диагональное преобладание присутствует).
*  Количество итераций ограничено 1000, т.к. предполагается,
*  что если решение не сошлось за 1000 итераций то оно и не сойдётся.
*  Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
*/
doublereal *SoprGrad(doublereal **A, doublereal *dV, doublereal *x, integer n){
	printf("Reshenie metodom sopryjennyh gradientov:\n");
	integer k=0;
	integer i; // счётчик
	doublereal *ap=nullptr,
		 *z=nullptr, *p=nullptr;

	ap = new doublereal[n];
	z = new doublereal[n];
	p = new doublereal[n];

	doublereal a=0.0, b=0.0, nz=0.0;

	// шаг 1.1
	//X0==
	if (x==nullptr) {
        x=new doublereal[n];
		for(i=0;i<n;++i) x[i] = 0.0;
	}

	// пороговое значение невязки
	doublereal e = dterminatedTResudual;
	
	// шаг 1.2
    // вычисление z - невязки начального приближения
	if (ap != nullptr) {
		delete[] ap;
		ap = nullptr;
	}
	ap=MatrixByVector(A,x,n);
	for (i=0; i<n; ++i) z[i]=dV[i]-ap[i];

	if (Scal(z,z,n)!=0){
		// шаг 1.3
	   for (i=0; i<n; ++i)	p[i]=z[i];
	   nz=1000.;
	   while ((nz>e) && (k<1000)) {
		   // шаг 2.1
		  if (ap != nullptr) {
			   delete[] ap;
			   ap = nullptr;
		  }
	 	  ap=MatrixByVector(A,p,n);
		  if (ap != nullptr) {
			  // шаг 2.2
			  //a=Scal(z,p,n)/Scal(z,ap,n);
			  a = Scal(z, p, n) / Scal(ap, p, n); // шаговый множитель
			  // шаг 2.3 и 2.4
			  for (i = 0; i < n; ++i) {
				  x[i] += a * p[i]; // очередное приближение
				  z[i] -= a * ap[i]; // невязка k+1-го приближения
			  }
			  // шаг 2.5
			  nz = NormaV(z, n);
			  if (k % 10 == 0) printf("iter residual\n");
#if doubleintprecision == 1
			  printf(" %lld %e\n", k, nz);
#else
			  printf(" %d %e\n", k, nz);
#endif

			  // шаг 3.1
			  b = Scal(z, ap, n) / Scal(p, ap, n);
			  // шаг 3.2
			  for (i = 0; i < n; ++i) {
				  p[i] = z[i] - b * p[i]; // новое направление минимизации
			  }
		  }
		  else {
			  printf("ERROR!!! MatrixByVector(A,p,n); return nullptr pointer.\n");
			  system("PAUSE");
			  exit(1);
		  }
          // шаг 3.3 
		  k++;
	   } // while

	   // Освобождение памяти
	   if (ap != nullptr) {
		   delete[] ap;
	   }
	   if (z != nullptr) {
		   delete[] z;
	   }
	   if (p != nullptr) {
		   delete[] p;
	   }

	   return x;
	}
	else {
		// Освобождение памяти
		if (ap != nullptr) {
			delete[] ap;
		}
		if (z != nullptr) {
			delete[] z;
		}
		if (p != nullptr) {
			delete[] p;
		}

		return x;
	}
} // SoprGrad


  // умножение матрицы на вектор
  // используется формат хранения CRS
  // Разреженная матрица A (val, col_ind, row_ptr) квадратная размером n*n.
  // Число уравнений равно числу неизвестных и равно n.
  // Уданной функции три эквивалентных представления отличающихся лишь типом аргументов. 13.января.2018
// Реализация перенесена в файл my_agregat_amg.cu т.к. он первым использует эту функцию.
//void MatrixCRSByVector(doublereal* val, integer* col_ind, integer* row_ptr, doublereal* V, doublereal* &tmp, integer n);

// умножение матрицы на вектор (отладочный вариант для поиска ошибок).
// используется формат хранения CRS
// Разреженная матрица A (val, col_ind, row_ptr) квадратная размером n*n.
// Число уравнений равно числу неизвестных и равно n.
void MatrixCRSByVectordebug(doublereal* val, integer* col_ind, integer* row_ptr, doublereal* V, doublereal* &tmp, integer n)
{
	integer i,j; // Счётчики цикла

    // вектор tmp индексируется начиная с нуля так же как и вектор V
	for (i=0; i<n; ++i) tmp[i]=0.0;
	/*
	// В целях увеличения быстродействия 
	// вся необходимая память выделяется заранее.
	if (tmp == nullptr)
	{
		printf("malloc: out of memory for vector tmp in MatrixCRSByVector\n"); // нехватка памяти
		system("pause");
		exit(0);  // завершение программы
	}*/
	
	doublereal sum;
	integer rowend, rowbeg;
    
/*
#ifdef _OPENMP
	omp_set_num_threads(inumcore);
#endif
*/

    //#pragma omp parallel for shared(row_ptr, val, col_ind, V, tmp) private(sum, rowend, rowbeg, i, j) schedule (guided)
	for (i=0; i<n; ++i) {
	    sum = 0.0;
		if (i > 23897) {

#if doubleintprecision == 1
		printf("diagnostic message node %lld\n", i);
		printf("start=%lld, end=%lld\n", row_ptr[i], row_ptr[i + 1]);
		for (j = rowbeg; j < rowend; ++j)
		{
			printf("val=%e, V=%e, col_ind=%lld, j=%lld\n", val[j], V[col_ind[j]], col_ind[j], j);
		}
#else
			printf("diagnostic message node %d\n", i);
			printf("start=%d, end=%d\n", row_ptr[i], row_ptr[i + 1]);
			for (j = rowbeg; j < rowend; ++j)
			{
				printf("val=%e, V=%e, col_ind=%d, j=%d\n", val[j], V[col_ind[j]], col_ind[j], j);
			}
#endif
		 	    
			
#if doubleintprecision == 1
			printf("diagnostic message node %lld\n", i);
#else
			printf("diagnostic message node %d\n", i);
#endif
			
			//system("pause");
			system("pause");
		}
		rowend=row_ptr[i+1];
		rowbeg=row_ptr[i];
	    for (j = rowbeg; j<rowend; ++j)
		{
		    	sum += val[j]*V[col_ind[j]];
		}
		tmp[i] = sum;
	}
	
	//return tmp;
} // MatrixCRSByVectordebug

// умножение транспонированной матрицы на вектор
// (используется, например, в методе BiCG - бисопряжённых градиентов)
// для исходной (не транспонированной матрицы) используется формат хранения CRS
// Разреженная матрица A (val, col_ind, row_ptr) квадратная размером n*n.
// Число уравнений равно числу неизвестных и равно n.
doublereal* MatrixTransposeCRSByVector(doublereal* val, integer* col_ind, integer* row_ptr, doublereal* V, integer n)
{
	
	doublereal* tmp=new doublereal[n]; // вектор индексируется начиная с нуля так же как и вектор V
	if (tmp == nullptr)
	{
		printf("malloc: out of memory for vector tmp in MatrixTransposeCRSByVector\n"); // нехватка памяти
		//system("pause");
		system("pause");
		exit(0);
		return nullptr; // завершение программы
	}
	
	
    integer i,j; // Счётчики цикла
	integer rowend, rowbeg;
    
	for (i=0; i<n; ++i) tmp[i]=0.0;

	for (j=0; j<n; ++j) {
		rowend=row_ptr[j+1];
		rowbeg=row_ptr[j];
	    for (i = rowbeg; i<rowend; ++i)
		{
		    	tmp[col_ind[i]] += val[i]*V[j];
		}
	}
	
	return tmp;
} // MatrixTransposeCRSByVector


/* Метод сопряжённых градиентов Хестенса и Штифеля [1952]
*  Входные параметры:
*  val, col_ind, row_ptr - разреженная матрица СЛАУ в формате CRS,
*  dV - вектор правой части, 
*  x - начальное приближение к решению или nullptr.
*  n - размерность СЛАУ An*n.
*  Разреженная матрица A (val, col_ind, row_ptr) квадратная размером n*n.
*  Число уравнений равно числу неизвестных и равно n.
*  Матрица A полагается положительно определённой и 
*  симметричной (диагональное преобладание присутствует).
*  Количество итераций ограничено 1000, т.к. предполагается,
*  что если решение не сошлось за 1000 итераций то оно и не сойдётся.
*  Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
*/
doublereal *SoprGradCRS(doublereal *val, integer* col_ind, integer* row_ptr, doublereal *dV, doublereal *x, integer n){
	printf("Conjugate Gradients Method...:\n");
	integer k=0;
	
	doublereal *ap=nullptr,
		 *z=nullptr, *p=nullptr;

	// Выделение оперативной памяти.
	ap = new doublereal[n];
	z = new doublereal[n];
	p = new doublereal[n];

	doublereal a=0.0, b=0.0, nz=0.0;

//#ifdef _OPENMP
  //  omp_set_num_threads(inumcore);
//#endif

	// шаг 1.1
	//X0==
	if (x==nullptr) {
        x=new doublereal[n];
#pragma omp parallel for
		for(integer i=0;i<n;++i) x[i] = 0.0;
	}

	// пороговое значение невязки
	doublereal e = dterminatedTResudual;
	
	// шаг 1.2
    // вычисление z - невязки начального приближения
	MatrixCRSByVector(val,col_ind,row_ptr,x,ap,n);
	
    #pragma omp parallel for shared(z,dV,ap) schedule (guided)
	for (integer i=0; i<n; ++i) z[i]=dV[i]-ap[i];

	if (Scal(z,z,n)!=0){
		// шаг 1.3
       #pragma omp parallel for shared(p,z) schedule (guided)
	   for (integer i=0; i<n; ++i)	p[i]=z[i];

	   nz=1000.;
	   while ((nz>e) && (k<2*n)) {
		   // шаг 2.1
		  // чтобы избежать утечки памяти
	 	  MatrixCRSByVector(val,col_ind,row_ptr,p,ap,n);
		  // шаг 2.2
		  //a=Scal(z,p,n)/Scal(z,ap,n);
		  a=Scal(z,p,n)/Scal(ap,p,n); // шаговый множитель
		  // шаг 2.3 и 2.4
		  #pragma omp parallel for shared(x,z,p,ap,a) schedule (guided)
		  for (integer i=0; i<n; ++i) {
		      x[i]+=a*p[i]; // очередное приближение
			  z[i]-=a*ap[i]; // невязка k+1-го приближения
		  }
		  // шаг 2.5
		  nz=NormaV(z,n);
		  if (k%10==0) printf("iter residual\n");
#if doubleintprecision == 1
		  printf(" %lld %e\n", k, nz);
#else
		  printf(" %d %e\n", k, nz);
#endif
		  
		  // шаг 3.1
		  b=Scal(z,ap,n)/Scal(p,ap,n);
		  // шаг 3.2
		  #pragma omp parallel for shared(p,z,b) schedule (guided)
		  for (integer i=0; i<n; ++i) {
		     p[i]=z[i]-b*p[i]; // новое направление минимизации
		  }
          // шаг 3.3 
		  k++;
	   } // while

	   // Освобождение памяти
	   if (ap != nullptr) {
		   delete[] ap;
	   }
	   if (z != nullptr) {
		   delete[] z;
	   }
	   if (p != nullptr) {
		   delete[] p;
	   }

	   return x;
	}
	else {
		// Освобождение памяти
		if (ap != nullptr) {
			delete[] ap;
		}
		if (z != nullptr) {
			delete[] z;
		}
		if (p != nullptr) {
			delete[] p;
		}

		return x;
	}
} // SoprGradCRS

// Метод бисопряжённых градиентов
// для возможно несимметричной матрицы А (val, col_ind, row_ptr).
// Запрограммировано по книжке Баландин, Шурина: "Методы
// решения СЛАУ большой размерности".
// dV - правая часть СЛАУ,
// x - начальное приближение к решению или nullptr.
// n - размерность А n*n.
// Количество итераций ограничено 2000.
// Точность выхода по невязке задаётся в глобальной константе:
//  dterminatedTResudual.
void BiSoprGradCRS(doublereal *val, integer* col_ind, integer* row_ptr, doublereal *dV, doublereal* &x, integer n, integer maxit){
	printf("BiConjugate Gradients Method...:\n");

	doublereal *r=new doublereal[n], *r_tilda=new doublereal[n];
	doublereal *p=new doublereal[n], *p_tilda=new doublereal[n];
	doublereal nz=0.0; // невязка
	doublereal *ap=new doublereal[n];
	doublereal a=0.0,b=0.0, dold=0.0, dnew=0.0;

	integer i=0; // счётчик цикла for
	integer k=0; // номер итерации.

	// Начальное приближение:
    //X0==
	if (x==nullptr) {
        x=new doublereal[n];
		for(i=0;i<n;++i) x[i] = 0.0;
	}

	// пороговое значение невязки
	doublereal e = 1e-10;// dterminatedTResudual;

	MatrixCRSByVector(val,col_ind,row_ptr,x,ap,n);
	for (i=0; i<n; ++i) {
		r[i]=dV[i]-ap[i];
		r_tilda[i]=r[i];
		p[i]=r[i];
		p_tilda[i]=r_tilda[i];
	}

	nz=NormaV(r,n); // начальное значение невязки
	dold=Scal(r,r_tilda,n);

    while ((nz>e) && (k<maxit)) {
		MatrixCRSByVector(val,col_ind,row_ptr,p,ap,n);

		a=dold/Scal(ap,p_tilda,n);
		for (i=0; i<n; ++i) {
           x[i]+=a*p[i];
		   r[i]-=a*ap[i];
		}
		if (ap != nullptr) {
			delete[] ap;
			ap = nullptr;
		}
		ap=MatrixTransposeCRSByVector(val,col_ind,row_ptr,p_tilda,n);
        for (i=0; i<n; ++i) {
			r_tilda[i]-=a*ap[i];
		}
		dnew=Scal(r,r_tilda,n);
		b=dnew/dold;
		dold=dnew;
		// вычисление невязки.
        nz=NormaV(r,n);
		if (k%10==0) printf("iter residual\n");
#if doubleintprecision == 1
		printf(" %lld %e\n", k, nz);
#else
		printf(" %d %e\n", k, nz);
#endif
		

		if (fabs(b) < 1e-270) {
			printf("\nBiCG divergence detected...\n");
            //system("pause");
			system("pause");
			exit(0); // выход из приложения.
			break; // выход из цикла while
		}

        for (i=0; i<n; ++i) {
			p[i]=r[i]+b*p[i];
			p_tilda[i]=r_tilda[i]+b*p_tilda[i];
		}

		k++; // переход к следующей итерации.
	}

	// Освобождение памяти
	if (r != nullptr) {
		delete[] r;
		r = nullptr;
	}
	if (r_tilda != nullptr) {
		delete[] r_tilda;
		r_tilda = nullptr;
	}
	if (p != nullptr) {
		delete[] p;
		p = nullptr;
	}
	if (p_tilda != nullptr) {
		delete[] p_tilda;
		p_tilda = nullptr;
	}
	if (ap != nullptr) {
		delete[] ap;
		ap = nullptr;
	}

	//system("pause");
	system("pause");

	//return x;

} // BiSoprGradCRS

// Прямой ход по разреженной нижнетреугольной матрице L.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// L - хранится в следующем виде:
// 1. ldiag - диагональные элементы L.
// 2. lltr - поддиагональные элементы в строчном формате,
// т.е. хранение построчное. 
// 3. jptr - соответствующие номера столбцов для lltr, 
// 4. iptr - информация о начале следующей строки для lltr.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(L)*f;
// Вектор f портится.
// пример (CSIR - формат):
//  L = 
//  9.0   0.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   11.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   2.0   10.0   0.0   0.0   0.0   0.0   
//  3.0   1.0   2.0   9.0   0.0   0.0   0.0   
//  1.0   0.0   0.0   1.0   12.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  1.0   2.0   0.0   0.0   1.0   0.0   8.0   
// ------------------------------------------
// ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
doublereal* inverseL(doublereal* f, doublereal* ldiag, doublereal* lltr, integer* jptr, integer* iptr, integer n) {
	doublereal *z=new doublereal[n];
	for (integer ii = 0; ii < n; ++ii) z[ii] = 0.0; // initialization

    if (z == nullptr)
	{
		printf("malloc: out of memory for vector z in inverse(L)*f \n"); // нехватка памяти
		//system("pause");
		system("pause");
		exit(0);
		return nullptr; // завершение программы
	}
	else {

		integer i = 0, j = 0;
		for (i = 0; i < n; ++i) {
			for (j = iptr[i]; j < iptr[i + 1]; ++j) {
				f[i] -= z[jptr[j]] * lltr[j];
			}
			z[i] = f[i] / ldiag[i];
		}
		return z;
	}

	return z;

}//inverseL

// Прямой ход по разреженной нижнетреугольной матрице L.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// L - хранится в следующем виде:
// 1. val - диагональные и поддиагональные элементы L.
// в столбцовом порядке. 
// 3. indx - соответствующие номера строк для val, 
// 4. pntr - информация о начале следующего столбца.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(L)*f;
// Вектор f портится.
// пример (CSIR - формат):
//  L = 
//  9.0   0.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   11.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   2.0   10.0   0.0   0.0   0.0   0.0   
//  3.0   1.0   2.0   9.0   0.0   0.0   0.0   
//  1.0   0.0   0.0   1.0   12.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  1.0   2.0   0.0   0.0   1.0   0.0   8.0   
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseL_ITL(doublereal* f, doublereal* val, integer* indx, integer* pntr, doublereal* &z, integer n, int* &inumerate) {
	
	// doublereal **fbuf;
	// набор векторов fbuf нужен только в параллельной версии, в серийной версии можно просто передавать nullptr.
	// количество векторов в fbuf равно количеству потоков.

    if (z == nullptr)
	{
		// Попробуем выделить память. 23.03.2019
		z=new doublereal[n];
		if (z==nullptr) {
			printf("malloc: out of memory for vector z in inverse(L)*f \n"); // нехватка памяти
		   // system("pause");
			system("pause");
		    exit(0); // завершение программы
		}
	}

	
	if (bStableVersion) {

		for (integer i = 0; i < n; ++i) {
			z[i] = f[i] / val[pntr[i]];

			// обработка i-го столбца
			// эта часть не поддаётся распараллеливанию.
			// из за зависимостей по данным для f.
			for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
				f[indx[j]] -= z[i] * val[j];
			}

		}

	}
	else {

		//bool bserial=true;

		//if (bserial) {


		if (0) {

			for (integer i = 0; i < n; ++i) {
				z[i] = f[i] / val[pntr[i]];

				// обработка i-го столбца
				// эта часть не поддаётся распараллеливанию.
				// из за зависимостей по данным для f.
				for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
					f[indx[j]] -= z[i] * val[j];
				}

			}

		}
		else {

			// однопоточное исполнение.
			if (number_cores() == 1) {

				for (integer i = 0; i < n; ++i) {
					z[i] = f[i] / val[pntr[i]];

					// обработка i-го столбца
					// эта часть не поддаётся распараллеливанию.
					// из за зависимостей по данным для f.
					for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
						f[indx[j]] -= z[i] * val[j];
					}

				}

			}
			if (number_cores() == 2) {

#ifdef _OPENMP 
				omp_set_num_threads(number_cores()); // установка числа потоков
#endif


				


#pragma omp parallel sections
				{
#pragma omp section
					{
						for (integer i = s_par[1].s; i < s_par[1].e; ++i) {

							//if (inumerate[i] == 1)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (integer i = s_par[2].s; i < s_par[2].e; ++i) {

							//if (inumerate[i] == 2)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}


				for (integer i = s_par[3].s; i < s_par[3].e; ++i) {

					//if (inumerate[i] == 3)
					{

						z[i] = f[i] / val[pntr[i]];

						// обработка i-го столбца
						// эта часть не поддаётся распараллеливанию.
						// из за зависимостей по данным для f.
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[indx[j]] -= z[i] * val[j];
						}
					}

				}

			}

			if (number_cores() == 3)
			{
#pragma omp parallel sections
				{
#pragma omp section
					{
						for (integer i = s_par[1].s; i < s_par[1].e; ++i) {

							//if (inumerate[i] == 1)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (integer i = s_par[2].s; i < s_par[2].e; ++i) {

							//if (inumerate[i] == 2)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (integer i = s_par[3].s; i < s_par[3].e; ++i) {

							//if (inumerate[i] == 3)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}

#pragma omp parallel sections
				{
#pragma omp section
					{
						for (integer i = s_par[4].s; i < s_par[4].e; ++i) {

							//if (inumerate[i] == 4)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (integer i = s_par[5].s; i < s_par[5].e; ++i) {

							//if (inumerate[i] == 5)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}

			}

			if ((number_cores() == 4)|| (number_cores() == 5)) {

#ifdef _OPENMP 
				omp_set_num_threads(number_cores()); // установка числа потоков
#endif





#pragma omp parallel sections
				{
#pragma omp section
					{
						for (integer i = s_par[4].s; i < s_par[4].e; ++i) {

							//if (inumerate[i] == 4)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (integer i = s_par[6].s; i < s_par[6].e; ++i) {

							//if (inumerate[i] == 6)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (integer i = s_par[7].s; i < s_par[7].e; ++i) {

							//if (inumerate[i] == 7)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (integer i = s_par[9].s; i < s_par[9].e; ++i) {

							//if (inumerate[i] == 9)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}


#pragma omp parallel sections
				{
#pragma omp section
					{
						for (integer i = s_par[5].s; i < s_par[5].e; ++i) {

							//if (inumerate[i] == 5)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
#pragma omp section
					{
						for (integer i = s_par[8].s; i < s_par[8].e; ++i) {

							//if (inumerate[i] == 8)
							{

								z[i] = f[i] / val[pntr[i]];

								// обработка i-го столбца
								// эта часть не поддаётся распараллеливанию.
								// из за зависимостей по данным для f.
								for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
									f[indx[j]] -= z[i] * val[j];
								}
							}

						}
					}
				}


				for (integer i = s_par[3].s; i < s_par[3].e; ++i) {

					//if (inumerate[i] == 3)
					{

						z[i] = f[i] / val[pntr[i]];

						// обработка i-го столбца
						// эта часть не поддаётся распараллеливанию.
						// из за зависимостей по данным для f.
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[indx[j]] -= z[i] * val[j];
						}
					}

				}

			}


		}
	}
		/*

	}
	else {
		// параллельное исполнение.
		// параллельный код требует правильного разрешения зависимостей по данным.

		// Нам понадобиться 
		// n=omp_get_max_threads(); 
		// дополнительных векторов.

		integer nt=0;
#pragma omp parallel shared(nt)
		{
			// число нитей.
			nt=omp_get_max_threads();
		}

		

		for (integer i=0; i<nt; ++i) {
			for (integer j=0; j<n; ++j) {
				fbuf[i][j]=0.0; // инициализация.
			}
		}

#pragma omp for  shared(n, z, val, f, fbuf, pntr, indx, fbuf) 
		for (integer i=0; i<n; ++i) {
		   // Проблема в том что здесь используется f[i], а оно может быть обновлённым, что здесь не учитывается !!!
            z[i]=f[i]/val[pntr[i]];
		    // обработка i-го столбца
		    // эта часть не поддаётся распараллеливанию.
            // из за зависимостей по данным для f.
		    for (integer j=pntr[i]+1; j<pntr[i+1]; ++j) {
			    fbuf[omp_get_thread_num()][indx[j]]-=z[i]*val[j];
		    }
		
	    }

	}
	*/
}//inverseL_ITL

// Прямой ход по разреженной нижнетреугольной матрице L.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// L - хранится в следующем виде:
// 1. val - диагональные и поддиагональные элементы L.
// в столбцовом порядке. 
// 3. indx - соответствующие номера строк для val, 
// 4. pntr - информация о начале следующего столбца.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(L)*f;
// Вектор f портится.
// пример (CSIR - формат):
//  L = 
//  9.0   0.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   11.0   0.0   0.0   0.0   0.0   0.0   
//  0.0   2.0   10.0   0.0   0.0   0.0   0.0   
//  3.0   1.0   2.0   9.0   0.0   0.0   0.0   
//  1.0   0.0   0.0   1.0   12.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  1.0   2.0   0.0   0.0   1.0   0.0   8.0   
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseL_ITL(doublereal* f, doublereal* val, integer* indx, integer* pntr, doublereal*& z, integer n) {

	// doublereal **fbuf;
	// набор векторов fbuf нужен только в параллельной версии, в серийной версии можно просто передавать nullptr.
	// количество векторов в fbuf равно количеству потоков.

	if (z == nullptr)
	{
		// Попробуем выделить память. 23.03.2019
		z = new doublereal[n];
		if (z == nullptr) {
			printf("malloc: out of memory for vector z in inverse(L)*f \n"); // нехватка памяти
		   // system("pause");
			system("pause");
			exit(0); // завершение программы
		}
	}


	

		for (integer i = 0; i < n; ++i) {
			z[i] = f[i] / val[pntr[i]];

			// обработка i-го столбца
			// эта часть не поддаётся распараллеливанию.
			// из за зависимостей по данным для f.
			for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
				f[indx[j]] -= z[i] * val[j];
			}

		}


		/*

	}
	else {
		// параллельное исполнение.
		// параллельный код требует правильного разрешения зависимостей по данным.

		// Нам понадобиться
		// n=omp_get_max_threads();
		// дополнительных векторов.

		integer nt=0;
#pragma omp parallel shared(nt)
		{
			// число нитей.
			nt=omp_get_max_threads();
		}



		for (integer i=0; i<nt; ++i) {
			for (integer j=0; j<n; ++j) {
				fbuf[i][j]=0.0; // инициализация.
			}
		}

#pragma omp for  shared(n, z, val, f, fbuf, pntr, indx, fbuf)
		for (integer i=0; i<n; ++i) {
		   // Проблема в том что здесь используется f[i], а оно может быть обновлённым, что здесь не учитывается !!!
			z[i]=f[i]/val[pntr[i]];
			// обработка i-го столбца
			// эта часть не поддаётся распараллеливанию.
			// из за зависимостей по данным для f.
			for (integer j=pntr[i]+1; j<pntr[i+1]; ++j) {
				fbuf[omp_get_thread_num()][indx[j]]-=z[i]*val[j];
			}

		}

	}
	*/

} //inverseL_ITL

// Обратный ход по разреженной верхнетреугольной матрице U.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// U=transpose(L);
// U - хранится в следующем виде:
// 1. udiag - диагональные элементы U.
// 2. uutr - наддиагональные элементы в столбцовом формате,
// т.е. хранение постолбцовое. 
// Так портрет симметричен, то:
// 3. jptr - соответствующие номера столбцов для lltr, 
// 4. iptr - информация о начале следующей строки для lltr.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(U)*f;
// Вектор f портится.
// пример (CSIR - формат):
//  U=transpose(L) = 
//  9.0   0.0   0.0   3.0   1.0   0.0   1.0   
//  0.0   11.0   2.0   1.0   0.0   0.0   2.0   
//  0.0   0.0   10.0   2.0   0.0   0.0   0.0   
//  0.0   0.0   0.0   9.0   1.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   12.0   0.0   1.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   0.0   8.0 
// ------------------------------------------
// udiag==ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// uutr==lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
doublereal* inverseU(doublereal* f, doublereal* udiag, doublereal* uutr, integer* jptr, integer* iptr, integer n) {
	doublereal *z=new doublereal[n];

    if (z == nullptr)
	{
		printf("malloc: out of memory for vector z in inverse(U)*f \n"); // нехватка памяти
		//system("pause");
		system("pause");
		exit(0);
		return nullptr; // завершение программы
	}

	integer i,j;
	for (i=(n-1); i>=0; i--) {
        z[i]=f[i]/udiag[i];
		// Обработка i-го столбца над диагональю:
		for (j=iptr[i]; j<iptr[i+1]; ++j) {
			f[jptr[j]]-=z[i]*uutr[j];
		}
		
	}
	return z;
}//inverseU



// Обратный ход по разреженной верхнетреугольной матрице U.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// U=transpose(L); - верхняя треугольная матрица.
// U - хранится в следующем виде:
// 1. val - диагональные и наддиагональные элементы U (в строковом формате).
// 2. indx - соответствующие номера столбцов, 
// 3. pntr - информация о начале следующей строки для val.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(U)*f;
// Вектор f портится.
// пример (CSIR_ITL - формат):
//  U=transpose(L) = 
//  9.0   0.0   0.0   3.0   1.0   0.0   1.0   
//  0.0   11.0   2.0   1.0   0.0   0.0   2.0   
//  0.0   0.0   10.0   2.0   0.0   0.0   0.0   
//  0.0   0.0   0.0   9.0   1.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   12.0   0.0   1.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   0.0   8.0 
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseU_ITL(doublereal* f, doublereal* val, integer* indx, integer* pntr, doublereal*& z, integer n, int* &inumerate) {

	if (z == nullptr)
	{
		z = new doublereal[n];
		if (z == nullptr) {
			printf("malloc: out of memory for vector z in inverse(U)*f \n"); // нехватка памяти
			//system("pause");
			system("pause");
			exit(0); // завершение программы
		}
	}




	// 09.08.2021

	if (number_cores() == 1) {

		for (integer i = (n - 1); i >= 0; --i) {

			// Обработка i-ой строки:
			// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
			//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
			for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
				f[i] -= z[indx[j]] * val[j];
			}
			// делим на диагональный элемент:
			z[i] = f[i] / val[pntr[i]];

		}

	}

	if (number_cores() == 2) {

		for (integer i = (s_par[3].e - 1); i >= s_par[3].s; --i) {

			//if (inumerate[i] == 3)
			{

				// Обработка i-ой строки:
				// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
				//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
				for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
					f[i] -= z[indx[j]] * val[j];
				}
				// делим на диагональный элемент:
				z[i] = f[i] / val[pntr[i]];

			}
		}

#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (integer i = (s_par[1].e - 1); i >= s_par[1].s; --i) {

					//if (inumerate[i] == 1) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (integer i = (s_par[2].e - 1); i >= s_par[2].s; --i) {

					//if (inumerate[i] == 2) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}

	}


	if (number_cores() == 3) {

#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (integer i = (s_par[4].e - 1); i >= s_par[4].s; --i) {

					//if (inumerate[i] == 4) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (integer i = (s_par[5].e - 1); i >= s_par[5].s; --i) {

					//if (inumerate[i] == 5) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}

		

#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (integer i = (s_par[1].e - 1); i >= s_par[1].s; --i) {

					//if (inumerate[i] == 1) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (integer i = (s_par[2].e - 1); i >= s_par[2].s; --i) {

					//if (inumerate[i] == 2) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}

#pragma omp section 
			{
				for (integer i = (s_par[3].e - 1); i >= s_par[3].s; --i) {

					//if (inumerate[i] == 3)
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}

	}


	if ((number_cores() == 4)||(number_cores() == 5)) {

		for (integer i = (s_par[3].e - 1); i >= s_par[3].s; --i) {

			//if (inumerate[i] == 3) 
			{

				// Обработка i-ой строки:
				// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
				//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
				for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
					f[i] -= z[indx[j]] * val[j];
				}
				// делим на диагональный элемент:
				z[i] = f[i] / val[pntr[i]];

			}
		}

#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (integer i = (s_par[5].e - 1); i >= s_par[5].s; --i) {

					//if (inumerate[i] == 5) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (integer i = (s_par[8].e - 1); i >= s_par[8].s; --i) {

					//if (inumerate[i] == 8)
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}




#pragma omp parallel sections
		{
#pragma omp section 
			{
				for (integer i = (s_par[4].e - 1); i >= s_par[4].s; --i) {

					//if (inumerate[i] == 4)
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (integer i = (s_par[6].e - 1); i >= s_par[6].s; --i) {

					//if (inumerate[i] == 6) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}

#pragma omp section 
			{
				for (integer i = (s_par[7].e - 1); i >= s_par[7].s; --i) {

					//if (inumerate[i] == 7) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
#pragma omp section 
			{
				for (integer i = (s_par[9].e - 1); i >= s_par[9].s; --i) {

					//if (inumerate[i] == 9) 
					{

						// Обработка i-ой строки:
						// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
						//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
						for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
							f[i] -= z[indx[j]] * val[j];
						}
						// делим на диагональный элемент:
						z[i] = f[i] / val[pntr[i]];

					}
				}
			}
		}

	}


}//inverseU_ITL



// Обратный ход по разреженной верхнетреугольной матрице U.
// симметричная положительно определённая матрица
// СЛАУ A представлена неполным разложением Холецкого 
// A~=L*transpose(L); L - нижняя треугольная матрица.
// U=transpose(L); - верхняя треугольная матрица.
// U - хранится в следующем виде:
// 1. val - диагональные и наддиагональные элементы U (в строковом формате).
// 2. indx - соответствующие номера столбцов, 
// 3. pntr - информация о начале следующей строки для val.
// f - вектор правой части размером nodes.
// возвращает вектор z=inverse(U)*f;
// Вектор f портится.
// пример (CSIR_ITL - формат):
//  U=transpose(L) = 
//  9.0   0.0   0.0   3.0   1.0   0.0   1.0   
//  0.0   11.0   2.0   1.0   0.0   0.0   2.0   
//  0.0   0.0   10.0   2.0   0.0   0.0   0.0   
//  0.0   0.0   0.0   9.0   1.0   0.0   0.0   
//  0.0   0.0   0.0   0.0   12.0   0.0   1.0   
//  0.0   0.0   0.0   0.0   0.0   8.0   0.0   
//  0.0   0.0   0.0   0.0   0.0   0.0   8.0 
// ------------------------------------------
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6
// pntr: 0 4 8 10 12 14 15 16
//-------------------------------------------
void inverseU_ITL(doublereal* f, doublereal* val, integer* indx, integer* pntr, doublereal* &z, integer n) {

    if (z == nullptr)
	{
		z = new doublereal[n];
		if (z==nullptr) {
			printf("malloc: out of memory for vector z in inverse(U)*f \n"); // нехватка памяти
		    //system("pause");
			system("pause");
		    exit(0); // завершение программы
		}
		for (integer i = 0; i < n; ++i) z[i] = 0.0; // initialization
	}
	

		

		// 09.08.2021


		for (integer i = (n - 1); i >= 0; --i) {

			// Обработка i-ой строки:
			// эта часть не поддаётся распараллеливанию зависимость по данным по z[].
			//#pragma omp parallel for shared(f, indx, z, val, i, pntr) private(j).
			for (integer j = pntr[i] + 1; j < pntr[i + 1]; ++j) {
				f[i] -= z[indx[j]] * val[j];
			}
			// делим на диагональный элемент:
			z[i] = f[i] / val[pntr[i]];

		}
	
	
}//inverseU_ITL



// Переводит из формата CSIR в формат CSIR_ITL
// форматы:
// CSIR: ldiag, lltr, jptr, iptr
// CSIR_ITL: val, indx, pntr
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
// формат CSIR:
// ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
//Формируем разреженный формат CSIR_ITL
//val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
//indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
//pntr: 0 4 8 10 12 14 15 16 
//--------------------------------------------
void convertCSIRtoCSIR_ITL(doublereal *ldiag, doublereal *lltr, integer *jptr, integer *iptr, integer n, integer nz, doublereal* &val, integer* &indx, integer* &pntr, integer nnz) {
	integer i,j,k;
	//nnz=n+nz; // размер массивов val и indx
	// выделение оперативной памяти:
	val = new doublereal[nnz];
	indx = new integer[nnz];
	pntr = new integer[n+1];
	for (i=0; i<=n; ++i) pntr[i]=nnz;

	if ((val == nullptr) || (indx == nullptr) || (pntr == nullptr))
	{
		printf("malloc: out of memory in convertCSIRtoCSIR_ITL \n"); // нехватка памяти
		//system("pause");
		system("pause");
		exit(0); // завершение программы
	}

	// Алгоритм:
	// По порядку для всех столбцов формата CSIR_ITL
	integer ic=0; // счётчик ненулевых элементов
	for (k=0; k<n; ++k) {
		// добавление диагонального элемента k - го стобца
		val[ic]=ldiag[k];
		indx[ic]=k;
		pntr[k]= myi_min(ic,pntr[k]);
		ic++;

		// добавление остальных элементов k-го столбца
		// сканирование матрицы в CSIR формате:
		for (i=1; i<n; ++i) {
			for (j=iptr[i]; j<iptr[i+1]; ++j)
				if (jptr[j] == k) {
					// добавление элемента в k-ый столбец
					val[ic]=lltr[j];
					indx[ic]=i;
                    pntr[k]= myi_min(ic,pntr[k]);
					ic++;
				}
		}

	}

} // convertCSIRtoCSIR_ITL

// Неполное разложение Холецкого
// для положительно определённой симметричной
// матрицы А размером n*n.
// n - размерность матрицы СЛАУ
// Матрица val изменяется и в ней возвращается
// неполное разложение Холецкого IC(0):
// val == U верхняя треугольная матрица
// A = transpose(U)*U=L*transpose(L);
// L=transpose(U);
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
//формат CSIR_ITL (верхний треугольник хранится построчно).
// val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
// indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
// pntr: 0 4 8 10 12 14 15 16 
//--------------------------------------------
// Результат факторизации без заполнения:
// изменённый массив val (indx и pntr остались без изменений):
// val (factorization)= 
// 3.0
// 1.0
// 0.3333333333333333
// 0.3333333333333333
// 3.3166247903554
// 0.6030226891555273
// 0.30151134457776363
// 0.6030226891555273
// 3.1622776601683795
// 0.6324555320336759
// 2.932575659723036
// 0.34099716973523675
// 3.4472773213410837
// 0.2578524458667825
// 2.8284271247461903
// 2.7310738989293286
//-------------------------------------------
void IC0Factor_ITL(doublereal* val, integer* indx, integer* pntr, integer n)
{
  integer d, g, h, i, j, k;
  doublereal z;

  for (k = 0; k < n - 1; ++k) {
    d = pntr[k];
    z = val[d] = sqrt(val[d]);

    for (i = d + 1; i < pntr[k+1]; ++i)
      val[i] /= z;

    for (i = d + 1; i < pntr[k+1]; ++i) {
      z = val[i];
      h = indx[i];
      g = i;

      for (j = pntr[h] ; j < pntr[h+1]; ++j)
        for ( ; g < pntr[k+1] && indx[g+1] <= indx[j]; g++)
          if (indx[g] == indx[j])
             val[j] -= z * val[g];
    }
  }
  d = pntr[n-1];
  val[d] = sqrt(val[d]);
} // IC0Factor_ITL

// Модифицированное неполное разложение Холецкого.
void IC0FactorModify_ITL(doublereal* val, integer* indx, integer* pntr, integer n)
{
  integer d, g, h, i, j, k;
  doublereal z, accumulate_fill_in;

  for (k = 0; k < n - 1; ++k) {
    d = pntr[k];
    z = val[d] = sqrt(val[d]);

    for (i = d + 1; i < pntr[k+1]; ++i)
      val[i] /= z;

    for (i = d + 1; i < pntr[k+1]; ++i) {
      z = val[i];
      h = indx[i];
      g = i;

      accumulate_fill_in = 0.0;

      for (j = pntr[h] ; j < pntr[h+1]; ++j)
        for ( ; g < pntr[k+1] && indx[g+1] <= indx[j]; g++)
          if (indx[g] == indx[j]) // номера столбцов равны
             val[j] -= z * val[g];
	  else //index does not match accumulate the fill-in value
		  accumulate_fill_in += z * val[g];

	  val[pntr[h]] -= accumulate_fill_in;

    }
  }
  d = pntr[n-1];
  val[d] = sqrt(val[d]);
} // IC0FactorModify_ITL

// Переводит из формата CSIR_ITL в формат CSIR (обратное преобразование)
// Память под все массивы предполагается выделенной заранее!!!
// форматы:
// CSIR_ITL: val, indx, pntr
// CSIR: ldiag, lltr, jptr, iptr
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
//Формируем разреженный формат CSIR_ITL
//val: 9.0 3.0 1.0 1.0 11.0 2.0 1.0 2.0 10.0 2.0 9.0 1.0 12.0 1.0 8.0 8.0 
//indx: 0 3 4 6 1 2 3 6 2 3 3 4 4 6 5 6 
//pntr: 0 4 8 10 12 14 15 16 
//--------------------------------------------
// формат CSIR:
// ldiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// lltr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
void convertCSIR_ITLtoCSIR(doublereal* ldiag, doublereal* lltr, integer* jptr, integer* iptr, integer n, integer nz, doublereal* val, integer* indx, integer* pntr, integer nnz) {
	integer i,j,k;//,k1;
	//integer imin=1;
	//nz=nnz-n; // размер массивов lltr и jptr
	// память предполагается выделенной заранее!!!
	// jptr и iptr изменяться не будут
	for (i=0; i<n; ++i) ldiag[i]=0.0;
	for (i=0; i<nz; ++i) {
		lltr[i]=0.0;
		//jptr[i]=0;
	}
	//for (i=0; i<=n; ++i) iptr[i]=nz;


	// Алгоритм:
	// По порядку для всех строк формата CSIR
	integer ic=0; // счётчик ненулевых элементов
	for (k=0; k<n; ++k) {
		// добавление диагонального элемента k - ой строки
		ldiag[k]=val[pntr[k]];

		// добавление остальных элементов k-ой строки
		// сканирование матрицы в CSIR_ITL формате:
		for (i=0; i<n-1; ++i) {
			for (j=pntr[i]+1; j<pntr[i+1]; ++j)
				if (indx[j] == k) {
					// добавление элемента в k-ую строку
					lltr[ic]=val[j];
					//jptr[ic]=i;
					//imin=myi_min(ic,iptr[k]);
                    //iptr[k]=imin;
					//if (imin==0) {
					//	for (k1=0; k1<k; k1++) iptr[k1]=0;
					//}
					ic++;
				}
		}

	}

} // convertCSIR_ITLtoCSIR

// неполное разложение Холецкого IC(0).
// входные данные нижний треугольник симметричной матрицы в формате CSIR.
// Внутри программы идут преобразования к формату CSIR_ITL библиотеки шаблонов ITL.
void ICFactor0(doublereal* ldiag, doublereal* lltr, integer* jptr, integer* iptr, integer n, integer nz) {
    
	doublereal *val;
	integer *indx, *pntr;

	// внутри происходит выделение памяти
	// преобразование (прямое и обратное) ресурсоёмкая операция для больших матриц,
	// поэтому от неё нужно отказаться.
	convertCSIRtoCSIR_ITL(ldiag, lltr, jptr, iptr, n, nz, val, indx, pntr, n+nz);
	printf("Incoplete Cholesky 49.9%%...\n");
	IC0Factor_ITL(val, indx, pntr, n);
	printf("Incoplete Cholesky 50%%...\n");
    convertCSIR_ITLtoCSIR(ldiag, lltr, jptr, iptr, n, nz, val, indx, pntr, n+nz);
	printf("Incoplete Cholesky 100%%...\n");

	// освобождение памяти
	delete val; delete indx; delete pntr;
} // ICFactor0


// умножение симметричной положительно определённой  матрицы на вектор 
// используется формат хранения CSIR. В силу симметрии хранятся только поддиагональные элементы altr. 
// Разреженная SPD матрица A (adiag, altr, jptr, iptr) квадратная размером n*n.
// Число уравнений равно числу неизвестных и равно n.
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   2.0   10.0   2.0   0.0   0.0   0.0    
// 3.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 1.0   2.0   0.0   0.0   1.0   0.0   8.0 
// ------------------------------------------
// формат CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0 1.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
void  SPDMatrixCSIRByVector(doublereal* adiag, doublereal* altr, integer* jptr, integer* iptr, doublereal* V, doublereal* &tmp, integer n)
{
	
	// вектор tmp индексируется начиная с нуля так же как и вектор V
	if (tmp == nullptr)
	{
		printf("in SPDMatrixCSIRByVector tmp==nullptr\n");
		//system("pause");
		system("pause");
		tmp =new doublereal[n];
		if (tmp==nullptr) {
			printf("malloc: out of memory for vector tmp in SPDMatrixCSIRByVector\n"); // нехватка памяти
		    //system("pause");
			system("pause");
		    exit(0); // завершение программы
		}
	}
	
	
   // integer i,j; // Счётчики цикла
    
/*
#ifdef _OPENMP
	omp_set_num_threads(inumcore);
#endif
*/

    #pragma omp parallel for shared(tmp, V, adiag)  
	for (integer i=0; i<n; ++i) tmp[i]=V[i]*adiag[i];

    // Последовательная секция
	/*
	for (i=0; i<n; ++i) {
	    for (j = iptr[i]; j<iptr[i+1]; ++j)
		{
		    tmp[i] += V[jptr[j]]*altr[j];
		    tmp[jptr[j]] += V[i]*altr[j];
		}
	}
	*/

	// Часть первая из двух.
	#pragma omp parallel for shared(tmp, V, altr, iptr, jptr,n) 
    for (integer i=0; i<n; ++i) {
	    for (integer j = iptr[i]; j<iptr[i+1]; ++j)
		{
		    tmp[i] += V[jptr[j]]*altr[j];
		}
	}

	if (number_cores() == 1) {

		// Вторая часть не поддаётся распараллеливанию
		for (integer i = 0; i < n; ++i) {

			// эта часть не поддаётся распараллеливанию.
			//#pragma omp parallel for shared(tmp, V, altr, i, iptr, jptr) private(j)
			for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
			{
				tmp[jptr[j]] += V[i] * altr[j];
			}
		}
	}

	if (number_cores() == 2) {


#pragma omp parallel sections
		{
#pragma omp section
			{
				// Вторая часть не поддаётся распараллеливанию
				for (integer i = s_par[1].s; i < s_par[1].e; ++i) {

					// эта часть не поддаётся распараллеливанию.
					//#pragma omp parallel for shared(tmp, V, altr, i, iptr, jptr) private(j)
					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				// Вторая часть не поддаётся распараллеливанию
				for (integer i = s_par[2].s; i < s_par[2].e; ++i) {

					// эта часть не поддаётся распараллеливанию.
					//#pragma omp parallel for shared(tmp, V, altr, i, iptr, jptr) private(j)
					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}

		// Вторая часть не поддаётся распараллеливанию
		for (integer i = s_par[3].s; i < s_par[3].e; ++i) {

			// эта часть не поддаётся распараллеливанию.
			//#pragma omp parallel for shared(tmp, V, altr, i, iptr, jptr) private(j)
			for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
			{
				tmp[jptr[j]] += V[i] * altr[j];
			}
		}
	}


	if (number_cores() == 3) {


#pragma omp parallel sections
		{
#pragma omp section
			{
				for (integer i = s_par[1].s; i < s_par[1].e; ++i) {

					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (integer i = s_par[2].s; i < s_par[2].e; ++i) {

					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (integer i = s_par[3].s; i < s_par[3].e; ++i) {

					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}


#pragma omp parallel sections
		{
#pragma omp section
			{
				for (integer i = s_par[4].s; i < s_par[4].e; ++i) {

					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (integer i = s_par[5].s; i < s_par[5].e; ++i) {

					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}

		
		
	}

	if ((number_cores() == 4)||(number_cores() == 5)) {

#pragma omp parallel sections
		{
#pragma omp section
			{
				for (integer i = s_par[4].s; i < s_par[4].e; ++i) {
					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (integer i = s_par[6].s; i < s_par[6].e; ++i) {
					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (integer i = s_par[7].s; i < s_par[7].e; ++i) {
					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (integer i = s_par[9].s; i < s_par[9].e; ++i) {
					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}

#pragma omp parallel sections
		{
#pragma omp section
			{
				for (integer i = s_par[5].s; i < s_par[5].e; ++i) {
					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
#pragma omp section
			{
				for (integer i = s_par[8].s; i < s_par[8].e; ++i) {
					for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
					{
						tmp[jptr[j]] += V[i] * altr[j];
					}
				}
			}
		}
		
		for (integer i = s_par[3].s; i < s_par[3].e; ++i) {
			for (integer j = iptr[i]; j < iptr[i + 1]; ++j)
			{
				tmp[jptr[j]] += V[i] * altr[j];
			}
		}
	}

} // SPDMatrixCSIRByVector

// умножение несимметричной положительно определённой  матрицы на вектор 
// используется формат хранения CSIR.  
// Разреженная матрица A (adiag, altr, autr, jptr, iptr) квадратная размером n*n.
// Число уравнений равно числу неизвестных и равно n.
// Диагональ adiag хранится отдельно. Нижний треугольник altr хранится построчно.
// Верхний треугольник хранится по столбцам autr. Портрет матрицы (позиции ненулевых 
// элементов ) предполагается симметричным. Массив jptr - номера столбцов для нижнего 
// треугольника, массив iptr - показывает где начинаются новые строки для нижнего треугольника.
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   1.0   10.0   2.0   0.0   0.0   0.0    
// 2.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 2.0   2.0   0.0   0.0   3.0   0.0   8.0 
// ------------------------------------------
// формат CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 1.0  2.0 1.0 2.0  1.0 1.0  2.0 2.0 3.0
// autr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
doublereal* MatrixCSIRByVector(doublereal* adiag, doublereal* altr, doublereal* autr, integer* jptr, integer* iptr, doublereal* V, integer n)
{
	
	doublereal* tmp=new doublereal[n]; // вектор индексируется начиная с нуля так же как и вектор V
	if (tmp == nullptr)
	{
		printf("malloc: out of memory for vector tmp in SPDMatrixCSIRByVector\n"); // нехватка памяти
		//system("pause");
		system("pause");
		exit(0);
		return nullptr; // завершение программы
	}
	
	
    integer i,j; // Счётчики цикла

	for (i=0; i<n; ++i) tmp[i]=V[i]*adiag[i];

    
	for (i=0; i<n; ++i) {
	    for (j = iptr[i]; j<iptr[i+1]; ++j)
		{
		    	tmp[i] += V[jptr[j]]*altr[j];
		        tmp[jptr[j]] += V[i]*autr[j];
		}
	}
	
	return tmp;
} // MatrixCSIRByVector

// умножение транспонированной несимметричной положительно определённой  матрицы на вектор 
// используется формат хранения CSIR.  
// Разреженная матрица A (adiag, altr, autr, jptr, iptr) квадратная размером n*n. Хранится 
// именно исходная матрица, а умножается её транспонированный вариант.
// Число уравнений равно числу неизвестных и равно n.
// Диагональ adiag хранится отдельно. Нижний треугольник altr хранится построчно.
// Верхний треугольник хранится по столбцам autr. Портрет матрицы (позиции ненулевых 
// элементов ) предполагается симметричным. Массив jptr - номера столбцов для нижнего 
// треугольника, массив iptr - показывает где начинаются новые строки для нижнего треугольника.
// пример:
// A = 
// 9.0   0.0   0.0   3.0   1.0   0.0   1.0    
// 0.0   11.0   2.0   1.0   0.0   0.0   2.0    
// 0.0   1.0   10.0   2.0   0.0   0.0   0.0    
// 2.0   1.0   2.0   9.0   1.0   0.0   0.0    
// 1.0   0.0   0.0   1.0   12.0   0.0   1.0    
// 0.0   0.0   0.0   0.0   0.0   8.0   0.0    
// 2.0   2.0   0.0   0.0   3.0   0.0   8.0 
// ------------------------------------------
// формат CSIR:
// adiag: 9.0 11.0 10.0 9.0 12.0 8.0 8.0
// altr: 1.0  2.0 1.0 2.0  1.0 1.0  2.0 2.0 3.0
// autr: 2.0 3.0 1.0 2.0 1.0 1.0 1.0 2.0
// jptr: 1 0 1 2 0 3 0 1 4
// iptr: 0 0 0 1 4 6 6 9
//-------------------------------------------
doublereal* MatrixTransposeCSIRByVector(doublereal* adiag, doublereal* altr, doublereal* autr, integer* jptr, integer* iptr, doublereal* V, integer n)
{
	
	doublereal* tmp=new doublereal[n]; // вектор индексируется начиная с нуля так же как и вектор V
	if (tmp == nullptr)
	{
		printf("malloc: out of memory for vector tmp in SPDMatrixCSIRByVector\n"); // нехватка памяти
		//system("pause");
		system("pause");
		exit(0);
		return nullptr; // завершение программы
	}
	
	
    integer i,j; // Счётчики цикла

	for (i=0; i<n; ++i) tmp[i]=V[i]*adiag[i];

    
	for (i=0; i<n; ++i) {
	    for (j = iptr[i]; j<iptr[i+1]; ++j)
		{
		    	tmp[i] += V[jptr[j]]*autr[j];
		        tmp[jptr[j]] += V[i]*altr[j];
		}
	}
	
	return tmp;
} // MatrixTransposeCSIRByVector


/* Метод сопряжённых градиентов Хестенса и Штифеля [1952]
*  Входные параметры:
*  adiag, altr, jptr, iptr - разреженная матрица СЛАУ в формате CSIR,
*  dV - вектор правой части, 
*  x - начальное приближение к решению или nullptr.
*  n - размерность СЛАУ An*n.
*  nz - размерность массивов altr, jptr.
*  Разреженная матрица A (adiag, altr, jptr, iptr) квадратная размером n*n.
*  Число уравнений равно числу неизвестных и равно n.
*  Матрица A полагается положительно определённой и 
*  симметричной (диагональное преобладание присутствует).
*  Хранится только нижний треугольник с диагональю altr и adiag.
*  Количество итераций ограничено 1000, т.к. предполагается,
*  что если решение не сошлось за 1000 итераций то оно и не сойдётся.
*  Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
*  В качестве предобуславливателя работает неполное разложение Холецкого:
*  M^(-1)==transpose(L)^(-1)*L^(-1); // обращённый предобуславливатель.
*  
*/
doublereal *SoprGradCSIR(doublereal* adiag, doublereal* altr, integer* jptr, integer* iptr, doublereal *dV, doublereal *x, integer n, integer nz){

	printf("Reshenie metodom sopryjennyh gradientov:\n");
	integer k=0;
	integer i=0; // счётчик
	doublereal *ap=nullptr, *vcopy=nullptr,
		 *z=nullptr, *p=nullptr;

	ap = new doublereal[n];
	vcopy = new doublereal[n];
	z = new doublereal[n];
	p = new doublereal[n];

    doublereal a=0.0, b=0.0, res=0.0;
	
	// для неполного разложения Холецкого:
	doublereal  *ldiag = nullptr, *lltr = nullptr;
	integer *jptrsort = nullptr;
	doublereal *f = nullptr;

	ldiag = new doublereal[n];
	lltr = new doublereal[nz];
	jptrsort = new integer[nz];
	f = new doublereal[n];

	doublereal dold=0.0, dnew=0.0;
	

	
	// инициализация
	for (i=0; i<n; ++i) ldiag[i]=adiag[i];
	for (i=0; i<nz; ++i) lltr[i]=altr[i];
	// неполное разложение Холецкого:
	// Возвращает левый нижний треугольный сомножитель.
	printf("Incoplete Cholesky decomposition beginig...:\n");
    ICFactor0(ldiag, lltr, jptr, iptr, n, nz);
	printf("Incoplete Cholesky decomposition finish...:\n");//*/

    
	for (i=0; i<nz; ++i) jptrsort[i]=jptr[i];
	for (i=0; i<n; ++i) QuickSort(jptrsort, iptr[i], iptr[i+1]-1);
    //printf("jptrsort...\n");
#if doubleintprecision == 1
	//for (i=0; i<nz; ++i) printf("%lld ",jptrsort[i]); system("pause");
#else
	//for (i=0; i<nz; ++i) printf("%d ",jptrsort[i]); system("pause");
#endif
	



	// шаг 1.1
	//X0==
	if (x==nullptr) {
        x=new doublereal[n];
		for(i=0;i<n;i++) x[i] = 0.0;
	}

	// пороговое значение невязки
	doublereal e = dterminatedTResudual;
	
	// шаг 1.2
    // вычисление z - невязки начального приближения
	SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, x, ap, n);
	for (i=0; i<n; ++i) z[i]=dV[i]-ap[i];
	for (i=0; i<n; ++i) vcopy[i]=z[i];
    f=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);
    for (i=0; i<n; ++i) vcopy[i]=f[i];
	if (f != nullptr) {
		delete[] f;
	}
	f=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
    dnew=Scal(z,f,n);

	if (fabs(dnew)>1.0e-100){
		// шаг 1.3
	   for (i=0; i<n; ++i)	p[i]=f[i];
	   res=1000.;
	   while ((fabs(res)>e) && (k<1000)) {
		   // шаг 2.1
		  SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, p, ap, n);

		  // шаг 2.2
		  a=dnew/Scal(p,ap,n);// шаговый множитель
		  // шаг 2.3 и 2.4
		  for (i=0; i<n; ++i) {
		      x[i]+=a*p[i]; // очередное приближение 
              z[i]-=a*ap[i];// невязка k+1-го приближения
		  }
          for (i=0; i<n; ++i) vcopy[i]=z[i];
		  if (f != nullptr) {
			  delete[] f;
		  }
          f=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);
		  if (f != nullptr) {
			  for (i = 0; i < n; ++i) vcopy[i] = f[i];
		  }
		  if (f != nullptr) {
			  delete[] f;
		  }
	      f=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
		  if (f != nullptr) {
			  // шаг 2.5
			  dold = dnew;
			  dnew = Scal(z, f, n);


			  res = dnew;
			  if (k % 10 == 0) printf("iter residual\n");
#if doubleintprecision == 1
			  printf(" %lld %e\n", k, res);
#else
			  printf(" %d %e\n", k, res);
#endif
			  
			  // шаг 3.1
			  b = dnew / dold;
			  // шаг 3.2
			  for (i = 0; i < n; ++i) {
				  p[i] = f[i] + b*p[i]; // новое направление минимизации
			  }
		  }
          // шаг 3.3
		  k++;
	   } // while

	   // Освобождение памяти
	   if (ap != nullptr) {
		   delete[] ap;
	   }
	   if (vcopy != nullptr) {
		   delete[] vcopy;
	   }
	   if (z != nullptr) {
		   delete[] z;
	   }
	   if (p != nullptr) {
		   delete[] p;
	   }
	   if (f != nullptr) {
		   delete[] f;
	   }

	   return x;
	}
	else {
		// Освобождение 
		if (ap != nullptr) {
			delete[] ap;
		}
		if (vcopy != nullptr) {
			delete[] vcopy;
		}
		if (z != nullptr) {
			delete[] z;
		}
		if (p != nullptr) {
			delete[] p;
		}
		if (f != nullptr) {
			delete[] f;
		}

		return x;
	}
} // SoprGradCSIR


// простая реализация явно преобразующая матрицу СЛАУ А.
// Матрица СЛАУ А задаётся в CSIR формате: adiag, altr, jptr, iptr.
// Неполное разложение Холецкого для А представляет её приближённо в виде:
// A = L*transpose(L); с нулевым заполнением. Массивы jptr и  iptr остаются теми же.
// Тогда матрица: A~=inverse(L)*A*inverse(transpose(L)) тоже симметрична и положительно определена.
// Правая часть преобразованной системы имеет вид: dV~=inverse(L)*dV.
// Решение СЛАУ тогда равно A~*x~=dV~; => x~=transpose(L)*x; => x=inverse(transpose(L))*x~;
// Предобуславливание неполным разложением Холецкого уменьшает количество итераций при решении СЛАУ,
// улучшает спектральные характеристики матрицы СЛАУ.
doublereal *SoprGradCSIR2(doublereal* adiag, doublereal* altr, integer* jptr, integer* iptr, doublereal *dV, doublereal *x, integer n, integer nz0){
	printf("Reshenie metodom sopryjennyh gradientov:\n");
	integer k=0;
	integer i=0; // счётчик
	doublereal *ap = nullptr, *vcopy = nullptr,
		 *z=nullptr, *p=nullptr;
	doublereal a=0.0, b=0.0, nz=0.0; // инициализация.

							   // для неполного разложения Холецкого:
	doublereal  *ldiag = nullptr, *lltr = nullptr;
	integer *jptrsort = nullptr;


	// allocate memory
	vcopy = new doublereal[n];
	z = new doublereal[n];
	p = new doublereal[n];

   

	ldiag = new doublereal[n];
	lltr = new doublereal[nz0];
	jptrsort = new integer[nz0];

	if ((vcopy != nullptr) && (z != nullptr) && (p != nullptr) && (ldiag != nullptr) && (lltr != nullptr) && (jptrsort != nullptr)) {

		// инициализация
		for (i = 0; i < n; ++i) ldiag[i] = adiag[i];
		for (i = 0; i < nz0; ++i) lltr[i] = altr[i];
		// неполное разложение Холецкого:
		// Возвращает левый нижний треугольный сомножитель.
		printf("Incoplete Cholesky decomposition beginig...:\n");
		ICFactor0(ldiag, lltr, jptr, iptr, n, nz0);
		printf("Incoplete Cholesky decomposition finish...:\n");//*/



	   /*
		ldiag[0]=1.0; ldiag[1]=1.0;  ldiag[2]=1.838477; ldiag[3]=2.00055;
		ldiag[4]=0.590477; ldiag[5]=1.0;  ldiag[6]=1.0;
		lltr[0]=-1.22383866; lltr[1]=-0.5439282932;  lltr[2]=-1.33247070; //*/

		/* // переставлены элементы
		ldiag[0]=1.0; ldiag[1]=1.0;  ldiag[2]=1.838477; ldiag[3]=2.00055;
		ldiag[4]=0.590477; ldiag[5]=1.465913;  ldiag[6]=0.37585673;
		lltr[0]=-1.22383866; lltr[1]=-1.33247070;  lltr[2]=-0.5439282932; lltr[3]=-0.1457305633;
		lltr[4]=-0.4998613742; lltr[5]=-1.401073265;  lltr[6]=-0.06498197865;//*/

		for (i = 0; i < nz0; ++i) jptrsort[i] = jptr[i];
		for (i = 0; i < n; ++i) QuickSort(jptrsort, iptr[i], iptr[i + 1] - 1);

		// шаг 1.1
		//X0==
		if (x == nullptr) {
			x = new doublereal[n];
			for (i = 0; i < n; ++i) x[i] = 0.0;
		}

		// пороговое значение невязки
		doublereal e = dterminatedTResudual;

		// шаг 1.2
		// вычисление z - невязки начального приближения
		//ap=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, x, n);
		//for (i=0; i<n; ++i) z[i]=dV[i]-ap[i];

		for (i = 0; i < n; ++i) vcopy[i] = x[i];
		ap = inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
		for (i = 0; i < n; ++i) vcopy[i] = ap[i];
		SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, ap, n);
		for (i = 0; i < n; ++i) vcopy[i] = ap[i];
		delete[] ap;
		ap = inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);
		for (i = 0; i < n; ++i) vcopy[i] = dV[i];
		delete[] dV;
		dV = inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);

		for (i = 0; i < n; ++i) z[i] = dV[i] - ap[i];

		if (Scal(z, z, n) != 0) {
			// шаг 1.3
			for (i = 0; i < n; ++i)	p[i] = z[i];
			nz = 1000.;
			while ((nz > e) && (k < 1000)) {
				// шаг 2.1
			   //ap=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, p, n);

				delete ap; // освобождение памяти
				for (i = 0; i < n; ++i) vcopy[i] = p[i];
				ap = inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
				for (i = 0; i < n; ++i) vcopy[i] = ap[i];
				SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, ap, n);
				for (i = 0; i < n; ++i) vcopy[i] = ap[i]; delete ap;
				ap = inverseL(vcopy, ldiag, lltr, jptrsort, iptr, n);

				// шаг 2.2
				//a=Scal(z,p,n)/Scal(z,ap,n);
				a = Scal(z, p, n) / Scal(ap, p, n); // шаговый множитель
				// шаг 2.3 и 2.4
				for (i = 0; i < n; ++i) {
					x[i] += a*p[i]; // очередное приближение
					z[i] -= a*ap[i]; // невязка k+1-го приближения
				}
				// шаг 2.5
				nz = NormaV(z, n);
				if (k % 10 == 0) printf("iter residual\n");
#if doubleintprecision == 1
				printf(" %lld %e\n", k, nz);
#else
				printf(" %d %e\n", k, nz);
#endif
				
				// шаг 3.1
				b = Scal(z, ap, n) / Scal(p, ap, n);
				// шаг 3.2
				for (i = 0; i < n; ++i) {
					p[i] = z[i] - b*p[i]; // новое направление минимизации
				}
				// шаг 3.3 
				k++;
			} // while

			// Освобождение 
			if (ap != nullptr) {
				delete[] ap;
			}
			if (vcopy != nullptr) {
				delete[] vcopy;
				vcopy = nullptr;
			}
			if (z != nullptr) {
				delete[] z;
			}
			if (p != nullptr) {
				delete[] p;
			}

			vcopy = new doublereal[n];
			if (vcopy != nullptr) {
				for (i = 0; i < n; ++i) vcopy[i] = x[i];
			}
			else {
				printf("vcopy is nullptr see SoprGradCSIR2 in mylinalg.c file\n");
				system("pause");
				exit(1);
			}
			if (x != nullptr) {
				delete[] x;
			}
			x = inverseU(vcopy, ldiag, lltr, jptrsort, iptr, n);
			if (vcopy != nullptr) {
				delete[] vcopy;
				vcopy = nullptr;
			}
			return x;
		}
		else {
			// Освобождение памяти
			if (ap != nullptr) {
				delete[] ap;
			}
			if (vcopy != nullptr) {
				delete[] vcopy;
			}
			if (z != nullptr) {
				delete[] z;
			}
			if (p != nullptr) {
				delete[] p;
			}

			return x;
		}

	}
	else {
		printf("problem memory allocate in SoprGradCSIR2\n");
		system("pause");
		return nullptr;
	}
} // SoprGradCSIR2

/* Метод сопряжённых градиентов Хестенса и Штифеля [1952]
*  Входные параметры:
*  M - разреженная матрица СЛАУ в формате SIMPLESPARSE,
*  dV - вектор правой части, 
*  x - начальное приближение к решению или nullptr.
*  n - размерность СЛАУ An*n.
*
*  Разреженная матрица M квадратная размером n*n.
*  Число уравнений равно числу неизвестных и равно n.
*  Матрица M предполагается положительно определённой и 
*  симметричной (диагональное преобладание присутствует).
*  Хранятся только ненулевые элементы. 
*  Количество итераций ограничено 30000, т.к. предполагается,
*  что если решение не сошлось за 30000 итераций то оно и не сойдётся.
*  Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
*  В качестве предобуславливателя работает неполное разложение Холецкого:
*  K^(-1)==transpose(L)^(-1)*L^(-1); // обращённый предобуславливатель.
* 
*  Удалось частично распараллелить, так что 4 ядерный процессор загружен на 54%
*  К сожалению, некоторые операции не поддаются распараллеливанию.
* 
*  18,07,2021 Метод ICCG присоединен к решению задач механики твёрдого тела матрицы СЛАУ в которых симметричны.
*/
void ICCG(integer iVar, doublereal *dV, doublereal* &x, integer n, bool bprintmessage, bool bdistwall, integer maxiter, int* &inumerate)
{

	// если bdistwall==true то решается СЛАУ для нахождения кратчайшего расстояния до стенки.

	bool bdebug=false;
	if (bdebug) {
		printf("ICCG debug\n");
	}
	else {
		printf("ICCG solver\n");
		printf("1977 Meijerink, Van der Vorst\n");
		printf("for Structural Mechanical simulations.\n");
		printf("first implemented 18.07.2021 in AliceFlow solver.\n");
	}

	doublereal dsize= static_cast<doublereal>(1.0*n); // вещественная длина вектора

	if (bprintmessage) {
		printf("Reshenie metodom sopryjennyh gradientov:\n");
	}
    // матрица СЛАУ
	// в формате CSIR:
	doublereal *adiag=nullptr, *altr=nullptr;
	integer *jptr=nullptr, *iptr=nullptr;

	// предобуславливатель:
	// неполным разложением Холесского в
	// формате CSIR_ITL:
	doublereal *val=nullptr;
	integer *indx=nullptr, *pntr=nullptr;
	
	
	
	
	// инициализация
	// Память выделяется внутри:

	ell_to_CSIR(adiag, altr, jptr, iptr, n);
	ell_to_CSIR_ITLSPD(val, indx, pntr, n);

	std::cout << "n=" << n << "!= ell_global_number_of_strings_26_11_2021=" << ell_global_number_of_strings_26_11_2021 << std::endl;
	//system("pause");

	for (int i = 0; i < ell_global_number_of_strings_26_11_2021; ++i) {
		delete[] data_ell[i];
		delete[] coll_ell[i];
    }
	delete[] data_ell;
	delete[] coll_ell;

	data_ell = nullptr;
	coll_ell = nullptr;

	if (bdebug) {
    	isfinite_vec(pntr[n], val , "val");
	}

	//printf("max memory fics 2...\n"); // debug
	//system("pause");
	//simplesparsefree(M, n);

	integer k=0;
	//integer i; // счётчик
	doublereal *ap=new doublereal[n], *vcopy=new doublereal[n], *f=new doublereal[n],
		 *z=new doublereal[n], *p=new doublereal[n];
    doublereal a, b, res, dbuf;
	

	doublereal dold, dnew;

	// неполное разложение Холецкого:
	// Возвращает левый нижний треугольный сомножитель.
	if (bprintmessage) {
		printf("Incoplete Cholesky decomposition beginig...:\n");
	}
	
	IC0Factor_ITL(val, indx, pntr, n);
	//IC0FactorModify_ITL(val, indx, pntr, n);

	if (bprintmessage) {
		printf("Incoplete Cholesky decomposition finish...:\n");//*/
	}


	// шаг 1.1
	//X0==
	if (x==nullptr) {
        x=new doublereal[n];
		for(integer  i=0;i<n;i++) x[i] = 0.0;
	}

	// пороговое значение невязки
	doublereal e = dterminatedTResudual;
	
	// шаг 1.2
    // вычисление z - невязки начального приближения
	SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, x, ap, n); 
	if (bdebug) {
    	isfinite_vec(n, ap , "ap");
		isfinite_vec(n, dV , "dV");
	}
#pragma omp parallel for
	for (integer  i=0; i<n; ++i) z[i]=dV[i]-ap[i];
	if (bdebug) {
    	isfinite_vec(n, z , "z");
	}
	// передача информации о начальном уровне невязки в текстовый файл:
#pragma omp parallel for
	for (integer  i=0; i<n; ++i) vcopy[i]=z[i];
	if (bdebug) {
    	isfinite_vec(n, vcopy , "vcopy");
	}
    inverseL_ITL(vcopy, val, indx, pntr, f, n, inumerate);
	if (bdebug) {
    	isfinite_vec(n, f , "inverse L: f");
	}
#pragma omp parallel for
    for (integer  i=0; i<n; ++i) vcopy[i]=f[i];
	if (bdebug) {
    	isfinite_vec(n, vcopy , "vcopy");
	}
	inverseU_ITL(vcopy, val, indx, pntr, f, n, inumerate);
	if (bdebug) {
    	isfinite_vec(n, f , "inverse U: f");
	}
    dnew=Scal(z,f,n);
	//dnew=sqrt(dnew)/dsize; // среднеквадратическая эту трогат нельзя она завязана на вычислительный процесс.

	
	// терминаьная невязка всегда на точность аппроксимации меньше стартовой невязки.
	//if (e*dnew<e) e*=dnew;
	doublereal me=sqrt(dnew)/dsize;
	if (me<1.0) e*=me;
	//dterminatedTResudual=e;
	
	
	
	//if (fabs(dnew) > e) {
	if (fabs(me) > e) {

		// шаг 1.3
#pragma omp parallel for
		for (integer  i = 0; i < n; ++i)	p[i] = f[i];
	   res=1000.;
	   while ((fabs(res) > e) && (k < maxiter)) {
		   // шаг 2.1
		  SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, p, ap, n);

		  // шаг 2.2
		  a=dnew/Scal(p,ap,n);// шаговый множитель
		  // шаг 2.3 и 2.4
          #pragma omp parallel for shared(x,z,p,ap,a,n) 
		  for (integer  i=0; i<n; ++i) {
		      x[i]+=a*p[i]; // очередное приближение 
              z[i]-=a*ap[i];// невязка k+1-го приближения
		  }
          #pragma omp parallel for shared(vcopy,z,n) 
          for (integer  i=0; i<n; ++i) vcopy[i]=z[i];
          inverseL_ITL(vcopy, val, indx, pntr, f, n, inumerate);
          #pragma omp parallel for shared(vcopy,f,n) 
          for (integer  i=0; i<n; ++i) vcopy[i]=f[i];
	      inverseU_ITL(vcopy, val, indx, pntr, f, n, inumerate);
		  // шаг 2.5
          dold=dnew;
		  dnew=Scal(z,f,n);

		  // res=dnew; // исходный код.
		  res=sqrt(dnew)/dsize;
		  if (bprintmessage) {
			  
			  if (k % 10 == 0) {
				  if (k % 100 == 0) {
					  std::cout << "iter residual" << std::endl;
				  }
				  std::cout << " " << k << " " << res << std::endl;
			  }
		     
		  }
		  // шаг 3.1
		  b=dnew/dold;
		  // шаг 3.2

          #pragma omp parallel for shared(p,f,b,n) private(dbuf) 
		  for (integer i=0; i<n; ++i) {
			 dbuf=p[i];
		     p[i]=f[i]+b*dbuf; // новое направление минимизации
		  }
          // шаг 3.3
		  k++;
	   } // while

	  
	   // Освобождение памяти
        delete[] ap;
		delete[] vcopy;
		delete[] z;
		delete[] p;
		delete[] f;  
	}
	else {
		// Освобождение памяти
		printf("ICCG inform: residual of the initial approximation is too small me= %e, e=%e...\n",me,e);
		delete[] ap; 
		delete[] vcopy;
		delete[] z; 
		delete[] p;
		delete[] f;		
	}

	// Освобождение памяти
	delete[] adiag; 
	delete[] altr;
	delete[] jptr;
	delete[] iptr;

	delete[] val; 
	delete[] indx;
	delete[] pntr;

#if doubleintprecision == 1
	//printf("pam %lld  \n",k); // контроль количества итераций.
#else
	//printf("pam %d  \n",k); // контроль количества итераций.
#endif
	
	
} // ICCG


// Старая стабильная версия без ILU распараллеливания
/* Метод сопряжённых градиентов Хестенса и Штифеля [1952]
*  Входные параметры:
*  M - разреженная матрица СЛАУ в формате SIMPLESPARSE,
*  dV - вектор правой части,
*  x - начальное приближение к решению или nullptr.
*  n - размерность СЛАУ An*n.
*
*  Разреженная матрица M квадратная размером n*n.
*  Число уравнений равно числу неизвестных и равно n.
*  Матрица M предполагается положительно определённой и
*  симметричной (диагональное преобладание присутствует).
*  Хранятся только ненулевые элементы.
*  Количество итераций ограничено 30000, т.к. предполагается,
*  что если решение не сошлось за 30000 итераций то оно и не сойдётся.
*  Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
*  В качестве предобуславливателя работает неполное разложение Холецкого:
*  K^(-1)==transpose(L)^(-1)*L^(-1); // обращённый предобуславливатель.
*
*  Удалось частично распараллелить, так что 4 ядерный процессор загружен на 54%
*  К сожалению, некоторые операции не поддаются распараллеливанию.
*
*  18,07,2021 Метод ICCG присоединен к решению задач механики твёрдого тела матрицы СЛАУ в которых симметричны.
*/
void ICCG(integer iVar, SIMPLESPARSE& M, doublereal* dV, doublereal*& x, integer n, bool bprintmessage, bool bdistwall, integer maxiter)
{

	// если bdistwall==true то решается СЛАУ для нахождения кратчайшего расстояния до стенки.

	bool bdebug = false;
	if (bdebug) {
		printf("ICCG debug\n");
	}
	else {
		printf("ICCG solver\n");
		printf("1977 Meijerink, Van der Vorst\n");
		printf("for Structural Mechanical simulations.\n");
		printf("first implemented 18.07.2021 in AliceFlow solver.\n");
	}

	doublereal dsize = static_cast<doublereal>(1.0 * n); // вещественная длина вектора

	if (bprintmessage) {
		printf("Reshenie metodom sopryjennyh gradientov:\n");
	}
	// матрица СЛАУ
	// в формате CSIR:
	doublereal* adiag = nullptr, * altr = nullptr;
	integer* jptr = nullptr, * iptr = nullptr;

	// предобуславливатель:
	// неполным разложением Холесского в
	// формате CSIR_ITL:
	doublereal* val = nullptr;
	integer* indx = nullptr, * pntr = nullptr;




	// инициализация
	// Память выделяется внутри:
/*#if !bStableVersion
	ell_to_CSIR(adiag, altr, jptr, iptr, n);
	ell_to_CSIR_ITLSPD(val, indx, pntr, n);



	for (int i = 0; i < n; ++i) {
		delete[] data_ell[i];
		delete[] coll_ell[i];
	}
	delete[] data_ell;
	delete[] coll_ell;

#else*/
	simplesparsetoCSIR(M, adiag, altr, jptr, iptr, n);
	simplesparsetoCSIR_ITLSPD(M, val, indx, pntr, n);
//#endif

	if (bdebug) {
		isfinite_vec(pntr[n], val, "val");
	}

	//printf("max memory fics 2...\n"); // debug
	//system("pause");
	simplesparsefree(M, n);

	integer k = 0;

	doublereal* ap = new doublereal[n], * vcopy = new doublereal[n], * f = new doublereal[n],
		* z = new doublereal[n], * p = new doublereal[n];
	doublereal a, b, res, dbuf;


	doublereal dold, dnew;

	// неполное разложение Холецкого:
	// Возвращает левый нижний треугольный сомножитель.
	if (bprintmessage) {
		printf("Incoplete Cholesky decomposition beginig...:\n");
	}

	IC0Factor_ITL(val, indx, pntr, n);
	//IC0FactorModify_ITL(val, indx, pntr, n);

	if (bprintmessage) {
		printf("Incoplete Cholesky decomposition finish...:\n");//*/
	}


	// шаг 1.1
	//X0==
	if (x == nullptr) {
		x = new doublereal[n];
#pragma omp parallel for
		for (integer i = 0; i < n; ++i) x[i] = 0.0;
	}

	// пороговое значение невязки
	doublereal e = dterminatedTResudual;

	// шаг 1.2
	// вычисление z - невязки начального приближения
	SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, x, ap, n);
	if (bdebug) {
		isfinite_vec(n, ap, "ap");
		isfinite_vec(n, dV, "dV");
	}
#pragma omp parallel for
	for (integer i = 0; i < n; ++i) z[i] = dV[i] - ap[i];
	if (bdebug) {
		isfinite_vec(n, z, "z");
	}
	// передача информации о начальном уровне невязки в текстовый файл:
#pragma omp parallel for
	for (integer i = 0; i < n; ++i) vcopy[i] = z[i];
	if (bdebug) {
		isfinite_vec(n, vcopy, "vcopy");
	}
	inverseL_ITL(vcopy, val, indx, pntr, f, n);
	if (bdebug) {
		isfinite_vec(n, f, "inverse L: f");
	}
#pragma omp parallel for
	for (integer i = 0; i < n; ++i) vcopy[i] = f[i];
	if (bdebug) {
		isfinite_vec(n, vcopy, "vcopy");
	}
	inverseU_ITL(vcopy, val, indx, pntr, f, n);
	if (bdebug) {
		isfinite_vec(n, f, "inverse U: f");
	}
	dnew = Scal(z, f, n);
	//dnew=sqrt(dnew)/dsize; // среднеквадратическая эту трогат нельзя она завязана на вычислительный процесс.


	// терминаьная невязка всегда на точность аппроксимации меньше стартовой невязки.
	//if (e*dnew<e) e*=dnew;
	doublereal me = sqrt(dnew) / dsize;
	if (me < 1.0) e *= me;
	//dterminatedTResudual=e;



	//if (fabs(dnew) > e) {
	if (fabs(me) > e) {

		// шаг 1.3
#pragma omp parallel for
		for (integer i = 0; i < n; ++i)	p[i] = f[i];
		res = 1000.;
		while ((fabs(res) > e) && (k < maxiter)) {
			// шаг 2.1
			SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, p, ap, n);

			// шаг 2.2
			a = dnew / Scal(p, ap, n);// шаговый множитель
			// шаг 2.3 и 2.4
#pragma omp parallel for shared(x,z,p,ap,a,n) schedule (guided)
			for (integer i = 0; i < n; ++i) {
				x[i] += a * p[i]; // очередное приближение 
				z[i] -= a * ap[i];// невязка k+1-го приближения
			}
#pragma omp parallel for shared(vcopy,z,n)  schedule (guided)
			for (integer i = 0; i < n; ++i) vcopy[i] = z[i];
			inverseL_ITL(vcopy, val, indx, pntr, f, n);
#pragma omp parallel for shared(vcopy,f,n) schedule (guided)
			for (integer i = 0; i < n; ++i) vcopy[i] = f[i];
			inverseU_ITL(vcopy, val, indx, pntr, f, n);
			// шаг 2.5
			dold = dnew;
			dnew = Scal(z, f, n);

			// res=dnew; // исходный код.
			res = sqrt(dnew) / dsize;
			if (bprintmessage) {
				if (k % 10 == 0) {
					std::cout << "iter residual" << std::endl;
				}
				std::cout << " " << k << " " << res << std::endl;

			}
			// шаг 3.1
			b = dnew / dold;
			// шаг 3.2

#pragma omp parallel for shared(p,f,b,n) private(dbuf) schedule (guided)
			for (integer i = 0; i < n; ++i) {
				dbuf = p[i];
				p[i] = f[i] + b * dbuf; // новое направление минимизации
			}
			// шаг 3.3
			k++;
		} // while


		// Освобождение памяти
		delete[] ap;
		delete[] vcopy;
		delete[] z;
		delete[] p;
		delete[] f;
	}
	else {
		// Освобождение памяти
		printf("ICCG inform: residual of the initial approximation is too small me= %e, e=%e...\n", me, e);
		delete[] ap;
		delete[] vcopy;
		delete[] z;
		delete[] p;
		delete[] f;
	}

	// Освобождение памяти
	delete[] adiag;
	delete[] altr;
	delete[] jptr;
	delete[] iptr;

	delete[] val;
	delete[] indx;
	delete[] pntr;

#if doubleintprecision == 1
	//printf("pam %lld  \n",k); // контроль количества итераций.
#else
	//printf("pam %d  \n",k); // контроль количества итераций.
#endif


} // ICCG

// алгоритм Ю.Г. Соловейчика [1993]
// для возможно несимметричных матриц.
// Запрограммирован по практикуму
// "Численные методы решения систем уравнений" [2004]
// Новосибирского Государственного Технического Университета (НГТУ).
doublereal* SoloveichikAlgCSIR_SPD(integer isize, // размер квадратной матрицы
						doublereal* &adiag, doublereal* &altr, integer* &jptr, integer* &iptr, // матрица СЛАУ
                         doublereal* &dV,  // вектор правой части
                         const doublereal *dX0, // вектор начального приближения
                         bool bconsole_message) // выводить ли значения невязки на консоль ?
{

     integer i=0,k=0; // счётчики цикла for
     doublereal *dx=nullptr, *dax=nullptr, *dr=nullptr, *dz=nullptr, *dp=nullptr, *dar1=nullptr, *dres=nullptr;
     doublereal dar=0.0, dbr=0.0, dnz=0.0, dscalp=0.0;
	 doublereal kend=1000; // ограничение на максимальное число итераций
	 doublereal epsilon=dterminatedTResudual;  // точность вычисления
	 bool bweShouldContinue=true;


    // Выделение памяти под динамические массивы
    dx=new doublereal[isize]; dax=new doublereal[isize]; dr= new doublereal[isize];
    dz=new doublereal[isize]; dp=new doublereal[isize]; dar1=new doublereal[isize];
	dres=new doublereal[isize]; // вектор результата
   

   // начальное приближение
   // X0 ==
   // под X0 понимается вектор поля температур к примеру.
   if (dX0==nullptr) {
	   for (i=0; i<isize; ++i) dx[i]=0.0;
   }
   else {
	   for (i=0; i<isize; ++i) dx[i]=dX0[i];
   }

   SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dx, dax, isize); // результат занесён в  dax
   for (i=0; i<isize; ++i) dr[i]= dV[i] - dax[i];  // начальная невязка
   dnz=Scal(dr,dr,isize); // начальное значение невязки
   for (i=0; i<isize; ++i) dz[i]=dr[i];  // вектор спуска (сопряжённое направление поиска).
   SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dz, dp, isize); // результат занесён в dp

   if (fabs(Scal( dp, dp, isize))>1e-270) 
   {
      k=1; // итерации начинаются именно с 1
      // начальное значение невязки вычислено выше
      while ((bweShouldContinue) && (k <= kend) && (dnz > epsilon))
	  {
         dscalp=1.0/Scal( dp, dp, isize);
         dar=Scal(dp, dr,isize)*dscalp;
         for (i=0; i<isize; ++i)
		 {
            dx[i]=dx[i]+dar*dz[i];
            dr[i]=dr[i]-dar*dp[i];
		 }
         dnz=dnz-dar*dar/dscalp; // норма невязки
         
         if (bconsole_message) 
		 {
            // печать невязки на консоль
            if ((k % 10) == 0)  printf("iter  residual\n");
#if doubleintprecision == 1
			printf("%lld %e \n", k, dnz);
#else
			printf("%d %e \n", k, dnz);
#endif
           
		 } 
		 SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dr, dar1, isize);// результат занесён в dar1=A*dr
         dbr=-Scal(dp,dar1,isize)*dscalp;
         for (i=0; i<isize; ++i)
		 {
            dz[i]=dr[i]+dbr*dz[i];
            dp[i]=dar1[i]+dbr*dp[i];
		 }
         k++;
         // если процесс расходится то его надо остановить
         if (dnz > 1e7) 
		 {
            // восстановление начального приближения
            for (i=0; i<isize; ++i) if (dX0==nullptr) dx[i]=0.0; else dx[i]=dX0[i];
            printf("\n divergence Soloveichik solver \n");
            bweShouldContinue=false;
            break; // выход из цикла while
		 }
 
	  } // while
      // возвращение результата
      for (i=0; i<isize; ++i) dres[i]=dx[i];
   }
   else
   {
      // возвращает начальное приближение
	   if (dX0 != nullptr) {
		   for (i = 0; i < isize; ++i) dres[i] = dX0[i];
	   }
	   else {
		   printf("error dX0 is nullptr in SoloveichikAlgCSIR_SPD in my_linalg.c\n");
		   system("pause");
		   exit(1);
	   }
   }

   // освобождение памяти выделенной под динамические массивы
   if (dx != nullptr) {
	   delete[] dx;
   }
   if (dax != nullptr) {
	   delete[] dax;
   }
   if (dr != nullptr) {
	   delete[] dr;
   }
   if (dz != nullptr) {
	   delete[] dz;
   }
   if (dp != nullptr) {
	   delete[] dp;
   }
   if (dar1 != nullptr) {
	   delete[] dar1;
   }

   return dres; 

} // SoloveichikAlgCSIR_SPD

// алгоритм Ю.Г. Соловейчика [1993]
// для возможно несимметричных матриц.
// Запрограммирован по практикуму
// "Численные методы решения систем уравнений" [2004]
// Новосибирского Государственного технического университета.
doublereal* SoloveichikAlgCSIR_SPDgood(integer isize, integer nz0,// размер квадратной матрицы
						doublereal* adiag, doublereal* altr, integer* jptr, integer* iptr, // матрица СЛАУ
                         doublereal *dV,  // вектор правой части
                         const doublereal *dX0, // вектор начального приближения
                         bool bconsole_message) // выводить ли значения невязки на консоль ?
{

     integer i,k; // счётчики цикла for
     doublereal *dx, *dax, *dr, *dz, *dp, *dar1, *dres, *df, *vcopy;
     doublereal dar, dbr, dnz, dscalp;
	 doublereal kend=1000; // ограничение на максимальное число итераций
	 doublereal epsilon=dterminatedTResudual;  // точность вычисления
	 bool bweShouldContinue=true;


    // Выделение памяти под динамические массивы
    dx=new doublereal[isize]; dr= new doublereal[isize];
    dz=new doublereal[isize]; dp=new doublereal[isize]; dar1=new doublereal[isize];
	dres=new doublereal[isize]; vcopy=new doublereal[isize]; // вектор результата
	df=new doublereal[isize];
   


	// для неполного разложения Холецкого:
	doublereal  *ldiag=new doublereal[isize], *lltr=new doublereal[nz0];
	integer *jptrsort=new integer[nz0];


    // инициализация
	for (i=0; i<isize; ++i) ldiag[i]=adiag[i];
	for (i=0; i<nz0; ++i) lltr[i]=altr[i];
	// неполное разложение Холецкого:
	// Возвращает левый нижний треугольный сомножитель.
	printf("Incoplete Cholesky decomposition beginig...:\n");
    ICFactor0(ldiag, lltr, jptr, iptr, isize, nz0);
	printf("Incoplete Cholesky decomposition finish...:\n");
    

   /*
	ldiag[0]=1.0; ldiag[1]=1.0;  ldiag[2]=1.838477; ldiag[3]=2.00055;
    ldiag[4]=0.590477; ldiag[5]=1.465913;  ldiag[6]=0.37585673;
	lltr[0]=-1.22383866; lltr[1]=-0.5439282932;  lltr[2]=-1.33247070; lltr[3]=-0.4998613742;
    lltr[4]=-0.1457305633; lltr[5]=-0.06498197865;  lltr[6]=-1.401073265;//*/

    //lltr[0]=-1.22383866; lltr[1]=-1.33247070;  lltr[2]=-0.5439282932; lltr[3]=-0.1457305633;
    //lltr[4]=-0.4998613742; lltr[5]=-1.401073265;  lltr[6]=-0.06498197865;

    for (i=0; i<nz0; ++i) jptrsort[i]=jptr[i];
	//for (i=0; i<isize; ++i) QuickSort(jptrsort, iptr[i], iptr[i+1]-1);

   // начальное приближение
   // X0 ==
   // под X0 понимается вектор поля температур к примеру.
   if (dX0==nullptr) {
	   for (i=0; i<isize; ++i) dx[i]=0.0;
   }
   else {
	   for (i=0; i<isize; ++i) dx[i]=dX0[i];
   }

   //dax=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dx, isize); // результат занесён в  dax
   for (i=0; i<isize; ++i) vcopy[i]=dx[i];
   dax=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, isize);
   for (i=0; i<isize; ++i) vcopy[i]=dax[i];
   SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, dax, isize);
   for (i=0; i<isize; ++i) vcopy[i]=dax[i]; delete dax;
   dax=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, isize);

   for (i=0; i<isize; ++i) vcopy[i]=dV[i]; delete dV;
   dV=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, isize);

   for (i=0; i<isize; ++i) dr[i]= dV[i] - dax[i];  // начальная невязка
   dnz=Scal(dr,dr,isize); // начальное значение невязки
   for (i=0; i<isize; ++i) dz[i]=dr[i];  // вектор спуска (сопряжённое направление поиска).
   //dp=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dz, isize); // результат занесён в dp
   for (i=0; i<isize; ++i) vcopy[i]=dz[i]; 
   dp=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, isize);
   for (i=0; i<isize; ++i) vcopy[i]=dp[i]; 
   SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, dp, isize);
   for (i=0; i<isize; ++i) vcopy[i]=dp[i];
   if (dp != nullptr) {
	   delete[] dp;
	   dp = nullptr;
   }
   dp=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, isize);

   if (fabs(Scal( dp, dp, isize))>1e-270) 
   {
      k=1; // итерации начинаются именно с 1
      // начальное значение невязки вычислено выше
      while ((bweShouldContinue) && (k <= kend) && (fabs(dnz) > epsilon))
	  {
         dscalp=1.0/Scal( dp, dp, isize);
         dar=Scal(dp, dr,isize)*dscalp;
         for (i=0; i<isize; ++i)
		 {
            dx[i]=dx[i]+dar*dz[i];
            dr[i]=dr[i]-dar*dp[i];
		 }
         //dnz=dnz-dar*dar/dscalp; // норма невязки
		 dnz=Scal( dr, dr, isize);
         
         if (bconsole_message) 
		 {
            // печать невязки на консоль
            if ((k % 10) == 0)  printf("iter  residual\n");
#if doubleintprecision == 1
			printf("%lld %e \n", k, dnz);
#else
			printf("%d %e \n", k, dnz);
#endif
            
		 } 
		 //dar1=SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, dr, isize);// результат занесён в dar1=A*dr
          
         for (i=0; i<isize; ++i) vcopy[i]=dr[i];
         dar1=inverseU(vcopy, ldiag, lltr, jptrsort, iptr, isize);
		 for (i=0; i<isize; ++i) vcopy[i]=dar1[i]; 
         SPDMatrixCSIRByVector(adiag, altr, jptr, iptr, vcopy, dar1, isize); 
		 for (i=0; i<isize; ++i) vcopy[i]=dar1[i]; 
		 if (dar1 != nullptr) {
			 delete[] dar1;
			 dar1 = nullptr;
		 }
         dar1=inverseL(vcopy, ldiag, lltr, jptrsort, iptr, isize);


         dbr=-Scal(dp,dar1,isize)*dscalp;
         for (i=0; i<isize; ++i)
		 {
            dz[i]=dr[i]+dbr*dz[i];
            dp[i]=dar1[i]+dbr*dp[i];
		 }
         k++;
         // если процесс расходится то его надо остановить
         if (dnz > 1e7) 
		 {
            // восстановление начального приближения
            for (i=0; i<isize; ++i) if (dX0==nullptr) dx[i]=0.0; else dx[i]=dX0[i];
            printf("\n divergence Soloveichik solver \n");
            bweShouldContinue=false;
            break; // выход из цикла while
		 }
 
	  } // while
      // возвращение результата
      //for (i=0; i<isize; ++i) dres[i]=dx[i];
	  dres=inverseU(dx, ldiag, lltr, jptrsort, iptr, isize);
   }
   else
   {
      // возвращает начальное приближение
	   if (dX0 != nullptr) {
		   for (i = 0; i < isize; ++i) dres[i] = dX0[i];
	   }
   }

   // освобождение памяти выделенной под динамические массивы
   delete[] dx;
   delete[] dax;
   delete[] dr;
   delete[] dz;
   delete[] dp;
   delete[] dar1;
   delete[] vcopy;
   delete[] df;

   return dres; 

} // SoloveichikAlgCSIR_SPDgood

// алгоритм Ю.Г. Соловейчика [1993]
// для возможно несимметричных матриц.
// Запрограммирован по практикуму
// "Численные методы решения систем уравнений" [2004]
// Новосибирского Государственного технического университета.
void SoloveichikAlgCRS(integer isize, // размер квадратной матрицы
						 doublereal* &val, integer* &col_ind, integer* &row_ptr, // матрица СЛАУ
                         doublereal* &dV,  // вектор правой части
                         doublereal* &dX0, // вектор начального приближения
                         bool bconsole_message, integer maxit) // выводить ли значения невязки на консоль ?
{

     
     doublereal *dx=nullptr, *dax=nullptr, *dr=nullptr, *dz=nullptr, *dp=nullptr, *dar1=nullptr, *dres=nullptr, *dstart=nullptr;
     doublereal dar, dbr, dnz, dscalp;
	 doublereal kend= static_cast<doublereal>(maxit); // ограничение на максимальное число итераций
	 doublereal epsilon=dterminatedTResudual;  // точность вычисления
	 bool bweShouldContinue=true;


    // Выделение памяти под динамические массивы
    dx=new doublereal[isize]; dax=new doublereal[isize]; dr= new doublereal[isize];
    dz=new doublereal[isize]; dp=new doublereal[isize]; dar1=new doublereal[isize];
	dres=new doublereal[isize], dstart=new doublereal[isize]; // вектор результата
   

   // начальное приближение
   // X0 ==
   // под X0 понимается вектор поля температур к примеру.
   if (dX0==nullptr) {
#pragma omp parallel for
	   for (integer i=0; i<isize; ++i) {
		   dx[i]=0.0;
		   dstart[i]=0.0;
	   }
	   dX0=new doublereal[isize];

   }
   else {
#pragma omp parallel for
	   for (integer i=0; i<isize; ++i) {
		   dx[i]=dX0[i];
           dstart[i]=dX0[i];
	   }
   }

   
   MatrixCRSByVector(val,col_ind,row_ptr,dx, dax,isize); // результат занесён в  dax
#if doubleintprecision == 1
		//printf("dax=%e, rthsd=%e, %lld\n",Scal(dax,dax,isize),Scal(dV,dV,isize), isize); // debug
#else
		//printf("dax=%e, rthsd=%e, %d\n",Scal(dax,dax,isize),Scal(dV,dV,isize), isize); // debug
#endif
#pragma omp parallel for
   for (integer i=0; i<isize; ++i) dr[i]= dV[i] - dax[i];  // начальная невязка
   dnz=Scal(dr,dr,isize); // начальное значение невязки
   //printf("%e\n",dnz); // debug
#pragma omp parallel for
   for (integer i=0; i<isize; ++i) dz[i]=dr[i];  // вектор спуска (сопряжённое направление поиска).
   MatrixCRSByVector(val,col_ind,row_ptr,dz, dp, isize);// результат занесён в dp

   if (fabs(Scal( dp, dp, isize))>1e-270) 
   {
      integer k=1; // итерации начинаются именно с 1
      // начальное значение невязки вычислено выше

      if ((k==1) && (fabs(dnz) < epsilon)) {
		  printf("residual on a first iteration == %e zero...\n",dnz);
		  //system("pause");
	  }

      while ((bweShouldContinue) && (k <= kend) && (fabs(dnz) > epsilon))
	  {
		  
         dscalp=1.0/Scal( dp, dp, isize);
         dar=Scal(dp, dr,isize)*dscalp;
         #pragma omp parallel for shared(dx,dr,dz,dp,dar,isize) schedule (guided)
         for (integer i=0; i<isize; ++i)
		 {
            dx[i]=dx[i]+dar*dz[i];
            dr[i]=dr[i]-dar*dp[i];
		 }
         dnz=dnz-dar*dar/dscalp; // норма невязки
         
         if (bconsole_message) 
		 {
            // печать невязки на консоль
            if ((k % 10) == 0)  printf("iter  residual\n");
#if doubleintprecision == 1
			printf("%lld %e \n", k, dnz);
			//printf("%lld %e \n",k, NormaChebyshev(dr, isize));
#else
			printf("%d %e \n", k, dnz);
			//printf("%d %e \n",k, NormaChebyshev(dr, isize));
#endif
            
		 } 
		 
		 MatrixCRSByVector(val,col_ind,row_ptr,dr, dar1, isize);// результат занесён в dar1=A*dr
         dbr=-Scal(dp,dar1,isize)*dscalp;
         #pragma omp parallel for shared(isize,dz,dp,dr,dar1,dbr) schedule (guided)
         for (integer i=0; i<isize; ++i)
		 {
            dz[i]=dr[i]+dbr*dz[i];
            dp[i]=dar1[i]+dbr*dp[i];
		 }
         k++;
         // если процесс расходится то его надо остановить
         if (dnz > 1e14) 
		 {
            // восстановление начального приближения
#pragma omp parallel for 
            for (integer i=0; i<isize; ++i) if (dX0==nullptr) dx[i]=0.0; else dx[i]=dstart[i];
            printf("\n divergence Soloveichik solver \n");
			// В этот файл пишется статистика об успешности решения СЛАУ:
			bweShouldContinue=false;
            break; // выход из цикла while
		 }
 
	  } // while
	  
      // возвращение результата
#pragma omp parallel for 
      for (integer i=0; i<isize; ++i) dres[i]=dx[i];
   }
   else
   {
      // возвращает начальное приближение
#pragma omp parallel for
	  for (integer i=0; i<isize; ++i) dres[i]=dstart[i];
	  printf(" (fabs(Scal( dp, dp, isize))>1e-270)==false\n");
	  //system("pause");
   }

   // освобождение памяти выделенной под динамические массивы
   if (dx != nullptr) {
	   delete[] dx;
   }
   if (dax != nullptr) {
	   delete[] dax;
   }
   if (dr != nullptr) {
	   delete[] dr;
   }
   if (dz != nullptr) {
	   delete[] dz;
   }
   if (dp != nullptr) {
	   delete[] dp;
   }
   if (dar1 != nullptr) {
	   delete[] dar1;
   }

   //return dres;
#pragma omp parallel for
   for (integer i=0; i<isize; ++i) dX0[i]=dres[i];
   if (dres != nullptr) {
	   delete[] dres;
   }
   if (dstart != nullptr) {
	   delete[] dstart;
   }

} // SoloveichikAlgCRS


/* Реализация на диннамическом массиве
// инициализирует разреженную матрицу
void initsimplesparse(SIMPLESPARSE &M) {
	M.a=nullptr;
	M.n=0;
	M.incCLUSTER_SIZE=10;
	M.POOL_SIZE=0;
} // initsimplesparse
*/

// Реализация на связном списке
// инициализирует разреженную матрицу
void initsimplesparse(SIMPLESPARSE &M, integer nodes) {
	M.n=0; // изначально все элементы нулевые 
	M.root=new NONZEROELEM*[nodes];
	integer i; // номер строки, номер уравнения в СЛАУ
	for (i=0; i<nodes; ++i) M.root[i]=nullptr; 
} // initsimplesparse

/* Реализация на массиве.
// Добавляет ненулевой элемент в
// простейшую разряженную матрицу M
void addelmsimplesparse(SIMPLESPARSE &M, doublereal aij, integer i, integer j, bool bset) {
	if (M.n==0) {
		// первый элемент
		M.POOL_SIZE+=M.incCLUSTER_SIZE;
		M.n++;
		M.a=new NONZEROELEM[M.POOL_SIZE];
		M.a[0].aij=aij;
		M.a[0].i=i;
		M.a[0].j=j;
	}
	else if (M.n<M.POOL_SIZE) 
	{
		bool flag=false; // элемент не найден
		integer i1; // счётчик
		for (i1=0; i1<M.n; i1++) if ((M.a[i1].i==i) && (M.a[i1].j==j)) {
           flag=true;
           if (bset) M.a[i1].aij=aij;  // установка
		   else M.a[i1].aij+=aij; // добавление
		}
		if (!flag) {
			M.a[M.n].aij=aij;
		    M.a[M.n].i=i;
		    M.a[M.n].j=j;
            M.n++;
		} 
	}
	else // M.n==M.POOL_SIZE
	{
        bool flag=false; // элемент не найден
		integer i1; // счётчик
		for (i1=0; i1<M.n; i1++) if ((M.a[i1].i==i) && (M.a[i1].j==j)) {
           flag=true;
           if (bset) M.a[i1].aij=aij;  // установка
		   else M.a[i1].aij+=aij; // добавление
		}
		if (!flag) {
           NONZEROELEM* list=new NONZEROELEM[M.POOL_SIZE];
		   for (i1=0; i1<M.n; i1++) list[i1]=M.a[i1]; // копирование
		   delete M.a;
		   M.POOL_SIZE+=M.incCLUSTER_SIZE;
		   M.a=new NONZEROELEM[M.POOL_SIZE];
           for (i1=0; i1<M.n; i1++) M.a[i1]=list[i1]; // обратное копирование
           M.a[M.n].aij=aij;
		   M.a[M.n].i=i;
		   M.a[M.n].j=j;
		   M.n++;

		}
	}
} // addelmsimplesparse
*/

// Реализация на связном списке
// Добавляет ненулевой элемент в
// простейшую разряженную матрицу M
// Проверки на равенство добавляемого элемента нулю нет, поэтому
// может добавить и нулевой элемент.
void addelmsimplesparse(SIMPLESPARSE &M, doublereal aij, integer i, integer j, bool bset) {
    NONZEROELEM* p;
	p=M.root[i];
	// линейный поиск элемента с ключом key
	while ((p!=nullptr) && (p->key!=j)) p=p->next;
	if (p!=nullptr) {
		// элемент найден
		if (bset) p->aij=aij; // установка
		else p->aij+=aij; // добавление
	}
	else 
	{
		// если такого элемента нет в списке
		// то добавление элемента в начало списка.
        NONZEROELEM* q=new NONZEROELEM;
		q->aij=aij;
		q->key=j;
		q->next=M.root[i];
		M.root[i]=q;
		q=nullptr;
		M.n++; // количество ненулевых элементов увеличилось на 1. 
	}
} // addelmsimplesparse

// Полная очистка строки i матрицы M.
void addelmsimplesparse_Stress_clean_string(SIMPLESPARSE &M, integer i)
{
	NONZEROELEM* p;
	p = M.root[i];
	if (p != nullptr) {
		NONZEROELEM* q = nullptr;

		q = p->next;
		p->next = nullptr;

		while (q != nullptr) {
			p = q;

			//printf(" Dirichlet p-aij=%d\n",p->aij);
			//system("pause");
			q = p->next;
			p->next = nullptr;
			delete p;
			p = nullptr;
			M.n--;
		}
		delete M.root[i];
		M.root[i] = nullptr;
		M.n--;
	}
}

  // Реализация на связном списке
  // Добавляет ненулевой элемент в
  // простейшую разряженную матрицу M
  // Проверки на равенство добавляемого элемента нулю нет, поэтому
  // может добавить и нулевой элемент.
void addelmsimplesparse_Stress(SIMPLESPARSE &M, doublereal aij, integer i, integer j, bool bset,bool bsetD) {
	const doublereal MY_ZERO_TOLERANCE = 1.0e-300;

	NONZEROELEM* p;
	p = M.root[i];// Корневой элемент строки i
	// линейный поиск элемента с ключом key
	while ((p != nullptr) && (p->key != j)) p = p->next;
	if (p != nullptr) {
		// элемент найден
		if (bsetD) {
			// Удалить строку.
			NONZEROELEM* q = nullptr;
			p = nullptr;
			p = M.root[i];
			q = p->next;
			p->next = nullptr;
			
		
			

			while (q != nullptr) {
				p = q;

				printf(" Dirichlet p-aij=%e\n",p->aij);
				system("pause");
				q = p->next;
				p->next = nullptr;
				delete p;
				p = nullptr;
				M.n--;
			}
			// Установить условие Дирихле равное единице.
			p = M.root[i];
			if (fabs(aij) > MY_ZERO_TOLERANCE) {
				p->aij = aij;
				p->key = j;
			}
			p = nullptr;
		}
		else {
			if (fabs(aij) > MY_ZERO_TOLERANCE) {
				if (bset) p->aij = aij; // установка
				else {
					if (fabs(aij) > MY_ZERO_TOLERANCE) {
						//printf("%e\n", p->aij);
						p->aij += aij; // добавление
						//printf("%e %e\n", p->aij,aij);
						//if (i == 3) {
							//if (fabs(p->aij) < MY_ZERO_TOLERANCE) {
								//printf("i=%d j=%d\n", i, p->key);
								//system("pause");
							//}
						//}
					}
				}
			}
		}
	}
	else
	{
		// если такого элемента нет в списке
		// то добавление элемента в начало списка.
		if (fabs(aij) > MY_ZERO_TOLERANCE) {
			NONZEROELEM* q = new NONZEROELEM;
			q->aij = aij;
			q->key = j;
			q->next = M.root[i];
			M.root[i] = q;
			q = nullptr;
			M.n++; // количество ненулевых элементов увеличилось на 1. 
		}
	}
} // addelmsimplesparse

// освобождение памяти для матрицы SIMPLESPARSE
void simplesparsefree(SIMPLESPARSE &M, integer nodes) {
	integer i; // счётчик цикла for
	for (i=0; i<nodes; ++i) {
        NONZEROELEM* p9, *q9;
		if (M.root != nullptr) {
			p9 = M.root[i]; q9 = p9;
			M.root[i] = nullptr;
			while (p9 != nullptr) {
				p9 = p9->next;
				q9->next = nullptr;
				delete q9;
				q9 = p9;
			}
		}
	}
	if (M.root != nullptr) {
		delete[] M.root;
		M.root = nullptr;
	}
} // simplesparsefree 


// При обнаружении сбоя включить stable версия данной функции.
// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CRS. Всего nodes - уравнений.
template <typename MY_IND_TYPE>
void simplesparsetoCRS(SIMPLESPARSE& M, doublereal*& val, MY_IND_TYPE*& col_ind, MY_IND_TYPE*& row_ptr, integer nodes) {
	//bool flag = true;
	//integer k; // счётчик
	//for (MY_IND_TYPE k = 0; k < nodes; ++k) if (M.root[k] == nullptr) {
		//flag = false; break;
	//}

	//if (flag) {
		val = new doublereal[M.n];
		col_ind = new MY_IND_TYPE[M.n];
		row_ptr = new MY_IND_TYPE[nodes + 1];

		/*bool* bcheck = new bool[M.n];
		for (MY_IND_TYPE i_1 = 0; i_1 < M.n; ++i_1) {
			bcheck[i_1] = false;
		}
		NONZEROELEM* p_1 = nullptr;
		for (MY_IND_TYPE k = 0; k < nodes; ++k) {
			p_1 = M.root[k];
			while (p_1 != nullptr) {
				if (bcheck[p_1->key]) {
					printf("ERROR MATRIX CHECK duplicate ja index string=%lld col_ind=%lld\n", k, p_1->key);
					system("pause");
				}
				bcheck[p_1->key] = true;

				p_1 = p_1->next;
			}
			p_1 = M.root[k];
			// Сброс.
			while (p_1 != nullptr) {
				bcheck[p_1->key] = false;
				p_1 = p_1->next;
			}
		}
		delete[] bcheck;*/

		// инициализация
		for (MY_IND_TYPE k = 0; k < (M.n); ++k) {
			val[k] = 0.0;
			col_ind[k] = 0;
		}
		for (MY_IND_TYPE k = 0; k <= nodes; ++k) {
			row_ptr[k] = static_cast<MY_IND_TYPE>(M.n); // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}

		// Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (k=0; k<M.n; ++k) {
			val[k]=M.a[k].aij;
			col_ind[k]=M.a[k].j;
			row_ptr[M.a[k].i]=myi_min(k,row_ptr[M.a[k].i]);
		}
		*/


		MY_IND_TYPE ik = 0; // счётчик ненулевых элементов СЛАУ
		NONZEROELEM* p;
		MY_IND_TYPE nnz = 0;
		for (MY_IND_TYPE k = 0; k < nodes; ++k) {
			p = M.root[k];
			while (p != nullptr) {
				//if (ik < M.n) {
					// Защита от записи нулевого коэффициента.
					if (fabs(p->aij) > 1.0e-300)
					{
						nnz++;// 20.03.2021
						val[ik] = p->aij;
						/*if (p->key == k) {

							printf("diag=%e string %lld \n", val[ik], k);
							system("pause");
						}*/
						// p->key начинается с нуля
						col_ind[ik] = static_cast<MY_IND_TYPE>(p->key);
						/*if (p->key < 0) {
							printf("%lld\n", row_ptr[k]);
							system("pause");
						}*/
						row_ptr[k] = std::min(ik, row_ptr[k]);
						ik++;
					}
					/*else {
						if (col_ind[ik] == k) {

							printf("diag=%e string %lld \n", val[ik], k);
							system("pause");
						}
					}*/
				//}
				/*else {
					printf("error: ik>=M.n in simplesparsetoCRS\n");
					printf("see module my_linalg.c\n");
					system("pause");
					exit(1);
				}*/
				p = p->next;
			}
		}

		// Исправлена ошибка. Истинное количество ненулевых элементов меньше из за
		// того что берутся лишь ненулевые элементы выше некоторого порога.
		std::cout << "reserve: M.n =" << M.n << ", natural count : nnz = " << nnz  <<"\n";
		//system("pause"); // 20.03.2021
		row_ptr[nodes] = nnz;

		// в каждой строке элементы отсортированы по номерам столбцов:
		for (MY_IND_TYPE k = 0; k < nodes; ++k) {


			if (0) {
				// BubbleSort
				BubbleSortCSIR(col_ind, val, row_ptr[k], row_ptr[k + 1] - 1);
			}
			else {
				QuickSortCSIR(col_ind, val, row_ptr[k], row_ptr[k + 1] - 1);
			}

		}

	//}
} // simplesparsetoCRS

// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CRS. Всего nodes - уравнений.
template <typename MY_IND_TYPE>
void simplesparsetoCRS_stable(SIMPLESPARSE &M, doublereal* &val, MY_IND_TYPE* &col_ind, MY_IND_TYPE* &row_ptr, integer nodes) {
	bool flag=true;
    integer k; // счётчик
	for (k=0; k<nodes; ++k) if (M.root[k]==nullptr) {
		flag=false; break;
	}

	if (flag) {
		val = new doublereal[M.n];
		col_ind = new MY_IND_TYPE[M.n];
		row_ptr = new MY_IND_TYPE[nodes+1];

		bool* bcheck = new bool[M.n];
		for (integer i_1 = 0; i_1 < M.n; ++i_1) {
			bcheck[i_1] = false;
		}
		NONZEROELEM* p_1=nullptr;
		for (k = 0; k < nodes; ++k) {
			p_1 = M.root[k];
			while (p_1 != nullptr) {
				if (bcheck[p_1->key]) {
					printf("ERROR MATRIX CHECK duplicate ja index string=%lld col_ind=%lld\n",k, p_1->key);
					system("pause");
				}
				bcheck[p_1->key] = true;
			
				p_1 = p_1->next;
			}
			p_1 = M.root[k];
			// Сброс.
			while (p_1 != nullptr) {				
					bcheck[p_1->key] = false;
					p_1 = p_1->next;
			}
		}
		delete[] bcheck;
		
		// инициализация
        for (k=0; k<(M.n); ++k) {
		   val[k]=0.0;
		   col_ind[k]=0;
	    }
        for (k=0; k<=nodes; ++k) {
		    row_ptr[k]=M.n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	    }

        // Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (k=0; k<M.n; ++k) {
			val[k]=M.a[k].aij;
            col_ind[k]=M.a[k].j;
            row_ptr[M.a[k].i]=myi_min(k,row_ptr[M.a[k].i]);
		}
		*/


		integer ik=0; // счётчик ненулевых элементов СЛАУ
		NONZEROELEM* p;
		integer nnz = 0;
        for (k=0; k<nodes; ++k) {
			p=M.root[k];
			while (p!=nullptr) {
				if (ik < M.n) {
					// Защита от записи нулевого коэффициента.
					if (fabs(p->aij) > 1.0e-300)
					{
						nnz++;// 20.03.2021
						val[ik] = p->aij;
						/*if (p->key == k) {

							printf("diag=%e string %lld \n", val[ik], k);
							system("pause");
						}*/
						// p->key начинается с нуля
						col_ind[ik] = p->key;
						if (p->key < 0) {
							printf("%lld\n", row_ptr[k]);
							system("pause");
						}
						row_ptr[k] = myi_min(ik, row_ptr[k]);
						ik++;
					}
					else {
						if (col_ind[ik] == k) {

							printf("diag=%e string %lld \n", val[ik], k);
							system("pause");
						}
					}
				}
				else {
					printf("error: ik>=M.n in simplesparsetoCRS\n");
					printf("see module my_linalg.c\n");
					system("pause");
					exit(1);
				}
				p=p->next;
			}
		}

		// Исправлена ошибка. Истинное количество ненулевых элементов меньше из за
		// того что берутся лишь ненулевые элементы выше некоторого порога.
		printf("reserve: M.n =%lld,  natural count: nnz=%lld\n", M.n, nnz);
		//system("pause"); // 20.03.2021
		row_ptr[nodes] = nnz;

		// в каждой строке элементы отсортированы по номерам столбцов:
		for (k = 0; k < nodes; ++k) {

			
			if (0) {
				// BubbleSort
				BubbleSortCSIR(col_ind, val, row_ptr[k], row_ptr[k + 1] - 1);
			}
			else {
				QuickSortCSIR(col_ind, val, row_ptr[k], row_ptr[k + 1] - 1);
			}

		}

	}
} // simplesparsetoCRS_stable

// Преобразует equation3D  формат хранения в CRS формат.
// Цель написания этого преобразователя: экономия оперативной памяти компьютера.
// Т.к. формат SIMPLESPARSE требует слишком много памяти.
integer equation3DtoCRS(equation3D* &sl, equation3D_bon* &slb, doublereal* &val,
	integer* &col_ind, integer* &row_ptr, 
	integer maxelm, integer maxbound, doublereal alpharelax, bool ballocmemory
	, BLOCK* &b, int &lb, SOURCE* &s, int &ls) {

	//integer iproblem_nodes = 0;
	//integer ipatch_problem_nodes = 0;

	// Если ballocmemory равен true то происходит выделение памяти.
	
	bool flag=true;
    integer k; // счётчик
	integer n=0; // число ненулевых элементов

	const doublereal nonzeroEPS = 1.0e-300;// 1e-37; // для отделения вещественного нуля

	// подсчёт количества ненулевых элементов
	// во внутренних точках расчётной области.
	for (k=0; k<maxelm; ++k) {
		
		if (sl[k].ap > nonzeroEPS) n++;
		/*
		//if (fabs(sl[k].ap)> 1e10*nonzeroEPS) n++; // диагональный элемент
		if (sl[k].ap > 1e10*nonzeroEPS) n++; // Диагональный элемент.
			else {
				// 5 августа 2016. 
				iproblem_nodes++;

			flag = false;
			printf("internal zero diagonal element.\n");
#if doubleintprecision == 1
			printf("ap[%lld]=%e, maxelm=%lld\n", k, sl[k].ap, maxelm);
#else
			printf("ap[%d]=%e, maxelm=%d\n", k, sl[k].ap, maxelm);
#endif
			
			printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e sum_nb=%e\n", sl[k].ae, sl[k].aw, sl[k].an, sl[k].as, sl[k].at, sl[k].ab, sl[k].ae + sl[k].aw + sl[k].an + sl[k].as + sl[k].at + sl[k].ab);
			if (sl[k].ap < 0.0) {
				printf("found negativ diagonal coefficient=%e...\n",sl[k].ap);
			}
			printf("fatal error equation3DtoCRS...\n");
			system("pause");
			//09_02_2017
			if (sl[k].ae + sl[k].aw + sl[k].an + sl[k].as + sl[k].at + sl[k].ab > 0.0) {
				sl[k].ap = sl[k].ae + sl[k].aw + sl[k].an + sl[k].as + sl[k].at + sl[k].ab;
				ipatch_problem_nodes++;
			}
			else {
				sl[k].ap = 1.0;
			}
			n++;
			//system("PAUSE");
			//exit(1);
			//n++;
			//sl[k].ap = fabs(sl[k].ae) + fabs(sl[k].aw) + fabs(sl[k].an) + fabs(sl[k].as) + fabs(sl[k].at) + fabs(sl[k].ab);
			
		}
			*/

        if ((sl[k].iE>-1) && (fabs(sl[k].ae) > nonzeroEPS)) n++;
        if ((sl[k].iN>-1) && (fabs(sl[k].an) > nonzeroEPS)) n++;
        if ((sl[k].iT>-1) && (fabs(sl[k].at) > nonzeroEPS))	n++;		
        if ((sl[k].iS>-1) && (fabs(sl[k].as) > nonzeroEPS)) n++;
        if ((sl[k].iW>-1) && (fabs(sl[k].aw) > nonzeroEPS)) n++;
        if ((sl[k].iB>-1) && (fabs(sl[k].ab) > nonzeroEPS)) n++;

		if (sl[k].bE2) {
			if ((sl[k].iE2 > -1) && (fabs(sl[k].ae2) > nonzeroEPS)) n++;
		}
		if (sl[k].bN2) {
			if ((sl[k].iN2 > -1) && (fabs(sl[k].an2) > nonzeroEPS)) n++;
		}
		if (sl[k].bT2) {
			if ((sl[k].iT2 > -1) && (fabs(sl[k].at2) > nonzeroEPS))	n++;
		}
		if (sl[k].bS2) {
			if ((sl[k].iS2 > -1) && (fabs(sl[k].as2) > nonzeroEPS)) n++;
		}
		if (sl[k].bW2) {
			if ((sl[k].iW2 > -1) && (fabs(sl[k].aw2) > nonzeroEPS)) n++;
		}
		if (sl[k].bB2) {
			if ((sl[k].iB2 > -1) && (fabs(sl[k].ab2) > nonzeroEPS)) n++;
		}

		if (sl[k].bE3) {
			if ((sl[k].iE3 > -1) && (fabs(sl[k].ae3) > nonzeroEPS)) n++;
		}
		if (sl[k].bN3) {
			if ((sl[k].iN3 > -1) && (fabs(sl[k].an3) > nonzeroEPS)) n++;
		}
		if (sl[k].bT3) {
			if ((sl[k].iT3 > -1) && (fabs(sl[k].at3) > nonzeroEPS))	n++;
		}
		if (sl[k].bS3) {
			if ((sl[k].iS3 > -1) && (fabs(sl[k].as3) > nonzeroEPS)) n++;
		}
		if (sl[k].bW3) {
			if ((sl[k].iW3 > -1) && (fabs(sl[k].aw3) > nonzeroEPS)) n++;
		}
		if (sl[k].bB3) {
			if ((sl[k].iB3 > -1) && (fabs(sl[k].ab3) > nonzeroEPS)) n++;
		}

		if (sl[k].bE4) {
			if ((sl[k].iE4 > -1) && (fabs(sl[k].ae4) > nonzeroEPS)) n++;
		}
		if (sl[k].bN4) {
			if ((sl[k].iN4 > -1) && (fabs(sl[k].an4) > nonzeroEPS)) n++;
		}
		if (sl[k].bT4) {
			if ((sl[k].iT4 > -1) && (fabs(sl[k].at4) > nonzeroEPS))	n++;
		}
		if (sl[k].bS4) {
			if ((sl[k].iS4 > -1) && (fabs(sl[k].as4) > nonzeroEPS)) n++;
		}
		if (sl[k].bW4) {
			if ((sl[k].iW4 > -1) && (fabs(sl[k].aw4) > nonzeroEPS)) n++;
		}
		if (sl[k].bB4) {
			if ((sl[k].iB4 > -1) && (fabs(sl[k].ab4) > nonzeroEPS)) n++;
		}


	}

	// подсчёт количества ненулевых элементов
    // для граничных точек расчётной области.
	for (k=0; k<maxbound; ++k) {
		if (fabs(slb[k].aw)>nonzeroEPS) n++; // диагональный элемент
		else {
			flag = false;
			printf("boundary zero diagonal element.\n");
		}

		if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) n++;
	}

	if (flag) {
		// memory +15N
		// Теперь выделение памяти будет происходить централизованно, вне данного кода.
		// Это сделано для кода BICGSTAB_internal3. дата изменения 12 апреля 2013.
		// Другой код, использующий equation3dtoCRS может оказаться неработоспособным после этого изменения.
		if ( ballocmemory) {
			// Важно выделить память с запасом, т.к. одна и таже память используется и для компонент скорости и для поправки давления.
		  // Это выделение оперативной памяти было актуально на 7 точечном шаблоне.
			// val = new doublereal[7*(maxelm+maxbound)+2*maxbound+2];
		  // col_ind = new integer[7*(maxelm+maxbound)+2*maxbound+2];
		   //val = new doublereal[n+2];
		   //col_ind = new integer[n+2];
			// Выделение оперативной памяти пригодное и для АЛИС сетки.
			// 2 * maxbound + 2 - это запас.
			// 26.09.2016 Сентябрь 2016 года.
			val = new doublereal[n + 2 * maxbound + 2];
			col_ind = new integer[n + 2 * maxbound + 2];
		    row_ptr = new integer[(maxelm+maxbound)+1];
		    if ((val==nullptr)||(col_ind==nullptr)||(row_ptr==nullptr)) {
			     // недостаточно памяти на данном оборудовании.
			     printf("Problem: not enough memory on your equipment...\n");
				 printf("Please any key to exit...\n");
				 exit(1);
			}
		}

		
		// инициализация
        for (k=0; k<(n); ++k) {
		   val[k]=0.0;
		   col_ind[k]=-1;
	    }
        for (k=0; k<=(maxelm+maxbound); ++k) {
		    row_ptr[k]=n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	    }

		// n - это в данном контексте nnz , т.е. число ненулевых элементов в матрице.
		

        // Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (k=0; k<M.n; ++k) {
			val[k]=M.a[k].aij;
            col_ind[k]=M.a[k].j;
            row_ptr[M.a[k].i]=myi_min(k,row_ptr[M.a[k].i]);
		}
		*/
		integer ik=0; // счётчик ненулевых элементов СЛАУ
		
		// для внутренних узлов расчётной области:
        for (k=0; k<maxelm; ++k) {

			if (fabs(sl[k].ap) > nonzeroEPS) {
                val[ik]=sl[k].ap/alpharelax;
				col_ind[ik]=sl[k].iP;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iE>-1) && (fabs(sl[k].ae) > nonzeroEPS)) {
                val[ik]=-sl[k].ae;
				col_ind[ik]=sl[k].iE;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iN>-1) && (fabs(sl[k].an) > nonzeroEPS)) {
                val[ik]=-sl[k].an;
				col_ind[ik]=sl[k].iN;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iT>-1) && (fabs(sl[k].at) > nonzeroEPS)) {
                val[ik]=-sl[k].at;
				col_ind[ik]=sl[k].iT;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}		
			if ((sl[k].iS>-1) && (fabs(sl[k].as) > nonzeroEPS)) {
                val[ik]=-sl[k].as;
				col_ind[ik]=sl[k].iS;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iW>-1) && (fabs(sl[k].aw) > nonzeroEPS)) {
				val[ik]=-sl[k].aw;
				col_ind[ik]=sl[k].iW;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iB>-1) && (fabs(sl[k].ab) > nonzeroEPS)) {
				val[ik]=-sl[k].ab;
				col_ind[ik]=sl[k].iB;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}

			if (sl[k].bE2) {
				if ((sl[k].iE2 > -1) && (fabs(sl[k].ae2) > nonzeroEPS)) {
					val[ik] = -sl[k].ae2;
					col_ind[ik] = sl[k].iE2;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bN2) {
				if ((sl[k].iN2 > -1) && (fabs(sl[k].an2) > nonzeroEPS)) {
					val[ik] = -sl[k].an2;
					col_ind[ik] = sl[k].iN2;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bT2) {
				if ((sl[k].iT2 > -1) && (fabs(sl[k].at2) > nonzeroEPS)) {
					val[ik] = -sl[k].at2;
					col_ind[ik] = sl[k].iT2;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bS2) {
				if ((sl[k].iS2 > -1) && (fabs(sl[k].as2) > nonzeroEPS)) {
					val[ik] = -sl[k].as2;
					col_ind[ik] = sl[k].iS2;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bW2) {
				if ((sl[k].iW2 > -1) && (fabs(sl[k].aw2) > nonzeroEPS)) {
					val[ik] = -sl[k].aw2;
					col_ind[ik] = sl[k].iW2;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bB2) {
				if ((sl[k].iB2 > -1) && (fabs(sl[k].ab2) > nonzeroEPS)) {
					val[ik] = -sl[k].ab2;
					col_ind[ik] = sl[k].iB2;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}

			if (sl[k].bE3) {
				if ((sl[k].iE3 > -1) && (fabs(sl[k].ae3) > nonzeroEPS)) {
					val[ik] = -sl[k].ae3;
					col_ind[ik] = sl[k].iE3;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bN3) {
				if ((sl[k].iN3 > -1) && (fabs(sl[k].an3) > nonzeroEPS)) {
					val[ik] = -sl[k].an3;
					col_ind[ik] = sl[k].iN3;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bT3) {
				if ((sl[k].iT3 > -1) && (fabs(sl[k].at3) > nonzeroEPS)) {
					val[ik] = -sl[k].at3;
					col_ind[ik] = sl[k].iT3;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bS3) {
				if ((sl[k].iS3 > -1) && (fabs(sl[k].as3) > nonzeroEPS)) {
					val[ik] = -sl[k].as3;
					col_ind[ik] = sl[k].iS3;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bW3) {
				if ((sl[k].iW3 > -1) && (fabs(sl[k].aw3) > nonzeroEPS)) {
					val[ik] = -sl[k].aw3;
					col_ind[ik] = sl[k].iW3;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bB3) {
				if ((sl[k].iB3 > -1) && (fabs(sl[k].ab3) > nonzeroEPS)) {
					val[ik] = -sl[k].ab3;
					col_ind[ik] = sl[k].iB3;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}

			if (sl[k].bE4) {
				if ((sl[k].iE4 > -1) && (fabs(sl[k].ae4) > nonzeroEPS)) {
					val[ik] = -sl[k].ae4;
					col_ind[ik] = sl[k].iE4;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bN4) {
				if ((sl[k].iN4 > -1) && (fabs(sl[k].an4) > nonzeroEPS)) {
					val[ik] = -sl[k].an4;
					col_ind[ik] = sl[k].iN4;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bT4) {
				if ((sl[k].iT4 > -1) && (fabs(sl[k].at4) > nonzeroEPS)) {
					val[ik] = -sl[k].at4;
					col_ind[ik] = sl[k].iT4;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bS4) {
				if ((sl[k].iS4 > -1) && (fabs(sl[k].as4) > nonzeroEPS)) {
					val[ik] = -sl[k].as4;
					col_ind[ik] = sl[k].iS4;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bW4) {
				if ((sl[k].iW4 > -1) && (fabs(sl[k].aw4) > nonzeroEPS)) {
					val[ik] = -sl[k].aw4;
					col_ind[ik] = sl[k].iW4;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			if (sl[k].bB4) {
				if ((sl[k].iB4 > -1) && (fabs(sl[k].ab4) > nonzeroEPS)) {
					val[ik] = -sl[k].ab4;
					col_ind[ik] = sl[k].iB4;
					row_ptr[k] = myi_min(ik, row_ptr[k]);
					ik++;
				}
			}
			

		}

		// для внутренних узлов расчётной области:
        for (k=0; k<maxbound; ++k) {
			if (fabs(slb[k].aw) > nonzeroEPS) {
               // val[ik]=slb[k].aw/alpharelax;
				val[ik]=slb[k].aw; // релаксация для граничных узлов не применяется.
				/*if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				     // Внимание !!! было произведено тестирование: один вариант был с нижней релаксацией для граничных узлов,
					 // а второй вариант был без нижней релаксации на граничных узлах. Было выяснено, что для сходимости
					 // более благоприятен вариант без нижней релаксации на граничных узлах.
					 // Данное изменение согласовано с функцией solve.

					 val[ik]/=alpharelax; // Если условия Неймана то нижняя релаксация.
				}*/
				col_ind[ik]=slb[k].iW;
                row_ptr[maxelm+k]= myi_min(ik,row_ptr[maxelm+k]);
				ik++;
			}
			if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				val[ik]=-slb[k].ai;
				col_ind[ik]=slb[k].iI;
                row_ptr[maxelm+k]= myi_min(ik,row_ptr[maxelm+k]);
				// Это очень важный вопрос и он требует проверки !
				
				ik++;
			}

		}

		// в каждой строке элементы отсортированы по номерам столбцов:
        for (k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR(col_ind, val, row_ptr[k]+1, row_ptr[k+1]-1); 

		/*
		FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, "matr.txt", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {

		 // debug
		for (k=0; k<=maxelm+maxbound; ++k) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld ",row_ptr[k]);
		#else
			fprintf(fp,"%d ",row_ptr[k]);
		#endif
		    
		}
       fprintf(fp,"\n");
	   for (k=0; k<row_ptr[maxelm+maxbound]; ++k) {
	   #if doubleintprecision == 1
			fprintf(fp, "%e %lld\n",val[k],col_ind[k]);
	   #else
			fprintf(fp, "%e %d\n",val[k],col_ind[k]);
	   #endif
		   
	   }
		
		fclose(fp);
	}
	printf("ready");
	system("pause");
	*/

	}

	integer ierr = 0;

	/*
	if (!flag) {
		printf("Error equation 3D to CRS: zero diagonal element...\n");
#if doubleintprecision == 1
		printf("iproblem_nodes=%lld ipatch_problem_nodes=%lld\n", iproblem_nodes, ipatch_problem_nodes);
#else
		printf("iproblem_nodes=%d ipatch_problem_nodes=%d\n", iproblem_nodes, ipatch_problem_nodes);
#endif
		
		//system("pause");
		if (ipatch_problem_nodes != iproblem_nodes) {
			system("pause");
			ierr = 1;
		}		
	}
	*/

	for (k=0; k<n; ++k) if (col_ind[k]==(-1)) {
#if doubleintprecision == 1
		printf("Error equation3D to CRS in string %lld nnz=%lld.\n", k, row_ptr[n]);
#else
		printf("Error equation3D to CRS in string %d nnz=%d.\n", k, row_ptr[n]);
#endif
		
		//system("pause");
		system("pause");
		ierr = 2;
	}

	for (k=0; k<maxelm+maxbound; ++k) {
		if (val[row_ptr[k]]<nonzeroEPS) {
#if doubleintprecision == 1
			printf("negativ diagonal element equation3DtoCRS %lld\n", k);
#else
			printf("negativ diagonal element equation3DtoCRS %d\n", k);
#endif
			printf("ap[iP]=%e\n", val[row_ptr[k]]);
			printf("maxelm=%lld iP=%lld row_ptr[iP]=%lld\n", maxelm, k, row_ptr[k]);
			//system("pause");
			system("pause");
			ierr = 3;
		}
	}

	if (ierr > 0) {
		for (integer i_7 = 0; i_7 < ls; i_7++) {
			printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", s[i_7].g.xS, s[i_7].g.xE, s[i_7].g.yS, s[i_7].g.yE, s[i_7].g.zS, s[i_7].g.zE);
		}
		system("pause");
	}

	return ierr;
} // equation3DtoCRS

// Преобразует equation3D  формат хранения в CRS формат.
// Цель написания этого преобразователя: экономия оперативной памяти компьютера.
// Т.к. формат SIMPLESPARSE требует слишком много памяти.
// nested desection версия алгоритма.
integer equation3DtoCRSnd(equation3D* &sl, equation3D_bon* &slb, doublereal* &val, 
	integer* &col_ind, integer* &row_ptr, 
					 int maxelm, int maxbound, doublereal alpharelax, 
	bool ballocmemory, int* &ifrontregulationgl, int* &ibackregulationgl
	, BLOCK* &b, int &lb, SOURCE* &s, int &ls) {

	// Если ballocmemory равен true то происходит выделение памяти.
	
	bool flag=true;
    integer k; // счётчик
	integer n=0; // число ненулевых элементов

    const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

	// подсчёт количества ненулевых элементов
	// во внутренних точках расчётной области.
	for (k=0; k<maxelm; ++k) {
		
		if (fabs(sl[k].ap)> nonzeroEPS) n++; // диагональный элемент
		else flag=false;

		if (sl[k].ap != sl[k].ap) {
			printf("NAN or INF in iP=%lld\n", k);
			system("pause");
			exit(1);
		}

        if ((sl[k].iE>-1) && (fabs(sl[k].ae) > nonzeroEPS)) n++;
        if ((sl[k].iN>-1) && (fabs(sl[k].an) > nonzeroEPS)) n++;
        if ((sl[k].iT>-1) && (fabs(sl[k].at) > nonzeroEPS))	n++;		
        if ((sl[k].iS>-1) && (fabs(sl[k].as) > nonzeroEPS)) n++;
        if ((sl[k].iW>-1) && (fabs(sl[k].aw) > nonzeroEPS)) n++;
        if ((sl[k].iB>-1) && (fabs(sl[k].ab) > nonzeroEPS)) n++;

		if (b_on_adaptive_local_refinement_mesh) {
			if ((sl[k].iE2>-1) && (fabs(sl[k].ae2) > nonzeroEPS)) n++;
			if ((sl[k].iN2>-1) && (fabs(sl[k].an2) > nonzeroEPS)) n++;
			if ((sl[k].iT2>-1) && (fabs(sl[k].at2) > nonzeroEPS)) n++;
			if ((sl[k].iS2>-1) && (fabs(sl[k].as2) > nonzeroEPS)) n++;
			if ((sl[k].iW2>-1) && (fabs(sl[k].aw2) > nonzeroEPS)) n++;
			if ((sl[k].iB2>-1) && (fabs(sl[k].ab2) > nonzeroEPS)) n++;

			if ((sl[k].iE3>-1) && (fabs(sl[k].ae3) > nonzeroEPS)) n++;
			if ((sl[k].iN3>-1) && (fabs(sl[k].an3) > nonzeroEPS)) n++;
			if ((sl[k].iT3>-1) && (fabs(sl[k].at3) > nonzeroEPS)) n++;
			if ((sl[k].iS3>-1) && (fabs(sl[k].as3) > nonzeroEPS)) n++;
			if ((sl[k].iW3>-1) && (fabs(sl[k].aw3) > nonzeroEPS)) n++;
			if ((sl[k].iB3>-1) && (fabs(sl[k].ab3) > nonzeroEPS)) n++;

			if ((sl[k].iE4>-1) && (fabs(sl[k].ae4) > nonzeroEPS)) n++;
			if ((sl[k].iN4>-1) && (fabs(sl[k].an4) > nonzeroEPS)) n++;
			if ((sl[k].iT4>-1) && (fabs(sl[k].at4) > nonzeroEPS)) n++;
			if ((sl[k].iS4>-1) && (fabs(sl[k].as4) > nonzeroEPS)) n++;
			if ((sl[k].iW4>-1) && (fabs(sl[k].aw4) > nonzeroEPS)) n++;
			if ((sl[k].iB4>-1) && (fabs(sl[k].ab4) > nonzeroEPS)) n++;
		}
	}

	// подсчёт количества ненулевых элементов
    // для граничных точек расчётной области.
	for (k=0; k<maxbound; ++k) {
		if (fabs(slb[k].aw)>nonzeroEPS) n++; // диагональный элемент
		else flag=false;

		if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) n++;
	}

	if (flag) {
		// memory +15N
		// Теперь выделение памяти будет происходить централизованно, вне данного кода.
		// Это сделано для кода BICGSTAB_internal3. дата изменения 12 апреля 2013.
		// Другой код, использующий equation3dtoCRS может оказаться неработоспособным после этого изменения.
		if ( ballocmemory) {
			// Важно выделить память с запасом, т.к. одна и таже память используется и для компонент скорости и для поправки давления.
		   val = new doublereal[7*(maxelm+maxbound)+2*maxbound+2];
		   col_ind = new integer[7*(maxelm+maxbound)+2*maxbound+2];
		   //val = new doublereal[n+2];
		   //col_ind = new integer[n+2];
		   row_ptr = new integer[(maxelm+maxbound)+1];
		   if ((val==nullptr)||(col_ind==nullptr)||(row_ptr==nullptr)) {
			     // недостаточно памяти на данном оборудовании.
			     printf("Problem: not enough memory on your equipment...\n");
				 printf("Please any key to exit...\n");
				 exit(1);
			}
		}

		
		// инициализация
        for (k=0; k<(n); ++k) {
		   val[k]=0.0;
		   col_ind[k]=-1;
	    }
        for (k=0; k<=(maxelm+maxbound); ++k) {
		    row_ptr[k]=n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	    }

        // Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (k=0; k<M.n; ++k) {
			val[k]=M.a[k].aij;
            col_ind[k]=M.a[k].j;
            row_ptr[M.a[k].i]= myi_min(k,row_ptr[M.a[k].i]);
		}
		*/
		integer ik=0; // счётчик ненулевых элементов СЛАУ
		
		for (integer knew=0; knew<maxelm+maxbound; knew++) {
			k=ifrontregulationgl[knew];
			if (k<maxelm) {

				// для внутренних узлов расчётной области:

				if (fabs(sl[k].ap) > nonzeroEPS) {
                val[ik]=sl[k].ap/alpharelax;
				col_ind[ik]=ibackregulationgl[sl[k].iP];
                row_ptr[knew]= myi_min(ik,row_ptr[knew]);
				ik++;
			}


			if ((sl[k].iE>-1) && (fabs(sl[k].ae) > nonzeroEPS)) {
                val[ik]=-sl[k].ae;
				col_ind[ik]=ibackregulationgl[sl[k].iE];
                row_ptr[knew]= myi_min(ik,row_ptr[knew]);
				ik++;
			}
			if ((sl[k].iN>-1) && (fabs(sl[k].an) > nonzeroEPS)) {
                val[ik]=-sl[k].an;
				col_ind[ik]=ibackregulationgl[sl[k].iN];
                row_ptr[knew]= myi_min(ik,row_ptr[knew]);
				ik++;
			}
			if ((sl[k].iT>-1) && (fabs(sl[k].at) > nonzeroEPS)) {
                val[ik]=-sl[k].at;
				col_ind[ik]=ibackregulationgl[sl[k].iT];
                row_ptr[knew]= myi_min(ik,row_ptr[knew]);
				ik++;
			}		
			if ((sl[k].iS>-1) && (fabs(sl[k].as) > nonzeroEPS)) {
                val[ik]=-sl[k].as;
				col_ind[ik]=ibackregulationgl[sl[k].iS];
                row_ptr[knew]= myi_min(ik,row_ptr[knew]);
				ik++;
			}
			if ((sl[k].iW>-1) && (fabs(sl[k].aw) > nonzeroEPS)) {
				val[ik]=-sl[k].aw;
				col_ind[ik]=ibackregulationgl[sl[k].iW];
                row_ptr[knew]= myi_min(ik,row_ptr[knew]);
				ik++;
			}
			if ((sl[k].iB>-1) && (fabs(sl[k].ab) > nonzeroEPS)) {
				val[ik]=-sl[k].ab;
				col_ind[ik]=ibackregulationgl[sl[k].iB];
                row_ptr[knew]= myi_min(ik,row_ptr[knew]);
				ik++;
			}

			if (b_on_adaptive_local_refinement_mesh) {
				if ((sl[k].iE2>-1) && (fabs(sl[k].ae2) > nonzeroEPS)) {
					val[ik] = -sl[k].ae2;
					col_ind[ik] = ibackregulationgl[sl[k].iE2];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iN2>-1) && (fabs(sl[k].an2) > nonzeroEPS)) {
					val[ik] = -sl[k].an2;
					col_ind[ik] = ibackregulationgl[sl[k].iN2];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iT2>-1) && (fabs(sl[k].at2) > nonzeroEPS)) {
					val[ik] = -sl[k].at2;
					col_ind[ik] = ibackregulationgl[sl[k].iT2];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iS2>-1) && (fabs(sl[k].as2) > nonzeroEPS)) {
					val[ik] = -sl[k].as2;
					col_ind[ik] = ibackregulationgl[sl[k].iS2];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iW2>-1) && (fabs(sl[k].aw2) > nonzeroEPS)) {
					val[ik] = -sl[k].aw2;
					col_ind[ik] = ibackregulationgl[sl[k].iW2];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iB2>-1) && (fabs(sl[k].ab2) > nonzeroEPS)) {
					val[ik] = -sl[k].ab2;
					col_ind[ik] = ibackregulationgl[sl[k].iB2];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}

				if ((sl[k].iE3>-1) && (fabs(sl[k].ae3) > nonzeroEPS)) {
					val[ik] = -sl[k].ae3;
					col_ind[ik] = ibackregulationgl[sl[k].iE3];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iN3>-1) && (fabs(sl[k].an3) > nonzeroEPS)) {
					val[ik] = -sl[k].an3;
					col_ind[ik] = ibackregulationgl[sl[k].iN3];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iT3>-1) && (fabs(sl[k].at3) > nonzeroEPS)) {
					val[ik] = -sl[k].at3;
					col_ind[ik] = ibackregulationgl[sl[k].iT3];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iS3>-1) && (fabs(sl[k].as3) > nonzeroEPS)) {
					val[ik] = -sl[k].as3;
					col_ind[ik] = ibackregulationgl[sl[k].iS3];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iW3>-1) && (fabs(sl[k].aw3) > nonzeroEPS)) {
					val[ik] = -sl[k].aw3;
					col_ind[ik] = ibackregulationgl[sl[k].iW3];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iB3>-1) && (fabs(sl[k].ab3) > nonzeroEPS)) {
					val[ik] = -sl[k].ab3;
					col_ind[ik] = ibackregulationgl[sl[k].iB3];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}

				if ((sl[k].iE4>-1) && (fabs(sl[k].ae4) > nonzeroEPS)) {
					val[ik] = -sl[k].ae4;
					col_ind[ik] = ibackregulationgl[sl[k].iE4];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iN4>-1) && (fabs(sl[k].an4) > nonzeroEPS)) {
					val[ik] = -sl[k].an4;
					col_ind[ik] = ibackregulationgl[sl[k].iN4];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iT4>-1) && (fabs(sl[k].at4) > nonzeroEPS)) {
					val[ik] = -sl[k].at4;
					col_ind[ik] = ibackregulationgl[sl[k].iT4];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iS4>-1) && (fabs(sl[k].as4) > nonzeroEPS)) {
					val[ik] = -sl[k].as4;
					col_ind[ik] = ibackregulationgl[sl[k].iS4];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iW4>-1) && (fabs(sl[k].aw4) > nonzeroEPS)) {
					val[ik] = -sl[k].aw4;
					col_ind[ik] = ibackregulationgl[sl[k].iW4];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
				if ((sl[k].iB4>-1) && (fabs(sl[k].ab4) > nonzeroEPS)) {
					val[ik] = -sl[k].ab4;
					col_ind[ik] = ibackregulationgl[sl[k].iB4];
					row_ptr[knew] = myi_min(ik, row_ptr[knew]);
					ik++;
				}
			}

			}
			else {
				// граничный узел
				k-=maxelm;
// для внутренних узлов расчётной области:

				if (fabs(slb[k].aw) > nonzeroEPS) {
               // val[ik]=slb[k].aw/alpharelax;
				val[ik]=slb[k].aw; // релаксация для граничных узлов не применяется.
				/*if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				     // Внимание !!! было произведено тестирование: один вариант был с нижней релаксацией для граничных узлов,
					 // а второй вариант был без нижней релаксации на граничных узлах. Было выяснено, что для сходимости
					 // более благоприятен вариант без нижней релаксации на граничных узлах.
					 // Данное изменение согласовано с функцией solve.

					 val[ik]/=alpharelax; // Если условия Неймана то нижняя релаксация.
				}*/
				col_ind[ik]=ibackregulationgl[slb[k].iW];
                row_ptr[knew]= myi_min(ik,row_ptr[knew]);
				ik++;
			}
			if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				val[ik]=-slb[k].ai;
				col_ind[ik]=ibackregulationgl[slb[k].iI];
                row_ptr[knew]= myi_min(ik,row_ptr[knew]);
				// Это очень важный вопрос и он требует проверки !
				
				ik++;
			}

			}



		}

		
       

		

		// в каждой строке элементы отсортированы по номерам столбцов:
        for (k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR(col_ind, val, row_ptr[k]+1, row_ptr[k+1]-1); 

		/*
		FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, "matr.txt", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {

		 // debug
		for (k=0; k<=maxelm+maxbound; ++k) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld ",row_ptr[k]);
		#else
			fprintf(fp,"%d ",row_ptr[k]);
		#endif
		   
		}
       fprintf(fp,"\n");
	   for (k=0; k<row_ptr[maxelm+maxbound]; ++k) {
	   #if doubleintprecision == 1
			fprintf(fp, "%e %lld\n",val[k],col_ind[k]);
	   #else
			fprintf(fp, "%e %d\n",val[k],col_ind[k]);
	   #endif
		   
	   }
		
		fclose(fp);
	}
	printf("ready");
	system("pause");
	*/

	}

	integer ierr = 0;

	if (!flag) {
		printf("Error equation 3D to CRS: zero diagonal element...\n");
		//system("pause");
		system("pause");
		ierr = 1;
	}

	for (k=0; k<n; ++k) if (col_ind[k]==(-1)) {
#if doubleintprecision == 1
		printf("Error equation3D to CRS in string %lld nnz=%lld.\n", k, row_ptr[n]);
#else
		printf("Error equation3D to CRS in string %d nnz=%d.\n", k, row_ptr[n]);
#endif
		
		//system("pause");
		system("pause");
		ierr = 2;
	}

	for (k=0; k<maxelm+maxbound; ++k) {
		if (val[row_ptr[k]]<nonzeroEPS) {
#if doubleintprecision == 1
			printf("negativ diagonal element equation3DtoCRS %lld\n", k);
#else
			printf("negativ diagonal element equation3DtoCRS %d\n", k);
#endif
			
			//system("pause");
			system("pause");
			ierr = 3;
		}
	}

	if (ierr > 0) {
		for (integer i_7 = 0; i_7 < ls; i_7++) {
			printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", s[i_7].g.xS, s[i_7].g.xE, s[i_7].g.yS, s[i_7].g.yE, s[i_7].g.zS, s[i_7].g.zE);
			system("pause");
		}
	}

	return ierr;
} // equation3DtoCRSnd





// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CSIR. Всего nodes - уравнений.
// Это работает только для SPD матриц.
// Симметричный положительно определённый случай,
// хранится только нижний треугольник.
void simplesparsetoCSIR(SIMPLESPARSE &M, doublereal* &adiag, doublereal* &altr, integer* &jptr, integer* &iptr, integer nodes) {
	bool flag=true;
    integer k; // счётчик
	for (k=0; k<nodes; ++k) if (M.root[k]==nullptr) {
		flag=false; break;
	}

	if (flag) {
		// поддиагональные элементы в altr хранятся построчно
		integer nz=static_cast<integer>((M.n-nodes)/2); // число ненулевых элементов
		adiag = new doublereal[nodes]; // диагональные элементы
		altr = new doublereal[nz]; // поддиагональные элементы
		jptr = new integer[nz]; // номера столцов для нижнего треугольника
		iptr = new integer[nodes+1]; // указатели на следующую строку

		
		// инициализация
		for (k=0; k<nodes; ++k) adiag[k]=0.0;
        for (k=0; k<(nz); ++k) {
		   altr[k]=0.0;
		   jptr[k]=0;
	    }
        for (k=0; k<=nodes; ++k) {
		    iptr[k]=nz; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	    }

        // Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (k=0; k<M.n; ++k) {
			val[k]=M.a[k].aij;
            col_ind[k]=M.a[k].j;
            row_ptr[M.a[k].i]= myi_min(k,row_ptr[M.a[k].i]);
		}
		*/
		/*
		integer ik=0; // счётчик ненулевых элементов СЛАУ
		NONZEROELEM* p;
        for (k=0; k<nodes; ++k) {
			p=M.root[k];
			while (p!=nullptr) {
				val[ik]=p->aij;
				col_ind[ik]=p->key;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
				p=p->next;
			}
		}
		*/

		integer ik=0, imin=1,k1; // счётчик ненулевых поддиагональных элементов СЛАУ
		bool bvisit;
		NONZEROELEM* p;
        for (k=0; k<nodes; ++k) {
			bvisit=false;
			p=M.root[k];
			while (p!=nullptr) {
				if (p->key==k) {
					adiag[k]=p->aij;
				}
				else if (p->key<k) {
					if (ik<(nz)) {
						altr[ik]=p->aij; // ненулевое значение
					    jptr[ik]=p->key; // номер столбца
					}
					else {
						printf("non simmetric matrix ICCG. simplesparsetoCSIR\n");
						//system("pause");
						system("pause");
					}
					bvisit=true;			   
				}
				imin= myi_min(ik,iptr[k]);
#if doubleintprecision == 1
				//printf("imin=%lld\n",imin);
#else
				//printf("imin=%d\n",imin);
#endif
				
                iptr[k]=imin;
                if (imin==0) for (k1=0; k1<k; k1++) iptr[k1]=0;	
				if (bvisit) { 
					ik++;
					bvisit=false;
				}
				p=p->next;
			}
		}


		for (k=0; k<nodes; ++k) QuickSortCSIR(jptr, altr, iptr[k], iptr[k+1]-1);

	}
} // simplesparsetoCSIR


// печать матрицы в консоль
void printM_and_CSIR(SIMPLESPARSE &sparseM, integer  n) {

	FILE *fp=nullptr;
   
#ifdef MINGW_COMPILLER
	int err = 0;
	fp=fopen64("matrix.txt", "w");
	if (fp == NULL) err = 1;
#else
	errno_t err = 0;
	err = fopen_s(&fp, "matrix.txt", "w");
#endif

	if (((err ) != 0)||(fp==nullptr)) {
		printf("Create File temp Error function printM_and_CSIR in my_linalg.cpp\n");
		//system("pause");
		system("pause");

	}
	else {
		if (fp != nullptr) {

			integer i;
			// печать простейшей формы разреженной матрицы.
			for (i = 0; i < n; ++i) {
				NONZEROELEM* pelm = sparseM.root[i];
				while (pelm != nullptr) {
#if doubleintprecision == 1
					fprintf(fp, "a[%lld][%lld]=%e  ", i, pelm->key, pelm->aij);
#else
					fprintf(fp, "a[%d][%d]=%e  ", i, pelm->key, pelm->aij);
#endif
					
					pelm = pelm->next;
				}
				fprintf(fp, "\n");
			}//*/
			fclose(fp); // закрытие файла
		}

	}
	//system("pause");
	system("pause");
}

// Реализация на связном списке.
// Преобразует простейший формат хранения разреженной матрицы
// в формат CSIR_ITL. Всего nodes - уравнений.
// Это работает только для SPD матриц.
// Симметричный положительно определённый случай,
// хранится только верхний треугольник.
// Память выделяется внутри метода.
void simplesparsetoCSIR_ITLSPD(SIMPLESPARSE &M, doublereal* &val, integer* &indx, integer* &pntr, integer nodes) {
	bool flag=true;
    integer k; // счётчик
	for (k=0; k<nodes; ++k) if (M.root[k]==nullptr) {
		flag=false; break;
	}

	if (flag) {
 
		//printM_and_CSIR(M, nodes); // debug

		// поддиагональные элементы в altr хранятся построчно
		integer nz=static_cast<integer>((M.n-nodes)/2 + nodes); // число ненулевых элементов
		val = new doublereal[nz]; // диагональные элементы и наддиагональные элементы
		indx = new integer[nz]; // номера столцов для нижнего треугольника
		pntr = new integer[nodes+1]; // указатели на следующую строку

		
		// инициализация
        for (k=0; k<(nz); ++k) {
		   val[k]=0.0;
		   indx[k]=0;
	    }
        for (k=0; k<=nodes; ++k) {
		    pntr[k]=nz; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
	    }

        

		integer ik=0; // счётчик ненулевых поддиагональных элементов СЛАУ
		NONZEROELEM* p;
        for (k=0; k<nodes; ++k) {
			
			p=M.root[k];
			while (p!=nullptr) {

				// k - номер диагонального элемента
				if ((p->key>=k)&&(p->key<nodes)) {
					if (ik<(nz)) {
						val[ik]=p->aij; // ненулевое значение
					    indx[ik]=p->key; // номер столбца	
					}
					else {
						std::cout << "\nik=" << ik << " nz=" << nz << " nodes=" << nodes << " k=" << k << std::endl;
						printf(" Error non simmetric matrix ICCG. simplesparsetoCSIR_ITLSPD\n");
					    //system("pause");
						system("pause");
					}
					pntr[k]= myi_min(ik,pntr[k]);

					ik++;
				}

				p=p->next;
			}

		}

		for (k=0; k<nodes; ++k) QuickSortCSIR(indx, val, pntr[k], pntr[k+1]-1);


		/*
		FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, "matr.txt", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
	#if doubleintprecision == 1
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");
		// debug
		for (k=0; k<=nodes; ++k) {
			fprintf(fp,"%lld ",pntr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<pntr[nodes]; ++k) {
			fprintf(fp, "%e %lld\n",val[k],indx[k]);
		}
		fprintf(fp, "nz==%lld\n", nz);
	#else
		// запись заголовка
		fprintf(fp, "TITLE = \"ALICEFLOW0_03\"\n");
		// debug
		for (k=0; k<=nodes; ++k) {
			fprintf(fp,"%d ",pntr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<pntr[nodes]; ++k) {
			fprintf(fp, "%e %d\n",val[k],indx[k]);
		}
		fprintf(fp, "nz==%d\n", nz);
	#endif
		
		
		fclose(fp);
	}
	printf("ready");
	system("pause");
	*/
	}
} // simplesparsetoCSIR_ITLSPD

/* Неполное LU разложение для несимметричных матриц
*  Пример А n*n=
*    9.0 0.0 0.0 3.0 1.0 0.0 1.0
*    0.0 11.0 2.0 1.0 0.0 0.0 2.0 
*    0.0 1.0 10.0 2.0 0.0 0.0 0.0 
*    2.0 1.0 2.0 9.0 1.0 0.0 0.0 
*    1.0 0.0 0.0 1.0 12.0 0.0 1.0 
*    0.0 0.0 0.0 0.0 0.0 8.0 0.0
*    2.0  2.0 0.0 0.0 3.0 0.0 8.0
*-----------------------------------------
*  инициализация (в этом виде данные поступают на вход процедуре)
*  память предполагается выделенной заранее:
*  верхняя треугольная матрица хранится построчно, в каждой строке
*  элементы отсортированы по убыванию номеров столбцов.
*  U_val:   1.0, 1.0, 3.0, 9.0,   2.0, 1.0, 2.0, 11.0,   2.0, 10.0, 1.0, 9.0, 1.0,12.0, 8.0, 8.0
*  U_ind:   6, 4, 3, 0,  6, 3, 2, 1,  3,2, 4,3, 6,4, 5, 6
*  U_ptr:   0, 4, 8, 10, 12, 14, 15, 16
*  нижняя треугольная матрица хранится постолбцово, в каждом столбце
*  элементы отсортированы по убыванию номеров строк.
*  L_val:  2.0, 1.0, 2.0, 9.0,    2.0, 1.0, 1.0, 11.0,  2.0, 10.0, 1.0, 9.0,  3.0, 12.0, 8.0, 8.0
*  L_ind:  6, 4, 3, 0,  6, 3, 2, 1,   3, 2,  4,3,  6, 4, 5, 6
*  L_ptr:  0, 4, 8, 10, 12, 14, 15, 16
*----------------------------------------------
*  Результат ILU разложения:
*  U_val: 1.0, 1.0, 3.0, 9.0, 2.0, 1.0, 2.0, 11.0, 2.0, 10.0, 1.0, 9.0, 1.0, 12.0, 8.0, 8.0.
*  L_val: 0.222, 0.111, 0.222, 1.0, -1.273, 0.091, 0.091, 1.0, 0.2, 1.0, 0.111, 1.0, -0.417, 1.0, 1.0, 1.0.
*/
void ILU0_Decomp_ITL(doublereal* &U_val, integer* &U_ind, integer* &U_ptr, doublereal* &L_val, integer* &L_ind, integer* &L_ptr, integer n)
{
	/*
	// выделение памяти
	integer n=7;
	//doublereal U_val[16] = { 3.0, 1.0, 1.0, 9.0,  2.0, 1.0, 2.0, 11.0, 2.0, 10.0, 1.0, 9.0, 1.0,12.0, 8.0, 8.0};
	//integer U_ind[16] = { 3, 4, 6, 0,  2, 3, 6, 1,  3,2, 4,3, 6,4, 5, 6};
	//integer U_ptr[8] = {0, 4, 8, 10, 12, 14, 15, 16};

	// Отсортированы в порядке убывания по столбцам.

	// verno
	doublereal U_val[16] = {  1.0, 1.0, 3.0, 9.0,     2.0, 1.0, 2.0, 11.0,   2.0, 10.0,   1.0, 9.0, 1.0,12.0, 8.0, 8.0};
	integer U_ind[16] = { 6, 4, 3, 0,  6, 3, 2, 1,  3,2, 4,3, 6,4, 5, 6};
	integer U_ptr[8] = {0, 4, 8, 10, 12, 14, 15, 16};

	//doublereal U_val[16] = {  9.0, 11.0, 2.0, 10.0, 3.0, 1.0, 2.0, 9.0, 1.0, 1.0, 12.0, 8.0, 1.0, 2.0, 1.0, 8.0 };
	//integer U_ind[16] = { 0, 1, 1, 2, 0, 1, 2, 3, 0, 3, 4, 5, 0, 1, 4, 6};
	//integer U_ptr[8] = {0, 1, 2, 4, 8, 11, 15, 16};

	//doublereal L_val[16] = {2.0, 1.0, 2.0, 1.0, 1.0, 2.0, 2.0, 1.0,  3.0};
	//integer L_ind[16] = { 3, 4, 6, 2, 3, 6, 3, 4, 6};
	//integer L_ptr[8] = {0, 3, 6, 7, 8, 8, 8, 8};

	// verno
	doublereal L_val[16] = {2.0, 1.0, 2.0, 9.0,    2.0, 1.0, 1.0, 11.0,  2.0, 10.0, 1.0, 9.0,  3.0, 12.0, 8.0, 8.0};
	integer L_ind[16] = { 6, 4, 3, 0,  6, 3, 2, 1,   3, 2,  4,3,  6, 4, 5, 6};
	integer L_ptr[8] = {0, 4, 8, 10, 12, 14, 15, 16};

     //doublereal L_val[16] = {9.0, 11.0, 1.0, 10.0, 2.0, 1.0, 2.0, 9.0, 1.0, 1.0, 12.0, 8.0, 2.0, 2.0, 3.0, 8.0};
	//integer L_ind[16] = { 0, 1, 1, 2, 0, 1, 2, 3, 0, 3, 4, 5, 0, 1, 4, 6};
	//integer L_ptr[8] = {0, 1, 2, 4, 8, 11, 12, 16};
	*/

	// решение
	integer i, j, qn, pn, rn; 
      for (i = 0; i < n - 1; ++i) {
	     doublereal multiplier = U_val[U_ptr[i+1]-1];
    
	     for (j = L_ptr[i]; j < L_ptr[i+1]; ++j)
	          L_val[j] /= multiplier;
    
	     for (j = U_ptr[i+1]; j < U_ptr[i+2]-1; ++j) {
	         multiplier = U_val[j];
	         qn = j + 1;
	         rn = L_ptr[i+1];
	         for (pn = L_ptr[U_ind[j]]; (pn < L_ptr[U_ind[j] + 1])&&(L_ind[pn] <= i + 1); pn++) {
	              while ((qn < U_ptr[i + 2])&&(U_ind[qn] < L_ind[pn])) qn++;

	              if ((qn < U_ptr[i + 2])&&(L_ind[pn] == U_ind[qn]))
	                     U_val[qn] -= multiplier * L_val[pn];
	         }
	         for (; pn < L_ptr[U_ind[j]+1]; pn++) {
	             while ((rn < L_ptr[i + 2])&&(L_ind[rn] < L_ind[pn]))  rn++;

	             if ((rn < L_ptr[i + 2])&&(L_ind[pn] == L_ind[rn]))
	                    L_val[rn] -= multiplier * L_val[pn];
	         }
	      }
      }
	  L_val[L_ptr[n-1]]=1.0;

	  // сортировка по возрастанию
	  for (i = 0; i < n; ++i) {
          QuickSortCSIR(U_ind, U_val, U_ptr[i], U_ptr[i+1]-1);
          QuickSortCSIR(L_ind, L_val, L_ptr[i], L_ptr[i+1]-1);
	  }

	/*
	printf("Uval: ");
	for (i=0; i<16; ++i) printf("%.3f, ",U_val[i]);
	printf("\n\n Lval: ");
	for (i=0; i<16; ++i) printf("%.3f, ",L_val[i]);
	system("pause");
	exit(0);
	*/
} // ILU0_Decomp_ITL

// 31 марта 2013: данный код не является рабочим, он не согласуется со SPARSKIT2 и не проходит проверку на практике.
// На практике он даёт переполнение и расходимость.
// неполное LU разложение с нулевым заполнением из книги Й. Саада
// Iterative Methods for Sparse linear systems.
// Только для матриц с симметричным портретом.
// на вход подаётся матрица А в CRS формате.
// на выходе матрица luval, jlu в MSR формате в матрице L на диагонали 1.0,
// uptr - указатели на диагональные элементы.
// ja и ia наследуются от матрицы a. Это неверно. У сада есть ещё массив jlu.
void ilu0_Saad(integer n, doublereal* a, integer* ja, integer* ia, doublereal* &luval, integer* &uptr, integer &icode) {
//void ilu0_Saadtest() {
	//integer n=5; // число строк в матрице.
	//doublereal a[12] = {1.0, 2.0,   3.0, 4.0, 5.0,   6.0, 7.0, 8.0, 9.0,  10.0, 11.0,   12.0};
	//integer ja[12] = {0, 3,   0, 1, 3,   0, 2, 3, 4,   2, 3,   4};
	//integer ia[6] = {0, 2, 5, 9, 11, 12};

	//integer n=7; // число строк в матрице.
	//doublereal a[25] = {9.0, 3.0, 1.0, 1.0,   11.0, 2.0, 1.0, 2.0,   1.0, 10.0, 2.0,   2.0, 1.0, 2.0, 9.0, 1.0,   1.0, 1.0, 12.0, 1.0,   8.0,   2.0, 2.0, 3.0, 8.0};
	//integer ja[25] = {0, 3, 4, 6,   1, 2, 3, 6,   1, 2, 3,   0, 1, 2, 3, 4,   0, 3, 4, 6,   5,  0, 1, 4, 6};
	//integer ia[8] = {0, 4, 8, 11, 16, 20, 21, 25};

	// MSR format:
	//doublereal *luval;
	//integer *jlu;
	// указатели на диагональный элемент.
	//integer *uptr;
	integer *iw=nullptr; // рабочий массив длины n.

	//integer icode;

    // *******

	// INPUT:
	// n - dimension of matrix
	// a, ja, ia - sparse matrix in CRS format
	// iw - integer work array of length n
	// OUTPUT:
	// luval - L/U matrices stored together. On return luval,
	//         ja, ia is the combined CSR data structure for 
	//         the LU factors.
	// uptr  - pointer to the diagonal elements in the CRS
	//        data structure luval, ja, ia.
	// icode - integer indicating error code on return
	//         icode == -1: normal return
	//         icode == k: encountered a zero pivot at step k

	luval = new doublereal[ia[n]];
	iw = new integer[n];
	uptr = new integer[n];
	for (integer i87 = 0; i87 < n; i87++) {
		uptr[i87] = -1;
	}

	if ((luval != nullptr) && (iw != nullptr) && (uptr != nullptr)) {

		icode = -1; // Normal return

		integer i = 0;
		// initialize work array iw to zero and luval array to a
		for (i = 0; i < ia[n]; ++i) luval[i] = a[i];

		for (i = 0; i < n; ++i) iw[i] = -1;

		// Main loop
		integer k = 0;
		integer j1 = 0, j2 = 0;
		integer j = 0;
		integer jrow = 0;
		doublereal t1 = 0.0;
		integer jj = 0, jw = 0;
		bool bcont = true;

		k = 0;
		while ((icode == -1) && (k < n)) {

			j1 = ia[k];
			j2 = ia[k + 1] - 1;
			for (j = j1; j <= j2; ++j) iw[ja[j]] = j;
			j = j1;
			jrow = ja[j];

			do {

				bcont = true;
				if (jrow >= k) {// Exit if diagonal element is reached
					// Store pointer to diagonal element
					uptr[k] = j;
					if (j < ia[n]) {
						if ((jrow != k) || (fabs(luval[j]) < 1e-37)) {
							icode = k; // Error: zero pivot
						}
						else {
							if (j < ia[n]) {
								luval[j] = 1.0 / luval[j];
							}
							else {
								printf("error: j>=ia[n] in ilu0_Saad\n");
								system("pause");
								exit(1);
							}
						}
					}
					else {
						printf("error 2: j>=ia[n] in ilu0_Saad\n");
						system("pause");
						exit(1);
					}

					bcont = false; // выход из цикла do
				}
				else {
					// Compute the multiplier for jrow
					if (j < ia[n]) {
						t1 = luval[j] * luval[uptr[jrow]];
						luval[j] = t1;

						for (jj = uptr[jrow] + 1; jj < ia[jrow + 1]; jj++) {
							jw = iw[ja[jj]];
							if (jw != (-1)) luval[jw] -= t1*luval[jj];
						}
					}
					else {
						printf("memory problem in ilu0_Saad. luval[j] in j>=ia[n].");
						system("pause");
						exit(1);
					}
					j++;
					jrow = ja[j];
				}
			} while ((bcont) && (j <= j2));

			if (icode == (-1)) {
				// Refresh all entries of iw to zero.
				for (i = j1; i <= j2; ++i) iw[ja[i]] = -1;
				k++;
			}
		}

	}
	else {
		printf("problem memory allocate for ILU0 in ilu0_Saad in my_linalg.c module.\n");
		system("pause");
		exit(1);
	}

	if (iw != nullptr) {
		delete[] iw;
		iw = nullptr;
	}

    //********
	//for (i=0; i<ia[n]; ++i) printf("%e ",luval[i]);
	//printf("\n");
#if doubleintprecision == 1
	//for (i=0; i<n; ++i) printf("%lld ",uptr[i]);
#else
	//for (i=0; i<n; ++i) printf("%d ",uptr[i]);
#endif
    
	//system("pause");

	/*
	if (icode==(-1)) {
            // ILU предобуславливатель:
            doublereal *U_val, *L_val;
	        integer  *U_ind, *U_ptr, *L_ind, *L_ptr;
			IMatrix xO;
			
			initIMatrix(&xO, n); // инициализация
             
            convertCRStoIMatrix(n, luval, ja, ia, uptr, &xO);
			delete luval; 
			delete uptr;
			convertIMatrixtoCSIR_ILU_ITL(&xO, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);
            freeIMatrix(&xO);
			// сортировка по возрастанию
	        for (i = 0; i < n; ++i) {
                QuickSortCSIR(U_ind, U_val, U_ptr[i], U_ptr[i+1]-1);
                QuickSortCSIR(L_ind, L_val, L_ptr[i], L_ptr[i+1]-1);
				L_val[L_ptr[i]]=1.0; // единица на главной диагонали.
	        }

			  // распечатка получившихся матриц
	         
	         for (i=0; i<U_ptr[n]; ++i) {
		         printf("%e ",U_val[i]);
	         }
	         printf("\n");
			 #if doubleintprecision == 1
				for (i=0; i<U_ptr[n]; ++i) {
					printf("%lld ",U_ind[i]);
				}
				printf("\n");
				for (i=0; i<n+1; ++i) {
					printf("%lld ",U_ptr[i]);
				}
			 #else
				for (i=0; i<U_ptr[n]; ++i) {
					printf("%d ",U_ind[i]);
				}
				printf("\n");
				for (i=0; i<n+1; ++i) {
					printf("%d ",U_ptr[i]);
				}
			 #endif
            
	         printf("\n");
	         system("pause");

			 
	         for (i=0; i<L_ptr[n]; ++i) {
		         printf("%e ",L_val[i]);
	         }
	         printf("\n");
			 #if doubleintprecision == 1
				 for (i=0; i<L_ptr[n]; ++i) {
					printf("%lld ",L_ind[i]);
				 }
				 printf("\n");
				 for (i=0; i<n+1; ++i) {
					 printf("%lld ",L_ptr[i]);
				 }
			 #else
				 for (i=0; i<L_ptr[n]; ++i) {
					 printf("%d ",L_ind[i]);
				 }
				 printf("\n");
				 for (i=0; i<n+1; ++i) {
					 printf("%d ",L_ptr[i]);
				 }
			 #endif
             
	         printf("\n");
	         system("pause");
	         

		}
	*/

} // ilu0_Saad

/* Метод бисопряжённых градиентов
* для возможно несимметричной матрицы А (val, col_ind, row_ptr).
* Запрограммировано по книжке Баландин, Шурина: "Методы
* решения СЛАУ большой размерности".
* dV - правая часть СЛАУ,
* x - начальное приближение к решению или nullptr.
* n - размерность А n*n.
* Количество итераций ограничено для обычных задач 2000.
* для сеток в несколько миллионов узлов 8000 итераций.
* Максимальное число итераций передаётся в переменной maxiter.
* Точность выхода по невязке задаётся в глобальной константе:
*  dterminatedTResudual.
* Иногда метод расходится. Если выбрать другой вектор r_tilda, то 
* процесс может стать сходящимся. Ограничение на выбор вектора r_tilda:
* главное чтобы скалярное произведение Scal(r,r_tilda,n) != 0.0.
*/
void BiSoprGrad(IMatrix *xO, equation3D* &sl, equation3D_bon* &slb,
	            doublereal *dV, doublereal* &x, integer maxelm, integer maxbound,
				bool bSaad, doublereal alpharelax, integer  maxiter,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls)
{
	printf("\nBiConjugate Gradients Method...:\n");

    
	integer n=maxelm+maxbound;
	integer iflag=1; // нужно продолжать.

	// Разреженная матрица СЛАУ
	// в CRS формате.
    doublereal *val=nullptr;
    integer* col_ind=nullptr, *row_ptr=nullptr;
	doublereal dbuf=0.0;

	// преобразование из SIMPLESPARSE формата в CRS формат хранения.
	//simplesparsetoCRS(M, val, col_ind, row_ptr, n);
	equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax,true, b, lb, s_loc, ls);
	// Сначало нужно  проверить надо ли решать СЛАУ,
	// т.к. тривиальное решение может подходить и тогда 
	// последующие действия могут вызвать ошибки.
	doublereal *dax=new doublereal[n];
	doublereal *ri=new doublereal[n];
	MatrixCRSByVector(val,col_ind,row_ptr,x,dax,n);
#pragma omp parallel for
	for (integer i=0; i<n; ++i) ri[i]=dV[i]-dax[i];
	if (dax != nullptr) {
		delete[] dax;
	}
	if (fabs(NormaV(ri,n))<dterminatedTResudual) iflag=0;
	if (ri != nullptr) {
		delete[] ri;
	}
	 //if (iflag) Bi_CGStabCRS((maxelm+maxbound), val, col_ind, row_ptr, dV, x, 8000); // debug equation3DtoCRS
	 // printf("test equation3DtoCRS .../n");
	  // system("pause");
	if (iflag==0) {
		if (val != nullptr) {
			delete[] val;
		}
		if (col_ind != nullptr) {
			delete[] col_ind;
		}
		if (row_ptr != nullptr) {
			delete[] row_ptr;
		}
	}

	if (iflag==1) {

	// ILU предобуславливатель:
    doublereal *U_val=nullptr, *L_val=nullptr;
	integer  *U_ind=nullptr, *U_ptr=nullptr, *L_ind=nullptr, *L_ptr=nullptr;

	if (!bSaad) {
		
		printf("Incoplete LU Decomposition begin...\n");
        convertIMatrixtoCSIR_ILU_ITL(xO, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);// освобождение памяти
	    //printf("max memory BiSoprGrad...\n"); system("pause"); // debug
	    // освобождение оперативной памяти
	    freeIMatrix(xO);

    
	    /* // debug TODO 2
	    printf("TODO 2\n");
	    for (i=0; i<U_ptr[n]; ++i) {
		     printf("%e ",U_val[i]);
	    }
	    printf("\n");
		#if doubleintprecision == 1
			for (i=0; i<U_ptr[n]; ++i) {
				printf("%lld ",U_ind[i]);
			}
			printf("\n");
			for (i=0; i<n+1; ++i) {
				printf("%lld ",U_ptr[i]);
			}
		#else
			for (i=0; i<U_ptr[n]; ++i) {
				printf("%d ",U_ind[i]);
			}
			printf("\n");
			for (i=0; i<n+1; ++i) {
				printf("%d ",U_ptr[i]);
			}
		#endif
       
	    printf("\n");
	    system("pause");
	    */

	    ILU0_Decomp_ITL(U_val, U_ind, U_ptr, L_val, L_ind, L_ptr, n);
	    printf("Incoplete LU Decomposition finish...\n");

	}
	else {
		// ILU(0) разложение из книги Й. Саада
        printf("Incoplete LU Decomposition I.Saad begin...\n");
		freeIMatrix(xO);
		doublereal *luval=nullptr;
		integer *uptr=nullptr;
		integer icode=-1;
        ilu0_Saad(n, val, col_ind, row_ptr, luval, uptr, icode); // ILU(0) разложение
		if (icode==(-1)) {
			IMatrix xO1;
            initIMatrix(&xO1, n); // инициализация

            convertCRStoIMatrix(n, luval, col_ind, row_ptr, uptr, &xO1);
			delete luval; 
			delete uptr;
			convertIMatrixtoCSIR_ILU_ITL(&xO1, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);
            freeIMatrix(&xO1);
			// сортировка по возрастанию
	        for (integer i = 0; i < n; ++i) {
                QuickSortCSIR(U_ind, U_val, U_ptr[i], U_ptr[i+1]-1);
                QuickSortCSIR(L_ind, L_val, L_ptr[i], L_ptr[i+1]-1);
				L_val[L_ptr[i]]=1.0; // единица на главной диагонали.
	        }

		}
		else {
#if doubleintprecision == 1
			printf("Error!!! zero  diagonal elem in %lld string matrix.\n", icode);
#else
			printf("Error!!! zero  diagonal elem in %d string matrix.\n", icode);
#endif
			
			//system("pause");
			system("pause");
			exit(0); // выход из программы.
		}

		printf("Incoplete LU Decomposition I.Saad finish...\n");

	}


	doublereal *r=new doublereal[n], *r_tilda=new doublereal[n];
	doublereal *p=new doublereal[n], *f=new doublereal[n], *p_tilda=new doublereal[n];
	doublereal nz=0.0; // невязка
	doublereal *ap=new doublereal[n], *vcopy=new doublereal[n];
	doublereal a=0.0, b=0.0, dold=0.0, dnew=0.0;

	
	integer k=0; // номер итерации.

	// Начальное приближение:
    //X0==
	if (x==nullptr) {
        x=new doublereal[n];
#pragma omp parallel for
		for(integer i=0;i<n;++i) x[i] = 0.0;
	}

	// пороговое значение невязки
	doublereal e = dterminatedTResudual;

	MatrixCRSByVector(val,col_ind,row_ptr,x,ap,n);
#pragma omp parallel for
	for (integer i=0; i<n; ++i) {
		r[i]=dV[i]-ap[i];
		r_tilda[i]=r[i];
	}

	 // p==M^(-1)*r;
#pragma omp parallel for
    for (integer i=0; i<n; ++i) vcopy[i]=r[i];
    inverseL_ITL(vcopy, L_val, L_ind, L_ptr, p, n);
#pragma omp parallel for
    for (integer i=0; i<n; ++i) vcopy[i]=p[i];  
	inverseU_ITL(vcopy, U_val, U_ind, U_ptr, p, n);

    // p_tilda==M^(-T)*r_tilda;
#pragma omp parallel for
	for (integer i=0; i<n; ++i) vcopy[i]=r_tilda[i];
    inverseL_ITL(vcopy, U_val, U_ind, U_ptr, p_tilda, n);
#pragma omp parallel for
    for (integer i=0; i<n; ++i) vcopy[i]=p_tilda[i];  
	inverseU_ITL(vcopy, L_val, L_ind, L_ptr, p_tilda, n);
	   


	nz=NormaV(r,n); // начальное значение невязки

#pragma omp parallel for
	for (integer i=0; i<n; ++i) vcopy[i]=r[i];
    inverseL_ITL(vcopy, L_val, L_ind, L_ptr, f, n);
#pragma omp parallel for
    for (integer i=0; i<n; ++i) vcopy[i]=f[i];  
	inverseU_ITL(vcopy, U_val, U_ind, U_ptr, f, n);
	// f==M^(-1)*r;
	dold=Scal(f,r_tilda,n); 

    while ((nz>e) && (k<maxiter)) { 
		MatrixCRSByVector(val,col_ind,row_ptr,p,ap,n);

		a=dold/Scal(ap,p_tilda,n);
        #pragma omp parallel for shared(n,x,r,p,ap,a) schedule (guided)
		for (integer i=0; i<n; ++i) {
           x[i]+=a*p[i];
		   r[i]-=a*ap[i];
		}
		if (ap != nullptr) {
			delete[] ap;
		}
		ap=MatrixTransposeCRSByVector(val,col_ind,row_ptr,p_tilda,n);

        #pragma omp parallel for shared(n,r_tilda,ap,a)  schedule (guided)
        for (integer i=0; i<n; ++i) {
			r_tilda[i]-=a*ap[i];
		}

		#pragma omp parallel for shared(n,vcopy,r)  schedule (guided)
        for (integer i=0; i<n; ++i) vcopy[i]=r[i];
        inverseL_ITL(vcopy, L_val, L_ind, L_ptr, f, n);

        #pragma omp parallel for shared(n,vcopy,f)  schedule (guided)
        for (integer i=0; i<n; ++i) vcopy[i]=f[i];
	    inverseU_ITL(vcopy, U_val, U_ind, U_ptr, f, n);

	    // f==M^(-1)*r;
		dnew=Scal(f,r_tilda,n);
		b=dnew/dold;
		dold=dnew;
		// вычисление невязки.
        nz=NormaV(r,n);
		if (k%10==0) printf("iter residual\n");
#if doubleintprecision == 1
		printf(" %lld %e\n", k, nz);
#else
		printf(" %d %e\n", k, nz);
#endif
		

		if ((fabs(b) < 1e-60) || (fabs(nz)>1e10)) {
			// метод Бисопряжённых градиентов иногда расходится.
			printf("\nBiCG divergence detected...\n");
           // system("pause");
			system("pause");
			exit(0); // выход из приложения.
			break; // выход из цикла while
		}

        #pragma omp parallel for shared(n,p,f,b)  schedule (guided)
        for (integer i=0; i<n; ++i) {
			p[i]=f[i]+b*p[i];
		}

		#pragma omp parallel for shared(n,vcopy,r_tilda)  schedule (guided)
		for (integer i=0; i<n; ++i) vcopy[i]=r_tilda[i];
        inverseL_ITL(vcopy, U_val, U_ind, U_ptr, f, n);
        #pragma omp parallel for shared(n,vcopy,f)  schedule (guided)
        for (integer i=0; i<n; ++i) vcopy[i]=f[i];
	    inverseU_ITL(vcopy, L_val, L_ind, L_ptr, f, n);
	    // f==M^(-T)*r_tilda;
		#pragma omp parallel for shared(n,p_tilda,f,b,dbuf) schedule (guided)
        for (integer i=0; i<n; ++i) {
			dbuf=p_tilda[i];
		    p_tilda[i]=f[i]+b*dbuf;
		}

		
		k++; // переход к следующей итерации.
	}

	   // Освобождение памяти
	if (r != nullptr) {
		delete[] r;
	}
	if (r_tilda != nullptr) {
		delete[] r_tilda;
	}
	if (p != nullptr) {
		delete[] p;
	}
	if (p_tilda != nullptr) {
		delete[] p_tilda;
	}
	if (ap != nullptr) {
		delete[] ap;
	}
	if (f != nullptr) {
		delete[] f;
	}
	if (vcopy != nullptr) {
		delete[] vcopy;
	}
	if (U_val != nullptr) {
		delete[] U_val;
	}
	if (U_ind != nullptr) {
		delete[] U_ind;
	}
	if (U_ptr != nullptr) {
		delete[] U_ptr;
	}
	if (L_val != nullptr) {
		delete[] L_val;
	}
	if (L_ind != nullptr) {
		delete[] L_ind;
	}
	if (L_ptr != nullptr) {
		delete[] L_ptr;
	}
	if (val != nullptr) {
		delete[] val;
	}
	if (col_ind != nullptr) {
		delete[] col_ind;
	}
	if (row_ptr != nullptr) {
		delete[] row_ptr;
	}
	

	} // if (iflag) end


} // BiSoprGrad

// алгоритм Ю.Г. Соловейчика [1993]
// для возможно несимметричных матриц.
// Запрограммирован по практикуму
// "Численные методы решения систем уравнений" [2004]
// Новосибирского Государственного технического университета.
// Добавлен ILU0 предобуславливатель. Также есть выбор между ILU0 предобуславливателем
// из книги Й. Саада (bSaad==true) или ILU0 предобуславливателем из библиотеки ITL.
void SoloveichikAlg( IMatrix *xO, equation3D* &sl, equation3D_bon* &slb,// Разреженная матрица СЛАУ
					     integer maxelm, integer maxbound, // число внутренних и граничных КО
                         doublereal *dV,  // вектор правой части
                         doublereal* &dX0, // вектор начального приближения
                         bool bconsole_message, // выводить ли значения невязки на консоль ?
						 bool bSaad, // если bSaad==true то использовать ilu0 разложение из книги Й. Саада иначе использовать ITL ilu0 разложение. 
						 integer imaxiter,// максимально допустимое кол-во итераций
						 doublereal alpharelax,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls)
{
    
	integer isize = xO->n;// размер квадратной матрицы
	 // Разреженная матрица СЛАУ
	 // в CRS формате.
     doublereal *val=nullptr;
     integer* col_ind=nullptr, *row_ptr=nullptr;

	 // преобразование из SIMPLESPARSE формата в CRS формат хранения.
	 //simplesparsetoCRS(M, val, col_ind, row_ptr, isize);
	 equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax,true, b, lb, s_loc, ls);

	 // ILU предобуславливатель:
     doublereal *U_val=nullptr, *L_val=nullptr;
	 integer  *U_ind=nullptr, *U_ptr=nullptr, *L_ind=nullptr, *L_ptr=nullptr;

	 if (!bSaad) {
		
		printf("Incoplete LU Decomposition begin...\n");
        convertIMatrixtoCSIR_ILU_ITL(xO, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);// освобождение памяти
	    //printf("max memory BiSoprGrad...\n"); system("pause"); // debug
	    // освобождение оперативной памяти
	    freeIMatrix(xO);

    
	    /* // debug TODO 2
	    printf("TODO 2\n");
	    for (i=0; i<U_ptr[n]; ++i) {
		     printf("%e ",U_val[i]);
	    }
	    printf("\n");
		#if doubleintprecision == 1
			for (i=0; i<U_ptr[n]; ++i) {
				printf("%lld ",U_ind[i]);
			}
			printf("\n");
			for (i=0; i<n+1; ++i) {
				printf("%lld ",U_ptr[i]);
			}
		#else
			for (i=0; i<U_ptr[n]; ++i) {
				printf("%d ",U_ind[i]);
			}
			printf("\n");
			for (i=0; i<n+1; ++i) {
				printf("%d ",U_ptr[i]);
			}
		#endif
        
	    printf("\n");
	    system("pause");
	    */

	    ILU0_Decomp_ITL(U_val, U_ind, U_ptr, L_val, L_ind, L_ptr, isize);
	    printf("Incoplete LU Decomposition finish...\n");

	}
	else {
		// ILU(0) разложение из книги Й. Саада
        printf("Incoplete LU Decomposition I.Saad begin...\n");
		freeIMatrix(xO);
		doublereal *luval=nullptr;
		integer *uptr=nullptr;
		integer icode=-1;
        ilu0_Saad(isize, val, col_ind, row_ptr, luval, uptr, icode); // ILU(0) разложение
		if (icode==(-1)) {
			IMatrix xO1;
            initIMatrix(&xO1, isize); // инициализация

            convertCRStoIMatrix(isize, luval, col_ind, row_ptr, uptr, &xO1);
			delete[] luval; 
			delete[] uptr;
			convertIMatrixtoCSIR_ILU_ITL(&xO1, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);
            freeIMatrix(&xO1);
			// сортировка по возрастанию
	        for (integer i = 0; i < isize; ++i) {
                QuickSortCSIR(U_ind, U_val, U_ptr[i], U_ptr[i+1]-1);
                QuickSortCSIR(L_ind, L_val, L_ptr[i], L_ptr[i+1]-1);
				L_val[L_ptr[i]]=1.0; // единица на главной диагонали.
	        }

		}
		else {
#if doubleintprecision == 1
			printf("Error!!! zero  diagonal elem in %lld string matrix.\n", icode);
#else
			printf("Error!!! zero  diagonal elem in %d string matrix.\n", icode);
#endif
			
			//system("pause");
			system("pause");
			exit(0); // выход из программы.
		}

		printf("Incoplete LU Decomposition I.Saad finish...\n");

	}


     integer i,k; // счётчики цикла for
     doublereal *dx=nullptr, *dax=nullptr, *dr=nullptr, *dz=nullptr, *dp=nullptr, *dar1=nullptr, *dres=nullptr, *f=nullptr, *vcopy=nullptr;
     doublereal dar, dbr, dnz, dscalp;
	 doublereal kend= static_cast<doublereal>(imaxiter); // ограничение на максимальное число итераций
	 doublereal epsilon=dterminatedTResudual;  // точность вычисления
	 bool bweShouldContinue=true;


    // Выделение памяти под динамические массивы
    dx=new doublereal[isize]; dax=new doublereal[isize]; dr= new doublereal[isize];
    dar1=new doublereal[isize]; vcopy=new doublereal[isize];dp= new doublereal[isize];
	dres=new doublereal[isize]; f=new doublereal[isize]; dz=new doublereal[isize];// вектор результата
   

   // начальное приближение
   // X0 ==
   // под X0 понимается вектор поля температур к примеру.
   if (dX0==nullptr) {
	   dX0=new doublereal[isize];
	   for (i=0; i<isize; ++i) {
		   dx[i]=0.0;
		   dX0[i]=0.0;
	   }
   }
   else {
	   for (i=0; i<isize; ++i) dx[i]=dX0[i];
   }

   
   MatrixCRSByVector(val,col_ind,row_ptr,dx, dax, isize); // результат занесён в  dax
   for (i=0; i<isize; ++i) dr[i]= dV[i] - dax[i];  // начальная невязка
   // dr=L^(-1)*(dV-A*dx);
   for (i=0; i<isize; ++i) vcopy[i]=dr[i]; 
   inverseL_ITL(vcopy, L_val, L_ind, L_ptr, dr, isize);
   dnz=Scal(dr,dr,isize); // начальное значение невязки
   // dz=U^(-1)*dr;
   for (i=0; i<isize; ++i) vcopy[i]=dr[i];  // вектор спуска (сопряжённое направление поиска).
   inverseU_ITL(vcopy, U_val, U_ind, U_ptr, dz, isize);
   // dp=L^(-1)*A*dz;
   MatrixCRSByVector(val,col_ind,row_ptr,dz,dp, isize);// результат занесён в dp
   for (i=0; i<isize; ++i) vcopy[i]=dp[i]; 
   inverseL_ITL(vcopy, L_val, L_ind, L_ptr, dp, isize);

   if (fabs(Scal( dp, dp, isize))>1e-270) 
   {
      k=1; // итерации начинаются именно с 1
      // начальное значение невязки вычислено выше
      while ((bweShouldContinue) && (k <= kend) && (fabs(dnz) > epsilon))
	  {
         dscalp=1.0/Scal( dp, dp, isize);
         dar=Scal(dp, dr,isize)*dscalp;
         for (i=0; i<isize; ++i)
		 {
            dx[i]=dx[i]+dar*dz[i];
            dr[i]=dr[i]-dar*dp[i];
		 }
         dnz=dnz-dar*dar/dscalp; // норма невязки
         
         if (bconsole_message) 
		 {
            // печать невязки на консоль
            if ((k % 10) == 0)  printf("iter  residual\n");
#if doubleintprecision == 1
			printf("%lld %e \n", k, dnz);
#else
			printf("%d %e \n", k, dnz);
#endif
            
		 } 
		 
         // f=U^(-1)*dr;
         for (i=0; i<isize; ++i) vcopy[i]=dr[i];  
         inverseU_ITL(vcopy, U_val, U_ind, U_ptr, f, isize);
         for (i=0; i<isize; ++i) vcopy[i]=f[i]; 
		 MatrixCRSByVector(val,col_ind,row_ptr,vcopy, dar1, isize);// результат занесён в dar1=A*U^(-1)*dr
         for (i=0; i<isize; ++i) vcopy[i]=dar1[i]; 
		 // dar1=L^(-1)*A*U^(-1)*dr;
         inverseL_ITL(vcopy, L_val, L_ind, L_ptr, dar1, isize);

         dbr=-Scal(dp,dar1,isize)*dscalp;
         for (i=0; i<isize; ++i)
		 {
            dz[i]=f[i]+dbr*dz[i];
            dp[i]=dar1[i]+dbr*dp[i];
		 }

         k++;
         // если процесс расходится то его надо остановить
         if (dnz > 1e14) 
		 {
            // восстановление начального приближения
            for (i=0; i<isize; ++i) if (dX0==nullptr) dx[i]=0.0; else dx[i]=dX0[i];
            printf("\n divergence Soloveichik solver \n");
			bweShouldContinue=false;
            break; // выход из цикла while
		 }
 
	  } // while

      // возвращение результата
      for (i=0; i<isize; ++i) dres[i]=dx[i];
   }
   else
   {
      // возвращает начальное приближение
	  for (i=0; i<isize; ++i) dres[i]=dX0[i];
   }

   // освобождение памяти выделенной под динамические массивы
   if (dx != nullptr) {
	   delete[] dx;
   }
   if (dax != nullptr) {
	   delete[] dax;
   }
   if (dr != nullptr) {
	   delete[] dr;
   }
   if (vcopy != nullptr) {
	   delete[] vcopy;
   }
   if (dz != nullptr) {
	   delete[] dz;
   }
   if (dp != nullptr) {
	   delete[] dp;
   }
   if (dar1!=nullptr) {
	   delete[] dar1;
   }
   if (f != nullptr) {
	   delete[] f;
   }
   if (U_val != nullptr) {
	   delete[] U_val;
   }
   if (U_ind != nullptr) {
	   delete[] U_ind;
   }
   if (U_ptr != nullptr) {
	   delete[] U_ptr;
   }
   if (L_val != nullptr) {
	   delete[] L_val;
   }
   if (L_ind != nullptr) {
	   delete[] L_ind;
   }
   if (L_ptr != nullptr) {
	   delete[] L_ptr;
   }
   if (val != nullptr) {
	   delete[] val;
   }
   if (col_ind != nullptr) {
	   delete[] col_ind;
   }
   if (row_ptr != nullptr) {
	   delete[] row_ptr;
   }

   for (i=0; i<isize; ++i) {
		 dX0[i]=dres[i];
	   }
   if (dres != nullptr) {
	   delete[] dres;
   }

} // SoloveichikAlg



// Метод Ван Дер Ворста Bi-CGStabCRS
// работает для возможно несимметричных вещественных матриц.
// Несимметричная матрица СЛАУ передаётся в CRS формате
// A (val, col_ind, row_ptr).
// Метод является комбинацией методов BiCG и GMRES(1). 
void Bi_CGStabCRS(integer n, doublereal *val, integer* col_ind, integer* row_ptr,
	doublereal *dV, doublereal* &dX0, integer maxit)
{

	bool bdebug=false;
	if (bdebug) {
		// т.к. мы находимся в отладочном режиме то уменьшаем количество итераций.
		maxit=200;
		printf("presolve:\n");
	}

	
	if (bdebug) {
    	isfinite_vec(row_ptr[n], val , "val");
	}

	integer iflag=1, icount=0;
	doublereal delta0, deltai=1.0E16;
	doublereal bet, roi;
	doublereal roim1=1.0, al=1.0, wi=1.0;
	doublereal *ri, *roc, *s, *t, *vi, *pi, *dx, *dax;
	doublereal epsilon=dterminatedTResudual;  // точность вычисления
    //printf("%e\n",epsilon); // контрольное значение невязки по которой осуществляется выход из итерационного процесса.
    //system("pause");
	

	ri=new doublereal[n]; roc=new doublereal[n]; s=new doublereal[n]; t=new doublereal[n];
	vi=new doublereal[n]; pi=new doublereal[n]; dx=new doublereal[n]; dax=new doublereal[n];

#pragma omp parallel for
	for (integer i=0; i<n; ++i) {
		s[i]=0.0;
		t[i]=0.0;
		vi[i]=0.0;
		pi[i]=0.0;
	}

    // начальное приближение
    // X0 ==
    // под X0 понимается вектор поля температур к примеру.
    if (dX0==nullptr) {
	   dX0=new doublereal[n];
#pragma omp parallel for 
	   for (integer i=0; i<n; ++i) {
		   dx[i]=0.0;
		   dX0[i]=0.0;
	   }
    }
    else {
#pragma omp parallel for
		for (integer i = 0; i < n; ++i) {
			dx[i] = dX0[i];
		}
    }

    MatrixCRSByVector(val,col_ind,row_ptr,dx,dax, n); // результат занесён в  dax
	if (bdebug) {
    	isfinite_vec(n,dax," dax");
	}
	
#pragma omp parallel for
	for (integer i=0; i<n; ++i) {
		ri[i]=dV[i]-dax[i];
		//roc[i]=ri[i];
		roc[i]=1.0; // Так чтобы скалярное произведение было строго ненулевым
		
	}
	if (bdebug) {
    	isfinite_vec(n, dV, "dV");
		isfinite_vec(n, ri, "ri");
		isfinite_vec(n, roc, "roc");
	}
	
	delta0=NormaV(ri,n);
	/*
	if (bdebug) {
	printf("Norma ri=%e",delta0);
	system("pause");
	if (fabs(delta0)>1e-20) {
		for (i=0; i<n; ++i) {
	    	printf("%e \n",ri[i]);
		    system("pause");
		}
	}
	}*/
	// Если решение сразу хорошее то не считать:
	if (fabs(delta0)<epsilon) iflag=0;

	if (bdebug) {
    	printf("solve:");
	}

	

	while (iflag != 0 && icount < maxit) {

		icount++;

		roi = Scal(roc, ri, n);
		if (bdebug) {
			if (fabs(roi) < 1e-20) {
				if (0) {
					printf("Norma ri=%e", delta0);
					//system("pause");
					system("pause");
					if (fabs(delta0) > 1e-20) {
						for (integer i = 0; i < n; ++i) {
							printf("roc=%e ri=%e\n", roc[i], ri[i]);
							//system("pause");
							system("pause");
						}
					}
				}

				printf("neverno vjbran vector roi\n");
				//system("pause");
				system("pause");
			}
		}

			if (!isfinite(roi)) {
				printf("roi!=roi solution bug. \n");
				system("pause");
			}
			if (fabs(wi) < 1.0e-30) {
				if (fabs(roim1) < 1.0e-30) {
					bet = 1.0;
				}
				else {
					bet = (roi / roim1);
				}
			}
			else {
				if (fabs(roim1) < 1.0e-30) {
					bet = (al / wi);
				}
				else {
					bet = (roi / roim1)*(al / wi);
				}
			}
			if (!isfinite(bet)) {
				printf("bet!=bet solution bug. \n");
				printf("%e %e %e %e\n", roi, roim1, al, wi);
				system("pause");
			}


			if (bdebug) {
				if (!isfinite(bet)) {
					printf("bet is infinity");
					//system("pause");
					system("pause");
				}
			}
#pragma omp parallel for shared(n,pi,ri,vi,wi,bet) schedule (guided)
			for (integer i = 0; i < n; ++i) {
				pi[i] = ri[i] + (pi[i] - vi[i] * wi)*bet;
			}
			if (bdebug) {
				isfinite_vec(n, pi, "pi");
			}

			MatrixCRSByVector(val, col_ind, row_ptr, pi, vi, n);
			if (bdebug) {
				isfinite_vec(n, vi, " vi");
			}
			al = roi / Scal(roc, vi, n);
			if (bdebug) {
				if (al != al) {
					printf("al is infinity: roi=%e, Scal(roc,vi)=%e", roi, Scal(roc, vi, n));
					//system("pause");
					system("pause");
				}
			}
#pragma omp parallel for shared(n,s,ri,vi,al) schedule (guided)
			for (integer i = 0; i < n; ++i) {
				s[i] = ri[i] - al*vi[i];
			}
			if (bdebug) {
				isfinite_vec(n, s, "s");
			}

			MatrixCRSByVector(val, col_ind, row_ptr, s, t, n);
			wi = Scal(t, s, n) / Scal(t, t, n);
			if (bdebug) {
				if (wi != wi) {
					printf("wi is infinity");
					//system("pause");
					system("pause");
				}
			}
#pragma omp parallel for shared(n,dx,al,pi,wi,s,ri,t) schedule (guided)
			for (integer i = 0; i < n; ++i) {
				dx[i] += al*pi[i] + wi*s[i];
				ri[i] = s[i] - wi*t[i];
			}
			if (bdebug) {
				isfinite_vec(n, dx, "dx");
				isfinite_vec(n, ri, "ri");
			}
			deltai = NormaV(ri, n);
			if (bdebug) {
				if (deltai != deltai) {
					printf("deltai is infinity");
					//system("pause");
					system("pause");
				}
			}
			// печать невязки на консоль
			//if ((icount % 10) == 0)  std::cout << "iter  residual" << std::endl;

			//std::cout<< icount << " " << deltai << std::endl;

		//if ((icount % 100)== 0) system("pause");

			if (deltai < epsilon) iflag = 0; // конец вычисления
			else roim1 = roi;
		//}
	}


	if (icount >= maxit - 1) {
		printf("Error !!! problem convergence Bi_CGStabCRS solver...\n");
		printf("You mast change solver...\n");
		system("PAUSE");
	}

	if (icount == 0) {
		//printf("internal: number iterations = %lld \n", icount);
	}
	else {
		//printf("internal: number iterations = %lld , finish residual = %e \n", icount, deltai);
	}
	//system("pause");

    // освобождение памяти
	delete[] ri; delete[] roc; delete[] s; delete[] t;
	delete[] vi; delete[] pi; delete[] dax;

#pragma omp parallel for
	for (integer i=0; i<n; ++i) dX0[i]=dx[i];

	//doublereal tmax = -1.0e30;
	//doublereal tmin = 1.0e30;
	//for (integer i = 0; i < n; ++i) {
		//if (dx[i] > tmax) tmax = dx[i];
		//if (dx[i] < tmin) tmin = dx[i];
	//}
	// Паралельное нахождение максимума и минимума.
	//tmax = get_max_array_elm(dx, n);
	//tmin = get_min_array_elm(dx, n);

	//printf("tmax=%e tmin=%e deltai=%e\n",tmax,tmin,deltai);

	delete[] dx; 


} // Bi_CGStabCRS


  // Метод Ван Дер Ворста Bi-CGStabCRS
  // работает для возможно несимметричных вещественных матриц.
  // Несимметричная матрица СЛАУ передаётся в CRS формате
  // A (val, col_ind, row_ptr).
  // Метод является комбинацией методов BiCG и GMRES(1). 
void Bi_CGStabCRS_smoother(integer n, doublereal *val, integer* col_ind, integer* row_ptr, 
	doublereal *dV, doublereal* &dX0, integer maxit)
{

//#define BICGSTAB_CRS_SMOOTHER

	
#ifdef  BICGSTAB_CRS_SMOOTHER
	bool bdebug = false;

	if (bdebug)
	{
		// т.к. мы находимся в отладочном режиме то уменьшаем количество итераций.
		maxit = 200;
		printf("presolve:\n");
	}


	if (bdebug) {
		isfinite_vec(row_ptr[n], val, "val");
	}
#endif

	integer iflag = 1, icount = 0;
	doublereal delta0, deltai;
	doublereal bet, roi;
	doublereal roim1 = 1.0, al = 1.0, wi = 1.0;
	doublereal *ri, *roc, *s, *t, *vi, *pi, *dx, *dax, *y, *z, *y0;
	doublereal epsilon = dterminatedTResudual;  // точность вычисления
												//printf("%e\n",epsilon); // контрольное значение невязки по которой осуществляется выход из итерационного процесса.
												//system("pause");

	ri = new doublereal[n]; roc = new doublereal[n]; s = new doublereal[n]; t = new doublereal[n];
	vi = new doublereal[n]; pi = new doublereal[n]; dx = new doublereal[n]; dax = new doublereal[n];
	y = new doublereal[n]; z = new doublereal[n]; y0 = new doublereal[n];

#pragma omp parallel for
	for (integer i = 0; i<n; ++i) {
		s[i] = 0.0;
		t[i] = 0.0;
		vi[i] = 0.0;
		pi[i] = 0.0;
		y[i] = 0.0;
		z[i] = 0.0;
		y0[i] = 0.0;
	}

	// начальное приближение
	// X0 ==
	// под X0 понимается вектор поля температур к примеру.
	if (dX0 == nullptr) {
		dX0 = new doublereal[n];
#pragma omp parallel for
		for (integer i = 0; i<n; ++i) {
			dx[i] = 0.0;
			dX0[i] = 0.0;
		}
	}
	else {
#pragma omp parallel for
		for (integer i = 0; i<n; ++i) dx[i] = dX0[i];
	}

	MatrixCRSByVector(val, col_ind, row_ptr, dx, dax, n); // результат занесён в  dax
#ifdef  BICGSTAB_CRS_SMOOTHER
	if (bdebug) {
		isfinite_vec(n, dax, " dax");
	}
#endif

#pragma omp parallel for
	for (integer i = 0; i<n; ++i) {
		ri[i] = dV[i] - dax[i];
		//roc[i]=ri[i];
		roc[i] = 1.0; // Так чтобы скалярное произведение было строго ненулевым

	}

#ifdef  BICGSTAB_CRS_SMOOTHER
	if (bdebug) {
		isfinite_vec(n, dV, "dV");
		isfinite_vec(n, ri, "ri");
		isfinite_vec(n, roc, "roc");
	}
#endif

	delta0 = NormaV(ri, n);
	
#ifdef  BICGSTAB_CRS_SMOOTHER
	if (bdebug) {
		printf("Norma ri=%e",delta0);
		system("pause");
			if (fabs(delta0)>1e-20) {
				for (integer i=0; i<n; ++i) {
					printf("%e \n",ri[i]);
					system("pause");
				}
			}
	}
#endif
	// Если решение сразу хорошее то не считать:
	if (fabs(delta0)<epsilon) iflag = 0;

#ifdef  BICGSTAB_CRS_SMOOTHER
	if (bdebug) {
		printf("solve:");
	}
#endif

	doublereal delta0_ = fabs(delta0);

	while ((iflag != 0 && icount < maxit)||(icount<5)) {

		icount++;

		if (!isfinite(Scal(roc, ri, n))) {
			roi = 0.0;
		}
		else {
			roi = Scal(roc, ri, n);
		}
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			if (fabs(roi)<1e-20) {
				if (0) {
					printf("Norma ri=%e", delta0);
					//system("pause");
					system("pause");
					if (fabs(delta0)>1e-20) {
						for (integer i = 0; i<n; ++i) {
							printf("roc=%e ri=%e\n", roc[i], ri[i]);
							//system("pause");
							system("pause");
						}
					}
				}

				printf("neverno vjbran vector roi\n");
				//system("pause");
				system("pause");
			}
			if (roi!=roi) {
				printf("roi is infinity");
				//system("pause");
				system("pause");
			}
		}
#endif
		if (fabs(wi) < 1.0e-30) {
			if (fabs(roim1) < 1.0e-30) {
				bet = 1.0;
			}
			else {
				bet = (roi / roim1);
			}
		}
		else {
			if (fabs(roim1) < 1.0e-30) {
				bet = (al / wi);
			}
			else {
				bet = (roi / roim1) * (al / wi);
			}
		}
		if ((bet != bet) || (!isfinite(bet))) {
			printf("bet!=bet solution bug. \n");
			printf("%e %e %e %e\n", roi, roim1, al, wi);
			system("pause");
		}

		
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			if (bet != bet) {
				printf("bet is infinity");
				//system("pause");
				system("pause");
			}
		}
#endif
#pragma omp parallel for shared(n,pi,ri,vi,wi,bet) schedule (guided)
		for (integer i = 0; i<n; ++i) {
			pi[i] = ri[i] + (pi[i] - vi[i] * wi)*bet;
		}
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			isfinite_vec(n, pi, "pi");
		}
#endif

		// Ky=pi
#pragma omp parallel for  schedule (guided)
		for (integer i = 0; i < n; ++i) y0[i] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !

		Seidel(n, val, col_ind, row_ptr, y0, pi);

		MatrixCRSByVector(val, col_ind, row_ptr, y0, vi, n);
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			isfinite_vec(n, vi, " vi");
		}
#endif
		if ((fabs(roi) < 1e-30) && (fabs(Scal(roc, vi, n)) < 1e-30)) {
			al = 1.0;
		}
		else if (fabs(roi) < 1e-30) {
			al = 0.0;
		}
		else {
			al = roi / Scal(roc, vi, n);
		}
		if (al != al) {
			printf("roi!=roi solution bug. \n");
			system("pause");
		}

		

#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			if (al != al) {
				printf("al is infinity: roi=%e, Scal(roc,vi)=%e", roi, Scal(roc, vi, n));
				//system("pause");
				system("pause");
			}
		}
#endif
#pragma omp parallel for shared(n,s,ri,vi,al) schedule (guided)
		for (integer i = 0; i<n; ++i) {
			s[i] = ri[i] - al*vi[i];
		}

#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			isfinite_vec(n, s, "s");
		}
#endif
		// Ky=s
#pragma omp parallel for  schedule (guided)
		for (integer i = 0; i < n; ++i) z[i] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !

		Seidel(n, val, col_ind, row_ptr, z, s);

		MatrixCRSByVector(val, col_ind, row_ptr, z, t, n);


#pragma omp parallel for  schedule (guided)
		for (integer i = 0; i < n; ++i) y[i] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !

		Seidel(n, val, col_ind, row_ptr, y, t);

		/*if ((fabs(Scal(t, s, n)) < 1e-30) && (fabs(Scal(t, t, n)) < 1e-30)) {
			wi = 1.0;
		}
		else if (fabs(Scal(t, s, n)) < 1e-30) {
			wi = 0.0;
			}
		else {
			doublereal ts = Scal(t, s, n);
			doublereal tt = Scal(t, t, n);
			if ((ts != ts) || (tt != tt)) {
				wi = 0.0;
			}
			else {
				wi = Scal(t, s, n) / Scal(t, t, n);
			}
		}
		*/
		wi = Scal(y, z, n) / Scal(y, y, n);
		
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			if (wi!=wi) {
				printf("wi is infinity");
				//system("pause");
				system("pause");
			}
		}
#endif
#pragma omp parallel for shared(n,dx,al,pi,wi,s,ri,t) schedule (guided)
		for (integer i = 0; i<n; ++i) {
			dx[i] += al*y0[i] + wi*z[i];
			ri[i] = s[i] - wi*t[i];
		}
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			isfinite_vec(n, dx, "dx");
			isfinite_vec(n, ri, "ri");
		}
#endif
		deltai = NormaV(ri, n);

#ifdef  BICGSTAB_CRS_SMOOTHER

		if (bdebug) {
			if (deltai!=deltai) {
				printf("deltai is infinity");
				//system("pause");
				system("pause");
			}
		}
#endif
		// печать невязки на консоль
		//if ((icount % 10) == 0)  std::cout << "iter  residual" << std::endl;

		//std::cout<< icount <<" "<<deltai<<std::endl;

		//if ((icount % 100)== 0) system("pause");
		if (deltai / delta0_ < 0.8) iflag = 0;
		if (deltai <epsilon) iflag = 0; // конец вычисления
		else roim1 = roi;
	}

	//printf("internal: %lld %e \n", icount, deltai);
	//printf("%d", icount);
	//system("pause");

	// освобождение памяти
	delete[] ri; delete[] roc; delete[] s; delete[] t;
	delete[] vi; delete[] pi; delete[] dax;

#pragma omp parallel for
	for (integer i = 0; i<n; ++i) dX0[i] = dx[i];

	delete[] dx;


	//getchar();

} // Bi_CGStabCRS_smoother



// Метод Ван Дер Ворста Bi-CGStabCRS
  // работает для возможно несимметричных вещественных матриц.
  // Несимметричная матрица СЛАУ передаётся в CRS формате
  // A (val, col_ind, row_ptr).
  // Метод является комбинацией методов BiCG и GMRES(1). 
void Bi_CGStabCRS_smoother_naiv(integer n, doublereal* val, integer* col_ind, integer* row_ptr,
	doublereal* dV, doublereal*& dX0, integer maxit)
{

	//#define BICGSTAB_CRS_SMOOTHER


#ifdef  BICGSTAB_CRS_SMOOTHER
	bool bdebug = false;

	if (bdebug)
	{
		// т.к. мы находимся в отладочном режиме то уменьшаем количество итераций.
		maxit = 200;
		printf("presolve:\n");
	}


	if (bdebug) {
		isfinite_vec(row_ptr[n], val, "val");
	}
#endif

	integer iflag = 1, icount = 0;
	doublereal delta0, deltai;
	doublereal bet, roi;
	doublereal roim1 = 1.0, al = 1.0, wi = 1.0;
	doublereal* ri, * roc, * s, * t, * vi, * pi, * dx, * dax;
	doublereal epsilon = dterminatedTResudual;  // точность вычисления
												//printf("%e\n",epsilon); // контрольное значение невязки по которой осуществляется выход из итерационного процесса.
												//system("pause");

	ri = new doublereal[n]; roc = new doublereal[n]; s = new doublereal[n]; t = new doublereal[n];
	vi = new doublereal[n]; pi = new doublereal[n]; dx = new doublereal[n]; dax = new doublereal[n];

#pragma omp parallel for
	for (integer i = 0; i < n; ++i) {
		s[i] = 0.0;
		t[i] = 0.0;
		vi[i] = 0.0;
		pi[i] = 0.0;
	}

	// начальное приближение
	// X0 ==
	// под X0 понимается вектор поля температур к примеру.
	if (dX0 == nullptr) {
		dX0 = new doublereal[n];
#pragma omp parallel for
		for (integer i = 0; i < n; ++i) {
			dx[i] = 0.0;
			dX0[i] = 0.0;
		}
	}
	else {
#pragma omp parallel for
		for (integer i = 0; i < n; ++i) dx[i] = dX0[i];
	}

	MatrixCRSByVector(val, col_ind, row_ptr, dx, dax, n); // результат занесён в  dax
#ifdef  BICGSTAB_CRS_SMOOTHER
	if (bdebug) {
		isfinite_vec(n, dax, " dax");
	}
#endif

#pragma omp parallel for
	for (integer i = 0; i < n; ++i) {
		ri[i] = dV[i] - dax[i];
		//roc[i]=ri[i];
		roc[i] = 1.0; // Так чтобы скалярное произведение было строго ненулевым

	}

#ifdef  BICGSTAB_CRS_SMOOTHER
	if (bdebug) {
		isfinite_vec(n, dV, "dV");
		isfinite_vec(n, ri, "ri");
		isfinite_vec(n, roc, "roc");
	}
#endif

	delta0 = NormaV(ri, n);

#ifdef  BICGSTAB_CRS_SMOOTHER
	if (bdebug) {
		printf("Norma ri=%e", delta0);
		system("pause");
		if (fabs(delta0) > 1e-20) {
			for (integer i = 0; i < n; ++i) {
				printf("%e \n", ri[i]);
				system("pause");
			}
		}
	}
#endif
	// Если решение сразу хорошее то не считать:
	if (fabs(delta0) < epsilon) iflag = 0;

#ifdef  BICGSTAB_CRS_SMOOTHER
	if (bdebug) {
		printf("solve:");
	}
#endif

	doublereal delta0_ = fabs(delta0);

	while ((iflag != 0 && icount < maxit) || (icount < 5)) {

		icount++;

		roi = Scal(roc, ri, n);
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			if (fabs(roi) < 1e-20) {
				if (0) {
					printf("Norma ri=%e", delta0);
					//system("pause");
					system("pause");
					if (fabs(delta0) > 1e-20) {
						for (integer i = 0; i < n; ++i) {
							printf("roc=%e ri=%e\n", roc[i], ri[i]);
							//system("pause");
							system("pause");
						}
					}
				}

				printf("neverno vjbran vector roi\n");
				//system("pause");
				system("pause");
			}
			if (roi != roi) {
				printf("roi is infinity");
				//system("pause");
				system("pause");
			}
		}
#endif
		bet = (roi / roim1) * (al / wi);
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			if (bet != bet) {
				printf("bet is infinity");
				//system("pause");
				system("pause");
			}
		}
#endif
#pragma omp parallel for shared(n,pi,ri,vi,wi,bet) schedule (guided)
		for (integer i = 0; i < n; ++i) {
			pi[i] = ri[i] + (pi[i] - vi[i] * wi) * bet;
		}
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			isfinite_vec(n, pi, "pi");
		}
#endif

		MatrixCRSByVector(val, col_ind, row_ptr, pi, vi, n);
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			isfinite_vec(n, vi, " vi");
		}
#endif
		al = roi / Scal(roc, vi, n);

#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			if (al != al) {
				printf("al is infinity: roi=%e, Scal(roc,vi)=%e", roi, Scal(roc, vi, n));
				//system("pause");
				system("pause");
			}
		}
#endif
#pragma omp parallel for shared(n,s,ri,vi,al) schedule (guided)
		for (integer i = 0; i < n; ++i) {
			s[i] = ri[i] - al * vi[i];
		}

#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			isfinite_vec(n, s, "s");
		}
#endif


		MatrixCRSByVector(val, col_ind, row_ptr, s, t, n);
		wi = Scal(t, s, n) / Scal(t, t, n);
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			if (wi != wi) {
				printf("wi is infinity");
				//system("pause");
				system("pause");
			}
		}
#endif
#pragma omp parallel for shared(n,dx,al,pi,wi,s,ri,t) schedule (guided)
		for (integer i = 0; i < n; ++i) {
			dx[i] += al * pi[i] + wi * s[i];
			ri[i] = s[i] - wi * t[i];
		}
#ifdef  BICGSTAB_CRS_SMOOTHER
		if (bdebug) {
			isfinite_vec(n, dx, "dx");
			isfinite_vec(n, ri, "ri");
		}
#endif
		deltai = NormaV(ri, n);

#ifdef  BICGSTAB_CRS_SMOOTHER

		if (bdebug) {
			if (deltai != deltai) {
				printf("deltai is infinity");
				//system("pause");
				system("pause");
			}
		}
#endif
		// печать невязки на консоль
		if ((icount % 10) == 0)  std::cout << "iter  residual" << std::endl;

		std::cout << icount << " " << deltai << std::endl;

		//if ((icount % 100)== 0) system("pause");
		if (deltai / delta0_ < 0.8) iflag = 0;
		if (deltai < epsilon) iflag = 0; // конец вычисления
		else roim1 = roi;
	}

	//printf("internal: %lld %e \n", icount, deltai);
	//printf("%d", icount);
	//system("pause");

	// освобождение памяти
	delete[] ri; delete[] roc; delete[] s; delete[] t;
	delete[] vi; delete[] pi; delete[] dax;

#pragma omp parallel for
	for (integer i = 0; i < n; ++i) dX0[i] = dx[i];

	delete[] dx;


	system("pause");

} // Bi_CGStabCRS_smoother_naiv

// Метод Ван Дер Ворста Bi-CGStab
// работает для возможно несимметричных вещественных матриц.
// встроен предобуславливатель ILU(0).
// Метод является комбинацией методов BiCG и GMRES(1). 
// Экспериментально было выяснено что он не рекомендуется к использованию из-за
// серьёзных проблем, отсутствия сходимости. Это утверждение не относится собственно
// к алгоритму BiCGStab Х ван дер Ворста, а относится к данной конкретной реализации 
// на языке СИ данного алгоритма.
void Bi_CGStab_internal1(IMatrix *xO, equation3D* &sl, equation3D_bon* &slb,
			   integer maxelm, integer maxbound,
			   doublereal *dV, doublereal* &dX0, integer maxit, doublereal alpharelax,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls)
{

     printf("Bi_CGStab preconditioning by ILU(0)...\n"); 

	 integer i=0; // счётчик цикла for 
	 integer n = xO->n;// размер квадратной матрицы
	 // Разреженная матрица СЛАУ
	 // в CRS формате.
     doublereal *val;
     integer* col_ind, *row_ptr;

	 // преобразование из SIMPLESPARSE формата в CRS формат хранения.
	 //simplesparsetoCRS(M, val, col_ind, row_ptr, n);
	 equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax,true, b, lb, s_loc, ls);

	 // ILU предобуславливатель:
     doublereal *U_val, *L_val;
	 integer  *U_ind, *U_ptr, *L_ind, *L_ptr;

	 printf("Incoplete LU Decomposition begin...\n");
     convertIMatrixtoCSIR_ILU_ITL(xO, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);
	 // освобождение оперативной памяти
	 freeIMatrix(xO);
	 ILU0_Decomp_ITL(U_val, U_ind, U_ptr, L_val, L_ind, L_ptr, n);
	 printf("Incoplete LU Decomposition finish...\n");

	 doublereal *dx=new doublereal[n], *dr=new doublereal[n], *dr0=new doublereal[n],
		  *dax=new doublereal[n], *vcopy=new doublereal[n], *p=new doublereal[n],
          *s=new doublereal[n], *dap=new doublereal[n];

	 doublereal alpha, omega, dnz, scal1, beta;
	 doublereal epsilon=dterminatedTResudual;  // точность вычисления
	 integer icount;

	 // начальное приближение
     // X0 ==
     // под X0 понимается вектор поля температур к примеру.
     if (dX0==nullptr) {
	     dX0=new doublereal[n];
	     for (i=0; i<n; ++i) {
		     dx[i]=0.0;
		     dX0[i]=0.0;
	     }
     }
     else {
	     for (i=0; i<n; ++i) dx[i]=10.0;//dX0[i];
     }

	 MatrixCRSByVector(val,col_ind,row_ptr,dx, dax, n); // результат занесён в  dax
     for (i=0; i<n; ++i) dr[i]= dV[i] - dax[i];  // начальная невязка
     // dr=L^(-1)*(dV-A*dx);
     for (i=0; i<n; ++i) vcopy[i]=dr[i]; 
     inverseL_ITL(vcopy, L_val, L_ind, L_ptr, dr, n);
      // dr=U^(-1)*dr;
     for (i=0; i<n; ++i) vcopy[i]=dr[i];
     inverseU_ITL(vcopy, U_val, U_ind, U_ptr, dr, n);
	 // dr - невязка начального приближения.
	 for (i=0; i<n; ++i) {
		 //dr0[i]=dr[i];
		 dr0[i]=1.0;
		 p[i]=dr[i];
	 }

	 icount=0;
	 do {
		 icount++;
		 scal1=Scal(dr,dr0,n);

         MatrixCRSByVector(val,col_ind,row_ptr,p, vcopy, n); // результат занесён в  vcopy
         inverseL_ITL(vcopy, L_val, L_ind, L_ptr, dap, n);
         for (i=0; i<n; ++i) vcopy[i]=dap[i];
         inverseU_ITL(vcopy, U_val, U_ind, U_ptr, dap, n);

		 alpha=scal1/Scal(dap,dr,n);

		 for (i=0; i<n; ++i) {
			 s[i]=dr[i]-alpha*dap[i];
		 }

		 MatrixCRSByVector(val,col_ind,row_ptr,s, vcopy, n); // результат занесён в  vcopy
         inverseL_ITL(vcopy, L_val, L_ind, L_ptr, dap, n);
         for (i=0; i<n; ++i) vcopy[i]=dap[i];
         inverseU_ITL(vcopy, U_val, U_ind, U_ptr, dap, n);

         omega=Scal(dap,s,n)/Scal(dap,dap,n);

         for (i=0; i<n; ++i) {
			 dx[i]=dx[i]+alpha*p[i]+omega*s[i];
			 dr[i]=s[i]-omega*dap[i];
		 }

         dnz=NormaV(dr,n);
		 // печать невязки на консоль
         if ((icount % 10) == 0)  std::cout << "iter  residual" << std::endl;
		 std::cout << icount << " " << dnz << std::endl;
         

         beta=(Scal(dr,dr0,n)/scal1)*(alpha/omega);

         for (i=0; i<n; ++i) {
			 p[i]=dr[i]+beta*(p[i]-omega*dap[i]);
		 }

	 } while ((fabs(dnz)> epsilon) && (icount<maxit));



	 delete[] vcopy; delete[] dr0; delete[] dr;
	 delete[] dax; delete[] s; delete[] dap; delete[] p;
	 delete[] U_val; delete[] L_val;
	 delete[] U_ind; delete[] U_ptr; delete[] L_ind; delete[] L_ptr;
	 delete[] val; delete[] col_ind; delete[] row_ptr;
	 
	 for (i=0; i<n; ++i) dX0[i]=dx[i];
	 delete[] dx;

} // Bi_CGStab_internal1



// А.А.Фомин, Л.Н.Фомина 
// Ускорение полилинейного рекуррентного метода в подпространствах крылова.
// Вестник томского государственного университета. Математика и механика №2(14) 2011год.
// Алгоритм основан на прямом сочетании алгоритмов LR1 и Bi-CGStab P.
// LR1 - полилинейный метод предложенный еще в книге С. Патанкара: гибрид
// прямого метода прогонки (алгоритм Томаса) и метода Гаусса-Зейделя.
// Bi-CGStab P - алгоритм Ван Дер Ворста с предобуславливанием: гибрид Bi-CG и GMRES(1).
// начало написания, тестирования и использования в AliceFlow_v0_06 датируется 
// 24 октября 2011 года на основе предыдущих разработок.
// Метод Ван Дер Ворста Bi-CGStabCRS
// работает для возможно несимметричных вещественных матриц.
// Несимметричная матрица СЛАУ передаётся в CRS формате
// A (val, col_ind, row_ptr).
// Метод является комбинацией методов BiCG и GMRES(1).
//  begin 2: 5 мая 2012 года. Начало разработки системы балансировки.
// Система балансировки связана с регулированием количества итераций предобуславливателя,
// что должно приводить либо к усилению предобуславливания либо к ослаблению предобуславливания.
// 
// 25 марта 2013 года. Предобуславливание осуществляется через один проход методом прогонки по каждому
// координатному направлению. Для задач чистой теплопроводности одного вызова полинейного метода в двух местах
// на одной итерации BiCGStab  достаточно.
// Также обнаружено, что с этими настройками при отсутствии теплоотвода (всюду нейман или 3 род) метод
// Lr1sk - демонстрирует расходимость (переполнение).
void LR1sK(FLOW &f, equation3D* &sl, equation3D_bon* &slb,
	       doublereal *val, integer* col_ind, integer* row_ptr,
		   integer maxelm, integer maxbound, integer iVar,
		   doublereal *dV, doublereal* &dX0, integer maxit, bool &bprintmessage, bool bexporttecplot)
{
	bprintmessage = false;
	// если произошла ошибка (потеря сходимости) то будет произведён экспорт в программу tecplot 360
	// для диагностики полевых величин.
	//bexporttecplot=false; 
	
	//doublereal *val;
	//integer* col_ind;
	//integer* row_ptr;
	integer n=maxelm+maxbound;
	// преобразование из SIMPLESPARSE формата в CRS формат хранения.
	//equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound);

	bool bnorelax=true; // если bnorelax==false внутри предобуславливателя LR1sk применяется нижняя релаксация по компонентам скорости.
	if (iVar==PAM) bnorelax=true; // для поправки давления релаксация не применяется.

	// После того как программа была распараллелена эфективность полилинейного предобуславливателя упала.
	// Поэтому здесь вместо 3 циклов поставим 4 и для давления 5.
	integer imaxdubl=4; // 4 стартовое значение imaxdubl.
	if (iVar==PAM) imaxdubl=5; // 5
	bool bprintf=false; // если bprintf==false то значения невязок внутри LR1sk не выводятся.
	integer iflag=1, icount=0;
	doublereal delta0, deltai;
	doublereal bet, roi;
	doublereal roim1=1.0, al=1.0, wi=1.0;
	doublereal *ri = nullptr, *roc = nullptr, *s = nullptr, *t = nullptr;
	doublereal *vi = nullptr, *pi = nullptr, *dx = nullptr, *dax = nullptr;
	doublereal *y = nullptr, *z = nullptr; // результат предобуславливания
	doublereal epsilon=dterminatedTResudual;  // точность вычислений


	ri=new doublereal[n]; roc=new doublereal[n]; s=new doublereal[n]; t=new doublereal[n];
	vi=new doublereal[n]; pi=new doublereal[n]; dx=new doublereal[n]; dax=new doublereal[n];
	y=new doublereal[n]; z=new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания

    #pragma omp parallel for shared(s, t, vi, pi, y, z) schedule (guided)
	for (integer i=0; i<n; ++i) {
		s[i]=0.0;
		t[i]=0.0;
		vi[i]=0.0;
		pi[i]=0.0;
		// инициализатор массивов для предобуславливания
		y[i]=0.0;
		z[i]=0.0;
	}

    // начальное приближение
    // X0 ==
    // под X0 понимается вектор поля температур к примеру.
    if (dX0==nullptr) {
	   dX0=new doublereal[n];

	   #pragma omp parallel for shared(dx, dX0)  schedule (guided)
	   for (integer i=0; i<n; ++i) {
		   dx[i]=0.0;
		   dX0[i]=0.0;
	   }
    }
    else {
       #pragma omp parallel for shared(dx, dX0)  schedule (guided)
	   for (integer i=0; i<n; ++i) dx[i]=dX0[i];
    }

    MatrixCRSByVector(val,col_ind,row_ptr,dx,dax, n); // результат занесён в  dax

	#pragma omp parallel for shared(ri, dV, dax, roc) schedule (guided)
	for (integer i=0; i<n; ++i) {
		ri[i]=dV[i]-dax[i];
		roc[i]=ri[i];
	}
	delta0=NormaV(ri,n);

	// Если решение сразу хорошее то не считать:
	if (fabs(delta0)<dterminatedTResudual) iflag=0;

	
	if (iflag !=0) {
		// конечная невязка всегда на точность аппроксимации меньше начальной невязки.
	    epsilon*=delta0;
	    dterminatedTResudual=epsilon;
	}
	
	
	// Идея алгоритма балансировки.
	// Есть возможность усиливать предобуславливание
	// увеличивая количество итераций полинейного метода.
	// Предложение если на двух последовательных итерациях
	// падение невязки составило менее res[iter+1]/res[iter]<0.05 раз то 
	// слишком хорошее падение и нужно уменьшить imaxdubl, если
	// 0.05 <= res[iter+1]/res[iter] < 0.23 то нормальная скорость
	// сходимости imaxdubl трогать не надо. А если 
	// res[iter+1]/res[iter] >= 0.23 то нужно усилить предобуславливание
	// увеличив imaxdubl.

	doublereal* resiter=new doublereal[2];
	const integer iNOW=1;
	const integer iOLD=0;
	resiter[iOLD]=delta0; resiter[iNOW]=delta0;

	// Меняя эти параметры можно управлять эффективностью солвера.
	// Оптимальные значения этих параметров нужно установить из 
	// вычислительного эксперимента.
	//const doublereal LBAR=0.05;
	//const doublereal RBAR=0.23;

	integer iflag1=1;
	if (fabs(delta0)<1e-14) iflag1=0;

	if (bprintmessage) {
		switch (iVar) {
		case VELOCITY_X_COMPONENT: printf("VX	"); break;
		case VELOCITY_Y_COMPONENT: printf("VY "); break;
		case VELOCITY_Z_COMPONENT: printf("VZ "); break;
		case PAM: printf("PAM "); break;
		}
	}

	// magic_const это пороговая константа по достижению
	// которой в алгоритме возникают проблемы вещественной арифметики,
	// что говорит о том что вычисления требуется прекратить.
	const doublereal magic_const = 1.0e-20;
	bool bsignal_out = false;

	// Мы обязательно должны сделать три итерации. (не менее 10).
	while (((icount < 10) && (iflag1 != 0)) || (iflag != 0 && icount < maxit)) {

		icount++;
		if (icount > maxit) break;
		if (bsignal_out) break;


		roi=Scal(roc,ri,n);
		bet=(roi/roim1)*(al/wi);
		#pragma omp parallel for shared(pi, ri, vi, wi, bet)  schedule (guided)
		for (integer i=0; i<n; ++i) {
			doublereal rbufpi=ri[i]+(pi[i]-vi[i]*wi)*bet;
			pi[i]=rbufpi;
		}
	
		// Ky=pi
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(y) schedule (guided)
		for (integer i=0; i<n; ++i) y[i]=0.0; // Если начинать не с нуля то небудет сходимости для PAM !.
		// My=pi;
		solveLRn(y, pi, n, iVar, imaxdubl, bprintf, bnorelax, f.neighbors_for_the_internal_node, f.maxelm, f.slau, f.slau_bon, f.iN, f.id, f.iWE, f.iSN, f.iBT, f.alpha, f.maxbound);

		MatrixCRSByVector(val,col_ind,row_ptr,y,vi, n); // vi=A*y;

		if ((fabs(roi)<magic_const) && (fabs(Scal(roc, vi, n))<magic_const)) {
			al = 1.0;
			bsignal_out = true;
		}
		else if (fabs(roi)<magic_const) {
			al = 0.0;
			bsignal_out = true;
		}
		else {
			al = roi / Scal(roc, vi, n);
		}


		//al=roi/Scal(roc,vi,n); // см. выше.
		#pragma omp parallel for shared(s,ri, al, vi) schedule (guided)
        for (integer i=0; i<n; ++i) {
			s[i]=ri[i]-al*vi[i];
		}

		// Kz=s
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(z) schedule (guided)
		for (integer i=0; i<n; ++i) z[i]=0.0; // Если начинать не с нуля то небудет сходимости для PAM !.
        solveLRn(z, s, n, iVar, imaxdubl, bprintf, bnorelax, f.neighbors_for_the_internal_node, f.maxelm, f.slau, f.slau_bon, f.iN, f.id, f.iWE, f.iSN, f.iBT, f.alpha, f.maxbound);
		
        MatrixCRSByVector(val,col_ind,row_ptr,z,t, n);


		
		if ((fabs(Scal(t, s, n)) < magic_const) && (fabs(Scal(t, t, n))<magic_const)) {
			wi = 1.0;
			bsignal_out = true;
		}
		else if (fabs(Scal(t, s, n)) < magic_const) {
			wi = 0.0;
			bsignal_out = true;
		}
		else {
			wi = Scal(t, s, n) / Scal(t, t, n);
		}
		// не бдем ждать пока dx испортится.
		if (bsignal_out) break;

		//wi = Scal(t, s, n) / Scal(t, t, n); // см. выше.
		#pragma omp parallel for shared(dx, al, y, wi, z, ri, s, t) schedule (guided)
		for (integer i=0; i<n; ++i) {
			//dx[i]+=al*pi[i]+wi*s[i]; // так было без предобуславливателя
			dx[i]+=al*y[i]+wi*z[i]; // так стало с предобуславливателем
			ri[i]=s[i]-wi*t[i];
		}
		deltai=NormaV(ri,n);

		// Балансировка.
		resiter[iOLD]=resiter[iNOW]; resiter[iNOW]=deltai;
		doublereal dres=resiter[iNOW]/resiter[iOLD];
		// 1.0e-30
		if (fabs(dres - 1.0) < magic_const) {
			printf(" stagnation LR1SK... ");
			//system("pause");
			break;
		}
		/*
		if (fabs(dres)>RBAR) {
			imaxdubl++;
		}
		else if (fabs(dres)<LBAR) {
			imaxdubl=my_imax(1,imaxdubl-1);
		}
		*/

		if (bprintmessage) {
			// печать невязки на консоль
            if ((icount % 10) == 0)  {
				std::cout << "iter  residual imaxdubl" << std::endl;				
			}

			std::cout << icount << " " << deltai << " " << imaxdubl << std::endl;			
            
			//system("pause");
			system("pause");
		}

		if (deltai <epsilon) iflag=0; // конец вычисления
		else roim1=roi;

		// Иногда возникают матрицы с которыми метод не может справится.
		// Графически это означает что невязка опустилась до некоторого 
		// горизонтального асимптотического предела. Определим этот момент 
		// значением imaxdubl>=100; Остаётся только выйти из цикла и всё.
		// Ну ещё как вариант можно попробовать другой солвер.
		// По моему в данных случаях без AMG не обойтись. AMG -
		// алгебраический мультигридовый солвер, который борется с 
		// жёсткостью системы производя вычисления на последовательности вложенных сеток.
		// С другой стороны усиление влияния предобуславливателя 
		// должно уменьшать спектральный радиус матрицы, что 
		// в итоге ведёт к тому что сходимость должна достигаться 
		// на матрицах любого вида. Просто для очень плохообусловленных матриц
		// нужен сильный предобуславливатель. Вывод такой: комбинация решающего
		// алгоритма BiCGStab и надлежащего предобуславливателя должно справляться 
		// с любыми даже сколь угодно плохо обусловленными задачами.
		//if (imaxdubl>=100) {
			// значение константы 100 также должно быть подобрано из вычислительного эксперимента.

			//iflag=0; // конец вычисления
			//printf("calculation can not cope with the stiffness of the problem...\n");
			//printf("Please, press any key to continue calculation...\n");
			//bexporttecplot=true; // проблемы со сходимотью => экспорт картинки для анализа в программу tecplot.
			//system("pause");
			//system("pause");
		//}

		// system("pause"); // debug 1 iteration LR1sk
	}

	/* // отладочный код.
	// поиск контрольного объёма с максимальным значением невязки:
	doublereal maxerr=-1e30;
	integer ierr=-1;
	for (i=0; i<n; ++i) {
		if (fabs(ri[i])>maxerr) {
			ierr=i;
			maxerr=fabs(ri[i]);
		}
	}
	#if doubleintprecision == 1
		printf("node number max residual is %lld, value residal is equal %e\n", ierr, maxerr);
	#else
		printf("node number max residual is %d, value residal is equal %e\n", ierr, maxerr);
	#endif
	
	system("pause");
	*/

    // освобождение памяти
	delete[] ri; delete[] roc; delete[] s; delete[] t;
	delete[] vi; delete[] pi; delete[] dax;
	delete[] y; delete[] z;

	#pragma omp parallel for shared(dx, dX0) schedule (guided)
	for (integer i=0; i<n; ++i) dX0[i]=dx[i];

	delete[] dx; 
	delete[] resiter;

#if doubleintprecision == 1
	switch (iVar) {
	case VELOCITY_X_COMPONENT:	printf("VX %lld  ", icount); // контроль количества итераций.
		break;
	case VELOCITY_Y_COMPONENT: printf("VY %lld  ", icount); // контроль количества итераций.
		break;
	case VELOCITY_Z_COMPONENT: printf("VZ %lld  ", icount); // контроль количества итераций.
		break;
	case PAM: if (eqin.itemper == 1) {
		// потом еще решается уравнение теплопередачи.
		printf("PAM %lld  ", icount); // контроль количества итераций.
	}
			  else {
				  printf("PAM %lld  \n", icount); // контроль количества итераций.
			  }
			  break;
	}
#else
	switch (iVar) {
	case VELOCITY_X_COMPONENT:	printf("VX %d  ", icount); // контроль количества итераций.
		break;
	case VELOCITY_Y_COMPONENT: printf("VY %d  ", icount); // контроль количества итераций.
		break;
	case VELOCITY_Z_COMPONENT: printf("VZ %d  ", icount); // контроль количества итераций.
		break;
	case PAM: if (eqin.itemper == 1) {
		// потом еще решается уравнение теплопередачи.
		printf("PAM %d  ", icount); // контроль количества итераций.
	}
			  else {
				  printf("PAM %d  \n", icount); // контроль количества итераций.
			  }
			  break;
	}
#endif
	

} // LR1sK



// А.А.Фомин, Л.Н.Фомина 
// Ускорение полинейного рекуррентного метода в подпространствах крылова.
// Вестник томского государственного университета. Математика и механика №2(14) 2011год.
// Алгоритм основан на прямом сочетании алгоритмов LR1 и Bi-CGStab P.
// LR1 - полинейный метод предложенный еще в книге С. Патанкара: гибрид
// прямого метода прогонки (алгоритм Томаса) и метода Гаусса-Зейделя.
// Bi-CGStab P - алгоритм Ван Дер Ворста с предобуславливанием: гибрид Bi-CG и GMRES(1).
// начало написания, тестирования и использования в AliceFlow_v0_06 датируется 
// 24 октября 2011 года на основе предыдущих разработок.
// Метод Ван Дер Ворста Bi-CGStabCRS
// работает для возможно несимметричных вещественных матриц.
// Несимметричная матрица СЛАУ передаётся в CRS формате
// A (val, col_ind, row_ptr).
// Метод является комбинацией методов BiCG и GMRES(1).
//
// Для задач теплопроводности с учётом конвекции не работает
// по отдельности ни Ван-Дер-Ворст ни полилинейный метод.
// Возможно из-за того что СЛАУ ещё более плохо обусловлена 
// чем для компонент скорости.
// Выход. Попробовать применить LR1sK метод к задаче теплопроводности 
// с учётом конвекции.
// Начало разработки 23 ноября 2011 года.
//
void LR1sK_temp(TEMPER &tGlobal, equation3D* &sl, equation3D_bon* &slb,
	       doublereal *val, integer* col_ind, integer* row_ptr,
		   integer maxelm, integer maxbound, 
		   doublereal *dV, doublereal* &dX0, integer maxit, integer inumiter, bool bprintmessage, bool &bexporttecplot)
{

	

	// inumiter - номер глобальной итерации (например номер итерации в стационарном алгоритме SIMPLE).
	// параметр inumiter - введён для того чтобы использоваться при отладке, когда нужно посмотреть
	// алгоритм решения СЛАУ на глобальной итерации (алгоритма SIMPLE) с номером большим чем inumiter.

	if (0) {
		if (inumiter>82) {
			printf("debug LR1sk for temperature solver...\n");
		}
	}

	bexporttecplot=false; // экспорт в tecplot делается лишь в случае проблем со сходимостью.

	//doublereal *val;
	//integer* col_ind;
	//integer* row_ptr;
	integer n=maxelm+maxbound;
	// преобразование из SIMPLESPARSE формата в CRS формат хранения.
	//equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound);

	//bool bnorelax=true; // Для уравнения теплопроводности не используется релаксация.
	
	integer imaxdubl=1; // стартовое количество итераций полинейного метода 3
	
	bool bprintf=false; // если bprintf==false то значения невязок внутри LR1sk не выводятся.
	integer iflag=1, icount=0;
	doublereal delta0, deltai;
	doublereal bet, roi;
	doublereal roim1=1.0, al=1.0, wi=1.0;
	doublereal *ri, *roc, *s, *t, *vi, *pi, *dx, *dax;
	doublereal *y, *z; // результат предобуславливания
	doublereal epsilon=dterminatedTResudual;  // точность вычисления
	

	ri=new doublereal[n]; roc=new doublereal[n]; s=new doublereal[n]; t=new doublereal[n];
	vi=new doublereal[n]; pi=new doublereal[n]; dx=new doublereal[n]; dax=new doublereal[n];
	y=new doublereal[n]; z=new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания

	#pragma omp parallel for shared(s,t,vi,pi,y,z,dax)  schedule (guided)
	for (integer i=0; i<n; ++i) {
		s[i]=0.0;
		t[i]=0.0;
		vi[i]=0.0;
		pi[i]=0.0;
		// инициализатор массивов для предобуславливания
		y[i]=0.0;
		z[i]=0.0;
		// результат умножения матрицы на вектор.
		dax[i]=0.0;
	}

    // начальное приближение
    // X0 ==
    // под X0 понимается вектор поля температур к примеру.
    if (dX0==nullptr) {
	   dX0=new doublereal[n];
	   #pragma omp parallel for shared(dx, dX0)  schedule (guided)
	   for (integer i=0; i<n; ++i) {
		   dx[i]=0.0;
		   dX0[i]=0.0;
	   }
    }
    else {
	   #pragma omp parallel for shared(dx, dX0)  schedule (guided)
	   for (integer i=0; i<n; ++i) dx[i]=dX0[i];
    }

	MatrixCRSByVector(val,col_ind,row_ptr,dx,dax, n); // результат занесён в  dax

	#pragma omp parallel for shared(ri,dV,dax,roc) schedule (guided)
	for (integer i=0; i<n; ++i) {
		ri[i]=dV[i]-dax[i];
		roc[i]=ri[i];
	}
	delta0=NormaV(ri,n);
	
	//printf("debug %e\n",NormaV(dax,n)); // проверка на корректное составление СЛАУ
	//system("pause");
	// Если решение сразу хорошее то не считать:
	if (fabs(delta0)<dterminatedTResudual) iflag=0; 

	//printf("delta0=%e\n",delta0);
	//system("pause");

	
	/*if (iflag != 0) {
       // конечная невязка всегда на точность аппроксимации меньше начальной невязки.
	   epsilon*=delta0; 
	   dterminatedTResudual=epsilon;
	}
	*/

	// Идея алгоритма балансировки.
	// Есть возможность усиливать предобуславливание
	// увеличивая количество итераций полинейного метода.
	// Предложение если на двух последовательных итерациях
	// падение невязки составило менее res[iter+1]/res[iter]<0.05 раз то 
	// слишком хорошее падение и нужно уменьшить imaxdubl, если
	// 0.05 <= res[iter+1]/res[iter] < 0.23 то нормальная скорость
	// сходимости imaxdubl трогать не надо. А если 
	// res[iter+1]/res[iter] >= 0.23 то нужно усилить предобуславливание
	// увеличив imaxdubl.

	doublereal* resiter=new doublereal[2];
	const integer iNOW=1;
	const integer iOLD=0;
	resiter[iOLD]=delta0; resiter[iNOW]=delta0;

	// Меняя эти параметры можно управлять эффективностью солвера.
	// Оптимальные значения этих параметров нужно установить из 
	// вычислительного эксперимента.
	//const doublereal LBAR=0.05;
	//const doublereal RBAR=0.23;
	

	while ( iflag != 0 && icount < maxit) {

		icount++;

		roi=Scal(roc,ri,n);
		bet=(roi/roim1)*(al/wi);

		#pragma omp parallel for shared(pi,ri,vi,wi,bet)  schedule (guided)
		for (integer i=0; i<n; ++i) {
			doublereal pibuf=ri[i]+(pi[i]-vi[i]*wi)*bet;
			pi[i]=pibuf;
		}
	
		// Ky=pi
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(y)  schedule (guided)
		for (integer i=0; i<n; ++i) y[i]=0.0; // Если начинать не с нуля то не будет сходимости для PAM !.
		solveLRn_temp(tGlobal, y, pi, n, imaxdubl, bprintf);

		MatrixCRSByVector(val,col_ind,row_ptr,y,vi, n); // vi==A*y;
		al=roi/Scal(roc,vi,n);

		#pragma omp parallel for shared(s,ri,al,vi)  schedule (guided)
        for (integer i=0; i<n; ++i) {
			s[i]=ri[i]-al*vi[i];
		}

		// Kz=s
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(z)  schedule (guided)
		for (integer i=0; i<n; ++i) z[i]=0.0; // Если начинать не с нуля то не будет сходимости для PAM !.
        solveLRn_temp(tGlobal, z, s, n, imaxdubl, bprintf);
		
        MatrixCRSByVector(val,col_ind,row_ptr,z,t, n); // t==A*z;
		wi=Scal(t,s,n)/Scal(t,t,n);

		#pragma omp parallel for shared(dx, al, y, wi, z, ri, s, t)  schedule (guided)
		for (integer i=0; i<n; ++i) {
			//dx[i]+=al*pi[i]+wi*s[i]; // так было без предобуславливателя
			dx[i]+=al*y[i]+wi*z[i]; // так стало с предобуславливателем
			ri[i]=s[i]-wi*t[i];
		}
		deltai=NormaV(ri,n);


		// Балансировка.
		resiter[iOLD]=resiter[iNOW]; resiter[iNOW]=deltai;
		//doublereal dres=resiter[iNOW]/resiter[iOLD];
		/*
		if (fabs(dres)>RBAR) {
			imaxdubl++;
		}
		else if (fabs(dres)<LBAR) {
			imaxdubl=my_imax(1,imaxdubl-1);
		}
		*/

		// печать невязки на консоль
		if (bprintmessage) {
            if ((icount % 10) == 0)  {
				std::cout << "iter  residual imaxdubl" <<std::endl;				
			}
			std::cout << icount<< " " << deltai << " " << imaxdubl << std::endl;           
		}

		if (deltai <epsilon) iflag=0; // конец вычисления
		else roim1=roi;

		// Иногда возникают матрицы с которыми метод не может справится.
		// Графически это означает что невязка опустилась до некоторого 
		// горизонтального асимптотического предела. Определим этот момент 
		// значением imaxdubl>=100; Остаётся только выйти из цикла и всё.
		// Ну ещё как вариант можно попробовать другой солвер.
		// По моему в данных случаях без AMG не обойтись. AMG -
		// алгебраический мультигридовый солвер, который борется с 
		// жёсткостью системы производя вычисления на последовательности вложенных сеток.
		// С другой стороны усиление влияния предобуславливателя 
		// должно уменьшать спектральный радиус матрицы, что 
		// в итоге ведёт к тому что сходимость должна достигаться 
		// на матрицах любого вида. Просто для очень плохообусловленных матриц
		// нужен сильный предобуславливатель. Вывод такой: комбинация решающего
		// алгоритма BiCGStab и надлежащего предобуславливателя должно справляться 
		// с любыми даже сколь угодно плохо обусловленными задачами.
		//if (imaxdubl>=100) {
			// значение константы 100 также должно быть подобрано из вычислительного эксперимента.

			//iflag=0; // конец вычисления
			//printf("calculation can not cope with the stiffness of the problem...\n");
			//printf("Please, press any key to continue calculation...\n");
			//bexporttecplot=true; // проблемы со сходимотью => экспорт картинки для анализа в программу tecplot.
			//system("pause");
			//system("pause");
		//}

		//system("pause"); // debug 1 iteration LR1sk
	}

    // освобождение памяти
	delete[] ri; delete[] roc; delete[] s; delete[] t;
	delete[] vi; delete[] pi; delete[] dax;
	delete[] y; delete[] z;

	#pragma omp parallel for shared(dX0, dx)  schedule (guided)
	for (integer  i=0; i<n; ++i) dX0[i]=dx[i];

	delete[] dx; 
	delete[] resiter;

#if doubleintprecision == 1
	printf("TEMP %lld  \n", icount); // контроль количества итераций.
#else
	printf("TEMP %d  \n", icount); // контроль количества итераций.
#endif
	

} // LR1sK_temp


// этот метод показывает значительно более лучшую сходимость, чем простой BiCGStabCRS,
// а также он гораздо лучше (и повидимому правильней) чем Bi_CGStab_internal1.
// дата написания Bi_CGStab_internal2: 31.03.2013. 
void Bi_CGStab_internal2(IMatrix *xO, equation3D* &sl, equation3D_bon* &slb,
			   integer maxelm, integer maxbound,
			   doublereal *dV, doublereal* &dX0, integer maxit, doublereal alpharelax,
			   bool bprintmessage, BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls)
{

	// inumiter - номер глобальной итерации (например номер итерации в стационарном алгоритме SIMPLE).
	// параметр inumiter - введён для того чтобы использоваться при отладке, когда нужно посмотреть
	// алгоритм решения СЛАУ на глобальной итерации (алгоритма SIMPLE) с номером большим чем inumiter.

	//bool bexporttecplot=false; // экспорт в tecplot делается лишь в случае проблем со сходимостью.

	// размерность квадратной матрицы.
	integer n=maxelm+maxbound;
	

	// Разреженная матрица СЛАУ
	 // в CRS формате.
     doublereal *val;
     integer *col_ind, *row_ptr;

	 // преобразование из SIMPLESPARSE формата в CRS формат хранения.
	 //simplesparsetoCRS(M, val, col_ind, row_ptr, n);
	 equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax,true, b, lb, s_loc, ls);

	 // ILU предобуславливатель:
     doublereal *U_val, *L_val;
	 integer  *U_ind, *U_ptr, *L_ind, *L_ptr;

	 bool bSaad=false; // внимание не используйте true т.к. оно содержит ошибки которые надо искать.

	 if (bSaad) {
		   // ILU(0) разложение из книги Й. Саада
           printf("Incoplete LU Decomposition I.Saad begin...\n");
		   freeIMatrix(xO);
		   doublereal *luval;
		   integer *uptr;
		   integer icode=-1;
           ilu0_Saad(n, val, col_ind, row_ptr, luval, uptr, icode); // ILU(0) разложение
		   if (icode==(-1)) {
			   IMatrix xO1;
               initIMatrix(&xO1, n); // инициализация

               convertCRStoIMatrix(n, luval, col_ind, row_ptr, uptr, &xO1);
			   delete luval; 
			   delete uptr;
			   convertIMatrixtoCSIR_ILU_ITL(&xO1, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);
               freeIMatrix(&xO1);
			   // сортировка по возрастанию
	           for (integer i = 0; i < n; ++i) {
                   QuickSortCSIR(U_ind, U_val, U_ptr[i], U_ptr[i+1]-1);
                   QuickSortCSIR(L_ind, L_val, L_ptr[i], L_ptr[i+1]-1);
				   L_val[L_ptr[i]]=1.0; // единица на главной диагонали.
	           }

		   }
		   else {
#if doubleintprecision == 1
			   printf("Error!!! zero  diagonal elem in %lld string matrix.\n", icode);
#else
			   printf("Error!!! zero  diagonal elem in %d string matrix.\n", icode);
#endif
			   
			   //system("pause");
			   system("pause");
			   exit(0); // выход из программы.
		   }

		   printf("Incoplete LU Decomposition I.Saad finish...\n");
	 }
	 else {

	      printf("Incoplete LU Decomposition begin...\n");
          convertIMatrixtoCSIR_ILU_ITL(xO, U_val, U_ind, U_ptr, L_val, L_ind, L_ptr);
	      // освобождение оперативной памяти
	      freeIMatrix(xO);
	      ILU0_Decomp_ITL(U_val, U_ind, U_ptr, L_val, L_ind, L_ptr, n);
	      printf("Incoplete LU Decomposition finish...\n");
	 }

	//bool bnorelax=true; // Для уравнения теплопроводности не используется релаксация.
	
	
	//bool bprintf=false; // если bprintf==false то значения невязок внутри LR1sk не выводятся.
	integer iflag=1, icount=0;
	doublereal delta0, deltai;
	doublereal bet, roi;
	doublereal roim1=1.0, al=1.0, wi=1.0;
	doublereal *ri, *roc, *s, *t, *vi, *pi, *dx, *dax, *vcopy;
	doublereal *y, *z; // результат предобуславливания
	doublereal epsilon=dterminatedTResudual;  // точность вычисления
	

	ri=new doublereal[n]; roc=new doublereal[n]; s=new doublereal[n]; t=new doublereal[n];
	vi=new doublereal[n]; pi=new doublereal[n]; dx=new doublereal[n]; dax=new doublereal[n];
	vcopy=new doublereal[n]; // временный вектор.
	y=new doublereal[n]; z=new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания

	#pragma omp parallel for shared(s,t,vi,pi,y,z,dax,vcopy)  schedule (guided)
	for (integer i=0; i<n; ++i) {
		s[i]=0.0;
		t[i]=0.0;
		vi[i]=0.0;
		pi[i]=0.0;
		// инициализатор массивов для предобуславливания
		y[i]=0.0;
		z[i]=0.0;
		// результат умножения матрицы на вектор.
		dax[i]=0.0;
		vcopy[i]=0.0; 
	}

    // начальное приближение
    // X0 ==
    // под X0 понимается вектор поля температур к примеру.
    if (dX0==nullptr) {
	   dX0=new doublereal[n];
	   #pragma omp parallel for shared(dx, dX0)  schedule (guided)
	   for (integer i=0; i<n; ++i) {
		   dx[i]=0.0;
		   dX0[i]=0.0;
	   }
    }
    else {
	   #pragma omp parallel for shared(dx, dX0)  schedule (guided)
	   for (integer i=0; i<n; ++i) dx[i]=dX0[i];
    }

	MatrixCRSByVector(val,col_ind,row_ptr,dx,dax, n); // результат занесён в  dax

	#pragma omp parallel for shared(ri,dV,dax,roc) schedule (guided)
	for (integer i=0; i<n; ++i) {
		ri[i]=dV[i]-dax[i];
		//roc[i]=ri[i];
		roc[i]=1.0;
	}
	delta0=NormaV(ri,n);
	
	//printf("debug %e\n",NormaV(dax,n)); // проверка на корректное составление СЛАУ
	//system("pause");
	// Если решение сразу хорошее то не считать:
	if (fabs(delta0)<dterminatedTResudual) iflag=0; 

	//printf("delta0=%e\n",delta0);
	//system("pause");

	
	/*if (iflag != 0) {
       // конечная невязка всегда на точность аппроксимации меньше начальной невязки.
	   epsilon*=delta0; 
	   dterminatedTResudual=epsilon;
	}
	*/
	

	while ( iflag != 0 && icount < maxit) {

		icount++;

		roi=Scal(roc,ri,n);
		bet=(roi/roim1)*(al/wi);

		#pragma omp parallel for shared(pi,ri,vi,wi,bet) schedule (guided)
		for (integer i=0; i<n; ++i) {
			doublereal pibuf=ri[i]+(pi[i]-vi[i]*wi)*bet;
			pi[i]=pibuf;
		}
	
		// Ky=pi
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(y) schedule (guided)
		for (integer i=0; i<n; ++i) y[i]=0.0; // Если начинать не с нуля то небудет сходимости для PAM !.
		//solveLRn_temp(tGlobal, y, pi, n, imaxdubl, bprintf);
		inverseL_ITL(pi, L_val, L_ind, L_ptr, y, n); // y=inverse(L)*pi;
		#pragma omp parallel for shared(y,vcopy)  schedule (guided)
        for (integer i=0; i<n; ++i) vcopy[i]=y[i]; // vcopy=inverse(L)*pi;
        inverseU_ITL(vcopy, U_val, U_ind, U_ptr, y, n); // y=inverse(U)*vcopy;

		MatrixCRSByVector(val,col_ind,row_ptr,y,vi, n); // vi==A*y;
		al=roi/Scal(roc,vi,n);

		#pragma omp parallel for shared(s,ri,al,vi)  schedule (guided)
        for (integer i=0; i<n; ++i) {
			s[i]=ri[i]-al*vi[i];
		}

		// Kz=s
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(z)  schedule (guided)
		for (integer i=0; i<n; ++i) z[i]=0.0; // Если начинать не с нуля то небудет сходимости для PAM !.
        //solveLRn_temp(tGlobal, z, s, n, imaxdubl, bprintf);
		inverseL_ITL(s, L_val, L_ind, L_ptr, z, n); // z=inverse(L)*s;
		#pragma omp parallel for shared(z,vcopy)  schedule (guided)
        for (integer i=0; i<n; ++i) vcopy[i]=z[i]; // vcopy=inverse(L)*s;
        inverseU_ITL(vcopy, U_val, U_ind, U_ptr, z, n); // z=inverse(U)*vcopy;
		
        MatrixCRSByVector(val,col_ind,row_ptr,z,t, n); // t==A*z;
		wi=Scal(t,s,n)/Scal(t,t,n);

		#pragma omp parallel for shared(dx, al, y, wi, z, ri, s, t)  schedule (guided)
		for (integer i=0; i<n; ++i) {
			//dx[i]+=al*pi[i]+wi*s[i]; // так было без предобуславливателя
			dx[i]+=al*y[i]+wi*z[i]; // так стало с предобуславливателем
			ri[i]=s[i]-wi*t[i];
		}
		deltai=NormaV(ri,n);
		
		

		// печать невязки на консоль
		if (bprintmessage) {
            if ((icount % 10) == 0)  {
				std::cout << "iter  residual" << std::endl;
			}

			std::cout << icount << " " << deltai << std::endl;
            
		}

		if (deltai <epsilon) iflag=0; // конец вычисления
		else roim1=roi;

	}

    // освобождение памяти
	delete[] ri; delete[] roc; delete[] s; delete[] t;
	delete[] vi; delete[] pi; delete[] dax;
	delete[] y; delete[] z;
	delete[] vcopy; 
	delete[] U_val; delete[] L_val;
	delete[] U_ind; delete[] U_ptr; delete[] L_ind; delete[] L_ptr;
	delete[] val; delete[] col_ind; delete[] row_ptr;

	#pragma omp parallel for shared(dX0, dx)  schedule (guided)
	for (integer i=0; i<n; ++i) dX0[i]=dx[i];

	delete[] dx; 

} // Bi_CGStab_internal2

// Возвращает позицию в векторе переменных.
integer iposfunc(integer iVar) {
	switch (iVar) {
	case VELOCITY_X_COMPONENT: return 0; break;
	case VELOCITY_Y_COMPONENT: return 1; break;
	case VELOCITY_Z_COMPONENT: return 2; break;
	case PAM: return 3; break;
	case TEMP: return 4; break;
	default: return -1; break;
	}
} // iposfunc

/* 
13 апреля 2013 года вся память выделяемая и уничтожаемая внутри
BICGSTAB_internal3 вынесена наружу, с целью препятствовать частым выделениям и уничтожениям памяти.
Это должно положительным образом сказаться на скорости работы BiCGStab.
// Объявление этой структуры данные перенесено выше по коду в модуль amg1r5.cpp.
*//*
typedef struct TQuickMemVorst {

	 //doublereal *rthdsd; // правая часть системы уравнений

	 // Исходная матрица в формате CRS.
	 doublereal *val;
     integer *col_ind;
	 integer *row_ptr;

	 // Рабочие вектора.
	 doublereal *ri, *roc, *s, *t, *vi, *pi, *dx, *dax;
	 doublereal *y, *z; // результат предобуславливания

	 doublereal *a; // CRS
	 integer *ja;
	 integer *ia;
	 // ILU предобуславливатель:
	 // Будем сразу хранить в MSR формате, который заведём с запасом по памяти.
	 doublereal *alurc;
	 integer *jlurc;
	 integer *jurc;
	 doublereal *vec;
	 doublereal *alu;
	 integer *jlu;
	 integer *ju;
	 // Копия для распараллеливания lusol_2
	 // Копия матрицы разложения нужна для распараллеливания обратного хода по U матрице.
	 doublereal *x1; // Копия искомого решения.
	 doublereal *alu1; // Копия матрицы 
	 integer *jlu1; // ILU2 разложения.
	 integer *ju1; // это тоже относится к копии матрицы.
	 integer *iw; // для ILU(0)
	 // для ILU(lfil)
	 integer* levs; 
	 doublereal* w; 
	 integer* jw;

	 integer iwk; // размерность памяти под матрицу предобуславливания.

	 //doublereal *trthdsd; // правая часть системы уравнений

	 // Исходная матрица в формате CRS.
	 doublereal *tval;
     integer *tcol_ind;
	 integer *trow_ptr;

	  // Рабочие вектора.
	 doublereal *tri, *troc, *ts, *tt, *tvi, *tpi, *tdx, *tdax;
	 doublereal *ty, *tz;

	 bool ballocCRScfd;
	 bool bsignalfreeCRScfd; // сигнал к освобождению памяти из под a,ja,ia. если true
	 doublereal *ta; // CRS
	 integer *tja;
	 integer *tia;
	  // ILU предобуславливатель:
	 // Будем сразу хранить в MSR формате, который заведём с запасом по памяти.
	 doublereal *talu;
	 integer *tjlu;
	 integer *tju;
	 integer *tiw; // для ILU(0)
	 // для ILU(lfil)
	 integer* tlevs; // для ILU(lfil)
	 doublereal* tw; 
	 integer* tjw;

	 integer tiwk; // размерность памяти под матрицу предобуславливания.

	 bool ballocCRSt;
	 bool bsignalfreeCRSt; // сигнал к освобождению памяти из под a,ja,ia. если true

	 // Идея запомнить сколько итераций было сделано для компонент скорости 
	 // чтобы использовать такое-же количество итераций для температуры в задачах с естественой конвекцией.
	 integer icount_vel; 
} QuickMemVorst;
*/


// 18.01.2019
// метод получает будущий вектор цвета не инициализированный, но память выделена.
// 1 - первая область, 2 вторая область, 3 - граница областей. Для number_part=2.
// 1 певая часть, 2 часть, 3 - сшивка, 4 третья часть, 5 - четвертая часть, 6 - сшивка, 7 сшивка. Для number_part=4.
// связи хранятся в матрице в координатном формате crs.
// length_separator = 1, 2, 3, ... клетки.
// number_part - число частей, степень двойки 1, 2, 4, 8, 16, 32,...
// Нумерация начинается с нуля.
void nested_desection_crs(integer*& col_ind, integer*& row_ptr, integer n, integer*& color, integer& dist_max
/*integer number_part, integer length_separator*/);

// Переупорядочивание в соответствии с алгоритмом вложенных сечений.
// Переупорядочивание векторов правой части и решения, а также матрицы СЛАУ.
void reorder_direct(integer*& col_ind, integer*& row_ptr, doublereal* &val, integer n,
	integer*& new_order, doublereal*& dV, doublereal*& dX0) 
{
	doublereal* tmp = new doublereal[n];

	for (integer i_1 = 0; i_1 < n; ++i_1) {
		tmp[new_order[i_1]] = dV[i_1];
	}
	for (integer i_1 = 0; i_1 < n; ++i_1) {
		dV[i_1] = tmp[i_1];
	}
	for (integer i_1 = 0; i_1 < n; ++i_1) {
		tmp[new_order[i_1]] = dX0[i_1];
	}
	for (integer i_1 = 0; i_1 < n; ++i_1) {
		dX0[i_1] = tmp[i_1];
	}

	delete[] tmp;


	integer nnz = row_ptr[n];
	integer *itmp = new integer[nnz];
	integer *ia= new integer[nnz];
	for (integer i_1 = 0; i_1 < nnz; ++i_1) {
		ia[i_1] = -100;// Заглушка вызывающая ошибку памяти. Её сразу можно отловить.
	}
	for (integer i_1 = 0; i_1 < n; ++i_1) {
		for (integer j_1 = row_ptr[i_1]; j_1 <= row_ptr[i_1 + 1] - 1; ++j_1) {
			//printf("%lld %lld\n", col_ind[row_ptr[i_1]],i_1);
			//system("pause");
			// Первое значение всегда диагональ.
			ia[j_1] = i_1;
		}
	}
	for (integer i_1 = 0; i_1 < nnz; ++i_1) {
		itmp[i_1]=new_order[col_ind[i_1]];
	}
	for (integer i_1 = 0; i_1 < nnz; ++i_1) {
		col_ind[i_1] = itmp[i_1];
	}
	for (integer i_1 = 0; i_1 < nnz; ++i_1) {
		itmp[i_1] = new_order[ia[i_1]];
	}
	for (integer i_1 = 0; i_1 < nnz; ++i_1) {
		ia[i_1] = itmp[i_1];
	}

	// val остаётся без изменений, изменились только индексы.

	delete[] itmp;

	// CountingSort по строкам.
	Ak2 Amat;
	Amat.i= new integer_mix_precision[nnz];
	Amat.j= new integer_mix_precision[nnz];
	Amat.aij = new real_mix_precision[nnz];
	Amat.abs_aij = new real_mix_precision[nnz];

	// Прямое копирование.
	for (integer i_1 = 0; i_1 < nnz; ++i_1) {
		Amat.i[i_1] = static_cast<integer_mix_precision>(ia[i_1]);
		Amat.j[i_1] = static_cast<integer_mix_precision>(col_ind[i_1]);
		Amat.aij[i_1] = static_cast<real_mix_precision>(val[i_1]);
		Amat.abs_aij[i_1] = static_cast<real_mix_precision>(abs(val[i_1]));
	}

	// Сортировка Тима Петерсона по строкам.
	// 18.01.2020
	timSort_amg(Amat, 0, nnz - 1);
	//HeapSort(Amat, 0, nnz - 1);
	//qs(Amat, 0, nnz - 1);

	// Обратное копирование.
	for (integer i_1 = 0; i_1 < nnz; ++i_1) {
		ia[i_1] = Amat.i[i_1];
		col_ind[i_1]= Amat.j[i_1];
		val[i_1] = Amat.aij[i_1];
	}

	delete[] Amat.i;
	delete[] Amat.j;
	delete[] Amat.aij;
	delete[] Amat.abs_aij;

	// coo -> crs
	bool* flag = new bool[n];
	for (integer i_1 = 0; i_1 < n; ++i_1) {
		flag[i_1] = false;
	}
	for (integer ii = 0; ii < nnz; ++ii) {
		if (flag[ia[ii]] == false) {
			row_ptr[ia[ii]] = ii;
			flag[ia[ii]] = true;
		}
	}
	row_ptr[n] = nnz;

	for (integer i_1 = 0; i_1 < n; ++i_1) {
		for (integer j_1 = row_ptr[i_1]; j_1 <= row_ptr[i_1 + 1] - 1; ++j_1) {
			if (col_ind[row_ptr[i_1]] != i_1) {
				printf("first index no diagonal\n");
				printf("%lld %lld\n", col_ind[row_ptr[i_1]], i_1);
				system("PAUSE");
				// Первое значение всегда диагональ.
			}
		}
	}

	delete[] flag;	
	delete[] ia;

} // reorder_direct

//integer icount_first_PAM_iteration_SIMPLE_algorithm_global = 0;

// этот метод показывает значительно более лучшую сходимость, чем простой BiCGStabCRS,
// а также он гораздо лучше (и повидимому правильней) чем Bi_CGStab_internal1.
// Bi_CGStab_internal3 использует предобуславливание из библиотеки Ю.Саада.
// дата написания Bi_CGStab_internal3: 31.03.2013. 
// раскурочил 9 августа 2015.
// 26 сентября 2016 Теперь метод пригоден и для АЛИС сетки.
// 18.01.2020 начало распараллеливания на два потока с помощью алгоритма вложенных сечений nested desection.
void Bi_CGStab_internal3(equation3D* &sl, equation3D_bon* &slb,
			   int maxelm, int maxbound,
			   doublereal *dV, doublereal* &dX0, integer maxit, doublereal alpharelax,
			   bool bprintmessage, integer iVar, QuickMemVorst& m,
	           int* &ifrontregulationgl, int* &ibackregulationgl,
	           BLOCK* &b, int &lb, SOURCE* &s, int &ls, integer inumber_iteration_SIMPLE,
               integer* &color, integer dist_max, bool breordering_for_parallel,
			   WALL* &w, int &lw)
{
	

	//printf("1. alpharelax=%e \n", alpharelax);

	if (breordering_for_parallel) {
		BiCGStab_internal3_incomming_now = true;
	}
	else {
		BiCGStab_internal3_incomming_now = false;
	}

	// Если память под решение еще не выделена.
	if (dX0 == nullptr) {
		dX0 = new doublereal[maxelm+maxbound];
		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
#pragma omp parallel for shared(m, dX0) schedule (guided)
			for (integer i_37 = 0; i_37 < maxelm + maxbound; i_37++) {
				dX0[i_37] = 0.0;
			}

		}
		if (iVar == TEMP) {
#pragma omp parallel for shared(m, dX0) schedule (guided)
			for (integer i_37 = 0; i_37 < maxelm + maxbound; i_37++) {
				dX0[i_37] = 0.0;
			}
		}

	}

	// inumiter - номер глобальной итерации (например номер итерации в стационарном алгоритме SIMPLE).
	// параметр inumiter - введён для того чтобы использоваться при отладке, когда нужно посмотреть
	// алгоритм решения СЛАУ на глобальной итерации (алгоритма SIMPLE) с номером большим чем inumiter.

	//bool bexporttecplot=false; // экспорт в tecplot делается лишь в случае проблем со сходимостью.
	bool brc=false;
	bool bpam_gsp=false; // предобуславливание с помощью нескольких итераций метода Гаусса-Зейделя.

	// размерность квадратной матрицы.
	integer n=maxelm+maxbound;
	

	// Разреженная матрица СЛАУ
	 // в CRS формате.
	bool bNORMIROVKA = true;
	for (integer k = 0; k < maxelm; ++k) {
		if (fabs(sl[k].ap - 1.0) > 1.0e-30) bNORMIROVKA = false;

		if (dV[k] != dV[k]) {
			printf("NAN or INF in iP=%lld in dV rthdsd internal\n", k);
			switch (iVar) {
			case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
			case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
			case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
			case NUSHA: printf("NU equation problem.\n"); break;
			case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY equation problem.\n"); break;
			case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA equation problem.\n"); break;
			case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS equation problem.\n"); break;
			case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS equation problem.\n"); break;
			case RE_THETA_LANGTRY_MENTER: printf("TURBULENT_RE_THETA_LANGTRY_MENTER equation problem\n"); break;
			case GAMMA_LANGTRY_MENTER: printf("TURBULENT_GAMMA_LANGTRY_MENTER equation problem\n"); break;
			case PAM: printf("PAM equation problem.\n"); break;
			}
			system("pause");
			exit(1);
		}
	}
	if (bNORMIROVKA) {
		//printf("bNORMIROVKA Ok\n");
		//system("pause");
	}

#pragma omp parallel
	for (integer k = maxelm; k < maxelm+maxbound; ++k) {
		if (dV[k] != dV[k]) {
			printf("NAN or INF in iP=%lld in dV rthdsd boundary\n", k);
			switch (iVar) {
			case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
			case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
			case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
			case NUSHA: printf("NU equation problem.\n"); break;
			case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY equation problem.\n"); break;
			case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA equation problem.\n"); break;
			case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS equation problem.\n"); break;
			case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS equation problem.\n"); break;
			case RE_THETA_LANGTRY_MENTER: printf("TURBULENT_RE_THETA_LANGTRY_MENTER equation problem\n"); break;
			case GAMMA_LANGTRY_MENTER: printf("TURBULENT_GAMMA_LANGTRY_MENTER equation problem\n"); break;
			case PAM: printf("PAM equation problem.\n"); break;
			}
			system("pause");
			exit(1);
		}
	}

	

	 if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM)||(iVar==NUSHA)||
		 (iVar== TURBULENT_KINETIK_ENERGY)||(iVar== TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA)||
		 (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		 (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		 if (ibackregulationgl!=nullptr) {
			 printf(" if (ibackregulationgl!=nullptr)\n"); // сюда он не должен заходить.
			 system("PAUSE");
			 // nested desection версия алгоритма.
			 integer ierr=equation3DtoCRSnd(sl, slb, m.val, m.col_ind, m.row_ptr, maxelm, maxbound, alpharelax,!m.ballocCRScfd, ifrontregulationgl, ibackregulationgl,b,lb,s,ls);
			 if (ierr > 0) {
				 switch (iVar) {
				 case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				 case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				 case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				 case NUSHA: printf("NU equation problem.\n"); break;
				 case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY equation problem.\n"); break;
				 case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA equation problem.\n"); break;
				 case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS equation problem.\n"); break;
				 case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS equation problem.\n"); break;
				 case RE_THETA_LANGTRY_MENTER: printf("TURBULENT_RE_THETA_LANGTRY_MENTER equation problem\n"); break;
				 case GAMMA_LANGTRY_MENTER: printf("TURBULENT_GAMMA_LANGTRY_MENTER equation problem\n"); break;
				 case PAM: printf("PAM equation problem.\n"); break;
				 }
			 }
		 }
		 else {
		     integer ierr=equation3DtoCRS(sl, slb, m.val, m.col_ind, m.row_ptr, maxelm, maxbound, alpharelax,!m.ballocCRScfd, b, lb, s, ls);
			 if (ierr > 0) {
				 switch (iVar) {
				 case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				 case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				 case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				 case NUSHA: printf("NU equation problem.\n"); break;
				 case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY equation problem.\n"); break;
				 case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA equation problem.\n"); break;
				 case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS equation problem.\n"); break;
				 case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS equation problem.\n"); break;
				 case RE_THETA_LANGTRY_MENTER: printf("TURBULENT_RE_THETA_LANGTRY_MENTER equation problem\n"); break;
				 case GAMMA_LANGTRY_MENTER: printf("TURBULENT_GAMMA_LANGTRY_MENTER equation problem\n"); break;
				 case PAM: printf("PAM equation problem.\n"); break;
				 }
			 }
		 }
	 }
	 if (iVar==TEMP) {
		 integer ierr=equation3DtoCRS(sl, slb, m.tval, m.tcol_ind, m.trow_ptr, maxelm, maxbound, alpharelax,!m.ballocCRSt, b, lb, s, ls);
		
		 doublereal alpharelax = 1.0;

		 bool bnonlinear18 = false;

		 if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE) {

			 // 02.12.2021
			 // Ахтунг, ахтунг !!!
			 // При нестационарном моделировании ни в коем случае нельзя присваивать
			 // alpharelax = 0.99999; это искажает нестационарный член уравнения теплопередачи.
			 // Проверено на amg1r5 и BiCGStab+ILU(k).

			  // Это не специальная нелинейная версия кода amg1r5 CAMG.
			 for (integer k = 0; k < lw; ++k) {
				 if ((w[k].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
					 (w[k].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY)) {
					// alpharelax = 0.99999; // Для того чтобы СЛАУ сходилась.
					 // 0.9999 - недостаточное значение, температуры не те получаются.
					 bnonlinear18 = true;
				 }
			 }

			 if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
				 (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
				 (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {
				 //alpharelax = 0.99999;
				 bnonlinear18 = true;
			 }

			 //std::cout << "unsteady thermal calculation" << std::endl;
			 //system("pause");
		 }
		 else {

			 // Это не специальная нелинейная версия кода amg1r5 CAMG.
			 for (integer k = 0; k < lw; ++k) {
				 if ((w[k].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
					 (w[k].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY)) {
					 alpharelax = 0.99999; // Для того чтобы СЛАУ сходилась.
					 // 0.9999 - недостаточное значение, температуры не те получаются.
					 bnonlinear18 = true;
				 }
			 }

			 if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
				 (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
				 (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {
				 alpharelax = 0.99999;
				 bnonlinear18 = true;
			 }
		 }
		 //if (adiabatic_vs_heat_transfer_coeff == ADIABATIC_WALL_BC) {
		 //printf("ADIABATIC WALL BC"); system("pause");
		 //}
		 
		 //printf("2. alpharelax=%e \n", alpharelax);

		 if (bnonlinear18) {
			 // Разреженная матрица СЛАУ
			 // в CRS формате.

			 {
				 integer nsize = maxelm + maxbound;
				// integer nnzsize = m.trow_ptr[maxelm+maxbound];

				 integer ierr = 0;
				 
				 for (integer i_5 = 0; i_5 < nsize; i_5++) {

					 for (integer i_6 = m.trow_ptr[i_5]; i_6 <= m.trow_ptr[i_5 + 1] - 1; i_6++) {
						 if (m.tcol_ind[i_6] == i_5) {
							 if (m.trow_ptr[i_5 + 1] > m.trow_ptr[i_5] + 1) {
								 dV[i_5] += (1.0 - alpharelax)*dX0[i_5]* m.tval[i_6] / alpharelax;
								 m.tval[i_6] /= alpharelax;
							 }
						 }
					 }
					 /*if (i_5 <= 1) {
						 printf("row_ptr75[%lld]=%lld \n",i_5,ia[i_5]);
						 system("pause");
					 }
					 */
				 }

				 if (ierr > 0) {
					 printf("Temperature equation problem.\n");
				 }
			 }

		 }
		 
		 if (ierr > 0) {
			 printf("Temperature equation problem.\n");
		 }
	 }
	 
	 // преобразование из SIMPLESPARSE формата в CRS формат хранения.
	 //simplesparsetoCRS(M, val, col_ind, row_ptr, n);
	

	 const integer ILU0=0;
	 const integer ILU_lfil=1;

	 const integer itype_ilu=ILU_lfil;//ILU_lfil;

	
	 if (b_on_adaptive_local_refinement_mesh) {
		 //integer* 
		 color = new integer[n];
		 for (integer i_1 = 0; i_1 < n; ++i_1) color[i_1] = 0; // initialization
		//integer 
		 dist_max = -1;
	 }
	 integer* new_order = new integer[n];	 
	// integer length_separator = my_amg_manager.lfil + 1;//только единица.
	 doublereal* valcopy=nullptr;
	 integer* col_indcopy = nullptr;
	 integer* row_ptrcopy = nullptr;

	 if (breordering_for_parallel) {
		 if ((my_amg_manager.lfil < 3) && (number_cores() > 1)) {
			 // nested_desection
			 
			 if (b_on_adaptive_local_refinement_mesh) {
				 if (iVar == TEMP) {
					 nested_desection_crs(m.tcol_ind, m.trow_ptr, n, color, dist_max);
				 }
				 if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
					 (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
					 (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
					 (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
					 nested_desection_crs(m.col_ind, m.row_ptr, n, color, dist_max);
				 }
			 }


			 // new_order[исходный индекс]=новый индекс.
			 //перенумеруется вектор правой части, вектор решения и стобец col_ind и row_ind. Снова запаковать в crs.
			 if (number_cores() == 2) {
				 // две части.
				 //PARDATA nd;
				 nd.ncore = 2; // два ядра.
				 // По умолчанию все разбиения неактивны.
				 nd.b0.active = true;
				 nd.b00.active = false;
				 nd.b01.active = false;
				 nd.b000.active = false;
				 nd.b001.active = false;
				 nd.b010.active = false;
				 nd.b011.active = false;
				 nd.b0.ileft_start = 0;
				 integer icount_new_order = 0;
				 integer isep = static_cast<integer>(0.5 * (dist_max));//2
				 if (!b_on_adaptive_local_refinement_mesh) {
					 isep = 2;
				 }
				 integer bar = 0;// 0; 1; 20;
				 if (my_amg_manager.lfil > 0) bar = 1;
				 integer il = 0, ic = 0, ir = n;
				 if (b_on_adaptive_local_refinement_mesh) {
					 while (abs(ir - il) > 1.4 * ic) {
						 icount_new_order = 0;
						 il = 0; ir = 0; ic = 0;// инициализация.
						 //if (length_separator == 1) 
						 {
							 // Разделитель одна клетка.
							 for (integer i_1 = 0; i_1 < n; ++i_1) {
								 if (color[i_1] < isep - bar) {
									 // первая часть.
									 new_order[i_1] = icount_new_order;
									 icount_new_order++;
									 il++;
								 }
							 }
							 nd.b0.ileft_finish = icount_new_order - 1;
							 nd.b0.iright_start = icount_new_order;
							 for (integer i_1 = 0; i_1 < n; ++i_1) {
								 if (color[i_1] > isep + bar) {
									 // вторая часть.
									 new_order[i_1] = icount_new_order;
									 icount_new_order++;
									 ir++;
								 }
							 }
							 nd.b0.iright_finish = icount_new_order - 1;
							 nd.b0.iseparate_start = icount_new_order;
							 for (integer i_1 = 0; i_1 < n; ++i_1) {
								 if ((color[i_1] >= isep - bar) && (color[i_1] <= isep + bar)) {
									 // разделитель.
									 new_order[i_1] = icount_new_order;
									 icount_new_order++;
									 ic++;
								 }
							 }
							 nd.b0.iseparate_finish = icount_new_order - 1;
							 printf("ileft=%lld center=%lld right=%lld\n", il, ic, ir);
							 if (ir > il) {
								 isep++;
							 }
							 else if (ir < il) {
								 isep--;
							 }
							 //system("pause");
						 }


					 }
				 }
				 else {
					 bar = 0; isep = 2; // только так.
					 icount_new_order = 0;
					 il = 0; ir = 0; ic = 0;// инициализация.
											//if (length_separator == 1) 
					 {
						 // Разделитель одна клетка.
						 for (integer i_1 = 0; i_1 < n; ++i_1) {
							 if (color[i_1] < isep - bar) {
								 // первая часть.
								 new_order[i_1] = icount_new_order;
								 icount_new_order++;
								 il++;
							 }
						 }
						 nd.b0.ileft_finish = icount_new_order - 1;
						 nd.b0.iright_start = icount_new_order;
						 for (integer i_1 = 0; i_1 < n; ++i_1) {
							 if (color[i_1] > isep + bar) {
								 // вторая часть.
								 new_order[i_1] = icount_new_order;
								 icount_new_order++;
								 ir++;
							 }
						 }
						 nd.b0.iright_finish = icount_new_order - 1;
						 nd.b0.iseparate_start = icount_new_order;
						 for (integer i_1 = 0; i_1 < n; ++i_1) {
							 if ((color[i_1] >= isep - bar) && (color[i_1] <= isep + bar)) {
								 // разделитель.
								 new_order[i_1] = icount_new_order;
								 icount_new_order++;
								 ic++;
							 }
						 }
						 nd.b0.iseparate_finish = icount_new_order - 1;
						 printf("ileft=%lld center=%lld right=%lld\n", il, ic, ir);
					 }
				 }

			 }

			 if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) ||
				 (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {

				 delete[] valcopy;
				 delete[] col_indcopy;
				 delete[] row_ptrcopy;

				 valcopy = new doublereal[m.row_ptr[n]];
				 col_indcopy = new integer[m.row_ptr[n]];
				 for (integer i_1 = 0; i_1 < m.row_ptr[n]; ++i_1) {
					 valcopy[i_1] = m.val[i_1];
					 col_indcopy[i_1] = m.col_ind[i_1];
				 }
				 row_ptrcopy = new integer[n + 1];
				 for (integer i_1 = 0; i_1 < n + 1; ++i_1) {
					 row_ptrcopy[i_1] = m.row_ptr[i_1];
				 }

				 reorder_direct(m.col_ind, m.row_ptr, m.val, n, new_order, dV, dX0);
			 }
			 if (iVar == TEMP) {

				 delete[] valcopy;
				 delete[] col_indcopy;
				 delete[] row_ptrcopy;

				 valcopy = new doublereal[m.trow_ptr[n]];
				 col_indcopy = new integer[m.trow_ptr[n]];
				 for (integer i_1 = 0; i_1 < m.trow_ptr[n]; ++i_1) {
					 valcopy[i_1] = m.tval[i_1];
					 col_indcopy[i_1] = m.tcol_ind[i_1];
				 }
				 row_ptrcopy = new integer[n + 1];
				 for (integer i_1 = 0; i_1 < n + 1; ++i_1) {
					 row_ptrcopy[i_1] = m.trow_ptr[i_1];
				 }

				 reorder_direct(m.tcol_ind, m.trow_ptr, m.tval, n, new_order, dV, dX0);
			 }
		 }
	 }
	 
     // Исходная матрица.
	 if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM)||(iVar==NUSHA)||
		 (iVar== TURBULENT_KINETIK_ENERGY)||(iVar== TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		 (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		 (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
	    if (!m.ballocCRScfd) {
	        // m.a=new doublereal[7*n+2]; // CRS
	        // m.ja=new integer[7*n+2];
			// 26 сентября 2016.
			m.a = new doublereal[m.row_ptr[n] + 2 * maxbound + 2];
			m.ja = new integer[m.row_ptr[n] + 2 * maxbound + 2];
	        m.ia=new integer[n+2];
		    // Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((m.a==nullptr)||(m.ja==nullptr)||(m.ia==nullptr)) {
			     // недостаточно памяти на данном оборудовании.
			     printf("Problem: not enough memory on your equipment...\n");
				 printf("Please any key to exit...\n");
				 exit(1);
			}
	    }
	 }
	 if (iVar==TEMP) {
		 if (!m.ballocCRSt) {
	       // m.ta=new doublereal[7*n+2]; // CRS
	       // m.tja=new integer[7*n+2];
			 // 26 сентября 2016.
			m.ta = new doublereal[m.trow_ptr[n] + 2 * maxbound + 2];
			m.tja = new integer[m.trow_ptr[n] + 2 * maxbound + 2];
	        m.tia=new integer[n+2];
		    // Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((m.ta==nullptr)||(m.tja==nullptr)||(m.tia==nullptr)) {
			     // недостаточно памяти на данном оборудовании.
			     printf("Problem: not enough memory on your equipment...\n");
				 printf("Please any key to exit...\n");
				 exit(1);
			}
	    }
	 }
	 // Эти структуры данных используются в библиотеке Ю.Саада SPARSKIT2.
	 // Основной момент который следует уяснить это:
	 // мы используем индексацию массивов в СИ начиная с нуля и до size не включая size;
	 // В библиотеке Sparskit2 нумерация элементов массива начинается с единицы и до size включая size.
	 // Главное что следует понять не следует пытаться переписать SPARSKIT2 так чтобы она соответствовала нумерации с нуля.
	 // оставим нумерацию с единицы иначе мы можем запутаться в логике отлаженного кода SPARSKIT2.
	 // Но мы в проекте AliceFlowv0_07 работаем с нумерации начинающейся с нуля. 
	 // Код Sparskit2 работает только для предобуславливания. Иными словами схема следующая:
	 // На входе матрица в CRS формате с нумерацией с нуля. Преобразуем её в матрицу в которой элементы нумеруются с единицы,
	 // для этого очевидно нужно к значениям элементов  col_ind и row_ptr прибавить единицу. Индексация же в матрице (преобразованной пусть начинается с нуля)
	 // для этого мы внутри кода SPARSKIT2 декременируем ссылку на массив тогда элемент с номером ноль станет первым тоесть будет иметь номер ноль. Потом в конце использования
	 // кода SPARSKIT2 мы увеливаем ссылки на массивы на единицу (откат назад). Получая на вход данную преобразованную матрицу CRS формата, мы делаем с помощью неё матрицу предобуславливателя
	 // в MSR формате, как обычно используя код SPARSKIT2 без каких либо изменений из вне (он получен с помощью f2c.exe). Используя матрицу предобуславливания в MSR формате
	 // мы пихаем её в функию lusol_ из SPARSKIT2 и получаем необходимый вектор x: (LU)x=y; по вектору y и матрице LU в формате MSR. Код в lusol_ содержит преобразование декремента
	 // указателей x, y что позволяет использовать обычные векторы x , y в которых нумерация начинается с нуля. В общем x,y менять не надо, они такие же как обычно в AliceFlowv0_07.
	 // В lusol_ указатели сначала декримируются --a; а в конце инкремируются ++a; Так что можно использовать код Sparskit2 без изменений !!!


	 
    if (bprintmessage) {
	    printf("Incoplete LU Decomposition begin...\n");
    }
    
	
	integer ierr=0;
	if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) 
	{
#pragma omp parallel
	   for (integer i=0; i<m.row_ptr[n]; ++i) {
		   m.a[i]=m.val[i];
		   m.ja[i]=m.col_ind[i]+1;
	   }
#pragma omp parallel
	   for (integer i=0; i<n+1; ++i) {
		  m.ia[i]=m.row_ptr[i]+1;
	   }
	}
	if (iVar==TEMP) 
	{
#pragma omp parallel
		for (integer i=0; i<m.trow_ptr[n]; ++i) {
		    m.ta[i]=m.tval[i];
		    m.tja[i]=m.tcol_ind[i]+1;
	    }
#pragma omp parallel
	    for (integer i=0; i<n+1; ++i) {
			m.tia[i]=m.trow_ptr[i]+1;
	    }
	}

	 if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
		 (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		 (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		 (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		 if (!m.ballocCRScfd) {
			 m.ri=new doublereal[n]; m.roc=new doublereal[n]; m.s=new doublereal[n]; m.t=new doublereal[n]; m.vec=new doublereal[n];
	         m.vi=new doublereal[n]; m.pi=new doublereal[n]; m.dx=new doublereal[n]; m.dax=new doublereal[n];
	         m.y=new doublereal[n]; m.z=new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
			 if ((m.ri==nullptr)||(m.roc==nullptr)||(m.s==nullptr)||(m.t==nullptr)||(m.vi==nullptr)||(m.pi==nullptr)||(m.dx==nullptr)||(m.dax==nullptr)||(m.y==nullptr)||(m.z==nullptr)) {
				  // недостаточно памяти на данном оборудовании.
			     printf("Problem: not enough memory on your equipment...\n");
				 printf("Please any key to exit...\n");
				 exit(1);
			 }
		 }
	 }
	 if (iVar==TEMP) {
         if (!m.ballocCRSt) {
             m.tri=new doublereal[n]; m.troc=new doublereal[n]; m.ts=new doublereal[n]; m.tt=new doublereal[n];
	         m.tvi=new doublereal[n]; m.tpi=new doublereal[n]; m.tdx=new doublereal[n]; m.tdax=new doublereal[n];
	         m.ty=new doublereal[n]; m.tz=new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
			 if ((m.tri==nullptr)||(m.troc==nullptr)||(m.ts==nullptr)||(m.tt==nullptr)||(m.tvi==nullptr)||(m.tpi==nullptr)||(m.tdx==nullptr)||(m.tdax==nullptr)||(m.ty==nullptr)||(m.tz==nullptr)) {
				  // недостаточно памяти на данном оборудовании.
			     printf("Problem: not enough memory on your equipment...\n");
				 printf("Please any key to exit...\n");
				 exit(1);
			 }
         }
	 }

	 if (0&&(my_amg_manager.lfil == 0)) {

		 // Почемуто не работает ilu0 на больших ренольдсах.

		 //if (itype_ilu==ILU0)
		 {

			 if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) ||
				 (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {

				 if (!m.ballocCRScfd) {
					 //m.alu=new doublereal[7*n+2]; // +2 запас по памяти.
					 //m.jlu=new integer[7*n+2];
					 // 26 сентября 2016.
					 m.alu = new doublereal[m.row_ptr[n] + 2 * maxbound + 2];
					 m.jlu = new integer[m.row_ptr[n] + 2 * maxbound + 2];

					 m.ju = new integer[n + 2];
					 if (ibackregulationgl != nullptr) {
						 // m.alu1=new doublereal[7*n+2]; // +2 запас по памяти.
						 // m.jlu1=new integer[7*n+2];
						 // m.ju1=new integer[n+2];
						 m.x1 = new doublereal[n + 2];
					 }
					 //m.alurc=new doublereal[7*n+2]; // +2 запас по памяти.
					 //m.jlurc=new integer[7*n+2];
					 // 26 сентября 2016.
					 m.alurc = new doublereal[m.row_ptr[n] + 2 * maxbound + 2];
					 m.jlurc = new integer[m.row_ptr[n] + 2 * maxbound + 2];

					 m.jurc = new integer[n + 2];
					 m.iw = new integer[n + 2]; // рабочий массив.
					 m.ballocCRScfd = true; // память выделена.

					 if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.ju == nullptr) || (m.iw == nullptr)) {
						 // недостаточно памяти на данном оборудовании.
						 printf("Problem: not enough memory on your equipment...\n");
						 printf("Please any key to exit...\n");
						 exit(1);
					 }
					 if (breordering_for_parallel) {
						 if ((m.alu1 == nullptr) || (m.jlu1 == nullptr) || (m.ju1 == nullptr) || (m.x1 == nullptr)) {
							 // недостаточно памяти на данном оборудовании.
							 printf("Problem: not enough memory on your equipment...\n");
							 printf("Please any key to exit...\n");
							 exit(1);
						 }
					 }

				 }
			 }
			 if (iVar == TEMP) {
				 if (!m.ballocCRSt) {
					 //m.talu=new doublereal[7*n+2]; // +2 запас по памяти.
					 //m.tjlu=new integer[7*n+2];
					 // 26 сентября 2016.
					 m.talu = new doublereal[m.trow_ptr[n] + 2 * maxbound + 2];
					 m.tjlu = new integer[m.trow_ptr[n] + 2 * maxbound + 2];

					 m.tju = new integer[n + 2];
					 m.tiw = new integer[n + 2]; // рабочий массив.
					 m.ballocCRSt = true; // память выделена.

					 if ((m.talu == nullptr) || (m.tjlu == nullptr) || (m.tju == nullptr) || (m.tiw == nullptr)) {
						 // недостаточно памяти на данном оборудовании.
						 printf("Problem: not enough memory on your equipment...\n");
						 printf("Please any key to exit...\n");
						 exit(1);
					 }

				 }
			 }


			 if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) ||
				 (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
				 ilu0_(n, m.a, m.ja, m.ia, m.alu, m.jlu, m.ju, m.iw, ierr);
				 /* if (ibackregulationgl!=nullptr) {
					  for (integer i87=0; i87<7*n+2; i87++) {
						  m.alu1[i87]= m.alu[i87];
						  m.jlu1[i87]=m.jlu[i87];
					  }
					  for (integer i87=0; i87<n+2; i87++) {
						  m.ju1[i87]=m.ju[i87];
					  }
				  }*/
			 }
			 if (iVar == TEMP) {
				 ilu0_(n, m.ta, m.tja, m.tia, m.talu, m.tjlu, m.tju, m.tiw, ierr);
			 }

			 if (ierr > 0) {
#if doubleintprecision == 1
				 printf("%lld string in matrix is zero diagonal element...\n", ierr - 1);
#else
				 printf("%d string in matrix is zero diagonal element...\n", ierr - 1);
#endif

				 //system("pause");
				 system("pause");
				 exit(0);
			 }
		 }
	 }
	 else {
		 if (itype_ilu == ILU_lfil) {

			 //bool btemp_quick = m.ballocCRSt;

			 integer lfil = 2; // 2 уровня (0, 1, 2)

			 lfil = my_amg_manager.lfil;

			 if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) ||
				 (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
				 if (!m.ballocCRScfd) {

					 // инициализация.
					 m.alu = nullptr;
					 m.jlu = nullptr;
					 m.ju = nullptr;
					 m.alu1 = nullptr;
					 m.jlu1 = nullptr;
					 m.ju1 = nullptr;
					 m.x1 = nullptr;
					 m.alurc = nullptr;
					 m.jlurc = nullptr;
					 m.jurc = nullptr;
					 m.levs = nullptr;
					 m.w = nullptr;
					 m.jw = nullptr;
					 m.w_dubl = nullptr;
					 m.jw_dubl = nullptr;

					 //m.iwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
					 // 26 сентября 2016.
					 if (lfil <= 2) {
						 m.iwk = (lfil + 1) * (m.row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
					 }
					 else if (lfil == 3) {
						 m.iwk = (2 * lfil + 1) * (m.row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
					 }
					 else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
						 m.iwk = (3 * lfil + 1) * (m.row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
					 }
					 else if (lfil >= 6) {
						 m.iwk = (5 * lfil + 1) * (m.row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
					 }

					 m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
					 m.jlu = new integer[m.iwk + 2];
					 m.ju = new integer[2 * n + 2];
					 if (ibackregulationgl != nullptr) {
						 //m.alu1=new doublereal[m.iwk+2]; // +2 запас по памяти.
						 //m.jlu1=new integer[m.iwk+2];
						 //m.ju1=new integer[n+2];
						 m.x1 = new doublereal[n + 2];
					 }
					 m.alurc = new doublereal[m.iwk + 2]; // +2 запас по памяти.
					 m.jlurc = new integer[m.iwk + 2];
					 m.jurc = new integer[2 * n + 2];
					 m.levs = new integer[m.iwk + 2]; // уровень.
					 m.w = new doublereal[static_cast<integer>(2 * n) + 2]; // +2 запас по памяти.
					 m.w_dubl = new doublereal[static_cast<integer>(2 * n) + 2]; // +2 запас по памяти.
					 if (lfil <= 2) {
						 m.jw = new integer[3 * n + 2]; // +2 запас по памяти.				
						 m.jw_dubl = new integer[3 * n + 2]; // +2 запас по памяти.
					 }
					 else if (lfil == 3) {
						 m.jw = new integer[3 * lfil * n + 2]; // +2 запас по памяти.				
						 m.jw_dubl = new integer[3 * lfil * n + 2]; // +2 запас по памяти.
					 }
					 else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
						 m.jw = new integer[4 * lfil * n + 2]; // +2 запас по памяти.				
						 m.jw_dubl = new integer[4 * lfil * n + 2]; // +2 запас по памяти.
					 }
					 else if (lfil >= 6) {
						 m.jw = new integer[6 * lfil * n + 2]; // +2 запас по памяти.				
						 m.jw_dubl = new integer[6 * lfil * n + 2]; // +2 запас по памяти.
					 }
					 m.ballocCRScfd = true; // память выделена.

					 if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.levs == nullptr) || (m.ju == nullptr) || (m.w == nullptr) || (m.jw == nullptr) || (m.w_dubl == nullptr) || (m.jw_dubl == nullptr)) {
						 // недостаточно памяти на данном оборудовании.
						 printf("Problem: not enough memory on your equipment...\n");
						 printf("Please any key to exit...\n");
						 exit(1);
					 }

				 }
			 }
			 if (iVar == TEMP) {
				 if (!m.ballocCRSt) {

					 // инициализация.
					 m.talu = nullptr;
					 m.tjlu = nullptr;
					 m.tju = nullptr;
					 m.tlevs = nullptr;
					 m.tw = nullptr;
					 m.tjw = nullptr;


					 //m.tiwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
					 // 26 сентября 2016.
					 if (lfil <= 2) {
						 m.tiwk = (lfil + 1) * (m.trow_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
					 }
					 else if (lfil == 3) {
						 m.tiwk = (2 * lfil + 1) * (m.trow_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
					 }
					 else if (/*(lfil >= 4) && */(lfil <= 5)) {
						 m.tiwk = (3 * lfil + 1) * (m.trow_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
					 }
					 else if (lfil >= 6) {
						 m.tiwk = (5 * lfil + 1) * (m.trow_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
					 }
					 m.talu = new doublereal[m.tiwk + 2]; // +2 запас по памяти.
					 m.tjlu = new integer[m.tiwk + 2];
					 m.tju = new integer[2 * n + 2];
					 m.tlevs = new integer[m.tiwk + 2]; // уровень.
					 m.tw = new doublereal[2 * n + 2]; // +2 запас по памяти.
					 if (lfil <= 2) {
						 m.tjw = new integer[3 * n + 2]; // +2 запас по памяти.
					 }
					 else if (lfil == 3) {
						 m.tjw = new integer[3 * lfil * n + 2]; // +2 запас по памяти.
					 }
					 else if (/*(lfil >= 4) && */(lfil <= 5)) {
						 m.tjw = new integer[4 * lfil * n + 2]; // +2 запас по памяти.
					 }
					 else if (lfil >= 6) {
						 m.tjw = new integer[6 * lfil * n + 2]; // +2 запас по памяти.
					 }
					 m.ballocCRSt = true; // память выделена.

					 if ((m.talu == nullptr) || (m.tjlu == nullptr) || (m.tlevs == nullptr) || (m.tju == nullptr) || (m.tw == nullptr) || (m.tjw == nullptr)) {
						 // недостаточно памяти на данном оборудовании.
						 printf("Problem: not enough memory on your equipment...\n");
						 printf("Please any key to exit...\n");
						 exit(1);
					 }
				 }
			 }


			 if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				 (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) ||
				 (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
				 // iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
				 iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);

				 if ((ierr == -2) || (ierr == -3)) {

					 integer ipassage = 1; // 4 января 2016.
					 do {
						 printf("\nPlease WAIT... ... ...\n");

						 // задаче не хватило памяти, значит нужно перевыделить !
						 if (m.alu != nullptr) delete[] m.alu;
						 if (m.jlu != nullptr) delete[] m.jlu;
						 /* if (ibackregulationgl!=nullptr) {
							   if (m.alu1!=nullptr) delete m.alu1;
							   if (m.jlu1!=nullptr) delete m.jlu1;
						  }*/
						 if (m.alurc != nullptr) delete[] m.alurc;
						 if (m.jlurc != nullptr) delete[] m.jlurc;
						 if (m.levs != nullptr) delete[] m.levs;

						 // инициализация !
						 m.alu = nullptr;
						 m.jlu = nullptr;
						 /* if (ibackregulationgl!=nullptr) {
							  m.alu1=nullptr;
							  m.jlu1=nullptr;
						  }*/
						 m.levs = nullptr;

						 // инициализация !
						 m.alurc = nullptr;
						 m.jlurc = nullptr;

						 //m.iwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
						 // 26 сентября 2016.
						 if (lfil <= 2) {
							 m.iwk = (lfil + 1) * (m.row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage) * n);
						 }
						 else if (lfil == 3) {
							 m.iwk = (2 * lfil + 1) * (m.row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage) * n);
						 }
						 else if (/*(lfil >= 4) && */(lfil <= 5)) {
							 m.iwk = (3 * lfil + 1) * (m.row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage) * n);
						 }
						 else if (lfil >= 6) {
							 m.iwk = (6 * lfil + 1) * (m.row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage) * n);
						 }

						 m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
						 m.jlu = new integer[m.iwk + 2];
						 /* (ibackregulationgl!=nullptr) {
							  m.alu1=new doublereal[m.iwk+2]; // +2 запас по памяти.
							  m.jlu1=new integer[m.iwk+2];
						 }*/
						 m.levs = new integer[m.iwk + 2]; // уровень.

						 if ((m.alu != nullptr) && (m.jlu != nullptr) && (m.levs != nullptr)) {
							 // iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
							 iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);
							 /*
							if (ibackregulationgl!=nullptr) {
								for (integer i87=0; i87<m.iwk+2; i87++) {
									m.alu1[i87]= m.alu[i87];
									m.jlu1[i87]=m.jlu[i87];
								}
								for (integer i87=0; i87<n+2; i87++) {
									m.ju1[i87]=m.ju[i87];
								}
						 }*/

						 }
						 else {
							 // недостаточно памяти на данном оборудовании.
							 ipassage = 4;
							 printf("Problem: not enough memory on your equipment...\n");
							 printf("Please any key to exit...\n");
							 exit(1);

						 }

						 ipassage++;
					 } while ((ierr != 0) && (ipassage < 4));

					 if (ipassage == 4) {
						 printf("Error memory alloc !!!\n");
						 printf("failed to obtain an expansion for the 4 approaches...\n");
						 printf("Please any key to exit...\n");
						 exit(1);
					 }
				 }
				 else {
					 /*
					  if (ibackregulationgl!=nullptr) {
						  for (integer i87=0; i87<m.iwk+2; i87++) {
							  m.alu1[i87]= m.alu[i87];
							  m.jlu1[i87]=m.jlu[i87];
						  }
						  for (integer i87=0; i87<n+2; i87++) {
							  m.ju1[i87]=m.ju[i87];
						  }
					  }*/
				 }

			 }
			 else if (iVar == TEMP) {

				 /*
				 if (0&&bglobal_unsteady_temperature_determinant) {
					 // модификация 20_10_2016.
					 // При нестационарном моделировании в твёрдом теле будем строить предобуславливатель лишь единожды на первом шаге.
					 //if (!btemp_quick) {
					 // iluk_call speed_up
					 // 10%=2 9%
					 // 15%=3 9%  19.25
					 // 20%=4 9%
					 // 30%=6 3%
					 if (rand()%20<3) {
						 iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
					 }
				 }
				 else {
				 */
				 iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);

				 //}

				 if ((ierr == -2) || (ierr == -3)) {

					 integer ipassage = 1;
					 do {
						 printf("\nPlease WAIT... ... ...\n");

						 // задаче не хватило памяти, значит нужно перевыделить !
						 if (m.talu != nullptr) delete[] m.talu;
						 if (m.tjlu != nullptr) delete[] m.tjlu;
						 if (m.tlevs != nullptr) delete[] m.tlevs;

						 // инициализация !
						 m.talu = nullptr;
						 m.tjlu = nullptr;
						 m.tlevs = nullptr;

						 //m.tiwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
						 // 26 сентября 2016.
						 if (lfil <= 2) {
							 m.tiwk = (lfil + 1) * (m.trow_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage) * n);
						 }
						 else if (lfil == 3) {
							 m.tiwk = (2 * lfil + 1) * (m.trow_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage) * n);
						 }
						 else if (/*(lfil>=4)&&*/(lfil <= 5)) {
							 m.tiwk = (3 * lfil + 1) * (m.trow_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage) * n);
						 }
						 else if (lfil >= 6) {
							 m.tiwk = (6 * lfil + 1) * (m.trow_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage) * n);
						 }

						 m.talu = new doublereal[m.tiwk + 2]; // +2 запас по памяти.
						 m.tjlu = new integer[m.tiwk + 2];
						 m.tlevs = new integer[m.tiwk + 2]; // уровень.

						 if ((m.talu != nullptr) && (m.tjlu != nullptr) && (m.tlevs != nullptr)) {
							 iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
						 }
						 else {
							 // недостаточно памяти на данном оборудовании.
							 ipassage = 4;
							 printf("Problem: not enough memory on your equipment...\n");
							 printf("Please any key to exit...\n");
							 exit(1);

						 }

						 ipassage++;
					 } while ((ierr != 0) && (ipassage < 4));

					 if (ipassage == 4) {
						 printf("Error memory alloc !!!\n");
						 printf("failed to obtain an expansion for the 4 approaches...\n");
						 printf("Please any key to exit...\n");
						 exit(1);
					 }
				 }
			 }



			 if (ierr != 0) {
#if doubleintprecision == 1
				 printf("error memory in iluk ierr=%lld\n", ierr);
#else
				 printf("error memory in iluk ierr=%d\n", ierr);
#endif

				 //system("pause");
				 system("pause");
				 exit(0);
			 }
		 }
	 }
	if (bprintmessage) {
	    printf("Incoplete LU Decomposition finish...\n");
	}
	

	//bool bnorelax=true; // Для уравнения теплопроводности не используется релаксация.
	
	
	//bool bprintf=false; // если bprintf==false то значения невязок внутри LR1sk не выводятся.
	integer iflag=1, icount=0;
	doublereal delta0=1.0e30, deltai=1.0e30;
	doublereal bet=0.0, roi=0.0;
	doublereal roim1=1.0, al=1.0, wi=1.0;
	//doublereal *ri, *roc, *s, *t, *vi, *pi, *dx, *dax;
	//doublereal *y, *z; // результат предобуславливания
	doublereal epsilon=dterminatedTResudual;  // точность вычисления
	if (iVar==TEMP) {
		epsilon*=1.0e-4; // 1.0e-4
	}


	

	#pragma omp parallel for shared(m,iVar)  schedule (guided)
	for (integer i=0; i<n; ++i) {
		if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		   m.s[i]=0.0;
		   m.t[i]=0.0;
		   m.vi[i]=0.0;
		   m.pi[i]=0.0;
		   // инициализатор массивов для предобуславливания
		   m.y[i]=0.0;
		   m.z[i]=0.0;
		   // результат умножения матрицы на вектор.
		   m.dax[i]=0.0;
		}
		if (iVar==TEMP) {
			m.ts[i]=0.0;
		    m.tt[i]=0.0;
		    m.tvi[i]=0.0;
		    m.tpi[i]=0.0;
		    // инициализатор массивов для предобуславливания
		    m.ty[i]=0.0;
		    m.tz[i]=0.0;
		    // результат умножения матрицы на вектор.
		    m.tdax[i]=0.0;
		}
	}

    // начальное приближение
    // X0 ==
    // под X0 понимается вектор поля температур к примеру.
    if (dX0==nullptr) {
	   dX0=new doublereal[n];
	   if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
		   (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		   (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		   (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
#pragma omp parallel for shared(m, dX0) schedule (guided)
		   for (integer i_37 = 0; i_37<n; i_37++) {
			   m.dx[i_37] = 0.0;
			   dX0[i_37] = 0.0;
		   }
		   
	   }
	   if (iVar==TEMP) {
		    #pragma omp parallel for shared(m, dX0) schedule (guided)
		   for (integer i_37 = 0; i_37<n; i_37++) {
		       m.tdx[i_37]=0.0;
		       dX0[i_37]=0.0;
	        }
	   }
	  
    }
    else {
      if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
		  (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		  (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		  (iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		  if (ibackregulationgl!=nullptr) {
               #pragma omp parallel for shared(m, dX0, ifrontregulationgl) schedule (guided)
	           for (integer i=0; i<n; ++i) {
				   // по новой нумерации с индексом i получает индекс старой нумерации iP
				   //iP=ifrontregulation[i];
		           m.dx[i]=dX0[ifrontregulationgl[i]];
	           }
		   }
		   else {
		    #pragma omp parallel for shared(m, dX0)  schedule (guided)
	        for (integer i=0; i<n; ++i) m.dx[i]=dX0[i];
		   }
	  }
	  if (iVar==TEMP) {
		   #pragma omp parallel for shared(m, dX0)  schedule (guided)
	       for (integer i=0; i<n; ++i) m.tdx[i]=dX0[i];
	  }
	  
    }

	if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		MatrixCRSByVector(m.val,m.col_ind,m.row_ptr, m.dx, m.dax, n); // результат занесён в  dax
	}
	if (iVar==TEMP) {
		MatrixCRSByVector(m.tval, m.tcol_ind,m.trow_ptr, m.tdx, m.tdax, n); // результат занесён в  dax
	}
	
	bool bOk_ZERO_rthdsd = false;
	bool bOk_Dirichlet = false;
	integer inumber_Dirichlet_node = 0;
	integer inumber_ZERO_rthdsd_node = 0;
	//bool bCheck_matrix = true;
	#pragma omp parallel for shared(dV,m,iVar,ifrontregulationgl)  schedule (guided)
	for (integer i=0; i<n; ++i) {
		if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			 if (ibackregulationgl!=nullptr) { 

				   // по новой нумерации с индексом i получает индекс старой нумерации iP
				   //iP=ifrontregulation[i];
				   m.ri[i]=dV[ifrontregulationgl[i]]-m.dax[i];
		           //m.roc[i]=m.ri[i];
                   m.roc[i]=1.0;
		     }
			 else {

				 if (b_on_adaptive_local_refinement_mesh) {
					 if ((iVar == PAM) && (fabs(dV[i]) < 1.0e-36)) {
						 if ((m.row_ptr[i + 1] - m.row_ptr[i] == 1) && (i == m.col_ind[m.row_ptr[i]])) {
							 bOk_Dirichlet = true;
							 inumber_Dirichlet_node++;
						 }
						 bOk_ZERO_rthdsd = true;
						 inumber_ZERO_rthdsd_node++;
					 }
					// if (iVar == PAM) {
						 //printf("ia=%lld\n",i);
						// for (integer j_7 = m.row_ptr[i]; j_7 <= m.row_ptr[i + 1] - 1; j_7++) {
							// printf("%lld ", m.col_ind[j_7]);
						 //}
						 //printf("\n");
						 //for (integer j_7 = m.row_ptr[i]; j_7 <= m.row_ptr[i + 1] - 1; j_7++) {
							// printf("%e ", m.val[j_7]);
							// if (m.col_ind[j_7] == i) {
								// if (m.val[j_7] <= 0.0) {
									// bCheck_matrix = false;
									 //printf("ERROR diagonal %e\n",m.val[j_7]);
								 //}								 
							 //}
							 //else {
								// if (m.val[j_7] >= 0.0) {
									// bCheck_matrix = false;
									 //printf("ERROR NO diagonal %e\n", m.val[j_7]);
								 //}
							 //}
						 //}
						// printf("\n");
						 //system("pause");
					 //}
				 }
		           m.ri[i]=dV[i]-m.dax[i];
				  // if (m.ri[i] != m.ri[i]) {
					//   printf(" m.ri[i]!= m.ri[i] solution bug. i=%lld\n",i);
					  // system("pause");
				   //}
		           //m.roc[i]=m.ri[i];
                   m.roc[i]=1.0;
				   //if (m.roc[i] != m.roc[i]) {
					 //  printf("m.roc[i]!=m.roc[i] solution bug.i=%lld \n",i);
					   //system("pause");
				   //}
			 }
		}
		if (iVar==TEMP) {
           m.tri[i]=dV[i]-m.tdax[i];
		   //m.troc[i]=m.tri[i];
		   m.troc[i]=1.0;
		}
	}

	doublereal norma_b= NormaV_for_gmres(dV, n);

	if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
	   delta0=NormaV(m.ri,n);
	}
	if (iVar==TEMP) {
		delta0=NormaV(m.tri,n);
	}
	if (0&&b_on_adaptive_local_refinement_mesh) {
		if (bOk_Dirichlet) {
			printf("b Diriclet Ok in %lld nodes\n", inumber_Dirichlet_node);
			system("pause");
		}
		if (bOk_ZERO_rthdsd) {
			printf("b ZERO RTHDSD  Ok in %lld nodes\n", inumber_ZERO_rthdsd_node);
			system("pause");
		}
	}
	//printf("debug %e\n",NormaV(dax,n)); // проверка на корректное составление СЛАУ
	//printf("%e \n",delta0); system("pause");
	//system("pause");
	// Если решение сразу хорошее то не считать:
	if (iVar==TEMP) {
		

		 if (fabs(delta0)<1.0e-4*dterminatedTResudual) iflag=0; 
	}
	else {
	   if (fabs(delta0)<dterminatedTResudual) iflag=0; 
	}
	integer iflag1=1;
	if (fabs(delta0)<1e-23) iflag1=0;
	if ((iVar == TEMP) && (iflag == 0) && (iflag1 == 0)) {
#if doubleintprecision == 1
		//printf("iflag=%lld, iflag1=%lld, delta0=%e\n", iflag, iflag1, delta0);
		std::cout << "iflag=" << iflag << ", iflag1=" << iflag1 << ", delta0=" << delta0 << std::endl;
#else
		//printf("iflag=%d, iflag1=%d, delta0=%e\n", iflag, iflag1, delta0);
		std::cout << "iflag=" << iflag << ", iflag1=" << iflag1 << ", delta0=" << delta0 << std::endl;
#endif
		 
		//system("pause");
		system("PAUSE");
	}
#if doubleintprecision == 1
	//printf("iflag=%lld, iflag1=%lld, delta0=%e\n",iflag,iflag1,delta0); system("pause");
#else
	//printf("iflag=%d, iflag1=%d, delta0=%e\n",iflag,iflag1,delta0); system("pause");
#endif
	

	/*if (iVar==PAM) {
	   // терминаьная невязка всегда на точность аппроксимации меньше стартовой невязки.
	   //if (e*dnew<e) e*=dnew;
	   doublereal me=sqrt(fabs(delta0))/n;
	   if (epsilon*me<epsilon) epsilon*=me;
	   //dterminatedTResudual=e;
	}*/
	
	

	//printf("delta0=%e\n",delta0);
	//system("pause");

	
	/*if (iflag != 0) {
       // конечная невязка всегда на точность аппроксимации меньше начальной невязки.
	   epsilon*=delta0; 
	   dterminatedTResudual=epsilon;
	}
	*/
	integer iN=10;
	if (n<=15000) {
		// задача очень малой размерности !
		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			iN = 1; // обязательно нужна хотя бы одна итерация.
					// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			//printf("%e\n",epsilon);
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon = 1.0e-3*fabs(delta0);
				//printf("%e\n", epsilon);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
			
		}
		if (iVar == TEMP) {
			iN = 2;
			//printf("%e\n", epsilon);
			epsilon = fmin(0.1*fabs(delta0), epsilon);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon=%e \n",epsilon);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016 20.05.2017
				//if (1.0e-3*fabs(delta0) < epsilon) {// не было.
					epsilon *= 1e-10;//1e-10
					iN = 20;
				//}
				//printf("%e\n", epsilon);
				//epsilon *= 1e-16;
				//iN = 30;
			}
		}
		if (iVar == PAM) {
			iN = 3; // решение для поправки давления должно быть получено точно.
			//printf("%e\n", epsilon);
			//1e-3->1e-4 30.03.2019
			if (1.0e-4*fabs(delta0)<epsilon) {
				epsilon = 1.0e-4*fabs(delta0);
				//printf("%e\n", epsilon);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
			//printf("%e",epsilon); system("pause");
			epsilon *= 1e-3;//1.0e-2
			
		}
	}
	else if ((n>15000)&&(n<30000)) {
		// задача очень малой размерности !
		if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		    iN=1; // обязательно нужна хотя бы одна итерация.
			// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon=1.0e-3*fabs(delta0);
			}
			if (iflag1==1) {
				iflag=1;
			}
		}
		if (iVar==TEMP) {
			iN=2;
			epsilon=fmin(0.1*fabs(delta0),epsilon);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon=%e \n",epsilon);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon *= 1e-10;
				iN = 120;
				//epsilon *= 1e-16;
				//iN = 30;
			}
		}
		if(iVar==PAM) {
			iN=3; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta0)<epsilon) {
				epsilon=1.0e-4*fabs(delta0);
			}
			if (iflag1==1) {
			   iflag=1;
			}
			epsilon *= 1e-3;//1.0e-2
			//printf("%e",epsilon); system("pause");
		}
	}
	else if ((n>=30000)&&(n<70000)) {
		// Здесь я немного увеличил число итераций и 
		// скорректировал условие окончания чтобы считало 
		// поточнее, но это не повлияло.
		// Главный вопрос в том что невязка по температуре почему-то не меняется.
		// задача небольшой размерности.
		if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		    iN=3; // обязательно нужна хотя бы одна итерация.
			// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon=1.0e-3*fabs(delta0);
			}
			if (iflag1==1) {
				iflag=1;
			}
			// 27.07.2016
			iN = 12;
			epsilon *= 1e-2;
		}
		if (iVar==TEMP) {
			iN=4;
			epsilon=fmin(0.1*fabs(delta0),epsilon);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon=%e \n",epsilon);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon *= 1e-10;
				iN = 120;
				//epsilon *= 1e-16;
				//iN = 30;
			}
		}
		if(iVar==PAM) {
			iN=6; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta0)<epsilon) {
				epsilon=1.0e-4*fabs(delta0);
			}
			if (iflag1==1) {
			   iflag=1;
			}
			//printf("%e",epsilon); system("pause");
			// 27.07.2016.
			//epsilon *= 1e-2;
			epsilon *= 1e-3;//1.0e-2
			iN = 20;
		}
	}
	else if ((n >= 70000) && (n<100000)) {
		// Здесь я немного увеличил число итераций и 
		// скорректировал условие окончания чтобы считало 
		// поточнее, но это не повлияло.
		// Главный вопрос в том что невязка по температуре почему-то не меняется.
		// задача небольшой размерности.
		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			iN = 3; // обязательно нужна хотя бы одна итерация.
					// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon = 1.0e-3*fabs(delta0);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
			// 27.07.2016
			iN = 12;
			epsilon *= 1e-2;
		}
		if (iVar == TEMP) {
			iN = 4;
			epsilon = fmin(0.1*fabs(delta0), epsilon);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon=%e \n",epsilon);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon *= 1e-10;
				iN = 20;
				//epsilon *= 1e-16;
				//iN = 30;
			}
		}
		if (iVar == PAM) {
			iN = 6; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta0)<epsilon) {
				epsilon = 1.0e-4*fabs(delta0);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
			//printf("%e",epsilon); system("pause");
			// 27.07.2016.
			epsilon *= 1e-2;
			iN = 20;
		}
	}
	else if ((n>=100000)&&(n<300000)) {
		// задача небольшой средней размерности.
		if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		    iN=3; // обязательно нужна хотя бы одна итерация.
			// Вообще говоря невязка для скоростей падает очень быстро поэтому всегда достаточно iN итераций для скорости.
			// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon=1.0e-3*fabs(delta0);
			}
			if (iflag1==1) {
				iflag=1;
			}
		}
		if (iVar==TEMP) {
			iN=4;
			epsilon=fmin(0.1*fabs(delta0),epsilon);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon=%e \n",epsilon);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon *= 1e-10;
				iN = 20;
				//epsilon *= 1e-16;
				//iN = 30;
			}
		}
		if(iVar==PAM) {
			iN=8; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta0)<epsilon) {
				epsilon=1.0e-4*fabs(delta0);
			}
			if (iflag1==1) {
			   iflag=1;
			}
			//printf("%e",epsilon); system("pause");
			// 27.07.2016.
			epsilon *= 1e-2;
		}
	}
	else if ((n>=300000)&&(n<1000000)) {
		// задача истинно средней размерности.
		if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		    iN=3; // обязательно нужна хотя бы одна итерация.
			// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon=1.0e-3*fabs(delta0);
			}
			if (iflag1==1) {
				iflag=1;
			}
		}
		if (iVar==TEMP) {
			iN=4;
			epsilon=fmin(0.1*fabs(delta0),epsilon);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon=%e \n",epsilon);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon *= 1e-10;
				iN = 20;
				//epsilon *= 1e-16;
				//iN = 30;
			}
		}
		if(iVar==PAM) {
			iN=16; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta0)<epsilon) {
				epsilon=1.0e-4*fabs(delta0);
			}
			if (iflag1==1) {
			   iflag=1;
			}
			//printf("%e",epsilon); system("pause");
			// 27.07.2016.
			epsilon *= 1e-2;
		}
	}
	else if ((n>=1000000)&&(n<3000000)) {
		// задача достаточно большой размерности.
		if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		    iN=6; // обязательно нужна хотя бы одна итерация.
			// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon=1.0e-3*fabs(delta0);
			}
			if (iflag1==1) {
				iflag=1;
			}
		}
		if (iVar==TEMP) {
			iN=8;
			epsilon=fmin(0.1*fabs(delta0),epsilon);
			if (bSIMPLErun_now_for_temperature  ) {
				//printf("epsilon=%e \n",epsilon);
				//system("pause");
				// Экспериментальным образом обнаружена недоэтерированость по температуре для гидродинамического решателя.
				// поэтому точность было решено увеличить на 5 порядков.
				// 27.07.2016
				epsilon *= 1e-10;
				iN = 20;
				//epsilon *= 1e-16;
				//iN = 30;
			}
		}
		if(iVar==PAM) {
			iN=23; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta0)<epsilon) {
				epsilon=1.0e-4*fabs(delta0);
			}
			if (iflag1==1) {
			   iflag=1;
			}
			//printf("%e",epsilon); system("pause");
			// 27.07.2016.
			epsilon *= 1e-2;
		}
	}
	else if (n>=3000000) {
		// задача очень большой размерности.
		if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		    iN=6; // обязательно нужна хотя бы одна итерация.
			// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon=1.0e-3*fabs(delta0);
			}
			if (iflag1==1) {
				iflag=1;
			}
		}
		if (iVar==TEMP) {
			iN=8;
			epsilon=fmin(0.1*fabs(delta0),epsilon);
		}
		if(iVar==PAM) {
			iN=36; // решение для поправки давления должно быть получено точно.
			if (1.0e-4*fabs(delta0)<epsilon) {
				epsilon=1.0e-4*fabs(delta0);
			}
			if (iflag1==1) {
			   iflag=1;
			}
			//printf("%e",epsilon); system("pause");
		}
	}
	
	
	if (iVar==TEMP) {
		maxit=3*4*m.icount_vel;
		//if (m.icount_vel<iN) {
		  //iN=m.icount_vel;
	//	}
		// 2000
		maxit=2000; // также увеличим максимальное число итераций при нахождении поля температур, т.к. мы его
		// решаем без нижней релаксации в общем случае в случае задачи сопряжённого теплообмена на огромной по размерности расчётной сетке.
#if doubleintprecision == 1
		//printf("iN=%lld\n",iN);
#else
		//printf("iN=%d\n",iN);
#endif
		
		//system("pause");
	}
	if (iVar==PAM) {
		// 90 процентов времени тратится на решение уравнения Пуассона.
		maxit=3*9*m.icount_vel;
		if (3*9*m.icount_vel<iN) {
		  //iN=3*9*m.icount_vel;
		}
		// Для течений с большими числами Рейнольдса или большим сеточным разрешением требуется большое количество итераций
		// для сходимости при решении уравнения на поправку давления.
		maxit=2000; // 2000
		//iN = 2000;
		// 23.01.2021 закоментировал все операторы break;
		//epsilon = 1.0e-10;
		//dterminatedTResudual = 1.0e-10;

		//if (icount_first_PAM_iteration_SIMPLE_algorithm_global < 6) {
			//maxit = 80;
			//icount_first_PAM_iteration_SIMPLE_algorithm_global++;
		//}
		if ((maxit==0)&&(iN==0)) {
			//maxit=iN=81;
			//maxit=iN=1000;
			/*
			if (iVar==PAM) {
			    for (icount=0; icount<20; icount++) {
			        PAMGSP(sl, slb, dX0, dV, maxelm, maxbound);
			    }
			}
			*/
		/*
			for (integer i9=m.row_ptr[maxelm-2]; i9<m.row_ptr[n]; i9++) {
			#if doubleintprecision == 1
				if (m.col_ind[i9]>=maxelm) {
					printf("boundary val=%e col_ind=%lld\n",m.val[i9],m.col_ind[i9]);
				}
				else {
					printf("internal val=%e col_ind=%lld\n",m.val[i9],m.col_ind[i9]);
				}
			#else
				if (m.col_ind[i9]>=maxelm) {
					printf("boundary val=%e col_ind=%d\n",m.val[i9],m.col_ind[i9]);
				}
				else {
					printf("internal val=%e col_ind=%d\n",m.val[i9],m.col_ind[i9]);
				}
			#endif
				
				system("pause");
			}
			
			for (integer i9=0; i9<=n; i9++) {
			#if doubleintprecision == 1
				printf("%lld %lld\n",m.row_ptr[i9],i9);
			#else
				printf("%d %d\n",m.row_ptr[i9],i9);
			#endif
			
			system("pause");
			}
		*/
			//Bi_CGStabCRS(n, m.val, m.col_ind, m.row_ptr, dV, dX0, 200);

			//BiSoprGradCRS( m.val, m.col_ind, m.row_ptr,dV,dX0,n,200);
		}
	}
	if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		maxit=100;//100
	}

	
	/*if (iVar==PAM) {
		printf(" %1.2e ",delta0);
		#if doubleintprecision == 1
			printf("icount=%lld iN=%lld iflag1=%lld iflag=%lld maxit=%lld",icount,iN,iflag1,iflag,maxit);
		#else
			printf("icount=%d iN=%d iflag1=%d iflag=%d maxit=%d",icount,iN,iflag1,iflag,maxit);
		#endif
		
		system("pause");
	}*/

	/*
	if (iVar==TEMP) {
		//BiSoprGradCRS( m.tval, m.tcol_ind, m.trow_ptr,dV,dX0,n,200);
    	Bi_CGStabCRS(n, m.tval, m.tcol_ind, m.trow_ptr, dV, dX0, 200);
	}
	else {
    	//BiSoprGradCRS( m.val, m.col_ind, m.row_ptr,dV,dX0,n,200);
    	Bi_CGStabCRS(n, m.val, m.col_ind, m.row_ptr, dV, dX0, 200);
	}
	*/
	
	// диагностическое сообщение какую переменную мы решаем.
	//switch (iVar) {
	//case PAM: printf("PAM\n");  break;
	//case VX:  printf("VX\n"); break;
	//case VY:  printf("VY\n"); break;
	//case VZ:  printf("VZ\n"); break;
	//case NUSHA:  printf("NU\n"); break;
	//case TURBULENT_KINETIK_ENERGY:  printf("TURBULENT_KINETIK_ENERGY\n"); break;
	//case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:  printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA\n"); break;
	//case TURBULENT_KINETIK_ENERGY_STD_K_EPS:  printf(" TURBULENT_KINETIK_ENERGY_STD_K_EPS \n"); break;
	//case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:  printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS\n"); break;
	//case RE_THETA_LANGTRY_MENTER: printf("TURBULENT_RE_THETA_LANGTRY_MENTER \n"); break;
	//case GAMMA_LANGTRY_MENTER: printf("TURBULENT_GAMMA_LANGTRY_MENTER \n");	break;
	//case TEMP:  printf("TEMP\n"); break;
	//}

	// Если число расходимостей превысит оговорённую константу то произойдёт выход из алгоритма.
	integer i_signal_break_pam_opening = 0;
	// x хорошее значение.
	const integer i_limit_signal_pam_break_opening = 4000;//20
	doublereal delta_old_iter = 1.0e10;

	integer count_iter_for_film_coef = 0;

	int ipovtor = 0;
	doublereal previos_res = 1.0;


	doublereal dres00 = 0.0;

	if (iVar == TEMP) {
		dres00 = (NormaV_for_gmres(m.tri, n) / norma_b);
	}
	else {
		dres00 = (NormaV_for_gmres(m.ri, n) / norma_b);
	}

	if (iVar == TEMP) {
		//printf("epsilon=%e delta0=%e\n", epsilon, delta0);
		//system("pause");
	}
	//printf("while (((icount < iN=%lld) && (iflag1=%lld != 0)) || (iflag=%lld != 0 && icount < maxit=%lld))\n",iN,iflag1,iflag,maxit);

	// Мы обязательно должны сделать несколько итераций. (не менее 10).
	// Если только решение не удовлетворяет уравнению тождественно.
	while (((icount < iN) && (iflag1 != 0)) || (iflag != 0 && icount < maxit)) {

		icount++;

		count_iter_for_film_coef++;
		// В случае задачи Ньютона - Рихмана, Стефана-Больцмана и миксового условия не итерируем до конца обрываем, 
		// т.к. нам требуется частая пересборка матрицы. 13 марта 2016.
		//if (((adiabatic_vs_heat_transfer_coeff > ADIABATIC_WALL_BC) || (breakRUMBAcalc_for_nonlinear_boundary_condition)) && (count_iter_for_film_coef>5)) break;


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			if (!isfinite(Scal(m.roc, m.ri, n))) {
				roi = 0.0;
			}
			else {
				roi = Scal(m.roc, m.ri, n);
			}
			if (roi != roi) {
				printf("roi!=roi solution bug. \n");
				system("pause");
			}
			if (fabs(wi) < 1.0e-30) {
				if (fabs(roim1) < 1.0e-30) {
					bet = 1.0;
				}
				else {
					bet = (roi / roim1);
				}
			}
			else {
				if (fabs(roim1) < 1.0e-30) {
					bet = (al / wi);
				}
				else {
					bet = (roi / roim1)*(al / wi);
				}
			}
			if ((bet != bet)||(!isfinite(bet))) {
				printf("bet!=bet solution bug. \n");
				printf("%e %e %e %e\n", roi, roim1, al, wi);
				system("pause");
			}

			//printf("%e %e %e %e\n",roi,roim1,al,wi);
			//system("pause");

#pragma omp parallel for shared(m,wi,bet)  schedule (guided)
			for (integer i = 0; i < n; ++i) {
				doublereal pibuf = m.ri[i] + (m.pi[i] - m.vi[i] * wi)*bet;
				if (pibuf != pibuf) {
					printf("pibuf!=pibuf solution bug. \n");
					printf("ri=%e pi=%e vi=%e wi=%e bet=%e\n", m.ri[i], m.pi[i], m.vi[i], wi, bet);
					system("pause");
				}
				m.pi[i] = pibuf;
			}
		}
		if (iVar == TEMP) {
			roi = Scal(m.troc, m.tri, n);
			bet = (roi / roim1)*(al / wi);

#pragma omp parallel for shared(m,wi,bet) schedule (guided)
			for (integer i = 0; i < n; ++i) {
				doublereal pibuf = m.tri[i] + (m.tpi[i] - m.tvi[i] * wi)*bet;
				m.tpi[i] = pibuf;
			}
		}

		// Ky=pi

		// (LU)y=pi; 
		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for shared(m)  schedule (guided)
			for (integer i = 0; i < n; ++i) m.y[i] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

			//  9 августа 2015 при внедрении перенумерации узлов nested desection
			if (bpam_gsp && (iVar == PAM)) {
				if (ibackregulationgl != nullptr) {
					PAMGSPnd(sl, slb, m.y, m.pi, maxelm, maxbound, ifrontregulationgl);
				}
				else {
					PAMGSP(sl, slb, m.y, m.pi, maxelm, maxbound);
				}
			}
			else {

				if (brc) {
#pragma omp parallel for
					for (integer i7 = 0; i7 < n; ++i7)
					{
						m.vec[i7] = m.pi[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < m.iwk + 2; ++i7)
					{
						m.alurc[i7] = m.alu[i7];
						m.jlurc[i7] = m.jlu[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < n + 2; ++i7) {
						m.jurc[i7] = m.ju[i7];

					}
				}

				if (ibackregulationgl != nullptr) {
					printf("ibackregulationgl != nullptr\n");
					system("pause");
					//lusol_2(n, m.pi, m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=pi;
					lusol_3(n, m.pi, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=pi;
				}
				else {
					lusol_(n, m.pi, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=pi;

				}

				if (brc) {
#pragma omp parallel for
					for (integer i7 = 0; i7 < n; ++i7)
					{
						m.pi[i7] = m.vec[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < m.iwk + 2; ++i7)
					{
						m.alu[i7] = m.alurc[i7];
						m.jlu[i7] = m.jlurc[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < n + 2; ++i7)
					{
						m.ju[i7] = m.jurc[i7];
					}
				}

			}

			MatrixCRSByVector(m.val, m.col_ind, m.row_ptr, m.y, m.vi, n); // vi==A*y;
		}
		if (iVar == TEMP) {
			// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for shared(m)  schedule (guided)
			for (integer i = 0; i < n; ++i) m.ty[i] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

			lusol_(n, m.tpi, m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=tpi;
			MatrixCRSByVector(m.tval, m.tcol_ind, m.trow_ptr, m.ty, m.tvi, n); // vi==A*y;
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {

			if ((fabs(roi) < 1e-30) && (fabs(Scal(m.roc, m.vi, n)) < 1e-30)) {
				al = 1.0;
			}
			else if (fabs(roi) < 1e-30) {
				al = 0.0;
			}
			else {
				al = roi / Scal(m.roc, m.vi, n);
			}
			if (al != al) {
				printf("roi!=roi solution bug. \n");
				system("pause");
			}

#pragma omp parallel for shared(m,al)  schedule (guided)
			for (integer i = 0; i < n; ++i) {
				m.s[i] = m.ri[i] - al * m.vi[i];
				if (m.s[i] != m.s[i]) {
					printf("m.s[i]!=m.s[i] solution bug. i==%lld \n",i);
					system("pause");
				}
			}
		}
		if (iVar == TEMP) {

			if ((fabs(roi) < 1e-30) && (fabs(Scal(m.troc, m.tvi, n)) < 1e-30)) {
				al = 1.0;
			}
			else if (fabs(roi) < 1e-30) {
				al = 0.0;
			}
			else {
				al = roi / Scal(m.troc, m.tvi, n);
			}

#pragma omp parallel for shared(m,al)  schedule (guided)
			for (integer i = 0; i < n; ++i) {
				m.ts[i] = m.tri[i] - al * m.tvi[i];
			}
		}
		// Kz=s

		// (LU)z=s; 
		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for shared(m)  schedule (guided)
			for (integer i = 0; i < n; ++i) m.z[i] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

			if (bpam_gsp && (iVar == PAM)) {
				printf("neponqtnoe perekluchenie na GZ.\n");
				printf("in BiCGStab_internal3 27.07.2016.");
				//system("pause");
				system("PAUSE");
				if (ibackregulationgl != nullptr) {
					PAMGSPnd(sl, slb, m.z, m.s, maxelm, maxbound, ifrontregulationgl);
				}
				else {
					PAMGSP(sl, slb, m.z, m.s, maxelm, maxbound);
				}
			}
			else {

				if (brc) {
#pragma omp parallel for
					for (integer i7 = 0; i7 < n; ++i7)
					{
						m.vec[i7] = m.s[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < m.iwk + 2; ++i7)
					{
						m.alurc[i7] = m.alu[i7];
#ifdef MY_DEBUG_NOT_NUMBER
						if (m.alurc[i7] != m.alurc[i7]) {
							printf("m.alurc[i7]!=m.alurc[i7] solution bug. i7=%lld\n",i7);
							system("pause");
						}
#endif
						m.jlurc[i7] = m.jlu[i7];
#ifdef MY_DEBUG_NOT_NUMBER
						if (m.jlurc[i7] != m.jlurc[i7]) {
							printf("m.jlurc[i7]!=m.jlurc[i7] solution bug. i7=%lld\n",i7);
							system("pause");
						}
#endif
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < n + 2; ++i7) 
					{
						m.jurc[i7] = m.ju[i7];
#ifdef MY_DEBUG_NOT_NUMBER
						if (m.jurc[i7] != m.jurc[i7]) {
							printf("m.jurc[i7]!=m.jurc[i7] solution bug. i7=%lld\n",i7);
							system("pause");
						}
#endif
					}
				}

				if (ibackregulationgl != nullptr) {
					//lusol_2(n, m.s, m.z, m.alu, m.jlu, m.ju,  m.x1, maxelm); // M*y=pi;
					lusol_3(n, m.s, m.z, m.alu, m.jlu, m.ju, maxelm); // M*y=pi;
				}
				else {
					lusol_(n, m.s, m.z, m.alu, m.jlu, m.ju, maxelm); // Mz=s;
				}

				if (brc) {
#pragma omp parallel for
					for (integer i7 = 0; i7 < n; ++i7)
					{
						m.s[i7] = m.vec[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < m.iwk + 2; ++i7) {
						m.alu[i7] = m.alurc[i7];
						m.jlu[i7] = m.jlurc[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < n + 2; ++i7)
					{
						m.ju[i7] = m.jurc[i7];
					}
				}
			}
			MatrixCRSByVector(m.val, m.col_ind, m.row_ptr, m.z, m.t, n); // t==A*z;
#ifdef MY_DEBUG_NOT_NUMBER
			for (integer i7 = 0; i7 < n; ++i7) {
				if (m.t[i7] != m.t[i7]) {
					printf("m.t[%lld]=%e  iVar=%lld\n",i7,m.t[i7],iVar);
					system("pause");
				}
			}
#endif
		}
		if (iVar == TEMP) {
			// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for shared(m)  schedule (guided)
			for (integer i = 0; i < n; ++i) m.tz[i] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

			lusol_(n, m.ts, m.tz, m.talu, m.tjlu, m.tju, maxelm);
			MatrixCRSByVector(m.tval, m.tcol_ind, m.trow_ptr, m.tz, m.tt, n); // t==A*z;
		}

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {

			//wi = Scal(m.t, m.s, n) / Scal(m.t, m.t, n);
			if ((fabs(Scal(m.t, m.s, n)) < 1e-30) && (fabs(Scal(m.t, m.t, n)) < 1e-30)) {
				wi = 1.0;
			}
			else if (fabs(Scal(m.t, m.s, n)) < 1e-30) {
				wi = 0.0;
			}
			else {
				doublereal ts = Scal(m.t, m.s, n);
				doublereal tt = Scal(m.t, m.t, n);
				if ((ts!=ts)||(tt!=tt)) {
					wi = 0.0;					
				}
				else {
					wi = Scal(m.t, m.s, n) / Scal(m.t, m.t, n);
				}
			}

#ifdef MY_DEBUG_NOT_NUMBER
			if (wi != wi) {
				printf("wi!=wi solution bug. \n");
				printf("(t,s)=%e (t,t)=%e", Scal(m.t, m.s, n), Scal(m.t, m.t, n));
				printf("iVar=%lld\n",iVar);
				for (integer i87 = 0; i87 < n; i87++) {
					if (m.t[i87] != m.t[i87]) {
						printf("m.t[%lld]=%e\n",i87,m.t[i87]);
						system("pause");
				}
					if (m.s[i87] != m.s[i87]) {
						printf("m.s[%lld]=%e\n", i87, m.s[i87]);
						system("pause");
					}
				}
				system("pause");
			}
#endif

			// printf("%e %e",Scal(m.t,m.s,n),Scal(m.t,m.t,n));

#pragma omp parallel for shared(m, al, wi)  schedule (guided)
			for (integer i = 0; i < n; ++i) {
				//dx[i]+=al*pi[i]+wi*s[i]; // так было без предобуславливателя
				m.dx[i] += al * m.y[i] + wi * m.z[i]; // так стало с предобуславливателем
				m.ri[i] = m.s[i] - wi * m.t[i];
			}
			deltai = NormaV(m.ri, n);

			if (b_on_adaptive_local_refinement_mesh) {
				//if (icount >= 96) system("pause");
				//if ((deltai > delta_old_iter)&&(icount>=3)) break;
				//system("pause");
			}

			for (integer i = 0; i < n; ++i) {
				// Обязательно делаем присваивание после того как проверена сходимость.
				//m.dx[i] += al * m.y[i] + wi * m.z[i]; // так стало с предобуславливателем
			}
		}
		if (iVar == TEMP) {
			wi = Scal(m.tt, m.ts, n) / Scal(m.tt, m.tt, n);

#pragma omp parallel for shared(m, al, wi)  schedule (guided)
			for (integer i = 0; i < n; ++i) {
				//dx[i]+=al*pi[i]+wi*s[i]; // так было без предобуславливателя
				m.tdx[i] += al * m.ty[i] + wi * m.tz[i]; // так стало с предобуславливателем
				m.tri[i] = m.ts[i] - wi * m.tt[i];
			}
			deltai = NormaV(m.tri, n);
		}
		//printf("icount=%lld deltai=%e\n",icount, deltai); system("pause");

		//if (breakRUMBAcalc_for_nonlinear_boundary_condition && (icount>10) && (bonly_solid_calculation)&& (iVar == TEMP)) iflag = 0;

		if ((iVar == TEMP) && (bonly_solid_calculation))
		{
			//printf("%lld %e\n", icount, deltai);
			std::cout << icount << " " << deltai << std::endl;
			
			//system("pause");
		}

		if (iVar == TEMP) {
			if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE) {

				// 02.12.2021
			    // Ахтунг, ахтунг !!!
			    // При нестационарном моделировании ни в коем случае нельзя обрывать вычисления.

				// Нестационарное моделирование. Считаем до наступления сходимости.
				previos_res = deltai;
			}
			else {
				// против стагнации в надежде получить решение.
				if (deltai > previos_res) {
					ipovtor++;
				}
				//else ipovtor = 0;
				previos_res = deltai;
				if (ipovtor > 10) break;
			}
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) {
			doublereal dres0 = (NormaV_for_gmres(m.dx, n) / norma_b);
			if (dres0 / dres00 < 0.1) {
				break;
			}
		}

		// печать невязки на консоль
		if (0&&bprintmessage)
		{
			/*
			switch (iVar) {
			case PAM: printf("\nPAM\n"); break;
			case VX: printf("\nVX\n"); break;
			case VY: printf("\nVY\n"); break;
			case VZ: printf("\nVZ\n"); break;
			case TEMP: printf("\nTEMP\n"); break;
			}
			*/
            if ((icount % 10) == 0)  {
				printf("iter  residual\n");
			}
			
			std::cout << icount << " " << deltai << std::endl;

			//system("pause");
		}
		// 28.07.2016.
		//std::cout << icount << " " << deltai << std::endl;
		
		//system("pause");
		

		if (deltai > delta_old_iter) i_signal_break_pam_opening++;

		// Если номер итерации больше третьей и процесс расходится и мы решаем поправку давления то стоп.
		// 31.10.2019
		if ((iVar == PAM)&&(icount>3)&&(deltai > delta_old_iter)) {
			// Естественная конвекция на АЛИС сетке в случае кода заданы только температуры а не мощности.
			if (inumber_iteration_SIMPLE < 10) {
				break;
			}
		}

		delta_old_iter = deltai;
		if (iVar == PAM) {
			if (i_signal_break_pam_opening > i_limit_signal_pam_break_opening) {
				// досрочный выход из цикла.
#if doubleintprecision == 1
				printf("icount PAM=%lld\n", icount);
#else
				printf("icount PAM=%d\n", icount);
#endif
				if (!b_on_adaptive_local_refinement_mesh) {
					break;
				}
				//system("pause");
			}
		}

		// Досрочный выход из итерационного процесса по опыту алгоритма FGMRES
		// Ю. Саада и М. Шульца.
		if (0&&((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER))) {
			// Нужно точнее, этой точности недостаточно
			if ((NormaV_for_gmres(m.ri, n) / norma_b) <= dterminatedTResudual) {
				iflag = 0; // конец вычисления
						   // 20.05.2017
				iflag1 = 0; // Выход по второму критерию чтоб не делать тысячи итераций.
				//printf("dosrochnji vjhod\n");
			}
		}
		if (iVar == TEMP) {
			if ((bSIMPLErun_now_for_temperature) && (bSIMPLErun_now_for_natural_convection)) {
				//printf("break %e %e\n", (NormaV_for_gmres(m.tri, n) / norma_b), dterminatedTResudual);
			
				// Внимание при естественной конвекции не менее 1.0e-12 иначе недоитерированность сильнейшая.
				// При недоитерированности решение совершенно нефизично.
				if ((NormaV_for_gmres(m.tri, n) / norma_b) <= 1.0e-12) {
					iflag = 0; // конец вычисления
							   // 20.05.2017
					iflag1 = 0; // Выход по второму критерию чтоб не делать тысячи итераций.
					//printf("dosrochnji vjhod\n");
				}
			}
			else if (bSIMPLErun_now_for_temperature) {
				// Внимание при вынужденной конвекции не менее 1.0e-12 иначе недоитерированность.
				// При недоитерированности решение может быть нефизично.
				if ((NormaV_for_gmres(m.tri, n) / norma_b) <= 1.0e-12) {
					iflag = 0; // конец вычисления
							   // 20.05.2017
					iflag1 = 0; // Выход по второму критерию чтоб не делать тысячи итераций.
					//printf("dosrochnji vjhod\n");
				}
			}
			else {
				if ((NormaV_for_gmres(m.tri, n) / norma_b) <= dterminatedTResudual) {
					iflag = 0; // конец вычисления
							   // 20.05.2017
					iflag1 = 0; // Выход по второму критерию чтоб не делать тысячи итераций.
					//printf("dosrochnji vjhod\n");
				}
			}
		}

		if (deltai < epsilon) {
			iflag = 0; // конец вычисления
			// 20.05.2017
			iflag1 = 0; // Выход по второму критерию чтоб не делать тысячи итераций.
		}
		else roim1=roi;

		if (iVar == TEMP) {
#if doubleintprecision == 1
			//printf("epsilon=%e deltai=%e icount=%lld\n",epsilon,deltai, icount);
#else
			//printf("epsilon=%e deltai=%e icount=%d\n",epsilon,deltai, icount);
#endif
			
			//system("pause");
		}

	}

    if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		if (!((maxit==0)&&(iN==0))) {
			if (ibackregulationgl!=nullptr) {
				#pragma omp parallel for shared(dX0, m) schedule (guided)
	            for (integer i=0; i<n; ++i) dX0[ifrontregulationgl[i]]=m.dx[i];
			}
			else {
	            #pragma omp parallel for shared(dX0, m)  schedule (guided)
	            for (integer i=0; i<n; ++i) dX0[i]=m.dx[i];
			}
		}
    }
    if (iVar==TEMP) {
	    #pragma omp parallel for shared(dX0, m) schedule (guided)
	    for (integer i=0; i<n; ++i) dX0[i]=m.tdx[i];
    }
	
	// Восстанавливаем обычный порядок значений.
	if (breordering_for_parallel) {
		if ((my_amg_manager.lfil < 3) && (number_cores() > 1)) {

			doublereal* tmp = new doublereal[n];
			integer* old_number = new integer[n];
			for (integer i_1 = 0; i_1 < n; ++i_1) {
				old_number[new_order[i_1]] = i_1;
			}

			for (integer i_1 = 0; i_1 < n; ++i_1) {
				tmp[old_number[i_1]] = dV[i_1];
			}
			for (integer i_1 = 0; i_1 < n; ++i_1) {
				dV[i_1] = tmp[i_1];
			}

			for (integer i_1 = 0; i_1 < n; ++i_1) {
				tmp[old_number[i_1]] = dX0[i_1];
			}
			for (integer i_1 = 0; i_1 < n; ++i_1) {
				dX0[i_1] = tmp[i_1];
			}

			delete[] old_number;
			delete[] tmp;


			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
				(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
				(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {


				for (integer i_1 = 0; i_1 < m.row_ptr[n]; ++i_1) {
					m.val[i_1] = valcopy[i_1];
					m.col_ind[i_1] = col_indcopy[i_1];
				}

				for (integer i_1 = 0; i_1 < n + 1; ++i_1) {
					m.row_ptr[i_1] = row_ptrcopy[i_1];
				}


			}
			if (iVar == TEMP) {

				for (integer i_1 = 0; i_1 < m.trow_ptr[n]; ++i_1) {
					m.tval[i_1] = valcopy[i_1];
					m.tcol_ind[i_1] = col_indcopy[i_1];
				}

				for (integer i_1 = 0; i_1 < n + 1; ++i_1) {
					m.trow_ptr[i_1] = row_ptrcopy[i_1];
				}
			}


			delete[] valcopy;
			delete[] col_indcopy;
			delete[] row_ptrcopy;
		}
	}

	// Это матрица в котрой нумерация (а индексация элементов с нуля) начинается с единицы. Она используется в библиотеке SPARSKIT2.
	if ((iVar==VELOCITY_X_COMPONENT)||(iVar==VELOCITY_Y_COMPONENT)||(iVar==VELOCITY_Z_COMPONENT)||(iVar==PAM) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
	   if (m.bsignalfreeCRScfd) {
		   // Это таже CRS матрица что и a,ja, ia только элементы в ней нумеруются также как и индексируются с нуля.
	       if (m.val!=nullptr) delete[] m.val;
		   if (m.col_ind!=nullptr) delete[] m.col_ind;
		   if (m.row_ptr!=nullptr) delete[] m.row_ptr; 
	       if (m.a!=nullptr) delete[] m.a;
		   if (m.ja!=nullptr) delete[] m.ja;
		   if (m.ia!=nullptr) delete[] m.ia; // уничтожаем матрицу в CRS формате.
		    // освобождение памяти
	       if (m.ri!=nullptr) delete[] m.ri;
		   if (m.roc!=nullptr) delete[] m.roc;
		   if (m.s!=nullptr) delete[] m.s;
		   if (m.t!=nullptr) delete[] m.t;
	       if (m.vi!=nullptr) delete[] m.vi;
		   if (m.pi!=nullptr) delete[] m.pi;
		   if (m.dax!=nullptr) delete[] m.dax;
	       if (m.y!=nullptr) delete[] m.y;
		   if (m.z!=nullptr) delete[] m.z;
		   if (m.dx!=nullptr) delete[] m.dx;
		   if (m.vec!=nullptr) delete[] m.vec;
		   // alu, jlu - MSR матрица неполного ILU разложения.
	       // ju - указатель на диагональные элементы, iw - вспомогательный вектор.
		   if (m.alu!=nullptr) delete[] m.alu; 
		   if (m.jlu!=nullptr) delete[] m.jlu;
		   if (m.ju!=nullptr) delete[] m.ju;
		   if (ibackregulationgl!=nullptr) {
		       if (m.alu1!=nullptr) delete[] m.alu1; 
		       if (m.jlu1!=nullptr) delete[] m.jlu1;
		       if (m.ju1!=nullptr) delete[] m.ju1;
			   if (m.x1!=nullptr) delete[] m.x1;
		   }
		   if (m.alurc!=nullptr) delete[] m.alurc; 
		   if (m.jlurc!=nullptr) delete[] m.jlurc;
		   if (m.jurc!=nullptr) delete[] m.jurc;
		   if (itype_ilu == ILU0) {
			   if (m.iw != nullptr) delete[] m.iw; // удаляем рабочий массив.
		   }
		   // Освобождение памяти.
		   if (itype_ilu == ILU_lfil) {
			   if (m.w != nullptr) delete[] m.w;
			   if (m.jw != nullptr) delete[] m.jw;
			   if (m.w != nullptr) delete[] m.w_dubl;
			   if (m.jw != nullptr) delete[] m.jw_dubl;
			   if (m.levs != nullptr) delete[] m.levs;
		   }
	       m.bsignalfreeCRScfd=false; // память корректно освобождена.
	   }
	}
	 if (iVar==TEMP) {
       if (m.bsignalfreeCRSt) {
		   // Это таже CRS матрица что и a,ja, ia только элементы в ней нумеруются также как и индексируются с нуля.
	       if (m.tval!=nullptr) delete[] m.tval;
		   if (m.tcol_ind!=nullptr) delete[] m.tcol_ind;
		   if (m.trow_ptr!=nullptr) delete[] m.trow_ptr; 
	       if (m.ta!=nullptr) delete[] m.ta;
		   if (m.tja!=nullptr) delete[] m.tja;
		   if (m.tia!=nullptr) delete[] m.tia; // уничтожаем матрицу в CRS формате.
		    // освобождение памяти
	       if (m.tri!=nullptr) delete[] m.tri;
		   if (m.troc!=nullptr) delete[] m.troc;
		   if (m.ts!=nullptr) delete[] m.ts;
		   if (m.tt!=nullptr) delete[] m.tt;
	       if (m.tvi!=nullptr) delete[] m.tvi; 
		   if (m.tpi!=nullptr) delete[] m.tpi;
		   if (m.tdax!=nullptr) delete[] m.tdax;
	       if (m.ty!=nullptr) delete[] m.ty; 
		   if (m.tz!=nullptr) delete[] m.tz;
		   if (m.tdx!=nullptr) delete[] m.tdx;
		   // alu, jlu - MSR матрица неполного ILU разложения.
	       // ju - указатель на диагональные элементы, iw - вспомогательный вектор.
		   if (m.talu!=nullptr) delete[] m.talu;
		   if (m.tjlu!=nullptr) delete[] m.tjlu; 
		   if (m.tju!=nullptr) delete[] m.tju;
		   if (itype_ilu == ILU0) {
			   if (m.tiw != nullptr) delete[] m.tiw; // удаляем рабочий массив.
		   }
		   // Освобождение памяти.
			   if (itype_ilu == ILU_lfil) {
				   if (m.tw != nullptr) delete[] m.tw;
				   if (m.tjw != nullptr) delete[] m.tjw;
				   if (m.tlevs != nullptr) delete[] m.tlevs;
			   }
	       m.bsignalfreeCRSt=false; // память корректно освобождена.
	   }
	 }	 

	 if (b_on_adaptive_local_refinement_mesh) {
		// delete[] color;
	 }
	 delete[] new_order;

	 delete[] valcopy;
	 delete[] col_indcopy;
	 delete[] row_ptrcopy;

	 BiCGStab_internal3_incomming_now = false;

	 // Шаманство.
	 //if (iVar==VX) {
		// m.icount_vel=icount;
	// }
	 //if ((iVar==VY)||(iVar==VZ)) {
	//	 if (icount>m.icount_vel) {
		//	 m.icount_vel=icount;
		 //}
	 //}
	 
	// отладочная печать количества сделанных итераций.
	 /*
	 #if doubleintprecision == 1
		printf("%lld ",icount);
	 #else
		printf("%d ",icount);
	 #endif
	
	if(iVar==PAM) {
		printf("%e ",deltai/delta0);
		//printf("%e %e\n",deltai,delta0);
		//system("pause");
	}
	*/
	//if (iVar==TEMP) printf(" %e ",deltai/delta0);
	 //system("pause");
} // Bi_CGStab_internal3


  // Здесь содержится обвязка вызывающая amg1r5.
  // локальное выделение памяти:всё внутри, многократные alloc и free.
// amg1r5 для решения задачи напряженно-деформированнного состояния.
// Нет сходимости. Дата успешного подключения 24 сентября 2017.
void amg_loc_memory_Stress(SIMPLESPARSE &sparseM, integer n,
	doublereal *dV, doublereal* &dX0,
	QuickMemVorst& m)
{
	// Замер времени.
	unsigned int calculation_main_start_time; // начало счёта мс.
	unsigned int calculation_main_end_time; // окончание счёта мс.

	calculation_main_start_time = clock(); // момент начала счёта.

		
	doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля
										   
										   // На случай если память не была выделена.
	if (dX0 == nullptr) {
		dX0 = new doublereal[n];
		for (integer i = 0; i<n; ++i) {
			dX0[i] = 0.0;
		}
	}

	const integer id = 0;

	simplesparsetoCRS(sparseM, m.val, m.col_ind, m.row_ptr, n); // преобразование матрицы из одного формата хранения в другой.
																//m.ballocCRScfd = true;
	simplesparsefree(sparseM, n);

	integer ierr = 0;
	doublereal eps = 1.0e-12;

	ierr = 0; // изначальное состояние безошибочное.
			  // Порог точности решения СЛАУ. Значение 1.0E-12 достаточно что проверено в ANSYS icepak.
	eps = 1.0e-12; // рекомендуемое значение которого достаточно. 

				   // Требования к оперативной памяти.
				   /*     VECTOR         NEEDED LENGTH (GUESS) */
				   /*       A               3*NNA + 5*NNU */
				   /*       JA              3*NNA + 5*NNU */
				   /*       IA              2.2*NNU */
				   /*       U               2.2*NNU */
				   /*       F               2.2*NNU */
				   /*       IG              5.4*NNU */

	integer nna = 0;
	for (integer k = 0; k < m.row_ptr[n]; ++k) {

		if (fabs(m.val[k]) > nonzeroEPS) {
			nna++;
		}
	}
	//printf("nna=%d istinnoe=%d\n", nna, m.row_ptr[n]);
	//system("pause");
	//integer nna = m.row_ptr[n]; // количество ненулевых элементов в матрице СЛАУ.



	integer nnu = n; // число неизвестных.

					 // данная константа работоспособна вплоть до размерностей сетки равных 34млн 463тысячи 250узлов.
					 //doublereal rsize=1.51; // 1048416
					 // Вынужденные течения достаточно 2.5.
					 // значения 3.5 недостаточно для 8 модулей Пионер. 
	doublereal rsize = 4.5; // на задаче Концевого Ю.А. Электростатика со столбиком в случае сетки со сгущением достаточно 2.0.

	integer nda = 0; // память под вектор значений матрицы слау.
	nda = static_cast<integer>(rsize*(3 * (nna)+5 * (nnu)));
	integer ndia = 0;
	ndia = static_cast<integer>(rsize*2.2*(nnu));
	integer ndja = 0;
	ndja = static_cast<integer>(rsize*(3 * (nna)+5 * (nnu)));
	integer ndu = 0;
	ndu = static_cast<integer>(rsize*2.2*(nnu));
	integer ndf = 0;
	ndf = static_cast<integer>(rsize*2.2*(nnu));
	integer ndig = 0;
	ndig = static_cast<integer>(rsize*5.4*(nnu));

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



	// рекомедуемые параметры по дефолту.

	integer iswtch = 0;
	iswtch = 4;
	integer iout = 0;
	iout = 13; // 13 обеспечивает печать изменения невязки в процессе счёта.
	integer iprint = 0;
	iprint = 10606;
	integer levelx = 0;
	levelx = 25;
	integer ifirst = 0;
	// начальное приближение:
	// 0 - используется из вне.
	// 1 - нулевое.
	// 2 - единицы.
	// 3 - случайная последовательность.
	ifirst = 13;//13 по умолчанию.
				//ifirst=11; // нулевое начальное приближение.
				//ifirst=10; // вроде как начальное приближение берётся из dX0.
				// но 10 никоим образом не улучшает сходимость.
	integer ncyc = 0;
	ncyc = 10110;
	integer madapt = 0;
	madapt = 27;
	integer nrd = 0;
	nrd = 1131;
	integer nsolco = 0;
	nsolco = 110;
	integer nru = 0;
	nru = 1131;
	doublereal ecg1 = 0.0;
	ecg1 = 0.0;
	doublereal ecg2 = 0.0;
	ecg2 = 0.25;
	doublereal ewt2 = 0.0;
	ewt2 = 0.35;
	integer nwt = 0;
	nwt = 2;
	integer ntr = 0;
	ntr = 0;

	integer matrix = 0;
	matrix=11; // symmetric SPD.
	//matrix = 22;
	ncyc = 10199;


	// allocate memory.
	doublereal *a = nullptr;
	//a=new doublereal[nda+1];
	// 15 jan 2016
	a = (doublereal*)malloc((static_cast<integer>(nda)+1) * sizeof(doublereal));
	if (a == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for a matrix in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	integer *ia = nullptr;
	//ia=new integer[ndia+1];
	ia = (integer*)malloc((static_cast<integer>(ndia)+1) * sizeof(integer));
	if (ia == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for ia matrix in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	integer *ja = nullptr;
	//ja=new integer[ndja+1];
	ja = (integer*)malloc((static_cast<integer>(ndja)+1) * sizeof(integer));
	if (ja == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for ja matrix in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	doublereal *u = nullptr;
	//u = new doublereal[ndu + 1];
	u = (doublereal*)malloc((static_cast<integer>(ndu)+1) * sizeof(doublereal));
	if (u == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for u vector in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	doublereal *f = nullptr;
	//f=new doublereal[ndf+1];
	f = (doublereal*)malloc((static_cast<integer>(ndf)+1) * sizeof(doublereal));
	if (f == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for f vector in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}
	integer *ig = nullptr;
	//ig=new integer[ndig+1];
	ig = (integer*)malloc((static_cast<integer>(ndig)+1) * sizeof(integer));
	if (ig == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for ig vector in amg1r5 algorithm...\n");
		printf("Please any key to exit...\n");
		//system("pause");
		system("pause");
		exit(1);
	}

	// Блок инициализации нулём, возможно будет работоспособно и без него.

	for (integer k = 0; k <= nda; ++k) {
		a[k] = 0.0;
	}
	for (integer k = 0; k <= ndia; ++k) {
		ia[k] = 0;
	}
	for (integer k = 0; k <= ndja; ++k) {
		ja[k] = 0;
	}
	for (integer k = 0; k <= ndu; ++k) {
		u[k] = 0.0;
	}
	for (integer k = 0; k <= ndf; ++k) {
		f[k] = 0.0;
	}
	for (integer k = 0; k <= ndig; ++k) {
		ig[k] = 0;
	}


	// обязателная инициализация.
	for (integer k = 0; k <= nnu + 1; ++k) ia[k + id] = nna + 1; // инициализация.
	if (id == 1) {
		ia[nnu + 2] = 0;
	}






	// начальное приближение.
	for (integer i = 0; i <= ndu; ++i) {
		u[i] = 0.0;
		if (i<n) {
			// обязательно нужно проверить была ли выделена оперативная память. 
			u[i + id] = dX0[i];
		}
	}

	// правая часть.
	for (integer i = 0; i <= ndf; ++i) {
		f[i] = 0.0;
		if (i<n) {
			// обязательно нужно проверить была ли выделена оперативная память. 
			f[i + id] = dV[i];
		}
	}

	// см. equation3DtoCRS.

	integer ik = 0; // счётчик ненулевых элементов СЛАУ

					// для внутренних узлов расчётной области:
	for (integer k = 0; k < n; ++k) {

		integer idiagonal_first_ik = ik;

		//сканируем строку.
		for (integer k1 = m.row_ptr[k]; k1 < m.row_ptr[k + 1]; k1++) {

			if (fabs(m.val[k1]) > nonzeroEPS) {
				if (m.col_ind[k1] != k) {
					// Внедиагональный элемент
					a[ik+1 + id] = m.val[k1];
					ja[ik+1 + id] = m.col_ind[k1] + 1;
					ia[k + id] = my_imin(ik + 1+1, ia[k + id]);
					ik++;
				}
				else {
					// диагональный элемент
					a[idiagonal_first_ik + id] = m.val[k1];
					ja[idiagonal_first_ik + id] = m.col_ind[k1] + 1;
					ia[k + id] = my_imin(idiagonal_first_ik + 1, ia[k + id]);
					//ik++;
				}
			}
		}
		ik++;

	}

	/*
	for (integer k = 0; k < n; ++k) {
		//сканируем строку.
		for (integer k1 = m.row_ptr[k]; k1 < m.row_ptr[k + 1]; k1++) {
			printf("val=%e col_ind=%d row_ptr=%d\n",a[k1+id],ja[k1+id],ia[k]);
		}
		system("pause");
	}
	*/

	// в каждой строке элементы отсортированы по номерам столбцов:
	// Но диагональный элемент всегда на первом месте в строке матрицы.
	integer imove = 0;
	if (id == 0) {
		imove = -1;
	}


	// сортировка ненужна порядок следования любой, но главное чтобы первый в строке был имено диагональный элемент.
	//for (integer k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR_amg(ja, a, ia[k+1]+1+imove, ia[k+2]-1+imove); // первый элемент всегда диагональный.
	//for (integer k=0; k<(maxelm+maxbound); ++k) QuickSortCSIR_amg(ja, a, ia[k+1]+imove, ia[k+2]-1+imove); 

	for (integer k = 1; k <= nnu; ++k) ig[k + imove] = ia[k + 1 + imove]; // инициализация.



	//printf("getready ...");
	//system("pause");

	// amg - особенно хорош для поправки давления в SIMPLE алгоритме.
	// алгоритм 1985 года.
	amg1r5_(a, ia, ja,
		u, f, ig, &nda, &ndia,
		&ndja, &ndu, &ndf, &ndig,
		&nnu, &matrix, &iswtch, &iout,
		&iprint, &levelx, &ifirst, &ncyc,
		&eps, &madapt, &nrd, &nsolco,
		&nru, &ecg1, &ecg2, &ewt2,
		&nwt, &ntr, &ierr);

	switch (ierr) {
	case 1: printf("dimension A small\n.");
		//system("pause");
		system("pause");
		break;
	case 2: printf("dimension IA small\n.");
		//system("pause");
		system("pause");
		break;
	case 3: printf("dimension JA small\n.");
		//system("pause");
		system("pause");
		break;
	case 4: printf("dimension U small\n.");
		//system("pause");
		system("pause");
		break;
	case 5: printf("dimension F small\n.");
		//system("pause");
		system("pause");
		break;
	case 6: printf("dimension IG small\n.");
		//system("pause");
		system("pause");
		break;
	}

	// возвращаем решение СЛАУ.
	for (integer i = 0; i<n; ++i) {
		// обратное копирование.
		dX0[i] = u[i + 1 + imove];
	}


	// освобождение памяти.
	if (a != nullptr) {
		// delete[] a;
		free(a);
		a = nullptr;
	}
	if (ia != nullptr) {
		// delete[] ia;
		free(ia);
		ia = nullptr;
	}
	if (ja != nullptr) {
		//delete[] ja;
		free(ja);
		ja = nullptr;
	}
	if (u != nullptr) {
		//delete[] u;
		free(u);
		u = nullptr;
	}
	if (f != nullptr) {
		//delete[] f;
		free(f);
		f = nullptr;
	}
	if (ig != nullptr) {
		//delete[] ig;
		free(ig);
		ig = nullptr;
	}

	calculation_main_end_time = clock();
	calculation_vorst_seach_time += calculation_main_end_time - calculation_main_start_time;

}




  // этот метод показывает значительно более лучшую сходимость, чем простой BiCGStabCRS,
  // а также он гораздо лучше (и повидимому правильней) чем Bi_CGStab_internal1.
  // Bi_CGStab_internal3 использует предобуславливание из библиотеки Ю.Саада.
  // дата написания Bi_CGStab_internal3: 31.03.2013. 
  // раскурочил 9 августа 2015.
  // 26 сентября 2016 Теперь метод пригоден и для АЛИС сетки.
// Stress
void Bi_CGStab_internal4(SIMPLESPARSE &sparseM,	integer n,
	doublereal *dV, doublereal* &dX0, integer maxit, 
	bool bprintmessage,  QuickMemVorst& m, WALL* &w, int &lw,
	bool* &bondary, unsigned char iVar)
{

	

	// inumiter - номер глобальной итерации (например номер итерации в стационарном алгоритме SIMPLE).
	// параметр inumiter - введён для того чтобы использоваться при отладке, когда нужно посмотреть
	// алгоритм решения СЛАУ на глобальной итерации (алгоритма SIMPLE) с номером большим чем inumiter.

	//bool bexporttecplot = false; // экспорт в tecplot делается лишь в случае проблем со сходимостью.
	bool brc = false;
	//bool bpam_gsp = false; // предобуславливание с помощью нескольких итераций метода Гаусса-Зейделя.

	// размерность квадратной матрицы.
	// n


	// Разреженная матрица СЛАУ
	// в CRS формате.

	if (bStableVersion)
	{

		simplesparsetoCRS(sparseM, m.val, m.col_ind, m.row_ptr, n); // преобразование матрицы из одного формата хранения в другой.
		//m.ballocCRScfd = true;
		simplesparsefree(sparseM, n);

	}
	else {
		// 02.11.2021
		if (m.val != nullptr) {
			delete[] m.val;
			m.val = nullptr;
		}
		if (m.col_ind != nullptr) {
			delete[] m.col_ind;
			m.col_ind = nullptr;
		}
		if (m.row_ptr != nullptr) {
			delete[] m.row_ptr;
			m.row_ptr = nullptr;
		}


		ell_to_CRS(m.val, m.col_ind, m.row_ptr, n); // преобразование матрицы из одного формата хранения в другой.

		std::cout << "n=" << n << "!= ell_global_number_of_strings_26_11_2021=" << ell_global_number_of_strings_26_11_2021 << std::endl;
		//system("pause");

		for (int i = 0; i < ell_global_number_of_strings_26_11_2021; ++i) {
			delete[] data_ell[i];
			delete[] coll_ell[i];
		}
		delete[] data_ell;
		delete[] coll_ell;
		delete[] coll_size;

		data_ell = nullptr;
		coll_ell = nullptr;
	}


	if (dX0 == nullptr) {
		dX0 = new doublereal[n];

#pragma omp parallel for shared(m, dX0) schedule (guided)
		for (integer i38 = 0; i38 < n; i38++) {
			dX0[i38] = 0.0;
		}

	}

	if (iVar == TEMP) {
		doublereal alpharelax = 1.0;

		if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE) {

			// 02.12.2021
			// Ахтунг, ахтунг !!!
			// При нестационарном моделировании ни в коем случае нельзя присваивать
			// alpharelax = 0.99999; это искажает нестационарный член уравнения теплопередачи.


			//std::cout << "unsteady thermal calculation" << std::endl;
			//system("pause");


			

		}
		else {

			// Это не специальная нелинейная версия кода amg1r5 CAMG.
			for (integer k = 0; k < lw; ++k) {
				if ((w[k].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
					(w[k].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY)) {
					alpharelax = 0.99999; // Для того чтобы решение СЛАУ сходилось.
					// 0.9999 - недостаточное значение, температуры не те получаются.
				}
			}

			if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
				(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
				(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {
				alpharelax = 0.99999; // Для того чтобы решение СЛАУ сходилось.
				// 0.9999 - недостаточное значение, температуры не те получаются.
			}

		}

		for (integer i = 0; i < n; ++i) {

			for (integer j = m.row_ptr[i]; j < m.row_ptr[i + 1]; ++j) {
				if (m.col_ind[j] == i) {
					// диагональ.
					if (((bondary != nullptr) && (!bondary[i])) &&
						(m.row_ptr[i + 1] > m.row_ptr[i] + 1)) {
						// Релаксация к предыдущему значению.

						dV[i] += (1.0 - alpharelax) * m.val[j] * dX0[i] / alpharelax;

						m.val[j] = (m.val[j] / alpharelax);
					}

					//debug
					//std::cout << "col_buffer[j]=" << col_ind[j] << std::endl;
					//system("pause");
				}
			}
		}
	}
	// debug
	//for (integer ii = m.row_ptr[650]; ii <= m.row_ptr[651] - 1; ++ii) {
		//printf("a[%d][%d]=%e\n",650, m.col_ind[ii],m.val[ii]);
	//}
	//for (integer ii = m.row_ptr[651]; ii <= m.row_ptr[652] - 1; ++ii) {
		//printf("a[%d][%d]=%e\n", 651, m.col_ind[ii], m.val[ii]);
	//}
	//system("pause");

	// преобразование из SIMPLESPARSE формата в CRS формат хранения.
	//simplesparsetoCRS(M, val, col_ind, row_ptr, n);


	const integer ILU0 = 0;
	const integer ILU_lfil = 1;

	const integer itype_ilu = ILU_lfil;//ILU_lfil;

	// Исходная матрица.
	// m.a=new doublereal[7*n+2]; // CRS
	// m.ja=new integer[7*n+2];
	// 26 сентября 2016.
	m.a = new doublereal[m.row_ptr[n] +  n + 2];
	m.ja = new integer[m.row_ptr[n] + n + 2];
	m.ia = new integer[n + 2];
	// Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
	if ((m.a == nullptr) || (m.ja == nullptr) || (m.ia == nullptr)) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

								
	
	// Эти структуры данных используются в библиотеке Ю.Саада SPARSKIT2.
	// Основной момент который следует уяснить это:
	// мы используем индексацию массивов в СИ начиная с нуля и до size не включая size;
	// В библиотеке Sparskit2 нумерация элементов массива начинается с единицы и до size включая size.
	// Главное что следует понять не следует пытаться переписать SPARSKIT2 так чтобы она соответствовала нумерации с нуля.
	// оставим нумерацию с единицы иначе мы можем запутаться в логике отлаженного кода SPARSKIT2.
	// Но мы в проекте AliceFlowv0_07 работаем с нумерации начинающейся с нуля. 
	// Код Sparskit2 работает только для предобуславливания. Иными словами схема следующая:
	// На входе матрица в CRS формате с нумерацией с нуля. Преобразуем её в матрицу в которой элементы нумеруются с единицы,
	// для этого очевидно нужно к значениям элементов  col_ind и row_ptr прибавить единицу. Индексация же в матрице (преобразованной пусть начинается с нуля)
	// для этого мы внутри кода SPARSKIT2 декременируем ссылку на массив тогда элемент с номером ноль станет первым тоесть будет иметь номер ноль. Потом в конце использования
	// кода SPARSKIT2 мы увеливаем ссылки на массивы на единицу (откат назад). Получая на вход данную преобразованную матрицу CRS формата, мы делаем с помощью неё матрицу предобуславливателя
	// в MSR формате, как обычно используя код SPARSKIT2 без каких либо изменений из вне (он получен с помощью f2c.exe). Используя матрицу предобуславливания в MSR формате
	// мы пихаем её в функию lusol_ из SPARSKIT2 и получаем необходимый вектор x: (LU)x=y; по вектору y и матрице LU в формате MSR. Код в lusol_ содержит преобразование декремента
	// указателей x, y что позволяет использовать обычные векторы x , y в которых нумерация начинается с нуля. В общем x,y менять не надо, они такие же как обычно в AliceFlowv0_07.
	// В lusol_ указатели сначала декримируются --a; а в конце инкремируются ++a; Так что можно использовать код Sparskit2 без изменений !!!



	if (bprintmessage) {
		printf("Incoplete LU Decomposition begin...\n");
	}


	integer ierr = 0;
	
#pragma omp parallel for
		for (integer i = 0; i<m.row_ptr[n]; ++i) 
		{
			m.a[i] = m.val[i];
			m.ja[i] = m.col_ind[i] + 1;
		}
#pragma omp parallel for
		for (integer i = 0; i<n + 1; ++i)
		{
			m.ia[i] = m.row_ptr[i] + 1;
		}
	

	
			m.ri = new doublereal[n]; m.roc = new doublereal[n]; m.s = new doublereal[n]; m.t = new doublereal[n]; m.vec = new doublereal[n];
			m.vi = new doublereal[n]; m.pi = new doublereal[n]; m.dx = new doublereal[n]; m.dax = new doublereal[n];
			m.y = new doublereal[n]; m.z = new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
			if ((m.ri == nullptr) || (m.roc == nullptr) || (m.s == nullptr) || (m.t == nullptr) || (m.vi == nullptr) || (m.pi == nullptr) || (m.dx == nullptr) || (m.dax == nullptr) || (m.y == nullptr) || (m.z == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
	
	

if (itype_ilu == ILU0) 
			{

		

			
				//m.alu=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.jlu=new integer[7*n+2];
				// 26 сентября 2016.
				m.alu = new doublereal[m.row_ptr[n] + n + 2];
				m.jlu = new integer[m.row_ptr[n] + n + 2];

				m.ju = new integer[n + 2];
				
				//m.alurc=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.jlurc=new integer[7*n+2];
				// 26 сентября 2016.
				m.alurc = new doublereal[m.row_ptr[n] + n + 2];
				m.jlurc = new integer[m.row_ptr[n] + n + 2];

				m.jurc = new integer[n + 2];
				m.iw = new integer[n + 2]; // рабочий массив.
				m.ballocCRScfd = true; // память выделена.

				if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.ju == nullptr) || (m.iw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
				//if ((m.alu1 == nullptr) || (m.jlu1 == nullptr) || (m.ju1 == nullptr) || (m.x1 == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					//printf("Problem: not enough memory on your equipment...\n");
					//printf("Please any key to exit...\n");
					//exit(1);
				//}

			
		

		
			ilu0_(n, m.a, m.ja, m.ia, m.alu, m.jlu, m.ju, m.iw, ierr);
			/* if (ibackregulationgl!=nullptr) {
			for (integer i87=0; i87<7*n+2; i87++) {
			m.alu1[i87]= m.alu[i87];
			m.jlu1[i87]=m.jlu[i87];
			}
			for (integer i87=0; i87<n+2; i87++) {
			m.ju1[i87]=m.ju[i87];
			}
			}*/
		

		if (ierr>0) {
#if doubleintprecision == 1
			printf("%lld string in matrix is zero diagonal element...\n", ierr - 1);
#else
			printf("%d string in matrix is zero diagonal element...\n", ierr - 1);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}
	}


if (itype_ilu == ILU_lfil) 
	{

		//bool btemp_quick = m.ballocCRSt;

		// 21 диагональ вместо 7 в расчётах напряженно-деформированного состояния.
		integer lfil = 6; // 2 уровня (0, 1, 2)
		lfil = my_amg_manager.lfil;// 13.10.2018


				// инициализация.
		m.alu = nullptr;
		m.jlu = nullptr;
		m.ju = nullptr;
		m.alu1 = nullptr;
		m.jlu1 = nullptr;
		m.ju1 = nullptr;
		m.x1 = nullptr;
		m.alurc = nullptr;
		m.jlurc = nullptr;
		m.jurc = nullptr;
		m.levs = nullptr;
		m.w = nullptr;
		m.jw = nullptr;
		m.w_dubl = nullptr;
		m.jw_dubl = nullptr;

		//m.iwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
		// 26 сентября 2016.
#if doubleintprecision == 0
		uint64_t size_now= (lfil + 12) * (m.row_ptr[n] + 26 * n + 2) + 4 * n;
		if (size_now > 2147463646ULL) {
			// 17.04.2021
			std::cout << "BiCGStab_internal4 memory size nnz limit reached!!!\n";
			std::cout << "nnz limit ratio = " << ((1.0 * size_now) / (1.0*2147463646ULL)) << std::endl;
			//system("PAUSE");
			//exit(1);
		}
#endif

		m.iwk = (int64_t)((lfil + 12))*((int64_t)((m.row_ptr[n] + 26*n + 2))) + 4 * n; // размерность памяти под матрицу предобуславливания.

		std::cout << m.iwk + 1 << std::endl;
		//system("pause");

		m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
		m.jlu = new integer[m.iwk + 2];
		m.ju = new integer[n + 2];

		m.alurc = new doublereal[m.iwk + 2]; // +2 запас по памяти.
		m.jlurc = new integer[m.iwk + 2];
		m.jurc = new integer[n + 2];
		m.levs = new integer[m.iwk + 2]; // уровень.
		m.w = new doublereal[n + 2]; // +2 запас по памяти.
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			m.jw = new integer[25 * n + 2]; // +2 запас по памяти.
			m.jw_dubl = new integer[25 * n + 2]; // +2 запас по памяти.
		}
		else {
			m.jw = new integer[25 * n + 2]; // +2 запас по памяти.
			m.jw_dubl = new integer[25 * n + 2]; // +2 запас по памяти.
		}
		
		m.w_dubl = new doublereal[n + 2]; // +2 запас по памяти.
		m.ballocCRScfd = true; // память выделена.

		if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.levs == nullptr) || (m.ju == nullptr) || (m.w == nullptr) || (m.jw == nullptr) || (m.w_dubl == nullptr) || (m.jw_dubl == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}





		// iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
		// recomended
		//iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);
		iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);

		if ((ierr == -2) || (ierr == -3)) {

			integer ipassage = 1; // 4 января 2016.
			do {
				printf("\nPlease WAIT... ... ...ierr=");
				std::cout << ierr << std::endl;

				// задаче не хватило памяти, значит нужно перевыделить !
				if (m.alu != nullptr) delete m.alu;
				if (m.jlu != nullptr) delete m.jlu;
				/* if (ibackregulationgl!=nullptr) {
				if (m.alu1!=nullptr) delete m.alu1;
				if (m.jlu1!=nullptr) delete m.jlu1;
				}*/
				if (m.alurc != nullptr) delete m.alurc;
				if (m.jlurc != nullptr) delete m.jlurc;
				if (m.levs != nullptr) delete m.levs;

				// инициализация !
				m.alu = nullptr;
				m.jlu = nullptr;
				/* if (ibackregulationgl!=nullptr) {
				m.alu1=nullptr;
				m.jlu1=nullptr;
				}*/
				m.levs = nullptr;

				// инициализация !
				m.alurc = nullptr;
				m.jlurc = nullptr;

				//m.iwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
				// 26 сентября 2016.
				if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
				{
					m.iwk = ((int64_t)(lfil + 11)) * ((int64_t)(m.row_ptr[n] + 26 * n + 2)) + (int64_t)(((1 + 3 + 3 * ipassage) * n));
				}
				else {
					m.iwk = ((int64_t)(lfil + 11)) * ((int64_t)(m.row_ptr[n] + 26 * n + 2)) + (int64_t)(((1 + 3 + 3 * ipassage) * n));
				}

				m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
				m.jlu = new integer[m.iwk + 2];
				/* (ibackregulationgl!=nullptr) {
				m.alu1=new doublereal[m.iwk+2]; // +2 запас по памяти.
				m.jlu1=new integer[m.iwk+2];
				}*/
				m.levs = new integer[m.iwk + 2]; // уровень.

				if ((m.alu != nullptr) && (m.jlu != nullptr) && (m.levs != nullptr)) {
					// iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
					// рекомендуется
					//iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);
					iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
					/*
					if (ibackregulationgl!=nullptr) {
					for (integer i87=0; i87<m.iwk+2; i87++) {
					m.alu1[i87]= m.alu[i87];
					m.jlu1[i87]=m.jlu[i87];
					}
					for (integer i87=0; i87<n+2; i87++) {
					m.ju1[i87]=m.ju[i87];
					}
					}*/

				}
				else {
					// недостаточно памяти на данном оборудовании.
					ipassage = 4;
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);

				}

				ipassage++;
			} while ((ierr != 0) && (ipassage < 4));

			if (ipassage == 4) {
				printf("Error memory alloc !!!\n");
				printf("failed to obtain an expansion for the 4 approaches...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}
		else {
			/*
			if (ibackregulationgl!=nullptr) {
			for (integer i87=0; i87<m.iwk+2; i87++) {
			m.alu1[i87]= m.alu[i87];
			m.jlu1[i87]=m.jlu[i87];
			}
			for (integer i87=0; i87<n+2; i87++) {
			m.ju1[i87]=m.ju[i87];
			}
			}*/
		}





		if (ierr != 0) {
#if doubleintprecision == 1
			printf("error memory in iluk ierr=%lld\n", ierr);
#else
			printf("error memory in iluk ierr=%d\n", ierr);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}


		if (bprintmessage) {
			printf("Incoplete LU Decomposition finish...\n");
		}
	}


	//bool bnorelax = true; // Для уравнения теплопроводности не используется релаксация.


	//bool bprintf = false; // если bprintf==false то значения невязок внутри LR1sk не выводятся.
	integer iflag = 1, icount = 0;
	doublereal delta0 = 1.0e30, deltai = 1.0e30;
	doublereal bet = 0.0, roi = 0.0;
	doublereal roim1 = 1.0, al = 1.0, wi = 1.0;
	//doublereal *ri, *roc, *s, *t, *vi, *pi, *dx, *dax;
	//doublereal *y, *z; // результат предобуславливания
	doublereal epsilon = dterminatedTResudual;  // точность вычисления
	
	epsilon *= 1.0e-4; // 1.0e-4
	
	



#pragma omp parallel for shared(m)  schedule (guided)
	for (integer i = 0; i<n; ++i) {
		
			m.s[i] = 0.0;
			m.t[i] = 0.0;
			m.vi[i] = 0.0;
			m.pi[i] = 0.0;
			// инициализатор массивов для предобуславливания
			m.y[i] = 0.0;
			m.z[i] = 0.0;
			// результат умножения матрицы на вектор.
			m.dax[i] = 0.0;
		
	}

	// начальное приближение
	// X0 ==
	// под X0 понимается вектор поля температур к примеру.
	if (dX0 == nullptr) {
		dX0 = new doublereal[n];
		
			
#pragma omp parallel for shared(m, dX0)  schedule (guided)
				for (integer i = 0; i<n; ++i) {
					m.dx[i] = 0.0;
					dX0[i] = 0.0;
				}		

	}
	else {
		
			
#pragma omp parallel for shared(m, dX0)  schedule (guided)
				for (integer i = 0; i<n; ++i) m.dx[i] = dX0[i];
			
		

	}

	
		MatrixCRSByVector(m.val, m.col_ind, m.row_ptr, m.dx, m.dax, n); // результат занесён в  dax
	


#pragma omp parallel for shared(dV,m)  schedule (guided)
	for (integer i = 0; i<n; ++i) {

				m.ri[i] = dV[i] - m.dax[i];
				//m.roc[i]=m.ri[i];
				m.roc[i] = 1.0;
			
		
	}
	
		delta0 = NormaV(m.ri, n);
	

	//printf("debug %e\n",NormaV(dax,n)); // проверка на корректное составление СЛАУ
	//printf("%e \n",delta0); system("pause");
	//system("pause");
	// Если решение сразу хорошее то не считать:
	
	if (fabs(delta0)<dterminatedTResudual) iflag = 0;
	
	integer iflag1 = 1;
	if (fabs(delta0)<1e-23) iflag1 = 0;
	/*
	if ((iVar == TEMP) && (iflag == 0) && (iflag1 == 0)) {
#if doubleintprecision == 1
		printf("iflag=%lld, iflag1=%lld, delta0=%e\n", iflag, iflag1, delta0);
#else
		printf("iflag=%d, iflag1=%d, delta0=%e\n", iflag, iflag1, delta0);
#endif

		//system("pause");
		system("PAUSE");
	}
	*/
#if doubleintprecision == 1
	//printf("iflag=%lld, iflag1=%lld, delta0=%e\n",iflag,iflag1,delta0); system("pause");
#else
	//printf("iflag=%d, iflag1=%d, delta0=%e\n",iflag,iflag1,delta0); system("pause");
#endif


	/*if (iVar==PAM) {
	// терминаьная невязка всегда на точность аппроксимации меньше стартовой невязки.
	//if (e*dnew<e) e*=dnew;
	doublereal me=sqrt(fabs(delta0))/n;
	if (epsilon*me<epsilon) epsilon*=me;
	//dterminatedTResudual=e;
	}*/



	//printf("delta0=%e\n",delta0);
	//system("pause");


	/*if (iflag != 0) {
	// конечная невязка всегда на точность аппроксимации меньше начальной невязки.
	epsilon*=delta0;
	dterminatedTResudual=epsilon;
	}
	*/
	integer iN = 10;
	if (n <= 15000) {
		// задача очень малой размерности !
		
			iN = 1; // обязательно нужна хотя бы одна итерация.
					// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
					//printf("%e\n",epsilon);
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon = 1.0e-3*fabs(delta0);
				//printf("%e\n", epsilon);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
		
	}
	else if ((n>15000) && (n<30000)) {
		// задача очень малой размерности !
		
			iN = 1; // обязательно нужна хотя бы одна итерация.
					// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon = 1.0e-3*fabs(delta0);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
		
	}
	else if ((n >= 30000) && (n<100000)) {
		// Здесь я немного увеличил число итераций и 
		// скорректировал условие окончания чтобы считало 
		// поточнее, но это не повлияло.
		// Главный вопрос в том что невязка по температуре почему-то не меняется.
		// задача небольшой размерности.
		
			iN = 3; // обязательно нужна хотя бы одна итерация.
					// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon = 1.0e-3*fabs(delta0);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
			// 27.07.2016
			iN = 12;
			epsilon *= 1e-2;
		
	}
	else if ((n >= 100000) && (n<300000)) {
		// задача небольшой средней размерности.
		
			iN = 3; // обязательно нужна хотя бы одна итерация.
					// Вообще говоря невязка для скоростей падает очень быстро поэтому всегда достаточно iN итераций для скорости.
					// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon = 1.0e-3*fabs(delta0);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
		
	}
	else if ((n >= 300000) && (n<1000000)) {
		// задача истинно средней размерности.
		
			iN = 3; // обязательно нужна хотя бы одна итерация.
					// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon = 1.0e-3*fabs(delta0);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
		
	}
	else if ((n >= 1000000) && (n<3000000)) {
		// задача достаточно большой размерности.
		
			iN = 6; // обязательно нужна хотя бы одна итерация.
					// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon = 1.0e-3*fabs(delta0);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
		
	}
	else if (n >= 3000000) {
		// задача очень большой размерности.
		
			iN = 6; // обязательно нужна хотя бы одна итерация.
					// если этого будет недостаточно то мы всё равно будем итерировать до тех пор пока невязка не станет меньше epsilon.
			if (1.0e-3*fabs(delta0)<epsilon) {
				epsilon = 1.0e-3*fabs(delta0);
			}
			if (iflag1 == 1) {
				iflag = 1;
			}
		
	}

	epsilon *= 1.0e-7;

	epsilon = 1.0e-10; // По рекомендации Алексея Ивановича. 20.02.2021
	printf("epsilon=%e \n",epsilon);
	
		//maxit = 1000;//2000
		maxit = 2500;


	/*if (iVar==PAM) {
	printf(" %1.2e ",delta0);
	#if doubleintprecision == 1
	printf("icount=%lld iN=%lld iflag1=%lld iflag=%lld maxit=%lld",icount,iN,iflag1,iflag,maxit);
	#else
	printf("icount=%d iN=%d iflag1=%d iflag=%d maxit=%d",icount,iN,iflag1,iflag,maxit);
	#endif

	system("pause");
	}*/

	/*
	if (iVar==TEMP) {
	//BiSoprGradCRS( m.tval, m.tcol_ind, m.trow_ptr,dV,dX0,n,200);
	Bi_CGStabCRS(n, m.tval, m.tcol_ind, m.trow_ptr, dV, dX0, 200);
	}
	else {
	//BiSoprGradCRS( m.val, m.col_ind, m.row_ptr,dV,dX0,n,200);
	Bi_CGStabCRS(n, m.val, m.col_ind, m.row_ptr, dV, dX0, 200);
	}
	*/

	// диагностическое сообщение какую переменную мы решаем.
	//switch (iVar) {
	//case PAM: printf("PAM\n");  break;
	//case VX:  printf("VX\n"); break;
	//case VY:  printf("VY\n"); break;
	//case VZ:  printf("VZ\n"); break;
	//case TEMP:  printf("TEMP\n"); break;
	//}

	// Если число расходимостей превысит оговорённую константу то произойдёт выход из алгоритма.
	integer i_signal_break_pam_opening = 0;
	// x хорошее значение.
	//const integer i_limit_signal_pam_break_opening = 4000;//20
	doublereal delta_old_iter = 1.0e10;

	integer count_iter_for_film_coef = 0;



	

	// Мы обязательно должны сделать несколько итераций. (не менее 10).
	// Если только решение не удовлетворяет уравнению тождественно.
	while (((icount < iN) && (iflag1 != 0)) || (iflag != 0 && icount < maxit)) {

		icount++;

		count_iter_for_film_coef++;
		// В случае задачи Ньютона - Рихмана, Стефана-Больцмана и миксового условия не итерируем до конца обрываем, 
		// т.к. нам требуется частая пересборка матрицы. 13 марта 2016.
		//if (((adiabatic_vs_heat_transfer_coeff > ADIABATIC_WALL_BC) || (breakRUMBAcalc_for_nonlinear_boundary_condition)) && (count_iter_for_film_coef>5)) break;


		
			roi = Scal(m.roc, m.ri, n);
			if (!isfinite(roi)) {
				printf("roi!=roi solution bug. \n");
				system("pause");
			}
			if (fabs(wi) < 1.0e-30) {
				if (fabs(roim1) < 1.0e-30) {
					bet = 1.0;
				}
				else {
					bet = (roi / roim1);
				}
			}
			else {
				if (fabs(roim1) < 1.0e-30) {
					bet = (al / wi);
				}
				else {
					bet = (roi / roim1)*(al / wi);
				}
			}
			if  (!isfinite(bet)) {
				printf("bet!=bet solution bug. \n");
				printf("%e %e %e %e\n", roi, roim1, al, wi);
				system("pause");
			}


			//bet = (roi / roim1)*(al / wi);

			//printf("%e %e %e %e\n",roi,roim1,al,wi);
			//system("pause");

#pragma omp parallel for shared(m,wi,bet)  schedule (guided)
			for (integer i = 0; i<n; ++i) {
				doublereal pibuf = m.ri[i] + (m.pi[i] - m.vi[i] * wi)*bet;
				m.pi[i] = pibuf;
			}
		

		// Ky=pi

		// (LU)y=pi; 
		
			// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for shared(m)  schedule (guided)
			for (integer i = 0; i<n; ++i) m.y[i] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

				//  9 августа 2015 при внедрении перенумерации узлов nested desection
			

				if (brc) {
#pragma omp parallel for
					for (integer i7 = 0; i7 < n; ++i7)
					{
						m.vec[i7] = m.pi[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7<m.iwk + 2; ++i7)
					{
						m.alurc[i7] = m.alu[i7];
						m.jlurc[i7] = m.jlu[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < n + 2; ++i7)
					{
						m.jurc[i7] = m.ju[i7];
					}
				}

				
					lusol_(n, m.pi, m.y, m.alu, m.jlu, m.ju, n); // M*y=pi;

				

				if (brc) {
#pragma omp parallel for
					for (integer i7 = 0; i7 < n; ++i7)
					{
						m.pi[i7] = m.vec[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7<m.iwk + 2; ++i7) 
					{
						m.alu[i7] = m.alurc[i7];
						m.jlu[i7] = m.jlurc[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < n + 2; ++i7)
					{
						m.ju[i7] = m.jurc[i7];
					}
				}

			

			MatrixCRSByVector(m.val, m.col_ind, m.row_ptr, m.y, m.vi, n); // vi==A*y;
		


		

			if ((fabs(roi)<1e-30) && (fabs(Scal(m.roc, m.vi, n))<1e-30)) {
				al = 1.0;
			}
			else if (fabs(roi)<1e-30) {
				al = 0.0;
			}
			else {
				al = roi / Scal(m.roc, m.vi, n);
			}


#pragma omp parallel for shared(m,al)  schedule (guided)
			for (integer i = 0; i<n; ++i) {
				m.s[i] = m.ri[i] - al*m.vi[i];
			}
		
		// Kz=s

		// (LU)z=s; 
		
			// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for shared(m)  schedule (guided)
			for (integer i = 0; i<n; ++i) m.z[i] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

			

				if (brc) {
#pragma omp parallel for
					for (integer i7 = 0; i7 < n; ++i7)
					{
						m.vec[i7] = m.s[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7<m.iwk + 2; ++i7)
					{
						m.alurc[i7] = m.alu[i7];
						m.jlurc[i7] = m.jlu[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < n + 2; ++i7)
					{
						m.jurc[i7] = m.ju[i7];
					}
				}

				
				lusol_(n, m.s, m.z, m.alu, m.jlu, m.ju, n); // Mz=s;
				

				if (brc) {
#pragma omp parallel for
					for (integer i7 = 0; i7 < n; ++i7)
					{
						m.s[i7] = m.vec[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
						m.alu[i7] = m.alurc[i7];
						m.jlu[i7] = m.jlurc[i7];
					}
#pragma omp parallel for
					for (integer i7 = 0; i7 < n + 2; ++i7)
					{
						m.ju[i7] = m.jurc[i7];
					}
				}
			
			MatrixCRSByVector(m.val, m.col_ind, m.row_ptr, m.z, m.t, n); // t==A*z;
		

		

			wi = Scal(m.t, m.s, n) / Scal(m.t, m.t, n);

			// printf("%e %e",Scal(m.t,m.s,n),Scal(m.t,m.t,n));

#pragma omp parallel for shared(m, al, wi) schedule (guided)
			for (integer i = 0; i<n; ++i) {
				//dx[i]+=al*pi[i]+wi*s[i]; // так было без предобуславливателя
				m.dx[i] += al*m.y[i] + wi*m.z[i]; // так стало с предобуславливателем
				m.ri[i] = m.s[i] - wi*m.t[i];
			}
			deltai = NormaV(m.ri, n);
		
			doublereal max_deformation = -1.0e+30;
			if (iVar == TOTALDEFORMATION) {
				for (integer i = 0; i < n; i=i+3) {
					doublereal total_deformation = sqrt(m.dx[i]* m.dx[i]+ m.dx[i+1] * m.dx[i+1]+ m.dx[i+2] * m.dx[i+2]);
					if (total_deformation > max_deformation) {
						max_deformation = total_deformation;
					}
				}
			}
			else {
				for (integer i = 0; i < n; ++i) {
					if (m.dx[i] > max_deformation) {
						max_deformation = m.dx[i];
					}
				}
			}


		//printf("deltai=%e\n",deltai); system("pause");

		// печать невязки на консоль
		if (bprintmessage) {
			if ((icount % 10) == 0) {
				if (iVar == TEMP) {
					std::cout << "iter  residual maximum_temperature" << std::endl;
				}
				if (iVar == TOTALDEFORMATION) {
					std::cout << "iter  residual maximum_deformation" << std::endl;
				}
			}
			//if ((icount % 1) == 0) {
				std::cout << icount << " " << deltai << " " << max_deformation << std::endl;
			//}

		}
		// 28.07.2016.
		//std::cout << icount << " " << deltai << std::endl;

		//system("pause");
		if (deltai > delta_old_iter) i_signal_break_pam_opening++;
		delta_old_iter = deltai;
		

		if (deltai < epsilon) {
			iflag = 0; // конец вычисления
					   // 20.05.2017
			iflag1 = 0; // Выход по второму критерию чтоб не делать тысячи итераций.
		}
		else roim1 = roi;

		
#if doubleintprecision == 1
			//printf("epsilon=%e deltai=%e icount=%lld\n",epsilon,deltai, icount);
#else
			//printf("epsilon=%e deltai=%e icount=%d\n",epsilon,deltai, icount);
#endif

			//system("pause");
		

	}

	
		if (!((maxit == 0) && (iN == 0))) {
			
#pragma omp parallel for shared(dX0, m)  schedule (guided)
				for (integer i = 0; i<n; ++i) dX0[i] = m.dx[i];
		}
	
	


	// Это матрица в котрой нумерация (а индексация элементов с нуля) начинается с единицы. Она используется в библиотеке SPARSKIT2.
	
			// Это таже CRS матрица что и a,ja, ia только элементы в ней нумеруются также как и индексируются с нуля.
		if (m.val != nullptr) {
			delete[] m.val;
			m.val = nullptr;
		}
		if (m.col_ind != nullptr) {
			delete[] m.col_ind;
			m.col_ind = nullptr;
		}
		if (m.row_ptr != nullptr) {
			delete[] m.row_ptr;
			m.row_ptr = nullptr;
		}
		if (m.a != nullptr) {
			delete[] m.a;
			m.a = nullptr;
		}
		if (m.ja != nullptr) {
			delete[] m.ja;
			m.ja = nullptr;
		}
		if (m.ia != nullptr) {
			delete[] m.ia; // уничтожаем матрицу в CRS формате.
			m.ia = nullptr;
		}
										   // освобождение памяти
		if (m.ri != nullptr) {
			delete[] m.ri;
			m.ri = nullptr;
		}
		if (m.roc != nullptr) {
			delete[] m.roc;
			m.roc = nullptr;
		}
		if (m.s != nullptr) {
			delete[] m.s;
			m.s = nullptr;
		}
		if (m.t != nullptr) {
			delete[] m.t;
			m.t = nullptr;
		}
		if (m.vi != nullptr) {
			delete[] m.vi;
			m.vi = nullptr;
		}
		if (m.pi != nullptr) {
			delete[] m.pi;
			m.pi = nullptr;
		}
		if (m.dax != nullptr) {
			delete[] m.dax;
			m.dax = nullptr;
		}
		if (m.y != nullptr) {
			delete[] m.y;
			m.y = nullptr;
		}
		if (m.z != nullptr) {
			delete[] m.z;
			m.z = nullptr;
		}
		if (m.dx != nullptr) {
			delete[] m.dx;
			m.dx = nullptr;
		}
		if (m.vec != nullptr) {
			delete[] m.vec;
			m.vec = nullptr;
		}
			// alu, jlu - MSR матрица неполного ILU разложения.
			// ju - указатель на диагональные элементы, iw - вспомогательный вектор.
		if (m.alu != nullptr) {
			delete[] m.alu;
			m.alu = nullptr;
		}
		if (m.jlu != nullptr) {
			delete[] m.jlu;
			m.jlu = nullptr;
		}
		if (m.ju != nullptr) {
			delete[] m.ju;
			m.ju = nullptr;
		}
			
		if (m.alurc != nullptr) {
			delete[] m.alurc;
			m.alurc = nullptr;
		}
		if (m.jlurc != nullptr) {
			delete[] m.jlurc;
			m.jlurc = nullptr;
		}
		if (m.jurc != nullptr) {
			delete[] m.jurc;
			m.jurc = nullptr;
		}
			if (itype_ilu == ILU0) {
				if (m.iw != nullptr) {
					delete[] m.iw; // удаляем рабочий массив.
					m.iw = nullptr;
				}
			}
			// Освобождение памяти.
if (itype_ilu == ILU_lfil) 
			{
	if (m.w != nullptr) {
		delete[] m.w;
		m.w = nullptr;
	}
	if (m.jw != nullptr) {
		delete[] m.jw;
		m.jw = nullptr;
	}
	if (m.w_dubl != nullptr) {
		delete[] m.w_dubl;
		m.w_dubl = nullptr;
	}
	if (m.jw_dubl != nullptr) {
		delete[] m.jw_dubl;
		m.jw_dubl = nullptr;
	}
	if (m.levs != nullptr) {
		delete[] m.levs;
		m.levs = nullptr;
	}
			}

			m.bsignalfreeCRScfd = false; // память корректно освобождена.
	

	// Шаманство.
	//if (iVar==VX) {
	// m.icount_vel=icount;
	// }
	//if ((iVar==VY)||(iVar==VZ)) {
	//	 if (icount>m.icount_vel) {
	//	 m.icount_vel=icount;
	//}
	//}

	// отладочная печать количества сделанных итераций.
	/*
	#if doubleintprecision == 1
	printf("%lld ",icount);
	#else
	printf("%d ",icount);
	#endif

	if(iVar==PAM) {
	printf("%e ",deltai/delta0);
	//printf("%e %e\n",deltai,delta0);
	//system("pause");
	}
	*/
	//if (iVar==TEMP) printf(" %e ",deltai/delta0);

} // Bi_CGStab_internal4

// Не работает.
void Direct_for_coef(doublereal** & Xmatr, doublereal* & koefmatr, doublereal* & bmatr, integer nodes) {
	/*IMatrix sparseS; // разреженная матрица в формате IMatrix
	initIMatrix(&sparseS, 10);
	for (int i_91 = 0; i_91 < 10; i_91++) {
		for (int i_92 = 0; i_92 < 10; i_92++) {
			setValueIMatrix(&sparseS, i_91, i_92, Xmatr[i_91][i_92]);
		}
	}
	// главный метод, возвращающий решение x,
	// принимает вектор свободных членов b и 
	// квадратную матрицу xO в специальном разреженном формате.
	// реализация без барьера и итерационного уточнения.
	calculateSPARSEgaussArray(&sparseS, koefmatr, bmatr);

	freeIMatrix(&sparseS);*/
	const doublereal epsilon = 1.0e-14;

	bool *vacant = new bool[nodes];
	integer* inum = new integer[nodes];
	for (integer k1 = 0; k1 < nodes; k1++) {
		vacant[k1] = true;
	}
	for (integer l1 = 0; l1 < nodes; l1++) {
		doublereal dmax = -1.0e30;
		integer m1 = -1;
		for (integer k1 = 0; k1 < nodes; k1++) {
			if (vacant[k1]) {
				if (fabs(Xmatr[k1][l1]) > dmax) {
					dmax = fabs(Xmatr[k1][l1]);
					m1 = k1;
				}
			}
		}
		inum[l1] = m1;
		vacant[m1] = true;
	}

	doublereal** A = new doublereal*[nodes];
	doublereal* b = new doublereal[nodes];
	for (integer k1 = 0; k1 < nodes; k1++) {
		A[k1] = new doublereal[nodes];
		for (integer l1 = 0; l1 < nodes; l1++) {
			A[k1][l1] = Xmatr[inum[k1]][l1];
		}
		b[k1] = bmatr[inum[k1]];
	}


	 { // если матрица ещё не приведена к верхнетреугольному виду
		doublereal M;
		// приведение к верхне треугольному виду:
		for (integer k1 = 0; k1 < nodes; k1++) {
			for (integer i1 = k1 + 1; i1 < nodes; i1++) {
				// Если на диагонали ноль:
				if (fabs(A[k1][k1]) < epsilon) {
					// решение не может быть получено, т.к.
					// на диагонали находится ноль.
					printf("\n inverse matrix simple ERROR !!! may be diagonal value is zero...\n");
					printf("\nSolution is not exist.\n");
					for (integer irow = 0; irow < nodes; irow++) {
						for (integer icol = 0; icol < nodes; icol++) {
							printf("%1.4e ", A[irow][icol]);
							//printf("%1.4e ", A[irow][icol]);

						}
						printf("\n");
					}
					//system("pause");
					system("pause");
					exit(0);
				}
				M = A[i1][k1] / A[k1][k1];
				for (integer j1 = k1; j1 < nodes; j1++) {
					A[i1][j1] -= M * A[k1][j1];
				}
				// преобразование правых частей:
				for (integer j1 = 0; j1 < nodes; j1++) b[j1] -= M * A[k1][j1];
			}
		}
	}
	doublereal sum = 0.0;


	// процесс обратного исключения
	for (integer i1 = nodes - 1; i1 >= 0; i1--) {
		// инициализация
		 sum = 0.0;

		 for (integer j1 = i1 + 1; j1 < nodes; j1++) {
				sum += A[i1][j1] * koefmatr[i1];
			}
		
		
		  koefmatr[i1] = (b[i1] - sum) / A[i1][i1];
		
	}

	for (integer j1 = 0; j1 < nodes; j1++) {
		b[j1] = koefmatr[j1];
	}
	for (integer j1 = 0; j1 < nodes; j1++) {
		koefmatr[j1]=b[inum[j1]];
	}

	delete[] vacant;
	delete[] inum;
	delete[] b;
	for (integer j1 = 0; j1 < nodes; j1++) {
		delete[] A[j1];
	}
	delete[] A;

}

// прямой метод для задач чистой теплопроводности.
void Direct(equation3D* &sl, equation3D_bon* &slb,
			   integer maxelm, integer maxbound,
			   doublereal *dV, doublereal* &dX0)
{
	IMatrix sparseS; // разреженная матрица в формате IMatrix
	initIMatrix(&sparseS, maxelm + maxbound);
 
    for (integer i=0; i<maxelm; ++i) {
        setValueIMatrix(&sparseS,sl[i].iP,sl[i].iP,sl[i].ap);
        const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

			
	    if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)){
               setValueIMatrix(&sparseS,sl[i].iP,sl[i].iE,-sl[i].ae);
		}
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) {
		       setValueIMatrix(&sparseS,sl[i].iP,sl[i].iN,-sl[i].an);
		}
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) {
               setValueIMatrix(&sparseS,sl[i].iP,sl[i].iT,-sl[i].at);
		}
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) {
               setValueIMatrix(&sparseS,sl[i].iP,sl[i].iS,-sl[i].as);
		}
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) {
               setValueIMatrix(&sparseS,sl[i].iP,sl[i].iW,-sl[i].aw);
		}
		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) {
               setValueIMatrix(&sparseS,sl[i].iP,sl[i].iB,-sl[i].ab);
		}	
    }

    // Запись уравнений для граничных узлов в матрицу:
    for (integer i=0; i<maxbound; ++i) {
         setValueIMatrix(&sparseS,slb[i].iW ,slb[i].iW, slb[i].aw);
						 
	     const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

	     if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) {
	         setValueIMatrix(&sparseS,slb[i].iW, slb[i].iI,-slb[i].ai);
	     }	 
    }


   // главный метод, возвращающий решение x,
   // принимает вектор свободных членов b и 
   // квадратную матрицу xO в специальном разреженном формате.
   // реализация без барьера и итерационного уточнения.
   calculateSPARSEgaussArray(&sparseS, dX0, dV);

	freeIMatrix(&sparseS);

}

/*
#include "hypre\hypre-2.0.0\src\utilities\_hypre_utilities.h"
#include "hypre\hypre-2.0.0\src\krylov\HYPRE_krylov.h"
#include "hypre\hypre-2.0.0\src\HYPRE.h"
#include "hypre\hypre-2.0.0\src\parcsr_ls\_hypre_parcsr_ls.h"

// решение на основе библиотеки hypre
void hypreSolve(equation3D* &sl, equation3D_bon* &slb,
			   integer maxelm, integer maxbound,
			   doublereal *dV, doublereal* &dX0)
{
}
*/


// Здесь содержится обвязка вызывающая amg1r5.
/*void amg(equation3D* , equation3D_bon* ,
			   integer , integer ,
			   doublereal *, doublereal* , integer ,
			   doublereal , integer );*/


// переменная с глобальной областью видимости необходима для 
// защиты от холостого рестарта, (перезапуск на сошедшмся решении).
// Это должно существенным образом экономить время пользователя.
// Нехолостой рестарт необходим на нелинейных задачах.
// 23 июля 2015.
doublereal finish_residual_lr1sk=0.0;

// Возрождение Lr1sk солвера.
// смысл сделать обвязку по новым требованиям.
// 23 июля 2015.
void Lr1sk_up(FLOW &f, TEMPER &t, equation3D* &sl, equation3D_bon* &slb,
			   integer maxelm, integer maxbound,
			   doublereal *dV, doublereal* &dX0, integer maxit, doublereal alpharelax, integer iVar,  bool bLRfree)
{
        // Замер времени.
	    unsigned int calculation_main_start_time; // начало счёта мс.
	    unsigned int calculation_main_end_time; // окончание счёта мс.

	    calculation_main_start_time=clock(); // момент начала счёта.

	integer i; // счётчик.

	 // На случай если память не была выделена.
	 if (dX0==nullptr) {
	    dX0=new doublereal[maxelm+maxbound];	    
	    for (i=0; i<maxelm+maxbound; ++i) {
	        dX0[i]=0.0;
	    }
	 }

	

	const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля
	doublereal res_sum=0.0;
	res_sum=0.0;
	for (i=0; i<maxelm; ++i) {
		// внутренность матрицы.
		doublereal buf=0.0;
		buf=(sl[i].ap*dX0[sl[i].iP]-dV[sl[i].iP]);
		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) buf-=sl[i].ab*dX0[sl[i].iB];
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) buf-=sl[i].ae*dX0[sl[i].iE];
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) buf-=sl[i].an*dX0[sl[i].iN];
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) buf-=sl[i].as*dX0[sl[i].iS];
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) buf-=sl[i].at*dX0[sl[i].iT];
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) buf-=sl[i].aw*dX0[sl[i].iW];
		buf*=buf;
		res_sum+=buf;
	}
	for (i=0; i<maxbound; ++i) {
		// граничные узлы.
		doublereal buf=0.0;
		buf=slb[i].aw*dX0[slb[i].iW]-dV[slb[i].iW];
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) buf-=slb[i].ai*dX0[slb[i].iI];
		buf*=buf;
		res_sum+=buf;
	}
	res_sum=sqrt(res_sum);

	
	//printf("residual start=%1.4e\n",res_sum);
	//system("pause");
	
	// результаты тестирования
	// задача, начальная невязка , значение евклидовой нормы невязки при которой решение является полученным.
	// tgf01 5.4357e-1 1.0209e-11
	// CGHV1J с метализацией 3.3667e-1 5.0712e-12
	// tgf02 7.6872e-11 1.434e-11
	// tgf05 1.0871e+0  2.2895e-11
	// резистор на 1мм поликоре 5.0e-2 4.9174e-14
	//Diamond ZUb 4 4.0016e-1  4.64444e-11
	// DiamondZUB 4.0016e-1 1.1443e-8
	// NXP100 4.3399e+0  7.8347e-11 (для решения хватило 8Гб ОЗУ.)

	if (bSIMPLErun_now_for_temperature) {
		// При решении CFD задач нам обязательно нужен перезапуск.
		finish_residual_lr1sk = 0.0;
	}
	//if (res_sum>1.0E-10) 
	if (res_sum>1.05*finish_residual_lr1sk) // защита от повторного холостого запуска экономит время конечного пользователя.
	{
	    // разреженная матрица в формате CRS
        doublereal *val=nullptr;
        integer *col_ind=nullptr, *row_ptr=nullptr;


		{

			// TODO получить val, col_ind, row_ptr
			integer nna=0; // количество ненулевых элементов в матрице СЛАУ.
	

	        // подсчёт числа ненулевых элементов в матрице.
	        nna=0;
	        for (i=0; i<maxelm; ++i) {
		        // внутренность матрицы.
		        if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
		        if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
		        if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
		        if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
		        if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
		        if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;
		        if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;
	        }
	        for (i=0; i<maxbound; ++i) {
		         // граничные узлы.
		         if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
		         if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
	        }

	        integer nnu=0; // число неизвестных.
	        nnu=maxelm+maxbound;

			

			// allocate memory.
	        val=new doublereal[nna];
	        if (val==nullptr) {
	           // недостаточно памяти на данном оборудовании.
		       printf("Problem: not enough memory on your equipment for val matrix in lr1sk_new algorithm...\n");
		       printf("Please any key to exit...\n");
		       //system("pause");
			   system("pause");
		       exit(1);
	        }
	        row_ptr=new integer[nnu+1];
	        if (row_ptr==nullptr) {
	              // недостаточно памяти на данном оборудовании.
		          printf("Problem: not enough memory on your equipment for row_ptr matrix in lr1sk_new algorithm...\n");
		          printf("Please any key to exit...\n");
		          //system("pause");
				  system("pause");
		          exit(1);
	        }
	        col_ind=new integer[nna];
	        if (col_ind==nullptr) {
	            // недостаточно памяти на данном оборудовании.
		        printf("Problem: not enough memory on your equipment for col_ind matrix in lr1sk_new algorithm...\n");
		        printf("Please any key to exit...\n");
		        //system("pause");
				system("pause");
		        exit(1);
        	}


			// Блок инициализации нулём, возможно будет работоспособно и без него.

	       for (integer k=0; k<nna; ++k) {
		       val[k]=0.0;
	       }
	       for (integer k=0; k<=nnu; ++k) {
		       row_ptr[k]=0;
	       }
	       for (integer k=0; k<nna; ++k) {
		       col_ind[k]=0;
	       }

		   // обязателная инициализация.
	       for (integer k=0; k<=nnu; ++k) row_ptr[k]=nna; // инициализация.

		   // см. equation3DtoCRS.

	    integer ik=0; // счётчик ненулевых элементов СЛАУ
		
		// для внутренних узлов расчётной области:
        for (integer k=0; k<maxelm; ++k) {

			if (fabs(sl[k].ap) > nonzeroEPS) {
                val[ik]=sl[k].ap/alpharelax;
				col_ind[ik]=sl[k].iP;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iE>-1) && (fabs(sl[k].ae) > nonzeroEPS)) {
                val[ik]=-sl[k].ae;
				col_ind[ik]=sl[k].iE;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iN>-1) && (fabs(sl[k].an) > nonzeroEPS)) {
                val[ik]=-sl[k].an;
				col_ind[ik]=sl[k].iN;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iT>-1) && (fabs(sl[k].at) > nonzeroEPS)) {
                val[ik]=-sl[k].at;
				col_ind[ik]=sl[k].iT;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}		
			if ((sl[k].iS>-1) && (fabs(sl[k].as) > nonzeroEPS)) {
                val[ik]=-sl[k].as;
				col_ind[ik]=sl[k].iS;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iW>-1) && (fabs(sl[k].aw) > nonzeroEPS)) {
				val[ik]=-sl[k].aw;
				col_ind[ik]=sl[k].iW;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}
			if ((sl[k].iB>-1) && (fabs(sl[k].ab) > nonzeroEPS)) {
				val[ik]=-sl[k].ab;
				col_ind[ik]=sl[k].iB;
                row_ptr[k]= myi_min(ik,row_ptr[k]);
				ik++;
			}


		}


		// для внутренних узлов расчётной области:
        for (integer k=0; k<maxbound; ++k) {
			if (fabs(slb[k].aw) > nonzeroEPS) {
				if (ik<nna) {
					// val[ik]=slb[k].aw/alpharelax;
					val[ik] = slb[k].aw; // релаксация для граничных узлов не применяется.
					/*if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
						 // Внимание !!! было произведено тестирование: один вариант был с нижней релаксацией для граничных узлов,
						 // а второй вариант был без нижней релаксации на граничных узлах. Было выяснено, что для сходимости
						 // более благоприятен вариант без нижней релаксации на граничных узлах.
						 // Данное изменение согласовано с функцией solve.

						 val[ik]/=alpharelax; // Если условия Неймана то нижняя релаксация.
					}*/
					col_ind[ik] = slb[k].iW;
					row_ptr[maxelm + k] = myi_min(ik, row_ptr[maxelm + k]);
					ik++;
				}
			}
			if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				val[ik]=-slb[k].ai;
				col_ind[ik]=slb[k].iI;
                row_ptr[maxelm+k]= myi_min(ik,row_ptr[maxelm+k]);
				// Это очень важный вопрос и он требует проверки !
				
				ik++;
			}

		}

		

		// А.А. Фомин, Л.Н. Фомина
		// Ускорение полилинейного рекуррентного метода в подпространствах крылова.
        // Вестник томского государственного университета. Математика и механика №2(14) 2011год.
		bool bprintmessage=false;
		bool bexporttecplot=false;

		if (iVar==TEMP) {
			integer inum_iter=500; // похоже это номер текущей итерации Simple алгоритма.
		    LR1sK_temp(t, sl, slb, val, col_ind, row_ptr, maxelm, maxbound,  dV, dX0, maxit,inum_iter, bprintmessage, bexporttecplot);
		}
		else {
			LR1sK(f, sl, slb, val, col_ind, row_ptr, maxelm, maxbound, iVar, dV, dX0, maxit,bprintmessage, bexporttecplot);
		}
			            


            if (val!=nullptr) {
		        delete[] val;
			}
			if (col_ind!=nullptr) {
		        delete[] col_ind;
			}
			if (row_ptr!=nullptr) {
		        delete[] row_ptr;
			}

            res_sum=0.0;
	        for (i=0; i<maxelm; ++i) {
		        // внутренность матрицы.
		        doublereal buf=0.0;
		        buf=(sl[i].ap*dX0[sl[i].iP]-dV[sl[i].iP]);
		        if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) buf-=sl[i].ab*dX0[sl[i].iB];
	            if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) buf-=sl[i].ae*dX0[sl[i].iE];
		        if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) buf-=sl[i].an*dX0[sl[i].iN];
		        if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) buf-=sl[i].as*dX0[sl[i].iS];
		        if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) buf-=sl[i].at*dX0[sl[i].iT];
		        if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) buf-=sl[i].aw*dX0[sl[i].iW];
	            buf*=buf;
		        res_sum+=buf;
	        }
	        for (i=0; i<maxbound; ++i) {
	    	   // граничные узлы.
		       doublereal buf=0.0;
		       buf=slb[i].aw*dX0[slb[i].iW]-dV[slb[i].iW];
		       if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) buf-=slb[i].ai*dX0[slb[i].iI];
		       buf*=buf;
		       res_sum+=buf;
	        }
#if doubleprecision == 1 
             	res_sum = sqrt(res_sum);
			
#else 
				res_sum = sqrtf(res_sum);
#endif
	       //printf("residual finish=%1.4e\n",res_sum);
	      // system("pause");
	       finish_residual_lr1sk=res_sum; // значение невязки решённой задачи.
		  // system("pause");

		}



	}

		 calculation_main_end_time=clock();
         calculation_vorst_seach_time+=calculation_main_end_time-calculation_main_start_time;
}

// Если использовать gmres с предобуславливателем AINV Bridson то нужно брать m_restart=32. 
// Так рекомендуют в статье для решения уравнений теплопередачи с числом неизвестных до 1.5М.
// Если использовать GMRES в качестве предобуславливателя то рекомендуют брать m_restart=2.
// Даный метод существенно замедляется при увеличении числа итераций. Например он не сходится на Диффузионно дрейфовой модели для 
// материала GaN в то время как BiCGStab сходится.
// Метод использует ilu2 предобуславливатель.
// Если использовать gmres с предобуславливателем AINV Bridson то нужно брать m_restart=32. 
// Так рекомендуют в статье для решения уравнений теплопередачи с числом неизвестных до 1.5М.
// Если использовать GMRES в качестве предобуславливателя то рекомендуют брать m_restart=2.
// Даный метод существенно замедляется при увеличении числа итераций. Например он не сходится на Диффузионно дрейфовой модели для 
// материала GaN в то время как BiCGStab сходится.
// Метод использует ilu2 предобуславливатель.
integer  fgmres(equation3D* &sl, equation3D_bon* &slb,
	int maxelm, int maxbound, doublereal *dV, doublereal* &dX0,
	integer maxit, integer &m_restart, doublereal alpharelax, bool bprintmessage, integer iVar,
	QuickMemVorst& m, int* &ifrontregulationgl, int* &ibackregulationgl,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls) {

	integer i_1 = 0;
	dterminatedTResudual = 1.0e-11; // Этого достаточно.

									// Мы используем m из BiCGStab_internal3 для хранения матрицы предобуславлдивания.
	//bool brc = false;
	//bool bpam_gsp = false; // предобуславливание с помощью нескольких итераций метода Гаусса-Зейделя.

	doublereal *val = nullptr;
	integer* col_ind = nullptr;
	integer* row_ptr = nullptr;
	integer n = maxelm + maxbound;

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (ibackregulationgl != nullptr) {
			// nested desection версия алгоритма.
			integer ierr = equation3DtoCRSnd(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
		else {
			integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
	}
	if (iVar == TEMP) {
		integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
		if (ierr > 0) {
			printf("Temperature equation problem.\n");
		}
	}

	const integer ILU0 = 0;
	const integer ILU_lfil = 1;

	const integer itype_ilu = ILU_lfil;//ILU_lfil;




								 // Исходная матрица.
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (!m.ballocCRScfd) {
			// m.a=new doublereal[7*n+2]; // CRS
			// m.ja=new integer[7*n+2];
			// 26 сентября 2016.
			m.a = new doublereal[row_ptr[n] + 2 * maxbound + 2];
			m.ja = new integer[row_ptr[n] + 2 * maxbound + 2];
			m.ia = new integer[n + 2];
			// Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((m.a == nullptr) || (m.ja == nullptr) || (m.ia == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}
	}
	if (iVar == TEMP) {
		if (!m.ballocCRSt) {
			// m.ta=new doublereal[7*n+2]; // CRS
			// m.tja=new integer[7*n+2];
			// 26 сентября 2016.
			m.ta = new doublereal[row_ptr[n] + 2 * maxbound + 2];
			m.tja = new integer[row_ptr[n] + 2 * maxbound + 2];
			m.tia = new integer[n + 2];
			// Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((m.ta == nullptr) || (m.tja == nullptr) || (m.tia == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}
	}
	// Эти структуры данных используются в библиотеке Ю.Саада SPARSKIT2.
	// Основной момент который следует уяснить это:
	// мы используем индексацию массивов в СИ начиная с нуля и до size не включая size;
	// В библиотеке Sparskit2 нумерация элементов массива начинается с единицы и до size включая size.
	// Главное что следует понять не следует пытаться переписать SPARSKIT2 так чтобы она соответствовала нумерации с нуля.
	// оставим нумерацию с единицы иначе мы можем запутаться в логике отлаженного кода SPARSKIT2.
	// Но мы в проекте AliceFlowv0_07 работаем с нумерации начинающейся с нуля. 
	// Код Sparskit2 работает только для предобуславливания. Иными словами схема следующая:
	// На входе матрица в CRS формате с нумерацией с нуля. Преобразуем её в матрицу в которой элементы нумеруются с единицы,
	// для этого очевидно нужно к значениям элементов  col_ind и row_ptr прибавить единицу. Индексация же в матрице (преобразованной пусть начинается с нуля)
	// для этого мы внутри кода SPARSKIT2 декременируем ссылку на массив тогда элемент с номером ноль станет первым тоесть будет иметь номер ноль. Потом в конце использования
	// кода SPARSKIT2 мы увеливаем ссылки на массивы на единицу (откат назад). Получая на вход данную преобразованную матрицу CRS формата, мы делаем с помощью неё матрицу предобуславливателя
	// в MSR формате, как обычно используя код SPARSKIT2 без каких либо изменений из вне (он получен с помощью f2c.exe). Используя матрицу предобуславливания в MSR формате
	// мы пихаем её в функию lusol_ из SPARSKIT2 и получаем необходимый вектор x: (LU)x=y; по вектору y и матрице LU в формате MSR. Код в lusol_ содержит преобразование декремента
	// указателей x, y что позволяет использовать обычные векторы x , y в которых нумерация начинается с нуля. В общем x,y менять не надо, они такие же как обычно в AliceFlowv0_07.
	// В lusol_ указатели сначала декримируются --a; а в конце инкремируются ++a; Так что можно использовать код Sparskit2 без изменений !!!



	if (bprintmessage) {
		printf("Incoplete LU Decomposition begin...\n");
	}


	integer ierr = 0;
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		for (integer i = 0; i<row_ptr[n]; ++i) {
			m.a[i] = val[i];
			m.ja[i] = col_ind[i] + 1;
		}
		for (integer i = 0; i<n + 1; ++i) {
			m.ia[i] = row_ptr[i] + 1;
		}
	}
	if (iVar == TEMP) {
		for (integer i = 0; i<row_ptr[n]; ++i) {
			m.ta[i] = val[i];
			m.tja[i] = col_ind[i] + 1;
		}
		for (integer i = 0; i<n + 1; ++i) {
			m.tia[i] = row_ptr[i] + 1;
		}
	}

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (!m.ballocCRScfd) {
			//m.ri = new doublereal[n]; m.roc = new doublereal[n]; m.s = new doublereal[n]; m.t = new doublereal[n]; m.vec = new doublereal[n];
			//m.vi = new doublereal[n]; m.pi = new doublereal[n]; m.dx = new doublereal[n]; m.dax = new doublereal[n];
			m.y = new doublereal[n];// m.z = new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
									/*
									if ((m.ri == nullptr) || (m.roc == nullptr) || (m.s == nullptr) || (m.t == nullptr) || (m.vi == nullptr) || (m.pi == nullptr) || (m.dx == nullptr) || (m.dax == nullptr) || (m.y == nullptr) || (m.z == nullptr)) {
									// недостаточно памяти на данном оборудовании.
									printf("Problem: not enough memory on your equipment...\n");
									printf("Please any key to exit...\n");
									exit(1);
									}
									*/
		}
	}
	if (iVar == TEMP) {
		if (!m.ballocCRSt) {
			//m.tri = new doublereal[n]; m.troc = new doublereal[n]; m.ts = new doublereal[n]; m.tt = new doublereal[n];
			//m.tvi = new doublereal[n]; m.tpi = new doublereal[n]; m.tdx = new doublereal[n]; m.tdax = new doublereal[n];
			m.ty = new doublereal[n]; //m.tz = new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
									  /*
									  if ((m.tri == nullptr) || (m.troc == nullptr) || (m.ts == nullptr) || (m.tt == nullptr) || (m.tvi == nullptr) || (m.tpi == nullptr) || (m.tdx == nullptr) || (m.tdax == nullptr) || (m.ty == nullptr) || (m.tz == nullptr)) {
									  // недостаточно памяти на данном оборудовании.
									  printf("Problem: not enough memory on your equipment...\n");
									  printf("Please any key to exit...\n");
									  exit(1);
									  }
									  */
		}
	}

if (itype_ilu == ILU0) 
	{

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {

			if (!m.ballocCRScfd) {
				//m.alu=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.jlu=new integer[7*n+2];
				// 26 сентября 2016.
				m.alu = new doublereal[row_ptr[n] + 2 * maxbound + 2];
				m.jlu = new integer[row_ptr[n] + 2 * maxbound + 2];

				m.ju = new integer[n + 2];
				if (ibackregulationgl != nullptr) {
					// m.alu1=new doublereal[7*n+2]; // +2 запас по памяти.
					// m.jlu1=new integer[7*n+2];
					// m.ju1=new integer[n+2];
					m.x1 = new doublereal[n + 2];
				}
				//m.alurc=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.jlurc=new integer[7*n+2];
				// 26 сентября 2016.
				m.alurc = new doublereal[row_ptr[n] + 2 * maxbound + 2];
				m.jlurc = new integer[row_ptr[n] + 2 * maxbound + 2];

				m.jurc = new integer[n + 2];
				m.iw = new integer[n + 2]; // рабочий массив.
				m.ballocCRScfd = true; // память выделена.

				if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.ju == nullptr) || (m.iw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
				if ((m.alu1 == nullptr) || (m.jlu1 == nullptr) || (m.ju1 == nullptr) || (m.x1 == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}
		if (iVar == TEMP) {
			if (!m.ballocCRSt) {
				//m.talu=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.tjlu=new integer[7*n+2];
				// 26 сентября 2016.
				m.talu = new doublereal[row_ptr[n] + 2 * maxbound + 2];
				m.tjlu = new integer[row_ptr[n] + 2 * maxbound + 2];

				m.tju = new integer[n + 2];
				m.tiw = new integer[n + 2]; // рабочий массив.
				m.ballocCRSt = true; // память выделена.

				if ((m.talu == nullptr) || (m.tjlu == nullptr) || (m.tju == nullptr) || (m.tiw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
			ilu0_(n, m.a, m.ja, m.ia, m.alu, m.jlu, m.ju, m.iw, ierr);
			/* if (ibackregulationgl!=nullptr) {
			for (integer i87=0; i87<7*n+2; i87++) {
			m.alu1[i87]= m.alu[i87];
			m.jlu1[i87]=m.jlu[i87];
			}
			for (integer i87=0; i87<n+2; i87++) {
			m.ju1[i87]=m.ju[i87];
			}
			}*/
		}
		if (iVar == TEMP) {
			ilu0_(n, m.ta, m.tja, m.tia, m.talu, m.tjlu, m.tju, m.tiw, ierr);
		}

		if (ierr>0) {
#if doubleintprecision == 1
			printf("%lld string in matrix is zero diagonal element...\n", ierr - 1);
#else
			printf("%d string in matrix is zero diagonal element...\n", ierr - 1);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}
	}


if (itype_ilu == ILU_lfil)
	{

		//bool btemp_quick = m.ballocCRSt;

		integer lfil = 3; // 2 уровня (0, 1, 2)
		lfil = my_amg_manager.lfil;

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
			if (!m.ballocCRScfd) {

				// инициализация.
				m.alu = nullptr;
				m.jlu = nullptr;
				m.ju = nullptr;
				m.alu1 = nullptr;
				m.jlu1 = nullptr;
				m.ju1 = nullptr;
				m.x1 = nullptr;
				m.alurc = nullptr;
				m.jlurc = nullptr;
				m.jurc = nullptr;
				m.levs = nullptr;
				m.w = nullptr;
				m.jw = nullptr;
				m.w_dubl = nullptr;
				m.jw_dubl = nullptr;

				//m.iwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
				// 26 сентября 2016.
				if (lfil <= 2) {
					m.iwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil == 3) {
					m.iwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (/*(lfil >= 4) && */(lfil <= 5)) {
					m.iwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil >= 6) {
					m.iwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}

				m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
				m.jlu = new integer[m.iwk + 2];
				m.ju = new integer[n + 2];
				if (ibackregulationgl != nullptr) {
					//m.alu1=new doublereal[m.iwk+2]; // +2 запас по памяти.
					//m.jlu1=new integer[m.iwk+2];
					//m.ju1=new integer[n+2];
					m.x1 = new doublereal[n + 2];
				}
				m.alurc = new doublereal[m.iwk + 2]; // +2 запас по памяти.
				m.jlurc = new integer[m.iwk + 2];
				m.jurc = new integer[n + 2];
				m.levs = new integer[m.iwk + 2]; // уровень.
				m.w = new doublereal[n + 2]; // +2 запас по памяти.
				m.w_dubl = new doublereal[n + 2]; // +2 запас по памяти.

				if (lfil <= 2) {
					m.jw = new integer[3 * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[3 * n + 2]; // +2 запас по памяти.
				}
				else if (lfil == 3) {
					m.jw = new integer[3 * lfil * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[3 * lfil * n + 2]; // +2 запас по памяти.
				}
				else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
					m.jw = new integer[4 * lfil * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[4 * lfil * n + 2]; // +2 запас по памяти.
				}
				else if (lfil >= 6) {
					m.jw = new integer[5 * lfil * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[5 * lfil * n + 2]; // +2 запас по памяти.
				}

				m.ballocCRScfd = true; // память выделена.

				if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.levs == nullptr) || (m.ju == nullptr) || (m.w == nullptr) || (m.jw == nullptr) || (m.w_dubl == nullptr) || (m.jw_dubl == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}
		if (iVar == TEMP) {
			if (!m.ballocCRSt) {

				// инициализация.
				m.talu = nullptr;
				m.tjlu = nullptr;
				m.tju = nullptr;
				m.tlevs = nullptr;
				m.tw = nullptr;
				m.tjw = nullptr;


				//m.tiwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
				// 26 сентября 2016.
				if (lfil <= 2) {
					m.tiwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil == 3) {
					m.tiwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
					m.tiwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil >= 6) {
					m.tiwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}

				m.talu = new doublereal[m.tiwk + 2]; // +2 запас по памяти.
				m.tjlu = new integer[m.tiwk + 2];
				m.tju = new integer[n + 2];
				m.tlevs = new integer[m.tiwk + 2]; // уровень.
				m.tw = new doublereal[n + 2]; // +2 запас по памяти.
				if (lfil <= 2) {
					m.tjw = new integer[3 * n + 2]; // +2 запас по памяти.
				}
				else if (lfil == 3) {
					m.tjw = new integer[3 * lfil* n + 2]; // +2 запас по памяти.
				}
				else if (/*(lfil >= 4) && */(lfil <= 5)) {
					m.tjw = new integer[4 * lfil* n + 2]; // +2 запас по памяти.
				}
				else if (lfil >= 6) {
					m.tjw = new integer[5 * lfil* n + 2]; // +2 запас по памяти.
				}


				m.ballocCRSt = true; // память выделена.

				if ((m.talu == nullptr) || (m.tjlu == nullptr) || (m.tlevs == nullptr) || (m.tju == nullptr) || (m.tw == nullptr) || (m.tjw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
			// iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
			iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);

			if ((ierr == -2) || (ierr == -3)) {

				integer ipassage = 1; // 4 января 2016.
				do {
					printf("\nPlease WAIT... ... ...\n");

					// задаче не хватило памяти, значит нужно перевыделить !
					if (m.alu != nullptr) delete m.alu;
					if (m.jlu != nullptr) delete m.jlu;
					/* if (ibackregulationgl!=nullptr) {
					if (m.alu1!=nullptr) delete m.alu1;
					if (m.jlu1!=nullptr) delete m.jlu1;
					}*/
					if (m.alurc != nullptr) delete m.alurc;
					if (m.jlurc != nullptr) delete m.jlurc;
					if (m.levs != nullptr) delete m.levs;

					// инициализация !
					m.alu = nullptr;
					m.jlu = nullptr;
					/* if (ibackregulationgl!=nullptr) {
					m.alu1=nullptr;
					m.jlu1=nullptr;
					}*/
					m.levs = nullptr;

					// инициализация !
					m.alurc = nullptr;
					m.jlurc = nullptr;

					//m.iwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
					// 26 сентября 2016.
					if (lfil <= 2) {
						m.iwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil == 3) {
						m.iwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
						m.iwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil >= 6) {
						m.iwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}


					m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
					m.jlu = new integer[m.iwk + 2];
					/* (ibackregulationgl!=nullptr) {
					m.alu1=new doublereal[m.iwk+2]; // +2 запас по памяти.
					m.jlu1=new integer[m.iwk+2];
					}*/
					m.levs = new integer[m.iwk + 2]; // уровень.

					if ((m.alu != nullptr) && (m.jlu != nullptr) && (m.levs != nullptr)) {
						// iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
						iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);
						/*
						if (ibackregulationgl!=nullptr) {
						for (integer i87=0; i87<m.iwk+2; i87++) {
						m.alu1[i87]= m.alu[i87];
						m.jlu1[i87]=m.jlu[i87];
						}
						for (integer i87=0; i87<n+2; i87++) {
						m.ju1[i87]=m.ju[i87];
						}
						}*/

					}
					else {
						// недостаточно памяти на данном оборудовании.
						ipassage = 4;
						printf("Problem: not enough memory on your equipment...\n");
						printf("Please any key to exit...\n");
						exit(1);

					}

					ipassage++;
				} while ((ierr != 0) && (ipassage<4));

				if (ipassage == 4) {
					printf("Error memory alloc !!!\n");
					printf("failed to obtain an expansion for the 4 approaches...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
			else {
				/*
				if (ibackregulationgl!=nullptr) {
				for (integer i87=0; i87<m.iwk+2; i87++) {
				m.alu1[i87]= m.alu[i87];
				m.jlu1[i87]=m.jlu[i87];
				}
				for (integer i87=0; i87<n+2; i87++) {
				m.ju1[i87]=m.ju[i87];
				}
				}*/
			}

		}
		else if (iVar == TEMP) {

			/*
			if (0&&bglobal_unsteady_temperature_determinant) {
			// модификация 20_10_2016.
			// При нестационарном моделировании в твёрдом теле будем строить предобуславливатель лишь единожды на первом шаге.
			//if (!btemp_quick) {
			// iluk_call speed_up
			// 10%=2 9%
			// 15%=3 9%  19.25
			// 20%=4 9%
			// 30%=6 3%
			if (rand()%20<3) {
			iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
			}
			}
			else {
			*/
			iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
			//}

			if ((ierr == -2) || (ierr == -3)) {

				integer ipassage = 1;
				do {
					printf("\nPlease WAIT... ... ...\n");

					// задаче не хватило памяти, значит нужно перевыделить !
					if (m.talu != nullptr) delete m.talu;
					if (m.tjlu != nullptr) delete m.tjlu;
					if (m.tlevs != nullptr) delete m.tlevs;

					// инициализация !
					m.talu = nullptr;
					m.tjlu = nullptr;
					m.tlevs = nullptr;

					//m.tiwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
					// 26 сентября 2016.
					if (lfil <= 2) {
						m.tiwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil == 3) {
						m.tiwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (/*(lfil >= 4) && */(lfil <= 5)) {
						m.tiwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil >= 6) {
						m.tiwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}

					m.talu = new doublereal[m.tiwk + 2]; // +2 запас по памяти.
					m.tjlu = new integer[m.tiwk + 2];
					m.tlevs = new integer[m.tiwk + 2]; // уровень.

					if ((m.talu != nullptr) && (m.tjlu != nullptr) && (m.tlevs != nullptr)) {
						iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
					}
					else {
						// недостаточно памяти на данном оборудовании.
						ipassage = 4;
						printf("Problem: not enough memory on your equipment...\n");
						printf("Please any key to exit...\n");
						exit(1);

					}

					ipassage++;
				} while ((ierr != 0) && (ipassage<4));

				if (ipassage == 4) {
					printf("Error memory alloc !!!\n");
					printf("failed to obtain an expansion for the 4 approaches...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
		}



		if (ierr != 0) {
#if doubleintprecision == 1
			printf("error memory in iluk ierr=%lld\n", ierr);
#else
			printf("error memory in iluk ierr=%d\n", ierr);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}
	}


	if (bprintmessage) {
		printf("Incoplete LU Decomposition finish...\n");
	}



	//bool bnorelax = true; // Для уравнения теплопроводности не используется релаксация.


	doublereal resid;
	integer i, j = 1, k;
	//Vector s(m + 1), cs(m + 1), sn(m + 1), w;
	doublereal* w = new doublereal[n];
	doublereal* s = new doublereal[m_restart + 2];
	doublereal* cs = new doublereal[m_restart + 2];
	doublereal* sn = new doublereal[m_restart + 2];

	doublereal *dx = new doublereal[n];
	doublereal *buffer = new doublereal[n];


	// начальное приближение
	// X0 ==
	// под X0 понимается вектор поля температур к примеру.
	if (dX0 == nullptr) {
		dX0 = new doublereal[n];
		for (i = 0; i<n; ++i) {
			dx[i] = 0.0;
			dX0[i] = 0.0;
		}
	}
	else {
		for (i = 0; i<n; ++i) dx[i] = dX0[i];
	}

	//doublereal normb = norm(M.solve(b));
	doublereal normb = 0.0;
	// здесь реализованы все три нормы
	// вообще говоря они все эквивалентны


	// Kbuffer=dV

	// (LU)buffer=dV; 
	/*
	if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
	// Очень важно начинать с нуля иначе не будет сходимости.
	#pragma omp parallel for shared(m) private(i_1) schedule (guided)
	for (i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

	//  9 августа 2015 при внедрении перенумерации узлов nested desection
	if (bpam_gsp && (iVar == PAM)) {
	if (ibackregulationgl != nullptr) {
	PAMGSPnd(sl, slb, m.y, w, maxelm, maxbound, ifrontregulationgl);
	}
	else {
	PAMGSP(sl, slb, m.y, w, maxelm, maxbound);
	}
	}
	else {

	if (brc) {
	for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
	for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
	m.alurc[i7] = m.alu[i7];
	m.jlurc[i7] = m.jlu[i7];
	}
	for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
	}

	if (ibackregulationgl != nullptr) {
	//lusol_2(n, w, m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=w;
	lusol_3(n, w, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=w;
	}
	else {
	lusol_(n, w, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=w;

	}

	if (brc) {
	for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
	for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
	m.alu[i7] = m.alurc[i7];
	m.jlu[i7] = m.jlurc[i7];
	}
	for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
	}

	}
	for (i_1 = 0; i_1 < n; ++i_1) w[i_1] = m.y[i_1];


	}
	*/
	/*
	if (iVar == TEMP) {
	// Очень важно начинать с нуля иначе не будет сходимости.
	#pragma omp parallel for shared(m) private(i) schedule (guided)
	for (i_1 = 0; i_1<n; ++i_1) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

	lusol_(n, dV, m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=w;
	for (i_1 = 0; i_1 < n; ++i_1) buffer[i_1] = m.ty[i_1];

	}
	*/
	normb = NormaV_for_gmres(dV, n);
	//normb = NormaV(buffer, n);

	//Vector r = M.solve(dV - A * x);
	doublereal *r = new doublereal[n];
	MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // результат занесён в  r
	for (i = 0; i < n; ++i) r[i] = dV[i] - r[i];

	//  calculate residual precontidioning;

	/*
	// Ky=r

	// (LU)y=r;
	if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
	// Очень важно начинать с нуля иначе не будет сходимости.
	#pragma omp parallel for shared(m) private(i_1) schedule (guided)
	for (i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

	//  9 августа 2015 при внедрении перенумерации узлов nested desection
	if (bpam_gsp && (iVar == PAM)) {
	if (ibackregulationgl != nullptr) {
	PAMGSPnd(sl, slb, m.y, r, maxelm, maxbound, ifrontregulationgl);
	}
	else {
	PAMGSP(sl, slb, m.y, r, maxelm, maxbound);
	}
	}
	else {

	if (brc) {
	for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
	for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
	m.alurc[i7] = m.alu[i7];
	m.jlurc[i7] = m.jlu[i7];
	}
	for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
	}

	if (ibackregulationgl != nullptr) {
	//lusol_2(n, v[0], m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=r;
	lusol_3(n, r, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=r;
	}
	else {
	lusol_(n, r, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=r;

	}

	if (brc) {
	for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
	for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
	m.alu[i7] = m.alurc[i7];
	m.jlu[i7] = m.jlurc[i7];
	}
	for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
	}

	}
	for (i_1 = 0; i_1 < n; ++i_1) r[i_1] = m.y[i_1];


	}
	if (iVar == TEMP) {
	// Очень важно начинать с нуля иначе не будет сходимости.
	#pragma omp parallel for shared(m) private(i) schedule (guided)
	for (i_1 = 0; i_1<n; ++i_1) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

	lusol_(n, r, m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=r;
	for (i_1 = 0; i_1 < n; ++i_1) r[i_1] = m.ty[i_1];

	}
	*/
	//doublereal beta = norm(r);
	doublereal beta = 0.0;



	beta = NormaV_for_gmres(r, n);

	if (fabs(normb) < 1.0e-30)
		normb = 1;

	doublereal norm_r = 0.0;


	norm_r = NormaV_for_gmres(r, n);

	if ((resid = norm_r / normb) <= dterminatedTResudual) {
		//tol = resid;
		maxit = 0;
		delete[] w;
		delete[] s;
		delete[] cs;
		delete[] sn;
		delete[] buffer;
		return 0;
	}

	doublereal** H = new doublereal*[m_restart + 2]; // Hessenberg
	for (i_1 = 0; i_1 < m_restart + 2; ++i_1) H[i_1] = new doublereal[m_restart + 2];
	for (i_1 = 0; i_1 < m_restart + 2; ++i_1)
	{
		for (integer j_1 = 0; j_1 < m_restart + 2; ++j_1)
		{
			H[i_1][j_1] = 0.0;
		}
	}

	//Vector *v = new Vector[m_restart + 1];
	doublereal** v = new doublereal*[m_restart + 2];
	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) v[i_1] = new doublereal[n];
	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			v[i_1][j_1] = 0.0;
		}
	}

	doublereal** Z = new doublereal*[m_restart + 2];
	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) Z[i_1] = new doublereal[n];
	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			Z[i_1][j_1] = 0.0;
		}
	}

	j = 1; // номер первой итерации
		   //doublereal delta = 1.0e-3;// DOPOLNENIE

	//integer i_copy;

	while (j <= maxit) {

		//v[0] = r * (1.0 / beta);    // ??? r / beta
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			v[0][j_1] = r[j_1] * (1.0 / beta);
		}

		//s = 0.0;
		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) s[i_1] = 0.0;
		//s[0] = beta;
		s[0] = 1.0;

		/*
		for (i_1 = 0; i_1 < m_restart + 2; ++i_1)
		{ // DOPOLNENIE
		for (integer j_1 = 0; j_1 < m_restart + 2; ++j_1)
		{
		H[i_1][j_1] = 0.0;
		}
		}
		*/

		// Ортогонализация Арнольди.
		for (i = 0; i < m_restart && j <= maxit; i++, j++) {

			//i_copy = i;


			// KZ[i]=v[i]

			// (LU)Z[i]=v[i];
			/*
			if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
			// Очень важно начинать с нуля иначе не будет сходимости.
			#pragma omp parallel for shared(m) private(i_1) schedule (guided)
			for (i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

			//  9 августа 2015 при внедрении перенумерации узлов nested desection
			if (bpam_gsp && (iVar == PAM)) {
			if (ibackregulationgl != nullptr) {
			PAMGSPnd(sl, slb, m.y, w, maxelm, maxbound, ifrontregulationgl);
			}
			else {
			PAMGSP(sl, slb, m.y, w, maxelm, maxbound);
			}
			}
			else {

			if (brc) {
			for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
			for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
			m.alurc[i7] = m.alu[i7];
			m.jlurc[i7] = m.jlu[i7];
			}
			for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
			}

			if (ibackregulationgl != nullptr) {
			//lusol_2(n, w, m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=w;
			lusol_3(n, w, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=w;
			}
			else {
			lusol_(n, w, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=w;

			}

			if (brc) {
			for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
			for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
			m.alu[i7] = m.alurc[i7];
			m.jlu[i7] = m.jlurc[i7];
			}
			for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
			}

			}
			//for (i_1 = 0; i_1 < n; ++i_1) w[i_1] = m.y[i_1];
			for (i_1 = 0; i_1 < n; ++i_1)  v[i + 1][i_1] = m.y[i_1];


			}
			*/

			if (iVar == TEMP) {
				// Очень важно начинать с нуля иначе не будет сходимости.
				//#pragma omp parallel for shared(m) private(i) schedule (guided)
				for (integer  i_1 = 0; i_1<n; ++i_1) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

				lusol_(n, v[i], m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=v[i];
				for (integer  i_1 = 0; i_1 < n; ++i_1) Z[i][i_1] = m.ty[i_1];
				//for (integer  i_1 = 0; i_1 < n; ++i_1) v[i + 1][i_1] = m.ty[i_1];

			}

			// Совсем без предобуславливателя.
			//for (i_1 = 0; i_1 < n; ++i_1) Z[i][i_1] = v[i][i_1];

			// закомментировано без предобуславливания.
			//w = M.solve(A * v[i]);
			MatrixCRSByVector(val, col_ind, row_ptr, Z[i], w, n); // результат занесён в  w


																  //doublereal av = sqrt(Scal(w,w,n)); // DOPOLNENIE
																  //doublereal av = sqrt(Scal(v[i + 1], v[i + 1], n)); // DOPOLNENIE

			for (k = 0; k <= i; ++k) {
				H[k][i] = Scal(w, v[k], n);
				//H[k][i] = Scal(v[i + 1], v[k], n);
				for (integer j_1 = 0; j_1 < n; ++j_1)
				{
					//v[i + 1][j_1] -= H[k][i] * v[k][j_1];
					w[j_1] -= H[k][i] * v[k][j_1];
				}
			}
			//H[i + 1][i] = norm(w);
			H[i + 1][i] = NormaV_for_gmres(w, n);
			//H[i + 1][i] = NormaV(v[i + 1], n);

			/*
			// DOPOLNENIE
			if ((av + delta * H[i+1][i]) == av)
			{
			for (k = 0; k <= i; ++k)//j
			{
			//doublereal htmp = Scal(w,v[k],n);
			doublereal htmp = Scal(v[i + 1], v[k], n);
			//htmp = r8vec_dot(n, v + k*n, v + (j - 1)*n);
			H[k][i] = H[k][i] + htmp;
			//h[(j - 1) + (k - 1)*(mr + 1)] = h[(j - 1) + (k - 1)*(mr + 1)] + htmp;
			for (integer j_1 = 0; j_1 < n; ++j_1)
			{
			v[i][j_1] = v[i][j_1] - htmp * v[k][j_1];
			}
			}
			H[i+1][i] = sqrt(Scal( v[i] , v[i],n));
			}

			if (H[i + 1][i] != 0.0) {
			for (integer j_1 = 0; j_1 < n; ++j_1)
			{
			//v[i + 1][j_1] = w[j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			v[i + 1][j_1] = v[i + 1][j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			}
			}
			*/
			for (integer j_1 = 0; j_1 < n; ++j_1)
			{
				v[i + 1][j_1] = w[j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
															  //v[i + 1][j_1] = v[i + 1][j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			}
			// Окончание ортогонализации Арнольди.
			// В v - хранится ортонормированный базис подпространства Крылова размерности m_restart.
			// H - Верхнетреугольная матрица Хессенберга - матрица коэффициентов ортогонализации.
			/*
			if (0 < i) {
			for (k = 0; k <= i + 1; ++k)
			{
			m.ty[k] = H[k][i];
			}
			for (k = 0; k <= i-1; ++k)
			{
			mult_givens(cs[k], sn[k], k, m.ty);
			}
			for (k = 0; k <= i + 1; ++k)
			{
			H[k][i] = m.ty[k];
			}

			}

			doublereal mu = sqrt(pow(H[i][i], 2)
			+ pow(H[i+1][i], 2));
			cs[i] = H[i][i] / mu;
			sn[i] = -H[i+1][i] / mu;
			H[i][i] = cs[i] * H[i][i] - sn[i] * H[i+1][i];
			H[i+1][i] = 0;
			mult_givens(cs[i], sn[i], i, s);
			*/
			//rho = fabs(s[i+1]);

			/*
			itr_used = itr_used + 1;

			if ( verbose )
			{
			cout << "  K =   " << k << "  Residual = " << rho << "\n";
			}

			if ( rho <= rho_tol && rho <= tol_abs )
			{
			break;
			}
			*/

			// 26.11.2017
			// Это проверенный и испытанный кусок кода.
			for (k = 0; k < i; ++k)
				ApplyPlaneRotation(H[k][i], H[k + 1][i], cs[k], sn[k]);

			GeneratePlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);

			/*
			// Это гораздо лучше.
			//step 1
			for (k = 0; k < i; ++k) {
			//ApplyPlaneRotation
			H[k][i] = cs[k] * H[k][i] + sn[k] * H[k + 1][i];
			H[k+1][i] = -sn[k] * H[k][i] + cs[k] * H[k + 1][i];
			}
			// step 2
			cs[i] = fabs(H[i][i]) / sqrt((H[i][i])*(H[i][i])+(H[i+1][i])*(H[i+1][i]));
			sn[i] = cs[i] * H[i + 1][i] / H[i][i];
			// step 3
			//ApplyPlaneRotation
			s[i] = cs[i] * s[i] + sn[i] * s[i+1];
			s[i+1] = -sn[i] * s[i] + cs[i] * s[i+1];
			// step 4
			H[i][i] = cs[i] * H[i][i] + sn[i] * H[i + 1][i];
			H[i + 1][i] = 0.0;
			*/

			// Вручную устраняем случай полного совпадения невязок на двух соседних итерациях,
			// т.к. иначе это приводит к развалу решения.
			//if (fabs(s[i] - s[i + 1]) < 1.0e-37) s[i + 1] = 1.05*s[i];

			//printf("%d %e \n", j, fabs(s[i + 1]) / normb);
			printf("%lld %e \n", j, beta*fabs(s[i + 1]));
			//system("pause");

			//resid = fabs(s[i + 1]) / normb;
			resid = beta*fabs(s[i + 1]);

			if ((resid) < dterminatedTResudual) {
				printf("dosrochnji vjhod\n");
				//system("pause");
				//Update(dx, i, n, H, s, v);
				Update_flexible(dx, i, n, H, s, v, beta, Z);// Внимание изменил на i-1.
															//tol = resid;
															//maxit = j;
				for (integer i_1 = 0; i_1<n; ++i_1) {
					dX0[i_1] = dx[i_1];

				}

				for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] v[i_1];
				delete[] v;
				for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] Z[i_1];
				delete[] Z;
				delete[] dx;
				delete[] buffer;
				delete[] r;
				delete[] w;
				delete[] s;
				delete[] cs;
				delete[] sn;
				for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) delete[] H[i_1];
				delete[] H;
				delete[] val;
				delete[] col_ind;
				delete[] row_ptr;
				return 0;

			}
		}



		//Update(dx, m_restart - 1, n, H, s, v);//i-1 //ERROR
		//--->Update(dx, i-1, n, H, s, v);//i-1 //ERROR
		Update_flexible(dx, i - 1, n, H, s, v, beta, Z);
		//r = M.solve(b - A * x);
		MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // Результат занесён в r
		for (integer  i_1 = 0; i_1 < n; ++i_1) r[i_1] = dV[i_1] - r[i_1];

		//  calculate residual precontidioning;

		// Ky=r
		/*
		// (LU)y=r;
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(m) private(i_1) schedule (guided)
		for (i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

		//  9 августа 2015 при внедрении перенумерации узлов nested desection
		if (bpam_gsp && (iVar == PAM)) {
		if (ibackregulationgl != nullptr) {
		PAMGSPnd(sl, slb, m.y, r, maxelm, maxbound, ifrontregulationgl);
		}
		else {
		PAMGSP(sl, slb, m.y, r, maxelm, maxbound);
		}
		}
		else {

		if (brc) {
		for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
		for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
		m.alurc[i7] = m.alu[i7];
		m.jlurc[i7] = m.jlu[i7];
		}
		for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
		}

		if (ibackregulationgl != nullptr) {
		//lusol_2(n, v[0], m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=r;
		lusol_3(n, r, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=r;
		}
		else {
		lusol_(n, r, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=r;

		}

		if (brc) {
		for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
		for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
		m.alu[i7] = m.alurc[i7];
		m.jlu[i7] = m.jlurc[i7];
		}
		for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
		}

		}
		for (i_1 = 0; i_1 < n; ++i_1) r[i_1] = m.y[i_1];


		}
		if (iVar == TEMP) {
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(m) private(i) schedule (guided)
		for (i_1 = 0; i_1<n; ++i_1) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

		lusol_(n, r, m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=r;
		for (i_1 = 0; i_1 < n; ++i_1) r[i_1] = m.ty[i_1];

		}
		*/

		/*
		i = i_copy - 1;
		m.ty[i] = s[i] / H[i][i];

		for (k = i; 0 <= k; k--)
		{
		m.ty[k] = s[k];
		for (j = k + 1; j <= i + 1; ++j)
		{
		m.ty[k] = m.ty[k] -H[k][j] * m.ty[j];
		}
		m.ty[k] = m.ty[k] / H[k][k];
		}

		for (k = 0; k < n; ++k)
		{
		for (j = 1; j <= k + 1; ++j)
		{
		dx[k] = dx[k] + v[j][k] * m.ty[j];
		}
		}
		*/
		/*
		if (rho <= rho_tol && rho <= tol_abs)
		{
		break;
		}
		*/

		//beta = norm(r);
		beta = NormaV_for_gmres(r, n);

		//resid = beta / normb;
		resid = beta;

		if ((resid) < dterminatedTResudual) {
			//tol = resid;
			//maxit = j;

			printf("end\n");
			//system("pause");

			for (integer i_1 = 0; i_1<n; ++i_1) {
				dX0[i_1] = dx[i_1];

			}
			for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] v[i_1];
			delete[] v;
			for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] Z[i_1];
			delete[] Z;
			delete[] dx;
			delete[] buffer;
			delete[] r;
			delete[] w;
			delete[] s;
			delete[] cs;
			delete[] sn;
			for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) delete[] H[i_1];
			delete[] H;
			delete[] val;
			delete[] col_ind;
			delete[] row_ptr;
			return 0;

		}
	}

	//tol = resid;
	for (i_1 = 0; i_1<n; ++i_1) {
		dX0[i_1] = dx[i_1];

	}
	for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] v[i_1];
	delete[] v;
	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] Z[i_1];
	delete[] Z;
	delete[] dx;
	delete[] buffer;
	delete[] r;
	delete[] w;
	delete[] s;
	delete[] cs;
	delete[] sn;
	for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) delete[] H[i_1];
	delete[] H;
	delete[] val;
	delete[] col_ind;
	delete[] row_ptr;
	return 1;

} //fgmres


// Централизованное выделение оперативной памяти совсем не увеличивает быстродействия 28.11.2017.
// Вариант соответствующий литературе.
// Если использовать gmres с предобуславливателем AINV Bridson то нужно брать m_restart=32. 
// Так рекомендуют в статье для решения уравнений теплопередачи с числом неизвестных до 1.5М.
// Если использовать GMRES в качестве предобуславливателя то рекомендуют брать m_restart=2.
// Даный метод существенно замедляется при увеличении числа итераций. Например он не сходится на Диффузионно дрейфовой модели для 
// материала GaN в то время как BiCGStab сходится.
// Метод использует ilu2 предобуславливатель.
integer  fgmres1(equation3D* &sl, equation3D_bon* &slb,
	int maxelm, int maxbound, doublereal *dV, doublereal* &dX0,
	integer maxit, int &m_restart, doublereal alpharelax, bool bprintmessage, integer iVar,
	QuickMemVorst& m, int* &ifrontregulationgl, int* &ibackregulationgl,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls) {

	//integer i_1 = 0;
	dterminatedTResudual = 1.0e-7; // recomended ( запас 1e-11. Этого достаточно.)

	// Мы используем m из BiCGStab_internal3 для хранения матрицы предобуславлдивания.
	bool brc = false;
	bool bpam_gsp = false; // предобуславливание с помощью нескольких итераций метода Гаусса-Зейделя.

	doublereal *val = nullptr;
	integer* col_ind = nullptr;
	integer* row_ptr = nullptr;
	integer n = maxelm + maxbound;

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar==NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		if (ibackregulationgl != nullptr) {
			// nested desection версия алгоритма.
			integer ierr = equation3DtoCRSnd(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				case NUSHA: printf("NU  equation problem.\n");  break;
				case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY  equation problem.\n");  break;
				case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA  equation problem.\n");  break;
				case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS  equation problem.\n"); break;
				case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS  equation problem.\n"); break;
				case RE_THETA_LANGTRY_MENTER: printf("TURBULENT_RE_THETA_LANGTRY_MENTER  equation  problem\n"); break;
				case GAMMA_LANGTRY_MENTER: printf("TURBULENT_GAMMA_LANGTRY_MENTER  equation  problem\n");	break;
				}
			}
		}
		else {
			integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
	}
	if (iVar == TEMP) {
		integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
		if (ierr > 0) {
			printf("Temperature equation problem.\n");
		}
	}

	const integer ILU0 = 0;
	const integer ILU_lfil = 1;

	const integer itype_ilu = ILU_lfil;//ILU_lfil;




								 // Исходная матрица.
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		if (!m.ballocCRScfd) {
			// m.a=new doublereal[7*n+2]; // CRS
			// m.ja=new integer[7*n+2];
			// 26 сентября 2016.
			m.a = new doublereal[row_ptr[n] + 2 * maxbound + 2];
			m.ja = new integer[row_ptr[n] + 2 * maxbound + 2];
			m.ia = new integer[n + 2];
			// Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((m.a == nullptr) || (m.ja == nullptr) || (m.ia == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}
	}
	if (iVar == TEMP) {
		if (!m.ballocCRSt) {
			// m.ta=new doublereal[7*n+2]; // CRS
			// m.tja=new integer[7*n+2];
			// 26 сентября 2016.
			m.ta = new doublereal[row_ptr[n] + 2 * maxbound + 2];
			m.tja = new integer[row_ptr[n] + 2 * maxbound + 2];
			m.tia = new integer[n + 2];
			// Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((m.ta == nullptr) || (m.tja == nullptr) || (m.tia == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}
	}
	// Эти структуры данных используются в библиотеке Ю.Саада SPARSKIT2.
	// Основной момент который следует уяснить это:
	// мы используем индексацию массивов в СИ начиная с нуля и до size не включая size;
	// В библиотеке Sparskit2 нумерация элементов массива начинается с единицы и до size включая size.
	// Главное что следует понять не следует пытаться переписать SPARSKIT2 так чтобы она соответствовала нумерации с нуля.
	// оставим нумерацию с единицы иначе мы можем запутаться в логике отлаженного кода SPARSKIT2.
	// Но мы в проекте AliceFlowv0_07 работаем с нумерации начинающейся с нуля. 
	// Код Sparskit2 работает только для предобуславливания. Иными словами схема следующая:
	// На входе матрица в CRS формате с нумерацией с нуля. Преобразуем её в матрицу в которой элементы нумеруются с единицы,
	// для этого очевидно нужно к значениям элементов  col_ind и row_ptr прибавить единицу. Индексация же в матрице (преобразованной пусть начинается с нуля)
	// для этого мы внутри кода SPARSKIT2 декременируем ссылку на массив тогда элемент с номером ноль станет первым тоесть будет иметь номер ноль. Потом в конце использования
	// кода SPARSKIT2 мы увеливаем ссылки на массивы на единицу (откат назад). Получая на вход данную преобразованную матрицу CRS формата, мы делаем с помощью неё матрицу предобуславливателя
	// в MSR формате, как обычно используя код SPARSKIT2 без каких либо изменений из вне (он получен с помощью f2c.exe). Используя матрицу предобуславливания в MSR формате
	// мы пихаем её в функию lusol_ из SPARSKIT2 и получаем необходимый вектор x: (LU)x=y; по вектору y и матрице LU в формате MSR. Код в lusol_ содержит преобразование декремента
	// указателей x, y что позволяет использовать обычные векторы x , y в которых нумерация начинается с нуля. В общем x,y менять не надо, они такие же как обычно в AliceFlowv0_07.
	// В lusol_ указатели сначала декримируются --a; а в конце инкремируются ++a; Так что можно использовать код Sparskit2 без изменений !!!



	if (bprintmessage) {
		printf("Incoplete LU Decomposition begin...\n");
	}


	integer ierr = 0;
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {

		for (integer i = 0; i<row_ptr[n]; ++i) {
			m.a[i] = val[i];
			m.ja[i] = col_ind[i] + 1;
		}
		for (integer i = 0; i<n + 1; ++i) {
			m.ia[i] = row_ptr[i] + 1;
		}
	}
	if (iVar == TEMP) {
		for (integer i = 0; i<row_ptr[n]; ++i) {
			m.ta[i] = val[i];
			m.tja[i] = col_ind[i] + 1;
		}
		for (integer i = 0; i<n + 1; ++i) {
			m.tia[i] = row_ptr[i] + 1;
		}
	}

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		if (!m.ballocCRScfd) {
			//m.ri = new doublereal[n]; m.roc = new doublereal[n]; m.s = new doublereal[n]; m.t = new doublereal[n]; m.vec = new doublereal[n];
			//m.vi = new doublereal[n]; m.pi = new doublereal[n]; m.dx = new doublereal[n]; m.dax = new doublereal[n];
			m.y = new doublereal[n];// m.z = new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
									/*
									if ((m.ri == nullptr) || (m.roc == nullptr) || (m.s == nullptr) || (m.t == nullptr) || (m.vi == nullptr) || (m.pi == nullptr) || (m.dx == nullptr) || (m.dax == nullptr) || (m.y == nullptr) || (m.z == nullptr)) {
									// недостаточно памяти на данном оборудовании.
									printf("Problem: not enough memory on your equipment...\n");
									printf("Please any key to exit...\n");
									exit(1);
									}
									*/
		}
	}
	if (iVar == TEMP) {
		if (!m.ballocCRSt) {
			//m.tri = new doublereal[n]; m.troc = new doublereal[n]; m.ts = new doublereal[n]; m.tt = new doublereal[n];
			//m.tvi = new doublereal[n]; m.tpi = new doublereal[n]; m.tdx = new doublereal[n]; m.tdax = new doublereal[n];
			m.ty = new doublereal[n]; //m.tz = new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
									  /*
									  if ((m.tri == nullptr) || (m.troc == nullptr) || (m.ts == nullptr) || (m.tt == nullptr) || (m.tvi == nullptr) || (m.tpi == nullptr) || (m.tdx == nullptr) || (m.tdax == nullptr) || (m.ty == nullptr) || (m.tz == nullptr)) {
									  // недостаточно памяти на данном оборудовании.
									  printf("Problem: not enough memory on your equipment...\n");
									  printf("Please any key to exit...\n");
									  exit(1);
									  }
									  */
		}
	}

if (itype_ilu == ILU0)
	{

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {

			if (!m.ballocCRScfd) {
				//m.alu=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.jlu=new integer[7*n+2];
				// 26 сентября 2016.
				m.alu = new doublereal[row_ptr[n] + 2 * maxbound + 2];
				m.jlu = new integer[row_ptr[n] + 2 * maxbound + 2];

				m.ju = new integer[n + 2];
				if (ibackregulationgl != nullptr) {
					// m.alu1=new doublereal[7*n+2]; // +2 запас по памяти.
					// m.jlu1=new integer[7*n+2];
					// m.ju1=new integer[n+2];
					m.x1 = new doublereal[n + 2];
				}
				//m.alurc=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.jlurc=new integer[7*n+2];
				// 26 сентября 2016.
				m.alurc = new doublereal[row_ptr[n] + 2 * maxbound + 2];
				m.jlurc = new integer[row_ptr[n] + 2 * maxbound + 2];

				m.jurc = new integer[n + 2];
				m.iw = new integer[n + 2]; // рабочий массив.
				m.ballocCRScfd = true; // память выделена.

				if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.ju == nullptr) || (m.iw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
				if ((m.alu1 == nullptr) || (m.jlu1 == nullptr) || (m.ju1 == nullptr) || (m.x1 == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}
		if (iVar == TEMP) {
			if (!m.ballocCRSt) {
				//m.talu=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.tjlu=new integer[7*n+2];
				// 26 сентября 2016.
				m.talu = new doublereal[row_ptr[n] + 2 * maxbound + 2];
				m.tjlu = new integer[row_ptr[n] + 2 * maxbound + 2];

				m.tju = new integer[n + 2];
				m.tiw = new integer[n + 2]; // рабочий массив.
				m.ballocCRSt = true; // память выделена.

				if ((m.talu == nullptr) || (m.tjlu == nullptr) || (m.tju == nullptr) || (m.tiw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			ilu0_(n, m.a, m.ja, m.ia, m.alu, m.jlu, m.ju, m.iw, ierr);
			/* if (ibackregulationgl!=nullptr) {
			for (integer i87=0; i87<7*n+2; i87++) {
			m.alu1[i87]= m.alu[i87];
			m.jlu1[i87]=m.jlu[i87];
			}
			for (integer i87=0; i87<n+2; i87++) {
			m.ju1[i87]=m.ju[i87];
			}
			}*/
		}
		if (iVar == TEMP) {
			ilu0_(n, m.ta, m.tja, m.tia, m.talu, m.tjlu, m.tju, m.tiw, ierr);
		}

		if (ierr>0) {
#if doubleintprecision == 1
			printf("%lld string in matrix is zero diagonal element...\n", ierr - 1);
#else
			printf("%d string in matrix is zero diagonal element...\n", ierr - 1);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}
	}


if (itype_ilu == ILU_lfil)
	{

		//bool btemp_quick = m.ballocCRSt;

		integer lfil = 3; // 2 уровня (0, 1, 2)

		lfil = my_amg_manager.lfil;

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			if (!m.ballocCRScfd) {

				// инициализация.
				m.alu = nullptr;
				m.jlu = nullptr;
				m.ju = nullptr;
				m.alu1 = nullptr;
				m.jlu1 = nullptr;
				m.ju1 = nullptr;
				m.x1 = nullptr;
				m.alurc = nullptr;
				m.jlurc = nullptr;
				m.jurc = nullptr;
				m.levs = nullptr;
				m.w = nullptr;
				m.jw = nullptr;
				m.w_dubl = nullptr;
				m.jw_dubl = nullptr;

				//m.iwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
				// 26 сентября 2016.
				if (lfil <= 2) {
					m.iwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil == 3) {
					m.iwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (/*(lfil >= 4) && */(lfil <= 5)) {
					m.iwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil >= 6) {
					m.iwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}

				m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
				m.jlu = new integer[m.iwk + 2];
				m.ju = new integer[n + 2];
				if (ibackregulationgl != nullptr) {
					//m.alu1=new doublereal[m.iwk+2]; // +2 запас по памяти.
					//m.jlu1=new integer[m.iwk+2];
					//m.ju1=new integer[n+2];
					m.x1 = new doublereal[n + 2];
				}
				m.alurc = new doublereal[m.iwk + 2]; // +2 запас по памяти.
				m.jlurc = new integer[m.iwk + 2];
				m.jurc = new integer[n + 2];
				m.levs = new integer[m.iwk + 2]; // уровень.
				m.w = new doublereal[n + 2]; // +2 запас по памяти.
				m.w_dubl = new doublereal[n + 2]; // +2 запас по памяти.

				if (lfil <= 2) {
					m.jw = new integer[3 * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[3 * n + 2]; // +2 запас по памяти.
				}
				else if (lfil == 3) {
					m.jw = new integer[3 * lfil * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[3 * lfil * n + 2]; // +2 запас по памяти.
				}
				else if (/*(lfil >= 4) && */(lfil <= 5)) {
					m.jw = new integer[4 * lfil * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[4 * lfil * n + 2]; // +2 запас по памяти.
				}
				else if (lfil >= 6) {
					m.jw = new integer[5 * lfil * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[5 * lfil * n + 2]; // +2 запас по памяти.
				}

				m.ballocCRScfd = true; // память выделена.

				if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.levs == nullptr) || (m.ju == nullptr) || (m.w == nullptr) || (m.jw == nullptr) || (m.w_dubl == nullptr) || (m.jw_dubl == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}
		if (iVar == TEMP) {
			if (!m.ballocCRSt) {

				// инициализация.
				m.talu = nullptr;
				m.tjlu = nullptr;
				m.tju = nullptr;
				m.tlevs = nullptr;
				m.tw = nullptr;
				m.tjw = nullptr;


				//m.tiwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
				// 26 сентября 2016.
				if (lfil <= 2) {
					m.tiwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil == 3) {
					m.tiwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (/*(lfil >= 4) && */(lfil <= 5)) {
					m.tiwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil >= 6) {
					m.tiwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}

				m.talu = new doublereal[m.tiwk + 2]; // +2 запас по памяти.
				m.tjlu = new integer[m.tiwk + 2];
				m.tju = new integer[n + 2];
				m.tlevs = new integer[m.tiwk + 2]; // уровень.
				m.tw = new doublereal[n + 2]; // +2 запас по памяти.
				if (lfil <= 2) {
					m.tjw = new integer[3 * n + 2]; // +2 запас по памяти.
				}
				else if (lfil == 3) {
					m.tjw = new integer[3 * lfil* n + 2]; // +2 запас по памяти.
				}
				else if (/*(lfil >= 4) && */(lfil <= 5)) {
					m.tjw = new integer[4 * lfil* n + 2]; // +2 запас по памяти.
				}
				else if (lfil >= 6) {
					m.tjw = new integer[5 * lfil* n + 2]; // +2 запас по памяти.
				}


				m.ballocCRSt = true; // память выделена.

				if ((m.talu == nullptr) || (m.tjlu == nullptr) || (m.tlevs == nullptr) || (m.tju == nullptr) || (m.tw == nullptr) || (m.tjw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
			(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			// iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
			iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);

			if ((ierr == -2) || (ierr == -3)) {

				integer ipassage = 1; // 4 января 2016.
				do {
					printf("\nPlease WAIT... ... ...\n");

					// задаче не хватило памяти, значит нужно перевыделить !
					if (m.alu != nullptr) delete m.alu;
					if (m.jlu != nullptr) delete m.jlu;
					/* if (ibackregulationgl!=nullptr) {
					if (m.alu1!=nullptr) delete m.alu1;
					if (m.jlu1!=nullptr) delete m.jlu1;
					}*/
					if (m.alurc != nullptr) delete m.alurc;
					if (m.jlurc != nullptr) delete m.jlurc;
					if (m.levs != nullptr) delete m.levs;

					// инициализация !
					m.alu = nullptr;
					m.jlu = nullptr;
					/* if (ibackregulationgl!=nullptr) {
					m.alu1=nullptr;
					m.jlu1=nullptr;
					}*/
					m.levs = nullptr;

					// инициализация !
					m.alurc = nullptr;
					m.jlurc = nullptr;

					//m.iwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
					// 26 сентября 2016.
					if (lfil <= 2) {
						m.iwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil == 3) {
						m.iwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (/*(lfil >= 4) && */(lfil <= 5)) {
						m.iwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil >= 6) {
						m.iwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}


					m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
					m.jlu = new integer[m.iwk + 2];
					/* (ibackregulationgl!=nullptr) {
					m.alu1=new doublereal[m.iwk+2]; // +2 запас по памяти.
					m.jlu1=new integer[m.iwk+2];
					}*/
					m.levs = new integer[m.iwk + 2]; // уровень.

					if ((m.alu != nullptr) && (m.jlu != nullptr) && (m.levs != nullptr)) {
						// iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
						iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);
						/*
						if (ibackregulationgl!=nullptr) {
						for (integer i87=0; i87<m.iwk+2; i87++) {
						m.alu1[i87]= m.alu[i87];
						m.jlu1[i87]=m.jlu[i87];
						}
						for (integer i87=0; i87<n+2; i87++) {
						m.ju1[i87]=m.ju[i87];
						}
						}*/

					}
					else {
						// недостаточно памяти на данном оборудовании.
						ipassage = 4;
						printf("Problem: not enough memory on your equipment...\n");
						printf("Please any key to exit...\n");
						exit(1);

					}

					ipassage++;
				} while ((ierr != 0) && (ipassage<4));

				if (ipassage == 4) {
					printf("Error memory alloc !!!\n");
					printf("failed to obtain an expansion for the 4 approaches...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
			else {
				/*
				if (ibackregulationgl!=nullptr) {
				for (integer i87=0; i87<m.iwk+2; i87++) {
				m.alu1[i87]= m.alu[i87];
				m.jlu1[i87]=m.jlu[i87];
				}
				for (integer i87=0; i87<n+2; i87++) {
				m.ju1[i87]=m.ju[i87];
				}
				}*/
			}

		}
		else if (iVar == TEMP) {

			/*
			if (0&&bglobal_unsteady_temperature_determinant) {
			// модификация 20_10_2016.
			// При нестационарном моделировании в твёрдом теле будем строить предобуславливатель лишь единожды на первом шаге.
			//if (!btemp_quick) {
			// iluk_call speed_up
			// 10%=2 9%
			// 15%=3 9%  19.25
			// 20%=4 9%
			// 30%=6 3%
			if (rand()%20<3) {
			iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
			}
			}
			else {
			*/
			iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
			//}

			if ((ierr == -2) || (ierr == -3)) {

				integer ipassage = 1;
				do {
					printf("\nPlease WAIT... ... ...\n");

					// задаче не хватило памяти, значит нужно перевыделить !
					if (m.talu != nullptr) delete m.talu;
					if (m.tjlu != nullptr) delete m.tjlu;
					if (m.tlevs != nullptr) delete m.tlevs;

					// инициализация !
					m.talu = nullptr;
					m.tjlu = nullptr;
					m.tlevs = nullptr;

					//m.tiwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
					// 26 сентября 2016.
					if (lfil <= 2) {
						m.tiwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil == 3) {
						m.tiwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (/*(lfil >= 4) && */(lfil <= 5)) {
						m.tiwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil >= 6) {
						m.tiwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}

					m.talu = new doublereal[m.tiwk + 2]; // +2 запас по памяти.
					m.tjlu = new integer[m.tiwk + 2];
					m.tlevs = new integer[m.tiwk + 2]; // уровень.

					if ((m.talu != nullptr) && (m.tjlu != nullptr) && (m.tlevs != nullptr)) {
						iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
					}
					else {
						// недостаточно памяти на данном оборудовании.
						ipassage = 4;
						printf("Problem: not enough memory on your equipment...\n");
						printf("Please any key to exit...\n");
						exit(1);

					}

					ipassage++;
				} while ((ierr != 0) && (ipassage<4));

				if (ipassage == 4) {
					printf("Error memory alloc !!!\n");
					printf("failed to obtain an expansion for the 4 approaches...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
		}



		if (ierr != 0) {
#if doubleintprecision == 1
			printf("error memory in iluk ierr=%lld\n", ierr);
#else
			printf("error memory in iluk ierr=%d\n", ierr);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}
	}


	if (bprintmessage) {
		printf("Incoplete LU Decomposition finish...\n");
	}



	//bool bnorelax = true; // Для уравнения теплопроводности не используется релаксация.


	doublereal resid;
	integer  j = 1, k;
	//Vector s(m + 1), cs(m + 1), sn(m + 1), w;
	doublereal* w = new doublereal[n];
	doublereal* s = new doublereal[m_restart + 2];
	doublereal* cs = new doublereal[m_restart + 2];
	doublereal* sn = new doublereal[m_restart + 2];

	doublereal *dx = new doublereal[n];
	doublereal *buffer = new doublereal[n];


	// начальное приближение
	// X0 ==
	// под X0 понимается вектор поля температур к примеру.
	if (dX0 == nullptr) {
		dX0 = new doublereal[n];

#pragma omp parallel for
		for (int i = 0; i<n; ++i) {
			dx[i] = 0.0;
			dX0[i] = 0.0;
		}
	}
	else {
#pragma omp parallel for
		for (int i = 0; i<n; ++i) dx[i] = dX0[i];
	}

	//doublereal normb = norm(M.solve(b));
	doublereal normb = 0.0;
	// здесь реализованы все три нормы
	// вообще говоря они все эквивалентны



	normb = NormaV_for_gmres(dV, n);
	//normb = NormaV(buffer, n);

	//Vector r = dV - A * x;
	doublereal *r = new doublereal[n];
	MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // результат занесён в  r

#pragma omp parallel for
	for (int i = 0; i < n; ++i) r[i] = dV[i] - r[i];

	//  calculate residual precontidioning;


	//doublereal beta = norm(r);
	doublereal beta = 0.0;



	beta = NormaV_for_gmres(r, n);

	if (fabs(normb) < 1.0e-30)
		normb = 1;

	doublereal norm_r = 0.0;


	norm_r = NormaV_for_gmres(r, n);

	if ((resid = norm_r / normb) <= dterminatedTResudual) {
		//tol = resid;
		maxit = 0;
		delete[] w;
		w = nullptr;
		delete[] s;
		s = nullptr;
		delete[] cs;
		cs = nullptr;
		delete[] sn;
		sn = nullptr;
		delete[] buffer;
		buffer = nullptr;

		return 0;
	}

	if ((iVar == TEMP) && (bSIMPLErun_now_for_temperature)) {
		// 07.02.2022
		// Понижаем допуск чтобы побыстрее считало.
		// Только в случае температуры и гидродинамического рамсчета.
		//0.01
		//dterminatedTResudual = 0.1 * fabs(norm_r / normb);


	}
	dterminatedTResudual = 1.0e-7;

	doublereal** H = new doublereal*[m_restart + 2]; // Hessenberg
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
	doublereal** v = new doublereal*[m_restart + 2];
	for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) v[i_1] = new doublereal[n];


#pragma omp parallel for
	for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			v[i_1][j_1] = 0.0;
		}
	}

	doublereal** Z = new doublereal*[m_restart + 2];
	for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) Z[i_1] = new doublereal[n];

#pragma omp parallel for
	for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			Z[i_1][j_1] = 0.0;
		}
	}

	j = 1; // номер первой итерации
		   //doublereal delta = 1.0e-3;// DOPOLNENIE

	//integer i_copy;

	while (j <= maxit) {

		//v[0] = r * (1.0 / beta);    // ??? r / beta
#pragma omp parallel for
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			v[0][j_1] = r[j_1] * (1.0 / beta);
		}

		//s = 0.0;
#pragma omp parallel for
		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) s[i_1] = 0.0;
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


		integer i = 0;

		// Ортогонализация Арнольди.
		for (i = 0; i < m_restart && j <= maxit; i++, j++) {

			//i_copy = i;


			// KZ[i]=v[i]

			// (LU)Z[i]=v[i];

			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
				(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
				(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
				(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
				// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for shared(m)  schedule (guided)
				for (integer  i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

															//  9 августа 2015 при внедрении перенумерации узлов nested desection
				if (bpam_gsp && (iVar == PAM)) {
					if (ibackregulationgl != nullptr) {
						PAMGSPnd(sl, slb, m.y, v[i], maxelm, maxbound, ifrontregulationgl);
					}
					else {
						PAMGSP(sl, slb, m.y, v[i], maxelm, maxbound);
					}
				}
				else {

					if (brc) {
						for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
						for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
							m.alurc[i7] = m.alu[i7];
							m.jlurc[i7] = m.jlu[i7];
						}
						for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
					}

					if (ibackregulationgl != nullptr) {
						//lusol_2(n, v[i], m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=v[i];
						lusol_3(n, v[i], m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=v[i];
					}
					else {
						lusol_(n, v[i], m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=v[i];

					}

					if (brc) {
						for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
						for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
							m.alu[i7] = m.alurc[i7];
							m.jlu[i7] = m.jlurc[i7];
						}
						for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
					}

				}
				//for (integer i_1 = 0; i_1 < n; ++i_1) w[i_1] = m.y[i_1];
				for (integer i_1 = 0; i_1 < n; ++i_1) Z[i][i_1] = m.y[i_1];
				//for (integer i_1 = 0; i_1 < n; ++i_1)  v[i + 1][i_1] = m.y[i_1];


			}


			if (iVar == TEMP) {
				// Очень важно начинать с нуля иначе не будет сходимости.
				//#pragma omp parallel for shared(m) private(i) schedule (guided)
				for (integer i_1 = 0; i_1<n; ++i_1) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

				lusol_(n, v[i], m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=v[i];
				for (integer i_1 = 0; i_1 < n; ++i_1) Z[i][i_1] = m.ty[i_1];
				//for (integer i_1 = 0; i_1 < n; ++i_1) v[i + 1][i_1] = m.ty[i_1];

			}

			// Совсем без предобуславливателя.
			//for (integer i_1 = 0; i_1 < n; ++i_1) Z[i][i_1] = v[i][i_1];

			// закомментировано без предобуславливания.
			//w = A * Z[i];
			MatrixCRSByVector(val, col_ind, row_ptr, Z[i], w, n); // результат занесён в  w

			for (k = 0; k <= i; ++k) {
				H[k][i] = Scal(w, v[k], n);

#pragma omp parallel for
				for (integer j_1 = 0; j_1 < n; ++j_1)
				{
					w[j_1] -= H[k][i] * v[k][j_1];
				}
			}
			H[i + 1][i] = NormaV_for_gmres(w, n);


#pragma omp parallel for
			for (integer j_1 = 0; j_1 < n; ++j_1)
			{
				v[i + 1][j_1] = w[j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			}
			// Окончание ортогонализации Арнольди.
			// В v - хранится ортонормированный базис подпространства Крылова размерности m_restart.
			// H - Верхнетреугольная матрица Хессенберга - матрица коэффициентов ортогонализации.


			// 26.11.2017
			// Это проверенный и испытанный кусок кода.
			for (k = 0; k < i; ++k)
				ApplyPlaneRotation(H[k][i], H[k + 1][i], cs[k], sn[k]);

			GeneratePlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);



			// Вручную устраняем случай полного совпадения невязок на двух соседних итерациях,
			// т.к. иначе это приводит к развалу решения.
			//if (fabs(s[i] - s[i + 1]) < 1.0e-37) s[i + 1] = 1.05*s[i];

			//----->printf("%lld %e \n", j, fabs(s[i + 1]) / normb);
			//printf("%d %e \n", j, beta*fabs(s[i + 1]));
			//system("pause");

			//resid = fabs(s[i + 1]) / normb;
			beta = NormaV_for_gmres(r, n);

			resid = beta / normb;

			if (j%30==0) printf("%lld %e \n", j, resid);
			//resid = beta*fabs(s[i + 1]);

			if ((resid) < dterminatedTResudual) {
				//------->printf("dosrochnji vjhod\n");
				//system("pause");				
				Update(dx, i, n, H, s, Z);
				//tol = resid;
				//maxit = j;

#pragma omp parallel for
				for (integer i_1 = 0; i_1<n; ++i_1) {
					dX0[i_1] = dx[i_1];
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
				delete[] val;
				delete[] col_ind;
				delete[] row_ptr;
				return 0;

			}
		}



		// i-1 -> m_restart-1
		Update(dx, i - 1, n, H, s, Z);//i-1 //ERROR

									  //r = M.solve(b - A * x);
		MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // Результат занесён в r

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < n; ++i_1) r[i_1] = dV[i_1] - r[i_1];

		//beta = norm(r);
		beta = NormaV_for_gmres(r, n);

		resid = beta / normb;
		//resid = beta;

		if ((resid) < dterminatedTResudual) {
			//tol = resid;
			//maxit = j;

			//--------->printf("end\n");
			//system("pause");

#pragma omp parallel for
			for (integer i_1 = 0; i_1<n; ++i_1) {
				dX0[i_1] = dx[i_1];

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
			delete[] val;
			delete[] col_ind;
			delete[] row_ptr;
			return 0;

		}
	}

	//tol = resid;

#pragma omp parallel for
	for (integer i_1 = 0; i_1<n; ++i_1) {
		dX0[i_1] = dx[i_1];

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
	delete[] val;
	delete[] col_ind;
	delete[] row_ptr;
	return 1;


} //fgmres1

  // Централизованное выделение оперативной памяти совсем не увеличивает быстродействия 28.11.2017.
  // Вариант соответствующий литературе.
  // Если использовать gmres с предобуславливателем AINV Bridson то нужно брать m_restart=32. 
  // Так рекомендуют в статье для решения уравнений теплопередачи с числом неизвестных до 1.5М.
  // Если использовать GMRES в качестве предобуславливателя то рекомендуют брать m_restart=2.
  // Даный метод существенно замедляется при увеличении числа итераций. Например он не сходится на Диффузионно дрейфовой модели для 
  // материала GaN в то время как BiCGStab сходится.
  // Метод использует ilu2 предобуславливатель.
  // fgmres + SSOR Хреновый метод, его явно недостаточно. Никогда не использовать.
integer  fgmres2(equation3D* &sl, equation3D_bon* &slb,
	int maxelm, int maxbound, doublereal *dV, doublereal* &dX0,
	integer maxit, integer &m_restart, doublereal alpharelax, bool bprintmessage, integer iVar,
	QuickMemVorst& m, int* &ifrontregulationgl, int* &ibackregulationgl,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls) {

	integer i_1 = 0;
	dterminatedTResudual = 1.0e-7; // recomended ( запас 1e-11. Этого достаточно.)

								   // Мы используем m из BiCGStab_internal3 для хранения матрицы предобуславлдивания.
	//bool brc = false;
	//bool bpam_gsp = false; // предобуславливание с помощью нескольких итераций метода Гаусса-Зейделя.

	doublereal *val = nullptr;
	integer* col_ind = nullptr;
	integer* row_ptr = nullptr;
	integer n = maxelm + maxbound;

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM) || (iVar == NUSHA) ||
		(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
		(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
		if (ibackregulationgl != nullptr) {
			// nested desection версия алгоритма.
			integer ierr = equation3DtoCRSnd(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				case NUSHA: printf("NU  equation problem.\n");  break;
				case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY  equation problem.\n");  break;
				case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA  equation problem.\n");  break;
				case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS  equation problem.\n"); break;
				case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS  equation problem.\n"); break;
				case RE_THETA_LANGTRY_MENTER: printf("TURBULENT_RE_THETA_LANGTRY_MENTER  equation problem\n"); break;
				case GAMMA_LANGTRY_MENTER: printf("TURBULENT_GAMMA_LANGTRY_MENTER  equation  problem\n");	break;
				}
			}
		}
		else {
			integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
	}
	else if (iVar == TEMP) {
		integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
		if (ierr > 0) {
			printf("Temperature equation problem.\n");
		}
	}
	else {
		printf("ERROR !!!: unknown iVar in function fgmres2 in my_linalg.cpp\n");
		system("PAUSE");
		exit(1);
	}



	//bool bnorelax = true; // Для уравнения теплопроводности не используется релаксация.


	doublereal resid;
	integer i, j = 1, k;
	//Vector s(m + 1), cs(m + 1), sn(m + 1), w;
	doublereal* w = new doublereal[n];
	doublereal* s = new doublereal[m_restart + 2];
	doublereal* cs = new doublereal[m_restart + 2];
	doublereal* sn = new doublereal[m_restart + 2];

	doublereal *dx = new doublereal[n];
	doublereal *buffer = new doublereal[n];


	// начальное приближение
	// X0 ==
	// под X0 понимается вектор поля температур к примеру.
	if (dX0 == nullptr) {
		dX0 = new doublereal[n];
		for (i = 0; i<n; ++i) {
			dx[i] = 0.0;
			dX0[i] = 0.0;
		}
	}
	else {
		for (i = 0; i<n; ++i) dx[i] = dX0[i];
	}

	//doublereal normb = norm(M.solve(b));
	doublereal normb = 0.0;
	// здесь реализованы все три нормы
	// вообще говоря они все эквивалентны



	normb = NormaV_for_gmres(dV, n);
	//normb = NormaV(buffer, n);

	//Vector r = dV - A * x;
	doublereal *r = new doublereal[n];
	MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // результат занесён в  r
	for (i = 0; i < n; ++i) r[i] = dV[i] - r[i];

	//  calculate residual precontidioning;


	//doublereal beta = norm(r);
	doublereal beta = 0.0;



	beta = NormaV_for_gmres(r, n);

	if (fabs(normb) < 1.0e-30) {
		normb = 1;
	}

	doublereal norm_r = 0.0;


	norm_r = NormaV_for_gmres(r, n);

	if ((resid = norm_r / normb) <= dterminatedTResudual) {
		//tol = resid;
		//maxit = 0;
		delete[] w;
		delete[] s;
		delete[] cs;
		delete[] sn;
		delete[] buffer;
		return 0;
	}



	doublereal** H = new doublereal*[m_restart + 2]; // Hessenberg
	for (i_1 = 0; i_1 < m_restart + 2; ++i_1) H[i_1] = new doublereal[m_restart + 2];


	for (i_1 = 0; i_1 < m_restart + 2; ++i_1)
	{
		for (integer j_1 = 0; j_1 < m_restart + 2; ++j_1)
		{
			H[i_1][j_1] = 0.0;
		}
	}

	//Vector *v = new Vector[m_restart + 1];
	doublereal** v = new doublereal*[m_restart + 2];
	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) v[i_1] = new doublereal[n];


	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			v[i_1][j_1] = 0.0;
		}
	}

	doublereal** Z = new doublereal*[m_restart + 2];
	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) Z[i_1] = new doublereal[n];

	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			Z[i_1][j_1] = 0.0;
		}
	}

	j = 1; // номер первой итерации
		   //doublereal delta = 1.0e-3;// DOPOLNENIE

	//integer i_copy;

	while (j <= maxit) {

		//v[0] = r * (1.0 / beta);    // ??? r / beta
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			v[0][j_1] = r[j_1] * (1.0 / beta);
		}

		//s = 0.0;
		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) s[i_1] = 0.0;
		s[0] = beta;
		//s[0] = 1.0;


		for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1)
		{ // DOPOLNENIE
			for (integer j_1 = 0; j_1 < m_restart + 2; ++j_1)
			{
				H[i_1][j_1] = 0.0;
			}
		}


		// Ортогонализация Арнольди.
		for (i = 0; i < m_restart && j <= maxit; i++, j++) {

			//i_copy = i;


			// KZ[i]=v[i]

			// (LU)Z[i]=v[i];
			/*
			if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)|| (iVar==NUSHA)||
		 (iVar== TURBULENT_KINETIK_ENERGY)||(iVar== TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA)||
			(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
			(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
			// Очень важно начинать с нуля иначе не будет сходимости.
			#pragma omp parallel for shared(m) private(i_1) schedule (guided)
			for (i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

			//  9 августа 2015 при внедрении перенумерации узлов nested desection
			if (bpam_gsp && (iVar == PAM)) {
			if (ibackregulationgl != nullptr) {
			PAMGSPnd(sl, slb, m.y, v[i], maxelm, maxbound, ifrontregulationgl);
			}
			else {
			PAMGSP(sl, slb, m.y, v[i], maxelm, maxbound);
			}
			}
			else {

			if (brc) {
			for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
			for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
			m.alurc[i7] = m.alu[i7];
			m.jlurc[i7] = m.jlu[i7];
			}
			for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
			}

			if (ibackregulationgl != nullptr) {
			//lusol_2(n, v[i], m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=v[i];
			lusol_3(n, v[i], m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=v[i];
			}
			else {
			lusol_(n, v[i], m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=v[i];

			}

			if (brc) {
			for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
			for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
			m.alu[i7] = m.alurc[i7];
			m.jlu[i7] = m.jlurc[i7];
			}
			for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
			}

			}
			//for (i_1 = 0; i_1 < n; ++i_1) w[i_1] = m.y[i_1];
			for (i_1 = 0; i_1 < n; ++i_1) Z[i][i_1] = m.y[i_1];
			//for (i_1 = 0; i_1 < n; ++i_1)  v[i + 1][i_1] = m.y[i_1];


			}
			*/

			if (iVar == TEMP) {

				if (1) {
					// Матрица col_ind, row_ptr, val собрана верно!!! 14 сентября 2017.

					// init
					for (integer i_1 = 0; i_1 < n; ++i_1) buffer[i_1] = v[i][i_1];

					doublereal omega = 1.855; // initialize.

											  // За подробностями смотри книгу Патрика Роуча стр. 183.
											  //doublereal rn = static_cast<doublereal>(n);
											  //optimal_omega(rn, omega);

					for (integer k_1 = 0; k_1 < 7; k_1++) {
						// направление возрастания.
						for (integer i_1 = 0; i_1 < n; ++i_1) {
							doublereal r_1 = v[i][i_1];
							doublereal ap = 0.0;
							for (integer j_1 = row_ptr[i_1]; j_1 <= row_ptr[i_1 + 1] - 1; ++j_1) {
								if (i_1 != col_ind[j_1]) {
									r_1 += -val[j_1] * buffer[col_ind[j_1]];
								}
								else {
									ap = val[j_1];
								}
							}
							//buffer[i_1] = r_1 / ap;
							buffer[i_1] = (1.0 - omega)*buffer[i_1] + ((omega)*(r_1)) / ap;
						}
						// направление убывания.
						for (integer i_1 = n - 1; i_1 >= 0; i_1--) {
							doublereal r_1 = v[i][i_1];
							doublereal ap = 0.0;
							for (integer j_1 = row_ptr[i_1]; j_1 <= row_ptr[i_1 + 1] - 1; ++j_1) {
								if (i_1 != col_ind[j_1]) {
									r_1 += -val[j_1] * buffer[col_ind[j_1]];
								}
								else {
									ap = val[j_1];
								}
							}
							//buffer[i_1] = r_1 / ap;
							buffer[i_1] = (1.0 - omega)*buffer[i_1] + ((omega)*(r_1)) / ap;
						}
					}

					for (integer i_1 = 0; i_1 < n; ++i_1) Z[i][i_1] = buffer[i_1];
				}

			}

			// Совсем без предобуславливателя.
			//for (i_1 = 0; i_1 < n; ++i_1) Z[i][i_1] = v[i][i_1];

			// закомментировано без предобуславливания.
			//w = A * Z[i];
			MatrixCRSByVector(val, col_ind, row_ptr, Z[i], w, n); // результат занесён в  w

			for (k = 0; k <= i; ++k) {
				H[k][i] = Scal(w, v[k], n);

				for (integer j_1 = 0; j_1 < n; ++j_1)
				{
					w[j_1] -= H[k][i] * v[k][j_1];
				}
			}
			H[i + 1][i] = NormaV_for_gmres(w, n);



			for (integer j_1 = 0; j_1 < n; ++j_1)
			{
				v[i + 1][j_1] = w[j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			}
			// Окончание ортогонализации Арнольди.
			// В v - хранится ортонормированный базис подпространства Крылова размерности m_restart.
			// H - Верхнетреугольная матрица Хессенберга - матрица коэффициентов ортогонализации.


			// 26.11.2017
			// Это проверенный и испытанный кусок кода.
			for (k = 0; k < i; ++k)
				ApplyPlaneRotation(H[k][i], H[k + 1][i], cs[k], sn[k]);

			GeneratePlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);



			// Вручную устраняем случай полного совпадения невязок на двух соседних итерациях,
			// т.к. иначе это приводит к развалу решения.
			//if (fabs(s[i] - s[i + 1]) < 1.0e-37) s[i + 1] = 1.05*s[i];

			printf("%lld %e \n", j, fabs(s[i + 1]) / normb);
			//printf("%d %e \n", j, beta*fabs(s[i + 1]));
			//system("pause");

			resid = fabs(s[i + 1]) / normb;
			//resid = beta*fabs(s[i + 1]);

			if ((resid) < dterminatedTResudual) {
				printf("dosrochnji vjhod\n");
				//system("pause");				
				Update(dx, i, n, H, s, Z);
				//tol = resid;
				//maxit = j;

				for (integer i_1 = 0; i_1<n; ++i_1) {
					dX0[i_1] = dx[i_1];
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
				delete[] val;
				delete[] col_ind;
				delete[] row_ptr;
				return 0;

			}
		}



		// i-1 -> m_restart-1
		Update(dx, i - 1, n, H, s, Z);//i-1 //ERROR

									  //r = M.solve(b - A * x);
		MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // Результат занесён в r
		for (integer i_1 = 0; i_1 < n; ++i_1) r[i_1] = dV[i_1] - r[i_1];

		//beta = norm(r);
		beta = NormaV_for_gmres(r, n);

		resid = beta / normb;
		//resid = beta;

		if ((resid) < dterminatedTResudual) {
			//tol = resid;
			//maxit = j;

			printf("end\n");
			//system("pause");

			for (integer i_1 = 0; i_1<n; ++i_1) {
				dX0[i_1] = dx[i_1];

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
			delete[] val;
			delete[] col_ind;
			delete[] row_ptr;
			return 0;

		}
	}

	//tol = resid;
	for (integer i_1 = 0; i_1<n; ++i_1) {
		dX0[i_1] = dx[i_1];

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
	delete[] val;
	delete[] col_ind;
	delete[] row_ptr;
	return 1;


} //fgmres2


// BiCGStabL, L=2 recomended
/* BiCGSTAB(L) algorithm for the n-by-n problem Ax = b */
// этот метод показывает значительно более лучшую сходимость, чем простой BiCGStabCRS,
// а также он гораздо лучше (и повидимому правильней) чем Bi_CGStab_internal1.
// Bi_CGStab_internal3 использует предобуславливание из библиотеки Ю.Саада.
// дата написания Bi_CGStab_internal3: 31.03.2013. 
// раскурочил 9 августа 2015.
// 26 сентября 2016 Теперь метод пригоден и для АЛИС сетки.
// BiCGStabL, L=2 recomended
/* BiCGSTAB(L) algorithm for the n-by-n problem Ax = b */
void Bi_CGStab_internal5(integer L, equation3D* &sl, equation3D_bon* &slb,
	int maxelm, int maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit, doublereal alpharelax,
	bool bprintmessage, integer iVar, QuickMemVorst& mstruct, 
	int* &ifrontregulationgl, int* &ibackregulationgl,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls)
{

	// Мы используем m из BiCGStab_internal3 для хранения матрицы предобуславлдивания.
	bool brc = false;
	bool bpam_gsp = false; // предобуславливание с помощью нескольких итераций метода Гаусса-Зейделя.

	
	integer n = maxelm + maxbound;


	const integer ILU0 = 0;
	const integer ILU_lfil = 1;

	const integer itype_ilu = ILU_lfil;//ILU_lfil;


	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (ibackregulationgl != nullptr) {
			// nested desection версия алгоритма.
			integer ierr = equation3DtoCRSnd(sl, slb, mstruct.val, mstruct.col_ind, mstruct.row_ptr, maxelm, maxbound, alpharelax, !mstruct.ballocCRScfd, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
		else {
			integer ierr = equation3DtoCRS(sl, slb, mstruct.val, mstruct.col_ind, mstruct.row_ptr, maxelm, maxbound, alpharelax, !mstruct.ballocCRScfd, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
	}
	if (iVar == TEMP) {
		integer ierr = equation3DtoCRS(sl, slb, mstruct.tval, mstruct.tcol_ind, mstruct.trow_ptr, maxelm, maxbound, alpharelax, !mstruct.ballocCRSt, b, lb, s_loc, ls);
		if (ierr > 0) {
			printf("Temperature equation problem.\n");
		}
	}


	// Исходная матрица.
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (!mstruct.ballocCRScfd) {
			// mstruct.a=new doublereal[7*n+2]; // CRS
			// mstruct.ja=new integer[7*n+2];
			// 26 сентября 2016.
			mstruct.a = new doublereal[mstruct.row_ptr[n] + 2 * maxbound + 2];
			mstruct.ja = new integer[mstruct.row_ptr[n] + 2 * maxbound + 2];
			mstruct.ia = new integer[n + 2];
			// Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((mstruct.a == nullptr) || (mstruct.ja == nullptr) || (mstruct.ia == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}
	}
	if (iVar == TEMP) {
		if (!mstruct.ballocCRSt) {
			// mstruct.ta=new doublereal[7*n+2]; // CRS
			// mstruct.tja=new integer[7*n+2];
			// 26 сентября 2016.
			mstruct.ta = new doublereal[mstruct.trow_ptr[n] + 2 * maxbound + 2];
			mstruct.tja = new integer[mstruct.trow_ptr[n] + 2 * maxbound + 2];
			mstruct.tia = new integer[n + 2];
			// Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((mstruct.ta == nullptr) || (mstruct.tja == nullptr) || (mstruct.tia == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}
	}
	// Эти структуры данных используются в библиотеке Ю.Саада SPARSKIT2.
	// Основной момент который следует уяснить это:
	// мы используем индексацию массивов в СИ начиная с нуля и до size не включая size;
	// В библиотеке Sparskit2 нумерация элементов массива начинается с единицы и до size включая size.
	// Главное что следует понять не следует пытаться переписать SPARSKIT2 так чтобы она соответствовала нумерации с нуля.
	// оставим нумерацию с единицы иначе мы можем запутаться в логике отлаженного кода SPARSKIT2.
	// Но мы в проекте AliceFlowv0_07 работаем с нумерации начинающейся с нуля. 
	// Код Sparskit2 работает только для предобуславливания. Иными словами схема следующая:
	// На входе матрица в CRS формате с нумерацией с нуля. Преобразуем её в матрицу в которой элементы нумеруются с единицы,
	// для этого очевидно нужно к значениям элементов  col_ind и row_ptr прибавить единицу. Индексация же в матрице (преобразованной пусть начинается с нуля)
	// для этого мы внутри кода SPARSKIT2 декременируем ссылку на массив тогда элемент с номером ноль станет первым тоесть будет иметь номер ноль. Потом в конце использования
	// кода SPARSKIT2 мы увеливаем ссылки на массивы на единицу (откат назад). Получая на вход данную преобразованную матрицу CRS формата, мы делаем с помощью неё матрицу предобуславливателя
	// в MSR формате, как обычно используя код SPARSKIT2 без каких либо изменений из вне (он получен с помощью f2c.exe). Используя матрицу предобуславливания в MSR формате
	// мы пихаем её в функию lusol_ из SPARSKIT2 и получаем необходимый вектор x: (LU)x=y; по вектору y и матрице LU в формате MSR. Код в lusol_ содержит преобразование декремента
	// указателей x, y что позволяет использовать обычные векторы x , y в которых нумерация начинается с нуля. В общем x,y менять не надо, они такие же как обычно в AliceFlowv0_07.
	// В lusol_ указатели сначала декримируются --a; а в конце инкремируются ++a; Так что можно использовать код Sparskit2 без изменений !!!



	if (bprintmessage) {
		printf("Incoplete LU Decomposition begin...\n");
	}


	integer ierr = 0;
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		for (integer i = 0; i<mstruct.row_ptr[n]; ++i) {
			mstruct.a[i] = mstruct.val[i];
			mstruct.ja[i] = mstruct.col_ind[i] + 1;
		}
		for (integer i = 0; i<n + 1; ++i) {
			mstruct.ia[i] = mstruct.row_ptr[i] + 1;
		}
	}
	if (iVar == TEMP) {
		for (integer i = 0; i<mstruct.trow_ptr[n]; ++i) {
			mstruct.ta[i] = mstruct.tval[i];
			mstruct.tja[i] = mstruct.tcol_ind[i] + 1;
		}
		for (integer i = 0; i<n + 1; ++i) {
			mstruct.tia[i] = mstruct.trow_ptr[i] + 1;
		}
	}

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (!mstruct.ballocCRScfd) {
			//mstruct.ri = new doublereal[n]; mstruct.roc = new doublereal[n]; mstruct.s = new doublereal[n]; mstruct.t = new doublereal[n]; mstruct.vec = new doublereal[n];
			//mstruct.vi = new doublereal[n]; mstruct.pi = new doublereal[n]; mstruct.dx = new doublereal[n]; mstruct.dax = new doublereal[n];
			mstruct.y = new doublereal[n];// mstruct.z = new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
										  /*
										  if ((mstruct.ri == nullptr) || (mstruct.roc == nullptr) || (mstruct.s == nullptr) || (mstruct.t == nullptr) || (mstruct.vi == nullptr) || (mstruct.pi == nullptr) || (mstruct.dx == nullptr) || (mstruct.dax == nullptr) || (mstruct.y == nullptr) || (mstruct.z == nullptr)) {
										  // недостаточно памяти на данном оборудовании.
										  printf("Problem: not enough memory on your equipment...\n");
										  printf("Please any key to exit...\n");
										  exit(1);
										  }
										  */
		}
	}
	if (iVar == TEMP) {
		if (!mstruct.ballocCRSt) {
			//mstruct.tri = new doublereal[n]; mstruct.troc = new doublereal[n]; mstruct.ts = new doublereal[n]; mstruct.tt = new doublereal[n];
			//mstruct.tvi = new doublereal[n]; mstruct.tpi = new doublereal[n]; mstruct.tdx = new doublereal[n]; mstruct.tdax = new doublereal[n];
			mstruct.ty = new doublereal[n]; //mstruct.tz = new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
											/*
											if ((mstruct.tri == nullptr) || (mstruct.troc == nullptr) || (mstruct.ts == nullptr) || (mstruct.tt == nullptr) || (mstruct.tvi == nullptr) || (mstruct.tpi == nullptr) || (mstruct.tdx == nullptr) || (mstruct.tdax == nullptr) || (mstruct.ty == nullptr) || (mstruct.tz == nullptr)) {
											// недостаточно памяти на данном оборудовании.
											printf("Problem: not enough memory on your equipment...\n");
											printf("Please any key to exit...\n");
											exit(1);
											}
											*/
		}
	}

if (itype_ilu == ILU0) 
	{

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {

			if (!mstruct.ballocCRScfd) {
				//mstruct.alu=new doublereal[7*n+2]; // +2 запас по памяти.
				//mstruct.jlu=new integer[7*n+2];
				// 26 сентября 2016.
				mstruct.alu = new doublereal[mstruct.row_ptr[n] + 2 * maxbound + 2];
				mstruct.jlu = new integer[mstruct.row_ptr[n] + 2 * maxbound + 2];

				mstruct.ju = new integer[n + 2];
				if (ibackregulationgl != nullptr) {
					// mstruct.alu1=new doublereal[7*n+2]; // +2 запас по памяти.
					// mstruct.jlu1=new integer[7*n+2];
					// mstruct.ju1=new integer[n+2];
					//mstruct.x1 = new doublereal[n + 2];
				}
				//mstruct.alurc=new doublereal[7*n+2]; // +2 запас по памяти.
				//mstruct.jlurc=new integer[7*n+2];
				// 26 сентября 2016.
				mstruct.alurc = new doublereal[mstruct.row_ptr[n] + 2 * maxbound + 2];
				mstruct.jlurc = new integer[mstruct.row_ptr[n] + 2 * maxbound + 2];

				mstruct.jurc = new integer[n + 2];
				mstruct.iw = new integer[n + 2]; // рабочий массив.
				mstruct.ballocCRScfd = true; // память выделена.

				if ((mstruct.alu == nullptr) || (mstruct.jlu == nullptr) || (mstruct.ju == nullptr) || (mstruct.iw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
				if ((mstruct.alu1 == nullptr) || (mstruct.jlu1 == nullptr) || (mstruct.ju1 == nullptr)) {// || (mstruct.x1 == nullptr)) {
																								// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}
		if (iVar == TEMP) {
			if (!mstruct.ballocCRSt) {
				//mstruct.talu=new doublereal[7*n+2]; // +2 запас по памяти.
				//mstruct.tjlu=new integer[7*n+2];
				// 26 сентября 2016.
				mstruct.talu = new doublereal[mstruct.trow_ptr[n] + 2 * maxbound + 2];
				mstruct.tjlu = new integer[mstruct.trow_ptr[n] + 2 * maxbound + 2];

				mstruct.tju = new integer[n + 2];
				mstruct.tiw = new integer[n + 2]; // рабочий массив.
				mstruct.ballocCRSt = true; // память выделена.

				if ((mstruct.talu == nullptr) || (mstruct.tjlu == nullptr) || (mstruct.tju == nullptr) || (mstruct.tiw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
			ilu0_(n, mstruct.a, mstruct.ja, mstruct.ia, mstruct.alu, mstruct.jlu, mstruct.ju, mstruct.iw, ierr);
			/* if (ibackregulationgl!=nullptr) {
			for (integer i87=0; i87<7*n+2; i87++) {
			mstruct.alu1[i87]= mstruct.alu[i87];
			mstruct.jlu1[i87]=mstruct.jlu[i87];
			}
			for (integer i87=0; i87<n+2; i87++) {
			mstruct.ju1[i87]=mstruct.ju[i87];
			}
			}*/
		}
		if (iVar == TEMP) {
			ilu0_(n, mstruct.ta, mstruct.tja, mstruct.tia, mstruct.talu, mstruct.tjlu, mstruct.tju, mstruct.tiw, ierr);
		}

		if (ierr>0) {
#if doubleintprecision == 1
			printf("%lld string in matrix is zero diagonal element...\n", ierr - 1);
#else
			printf("%d string in matrix is zero diagonal element...\n", ierr - 1);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}
	}

if (itype_ilu == ILU_lfil) 
	{

		//bool btemp_quick = mstruct.ballocCRSt;

		integer lfil = 2; // 2 уровня (0, 1, 2)
		lfil = my_amg_manager.lfil;

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
			if (!mstruct.ballocCRScfd) {

				// инициализация.
				mstruct.alu = nullptr;
				mstruct.jlu = nullptr;
				mstruct.ju = nullptr;
				mstruct.alu1 = nullptr;
				mstruct.jlu1 = nullptr;
				mstruct.ju1 = nullptr;
				//mstruct.x1 = nullptr;
				mstruct.alurc = nullptr;
				mstruct.jlurc = nullptr;
				mstruct.jurc = nullptr;
				mstruct.levs = nullptr;
				mstruct.w = nullptr;
				mstruct.jw = nullptr;
				mstruct.w_dubl = nullptr;
				mstruct.jw_dubl = nullptr;

				//mstruct.iwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
				// 26 сентября 2016.
				if (lfil <= 2) {
					mstruct.iwk = (lfil + 1) * (mstruct.row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil == 3) {
					mstruct.iwk = (2 * lfil + 1) * (mstruct.row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
					mstruct.iwk = (3 * lfil + 1) * (mstruct.row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil >= 6) {
					mstruct.iwk = (4 * lfil + 1) * (mstruct.row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}

				mstruct.alu = new doublereal[mstruct.iwk + 2]; // +2 запас по памяти.
				mstruct.jlu = new integer[mstruct.iwk + 2];
				mstruct.ju = new integer[n + 2];
				if (ibackregulationgl != nullptr) {
					//mstruct.alu1=new doublereal[mstruct.iwk+2]; // +2 запас по памяти.
					//mstruct.jlu1=new integer[mstruct.iwk+2];
					//mstruct.ju1=new integer[n+2];
					//mstruct.x1 = new doublereal[n + 2];
				}
				mstruct.alurc = new doublereal[mstruct.iwk + 2]; // +2 запас по памяти.
				mstruct.jlurc = new integer[mstruct.iwk + 2];
				mstruct.jurc = new integer[n + 2];
				mstruct.levs = new integer[mstruct.iwk + 2]; // уровень.
				mstruct.w = new doublereal[n + 2]; // +2 запас по памяти.
				mstruct.w_dubl = new doublereal[n + 2]; // +2 запас по памяти.

				if (lfil <= 2) {
					mstruct.jw = new integer[3 * n + 2]; // +2 запас по памяти.				
					mstruct.jw_dubl = new integer[3 * n + 2]; // +2 запас по памяти.
				}
				else if (lfil == 3) {
					mstruct.jw = new integer[3 * lfil * n + 2]; // +2 запас по памяти.				
					mstruct.jw_dubl = new integer[3 * lfil * n + 2]; // +2 запас по памяти.
				}
				else if (/*(lfil >= 4) && */(lfil <= 5)) {
					mstruct.jw = new integer[4 * lfil * n + 2]; // +2 запас по памяти.				
					mstruct.jw_dubl = new integer[4 * lfil * n + 2]; // +2 запас по памяти.
				}
				else if (lfil >= 6) {
					mstruct.jw = new integer[5 * lfil * n + 2]; // +2 запас по памяти.				
					mstruct.jw_dubl = new integer[5 * lfil * n + 2]; // +2 запас по памяти.
				}

				mstruct.ballocCRScfd = true; // память выделена.

				if ((mstruct.alu == nullptr) || (mstruct.jlu == nullptr) || (mstruct.levs == nullptr) || (mstruct.ju == nullptr) || (mstruct.w == nullptr) || (mstruct.jw == nullptr) || (mstruct.w_dubl == nullptr) || (mstruct.jw_dubl == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}
		if (iVar == TEMP) {
			if (!mstruct.ballocCRSt) {

				// инициализация.
				mstruct.talu = nullptr;
				mstruct.tjlu = nullptr;
				mstruct.tju = nullptr;
				mstruct.tlevs = nullptr;
				mstruct.tw = nullptr;
				mstruct.tjw = nullptr;


				//mstruct.tiwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
				// 26 сентября 2016.
				if (lfil <= 2) {
					mstruct.tiwk = (lfil + 1) * (mstruct.trow_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil == 3) {
					mstruct.tiwk = (2 * lfil + 1) * (mstruct.trow_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
					mstruct.tiwk = (3 * lfil + 1) * (mstruct.trow_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil >= 6) {
					mstruct.tiwk = (4 * lfil + 1) * (mstruct.trow_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}

				mstruct.talu = new doublereal[mstruct.tiwk + 2]; // +2 запас по памяти.
				mstruct.tjlu = new integer[mstruct.tiwk + 2];
				mstruct.tju = new integer[n + 2];
				mstruct.tlevs = new integer[mstruct.tiwk + 2]; // уровень.
				mstruct.tw = new doublereal[n + 2]; // +2 запас по памяти.
				if (lfil <= 2) {
					mstruct.tjw = new integer[3 * n + 2]; // +2 запас по памяти.
				}
				else if (lfil == 3) {
					mstruct.tjw = new integer[3 * lfil* n + 2]; // +2 запас по памяти.
				}
				else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
					mstruct.tjw = new integer[4 * lfil* n + 2]; // +2 запас по памяти.
				}
				else if (lfil >= 6) {
					mstruct.tjw = new integer[5 * lfil* n + 2]; // +2 запас по памяти.
				}


				mstruct.ballocCRSt = true; // память выделена.

				if ((mstruct.talu == nullptr) || (mstruct.tjlu == nullptr) || (mstruct.tlevs == nullptr) || (mstruct.tju == nullptr) || (mstruct.tw == nullptr) || (mstruct.tjw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
			// iluk_(n, mstruct.a, mstruct.ja, mstruct.ia, lfil, mstruct.alu, mstruct.jlu, mstruct.ju, mstruct.levs, mstruct.iwk, mstruct.w, mstruct.jw, ierr);
			iluk_2(n, mstruct.a, mstruct.ja, mstruct.ia, lfil, mstruct.alu, mstruct.jlu, mstruct.ju, mstruct.levs, mstruct.iwk, mstruct.w, mstruct.jw, mstruct.w_dubl, mstruct.jw_dubl, ierr);

			if ((ierr == -2) || (ierr == -3)) {

				integer ipassage = 1; // 4 января 2016.
				do {
					printf("\nPlease WAIT... ... ...\n");

					// задаче не хватило памяти, значит нужно перевыделить !
					if (mstruct.alu != nullptr) delete mstruct.alu;
					if (mstruct.jlu != nullptr) delete mstruct.jlu;
					/* if (ibackregulationgl!=nullptr) {
					if (mstruct.alu1!=nullptr) delete mstruct.alu1;
					if (mstruct.jlu1!=nullptr) delete mstruct.jlu1;
					}*/
					if (mstruct.alurc != nullptr) delete mstruct.alurc;
					if (mstruct.jlurc != nullptr) delete mstruct.jlurc;
					if (mstruct.levs != nullptr) delete mstruct.levs;

					// инициализация !
					mstruct.alu = nullptr;
					mstruct.jlu = nullptr;
					/* if (ibackregulationgl!=nullptr) {
					mstruct.alu1=nullptr;
					mstruct.jlu1=nullptr;
					}*/
					mstruct.levs = nullptr;

					// инициализация !
					mstruct.alurc = nullptr;
					mstruct.jlurc = nullptr;

					//mstruct.iwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
					// 26 сентября 2016.
					if (lfil <= 2) {
						mstruct.iwk = (lfil + 1) * (mstruct.row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil == 3) {
						mstruct.iwk = (2 * lfil + 1) * (mstruct.row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
						mstruct.iwk = (3 * lfil + 1) * (mstruct.row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil >= 6) {
						mstruct.iwk = (4 * lfil + 1) * (mstruct.row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}


					mstruct.alu = new doublereal[mstruct.iwk + 2]; // +2 запас по памяти.
					mstruct.jlu = new integer[mstruct.iwk + 2];
					/* (ibackregulationgl!=nullptr) {
					mstruct.alu1=new doublereal[mstruct.iwk+2]; // +2 запас по памяти.
					mstruct.jlu1=new integer[mstruct.iwk+2];
					}*/
					mstruct.levs = new integer[mstruct.iwk + 2]; // уровень.

					if ((mstruct.alu != nullptr) && (mstruct.jlu != nullptr) && (mstruct.levs != nullptr)) {
						// iluk_(n, mstruct.a, mstruct.ja, mstruct.ia, lfil, mstruct.alu, mstruct.jlu, mstruct.ju, mstruct.levs, mstruct.iwk, mstruct.w, mstruct.jw, ierr);
						iluk_2(n, mstruct.a, mstruct.ja, mstruct.ia, lfil, mstruct.alu, mstruct.jlu, mstruct.ju, mstruct.levs, mstruct.iwk, mstruct.w, mstruct.jw, mstruct.w_dubl, mstruct.jw_dubl, ierr);
						/*
						if (ibackregulationgl!=nullptr) {
						for (integer i87=0; i87<mstruct.iwk+2; i87++) {
						mstruct.alu1[i87]= mstruct.alu[i87];
						mstruct.jlu1[i87]=mstruct.jlu[i87];
						}
						for (integer i87=0; i87<n+2; i87++) {
						mstruct.ju1[i87]=mstruct.ju[i87];
						}
						}*/

					}
					else {
						// недостаточно памяти на данном оборудовании.
						ipassage = 4;
						printf("Problem: not enough memory on your equipment...\n");
						printf("Please any key to exit...\n");
						exit(1);

					}

					ipassage++;
				} while ((ierr != 0) && (ipassage<4));

				if (ipassage == 4) {
					printf("Error memory alloc !!!\n");
					printf("failed to obtain an expansion for the 4 approaches...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
			else {
				/*
				if (ibackregulationgl!=nullptr) {
				for (integer i87=0; i87<mstruct.iwk+2; i87++) {
				mstruct.alu1[i87]= mstruct.alu[i87];
				mstruct.jlu1[i87]=mstruct.jlu[i87];
				}
				for (integer i87=0; i87<n+2; i87++) {
				mstruct.ju1[i87]=mstruct.ju[i87];
				}
				}*/
			}

		}
		else if (iVar == TEMP) {

			/*
			if (0&&bglobal_unsteady_temperature_determinant) {
			// модификация 20_10_2016.
			// При нестационарном моделировании в твёрдом теле будем строить предобуславливатель лишь единожды на первом шаге.
			//if (!btemp_quick) {
			// iluk_call speed_up
			// 10%=2 9%
			// 15%=3 9%  19.25
			// 20%=4 9%
			// 30%=6 3%
			if (rand()%20<3) {
			iluk_(n, mstruct.ta, mstruct.tja, mstruct.tia, lfil, mstruct.talu, mstruct.tjlu, mstruct.tju, mstruct.tlevs, mstruct.tiwk, mstruct.tw, mstruct.tjw, ierr);
			}
			}
			else {
			*/
			iluk_(n, mstruct.ta, mstruct.tja, mstruct.tia, lfil, mstruct.talu, mstruct.tjlu, mstruct.tju, mstruct.tlevs, mstruct.tiwk, mstruct.tw, mstruct.tjw, ierr);
			//}

			if ((ierr == -2) || (ierr == -3)) {

				integer ipassage = 1;
				do {
					printf("\nPlease WAIT... ... ...\n");

					// задаче не хватило памяти, значит нужно перевыделить !
					if (mstruct.talu != nullptr) delete mstruct.talu;
					if (mstruct.tjlu != nullptr) delete mstruct.tjlu;
					if (mstruct.tlevs != nullptr) delete mstruct.tlevs;

					// инициализация !
					mstruct.talu = nullptr;
					mstruct.tjlu = nullptr;
					mstruct.tlevs = nullptr;

					//mstruct.tiwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
					// 26 сентября 2016.
					if (lfil <= 2) {
						mstruct.tiwk = (lfil + 1) * (mstruct.trow_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil == 3) {
						mstruct.tiwk = (2 * lfil + 1) * (mstruct.trow_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (/*(lfil >= 4) &&*/ (lfil <= 5)) {
						mstruct.tiwk = (3 * lfil + 1) * (mstruct.trow_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil >= 6) {
						mstruct.tiwk = (4 * lfil + 1) * (mstruct.trow_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}

					mstruct.talu = new doublereal[mstruct.tiwk + 2]; // +2 запас по памяти.
					mstruct.tjlu = new integer[mstruct.tiwk + 2];
					mstruct.tlevs = new integer[mstruct.tiwk + 2]; // уровень.

					if ((mstruct.talu != nullptr) && (mstruct.tjlu != nullptr) && (mstruct.tlevs != nullptr)) {
						iluk_(n, mstruct.ta, mstruct.tja, mstruct.tia, lfil, mstruct.talu, mstruct.tjlu, mstruct.tju, mstruct.tlevs, mstruct.tiwk, mstruct.tw, mstruct.tjw, ierr);
					}
					else {
						// недостаточно памяти на данном оборудовании.
						ipassage = 4;
						printf("Problem: not enough memory on your equipment...\n");
						printf("Please any key to exit...\n");
						exit(1);

					}

					ipassage++;
				} while ((ierr != 0) && (ipassage<4));

				if (ipassage == 4) {
					printf("Error memory alloc !!!\n");
					printf("failed to obtain an expansion for the 4 approaches...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
		}



		if (ierr != 0) {
#if doubleintprecision == 1
			printf("error memory in iluk ierr=%lld\n", ierr);
#else
			printf("error memory in iluk ierr=%d\n", ierr);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}
	}


	if (bprintmessage) {
		printf("Incoplete LU Decomposition finish...\n");
	}


	//bool bnorelax = true; // Для уравнения теплопроводности не используется релаксация.



						  //doublereal *work = new doublereal[(2*L+3)*n]; // required workspace

	doublereal  **r = new doublereal*[L + 1];
	doublereal  **u = new doublereal*[L + 1];
	// границы массивов
	for (integer i = 0; i <= L; ++i) {
		r[i] = new doublereal[n];
		u[i] = new doublereal[n];
	}


	doublereal *x = new doublereal[n];
	// начальное приближение
	// X0 ==
	// под X0 понимается вектор поля температур к примеру.
	if (dX0 == nullptr) {
		dX0 = new doublereal[n];
		for (integer i = 0; i<n; ++i) {
			x[i] = 0.0;
			dX0[i] = 0.0;
		}
	}
	else {
		for (integer i = 0; i<n; ++i) x[i] = dX0[i];
	}

	//doublereal normb = norm(M.solve(b));
	doublereal bnrm = 0.0;
	// здесь реализованы все три нормы
	// вообще говоря они все эквивалентны

	bnrm = NormaV_for_gmres(dV, n);
	if (bnrm == 0.0) bnrm = 1.0;

	integer iter = 0;

	doublereal *gamma = new doublereal[L + 1];
	doublereal *gamma_p = new doublereal[L + 1];
	doublereal *gamma_pp = new doublereal[L + 1];

	doublereal *tau = new doublereal[L * L];
	doublereal *sigma = new doublereal[L + 1];

	ierr = 0; // error code to return, if any
	//const doublereal breaktol = 1.0e-30;

	/**** FIXME: check for breakdown conditions(?) during iteration  ****/



	// rtilde = r[0] = b - Ax
	//doublereal *rtilde = work + (2 * L + 2) * n;
	doublereal *rtilde = new doublereal[(2 * L + 2) * n];
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		MatrixCRSByVector(mstruct.val, mstruct.col_ind, mstruct.row_ptr, x, r[0], n); // результат занесён в  r[0]
	}
	if (iVar == TEMP) {
		MatrixCRSByVector(mstruct.tval, mstruct.tcol_ind, mstruct.trow_ptr, x, r[0], n); // результат занесён в  r[0]
	}
	for (integer m = 0; m < n; ++m) rtilde[m] = r[0][m] = dV[m] - r[0][m];

	{ /* Sleipjen normalizes rtilde in his code; it seems to help slightly */
		doublereal s = 1.0 / NormaV_for_gmres(rtilde, n);
		for (integer m = 0; m < n; ++m) rtilde[m] *= s;
	}

	doublereal tol = 1.0e-16;

	for (integer m = 0; m < n; m++) u[0][m] = 0.0;

	doublereal rho = 1.0, alpha = 0, omega = 1;

	doublereal resid;

	while ((resid = NormaV_for_gmres(r[0], n)) > tol * bnrm) {
		resid = NormaV_for_gmres(r[0], n);
		if (resid / bnrm < 1.0e-7) break;
		++iter;

		if ((iVar == TEMP) && (bonly_solid_calculation)) {
			printf("%lld %e\n", iter, resid / bnrm);
			//system("pause");
		}

		rho = -omega * rho;
		for (integer j = 0; j < L; ++j) {
			//if (fabs(rho) < breaktol) { ierr = -1; goto finish; }
			doublereal rho1 = Scal(r[j], rtilde, n);
			//doublereal beta = alpha * rho1 / rho;
			doublereal beta = 0.0;

			if ((fabs(alpha * rho1)<1e-30) && (fabs(rho)<1e-30)) {
				beta = 1.0;
			}
			else if (fabs(alpha * rho1)<1e-30) {
				beta = 0.0;
			}
			else {
				beta = alpha * rho1 / rho;
			}

			rho = rho1;
			for (integer i = 0; i <= j; ++i)
				for (integer m = 0; m < n; ++m) u[i][m] = r[i][m] - beta * u[i][m];

			// Ky=u[j]

			// (LU)y=u[j]; 
			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
				// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for
				for (integer i = 0; i<n; ++i) mstruct.y[i] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

																  //  9 августа 2015 при внедрении перенумерации узлов nested desection
				if (bpam_gsp && (iVar == PAM)) {
					if (ibackregulationgl != nullptr) {
						PAMGSPnd(sl, slb, mstruct.y, u[j], maxelm, maxbound, ifrontregulationgl);
					}
					else {
						PAMGSP(sl, slb, mstruct.y, u[j], maxelm, maxbound);
					}
				}
				else {

					if (brc) {
						//for (integer i7 = 0; i7<n; ++i7) mstruct.vec[i7] = mstruct.pi[i7];
						for (integer i7 = 0; i7<mstruct.iwk + 2; ++i7) {
							mstruct.alurc[i7] = mstruct.alu[i7];
							mstruct.jlurc[i7] = mstruct.jlu[i7];
						}
						for (integer i7 = 0; i7<n + 2; ++i7) mstruct.jurc[i7] = mstruct.ju[i7];
					}

					if (ibackregulationgl != nullptr) {
						//lusol_2(n, u[j], mstruct.y, mstruct.alu, mstruct.jlu, mstruct.ju, mstruct.x1, maxelm); // M*y=u[j];
						lusol_3(n, u[j], mstruct.y, mstruct.alu, mstruct.jlu, mstruct.ju, maxelm); // M*y=u[j];
					}
					else {
						lusol_(n, u[j], mstruct.y, mstruct.alu, mstruct.jlu, mstruct.ju, maxelm); // M*y=u[j];

					}

					if (brc) {
						//for (integer i7 = 0; i7<n; ++i7) mstruct.pi[i7] = mstruct.vec[i7];
						for (integer i7 = 0; i7<mstruct.iwk + 2; ++i7) {
							mstruct.alu[i7] = mstruct.alurc[i7];
							mstruct.jlu[i7] = mstruct.jlurc[i7];
						}
						for (integer i7 = 0; i7<n + 2; ++i7) mstruct.ju[i7] = mstruct.jurc[i7];
					}

				}

				MatrixCRSByVector(mstruct.val, mstruct.col_ind, mstruct.row_ptr, mstruct.y, u[j + 1], n); // u[j + 1]==A*y;
			}
			if (iVar == TEMP) {
				// Очень важно начинать с нуля иначе не будет сходимости.
				//#pragma omp parallel for shared(m) schedule (guided)
				for (integer i = 0; i<n; ++i) mstruct.ty[i] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

				lusol_(n, u[j], mstruct.ty, mstruct.talu, mstruct.tjlu, mstruct.tju, maxelm); // M*ty=u[j];
				MatrixCRSByVector(mstruct.tval, mstruct.tcol_ind, mstruct.trow_ptr, mstruct.ty, u[j + 1], n); // u[j + 1]==A*ty;
			}

			//MatrixCRSByVector(val, col_ind, row_ptr, u[j], u[j + 1], n); // результат занесён в  u[j + 1]



			if ((fabs(rho)<1e-30) && (fabs(Scal(u[j + 1], rtilde, n))<1e-30)) {
				alpha = 1.0;
			}
			else if (fabs(rho)<1e-30) {
				alpha = 0.0;
			}
			else {
				alpha = rho / Scal(u[j + 1], rtilde, n);
			}
			//alpha = rho / Scal( u[j + 1], rtilde, n);


			for (integer i = 0; i <= j; ++i) {
				for (integer m = 0; m < n; ++m)  r[i][m] -= alpha * u[i + 1][m];
			}


			// Ky=r[j]

			// (LU)y=r[j]; 
			if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
				// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for
				for (integer i = 0; i<n; ++i) mstruct.y[i] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

																  //  9 августа 2015 при внедрении перенумерации узлов nested desection
				if (bpam_gsp && (iVar == PAM)) {
					if (ibackregulationgl != nullptr) {
						PAMGSPnd(sl, slb, mstruct.y, r[j], maxelm, maxbound, ifrontregulationgl);
					}
					else {
						PAMGSP(sl, slb, mstruct.y, r[j], maxelm, maxbound);
					}
				}
				else {

					if (brc) {
						//for (integer i7 = 0; i7<n; ++i7) mstruct.vec[i7] = mstruct.pi[i7];
						for (integer i7 = 0; i7<mstruct.iwk + 2; ++i7) {
							mstruct.alurc[i7] = mstruct.alu[i7];
							mstruct.jlurc[i7] = mstruct.jlu[i7];
						}
						for (integer i7 = 0; i7<n + 2; ++i7) mstruct.jurc[i7] = mstruct.ju[i7];
					}

					if (ibackregulationgl != nullptr) {
						//lusol_2(n, r[j], mstruct.y, mstruct.alu, mstruct.jlu, mstruct.ju, mstruct.x1, maxelm); // M*y=r[j];
						lusol_3(n, r[j], mstruct.y, mstruct.alu, mstruct.jlu, mstruct.ju, maxelm); // M*y=r[j];
					}
					else {
						lusol_(n, r[j], mstruct.y, mstruct.alu, mstruct.jlu, mstruct.ju, maxelm); // M*y=r[j];

					}

					if (brc) {
						//for (integer i7 = 0; i7<n; ++i7) mstruct.pi[i7] = mstruct.vec[i7];
						for (integer i7 = 0; i7<mstruct.iwk + 2; ++i7) {
							mstruct.alu[i7] = mstruct.alurc[i7];
							mstruct.jlu[i7] = mstruct.jlurc[i7];
						}
						for (integer i7 = 0; i7<n + 2; ++i7) mstruct.ju[i7] = mstruct.jurc[i7];
					}

				}

				MatrixCRSByVector(mstruct.val, mstruct.col_ind, mstruct.row_ptr, mstruct.y, r[j + 1], n); // r[j + 1]==A*y;
			}
			if (iVar == TEMP) {
				// Очень важно начинать с нуля иначе не будет сходимости.
#pragma omp parallel for 
				for (integer i = 0; i<n; ++i) mstruct.ty[i] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

				lusol_(n, r[j], mstruct.ty, mstruct.talu, mstruct.tjlu, mstruct.tju, maxelm); // M*ty=r[j];
				MatrixCRSByVector(mstruct.tval, mstruct.tcol_ind, mstruct.trow_ptr, mstruct.ty, r[j + 1], n); // r[j + 1]==A*ty;
			}


			//MatrixCRSByVector(val, col_ind, row_ptr, r[j], r[j + 1], n); // результат занесён в  r[j + 1]

			for (integer m = 0; m < n; ++m)  x[m] += alpha * u[0][m];



		}


		for (integer j = 1; j <= L; ++j) {
			for (integer i = 1; i < j; ++i) {
				integer ij = (j - 1)*L + (i - 1);
				tau[ij] = Scal(r[j], r[i], n) / sigma[i];

				for (integer m = 0; m < n; ++m)  r[j][m] -= tau[ij] * r[i][m];
			}
			sigma[j] = Scal(r[j], r[j], n);
			gamma_p[j] = Scal(r[0], r[j], n) / sigma[j];
		}



		omega = gamma[L] = gamma_p[L];
		for (integer j = L - 1; j >= 1; --j) {
			gamma[j] = gamma_p[j];
			for (integer i = j + 1; i <= L; ++i)
				gamma[j] -= tau[(i - 1)*L + (j - 1)] * gamma[i];
		}
		for (integer j = 1; j < L; ++j) {
			gamma_pp[j] = gamma[j + 1];
			for (integer i = j + 1; i < L; ++i)
				gamma_pp[j] += tau[(i - 1)*L + (j - 1)] * gamma[i + 1];
		}


		for (integer m = 0; m < n; ++m) x[m] += gamma[1] * r[0][m];
		for (integer m = 0; m < n; ++m) r[0][m] -= gamma_p[L] * r[L][m];
		for (integer m = 0; m < n; ++m) u[0][m] -= gamma[L] * u[L][m];

		for (integer j = 1; j < L; ++j) { /* TODO: use blas DGEMV (for L > 2) */

			for (integer m = 0; m < n; ++m) x[m] += gamma_pp[j] * r[j][m];
			for (integer m = 0; m < n; ++m)  r[0][m] -= gamma_p[j] * r[j][m];
			for (integer m = 0; m < n; ++m)  u[0][m] -= gamma[j] * u[j][m];
		}

		if (iter == maxit) { ierr = 1; break; }

		//printf("%d %e %e\n",iter, NormaV_for_gmres(r[0], n),  bnrm);
		//system("pause");

	}

	printf("final residual = %e\n", NormaV_for_gmres(r[0], n) / bnrm);
	//system("pause");

//finish:

	//delete[] val;
	//delete[] col_ind;
	//delete[] row_ptr;

	if (iVar == TEMP) {
		if (mstruct.bsignalfreeCRSt) {
			// Это таже CRS матрица что и a,ja, ia только элементы в ней нумеруются также как и индексируются с нуля.
			if (mstruct.tval != nullptr) delete mstruct.tval;
			if (mstruct.tcol_ind != nullptr) delete mstruct.tcol_ind;
			if (mstruct.trow_ptr != nullptr) delete mstruct.trow_ptr;

			if (mstruct.ta != nullptr) delete mstruct.ta;
			if (mstruct.tja != nullptr) delete mstruct.tja;
			if (mstruct.tia != nullptr) delete mstruct.tia; // уничтожаем матрицу в CRS формате.

			if (mstruct.ty != nullptr) delete mstruct.ty;

			// alu, jlu - MSR матрица неполного ILU разложения.
			// ju - указатель на диагональные элементы, iw - вспомогательный вектор.
			if (mstruct.talu != nullptr) delete mstruct.talu;
			if (mstruct.tjlu != nullptr) delete mstruct.tjlu;
			if (mstruct.tju != nullptr) delete mstruct.tju;
if (itype_ilu == ILU0)
			{
				if (mstruct.tiw != nullptr) delete mstruct.tiw; // удаляем рабочий массив.
			}

			// Освобождение памяти.
if (itype_ilu == ILU_lfil)
			{
				if (mstruct.tw != nullptr) delete mstruct.tw;
				if (mstruct.tjw != nullptr) delete mstruct.tjw;
				if (mstruct.tlevs != nullptr) delete mstruct.tlevs;
			}

			mstruct.bsignalfreeCRSt = false; // память корректно освобождена.
		}
	}

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (mstruct.bsignalfreeCRScfd) {
			// Это таже CRS матрица что и a,ja, ia только элементы в ней нумеруются также как и индексируются с нуля.
			if (mstruct.val != nullptr) delete mstruct.val;
			if (mstruct.col_ind != nullptr) delete mstruct.col_ind;
			if (mstruct.row_ptr != nullptr) delete mstruct.row_ptr;

			if (mstruct.a != nullptr) delete mstruct.a;
			if (mstruct.ja != nullptr) delete mstruct.ja;
			if (mstruct.ia != nullptr) delete mstruct.ia; // уничтожаем матрицу в CRS формате.

			if (mstruct.y != nullptr) delete mstruct.y;

			// alu, jlu - MSR матрица неполного ILU разложения.
			// ju - указатель на диагональные элементы, iw - вспомогательный вектор.
			if (mstruct.alu != nullptr) delete mstruct.alu;
			if (mstruct.jlu != nullptr) delete mstruct.jlu;
			if (mstruct.ju != nullptr) delete mstruct.ju;
if (itype_ilu == ILU0)
			{
				if (mstruct.iw != nullptr) delete mstruct.iw; // удаляем рабочий массив.
			}

			// Освобождение памяти.
if (itype_ilu == ILU_lfil) 
			{
				if (mstruct.w != nullptr) delete mstruct.w;
				if (mstruct.jw != nullptr) delete mstruct.jw;
				if (mstruct.levs != nullptr) delete mstruct.levs;
			}

			mstruct.bsignalfreeCRScfd = false; // память корректно освобождена.
		}
	}

	delete[] sigma;
	delete[] tau;
	delete[] gamma_pp;
	delete[] gamma_p;
	delete[] gamma;
	for (integer i = 0; i <= L; ++i) {
		delete[] r[i];
		delete[] u[i];
	}
	delete[] r;
	delete[] u;


	for (integer m = 0; m<n; m++) {
		dX0[m] = x[m];

	}

	delete[] x;

}

// Внимание алгоритму присуща стагнация. Его следует использовать только в качестве сглаживателя.
integer  gmres_internal1(equation3D* &sl, equation3D_bon* &slb,
	int maxelm, int maxbound,
	doublereal *dV, doublereal* &dX0, integer maxit, doublereal alpharelax,
	bool bprintmessage, integer iVar, integer &m_restart, 
	int* &ifrontregulationgl, int* &ibackregulationgl,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls)
{
	integer n = maxelm + maxbound;

	// Разреженная матрица СЛАУ
	// в CRS формате.

	doublereal* val = nullptr;
	integer* col_ind = nullptr, *row_ptr = nullptr;

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (ibackregulationgl != nullptr) {
			// nested desection версия алгоритма.
			integer ierr = equation3DtoCRSnd(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
		else {
			integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
	}
	if (iVar == TEMP) {
		integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
		if (ierr > 0) {
			printf("Temperature equation problem.\n");
		}
	}

	// GMRES Саад и Шульц. [1986]
	printf("Y.Saad and Shulc [1986]. General Minimum Residual Method (GMRES).\n");
	gmres(n, val, col_ind, row_ptr, dV, dX0, maxit, m_restart);

	delete[] val;
	delete[] col_ind;
	delete[] row_ptr;

	return 1;
}


// Если использовать gmres с предобуславливателем AINV Bridson то нужно брать m_restart=32. 
// Так рекомендуют в статье для решения уравнений теплопередачи с числом неизвестных до 1.5М.
// Если использовать GMRES в качестве предобуславливателя то рекомендуют брать m_restart=2.
// Даный метод существенно замедляется при увеличении числа итераций. Например он не сходится на Диффузионно дрейфовой модели для 
// материала GaN в то время как BiCGStab сходится.
// Метод использует ilu2 предобуславливатель.
integer  gmres_internal2_stable(equation3D* &sl, equation3D_bon* &slb,
	int maxelm, int maxbound, doublereal *dV, doublereal* &dX0,
	integer maxit, integer &m_restart, doublereal alpharelax, bool bprintmessage, integer iVar,
	QuickMemVorst& m, int* &ifrontregulationgl, int* &ibackregulationgl,
	BLOCK* &b, int &lb, SOURCE* &s_loc, int &ls) {

	integer i_1 = 0;

	// Мы используем m из BiCGStab_internal3 для хранения матрицы предобуславлдивания.
	//bool brc = false;
	//bool bpam_gsp = false; // предобуславливание с помощью нескольких итераций метода Гаусса-Зейделя.

	doublereal *val = nullptr;
	integer* col_ind = nullptr;
	integer* row_ptr = nullptr;
	integer n = maxelm + maxbound;

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (ibackregulationgl != nullptr) {
			// nested desection версия алгоритма.
			integer ierr = equation3DtoCRSnd(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
		else {
			integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
			if (ierr > 0) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX equation problem.\n"); break;
				case VELOCITY_Y_COMPONENT: printf("VY equation problem.\n"); break;
				case VELOCITY_Z_COMPONENT: printf("VZ equation problem.\n"); break;
				case PAM: printf("PAM equation problem.\n"); break;
				}
			}
		}
	}
	if (iVar == TEMP) {
		integer ierr = equation3DtoCRS(sl, slb, val, col_ind, row_ptr, maxelm, maxbound, alpharelax, true, b, lb, s_loc, ls);
		if (ierr > 0) {
			printf("Temperature equation problem.\n");
		}
	}

	const integer ILU0 = 0;
	const integer ILU_lfil = 1;

	const integer itype_ilu = ILU_lfil;//ILU_lfil;




								 // Исходная матрица.
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (!m.ballocCRScfd) {
			// m.a=new doublereal[7*n+2]; // CRS
			// m.ja=new integer[7*n+2];
			// 26 сентября 2016.
			m.a = new doublereal[row_ptr[n] + 2 * maxbound + 2];
			m.ja = new integer[row_ptr[n] + 2 * maxbound + 2];
			m.ia = new integer[n + 2];
			// Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((m.a == nullptr) || (m.ja == nullptr) || (m.ia == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}
	}
	if (iVar == TEMP) {
		if (!m.ballocCRSt) {
			// m.ta=new doublereal[7*n+2]; // CRS
			// m.tja=new integer[7*n+2];
			// 26 сентября 2016.
			m.ta = new doublereal[row_ptr[n] + 2 * maxbound + 2];
			m.tja = new integer[row_ptr[n] + 2 * maxbound + 2];
			m.tia = new integer[n + 2];
			// Флаг что память выделена ставить ещё рано, требуется ещё выделить память под матрицу ILU разложения.
			if ((m.ta == nullptr) || (m.tja == nullptr) || (m.tia == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}
	}
	// Эти структуры данных используются в библиотеке Ю.Саада SPARSKIT2.
	// Основной момент который следует уяснить это:
	// мы используем индексацию массивов в СИ начиная с нуля и до size не включая size;
	// В библиотеке Sparskit2 нумерация элементов массива начинается с единицы и до size включая size.
	// Главное что следует понять не следует пытаться переписать SPARSKIT2 так чтобы она соответствовала нумерации с нуля.
	// оставим нумерацию с единицы иначе мы можем запутаться в логике отлаженного кода SPARSKIT2.
	// Но мы в проекте AliceFlowv0_07 работаем с нумерации начинающейся с нуля. 
	// Код Sparskit2 работает только для предобуславливания. Иными словами схема следующая:
	// На входе матрица в CRS формате с нумерацией с нуля. Преобразуем её в матрицу в которой элементы нумеруются с единицы,
	// для этого очевидно нужно к значениям элементов  col_ind и row_ptr прибавить единицу. Индексация же в матрице (преобразованной пусть начинается с нуля)
	// для этого мы внутри кода SPARSKIT2 декременируем ссылку на массив тогда элемент с номером ноль станет первым тоесть будет иметь номер ноль. Потом в конце использования
	// кода SPARSKIT2 мы увеливаем ссылки на массивы на единицу (откат назад). Получая на вход данную преобразованную матрицу CRS формата, мы делаем с помощью неё матрицу предобуславливателя
	// в MSR формате, как обычно используя код SPARSKIT2 без каких либо изменений из вне (он получен с помощью f2c.exe). Используя матрицу предобуславливания в MSR формате
	// мы пихаем её в функию lusol_ из SPARSKIT2 и получаем необходимый вектор x: (LU)x=y; по вектору y и матрице LU в формате MSR. Код в lusol_ содержит преобразование декремента
	// указателей x, y что позволяет использовать обычные векторы x , y в которых нумерация начинается с нуля. В общем x,y менять не надо, они такие же как обычно в AliceFlowv0_07.
	// В lusol_ указатели сначала декримируются --a; а в конце инкремируются ++a; Так что можно использовать код Sparskit2 без изменений !!!



	if (bprintmessage) {
		printf("Incoplete LU Decomposition begin...\n");
	}


	integer ierr = 0;
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		for (integer i = 0; i<row_ptr[n]; ++i) {
			m.a[i] = val[i];
			m.ja[i] = col_ind[i] + 1;
		}
		for (integer i = 0; i<n + 1; ++i) {
			m.ia[i] = row_ptr[i] + 1;
		}
	}
	if (iVar == TEMP) {
		for (integer i = 0; i<row_ptr[n]; ++i) {
			m.ta[i] = val[i];
			m.tja[i] = col_ind[i] + 1;
		}
		for (integer i = 0; i<n + 1; ++i) {
			m.tia[i] = row_ptr[i] + 1;
		}
	}

	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
		if (!m.ballocCRScfd) {
			//m.ri = new doublereal[n]; m.roc = new doublereal[n]; m.s = new doublereal[n]; m.t = new doublereal[n]; m.vec = new doublereal[n];
			//m.vi = new doublereal[n]; m.pi = new doublereal[n]; m.dx = new doublereal[n]; m.dax = new doublereal[n];
			m.y = new doublereal[n];// m.z = new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
									/*
									if ((m.ri == nullptr) || (m.roc == nullptr) || (m.s == nullptr) || (m.t == nullptr) || (m.vi == nullptr) || (m.pi == nullptr) || (m.dx == nullptr) || (m.dax == nullptr) || (m.y == nullptr) || (m.z == nullptr)) {
									// недостаточно памяти на данном оборудовании.
									printf("Problem: not enough memory on your equipment...\n");
									printf("Please any key to exit...\n");
									exit(1);
									}
									*/
		}
	}
	if (iVar == TEMP) {
		if (!m.ballocCRSt) {
			//m.tri = new doublereal[n]; m.troc = new doublereal[n]; m.ts = new doublereal[n]; m.tt = new doublereal[n];
			//m.tvi = new doublereal[n]; m.tpi = new doublereal[n]; m.tdx = new doublereal[n]; m.tdax = new doublereal[n];
			m.ty = new doublereal[n]; //m.tz = new doublereal[n]; // выделение оперативной памяти для результатов предобуславливания
									  /*
									  if ((m.tri == nullptr) || (m.troc == nullptr) || (m.ts == nullptr) || (m.tt == nullptr) || (m.tvi == nullptr) || (m.tpi == nullptr) || (m.tdx == nullptr) || (m.tdax == nullptr) || (m.ty == nullptr) || (m.tz == nullptr)) {
									  // недостаточно памяти на данном оборудовании.
									  printf("Problem: not enough memory on your equipment...\n");
									  printf("Please any key to exit...\n");
									  exit(1);
									  }
									  */
		}
	}

if (itype_ilu == ILU0)
	{

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {

			if (!m.ballocCRScfd) {
				//m.alu=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.jlu=new integer[7*n+2];
				// 26 сентября 2016.
				m.alu = new doublereal[row_ptr[n] + 2 * maxbound + 2];
				m.jlu = new integer[row_ptr[n] + 2 * maxbound + 2];

				m.ju = new integer[n + 2];
				if (ibackregulationgl != nullptr) {
					// m.alu1=new doublereal[7*n+2]; // +2 запас по памяти.
					// m.jlu1=new integer[7*n+2];
					// m.ju1=new integer[n+2];
					m.x1 = new doublereal[n + 2];
				}
				//m.alurc=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.jlurc=new integer[7*n+2];
				// 26 сентября 2016.
				m.alurc = new doublereal[row_ptr[n] + 2 * maxbound + 2];
				m.jlurc = new integer[row_ptr[n] + 2 * maxbound + 2];

				m.jurc = new integer[n + 2];
				m.iw = new integer[n + 2]; // рабочий массив.
				m.ballocCRScfd = true; // память выделена.

				if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.ju == nullptr) || (m.iw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
				if ((m.alu1 == nullptr) || (m.jlu1 == nullptr) || (m.ju1 == nullptr) || (m.x1 == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}
		if (iVar == TEMP) {
			if (!m.ballocCRSt) {
				//m.talu=new doublereal[7*n+2]; // +2 запас по памяти.
				//m.tjlu=new integer[7*n+2];
				// 26 сентября 2016.
				m.talu = new doublereal[row_ptr[n] + 2 * maxbound + 2];
				m.tjlu = new integer[row_ptr[n] + 2 * maxbound + 2];

				m.tju = new integer[n + 2];
				m.tiw = new integer[n + 2]; // рабочий массив.
				m.ballocCRSt = true; // память выделена.

				if ((m.talu == nullptr) || (m.tjlu == nullptr) || (m.tju == nullptr) || (m.tiw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
			ilu0_(n, m.a, m.ja, m.ia, m.alu, m.jlu, m.ju, m.iw, ierr);
			/* if (ibackregulationgl!=nullptr) {
			for (integer i87=0; i87<7*n+2; i87++) {
			m.alu1[i87]= m.alu[i87];
			m.jlu1[i87]=m.jlu[i87];
			}
			for (integer i87=0; i87<n+2; i87++) {
			m.ju1[i87]=m.ju[i87];
			}
			}*/
		}
		if (iVar == TEMP) {
			ilu0_(n, m.ta, m.tja, m.tia, m.talu, m.tjlu, m.tju, m.tiw, ierr);
		}

		if (ierr>0) {
#if doubleintprecision == 1
			printf("%lld string in matrix is zero diagonal element...\n", ierr - 1);
#else
			printf("%d string in matrix is zero diagonal element...\n", ierr - 1);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}
	}


if (itype_ilu == ILU_lfil)
	{

		//bool btemp_quick = m.ballocCRSt;

		// Сделал параметр lfil передаваемым из интерфейса пользователя.
		integer lfil = my_amg_manager.lfil;// 3; // 2 уровня (0, 1, 2)
		

		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
			if (!m.ballocCRScfd) {

				// инициализация.
				m.alu = nullptr;
				m.jlu = nullptr;
				m.ju = nullptr;
				m.alu1 = nullptr;
				m.jlu1 = nullptr;
				m.ju1 = nullptr;
				m.x1 = nullptr;
				m.alurc = nullptr;
				m.jlurc = nullptr;
				m.jurc = nullptr;
				m.levs = nullptr;
				m.w = nullptr;
				m.jw = nullptr;
				m.w_dubl = nullptr;
				m.jw_dubl = nullptr;

				//m.iwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
				// 26 сентября 2016.
				if (lfil <= 2) {
					m.iwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil == 3) {
					m.iwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if ((lfil >= 4) && (lfil <= 5)) {
					m.iwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil >= 6) {
					m.iwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}

				m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
				m.jlu = new integer[m.iwk + 2];
				m.ju = new integer[n + 2];
				if (ibackregulationgl != nullptr) {
					//m.alu1=new doublereal[m.iwk+2]; // +2 запас по памяти.
					//m.jlu1=new integer[m.iwk+2];
					//m.ju1=new integer[n+2];
					m.x1 = new doublereal[n + 2];
				}
				m.alurc = new doublereal[m.iwk + 2]; // +2 запас по памяти.
				m.jlurc = new integer[m.iwk + 2];
				m.jurc = new integer[n + 2];
				m.levs = new integer[m.iwk + 2]; // уровень.
				m.w = new doublereal[n + 2]; // +2 запас по памяти.
				m.w_dubl = new doublereal[n + 2]; // +2 запас по памяти.

				if (lfil <= 2) {
					m.jw = new integer[3 * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[3 * n + 2]; // +2 запас по памяти.
				}
				else if (lfil == 3) {
					m.jw = new integer[3 * lfil * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[3 * lfil * n + 2]; // +2 запас по памяти.
				}
				else if ((lfil >= 4) && (lfil <= 5)) {
					m.jw = new integer[4 * lfil * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[4 * lfil * n + 2]; // +2 запас по памяти.
				}
				else if (lfil >= 6) {
					m.jw = new integer[5 * lfil * n + 2]; // +2 запас по памяти.				
					m.jw_dubl = new integer[5 * lfil * n + 2]; // +2 запас по памяти.
				}

				m.ballocCRScfd = true; // память выделена.

				if ((m.alu == nullptr) || (m.jlu == nullptr) || (m.levs == nullptr) || (m.ju == nullptr) || (m.w == nullptr) || (m.jw == nullptr) || (m.w_dubl == nullptr) || (m.jw_dubl == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

			}
		}
		if (iVar == TEMP) {
			if (!m.ballocCRSt) {

				// инициализация.
				m.talu = nullptr;
				m.tjlu = nullptr;
				m.tju = nullptr;
				m.tlevs = nullptr;
				m.tw = nullptr;
				m.tjw = nullptr;


				//m.tiwk=(lfil+1)*7*n+4*n; // размерность памяти под матрицу предобуславливания.
				// 26 сентября 2016.
				if (lfil <= 2) {
					m.tiwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil == 3) {
					m.tiwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if ((lfil >= 4) && (lfil <= 5)) {
					m.tiwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}
				else if (lfil >= 6) {
					m.tiwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + 4 * n; // размерность памяти под матрицу предобуславливания.
				}

				m.talu = new doublereal[m.tiwk + 2]; // +2 запас по памяти.
				m.tjlu = new integer[m.tiwk + 2];
				m.tju = new integer[n + 2];
				m.tlevs = new integer[m.tiwk + 2]; // уровень.
				m.tw = new doublereal[n + 2]; // +2 запас по памяти.
				if (lfil <= 2) {
					m.tjw = new integer[3 * n + 2]; // +2 запас по памяти.
				}
				else if (lfil == 3) {
					m.tjw = new integer[3 * lfil* n + 2]; // +2 запас по памяти.
				}
				else if ((lfil >= 4) && (lfil <= 5)) {
					m.tjw = new integer[4 * lfil* n + 2]; // +2 запас по памяти.
				}
				else if (lfil >= 6) {
					m.tjw = new integer[5 * lfil* n + 2]; // +2 запас по памяти.
				}


				m.ballocCRSt = true; // память выделена.

				if ((m.talu == nullptr) || (m.tjlu == nullptr) || (m.tlevs == nullptr) || (m.tju == nullptr) || (m.tw == nullptr) || (m.tjw == nullptr)) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
		}


		if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == PAM)) {
			// iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
			iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);

			if ((ierr == -2) || (ierr == -3)) {

				integer ipassage = 1; // 4 января 2016.
				do {
					printf("\nPlease WAIT... ... ...\n");

					// задаче не хватило памяти, значит нужно перевыделить !
					if (m.alu != nullptr) delete m.alu;
					if (m.jlu != nullptr) delete m.jlu;
					/* if (ibackregulationgl!=nullptr) {
					if (m.alu1!=nullptr) delete m.alu1;
					if (m.jlu1!=nullptr) delete m.jlu1;
					}*/
					if (m.alurc != nullptr) delete m.alurc;
					if (m.jlurc != nullptr) delete m.jlurc;
					if (m.levs != nullptr) delete m.levs;

					// инициализация !
					m.alu = nullptr;
					m.jlu = nullptr;
					/* if (ibackregulationgl!=nullptr) {
					m.alu1=nullptr;
					m.jlu1=nullptr;
					}*/
					m.levs = nullptr;

					// инициализация !
					m.alurc = nullptr;
					m.jlurc = nullptr;

					//m.iwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
					// 26 сентября 2016.
					if (lfil <= 2) {
						m.iwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil == 3) {
						m.iwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if ((lfil >= 4) && (lfil <= 5)) {
						m.iwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil >= 6) {
						m.iwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}


					m.alu = new doublereal[m.iwk + 2]; // +2 запас по памяти.
					m.jlu = new integer[m.iwk + 2];
					/* (ibackregulationgl!=nullptr) {
					m.alu1=new doublereal[m.iwk+2]; // +2 запас по памяти.
					m.jlu1=new integer[m.iwk+2];
					}*/
					m.levs = new integer[m.iwk + 2]; // уровень.

					if ((m.alu != nullptr) && (m.jlu != nullptr) && (m.levs != nullptr)) {
						// iluk_(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, ierr);
						iluk_2(n, m.a, m.ja, m.ia, lfil, m.alu, m.jlu, m.ju, m.levs, m.iwk, m.w, m.jw, m.w_dubl, m.jw_dubl, ierr);
						/*
						if (ibackregulationgl!=nullptr) {
						for (integer i87=0; i87<m.iwk+2; i87++) {
						m.alu1[i87]= m.alu[i87];
						m.jlu1[i87]=m.jlu[i87];
						}
						for (integer i87=0; i87<n+2; i87++) {
						m.ju1[i87]=m.ju[i87];
						}
						}*/

					}
					else {
						// недостаточно памяти на данном оборудовании.
						ipassage = 4;
						printf("Problem: not enough memory on your equipment...\n");
						printf("Please any key to exit...\n");
						exit(1);

					}

					ipassage++;
				} while ((ierr != 0) && (ipassage<4));

				if (ipassage == 4) {
					printf("Error memory alloc !!!\n");
					printf("failed to obtain an expansion for the 4 approaches...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
			else {
				/*
				if (ibackregulationgl!=nullptr) {
				for (integer i87=0; i87<m.iwk+2; i87++) {
				m.alu1[i87]= m.alu[i87];
				m.jlu1[i87]=m.jlu[i87];
				}
				for (integer i87=0; i87<n+2; i87++) {
				m.ju1[i87]=m.ju[i87];
				}
				}*/
			}

		}
		else if (iVar == TEMP) {

			/*
			if (0&&bglobal_unsteady_temperature_determinant) {
			// модификация 20_10_2016.
			// При нестационарном моделировании в твёрдом теле будем строить предобуславливатель лишь единожды на первом шаге.
			//if (!btemp_quick) {
			// iluk_call speed_up
			// 10%=2 9%
			// 15%=3 9%  19.25
			// 20%=4 9%
			// 30%=6 3%
			if (rand()%20<3) {
			iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
			}
			}
			else {
			*/
			iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
			//}

			if ((ierr == -2) || (ierr == -3)) {

				integer ipassage = 1;
				do {
					printf("\nPlease WAIT... ... ...\n");

					// задаче не хватило памяти, значит нужно перевыделить !
					if (m.talu != nullptr) delete m.talu;
					if (m.tjlu != nullptr) delete m.tjlu;
					if (m.tlevs != nullptr) delete m.tlevs;

					// инициализация !
					m.talu = nullptr;
					m.tjlu = nullptr;
					m.tlevs = nullptr;

					//m.tiwk=(lfil+1)*7*n+((1+3+3*ipassage)*n);
					// 26 сентября 2016.
					if (lfil <= 2) {
						m.tiwk = (lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil == 3) {
						m.tiwk = (2 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if ((lfil >= 4) && (lfil <= 5)) {
						m.tiwk = (3 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}
					else if (lfil >= 6) {
						m.tiwk = (4 * lfil + 1) * (row_ptr[n] + 2 * maxbound + 2) + ((1 + 3 + 3 * ipassage)*n);
					}

					m.talu = new doublereal[m.tiwk + 2]; // +2 запас по памяти.
					m.tjlu = new integer[m.tiwk + 2];
					m.tlevs = new integer[m.tiwk + 2]; // уровень.

					if ((m.talu != nullptr) && (m.tjlu != nullptr) && (m.tlevs != nullptr)) {
						iluk_(n, m.ta, m.tja, m.tia, lfil, m.talu, m.tjlu, m.tju, m.tlevs, m.tiwk, m.tw, m.tjw, ierr);
					}
					else {
						// недостаточно памяти на данном оборудовании.
						ipassage = 4;
						printf("Problem: not enough memory on your equipment...\n");
						printf("Please any key to exit...\n");
						exit(1);

					}

					ipassage++;
				} while ((ierr != 0) && (ipassage<4));

				if (ipassage == 4) {
					printf("Error memory alloc !!!\n");
					printf("failed to obtain an expansion for the 4 approaches...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}
		}



		if (ierr != 0) {
#if doubleintprecision == 1
			printf("error memory in iluk ierr=%lld\n", ierr);
#else
			printf("error memory in iluk ierr=%d\n", ierr);
#endif

			//system("pause");
			system("pause");
			exit(0);
		}
	}


	if (bprintmessage) {
		printf("Incoplete LU Decomposition finish...\n");
	}


	//bool bnorelax = true; // Для уравнения теплопроводности не используется релаксация.


	doublereal resid;
	integer i, j = 1, k;
	//Vector s(m + 1), cs(m + 1), sn(m + 1), w;
	doublereal* w = new doublereal[n];
	doublereal* s = new doublereal[m_restart + 2];
	doublereal* cs = new doublereal[m_restart + 2];
	doublereal* sn = new doublereal[m_restart + 2];

	doublereal *dx = new doublereal[n];
	doublereal *buffer = new doublereal[n];


	// начальное приближение
	// X0 ==
	// под X0 понимается вектор поля температур к примеру.
	if (dX0 == nullptr) {
		dX0 = new doublereal[n];
		for (i = 0; i<n; ++i) {
			dx[i] = 0.0;
			dX0[i] = 0.0;
		}
	}
	else {
		for (i = 0; i<n; ++i) dx[i] = dX0[i];
	}

	//doublereal normb = norm(M.solve(b));
	doublereal normb = 0.0;
	// здесь реализованы все три нормы
	// вообще говоря они все эквивалентны


	// Kbuffer=dV

	// (LU)buffer=dV; 
	/*
	if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
	// Очень важно начинать с нуля иначе не будет сходимости.
	#pragma omp parallel for shared(m) private(i_1) schedule (guided)
	for (i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

	//  9 августа 2015 при внедрении перенумерации узлов nested desection
	if (bpam_gsp && (iVar == PAM)) {
	if (ibackregulationgl != nullptr) {
	PAMGSPnd(sl, slb, m.y, w, maxelm, maxbound, ifrontregulationgl);
	}
	else {
	PAMGSP(sl, slb, m.y, w, maxelm, maxbound);
	}
	}
	else {

	if (brc) {
	for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
	for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
	m.alurc[i7] = m.alu[i7];
	m.jlurc[i7] = m.jlu[i7];
	}
	for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
	}

	if (ibackregulationgl != nullptr) {
	//lusol_2(n, w, m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=w;
	lusol_3(n, w, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=w;
	}
	else {
	lusol_(n, w, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=w;

	}

	if (brc) {
	for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
	for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
	m.alu[i7] = m.alurc[i7];
	m.jlu[i7] = m.jlurc[i7];
	}
	for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
	}

	}
	for (i_1 = 0; i_1 < n; ++i_1) w[i_1] = m.y[i_1];


	}
	*/
	/*
	if (iVar == TEMP) {
	// Очень важно начинать с нуля иначе не будет сходимости.
	#pragma omp parallel for shared(m) private(i) schedule (guided)
	for (i_1 = 0; i_1<n; ++i_1) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

	lusol_(n, dV, m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=w;
	for (i_1 = 0; i_1 < n; ++i_1) buffer[i_1] = m.ty[i_1];

	}
	*/
	normb = NormaV_for_gmres(dV, n);
	//normb = NormaV(buffer, n);

	//Vector r = M.solve(dV - A * x);
	doublereal *r = new doublereal[n];
	MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // результат занесён в  r
	for (i = 0; i < n; ++i) r[i] = dV[i] - r[i];

	//  calculate residual precontidioning;

	/*
	// Ky=r

	// (LU)y=r;
	if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
	// Очень важно начинать с нуля иначе не будет сходимости.
	#pragma omp parallel for shared(m) private(i_1) schedule (guided)
	for (i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

	//  9 августа 2015 при внедрении перенумерации узлов nested desection
	if (bpam_gsp && (iVar == PAM)) {
	if (ibackregulationgl != nullptr) {
	PAMGSPnd(sl, slb, m.y, r, maxelm, maxbound, ifrontregulationgl);
	}
	else {
	PAMGSP(sl, slb, m.y, r, maxelm, maxbound);
	}
	}
	else {

	if (brc) {
	for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
	for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
	m.alurc[i7] = m.alu[i7];
	m.jlurc[i7] = m.jlu[i7];
	}
	for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
	}

	if (ibackregulationgl != nullptr) {
	//lusol_2(n, v[0], m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=r;
	lusol_3(n, r, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=r;
	}
	else {
	lusol_(n, r, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=r;

	}

	if (brc) {
	for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
	for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
	m.alu[i7] = m.alurc[i7];
	m.jlu[i7] = m.jlurc[i7];
	}
	for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
	}

	}
	for (i_1 = 0; i_1 < n; ++i_1) r[i_1] = m.y[i_1];


	}
	if (iVar == TEMP) {
	// Очень важно начинать с нуля иначе не будет сходимости.
	#pragma omp parallel for shared(m) private(i) schedule (guided)
	for (i_1 = 0; i_1<n; ++i_1) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

	lusol_(n, r, m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=r;
	for (i_1 = 0; i_1 < n; ++i_1) r[i_1] = m.ty[i_1];

	}
	*/
	//doublereal beta = norm(r);
	doublereal beta = 0.0;



	beta = NormaV_for_gmres(r, n);

	if (fabs(normb) < 1.0e-30)
		normb = 1;

	doublereal norm_r = 0.0;


	norm_r = NormaV_for_gmres(r, n);

	if ((resid = norm_r / normb) <= dterminatedTResudual) {
		//tol = resid;
		//maxit = 0;
		delete[] w;
		delete[] s;
		delete[] cs;
		delete[] sn;
		delete[] buffer;
		return 0;
	}

	doublereal** H = new doublereal*[m_restart + 2]; // Hessenberg
	for (i_1 = 0; i_1 < m_restart + 2; ++i_1) H[i_1] = new doublereal[m_restart + 2];
	for (i_1 = 0; i_1 < m_restart + 2; ++i_1)
	{
		for (integer j_1 = 0; j_1 < m_restart + 2; ++j_1)
		{
			H[i_1][j_1] = 0.0;
		}
	}

	//Vector *v = new Vector[m_restart + 1];
	doublereal** v = new doublereal*[m_restart + 2];
	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) v[i_1] = new doublereal[n];
	for (i_1 = 0; i_1 <= m_restart + 1; ++i_1) {
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			v[i_1][j_1] = 0.0;
		}
	}

	j = 1; // номер первой итерации
		   //doublereal delta = 1.0e-3;// DOPOLNENIE

	//integer i_copy;

	while (j <= maxit) {

		//v[0] = r * (1.0 / beta);    // ??? r / beta
		for (integer j_1 = 0; j_1 < n; ++j_1)
		{
			v[0][j_1] = r[j_1] * (1.0 / beta);
		}

		//s = 0.0;
		for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) s[i_1] = 0.0;
		s[0] = beta;

		/*
		for (i_1 = 0; i_1 < m_restart + 2; ++i_1)
		{ // DOPOLNENIE
		for (integer j_1 = 0; j_1 < m_restart + 2; ++j_1)
		{
		H[i_1][j_1] = 0.0;
		}
		}
		*/

		// Ортогонализация Арнольди.
		for (i = 0; i < m_restart && j <= maxit; i++, j++) {

			//i_copy = i;

			// закомментировано без предобуславливания.
			//w = M.solve(A * v[i]);
			MatrixCRSByVector(val, col_ind, row_ptr, v[i], w, n); // результат занесён в  w

																  // Kw=A*v[i]

																  // (LU)w=A*v[i];
																  /*
																  if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
																  // Очень важно начинать с нуля иначе не будет сходимости.
																  #pragma omp parallel for shared(m) private(i_1) schedule (guided)
																  for (i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

																  //  9 августа 2015 при внедрении перенумерации узлов nested desection
																  if (bpam_gsp && (iVar == PAM)) {
																  if (ibackregulationgl != nullptr) {
																  PAMGSPnd(sl, slb, m.y, w, maxelm, maxbound, ifrontregulationgl);
																  }
																  else {
																  PAMGSP(sl, slb, m.y, w, maxelm, maxbound);
																  }
																  }
																  else {

																  if (brc) {
																  for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
																  for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
																  m.alurc[i7] = m.alu[i7];
																  m.jlurc[i7] = m.jlu[i7];
																  }
																  for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
																  }

																  if (ibackregulationgl != nullptr) {
																  //lusol_2(n, w, m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=w;
																  lusol_3(n, w, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=w;
																  }
																  else {
																  lusol_(n, w, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=w;

																  }

																  if (brc) {
																  for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
																  for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
																  m.alu[i7] = m.alurc[i7];
																  m.jlu[i7] = m.jlurc[i7];
																  }
																  for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
																  }

																  }
																  //for (i_1 = 0; i_1 < n; ++i_1) w[i_1] = m.y[i_1];
																  for (i_1 = 0; i_1 < n; ++i_1)  v[i + 1][i_1] = m.y[i_1];


																  }
																  */
																  /*
																  if (iVar == TEMP) {
																  // Очень важно начинать с нуля иначе не будет сходимости.
																  #pragma omp parallel for shared(m) private(i) schedule (guided)
																  for (i_1 = 0; i_1<n; ++i_1) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

																  lusol_(n, buffer, m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=A*v[i];
																  //for (i_1 = 0; i_1 < n; ++i_1) w[i_1] = m.ty[i_1];
																  for (i_1 = 0; i_1 < n; ++i_1) v[i + 1][i_1] = m.ty[i_1];

																  }
																  */
																  //doublereal av = sqrt(Scal(w,w,n)); // DOPOLNENIE
																  //doublereal av = sqrt(Scal(v[i + 1], v[i + 1], n)); // DOPOLNENIE

			for (k = 0; k <= i; ++k) {
				H[k][i] = Scal(w, v[k], n);
				//H[k][i] = Scal(v[i + 1], v[k], n);
				for (integer j_1 = 0; j_1 < n; ++j_1)
				{
					//v[i + 1][j_1] -= H[k][i] * v[k][j_1];
					w[j_1] -= H[k][i] * v[k][j_1];
				}
			}
			//H[i + 1][i] = norm(w);
			H[i + 1][i] = NormaV_for_gmres(w, n);
			//H[i + 1][i] = NormaV(v[i + 1], n);

			/*
			// DOPOLNENIE
			if ((av + delta * H[i+1][i]) == av)
			{
			for (k = 0; k <= i; ++k)//j
			{
			//doublereal htmp = Scal(w,v[k],n);
			doublereal htmp = Scal(v[i + 1], v[k], n);
			//htmp = r8vec_dot(n, v + k*n, v + (j - 1)*n);
			H[k][i] = H[k][i] + htmp;
			//h[(j - 1) + (k - 1)*(mr + 1)] = h[(j - 1) + (k - 1)*(mr + 1)] + htmp;
			for (integer j_1 = 0; j_1 < n; ++j_1)
			{
			v[i][j_1] = v[i][j_1] - htmp * v[k][j_1];
			}
			}
			H[i+1][i] = sqrt(Scal( v[i] , v[i],n));
			}

			if (H[i + 1][i] != 0.0) {
			for (integer j_1 = 0; j_1 < n; ++j_1)
			{
			//v[i + 1][j_1] = w[j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			v[i + 1][j_1] = v[i + 1][j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			}
			}
			*/
			for (integer j_1 = 0; j_1 < n; ++j_1)
			{
				v[i + 1][j_1] = w[j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
															  //v[i + 1][j_1] = v[i + 1][j_1] * (1.0 / H[i + 1][i]); // ??? w / H(i+1, i)
			}
			// Окончание ортогонализации Арнольди.
			// В v - хранится ортонормированный базис подпространства Крылова размерности m_restart.
			// H - Верхнетреугольная матрица Хессенберга - матрица коэффициентов ортогонализации.
			/*
			if (0 < i) {
			for (k = 0; k <= i + 1; ++k)
			{
			m.ty[k] = H[k][i];
			}
			for (k = 0; k <= i-1; ++k)
			{
			mult_givens(cs[k], sn[k], k, m.ty);
			}
			for (k = 0; k <= i + 1; ++k)
			{
			H[k][i] = m.ty[k];
			}

			}

			doublereal mu = sqrt(pow(H[i][i], 2)
			+ pow(H[i+1][i], 2));
			cs[i] = H[i][i] / mu;
			sn[i] = -H[i+1][i] / mu;
			H[i][i] = cs[i] * H[i][i] - sn[i] * H[i+1][i];
			H[i+1][i] = 0;
			mult_givens(cs[i], sn[i], i, s);
			*/
			//rho = fabs(s[i+1]);

			/*
			itr_used = itr_used + 1;

			if ( verbose )
			{
			cout << "  K =   " << k << "  Residual = " << rho << "\n";
			}

			if ( rho <= rho_tol && rho <= tol_abs )
			{
			break;
			}
			*/


			for (k = 0; k < i; ++k)
				ApplyPlaneRotation(H[k][i], H[k + 1][i], cs[k], sn[k]);

			GeneratePlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
			ApplyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);


			// Вручную устраняем случай полного совпадения невязок на двух соседних итерациях,
			// т.к. иначе это приводит к развалу решения.
			//if (fabs(s[i] - s[i + 1]) < 1.0e-37) s[i + 1] = 1.05*s[i];

			printf("%lld %e \n", j, fabs(s[i + 1]) / normb);
			system("pause");

			resid = fabs(s[i + 1]) / normb;

			if ((resid) < dterminatedTResudual) {
				Update(dx, i, n, H, s, v);
				//tol = resid;
				//maxit = j;
				for (integer i_1 = 0; i_1<n; ++i_1) {
					dX0[i_1] = dx[i_1];

				}
				for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] v[i_1];
				delete[] v;
				delete[] dx;
				delete[] buffer;
				delete[] r;
				delete[] w;
				delete[] s;
				delete[] cs;
				delete[] sn;
				for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) delete[] H[i_1];
				delete[] H;
				delete[] val;
				delete[] col_ind;
				delete[] row_ptr;
				return 0;
			}
		}



		//Update(dx, m_restart - 1, n, H, s, v);//i-1 //ERROR
		Update(dx, i - 1, n, H, s, v);//i-1 //ERROR
									  //r = M.solve(b - A * x);
		MatrixCRSByVector(val, col_ind, row_ptr, dx, r, n); // Результат занесён в r
		for (integer i_1 = 0; i_1 < n; ++i_1) r[i_1] = dV[i_1] - r[i_1];

		//  calculate residual precontidioning;

		// Ky=r
		/*
		// (LU)y=r;
		if ((iVar == VX) || (iVar == VY) || (iVar == VZ) || (iVar == PAM)) {
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(m) private(i_1) schedule (guided)
		for (i_1 = 0; i_1<n; ++i_1) m.y[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для PAM !.

		//  9 августа 2015 при внедрении перенумерации узлов nested desection
		if (bpam_gsp && (iVar == PAM)) {
		if (ibackregulationgl != nullptr) {
		PAMGSPnd(sl, slb, m.y, r, maxelm, maxbound, ifrontregulationgl);
		}
		else {
		PAMGSP(sl, slb, m.y, r, maxelm, maxbound);
		}
		}
		else {

		if (brc) {
		for (integer i7 = 0; i7<n; ++i7) m.vec[i7] = m.pi[i7];
		for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
		m.alurc[i7] = m.alu[i7];
		m.jlurc[i7] = m.jlu[i7];
		}
		for (integer i7 = 0; i7<n + 2; ++i7) m.jurc[i7] = m.ju[i7];
		}

		if (ibackregulationgl != nullptr) {
		//lusol_2(n, v[0], m.y, m.alu, m.jlu, m.ju, m.x1, maxelm); // M*y=r;
		lusol_3(n, r, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=r;
		}
		else {
		lusol_(n, r, m.y, m.alu, m.jlu, m.ju, maxelm); // M*y=r;

		}

		if (brc) {
		for (integer i7 = 0; i7<n; ++i7) m.pi[i7] = m.vec[i7];
		for (integer i7 = 0; i7<m.iwk + 2; ++i7) {
		m.alu[i7] = m.alurc[i7];
		m.jlu[i7] = m.jlurc[i7];
		}
		for (integer i7 = 0; i7<n + 2; ++i7) m.ju[i7] = m.jurc[i7];
		}

		}
		for (i_1 = 0; i_1 < n; ++i_1) r[i_1] = m.y[i_1];


		}
		if (iVar == TEMP) {
		// Очень важно начинать с нуля иначе не будет сходимости.
		#pragma omp parallel for shared(m) private(i) schedule (guided)
		for (i_1 = 0; i_1<n; ++i_1) m.ty[i_1] = 0.0; // Если начинать не с нуля то небудет сходимости для TEMP !.

		lusol_(n, r, m.ty, m.talu, m.tjlu, m.tju, maxelm); // M*ty=r;
		for (i_1 = 0; i_1 < n; ++i_1) r[i_1] = m.ty[i_1];

		}
		*/

		/*
		i = i_copy - 1;
		m.ty[i] = s[i] / H[i][i];

		for (k = i; 0 <= k; k--)
		{
		m.ty[k] = s[k];
		for (j = k + 1; j <= i + 1; ++j)
		{
		m.ty[k] = m.ty[k] -H[k][j] * m.ty[j];
		}
		m.ty[k] = m.ty[k] / H[k][k];
		}

		for (k = 0; k < n; ++k)
		{
		for (j = 1; j <= k + 1; ++j)
		{
		dx[k] = dx[k] + v[j][k] * m.ty[j];
		}
		}
		*/
		/*
		if (rho <= rho_tol && rho <= tol_abs)
		{
		break;
		}
		*/

		//beta = norm(r);
		beta = NormaV_for_gmres(r, n);

		resid = beta / normb;

		if ((resid) < dterminatedTResudual) {
			//tol = resid;
			//maxit = j;
			for (integer i_1 = 0; i_1<n; ++i_1) {
				dX0[i_1] = dx[i_1];

			}
			for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] v[i_1];
			delete[] v;

			delete[] dx;
			delete[] buffer;
			delete[] r;
			delete[] w;
			delete[] s;
			delete[] cs;
			delete[] sn;
			for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) delete[] H[i_1];
			delete[] H;
			delete[] val;
			delete[] col_ind;
			delete[] row_ptr;
			return 0;
		}
	}

	//tol = resid;
	for (integer i_1 = 0; i_1<n; ++i_1) {
		dX0[i_1] = dx[i_1];

	}
	for (integer i_1 = 0; i_1 <= m_restart + 1; ++i_1) delete[] v[i_1];
	delete[] v;

	delete[] dx;
	delete[] buffer;
	delete[] r;
	delete[] w;
	delete[] s;
	delete[] cs;
	delete[] sn;
	for (integer i_1 = 0; i_1 < m_restart + 2; ++i_1) delete[] H[i_1];
	delete[] H;
	delete[] val;
	delete[] col_ind;
	delete[] row_ptr;
	return 1;

} //gmres_internal2_stable




  // Детали решают всё. Экспериментально было выяснено, что метод
  // Bi_CGStab_internal1 не сходится должным образом. Это можно объяснить тем
  // что реализация вычислительной процедуры на языке СИ недостаточно правильна.
  // Поэтому 31 марта 2013 года мы попробуем вставить ILU предобуславливание в
  // BiCGStab на подобии алгоритма Lr1sk (см. версию Bi_CGStab_internal2).
void Bi_CGStab(IMatrix *xO, equation3D* &sl, equation3D_bon* &slb,
	int maxelm, int maxbound,
	doublereal *dV, doublereal* &dX0,
	integer maxit, doublereal alpharelax, integer iVar,
	QuickMemVorst& m, bool bLRfree, BLOCK* &b, int &lb, 
	int* &ifrontregulationgl,
	int* &ibackregulationgl,
	doublereal dgx, doublereal dgy, doublereal dgz,
	SOURCE* &s_loc, int &ls, integer inumber_iteration_SIMPLE,
	integer* &color, integer dist_max, WALL* &w, int &lw, 
	int * &whot_is_block, bool bprint_message)
{

	

	if (iVar == TURBULENT_KINETIK_ENERGY) {
		for (integer i_1 = 0; i_1 < maxelm + maxbound; ++i_1) {
			dV[i_1] *= 1.0e6;
		}
	}

	for (integer i_1 = 0; i_1 < maxelm + maxbound; ++i_1) {
		if (dV[i_1] != dV[i_1]) {
			switch (iVar) {
			case VELOCITY_X_COMPONENT: printf("VX rthdsd problem\n");
				break;
			case VELOCITY_Y_COMPONENT: printf("VY rthdsd problem\n");
				break;
			case VELOCITY_Z_COMPONENT: printf("VZ rthdsd problem\n");
				break;
			case PAM: printf("PAM rthdsd problem iP=%lld\n",i_1);
				break;
			case NUSHA: printf("NU rthdsd problem\n");
				break;
			case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY rthdsd problem\n");
				break;
			case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA rthdsd problem\n");
				break;
			case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS rthdsd problem\n");
				break;
			case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS rthdsd problem\n");
				break;
			case RE_THETA_LANGTRY_MENTER: printf("TURBULENT_RE_THETA_LANGTRY_MENTER  rthdsd problem\n"); break;
			case GAMMA_LANGTRY_MENTER: printf("TURBULENT_GAMMA_LANGTRY_MENTER  rthdsd problem\n");	break;
			case TEMP: printf("TEMP rthdsd problem\n");
				break;
			}
			printf("May be NAN or INF in premeshin.txt file. Power in control volume= %lld is undefined...\n", i_1);
			printf("ispolzuite poslednuu versiu Mesh generator AliceMesh. 04.05.2019.\n");
			system("pause");
			exit(1);
		}
	}


	// отладка.
#if doubleintprecision == 1
	//printf("iswitchsolveramg_vs_BiCGstab_plus_ILU2=%lld\n", iswitchsolveramg_vs_BiCGstab_plus_ILU2);
#else
	//printf("iswitchsolveramg_vs_BiCGstab_plus_ILU2=%d\n", iswitchsolveramg_vs_BiCGstab_plus_ILU2);
#endif


	//system("pause");

	// iVar==PAM && bLRfree определяет особый класс задач:
	// каверна и тест Валь Девиса.
	// Для промышлености это различные охладители выделяющие тепло в воздушном 
	// пространстве без контактов с теплоотводом. Это может встречаться как в реальных режимах эксплуатации
	// так и тестировании на тепловую надёжность прибора.

	// Замер времени.
	unsigned int calculation_main_start_time; // начало счёта мс.
	unsigned int calculation_main_end_time; // окончание счёта мс.

	calculation_main_start_time = clock(); // момент начала счёта.

	set_RUMBA_Classic_AMG_Setting(iVar);
	
	if (my_amg_manager.iprint_log == false) {
		print_log_message = false; // Не будем печатать лог на консоль.
	}


	// Метод Ван Дер Ворста Bi-CGStab
	// работает для возможно несимметричных вещественных матриц.
	// встроен предобуславливатель ILU(0).
	// Метод является комбинацией методов BiCG и GMRES(1). 
	// Экспериментально было выяснено что он не рекомендуется к использованию из-за
	// серьёзных проблем, отсутствия сходимости. Это утверждение не относится собственно
	// к алгоритму BiCGStab Х ван дер Ворста, а относится к данной конкретной реализации 
	// на языке СИ данного алгоритма.
	//if (!bBiCGStabSaad) { // см solve
	//Bi_CGStab_internal1(xO, sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax);
	bool bprintmessage = false; // печать невязки на консоль
								// данный метод со встроенным предобуславливателем ILU0 сходится гораздо лучше,
								// чем обычный BiCGStabCRS, но всё-таки хуже чем Lr1sk алгоритм.
								//Bi_CGStab_internal2(xO, sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage);
								//}
								//if (bBiCGStabSaad) { // см solve
								// internal3:
								// освобождение оперативной памяти
								//freeIMatrix(xO);
								// Надо проверить, что будет при встраивании предобуславливателя ILUT.

								//}
								/*
								if (iVar==PAM) {
								//PAMGSP(sl, slb, dX0, dV, maxelm, maxbound);
								ICCG(sparseM, dV, f.potent[PAM], f.maxelm + f.maxbound ,bprintmessage,false,2000);
								}
								else {
								Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,w,lw);
								}*/


	for (integer inumber = 0; inumber < maxelm + maxbound; inumber++) {
		if (dV[inumber] != dV[inumber]) {
			printf("dV!=dV assemble bug. inumber=%lld dV=%e\n", inumber, dV[inumber]);
			system("pause");
		}
		if (dX0[inumber] != dX0[inumber]) {
			printf("dX0!=dX0 assemble bug. inumber=%lld dX0=%e\n", inumber, dX0[inumber]);
			system("pause");
		}
		if (inumber < maxelm) {

		}
		else {
			integer iW = inumber - maxelm;
			if (slb[iW].b != slb[iW].b) {
				printf("slb.b!=slb.b assemble bug. inumber=%lld slb.b=%e\n", iW, slb[iW].b);
				system("pause");
			}
			if (slb[iW].aw != slb[iW].aw) {
				printf("slb.aw!=slb.aw assemble bug. inumber=%lld slb.aw=%e\n", iW, slb[iW].aw);
				system("pause");
			}
			if (slb[iW].ai != slb[iW].ai) {
				printf("slb.ai!=slb.ai assemble bug. inumber=%lld slb.ai=%e\n", iW, slb[iW].ai);
				system("pause");
			}
		}
	}




	if (!bdontstartsolver) {
		if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 0) {

			// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
			// BiCGStab + ILU(k). k=1 or 2 recomended.
			const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
			Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl,b,lb,s_loc,ls, inumber_iteration_SIMPLE, color, dist_max, breordering,w,lw);
			
			//integer L = 2;
			//Bi_CGStab_internal5(L, sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl);

		}
		else if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 1) {

			// здесь предложена реализация алгебраического многосеточного метода
			// под названием amg1r5 предложенная широкой публике в 1985 году.
			// Впервые в мировой истории многосеточный метод появился в статье 
			// Радия Петровича федоренко в 1961 году.

			// дата первого успешного запуска в коде AliceFlow_v0_07:
			// 15 июля 2015 года. Среда. 

			if (AMG1R5_OUT_ITERATOR::NONE_only_amg1r5 == AMG1r5Dat.stabilization_amg1r5_algorithm) {

				if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
					(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
					(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
					(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
					Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
						ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, 
						dist_max, breordering,w,lw);
				}
				else {
					bool worked_successfully = false;
					if ((iVar==PAM) && (inumber_iteration_SIMPLE < 10)) {
						// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
						const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
						Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, 
							ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color,
							dist_max, breordering, w, lw);

					}
					else {
						// amg1r5 realisation.
						amg(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m,
							ifrontregulationgl, ibackregulationgl, 0, worked_successfully, b, lb,
							s_loc, ls, w, lw, whot_is_block, bprint_message);
					}

					if (!bsolid_static_only) {
						if (!worked_successfully) {
							//30.03.2019
							// СБРОС обнуление.
							for (integer i_5 = 0; i_5 < maxelm + maxbound; i_5++) {
								if (i_5 < maxelm) {
									dX0[i_5] = 0.0;
								}
								else {
									if (slb[i_5 - maxelm].iI > -1) {
										// Однородное условие Неймана.
										dX0[i_5] = 0.0;
									}
								}
							}
							// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
							printf("Redirecting to BiCGStab + ILU2 solver.\n");
							const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
							Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
								ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, 
								dist_max, breordering, w, lw);
						}
					}
				}
			}

			if (AMG1R5_OUT_ITERATOR::BiCGStab_plus_amg1r5 == AMG1r5Dat.stabilization_amg1r5_algorithm) {
				// BiCGStab[1992] + amg1r5[1986]
			// Предобуславливание, Многосеточные технологии, Стабилизация.
			// 23-24 декабря 2017.

				if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
					(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
					(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
					(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					const bool breordering = false;
					Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
						ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, 
						color, dist_max, breordering, w, lw);
				}
				else {

					// H.A. VAN DER Vorst, BiCGStab, 1992.
					// Руге и Штубен, 1986.

					bool worked_successfully = false;
					const integer iHAVorstModification_id = 1;
					if ((iVar == PAM) && (inumber_iteration_SIMPLE < 10)) {
						// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
						const bool breordering = false;
						Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
							ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color,
							dist_max, breordering, w, lw);

					}
					else {
						amg(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m, ifrontregulationgl,
							ibackregulationgl, iHAVorstModification_id, worked_successfully, b, lb, s_loc, ls, 
							w, lw, whot_is_block,  bprint_message);
					}

					if (iVar == PAM) {
						if (!worked_successfully) {
							//30.03.2019
							printf("PAM equation divergence detected BiCGStab + amg1r5 solver.\n");
							// СБРОС огбнуление.
							for (integer i_5 = 0; i_5 < maxelm + maxbound; i_5++) {
								if (i_5 < maxelm) {
									dX0[i_5] = 0.0;
								}
								else {
									if (slb[i_5 - maxelm].iI > -1) {
										// Однородное условие Неймана.
										dX0[i_5] = 0.0;
									}
								}
							}
							printf("Redirecting to BiCGStab + ILU2 solver.\n");
							// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
							const bool breordering = false;
							Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
								ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE,
								color, dist_max, breordering,w,lw);
						}
					}

				}
			}

			if (AMG1R5_OUT_ITERATOR::FGMRes_plus_amg1r5 == AMG1r5Dat.stabilization_amg1r5_algorithm) {
				// FGMRes[1986] + amg1r5[1986]
			// Предобуславливание, Многосеточные технологии.
			// 31 декабря 2017.

				if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
					(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
					(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
					(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					const bool breordering = false;
					Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, 
						ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering,w,lw);
				}
				else {

					// Ю.Саад и Шульц, FGMRes, 1986.
					// Руге и Штубен, 1986.
					bool worked_successfully = false;
#ifdef _OPENMP
					omp_set_num_threads(6);
#endif
					const integer iHAVorstModification_id = 2;
					if ((iVar == PAM) && (inumber_iteration_SIMPLE < 10)) {
						// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
						const bool breordering = false;
						Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, 
							ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, 
							color, dist_max, breordering, w, lw);

					}
					else {
						amg(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m, ifrontregulationgl, 
							ibackregulationgl, iHAVorstModification_id, worked_successfully, b, lb, s_loc, ls,
							w, lw, whot_is_block, bprint_message);
					}
#ifdef _OPENMP
					omp_set_num_threads(1);
#endif
					if (iVar == PAM) {
						if (!worked_successfully) {
							//30.03.2019
							printf("PAM equation divergence detected FGMRES + amg1r5 solver.\n");
							// СБРОС огбнуление.
							for (integer i_5 = 0; i_5 < maxelm + maxbound; i_5++) {
								if (i_5 < maxelm) {
									dX0[i_5] = 0.0;
								}
								else {
									if (slb[i_5 - maxelm].iI > -1) {
										// Однородное условие Неймана.
										dX0[i_5] = 0.0;
									}
								}
							}
							printf("Redirecting to BiCGStab + ILU2 solver.\n");
							// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
							const bool breordering = false;
							Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering,w,lw);
						}
					}

				}
			}


			if (AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5 == AMG1r5Dat.stabilization_amg1r5_algorithm) {
				// BiCGStab[1992] + amg1r5[1986]
				// Предобуславливание, Многосеточные технологии, Стабилизация.
				// 23-24 декабря 2017.

				if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
					(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
					(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
					(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					const bool breordering = false;
					Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering,w,lw);
				}
				else {

					
					if (iVar == PAM) {

						// H.A. VAN DER Vorst, BiCGStab, 1992.
						// Руге и Штубен, 1986.

						bool worked_successfully = false;
						const integer iHAVorstModification_id = 1; // BiCGStab  for PAM
						amg(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m,
							ifrontregulationgl, ibackregulationgl, iHAVorstModification_id,
							worked_successfully, b, lb, s_loc, ls, w, lw, whot_is_block, bprint_message);


						if (!worked_successfully) {
							//30.03.2019
							printf("PAM equation divergence detected BiCGStab + amg1r5 solver.\n");
							// СБРОС огбнуление.
							for (integer i_5 = 0; i_5 < maxelm + maxbound; i_5++) {
								if (i_5 < maxelm) {
									dX0[i_5] = 0.0;
								}
								else {
									if (slb[i_5 - maxelm].iI > -1) {
										// Однородное условие Неймана.
										dX0[i_5] = 0.0;
									}
								}
							}
							printf("Redirecting to BiCGStab + ILU2 solver.\n");
							// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
							const bool breordering = false;
							Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);
						}
					}
					else {
						// H.A. VAN DER Vorst, BiCGStab, 1992.
						// Руге и Штубен, 1986.

						bool worked_successfully = false;
						const integer iHAVorstModification_id = 3; // non linear solver
						amg(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m, 
							ifrontregulationgl, ibackregulationgl, iHAVorstModification_id, 
							worked_successfully, b, lb, s_loc, ls, w, lw, whot_is_block, bprint_message);


					}

				}
			}


			//std::cout << "statistics...\n";
			//system("pause");

		}
		else if (2 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) {
			// LR1sK
			printf("ERROR !!! Call Lr1sk should be earlier in solver mysolverv0_03.c source code file.\n");
			printf("varialable is equal ");
			switch (iVar) {
			case VELOCITY_X_COMPONENT: printf("Vx \n");  break;
			case VELOCITY_Y_COMPONENT: printf("Vy \n");  break;
			case VELOCITY_Z_COMPONENT: printf("Vz \n");  break;
			case NUSHA: printf("NU \n");  break;
			case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY \n");  break;
			case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA \n");  break;
			case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS \n"); break;
			case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS \n"); break;
			case RE_THETA_LANGTRY_MENTER: printf("TURBULENT_RE_THETA_LANGTRY_MENTER \n"); break;
			case GAMMA_LANGTRY_MENTER: printf("TURBULENT_GAMMA_LANGTRY_MENTER \n");	break;
			case PAM: printf("PAM \n");  break;
			case TEMP: printf("TEMP \n"); break;
			}
			printf("Redirecting to BiCGStab + ILU2 solver.\n");
			system("PAUSE");
			const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
			Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);
			
			//system("pause");
			//system("PAUSE");
			//exit(1);
		}
		else if (4 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) {
			/*
	#if GPU_LIB_INCLUDE_MY_PROJECT == 1
				// Этот метод заимствован из библиотеки CUSP 0.5.1 распространяемой по
				// OpenSource Apache license 2.0.
				// В данном случае на многих ядрах видеокарты используется алгоритм
				// BiCGStab Хенка Ван дер Ворста и AINV (NS Brigson) в качестве предобуславливателя
				// GPU Accelerating FREE now!!!. Чем мощнее ваша видюха тем больше выигрыш в скорости вычисления.
				cusp_solver_GPU_AINV_Bridson(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar); // На GPU!!!

	#else
				printf("WARNING: CUSP 0.5.1 library is not connected\n");
				printf("Redirecting to BiCGStab + ILU2 solver.\n");
				// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
				const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
				Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering,w,lw);
	#endif
	*/
	// Будем использовать алгоритм ddemidova на CPU. 17.07.2021
	// Теперь в программе одновременно доступно вычисление как на CPU так и на GPU с помощью библиотеки amgcl ddemidova.
			if ((iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {

				// RE_THETA_LANGTRY_MENTER можно решать amgcl но будет медленней.

				// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
				const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
				Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);
			}
			else if (0 && ((iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA))) {
				// Обнаружена проблема работоспособности amgcl на уравнении для K в K-Omega модели турбулентности. 11,10,2019
				// Переключаемся на алгоритм Юсефа Саада.
				// в amgcl отключил проверку в bicgstab для rhs.

				// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
				const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
				Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);
			}
			else if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) {

				// 16.01.2022 Обнаружено что spai0-smoothed aggregation FGMRES(600) дал сбой на компонентах скорости. Решено перенаправить на стандартный BiCGStab+ILU(lfil).

				// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
				const bool breordering = false;
				Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
					ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);
			}
			else {
#ifdef AMGCL_INCLUDE_IN_MY_PROJECT
				const bool bprint_preconditioner_amgcl = bprint_message; // true;
				amgcl_solver_cpu(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar, bprint_preconditioner_amgcl, dgx, dgy, dgz, inumber_iteration_SIMPLE, w, lw, bprint_message);
#endif
		}

		}
		else if (11 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) {

#if GPU_LIB_INCLUDE_MY_PROJECT == 1
		// Этот метод заимствован из библиотеки CUSP 0.5.1 распространяемой по
		// OpenSource Apache license 2.0.
		// В данном случае на одном ядре центрального процессора используется алгоритм
		// BiCGStab Хенка Ван дер Ворста и AINV (NS Brigson) в качестве предобуславливателя 

		if (bglobal_unsteady_temperature_determinant) {
			cusp_solver_global_allocate(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
		}
		else {
			// 15_10_2016 GPU CUSP bicgstab + AINV (NS Bridson)
			//cusp_solver(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);// Рабочий.
			cusp_solver_host(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
		}
#else
		/*
		if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 13) {
			//fgmres2(sl, slb, maxelm, maxbound, dV, dX0, 2000, m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl);
			integer L = 1;
			Bi_CGStab_internal5(L, sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl);

		}
		*/
		// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
		//---->//printf("Redirecting to FGMRES(lfil) + ILU2 solver.\n");
		//--->//fgmres1(sl, slb, maxelm, maxbound, dV, dX0, 2000, my_amg_manager.m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);

			if (b_on_adaptive_local_refinement_mesh) {


				if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) {

					// 16.01.2022 Обнаружено что spai0-smoothed aggregation FGMRES(600) дал сбой на компонентах скорости. Решено перенаправить на стандартный BiCGStab+ILU(lfil).

					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					//const bool breordering = false;
					//Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
						//ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);

					// Работает только для структурированной расчётной сетки.
					// Быстрее в два раза. Достаточно относительной точности 0.1.
					specialSpeed_Alice(sl, slb, dV, dX0, maxelm, maxbound, iVar, false, alpharelax);

				}
				else if (iVar == PAM) {

					if (maxelm + maxbound < 100000) {

						// Работает только для структурированной расчётной сетки.
						// Метод П.Л. Чебышева.
						chebyshev3Dnow_Alice(sl, slb, dV, dX0, maxelm, maxbound, PAM, false);
						//getchar();
					}
					else {

						// Ю.Саад и Шульц, FGMRes, 1986.
						// Руге и Штубен, 1986.
						bool worked_successfully = false;
#ifdef _OPENMP
						omp_set_num_threads(6);
#endif
						const integer iHAVorstModification_id = 2;
						if (0 && (iVar == PAM) && (inumber_iteration_SIMPLE < 10)) {
							// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
							const bool breordering = false;
							Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
								ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE,
								color, dist_max, breordering, w, lw);

						}
						else {
							amg(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m, ifrontregulationgl,
								ibackregulationgl, iHAVorstModification_id, worked_successfully, b, lb, s_loc, ls,
								w, lw, whot_is_block, bprint_message);
						}
					}
				}
				else if (iVar==TEMP) {

					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					//printf("Redirecting to FGMRES(lfil) + ILU2 solver.\n");
					fgmres1(sl, slb, maxelm, maxbound, dV, dX0, 2000, my_amg_manager.m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
				}
				else {

					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					//printf("Redirecting to FGMRES(lfil) + ILU2 solver.\n");
					fgmres1(sl, slb, maxelm, maxbound, dV, dX0, 2000, my_amg_manager.m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
				}
			}
			else {

				if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) {

					// 16.01.2022 Обнаружено что spai0-smoothed aggregation FGMRES(600) дал сбой на компонентах скорости. Решено перенаправить на стандартный BiCGStab+ILU(lfil).

					

					//if (maxelm + maxbound < 6000000)
					{
						// Работает только для структурированной расчётной сетки.
						// Быстрее в два раза. Достаточно относительной точности 0.1.
						//if (maxelm + maxbound < 110000) 
						{

							// Вариант алгоритма на cpu. На небольших размерностях накладные расходы от передачи данных на видеокарту и обратно
							// превышат время решения на процессоре.
							specialSpeed(sl, slb, dV, dX0, maxelm, maxbound, iVar, false, alpharelax);
						}
						//else {
							// Алгоритм решения на видеокарте nvidia.
							//specialSpeed_cuda(sl, slb, dV, dX0, maxelm, maxbound, iVar, false, alpharelax);
						//}
					}
					//else {

						// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					  //  const bool breordering = false;
					    //Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
						  //       ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);
					//}

				}
				else if (iVar == PAM) {

					if (maxelm+maxbound < 110000) {
						// Работает только для структурированной расчётной сетки.
						// Метод П.Л. Чебышева.
						//chebyshev3Dnow(sl, slb, dV, dX0, maxelm, maxbound, PAM, false, false);
						bool bdivergence_chebyshev = chebyshev3Dnow_cuda(sl, slb, dV, dX0, maxelm, maxbound, PAM, false, false);
						//bicgstab_chebyshev_3Dnow_cuda(sl, slb, dV, dX0, maxelm, maxbound, PAM, false, false);
						//getchar();
						if (bdivergence_chebyshev) {
							// Переключение на алгебраический многосеточный метод в случае расходимости Чебышева.

							// Ю.Саад и Шульц, FGMRes, 1986.
							// Руге и Штубен, 1986.
							bool worked_successfully = false;

							const integer iHAVorstModification_id = 2;

							amg(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m, ifrontregulationgl,
								ibackregulationgl, iHAVorstModification_id, worked_successfully, b, lb, s_loc, ls,
								w, lw, whot_is_block, bprint_message);


							// Не всегда хорошо сходится.
							// Солвер Декниса Демидова на видеокарте.
							//const bool bprint_preconditioner_amgcl = false;
							//amgcl_solver(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar, bprint_preconditioner_amgcl, dgx, dgy, dgz, inumber_iteration_SIMPLE, w, lw, bprint_message);
						}

					}
					else {

						if (0) {
							// Ю.Саад и Шульц, FGMRes, 1986.
							// Руге и Штубен, 1986.
							bool worked_successfully = false;
#ifdef _OPENMP
							omp_set_num_threads(6);
#endif
							const integer iHAVorstModification_id = 2;
							if (0 && (iVar == PAM) && (inumber_iteration_SIMPLE < 10)) {
								// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
								const bool breordering = false;
								Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
									ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE,
									color, dist_max, breordering, w, lw);

							}
							else {
								amg(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m, ifrontregulationgl,
									ibackregulationgl, iHAVorstModification_id, worked_successfully, b, lb, s_loc, ls,
									w, lw, whot_is_block, bprint_message);
							}
						}
						else {
							

							//cg_chebyshev_3Dnow(sl, slb, dV, dX0, maxelm, maxbound, PAM, false, false);
							//bool bdivergence_chebyshev = bicgstab_chebyshev_3Dnow(sl, slb, dV, dX0, maxelm, maxbound, PAM, false, false);
							bool bdivergence_chebyshev = bicgstab_chebyshev_3Dnow_cuda(sl, slb, dV, dX0, maxelm, maxbound, PAM, false, false, lb);

							if (bdivergence_chebyshev) {
								// Переключение на алгебраический многосеточный метод в случае расходимости Чебышева.

								// Ю.Саад и Шульц, FGMRes, 1986.
							    // Руге и Штубен, 1986.
							 	bool worked_successfully = false;

								const integer iHAVorstModification_id = 2;

								amg(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m, ifrontregulationgl,
									ibackregulationgl, iHAVorstModification_id, worked_successfully, b, lb, s_loc, ls,
									w, lw, whot_is_block, bprint_message);


								// Не всегда хорошо сходится.
								// Солвер Декниса Демидова на видеокарте.
								//const bool bprint_preconditioner_amgcl = false;
								//amgcl_solver(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar, bprint_preconditioner_amgcl, dgx, dgy, dgz, inumber_iteration_SIMPLE, w, lw, bprint_message);
							}
						}
					}
				}
				else if (iVar == TEMP) {

					if (bSIMPLErun_now_for_temperature == false) {

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

						if (bmyconvective == false) {
							// Теплопередача внутри твёрдого тела.
							// Работает только для структурированной расчётной сетки.
					        // Метод П.Л. Чебышева.
							//std::cout << "alpharelax=" << alpharelax << std::endl;

							//chebyshev3Dnow(sl, slb, dV, dX0, maxelm, maxbound, TEMP, false, true);
							//cg_chebyshev_3Dnow(sl, slb, dV, dX0, maxelm, maxbound, TEMP, false, true);
							 bicgstab_chebyshev_3Dnow(sl, slb, dV, dX0, maxelm, maxbound, TEMP, false, true);
							// bicgstab_chebyshev_3Dnow_cache(sl, slb, dV, dX0, maxelm, maxbound, TEMP, false, true);
							//getchar();
						}
						else {
							// Теплопередача с конвекцией.

							// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					        //printf("Redirecting to FGMRES(lfil) + ILU2 solver.\n");
							fgmres1(sl, slb, maxelm, maxbound, dV, dX0, 2000, my_amg_manager.m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
						}
					}
					else {
						// Теплопередача в составе гидродинамики CFD.

						// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
						//printf("Redirecting to FGMRES(lfil) + ILU2 solver.\n");
						fgmres1(sl, slb, maxelm, maxbound, dV, dX0, 2000, my_amg_manager.m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
					}
				}
				else {
					// Турбулентные характеристики для моделей Спаларта Аллмареса или Ментора (Ментора-Лантгрии).

					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					printf("Redirecting to FGMRES(lfil) + ILU2 solver.\n");
					fgmres1(sl, slb, maxelm, maxbound, dV, dX0, 2000, my_amg_manager.m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
				}
			}

#endif

		}
		
		else if ((5 == iswitchsolveramg_vs_BiCGstab_plus_ILU2)
		|| (9 == iswitchsolveramg_vs_BiCGstab_plus_ILU2)
			|| (10 == iswitchsolveramg_vs_BiCGstab_plus_ILU2)) {

#if GPU_LIB_INCLUDE_MY_PROJECT == 1
			// Этот метод заимствован из библиотеки ViennaCL 1.7.1 распространяемой по
			// OpenSource MIT (X11) license.
			// В данном случае вызывается связка BiCGStab Хенка Ван Дер Ворста и алгебраический
			// мноосеточный метод (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 5).

			// Или (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 9) Связка bicgStab+ilu0.
			// На момент 2017 года amg алгоритмы в ViennaCL носят пробный или экспериментальный характер,
			// они широко предлагают связку Крыловского метода с неполным ilu0 разложением.

			// Или (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 10) Связка bicgStab+ilut (ilu treshold).
			// На момент 2017 года amg алгоритмы в ViennaCL носят пробный или экспериментальный характер,
			// они широко предлагают связку Крыловского метода с неполным ilut (ilu treshold) разложением.

			// переключение методов осуществляется напрямую в коде viennacl_solver с использованием
			// переменной iswitchsolveramg_vs_BiCGstab_plus_ILU2.
			if ((iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 5)||(iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 9)) {
				if (doubleintprecision == 1) {
					printf("ERROR ViennaCL Library!!! type int64_t is usage.\n");
					printf("Library ViennaCL 1.7.1 not supported type int64_t for long long int.\n");
					system("PAUSE");
                }
				else {
#ifdef VIENNA_CL_INCLUDE_IN_MY_PRJ
					if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 5) {
						//amg из ViennaCL так и не заработал. Он строит иерархию уровней сетки
						// за 44с на задаче в 1.1млн неизвестных. Потом долго что то итерирует и выдает 
						// переполнение inf в векторе результата. 
						viennacl_solver(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
					}
					if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 9) {
						//Для того чтобы заработало Vienna CL bicgstab+ilu0 решатель необходимо использовать тип int
						// а не int64t.
						viennacl_solver(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
						// serial - однопоточная версия.
						//viennacl_solver_serial(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
					}	
#endif
				}
			}
			else {
				//system("pause");
				//my_amg_manager.memory_size_Stress 300
				//integer m_restart = my_amg_manager.memory_size_Stress;
				integer m_restart = my_amg_manager.m_restart; // fgmres(20); // recomended
				//maxit
				//fgmres(sl, slb, maxelm, maxbound, dV, dX0, 2000, m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl);
				if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 10) {
				  fgmres1(sl, slb, maxelm, maxbound, dV, dX0, 2000, m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl,b,lb,s_loc,ls);
				}
				
				//gmres_internal1(sl, slb, maxelm, maxbound, dV, dX0, 2000, alpharelax, bprintmessage, TEMP, m_restart, ifrontregulationgl, ibackregulationgl);
				//integer L = 2;
				//maxit = 2000;
				//Bi_CGStab_internal5(L, sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl);
				printf("GMRES done.\n");
				//system("pause");
			}
#else
			// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
			//printf("Redirecting to FGMRES(20) + ILU2 solver.\n");
			//fgmres1(sl, slb, maxelm, maxbound, dV, dX0, 2000, my_amg_manager.m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);
			//Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl,w,lw);
			
			// Метод из библиотеки AMGCL работает и для гидродинамики и для температуры.
			// Дата присоединения к проекту 7.05.2019, 8.05.2019.
			// На задачу в 0.4млн неизвестных алгоритмом bicgstab+amgcl делал 
			// 63; 55; 35; 25 итераций на каждую матрицу и посчитал за 44с 690ms.
			// Для сравнения алгоритм bicgstab+amg1r5 считает эту задачу за 43-45s.
			// amg1r5+bicgstab делает 11; 12; 12; 11 итераций на каждую матрицу и посчитал за 43с 810ms.
			// Методы amg1r5 и samg amgcl дают примерно одинаковое время решения на размерности 0.4млн неизвестных.
			// Время bicgstab +amg1r5  на задаче в 1.5лн неизвестных равно 3m 9s 890ms.
			// Время bicgstab +samg amgcl  на задаче в 1.5лн неизвестных равно 3m 4s 870ms.
			// Методы amg1r5 и samg amgcl дают примерно одинаковое время решения на размерности 1.5млн неизвестных.
			
			if ((iVar== GAMMA_LANGTRY_MENTER)||(iVar== RE_THETA_LANGTRY_MENTER)) {

				// RE_THETA_LANGTRY_MENTER можно решать amgcl но будет медленней.

				// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
				const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
				Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);
			}
			else if (0&&((iVar == TURBULENT_KINETIK_ENERGY) || (iVar== TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA))) {
				// Обнаружена проблема работоспособности amgcl на уравнении для K в K-Omega модели турбулентности. 11,10,2019
				// Переключаемся на алгоритм Юсефа Саада.
				// в amgcl отключил проверку в bicgstab для rhs.

				// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
				const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
				Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering,w,lw);
			}
			else if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) {

				// 16.01.2022 Обнаружено что spai0-smoothed aggregation FGMRES(600) дал сбой на компонентах скорости. Решено перенаправить на стандартный BiCGStab+ILU(lfil).

				// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
				const bool breordering = false;
				Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m,
					ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);
			}
			else {

#ifdef AMGCL_INCLUDE_IN_MY_PROJECT 
				//if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 9) {
					//if (doubleintprecision == 1) {
						//printf("ERROR ViennaCL Library!!! type int64_t is usage.\n");
						//printf("Library ViennaCL 1.7.1 not supported type int64_t for long long int.\n");
						//system("PAUSE");
					//}
				//	else
				    //{
						//Для того чтобы заработало Vienna CL bicgstab+ilu0 решатель необходимо использовать тип int
						// а не int64t.
						//--->//viennacl_solver(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
						// serial - однопоточная версия.
						//viennacl_solver_serial(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
					//}
					
				//}
				//else if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 5) {
					//{
						// Vienna AMG
						//amg из ViennaCL так и не заработал. Он строит иерархию уровней сетки
						// за 44с на задаче в 1.1млн неизвестных. Потом долго что то итерирует и выдает 
						// переполнение inf в векторе результата. 
						//viennacl_solver(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
					//}
				//}
				//else 
				{
					if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) {
						//printf("*********Denis Demidov AMGCL...***********\n");
					}
					if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 10) {
						const bool bprint_preconditioner_amgcl = false;
						amgcl_solver(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar, bprint_preconditioner_amgcl, dgx, dgy, dgz, inumber_iteration_SIMPLE, w, lw, bprint_message);
					}
					else {
						const bool bprint_preconditioner_amgcl = true;
						amgcl_solver(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar, bprint_preconditioner_amgcl, dgx, dgy, dgz, inumber_iteration_SIMPLE, w, lw, bprint_message);
					}
				}
#else
				fgmres1(sl, slb, maxelm, maxbound, dV, dX0, 2000, my_amg_manager.m_restart, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls);

#endif
			}

#endif
		}
		else if (6 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) {
#if GPU_LIB_INCLUDE_MY_PROJECT == 1
			// Этот метод заимствован из библиотеки CUSP 0.5.1 распространяемой по
			// OpenSource Apache license 2.0.
			// В данном случае на одном ядре центрального процессора используется алгоритм
			// BiCGStab Хенка Ван дер Ворста и алгебраический многосеточный метод сглаженного аггрегирования 
			// в качестве предобуславливателя SAMG.

			cusp_solver_amghost(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
#else
		    printf("WARNING: CUSP 0.5.1 library is not connected\n");
		    printf("Redirecting to BiCGStab + ILU2 solver.\n");
			// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
			const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
			Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering,w,lw);
#endif
		}
		else if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 8) {
#if GPU_LIB_INCLUDE_MY_PROJECT == 1
			// Этот метод заимствован из библиотеки CUSP 0.5.1 распространяемой по
			// OpenSource Apache license 2.0.
			// В данном случае на всех ядрах графического ускорителя в FP64 точности используется алгоритм
			// BiCGStab Хенка Ван дер Ворста и алгебраический многосеточный метод сглаженного аггрегирования 
			// в качестве предобуславливателя SAMG.

			// Geforce GTX 1080 Ti имеет 0.388ТФЛОПС в FP64 точности.
			// Для сравнения один поток процессора core i7 6850K имеет всего 0.0144ТФЛОПС.

			cusp_solver_GPU_SAMG(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
#else
		    printf("WARNING: CUSP 0.5.1 library is not connected\n");
		    printf("Redirecting to BiCGStab + ILU2 solver.\n");
			// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
			const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
			Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering,w,lw);
#endif
		}
		else if ((iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 3) 
		|| (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 7)) {

		integer iswitchsolveramg_vs_BiCGstab_plus_ILU2_memo_loc = iswitchsolveramg_vs_BiCGstab_plus_ILU2;
		iswitchsolveramg_vs_BiCGstab_plus_ILU2 = 7;


			if (iswitchsolveramg_vs_BiCGstab_plus_ILU2_memo_loc==3) {

				if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == NUSHA) ||
					(iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
					(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)||
					(iVar == GAMMA_LANGTRY_MENTER) || (iVar == RE_THETA_LANGTRY_MENTER)) {
					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
					Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering,w,lw);
				}
				else {
					// Для температуры и поправки давления.
					// Температура вопрос.

					
					//РУМБАv0_14
					real_mix_precision ret74 = 0.0;
					my_agr_amg_loc_memory(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m,  ret74,
						b, lb, ifrontregulationgl, ibackregulationgl, s_loc, ls, inumber_iteration_SIMPLE,
						w,lw, whot_is_block);


				}
			}
			else {
				// Только РУМБАv0_14
				// 11 января 2016. классический агломеративный алгебраический многосеточный метод.
			    // Это моя собственная разработка РУМБА 0.14.
			    //if (iVar != PAM) {
			    //doublereal theta82 = 0.24;
			    //doublereal theta83 = 0.23;
			    //doublereal magic82 = 0.4;
			    //doublereal magic83 = 0.5; 
			    //doublereal ret74 = 0.0;
			    //-->doublereal theta82 = 0.24;// 0.25; //0.24
			    //-->doublereal theta83 = 0.23;// 0.25; // 0.23
			    // 0.3 0.4 0.44 0.45  0.5
			    // 16  25   38        17
			    //--->doublereal magic82 = 0.4; // 0.35; // 0.4 // 0.43
			    //----->doublereal magic83 = 0.4;// 0.35; // 0.42

				if ((iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA)) {
					// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
					const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
					Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w, lw);
				}
				else {

					//real_mix_precision theta82 = static_cast<real_mix_precision>(my_amg_manager.theta);
					//real_mix_precision theta83 = static_cast<real_mix_precision>(my_amg_manager.theta);
					//real_mix_precision magic82 = static_cast<real_mix_precision>(my_amg_manager.magic);
					//real_mix_precision magic83 = static_cast<real_mix_precision>(my_amg_manager.magic);

					real_mix_precision ret74 = 0.0;
					my_agr_amg_loc_memory(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m, ret74, 
						b, lb, ifrontregulationgl, ibackregulationgl, s_loc, ls, inumber_iteration_SIMPLE,
						w, lw, whot_is_block);
				
				}
				/*
				}
				else {
				//doublereal theta82=0.24;
				//doublereal theta83 = 0.23;
				//doublereal magic82 = 0.4;
				//doublereal magic83 = 0.5;
				//doublereal ret74 = 0.0;
				errno_t err_optimetric;
				FILE* fp_optimetric;
				err_optimetric = fopen_s(&fp_optimetric, "optimetric.txt", "a");
				if (err_optimetric != 0) {
				printf("Error open file log.txt\n");
				printf("Please, press any key to continue...\n");
				//system("pause");
				system("pause");
				exit(0);
				}
				for (doublereal theta82 = 0.21; theta82 < 0.26; theta82 += 0.01) {
				for (doublereal theta83 = 0.21; theta83 < 0.26; theta83 += 0.01) {
				for (doublereal magic82 = 0.3; magic82 < 0.35; magic82 += 0.01) {
				for (doublereal magic83 = 0.35; magic83 < 0.36; magic83 += 0.01) {
				doublereal ret74 = 0.0;
				for (integer i26 = 0; i26 < maxelm + maxbound; i26++) {
				dX0[i26] = 0.0;// init
				}
				my_agr_amg_loc_memory(sl, slb, maxelm, maxbound, dV, dX0, alpharelax, iVar, bLRfree, m, theta82, theta83, magic82, magic83, ret74);
				fprintf(fp_optimetric, "theta82=%e theta83=%e magic82=%e magic83=%e ret74=%e\n", theta82, theta83, magic82, magic83, ret74);
				}
				}
				}
				}
				fclose(fp_optimetric);
				printf("optimisation compleate\n");
				system("pause");
				}
				*/
			}
			iswitchsolveramg_vs_BiCGstab_plus_ILU2 = iswitchsolveramg_vs_BiCGstab_plus_ILU2_memo_loc;
			

		}
		else {		

			if (4 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) {
				// cusp call.
#if GPU_LIB_INCLUDE_MY_PROJECT == 1
				//cusp_solver_GPU_AINV_Bridson(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
				// smoothes aggregation algebraic multigrid
				cusp_solver_GPU_SAMG(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, iVar);
#else
				printf("Redirecting to BiCGStab + ILU2 solver.\n");
				// старый добрый проверенный метод Ю. Саада из SPARSKIT2.
				const bool breordering = false; // только false никакого reordering переупорядочивания 17,10,2020.
				Bi_CGStab_internal3(sl, slb, maxelm, maxbound, dV, dX0, maxit, alpharelax, bprintmessage, iVar, m, ifrontregulationgl, ibackregulationgl, b, lb, s_loc, ls, inumber_iteration_SIMPLE, color, dist_max, breordering, w,lw);
#endif
			}

		}
	}

	if (iVar == TURBULENT_KINETIK_ENERGY) {

		const integer isize = maxelm + maxbound;

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < isize; ++i_1) {
			dV[i_1] *= 1.0e-6;
			dX0[i_1] *= 1.0e-6;
		}
	}

	// Замечательный по своей простоте прямой метод.
	//Direct(sl, slb, maxelm, maxbound, dV, dX0);
	//hypreSolve(sl, slb, maxelm, maxbound, dV, dX0);

	/*
	switch (iVar) {
	case VX: SOR3Dnow(sl, slb, dX0,  maxelm, maxbound, VX);
	break;
	case VY: SOR3Dnow(sl, slb, dX0,  maxelm, maxbound, VY);
	break;
	case VZ: SOR3Dnow(sl, slb, dX0,  maxelm, maxbound, VZ);
	break;
	case PAM: SOR3Dnow(sl, slb, dX0, maxelm, maxbound, PAM);
	break;
	}
	*/
	calculation_main_end_time = clock();
	calculation_vorst_seach_time += calculation_main_end_time - calculation_main_start_time;
}

#endif // !MY_LINALG_CPP
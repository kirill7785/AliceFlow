// my_aggregat_amg.cpp: определяет точку входа для консольного приложения.
// classical aglomeration algebraic multigrid method.
// классический агломеративный алгебраический многосеточный метод.
// 31 - октября 2015. Версия 0.06 в которой выполнена работа по дальнейшему улучшению
// быстродействия алгоритма. Исправлен ряд ошибок. Задача Пуассона на сетке 800х800 решается за 9мин 32с.
// Данная версия алгоритма содержит очень большое количество оптимизирующих патчей и 
// поэтому код стало очень сложно сопровождать. Надо провести реальное тестирование в следующей версии кода.
// 22-23 октября 2015. Версия 0.05 в которой запрограммирована идея Писсанецки для
// разреженного матричного умножения. Реализация идеи Писсанецки ускорила код на 15%
// по сравнению с версией 0.04. Также в версии 0.05 реализован алгоритм Густавсона для
// разреженного матричного умножения gustavson sparse matrix multiplayer IBM 1978.
// Алгоритм Густавсона ускорил код по сравнению с версией 0.04 более чем вдвое.
// 18 октября 2015. Полностью работоспособный мультигрид.
// Тестировалось на условиях Дирихле но должно работать на любых 
// краевых задачах. 18 октября датируется версия 0.04. Версия 0.04 на треть
// быстрее версии 0_03. Были ускорены как операции построения C-F разбиения, 
// так и нахождение оператора Галёркина. При нахождении С-F разбиения 
// учитывается уже построеннная его часть и поэтому число сканирований на
// на поздних циклах сокращается охватывая только не построенную часть.
// При нахождении произведеия Галёркина получена самая оптимальная по быстродествию версия,
// Основанная на алгоритме слияния отсортированных списков.
// version 0.03 15 october 2015. (Успешное построение иерархии сеток, операторов restriction , interpolation.)
// По сравнению с версией 0.01 исправлены : бинарный поиск, нахождение оператора Галеркина. Внесены прочие мелкие исправления.
// По сравнению с версией 0.02 код стал быстрее, продуманнее. Найдена ошибка в MergeSort поэтому от неё полностью отказались в пользу HeapSort.
// Код работает на любых размерностях до 10млн ненулевых элементов но очень медленно на высоких размерностях.
// Код быстрее версии 0.02 на 50%.
/*    ###     ##    ###   ####
*    #  #    #  #  # ##  #   
*   #####   #    #   ##  #   ##  for AliceFlowv0_22 AliceMeshv_0_34 and anather products...
*  #    #  #         ##   ####
*/



#include "stdafx.h"
#include <math.h>
#include <stdlib.h>
#include <windows.h>
//#include <omp.h>

#define Real double
#define integer int

typedef struct TAk {
	integer i, j;
	Real aij;
	//в последней версии алгоритма ind не используется.
	integer ind; // позиция в первоначальной сортировке.
} Ak;

typedef struct TList {
	integer ii, i, countsosed;
	TList* next;
	TList* prev;
} List;

// Для генерации матрицы СЛАУ требуется в случае реализации
// на динамических массивах переупорядочивание элементов:
// сортировка. Здесь будет реализована быстрая сортировка.
// Брайан Керниган и Денис Ритчи "The C programming language".
// swap: Обмен местами A[i] и A[j]
void swap(Ak * &A, integer i, integer j)
{
	Ak A_temp;

	// change A[i] <-> A[j]
	A_temp = A[i];
	A[i] = A[j];
	A[j] = A_temp;

} // swap

/*
  // Вот алгоритм PivotList
integer PivotList(Ak * &A, integer n, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue = A[first].i*n+A[first].j;
	integer PivotPoint = first;

	for (integer index = (first + 1); index <= last; index++) {
		if (A[index].i*n+A[index].j<PivotValue) {
			PivotPoint++;
			swap(A, PivotPoint, index);
		}
	}

	swap(A, first, PivotPoint);

	return PivotPoint;
} // PivotListamg
*/
// Вот алгоритм PivotList
integer PivotList(Ak * &A, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue_j = A[first].j;
	integer PivotValue_i = A[first].i;
	integer PivotPoint = first;

	for (integer index = (first + 1); index <= last; index++) {
		if (A[index].i<PivotValue_i) {
			PivotPoint++;
			swap(A, PivotPoint, index);
		}
		else if ((A[index].i == PivotValue_i) && (A[index].j < PivotValue_j)) {
			PivotPoint++;
			swap(A, PivotPoint, index);
		}
	}

	swap(A, first, PivotPoint);

	return PivotPoint;
} // PivotList

  // Вот алгоритм PivotList
integer PivotList_j(Ak * &A, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue_j = A[first].j;
	integer PivotValue_i = A[first].i;
	integer PivotPoint = first;

	for (integer index = (first + 1); index <= last; index++) {
		if ( A[index].j<PivotValue_j) {
			PivotPoint++;
			swap(A, PivotPoint, index);
		}
		else if ((A[index].j == PivotValue_j) && (A[index].i < PivotValue_i)) {
			PivotPoint++;
			swap(A, PivotPoint, index);
		}
	}

	swap(A, first, PivotPoint);

	return PivotPoint;
} // PivotListamg_j


  // Пирамидальная сортировка

  // Переформировать пирамиду
void FixHeap_j(Ak* &A,
	integer root,
	Ak m,
	integer bound,
	/*integer n,*/ integer iadd)
{
	integer vacant;
	integer largerChild;

	// list сортируемый список пирамида
	// root номер корня пирамиды
	// m ключевое значение вставляемое в пирамиду
	// bound правая граница (номер) в пирамиде
	vacant = root;
	while (2 * vacant <= bound)
	{
		largerChild = 2 * vacant;
		integer lCadd = largerChild + iadd;
		integer lCadd1 = lCadd + 1;

		// поиск наибольшего из двух непосредственных потомков
		bool compare_result = false;
		if (A[lCadd1].j > A[lCadd].j) {
			compare_result = true;
		}
		else if (A[lCadd1].j == A[lCadd].j) {
			if (A[lCadd1].i > A[lCadd].i) {
				compare_result = true;
			}
		}
		if ((largerChild<bound) && ( compare_result /*A[largerChild + 1+iadd].j*n+ A[largerChild + 1 + iadd].i> A[largerChild + iadd].j*n+ A[largerChild + iadd].i*/))
		{
			largerChild = largerChild + 1;
		}

		lCadd = largerChild + iadd;
		// находится ли ключ выше текущего потомка ?
		compare_result = false;
		if (m.j > A[lCadd].j) {
			compare_result = true;
		}
		else if (m.j == A[lCadd].j) {
			if (m.i > A[lCadd].i) {
				compare_result = true;
			}
		}
		if ( compare_result /*m.j*n+m.i >  A[largerChild + iadd].j*n+ A[largerChild + iadd].i*/)
		{
			// да, цикл завершается
			break;
		}
		else
		{
			// нет, большего непосредственного потомка
			// следует поднять
			A[vacant + iadd] = A[lCadd];
			vacant = largerChild;
		}
	}
	A[vacant + iadd] = m;
} // FixHeap_j



  // Переформировать пирамиду
void FixHeap(Ak* &A,
	integer root,
	Ak m,
	integer bound,
	integer n, integer iadd)
{
	integer vacant;
	integer largerChild;

	// list сортируемый список пирамида
	// root номер корня пирамиды
	// m ключевое значение вставляемое в пирамиду
	// bound правая граница (номер) в пирамиде
	vacant = root;
	while (2 * vacant <= bound)
	{
		largerChild = 2 * vacant;
		integer lCadd = largerChild + iadd;
		integer lCadd1 = lCadd + 1;

		// поиск наибольшего из двух непосредственных потомков
		//integer key1 = A[largerChild + 1 + iadd].i*n + A[largerChild + 1 + iadd].j;
		//integer key2 = A[largerChild + iadd].i*n + A[largerChild + iadd].j;
		bool compare_result = false;
		if (A[lCadd1].i > A[lCadd].i) {
			compare_result = true;
		}
		else if (A[lCadd1].i == A[lCadd].i) {
			if (A[lCadd1].j > A[lCadd].j) {
				compare_result = true;
			}
		}
		if ((largerChild<bound) && compare_result/*(key1>key2)*/)
		{
			largerChild = largerChild + 1;
		}

		lCadd = largerChild + iadd;
		// находится ли ключ выше текущего потомка ?
		//integer key5 = m.i*n + m.j;
		//integer key6 = A[largerChild + iadd].i*n + A[largerChild + iadd].j;
		compare_result = false;
		if (m.i > A[lCadd].i) {
			compare_result = true;
		}
		else if (m.i == A[lCadd].i) {
			if (m.j > A[lCadd].j) {
				compare_result = true;
			}
		}
		if (compare_result/*key5 > key6*/)
		{
			// да, цикл завершается
			break;
		}
		else
		{
			// нет, большего непосредственного потомка
			// следует поднять
			A[vacant + iadd] = A[lCadd];
			vacant = largerChild;
		}
	}
	A[vacant+iadd] = m;
} // FixHeap

// HeapSort очень быстрая только на размерностях до 10^5 а после вплоть до 10^8 она может замедляться раза в три. 

  // Пирамидальная сортировка оптимальна как
  // по памяти, так и по быстродействию, к тому же её алгоритм
  // очень интересен.
  // Ограничение состоит в том, что нумерация массива должна начинаться с 1.
void HeapSort(Ak * &A, integer n, integer first, integer last)
{

	Ak maxelm; // элемент с наибольшим значением ключа

								 // конструирование пирамиды
	for (integer i =  ((last-first+1) / 2); i >= 1; i--)
	{
		FixHeap(A, i, A[i + first - 1], last-first+1,n,first-1);
	}
	for (integer i = last-first+1; i >= 2; i--)
	{
		// скопировать корень пирамиды в список
		// переформировать пирамиду
		maxelm = A[first];
		FixHeap(A, 1, A[i+first-1], i - 1,n,first-1);
		A[i+first-1] = maxelm;
	}
} // HeapSort



  // Пирамидальная сортировка оптимальна как
  // по памяти, так и по быстродействию, к тому же её алгоритм
  // очень интересен.
  // Ограничение состоит в том, что нумерация массива должна начинаться с 1.
void HeapSort_j(Ak * &A,  integer first, integer last)
{

	Ak maxelm; // элемент с наибольшим значением ключа

	integer iadd = first - 1;
			   // конструирование пирамиды
	for (integer i =  ((last - first + 1 )/ 2); i >= 1; i--)
	{
		FixHeap_j(A, i, A[i+iadd], last - first + 1,  iadd);
	}
	for (integer i = last-first+1; i >= 2; i--)
	{
		// скопировать корень пирамиды в список
		// переформировать пирамиду
		maxelm = A[1+iadd];
		FixHeap_j(A, 1, A[i+iadd], i - 1,  iadd);
		A[i + iadd] = maxelm;
	}
} // HeapSort_j

// Сортировка слиянием. 
// Требует дополнительной пямяти.
void MergeSort(Ak * &Aorig, integer size) {
	// предполагается индексация от нуля до size-1.
	// Массив А предполагается не менее двойного размера.
	// Двойная память это недостаток данного алгоритма.
	Ak* A = Aorig;
	Ak* B = A + size;
	Ak* C;

	// A - сортируемый массив, B - вспомогательная память (справа от основных данных в А такого-же размера что и основные данные в А).
	// C - указатель для обмена.
	for (int i = 1; i < size; i = i * 2) // размер объединяемых фрагментов
	{
		for (int j = 0; j < size; j = j + 2 * i) // начало первого из объединяемых 
			// фрагментов
		{
			int r = j + i; // начало второго из объединяемых фрагментов
			int n1=0, n2=0;
			if (i < size - j) { n1 = i; }
			else { n1 = size - j; };
			if (i < size - r) { n2 = i; }
			else { n2 = size - r; };

			if (n1 < 0) n1 = 0;
			if (n2 < 0) n2 = 0;

			// слияние упорядоченных фрагментов
			for (int ia = 0, ib = 0, k = 0; k < n1 + n2; k++)
			{
				if (ia >= n1) B[j + k] = A[r + ib++];
				else
					if (ib >= n2) B[j + k] = A[j + ia++];
					else {
						bool compare_result = false;
						int lCadd = j + ia;
						int lCadd1 = r + ib;
						if (A[lCadd1].i > A[lCadd].i) {
							compare_result = true;
						}
						else if (A[lCadd1].i == A[lCadd].i) {
							if (A[lCadd1].j > A[lCadd].j) {
								compare_result = true;
							}
						}
						if (compare_result) {
							B[j + k] = A[j + ia++];
						}
						else {
							B[j + k] = A[r+ib++];
						}
					}
			}



		}
		C = A; A = B; B = C;

	}

	C = A; A = B; B = C;

	// Копирование, если результат размещен не в основном а в вспомогательном массиве
	if (B != Aorig)
		memcpy(Aorig, B, size*sizeof(Ak));

	A = NULL; B = NULL; C = NULL;
} // MergeSort


// Сортировка слиянием. 
// Требует дополнительной пямяти.
void MergeSort_j(Ak * &Aorig, integer size) {
	// предполагается индексация от нуля до size-1.
	// Массив А предполагается не менее двойного размера.
	// Двойная память это недостаток данного алгоритма.
	Ak* A = Aorig;
	Ak* B = A + size;
	Ak* C;

	// A - сортируемый массив, B - вспомогательная память (справа от основных данных в А такого-же размера что и основные данные в А).
	// C - указатель для обмена.
	for (int i = 1; i < size; i = i * 2) // размер объединяемых фрагментов
	{
		for (int j = 0; j < size; j = j + 2 * i) // начало первого из объединяемых 
			// фрагментов
		{
			int r = j + i; // начало второго из объединяемых фрагментов
			int n1 = 0, n2 = 0;
			if (i < size - j) { n1 = i; }
			else { n1 = size - j; };
			if (i < size - r) { n2 = i; }
			else { n2 = size - r; };

			if (n1 < 0) n1 = 0;
			if (n2 < 0) n2 = 0;

			// слияние упорядоченных фрагментов
			for (int ia = 0, ib = 0, k = 0; k < n1 + n2; k++)
			{
				if (ia >= n1) B[j + k] = A[r + ib++];
				else
					if (ib >= n2) B[j + k] = A[j + ia++];
					else {
						bool compare_result = false;
						int lCadd = j + ia;
						int lCadd1 = r + ib;
						if (A[lCadd1].j > A[lCadd].j) {
							compare_result = true;
						}
						else if (A[lCadd1].j == A[lCadd].j) {
							if (A[lCadd1].i > A[lCadd].i) {
								compare_result = true;
							}
						}
						if (compare_result) {
							B[j + k] = A[j + ia++];
						}
						else {
							B[j + k] = A[r + ib++];
						}
					}
			}



		}
		C = A; A = B; B = C;

	}

	C = A; A = B; B = C;

	// Копирование, если результат размещен не в основном а в вспомогательном массиве
	if (B != Aorig)
		memcpy(Aorig, B, size*sizeof(Ak));

	A = NULL; B = NULL; C = NULL;
} // MergeSort_j



  // Быстрая сортировка Хоара.
  // Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
  // стр. 106.
void QuickSort(Ak * &A,  integer first, integer last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	if (0) {
		/*
		// BubbleSort
		integer numberOfPairs = last - first + 1;
		bool swappedElements = true;
		while (swappedElements) {
			numberOfPairs--;
			swappedElements = false;
			for (integer i = first; i <= first + numberOfPairs - 1; i++) {
				if (A[i].i*n+A[i].j>A[i + 1].i*n+A[i+1].j) {
					swap(A, i, i + 1);
					swappedElements = true;
				}
			}
		}
		*/
	}
	else
	{
		integer pivot;

		if (first < last) {
			pivot = PivotList(A, first, last);
			QuickSort(A, first, pivot - 1);
			QuickSort(A, pivot + 1, last);
		}
	}
} // QuickSort

  // Быстрая сортировка Хоара.
  // Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
  // стр. 106.
void QuickSort_j(Ak * &A, integer first, integer last) {
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
			for (integer i = first; i <= first + numberOfPairs - 1; i++) {
				if (A[i].j > A[i + 1].j) {
					swap(A, i, i + 1);
					swappedElements = true;
				}
			}
		}
	}
	else
	{
		integer pivot;

		if (first < last) {
			pivot = PivotList_j(A,  first, last);
			QuickSort_j(A, first, pivot - 1);
			QuickSort_j(A, pivot + 1, last);
		}
	}
} // QuickSort_j

/*
// Это исторически первоначальный код содержащий лишь построение последовательности 
// вложенных матриц и не иодержащий построения операторов restriction и prolongation. 
// создание этого кода завершено 1 сентября 2015 года.

integer aggregative_amg(Ak* &A, integer nnz, integer n, Real* &x, Real* &b) {


// нумерация начинается с единицы.

const integer maxlevel = 30;
integer ilevel = 1;
integer nnz_a[maxlevel];
integer n_a[maxlevel];
nnz_a[0] = nnz;
n_a[0] = n;
bool* flag = new bool[n + 1];
integer iadd = 0;

while ((ilevel<maxlevel-1)&&(n_a[ilevel - 1] > 50)) {
//heapsort(A, key=i*n_a[ilevel - 1] + j, 1, nnz_a[ilevel - 1]);
QuickSort(A, n_a[ilevel - 1], 1 + iadd, nnz_a[ilevel - 1] + iadd);

for (integer ii = 1+iadd; ii <= nnz_a[ilevel - 1]+iadd; ii++) {
//	printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n",ii,A[ii].aij,ii,A[ii].i,ii,A[ii].j);
//if (ii % 20 == 0) getchar();
}

for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
flag[ii] = false;
}
for (integer ii = 1+iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
A[ii].ind = ii-iadd;
}



// Copy(A) на nnz ячеек правее.
for (integer ii = nnz_a[ilevel - 1] + 1 + iadd; ii <= 2 * nnz_a[ilevel - 1] + iadd; ii++) {
A[ii] = A[ii - nnz_a[ilevel - 1]];
}

integer n_coarce = 1; // номер агрегата.
const integer max_sosed = 20;
const integer NULL_SOSED = -1;
integer vacant = NULL_SOSED;
for (integer ii = 1 + iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
if (flag[A[ii].i] == false) {
// Вычисляем по немодифиуцированной матрице А (хранящейся слева).

integer set[max_sosed]; // не более 20 узлов в одном агрегате.
for (integer js = 0; js < max_sosed; js++) {
set[js] = NULL_SOSED;
}
integer ic = 0;
set[ic] = A[ii].i;
ic++;

// если узел j ещё не был добавлен в агрегат.
if (flag[A[ii].j] == false) {
vacant = A[ii].j;
for (integer js = 0; js < ic; js++) {
if (vacant == set[js]) {
vacant = NULL_SOSED;
}
}
if (vacant != NULL_SOSED) {
set[ic] = vacant;
ic++;
}
}
integer iscan = ii + 1;
while ((iscan <= nnz_a[ilevel - 1] + iadd) && (A[iscan].i == set[0])) {
// если узел j ещё не был добавлен в агрегат.
if (flag[A[iscan].j] == false) {
vacant = A[iscan].j;
for (integer js = 0; js < ic; js++) {
if (vacant == set[js]) {
vacant = NULL_SOSED;
}
}
if (vacant != NULL_SOSED) {
set[ic] = vacant;
ic++;

}
}

iscan++;

} // while


// (i,j) -> (I,J)
// модифицируем копию A находящуюся справа.
for (integer k = nnz_a[ilevel - 1] + 1 + iadd; k <= 2 * nnz_a[ilevel - 1] + iadd; k++) {
bool found = false;
for (integer k1 = 0; k1 < ic; k1++) {
if (A[k - nnz_a[ilevel - 1]].i == set[k1]) found = true;
}
if (found) A[k].i = n_coarce;
found = false;
for (integer k1 = 0; k1 < ic; k1++) {
if (A[k - nnz_a[ilevel - 1]].j == set[k1]) found = true;
}
if (found) A[k].j = n_coarce;
}

// Помечаем узлы как включённые в агрегат.
for (integer js = 0; js < ic; js++) {
flag[set[js]] = true;
}

n_coarce++;

// Один агрегат создан.



} // узел не был ещё включён в агрегат.
} // агрегаты созданы.

//printf("%d %d\n",n,n_coarce-1);
//getchar();
n_a[ilevel] = n_coarce - 1;

// сортировка по новому ключу key=i*(iglcoarce_number-1)+j;
// в позиции ind сохранён индекс предыдущей позиции.
//heapsort(A, key = i*(n_coarce - 1) + j, nnz + 1, 2 * nnz);
QuickSort(A, n_coarce - 1, nnz_a[ilevel-1] + 1 + iadd, 2 * nnz_a[ilevel-1] + iadd);


for (integer ii = 1+nnz_a[ilevel-1]+iadd; ii <= 2*nnz_a[ilevel - 1]+iadd; ii++) {
//	printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n",ii,A[ii].aij,ii,A[ii].i,ii,A[ii].j);
//if (ii % 20 == 0) getchar();
}



// копирование в третью часть матрицы Слау на coarce сетке.
// инициализация.
// запас памяти 1 nnz на fine уровне.
for (integer ii = 2 * nnz_a[ilevel - 1]+1 + iadd; ii <= 3 * nnz_a[ilevel - 1] + iadd; ii++) {
A[ii].aij = 0.0;
A[ii].ind = NULL_SOSED;
A[ii].i = NULL_SOSED;
A[ii].j = NULL_SOSED;
}

for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) flag[ii] = false;


integer ic1 = 2 * nnz_a[ilevel - 1] + 1 + iadd;
integer im = 1;
integer im0 = 0;
for (integer ii = nnz_a[ilevel - 1] + 1 + iadd; ii <= 2 * nnz_a[ilevel - 1] + iadd; ii++) {
if (flag[A[ii].i] == false) {
integer istr = A[ii].i;
while ((ii + im0 <= 2 * nnz_a[ilevel - 1] + iadd)&&(istr==A[ii+im0].i)) {

if (ic1 <= 3 * nnz_a[ilevel - 1] + iadd) {
A[ic1].i = A[ii + im0].i;
A[ic1].j = A[ii + im0].j;
A[ic1].aij += A[ii + im0].aij;
while ((ii + im <= 2 * nnz_a[ilevel - 1] + iadd) && (A[ii + im0].i == A[ii + im].i) && (A[ii + im0].j == A[ii + im].j))
{
A[ic1].aij += A[ii + im].aij;
im++;
}
ic1++;
im0 = im;
im++;
}
else {
printf("error 1\n");
getchar();
}
}
flag[A[ii].i] = true;
im = 1;
im0 = 0;
}
}

nnz_a[ilevel] = ic1 - 1 - 2 * nnz_a[ilevel - 1] - iadd;
iadd += 2 * nnz_a[ilevel - 1];

printf("nnz : fine=%d, coarse=%d, operator complexity=%e. n : fine=%d, coarse=%d grid complexity=%e.\n", nnz_a[ilevel - 1], nnz_a[ilevel], (double)(nnz_a[ilevel])/(double)( nnz_a[ilevel - 1]), n_a[ilevel - 1], n_a[ilevel], (double)(n_a[ilevel])/ (double)(n_a[ilevel - 1]));
getchar();

ilevel++;// грубосеточная матрица построена.

} // иерархия сеток построена.

return 0;

} // aggregative_amg


*/

// seidel и residual вся информация  хранится в индексах i, j поэтому вектор х обрабатывается верно в плане нумерации,
// но для порядка обхода матрицы важно чтобы она была упорядочена по индексу i.


// smoother.
// 1 september 2015.
void seidel(Ak* &A, integer istart, integer iend, Real* &x, Real* &b, bool* &flag, integer n)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[A[ii].i] == false) {
			integer istr = A[ii].i;
			integer ic = ii;
			Real ap = 0.0;
			x[istr] = b[istr];
			while ((ic<=iend)&&(A[ic].i == istr)) {
				if (A[ic].j != istr) {
					x[istr] += -A[ic].aij*x[A[ic].j];
				}
				else ap = A[ic].aij;
				ic++;
			}
			/*
			if (fabs(ap) < 1.0e-30) {
				printf("zero diagonal elements in string %d",istr);
				getchar();
				exit(1);
			}
			else */{
				x[istr] /= ap;
			}
			flag[A[ii].i] = true;
		}
	}


} // seidel
/*
integer *row_ptr_start = new integer[4 * n_a[0] + 1];
integer *row_ptr_end = new integer[4 * n_a[0] + 1];
// istart - начальная позиция ненулевых элементов в матрице А.
// iend - конечная позиция ненулевых элементов в матрице А.
for (integer i = 1; i <= n; i++) {
	flag[i] = false;
}
for (integer ii = 1; ii <= nnz_a[0]; ii++) {
	if (flag[A[ii].i] == false) {
		integer istr = A[ii].i;
		integer ic = ii;
		integer icdiag = ii;
		row_ptr_start[istr] = ii;
		Real ap = 0.0;
		//x[istr] = b[istr];
		while ((ic <= nnz_a[0]) && (A[ic].i == istr)) {
			if (A[ic].j != istr) {
				//x[istr] += -A[ic].aij*x[A[ic].j];
			}
			else {
				ap = A[ic].aij;
				icdiag = ic;
			}
			ic++;
		}
		row_ptr_end[istr] = ic - 1;
		if (fabs(ap) < 1.0e-30) {
			printf("zero diagonal elements in string %d", istr);
			getchar();
			exit(1);
		}
		else {
			//x[istr] /= ap;
		}

		flag[A[ii].i] = true;
		Ak temp = A[ii];
		A[ii] = A[icdiag];
		A[icdiag] = temp;
		A[ii].aij = 1.0 / ap; // умножение быстрей деления.
	}
}
*/

// smoother.
// 9 september 2015.
// q - quick.
void seidelq(Ak* &A, integer istartq, integer iendq, Real* &x, Real* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq+iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		x[istr] = b[istr];
	
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++) 
		{
			x[istr] += -A[ii1].aij*x[A[ii1].j];
		}
		x[istr] *= A[row_ptr_start[ii]].aij;
	}


} // seidelq

/*
// smoother.
// 9 september 2015.
// q - quick.
// Якоби расходится.
void seidelq(Ak* &A, integer istartq, integer iendq, Real* &x, Real* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	//Real *sum=new Real[endpos-startpos+1];

#pragma omp parallel for
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		x[istr] = b[istr];
		//	Real sum = 0.0;
		//#pragma omp parallel for reduction(+:sum)
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++)
		{
			x[istr] += -A[ii1].aij*x[A[ii1].j];
			//sum = sum - A[ii1].aij*x[A[ii1].j];
		}
		//x[istr] += sum;
		x[istr] *= A[row_ptr_start[ii]].aij;
	}

	
	//for (integer ii = startpos; ii <= endpos; ii++) {
	//	integer istr = ii - iadd;
	//	x[istr] = sum[istr]*A[row_ptr_start[ii]].aij;
	//}
	

	//delete[] sum;
} // seidelq
*/


// residual.
// 2 september 2014
void residual(Ak* &A, integer istart, integer iend, Real* &x, Real* &b, bool* &flag, integer n, Real* &residual)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
		residual[i] = 0.0;
	}
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[A[ii].i] == false) {
			integer istr = A[ii].i;
			integer ic = ii;
			residual[istr] = b[istr];
			while ((ic <= iend) && (A[ic].i == istr)) {
				if (A[ic].j != istr)  residual[istr] += -A[ic].aij*x[A[ic].j];
				else residual[istr] += -A[ic].aij*x[istr];
				ic++;
			}
			flag[A[ii].i] = true;
		}
	}


} // residual

// smoother.
// 9 september 2015.
// q - quick.
void residualq(Ak* &A, integer istartq, integer iendq, Real* &x, Real* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, Real* &residual)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		residual[istr] = b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++) { residual[istr] += -A[ii1].aij*x[A[ii1].j]; }
		residual[istr] += (-1.0/A[row_ptr_start[ii]].aij)*x[istr];
	}


} // residualq

// smoother.
// 9 september 2015.
// q - quick.
void residualqspeshial(Ak* &A, integer istartq, integer iendq, Real* &x, Real* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, Real* &residual)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		residual[istr] = b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++) { residual[istr] += -A[ii1].aij*x[A[ii1].j]; }
		//residual[istr] -=  A[row_ptr_start[ii]].aij*x[istr];
		// 2.0 2.441 2.543 2.546
		//residual[istr] -= 2.0*A[row_ptr_start[ii]].aij*x[istr];
		Real omega = 1.855; // SOR
		//residual[istr] += ((-1.0 / A[row_ptr_start[ii]].aij)*x[istr]); // верный вариант.
		residual[istr] += ((-1.0 / A[row_ptr_start[ii]].aij)*x[istr])/16.0;
		//residual[istr] *= omega;
		// 2.423
		// 1.855 4.448
		// 1.0  5.074
		// 8.713 5.706
	}


} // residualqspeshial

// restriction
// 3 september 2015.
// 6 september 2015.
// должна быть отсортирована по i, порядок по j неважен.
void restriction(Ak* &R, integer istart, integer iend, bool* &flag, Real * &x_fine, Real * &x_coarse, integer n_fine, integer n_coarse) {
	// x_coarse[1..n_coarse]=R n_coarse x n_fine*x_fine[1..n_fine];
	for (integer i = 1; i <= n_coarse; i++) {
		flag[i] = false;
	}
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[R[ii].i] == false) {
			integer istr = R[ii].i;
			integer ic = ii;
			x_coarse[istr] = 0.0;
			while ((ic <= iend) && (R[ic].i == istr)) {
				x_coarse[istr] += R[ic].aij*x_fine[R[ic].j];
				ic++;
			}
			flag[R[ii].i] = true;
		}
	}

} // restriction

  // prolongation
  // 3 september 2015.
// 6 september 2015.
  // должна быть отсортирована по j, порядок по i неважен.
// P=transpose(R).
void prolongation(Ak* &P, integer istart, integer iend, bool* &flag, Real * &x_fine, Real * &x_coarse, integer n_fine, integer n_coarse) {
	// x_fine[1..n_fine]=P n_fine x n_coarse*x_coarse[1..n_coarse];
	// P=c*transpose(R); Pij=c*Rji;
	for (integer j = 1; j <= n_fine; j++) {
		flag[j] = false;
	}
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[P[ii].j] == false) {
			integer jstr = P[ii].j;
			integer ic = ii;
			x_fine[jstr] = 0.0;
			while ((ic <= iend) && (P[ic].j == jstr)) {
				x_fine[jstr] += P[ic].aij*x_coarse[P[ic].i];
				ic++;
			}
			
			flag[P[ii].j] = true;
		}
	}

} // prolongation

// Evklid norma of vector r in size n.
Real norma(Real * &r, integer n) {
	Real ret = 0.0;
	for (integer ii = 1; ii <= n; ii++) {
		ret += r[ii] * r[ii] / n;
	}
	ret = sqrt(ret);
	return ret;
}

// экспорт полевой величины u в программу tecplot 360.
void exporttecplot(Real* u, integer n_size) {
	FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s(&fp, "fedorenko1.PLT", "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		if (fp != NULL) {
			// запись имён переменных
			fprintf_s(fp, "VARIABLES = x y u\n");
			fprintf_s(fp, "zone\n");
			integer m = (int)(sqrt((double)(1.0*n_size)));
			integer n = m;
			Real h = 1.0 / (m - 1);
			fprintf_s(fp, "I=%d, J=%d, K=1, F=POINT\n", m, n);
			for (integer j = 0; j < n; j++) for (integer i = 0; i < m; i++)   fprintf_s(fp, "%e %e %e\n", i*h, j*h, u[i*m + j + 1]);
			fclose(fp);
			WinExec("C:\\Program Files\\Tecplot\\Tecplot 360 EX 2014 R1\\bin\\tec360.exe fedorenko1.PLT", SW_NORMAL);
			//WinExec("C:\\Program\ Files\ (x86)\\Tecplot\\Tec360\ 2008\\bin\\tec360.exe fedorenko1.PLT", SW_NORMAL);
		}
		else {
			printf("Create File Error\n");
		}
	}

	system("PAUSE");

} // exporttecplot

// Быстрая сортировка Ч. Хоара. Время n*log2(n).
void quickSort_set(integer* &ifrQ, integer left, integer right)
{
	integer i = left, j = right;
	integer tmp;
	integer pivot = ifrQ[(left + right) / 2];

	/* partition */
	while (i <= j) {
		while (ifrQ[i]<pivot)
			i++;
		while (ifrQ[j]>pivot)
			j--;
		if (i <= j) {
			tmp = ifrQ[i];
			ifrQ[i] = ifrQ[j];
			ifrQ[j] = tmp;
			i++;
			j--;
		}
	};

	/* recursion */
	if (left<j)
		quickSort_set(ifrQ, left, j);
	if (i<right)
		quickSort_set(ifrQ, i, right);


}

// Двоичный поиск.
integer BinarySearch(integer* &A, integer key, integer n)
{
	integer left = 0, right = n, mid;
	while (left <= right)
	{
		mid = left + (right - left) / 2;
		if (key<A[mid]) right = mid - 1;
		else if (key>A[mid]) left = mid + 1;
		else {
			// надо обязательно убедиться что mid самый левый представитель.
			while ((mid > 0) && (A[mid - 1] == A[mid])) {
				mid--;
			}
			return mid;
		}
	}
	return -1;
}
// Двоичный поиск.
integer BinarySearchAi(Ak* &A, integer key, integer istart, integer iend)
{
	integer left = istart, right = iend, mid;
	while (left <= right)
	{
		mid = left + (right - left) / 2;
		if (key<A[mid].i) right = mid - 1;
		else if (key>A[mid].i) left = mid + 1;
		else {
			// надо обязательно убедиться что mid самый левый представитель.
			while ((mid > istart) && (A[mid - 1].i == A[mid].i)) {
				mid--;
			}
			return mid;
		}
	}
	return -1;
}
// Двоичный поиск.
integer BinarySearchAj(Ak* &A, integer key, integer istart, integer iend)
{
	integer left = istart, right = iend, mid;
	while (left <= right)
	{
		mid = left + (right - left) / 2;
		if (key<A[mid].j) right = mid - 1;
		else if (key>A[mid].j) left = mid + 1;
		else {
			// надо обязательно убедиться что mid самый левый представитель.
			while ((mid > istart) && (A[mid - 1].j == A[mid].j)) {
				mid--;
			}
			return mid;
		};
	}
	return -1;
}



// 16 сентября 2015 года обнаружено что операции 
// сгрубления и интерполяции сделаны совершенно неверно,
// и если сгрубление еще в какой-то мере проецирует то интерполяция просто никакая.
// Операции сгрубления и интерполляции будут сделаны заново на основе статьи К.Н. Волкова в новой версии солвера.
// 3 september 2015 Villa Borgese.
integer aggregative_amg(Ak* &A, integer nnz,
	                integer n,
	                Ak* &R, // restriction
	                Ak* &P, // prolongation
	                Real* &x, Real* &b) {

	// количество рекурсивных вызовов ограничено, поэтому QuickSort не подходит.
	bool bquicktypesort = false;


	// x_coarse[1..n_coarse]=R n_coarse x n_fine*x_fine[1..n_fine];
	// x_fine[1..n_fine]=P n_fine x n_coarse*x_coarse[1..n_coarse];

	// нумерация начинается с единицы.

	const integer maxlevel = 30;
	integer ilevel = 1;
	integer nnz_a[maxlevel];
	integer n_a[maxlevel];
	nnz_a[0] = nnz;
	n_a[0] = n;
	bool* flag = new bool[n + 1];
	bool* flag_ = new bool[n + 1];
	integer iadd = 0;
	integer nnzR = 1;
	integer iaddR = 0;
	integer nnz_aRP[maxlevel];
	bool bcontinue = true;

	while ((ilevel<maxlevel-1)&&(n_a[ilevel - 1] > 50)&&(bcontinue)) {

		
		//heapsort(A, key=i*n_a[ilevel - 1] + j, 1, nnz_a[ilevel - 1]);
		if (bquicktypesort) {
			QuickSort(A,  1 + iadd, nnz_a[ilevel - 1] + iadd);
		}
		else {
			if (nnz_a[ilevel - 1] < 100000) {
				HeapSort(A, n_a[ilevel - 1], 1 + iadd, nnz_a[ilevel - 1] + iadd);
			}
			else {
				Ak* Aorig = &A[1 + iadd];
				MergeSort(Aorig, nnz_a[ilevel - 1]);
				Aorig = NULL;
			}
		}

		for (integer ii = 1+iadd; ii <= nnz_a[ilevel - 1]+iadd; ii++) {
		//	printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n",ii,A[ii].aij,ii,A[ii].i,ii,A[ii].j);
			//if (ii % 20 == 0) getchar();
		}

		for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
			flag[ii] = false;
		}
		for (integer ii = 1+iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
			A[ii].ind = ii-iadd;
		}



		// Copy(A) на nnz ячеек правее.
		for (integer ii = nnz_a[ilevel - 1] + 1 + iadd; ii <= 2 * nnz_a[ilevel - 1] + iadd; ii++) {
			A[ii] = A[ii - nnz_a[ilevel -1]];
		}

		// сделаем копию А упорядочим её по j 
		for (integer k = 1 + iadd; k <= nnz_a[ilevel - 1] + iadd; k++) {
			A[k + 2 * nnz_a[ilevel - 1]] = A[k]; // copy
			A[k + 2 * nnz_a[ilevel - 1]].ind = k; // запомиаем номер до сортировки.
		}
		// сортировка по j.
		if (nnz_a[ilevel - 1] < 100000) {
			HeapSort_j(A,/* n_a[ilevel - 1]*/ 1 + iadd + 2 * nnz_a[ilevel - 1], nnz_a[ilevel - 1] + iadd + 2 * nnz_a[ilevel - 1]);
		}
		else {
			Ak* Aorig = &A[1 + iadd + 2 * nnz_a[ilevel - 1]];
			MergeSort_j(Aorig, nnz_a[ilevel - 1]);
			Aorig = NULL;
		}


		integer n_coarce = 1; // номер агрегата.
		nnzR = 1;
		const integer max_sosed = 27850;
		const integer NULL_SOSED = -1;
		integer vacant = NULL_SOSED;
		for (integer ii = 1 + iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
			if (flag[A[ii].i] == false) {
				// Вычисляем по немодифиуцированной матрице А (хранящейся слева).

				Real sum = 0.0;
				integer nnzRl = nnzR + iaddR;

				integer set[max_sosed]; // не более 20 узлов в одном агрегате.
				// инициализация убрана потомучто она не нужна и она сильно тормозит быстродействие.
				//for (integer js = 0; js < max_sosed; js++) {
					//set[js] = NULL_SOSED;
				//}
				integer ic = 0;
				set[ic] = A[ii].i;
				Real theta = 0.25; // контроль числа сильных связей между переменными.
				Real max_vnediagonal = -1.0; // максимальное значение модуля вне диагонального элемента. 
				// добавляем диагональный элемент.
				for (integer is0 = ii; (is0 <= nnz_a[ilevel - 1] + iadd)&&(A[is0].i==set[0]); is0++) {
					if (A[is0].j == set[0]) {
						sum += fabs(A[is0].aij);
						R[nnzRl].aij = fabs(A[is0].aij);
						R[nnzRl].i = n_coarce; // индекс грубой сетки
						R[nnzRl].j = set[0]; //индекс на подробной сетке.
						nnzRl++;
						break;
					}
					else {
						if (fabs(A[is0].aij) > max_vnediagonal) {
							max_vnediagonal = fabs(A[is0].aij);
						}
					}
				}
				
				ic++;

				

				// если узел j ещё не был добавлен в агрегат.
				if (flag[A[ii].j] == false) {
					if ((A[ii].j != set[0]) && (fabs(A[ii].aij) >= theta*max_vnediagonal)) {
						vacant = A[ii].j;
						for (integer js = 0; js < ic; js++) {
							if (vacant == set[js]) {
								vacant = NULL_SOSED;
							}
						}
						if (vacant != NULL_SOSED) {
							set[ic] = vacant;
							sum += fabs(A[ii].aij);
							R[nnzRl].aij = fabs(A[ii].aij);
							R[nnzRl].i = n_coarce;
							R[nnzRl].j = vacant;
							nnzRl++;
							ic++;
						}
					}
				}
				integer iscan = ii + 1;
				while ((iscan <= nnz_a[ilevel - 1] + iadd) && (A[iscan].i == set[0])) {
					// если узел j ещё не был добавлен в агрегат.
					if (flag[A[iscan].j] == false) {
						if ((A[iscan].j != set[0]) && (fabs(A[iscan].aij) >= theta*max_vnediagonal)) {
							vacant = A[iscan].j;
							for (integer js = 0; js < ic; js++) {
								if (vacant == set[js]) {
									vacant = NULL_SOSED;
								}
							}
							if (vacant != NULL_SOSED) {
								set[ic] = vacant;
								sum += fabs(A[iscan].aij);
								R[nnzRl].aij = fabs(A[iscan].aij);
								R[nnzRl].i = n_coarce;
								R[nnzRl].j = vacant;
								nnzRl++;
								ic++;

							}
						}
					}

					iscan++;

				} // while

				// R restriction.
				for (integer k1 = 0; k1 < ic; k1++) {
					R[nnzR+ iaddR].aij /= sum;
					nnzR++;
				}
				
				
				//{
					//integer* Aset = new integer[ic];
					//for (integer ii3 = 0; ii3 < ic; ii3++) Aset[ii3] = set[ii3]; // copy
					// HeapSort(Aset,0,ic-1);
					//quickSort_set(Aset, 0, ic - 1);

					// bynarySearh 88.56%
					// aggregativeamg 8.89%
					// seidel 0.87%
					// (i,j) -> (I,J)
					// модифицируем копию A находящуюся справа.
					//for (integer k = nnz_a[ilevel - 1] + 1 + iadd; k <= 2 * nnz_a[ilevel - 1] + iadd; k++) {
						//bool found = false;

						//if (BinarySearch(Aset, A[k - nnz_a[ilevel - 1]].i, ic - 1) > -1) found = true;
						//for (integer k1 = 0; k1 < ic; k1++) {
						//	if (A[k - nnz_a[ilevel - 1]].i == set[k1]) found = true;
						//	}
						//if (found) A[k].i = n_coarce;
						//found = false;
						//if (BinarySearch(Aset, A[k - nnz_a[ilevel - 1]].j, ic - 1) > -1) found = true;
						
						//for (integer k1 = 0; k1 < ic; k1++) {
						//if (A[k - nnz_a[ilevel - 1]].j == set[k1]) found = true;
						//}
						
						//if (found) A[k].j = n_coarce;

					//}

					//delete[] Aset;
				//}
			
				
				{
					// 7 сентября 2015.
					// мы воспользовались тем что А упорядочена по i.
					for (integer k1 = 0; k1 < ic; k1++) {
						integer key = set[k1];
						integer ifound=BinarySearchAi(A, key, 1 + iadd, nnz_a[ilevel - 1] + iadd);
						if (ifound > -1) {
							integer if1 = ifound;
							while ((if1 <= nnz_a[ilevel - 1] + iadd) && (A[if1].i == key)){
								A[if1  + nnz_a[ilevel - 1]].i = n_coarce;
								if1++;
							}
							if1 = ifound-1;
							while ((if1 >= 1 + iadd) && (A[if1].i == key)){
								A[if1 +nnz_a[ilevel - 1]].i = n_coarce;
								if1--;
							}

						}
					}
					// копирование А на третью позицию и сортировка по j. Делается единожды в самом начале.
					// и сделаем тоже самое для j.
					// мы воспользовались тем что А упорядочена по j.
					for (integer k1 = 0; k1 < ic; k1++) {
						integer key = set[k1];
						integer ifound = BinarySearchAj(A, key, 1 + iadd + 2 * nnz_a[ilevel - 1], nnz_a[ilevel - 1] + iadd + 2 * nnz_a[ilevel - 1]);
						if (ifound > -1) {
							integer if1 = ifound;
							while ((if1 <= nnz_a[ilevel - 1] + iadd + 2 * nnz_a[ilevel - 1]) && (A[if1].j == key)) {
								A[A[if1].ind  + nnz_a[ilevel - 1]].j = n_coarce;
								if1++;
							}
							if1 = ifound - 1;
							while ((if1 >= 1 + iadd + 2 * nnz_a[ilevel - 1]) && (A[if1].j == key)) {
								A[A[if1].ind + nnz_a[ilevel - 1]].j = n_coarce;
								if1--;
							}

						}
					}
					// Временная память в разделе 1 + iadd + 2 * nnz_a[ilevel - 1]..nnz_a[ilevel - 1] + iadd + 2 * nnz_a[ilevel - 1]
					// больше ненужна.

				}
				

				// Помечаем узлы как включённые в агрегат.
				for (integer js = 0; js < ic; js++) {
					flag[set[js]] = true;
				}

				n_coarce++;
				
				// Один агрегат создан.



			} // узел не был ещё включён в агрегат.
		} // агрегаты созданы.


		// отладочная печать в рабочей версии требуется закоментировать.
		//for (integer ii = 1; ii <= n; ii++) {
			//flag_[ii] = false;
		//}
		//for (integer ii = 1 + iaddR; ii <= iaddR + nnzR - 1; ii++) {
			//if (flag_[R[ii].i] == false) {
				//integer istr = R[ii].i;
				//integer ic7 = ii;
				//while ((ic7 <= iaddR + nnzR - 1) && (R[ic7].i == istr)) {
					//printf("%e ", R[ic7].aij);
					//ic7++;
				//}
				//printf("\n");
				//system("pause");
				//flag_[R[ii].i] = true;
			//}
		//}


		// оператор restriction построен и он упорядочен по i.
		// число ненулевых элементов nnzR-1.
		// P=Copy(R);
		for (integer ii = 1 + iaddR; ii <= iaddR + nnzR - 1; ii++) {
			P[ii] = R[ii];
			//printf("ii=%d aij=%e, i=%d j=%d\n",ii,P[ii].aij,P[ii].i,P[ii].j);
			//getchar();
		}
		
		// heapsort(P,key==j,iaddR+1,iaddR+nnzR - 1);
		if ( bquicktypesort) {
			QuickSort_j(P, 1 + iaddR, iaddR + nnzR - 1);
		}
		else {
		    HeapSort_j(P, /*n_a[ilevel-1]*/ 1 + iaddR, iaddR + nnzR - 1);
		}

		// оператор интерполяции это не просто транспонированный оператор проекции а
		// а транспонированный оператор проекции умноженный на константу. Константа 
		// определяется из следующего соображения : если сумма элементов оператора рестрикции в стоке единица,
		// то соответственно в столбце у оператора интерполяции максимальный элемент равен единица.
		// этот код обязательно должен быть включён чтобы пара рестрикция-интерполяция была верна.
		
		for (integer ii = 1; ii <= n; ii++) {
			flag_[ii] = false;
		}
		Real mul = -1.e30;
		for (integer ii = 1 + iaddR; ii <= iaddR + nnzR - 1; ii++) {
			if (flag_[P[ii].j] == false) {
				integer jstr = P[ii].j;
				integer ic7 = ii;
				mul = -1.e30;
				while ((ic7 <= iaddR + nnzR - 1) && (P[ic7].j == jstr)) {
					if (fabs(P[ic7].aij) > mul) {
						mul = fabs(P[ic7].aij);
					}
					ic7++;
				}
				ic7 = ii;
				while ((ic7 <= iaddR + nnzR - 1) && (P[ic7].j == jstr)) {
					P[ic7].aij /= mul; // максимальный элемент единица.
					ic7++;
				}
				flag_[P[ii].j] = true;
			}
		}
		
		//for (integer ii = 1; ii <= n; ii++) {
			//flag_[ii] = false;
		//}
		//for (integer ii = 1 + iaddR; ii <= iaddR + nnzR - 1; ii++) {
			//if (flag_[P[ii].j] == false) {
				//integer jstr = P[ii].j;
				//integer ic7 = ii;
				//while ((ic7 <= iaddR + nnzR - 1) && (P[ic7].j == jstr)) {
					//printf("%e ",P[ic7].aij);
					//ic7++;
				//}
				//printf("\n");
				//system("pause");
				//flag_[P[ii].j] = true;
			//}
		//}
		

		
		
		nnz_aRP[ilevel - 1] = nnzR - 1;
		iaddR += nnzR - 1;

		//printf("%d %d\n",n,n_coarce-1);
		//getchar();
		n_a[ilevel] = n_coarce - 1;

		// сортировка по новому ключу key=i*(iglcoarce_number-1)+j;
		// в позиции ind сохранён индекс предыдущей позиции.
		//heapsort(A, key = i*(n_coarce - 1) + j, nnz + 1, 2 * nnz);
		if (bquicktypesort) {
			QuickSort(A, /*n_coarce - 1,*/ nnz_a[ilevel - 1] + 1 + iadd, 2 * nnz_a[ilevel - 1] + iadd);
		}
		else {
			if (nnz_a[ilevel - 1] < 100000) {
				HeapSort(A,  n_coarce - 1, nnz_a[ilevel - 1] + 1 + iadd, 2 * nnz_a[ilevel - 1] + iadd);
			}
			else {
				Ak* Aorig = &A[nnz_a[ilevel - 1] + 1 + iadd];
				MergeSort(Aorig, nnz_a[ilevel - 1]);
				Aorig = NULL;
			}
		}

		for (integer ii = 1+nnz_a[ilevel-1]+iadd; ii <= 2*nnz_a[ilevel - 1]+iadd; ii++) {
			//	printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n",ii,A[ii].aij,ii,A[ii].i,ii,A[ii].j);
			//if (ii % 20 == 0) getchar();
		}



		// копирование в третью часть матрицы Слау на coarce сетке.
		// инициализация.
		// запас памяти 1 nnz на fine уровне.
		for (integer ii = 2 * nnz_a[ilevel - 1]+1 + iadd; ii <= 3 * nnz_a[ilevel - 1] + iadd; ii++) {
			A[ii].aij = 0.0;
			A[ii].ind = NULL_SOSED;
			A[ii].i = NULL_SOSED;
			A[ii].j = NULL_SOSED;
		}

		for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) flag[ii] = false;


		integer ic1 = 2 * nnz_a[ilevel - 1] + 1 + iadd;
		integer im = 1;
		integer im0 = 0;
		for (integer ii = nnz_a[ilevel - 1] + 1 + iadd; ii <= 2 * nnz_a[ilevel - 1] + iadd; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				while ((ii + im0 <= 2 * nnz_a[ilevel - 1] + iadd)&&(istr==A[ii+im0].i)) {

					if (ic1 <= 3 * nnz_a[ilevel - 1] + iadd) {
						A[ic1].i = A[ii + im0].i;
						A[ic1].j = A[ii + im0].j;
						A[ic1].aij += A[ii + im0].aij;
						while ((ii + im <= 2 * nnz_a[ilevel - 1] + iadd) && (A[ii + im0].i == A[ii + im].i) && (A[ii + im0].j == A[ii + im].j))
				        {
					         A[ic1].aij += A[ii + im].aij;
					         im++;
				        }
						ic1++;
						im0 = im;
						im++;
					}
					else {
						printf("error 1\n");
						system("PAUSE");
					}
				}
				flag[A[ii].i] = true;
				im = 1;
				im0 = 0;
			}
		}
		
		nnz_a[ilevel] = ic1 - 1 - 2 * nnz_a[ilevel - 1] - iadd;
		iadd += 2 * nnz_a[ilevel - 1];

		printf("nnz : fine=%d, coarse=%d, operator complexity=%e. \n", nnz_a[ilevel - 1], nnz_a[ilevel], (double)(nnz_a[ilevel]) / (double)(nnz_a[ilevel - 1]));
	    printf("n : fine=%d, coarse=%d grid complexity=%e.\n",   n_a[ilevel - 1], n_a[ilevel], (double)(n_a[ilevel]) / (double)(n_a[ilevel - 1]));
		printf("nnz_aRP = %d\n",nnz_aRP[ilevel-1]);
		//getchar();

		ilevel++;// грубосеточная матрица построена.
		if (ilevel >= 2) {
			if (n_a[ilevel - 2] <= n_a[ilevel - 1]) bcontinue = false;
		}

	} // иерархия сеток построена.

	//for (integer ii = 1; ii <= nnz_aRP[0] + nnz_aRP[1]; ii++) {
		//printf("ii=%d aij=%e, i=%d j=%d\n", ii, R[ii].aij, R[ii].i, R[ii].j);
		//getchar();
	//}


	//exporttecplot(b,n);

	//Real *test_coarse = new Real[n_a[1] + 1];

	// restriction
	//restriction(R, 1, nnz_aRP[0], flag, b, test_coarse, n_a[0], n_a[1]);
	//for (integer ii = 1; ii <= n_a[0]; ii++) {
		//b[ii] = 0.0;
	//}

	//{
		//Real *test_coarse1 = new Real[n_a[2] + 1];

		// restriction
		//restriction(R, 1+nnz_aRP[0], nnz_aRP[0]+nnz_aRP[1], flag, test_coarse, test_coarse1, n_a[1], n_a[2]);
		//for (integer ii = 1; ii <= n_a[1]; ii++) {
			//test_coarse[ii] = 0.0;
		//}

		//{
		//	Real *test_coarse2 = new Real[n_a[3] + 1];

			// restriction
			//restriction(R, 1 + nnz_aRP[0]+nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1]+nnz_aRP[2], flag, test_coarse1, test_coarse2, n_a[2], n_a[3]);
			//for (integer ii = 1; ii <= n_a[2]; ii++) {
				//test_coarse1[ii] = 0.0;
			//}

			//prolongation(P, 1 + nnz_aRP[0]+ nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], flag, test_coarse1, test_coarse2, n_a[2], n_a[3]);
		//}

		//prolongation(P, 1 + nnz_aRP[0], nnz_aRP[0] + nnz_aRP[1], flag, test_coarse, test_coarse1, n_a[1], n_a[2]);
	//}

	//prolongation(P, 1, nnz_aRP[0], flag, b, test_coarse, n_a[0], n_a[1]);

	//exporttecplot(b, n);
	

	// подготовка матрицы к cycling:
	
	// smoother.
// 1 september 2015.
//void seidel(Ak* &A, integer istart, integer iend, Real* &x, Real* &b, bool* &flag, integer n)
//{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	//for (integer i = 1; i <= n; i++) {
		//flag[i] = false;
	//}
	//for (integer ii = istart; ii <= iend; ii++) {
		//if (flag[A[ii].i] == false) {
			//integer istr = A[ii].i;
			//integer ic = ii;
			//Real ap = 0.0;
			//x[istr] = b[istr];
			//while ((ic<=iend)&&(A[ic].i == istr)) {
				//if (A[ic].j != istr) {
					//x[istr] += -A[ic].aij*x[A[ic].j];
				//}
				//else ap = A[ic].aij;
				//ic++;
			//}
			//if (fabs(ap) < 1.0e-30) {
				//printf("zero diagonal elements in string %d",istr);
				//getchar();
				//exit(1);
			//}
			//else {
				//x[istr] /= ap;
			//}
			//flag[A[ii].i] = true;
		//}
	//}


//} // seidel

	
	integer *row_ptr_start = new integer[4*n_a[0]+1];
	integer *row_ptr_end = new integer[4 * n_a[0] + 1];
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	for (integer ii = 1; ii <= nnz_a[0]; ii++) {
		if (flag[A[ii].i] == false) {
			integer istr = A[ii].i;
			integer ic = ii;
			integer icdiag = ii;
			row_ptr_start[istr] = ii;
			Real ap = 0.0;
			//x[istr] = b[istr];
			while ((ic <= nnz_a[0]) && (A[ic].i == istr)) {
				if (A[ic].j != istr) {
					//x[istr] += -A[ic].aij*x[A[ic].j];
				}
				else {
					ap = A[ic].aij;
					icdiag = ic;
				}
				ic++;
			}
			row_ptr_end[istr] = ic - 1;
			if (fabs(ap) < 1.0e-30) {
				printf("zero diagonal elements in string %d", istr);
				system("PAUSE");
				exit(1);
			}
			else {
				//x[istr] /= ap;
			}
			
			flag[A[ii].i] = true;
			Ak temp = A[ii];
			A[ii] = A[icdiag];
			A[icdiag] = temp;
			A[ii].aij = 1.0 / ap; // умножение быстрей деления.
		}
	}
	// первый уровень вложенности.
	if (ilevel > 1) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		for (integer ii = 2*nnz_a[0] + 1; ii <= 2*nnz_a[0]+nnz_a[1]; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr+n_a[0]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <= 2 * nnz_a[0] + nnz_a[1]) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr+n_a[0]] = ic - 1;
				if (fabs(ap) < 1.0e-30) {
					printf("zero diagonal elements in string %d", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}

	// второй уровень вложенности.
	
	if (ilevel > 2) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		for (integer ii = 2 * nnz_a[0]+2*nnz_a[1] + 1; ii <= 2 * nnz_a[0]+ 2*nnz_a[1] + nnz_a[2]; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0]+n_a[1]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <= 2 * nnz_a[0] + 2*nnz_a[1]+nnz_a[2]) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0]+n_a[1]] = ic - 1;
				if (fabs(ap) < 1.0e-30) {
					printf("zero diagonal elements in string %d", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}
	

	// третий уровень вложенности.

	if (ilevel > 3) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		for (integer ii = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 1; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + nnz_a[3]; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0] + n_a[1]+n_a[2]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + nnz_a[3]) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0] + n_a[1] + n_a[2]] = ic - 1;
				if (fabs(ap) < 1.0e-30) {
					printf("zero diagonal elements in string %d", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}


	// 14 сентября 2015 понедельник
	// четвёртый уровень вложенности.

	if (ilevel > 4) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		integer ist = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 1;
		integer iend = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + nnz_a[4];
		for (integer ii = ist; ii <= iend; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <= iend ) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3]] = ic - 1;
				if (fabs(ap) < 1.0e-30) {
					printf("zero diagonal elements in string %d", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}


	// пятый уровень вложенности.

	if (ilevel > 5) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		integer ist = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 1;
		integer iend = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + nnz_a[5];
		for (integer ii = ist; ii <= iend; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <= iend) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3]+n_a[4]] = ic - 1;
				if (fabs(ap) < 1.0e-30) {
					printf("zero diagonal elements in string %d", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}

	// шестой уровень вложенности.

	if (ilevel > 6) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		integer ist = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 1;
		integer iend = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + nnz_a[6];
		for (integer ii = ist; ii <= iend; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <= iend) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5]] = ic - 1;
				if (fabs(ap) < 1.0e-30) {
					printf("zero diagonal elements in string %d", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}


	
	
	// smoother.
// 9 september 2015.
// q - quick.
// seidelq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);
//void seidelq(Ak* &A, integer istartq, integer iendq, Real* &x, Real* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
//{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	//integer startpos = istartq + iadd;
	//integer endpos = iendq+iadd;
	//for (integer ii = startpos; ii <= endpos; ii++) {
		//integer istr = ii - iadd;
		//x[istr] = b[istr];
		//for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++)
		//{
		//x[istr] += -A[ii1].aij*x[A[ii1].j]; 
		//}
		//x[istr] *= A[row_ptr_start[ii]].aij;
	//}


//} // seidelq
	


	printf("cycling: V cycle.\n");
	printf("level=%d\n", ilevel);
	printf("multigrid R.P.Fedorenko 1961.\n");
	printf("aggregative algebraic multigrid method.\n");
	//system("pause");

	// 10 11 21 multigrid tutorial Вильм Бригг.

	integer nu1 = 4;
	integer nu2 = 3;

	//ilevel = 1; //debug
	Real rho = 1.0;
	Real dres = 1.0;
	int iiter = 1;
	const Real tolerance = 1.0e-12;


	Real *residual_fine = new Real[n_a[0] + 1];
	Real *residual_coarse = NULL;
	Real* error_approx_coarse = NULL;
	Real *residual_fine1 = NULL;
	Real *residual_coarse1 = NULL;
	Real* error_approx_coarse1 = NULL;
	Real *error_approx_fine1 = NULL;
	Real *residual_fine2 = NULL;
	Real *residual_coarse2 = NULL;
	Real* error_approx_coarse2 = NULL;
	Real *error_approx_fine2 = NULL;
	Real *residual_fine3 = NULL;
	Real *residual_coarse3 = NULL;
	Real* error_approx_coarse3 = NULL;
	Real *error_approx_fine3 = NULL;
	Real *residual_fine4 = NULL;
	Real *residual_coarse4 = NULL;
	Real *error_approx_coarse4 = NULL;
	Real *error_approx_fine4 = NULL;
	Real *residual_fine5 = NULL;
	Real *residual_coarse5 = NULL;
	Real* error_approx_coarse5 = NULL;
	Real *error_approx_fine5 = NULL;
	Real *residual_fine6 = NULL;
	Real *residual_coarse6 = NULL;
	Real* error_approx_coarse6 = NULL;
	Real *error_approx_fine6 = NULL;

	if (ilevel > 1) {
		residual_coarse = new Real[n_a[1] + 1];
		error_approx_coarse = new Real[n_a[1] + 1];
		if (ilevel > 2) {
			// residual
			residual_fine1 = new Real[n_a[1] + 1];
			residual_coarse1 = new Real[n_a[2] + 1];
			error_approx_coarse1 = new Real[n_a[2] + 1];
			error_approx_fine1 = new Real[n_a[1] + 1];
			if (ilevel > 3) {
				// residual
				residual_fine2 = new Real[n_a[2] + 1];
				residual_coarse2 = new Real[n_a[3] + 1];
				error_approx_coarse2 = new Real[n_a[3] + 1];
				error_approx_fine2 = new Real[n_a[2] + 1];
				if (ilevel > 4) {
					// residual
					residual_fine3 = new Real[n_a[3] + 1];
					residual_coarse3 = new Real[n_a[4] + 1];
					error_approx_coarse3 = new Real[n_a[4] + 1];
					error_approx_fine3 = new Real[n_a[3] + 1];
					if (ilevel > 5) {
						// residual
						residual_fine4 = new Real[n_a[4] + 1];
						residual_coarse4 = new Real[n_a[5] + 1];
						error_approx_coarse4 = new Real[n_a[5] + 1];
						error_approx_fine4 = new Real[n_a[4] + 1];
						if (ilevel > 6) {
							// residual
							residual_fine5 = new Real[n_a[5] + 1];
							residual_coarse5 = new Real[n_a[6] + 1];
							error_approx_coarse5 = new Real[n_a[6] + 1];
							error_approx_fine5 = new Real[n_a[5] + 1];
							if (ilevel > 7) {
								// residual
								residual_fine6 = new Real[n_a[6] + 1];
								residual_coarse6 = new Real[n_a[7] + 1];
								error_approx_coarse6 = new Real[n_a[7] + 1];
								error_approx_fine6 = new Real[n_a[6] + 1];
							}
						}
					}
				}
			}
		}
	}
	Real *error_approx_fine = new Real[n_a[0] + 1];


	//for (integer iprohod = 0; iprohod < 20; iprohod++) {
	while (dres>tolerance) {

		// smother
		for (integer iter = 0; iter < nu1; iter++) {
			//seidel(A, 1, nnz_a[0], x, b, flag, n_a[0]);
			//quick seidel
			seidelq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);

		}
		
		//exporttecplot(x, n);

		// residual_r
		//Real *residual_fine = new Real[n_a[0] + 1];
		//residual(A, 1, nnz_a[0], x, b, flag, n_a[0], residual_fine);
		residualq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine);
		dres = norma(residual_fine, n_a[0]);
		printf("%d %e rho=%e\n",iiter, dres, dres/rho);
		iiter++;
		//rho=norma(residual_fine, n_a[0]);
		rho = dres;
		//if (iprohod%5==0) getchar();
		if (ilevel > 1) {

			//Real *residual_coarse = new Real[n_a[1] + 1];

			// restriction
			restriction(R, 1, nnz_aRP[0], flag, residual_fine, residual_coarse, n_a[0], n_a[1]);

			// A*e=r;
			//Real* error_approx_coarse = new Real[n_a[1] + 1];
			for (integer ii = 1; ii <= n_a[1]; ii++) {
				error_approx_coarse[ii] = 0.0;
			}
			// pre smothing
			for (integer iter = 0; iter < nu1; iter++) {
				//seidel(A, 1 + 2 * nnz_a[0], 2 * nnz_a[0] + nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1]);
				seidelq(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0]);
			}
			
			if (ilevel > 2) {
				// residual
				//Real *residual_fine1 = new Real[n_a[1] + 1];
				//residual(A, 1+2*nnz_a[0], 2*nnz_a[0]+nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1], residual_fine1);
				//residualq(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0], residual_fine1);
				residualqspeshial(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0], residual_fine1);


				//Real *residual_coarse1 = new Real[n_a[2] + 1];

				// restriction
				restriction(R, 1+nnz_aRP[0], nnz_aRP[0]+nnz_aRP[1], flag, residual_fine1, residual_coarse1, n_a[1], n_a[2]);
			
				// A*e=r;
				//Real* error_approx_coarse1 = new Real[n_a[2] + 1];
				for (integer ii = 1; ii <= n_a[2]; ii++) {
					error_approx_coarse1[ii] = 0.0;
				}
				// pre smothing
				for (integer iter = 0; iter < nu1; iter++) {
					//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1], 2 * nnz_a[0] + 2 * nnz_a[1] + nnz_a[2], error_approx_coarse1, residual_coarse1, flag, n_a[2]);
					seidelq(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0]+n_a[1]);
				}
				if (ilevel > 3) {
					// residual
					//Real *residual_fine2 = new Real[n_a[2] + 1];
					//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1], 2 * nnz_a[0] + 2 * nnz_a[1] + nnz_a[2], error_approx_coarse1, residual_coarse1, flag, n_a[2], residual_fine2);
					//residualq(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0]+n_a[1], residual_fine2);
					residualqspeshial(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0] + n_a[1], residual_fine2);

					//Real *residual_coarse2 = new Real[n_a[3] + 1];

					// restriction
					restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], flag, residual_fine2, residual_coarse2, n_a[2], n_a[3]);

					// A*e=r;
					//Real* error_approx_coarse2 = new Real[n_a[3] + 1];
					for (integer ii = 1; ii <= n_a[3]; ii++) {
						error_approx_coarse2[ii] = 0.0;
					}
					// pre smothing
					for (integer iter = 0; iter < nu1; iter++) {
						//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + nnz_a[3], error_approx_coarse2, residual_coarse2, flag, n_a[3]);
						seidelq(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1]+n_a[2]);

					}
					if (ilevel > 4) {
						// residual
						//Real *residual_fine3 = new Real[n_a[3] + 1];
						//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + nnz_a[3], error_approx_coarse2, residual_coarse2, flag, n_a[3], residual_fine3);
						//residualq(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2], residual_fine3);
						//speshial
						residualqspeshial(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2], residual_fine3);



						//Real *residual_coarse3 = new Real[n_a[4] + 1];

						// restriction
						restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], flag, residual_fine3, residual_coarse3, n_a[3], n_a[4]);

						// A*e=r;
						//Real* error_approx_coarse3 = new Real[n_a[4] + 1];
						for (integer ii = 1; ii <= n_a[4]; ii++) {
							error_approx_coarse3[ii] = 0.0;
						}
						// pre smothing
						for (integer iter = 0; iter < nu1; iter++) {
							//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + nnz_a[4], error_approx_coarse3, residual_coarse3, flag, n_a[4]);
							seidelq(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3]);
						}
						if (ilevel > 5) {
							// residual
							//Real *residual_fine4 = new Real[n_a[4] + 1];
							//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + nnz_a[4], error_approx_coarse3, residual_coarse3, flag, n_a[4], residual_fine4);
							//residualq(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3], residual_fine4);
							//speshial 14 september 2015.
							residualqspeshial(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3], residual_fine4);


							//Real *residual_coarse4 = new Real[n_a[5] + 1];

							// restriction
							restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3]+ nnz_aRP[4], flag, residual_fine4, residual_coarse4, n_a[4], n_a[5]);

							// A*e=r;
							//Real* error_approx_coarse4 = new Real[n_a[5] + 1];
							for (integer ii = 1; ii <= n_a[5]; ii++) {
								error_approx_coarse4[ii] = 0.0;
							}
							// pre smothing
							for (integer iter = 0; iter < nu1; iter++) {
								//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + nnz_a[5], error_approx_coarse4, residual_coarse4, flag, n_a[5]);
								seidelq(A, 1, n_a[5], error_approx_coarse4, residual_coarse4, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4]);
							}
							if (ilevel > 6) {
								// residual
								//Real *residual_fine5 = new Real[n_a[5] + 1];
								//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + nnz_a[5], error_approx_coarse4, residual_coarse4, flag, n_a[5], residual_fine5);
								//if (ilevel <= 15) {
									residualq(A, 1, n_a[5], error_approx_coarse4, residual_coarse4, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4], residual_fine5);
								//}
								//else {
									// приводит к расходимости.
									//speshial 14 september 2015.
									// это уже приводит к увеличению числа итераций на примере сетки в 1млн узлов. остановимся.
									//residualqspeshial(A, 1, n_a[5], error_approx_coarse4, residual_coarse4, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4], residual_fine5);
								//}

								//Real *residual_coarse5 = new Real[n_a[6] + 1];

								// restriction
								restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], flag, residual_fine5, residual_coarse5, n_a[5], n_a[6]);

								// A*e=r;
								//Real* error_approx_coarse5 = new Real[n_a[6] + 1];
								for (integer ii = 1; ii <= n_a[6]; ii++) {
									error_approx_coarse5[ii] = 0.0;
								}
								// pre smothing
								for (integer iter = 0; iter < nu1; iter++) {
									//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + nnz_a[6], error_approx_coarse5, residual_coarse5, flag, n_a[6]);
									seidelq(A, 1, n_a[6], error_approx_coarse5, residual_coarse5, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5]);
								}

								if (ilevel > 7) {
									// residual
									//Real *residual_fine6 = new Real[n_a[6] + 1];
									//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] +2*nnz_a[5], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + nnz_a[6], error_approx_coarse5, residual_coarse5, flag, n_a[6], residual_fine6);
									residualq(A, 1, n_a[6], error_approx_coarse5, residual_coarse5, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5], residual_fine6);

									//Real *residual_coarse6 = new Real[n_a[7] + 1];

									// restriction
									restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5]+nnz_aRP[6], flag, residual_fine6, residual_coarse6, n_a[6], n_a[7]);

									// A*e=r;
									//Real* error_approx_coarse6 = new Real[n_a[7] + 1];
									for (integer ii = 1; ii <= n_a[7]; ii++) {
										error_approx_coarse6[ii] = 0.0;
									}
									// pre smothing
									for (integer iter = 0; iter < nu1; iter++) {
										seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + nnz_a[7], error_approx_coarse6, residual_coarse6, flag, n_a[7]);
									}

									if (ilevel > 8) {
										// residual
										Real *residual_fine7 = new Real[n_a[7] + 1];
										residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5]+2*nnz_a[6], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] +2*nnz_a[6]+ nnz_a[7], error_approx_coarse6, residual_coarse6, flag, n_a[7], residual_fine7);


										Real *residual_coarse7 = new Real[n_a[8] + 1];

										// restriction
										restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6]+nnz_aRP[7], flag, residual_fine7, residual_coarse7, n_a[7], n_a[8]);

										// A*e=r;
										Real* error_approx_coarse7 = new Real[n_a[8] + 1];
										for (integer ii = 1; ii <= n_a[8]; ii++) {
											error_approx_coarse7[ii] = 0.0;
										}
										// pre smothing
										for (integer iter = 0; iter < nu1; iter++) {
											seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6]+2*nnz_a[7], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] +2*nnz_a[7]+ nnz_a[8], error_approx_coarse7, residual_coarse7, flag, n_a[8]);
										}


										if (ilevel > 9) {
											// residual
											Real *residual_fine8 = new Real[n_a[8] + 1];
											integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7];
											integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + nnz_a[8];
											residual(A, n1, n2, error_approx_coarse7, residual_coarse7, flag, n_a[8], residual_fine8);


											Real *residual_coarse8 = new Real[n_a[9] + 1];

											// restriction
											integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6]+nnz_aRP[7];
											integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7]+nnz_aRP[8];
											restriction(R,n3 ,n4 , flag, residual_fine8, residual_coarse8, n_a[8], n_a[9]);

											// A*e=r;
											Real* error_approx_coarse8 = new Real[n_a[9] + 1];
											for (integer ii = 1; ii <= n_a[9]; ii++) {
												error_approx_coarse8[ii] = 0.0;
											}
											// pre smothing
											for (integer iter = 0; iter < nu1; iter++) {
												integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8];
												integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + nnz_a[9];
												seidel(A, n5 , n6, error_approx_coarse8, residual_coarse8, flag, n_a[9]);
											}

											if (ilevel > 10) {
												// 8 сентября 2015 РИМИНИ пляж 

												// residual
												Real *residual_fine9 = new Real[n_a[9] + 1];
												integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7]+2*nnz_a[8];
												integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] +2*nnz_a[8]+ nnz_a[9];
												residual(A, n1, n2, error_approx_coarse8, residual_coarse8, flag, n_a[9], residual_fine9);


												Real *residual_coarse9 = new Real[n_a[10] + 1];

												// restriction
												integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7]+nnz_aRP[8];
												integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8]+nnz_aRP[9];
												restriction(R, n3, n4, flag, residual_fine9, residual_coarse9, n_a[9], n_a[10]);

												// A*e=r;
												Real* error_approx_coarse9 = new Real[n_a[10] + 1];
												for (integer ii = 1; ii <= n_a[10]; ii++) {
													error_approx_coarse9[ii] = 0.0;
												}
												// pre smothing
												for (integer iter = 0; iter < nu1; iter++) {
													integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8]+2*nnz_a[9];
													integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] +2*nnz_a[9]+ nnz_a[10];
													seidel(A, n5, n6, error_approx_coarse9, residual_coarse9, flag, n_a[10]);
												}

												if (ilevel > 11) {
													// 8 сентября 2015 РИМИНИ пляж 

													// residual
													Real *residual_fine10 = new Real[n_a[10] + 1];
													integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8]+2*nnz_a[9];
													integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] +2*nnz_a[9]+ nnz_a[10];
													residual(A, n1, n2, error_approx_coarse9, residual_coarse9, flag, n_a[10], residual_fine10);


													Real *residual_coarse10 = new Real[n_a[11] + 1];

													// restriction
													integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8]+nnz_aRP[9];
													integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9]+nnz_aRP[10];
													restriction(R, n3, n4, flag, residual_fine10, residual_coarse10, n_a[10], n_a[11]);

													// A*e=r;
													Real* error_approx_coarse10 = new Real[n_a[11] + 1];
													for (integer ii = 1; ii <= n_a[11]; ii++) {
														error_approx_coarse10[ii] = 0.0;
													}
													// pre smothing
													for (integer iter = 0; iter < nu1; iter++) {
														integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9]+2*nnz_a[10];
														integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2*nnz_a[10] + nnz_a[11];
														seidel(A, n5, n6, error_approx_coarse10, residual_coarse10, flag, n_a[11]);
													}

													if (ilevel > 12) {
														// 11 сентября 2015 РИМИНИ пляж 

														// residual
														Real *residual_fine11 = new Real[n_a[11] + 1];
														integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9]+2*nnz_a[10];
														integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] +2*nnz_a[10]+ nnz_a[11];
														residual(A, n1, n2, error_approx_coarse10, residual_coarse10, flag, n_a[11], residual_fine11);


														Real *residual_coarse11 = new Real[n_a[12] + 1];

														// restriction
														integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9]+nnz_aRP[10];
														integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] +nnz_aRP[10]+ nnz_aRP[11];
														restriction(R, n3, n4, flag, residual_fine11, residual_coarse11, n_a[11], n_a[12]);

														// A*e=r;
														Real* error_approx_coarse11 = new Real[n_a[12] + 1];
														for (integer ii = 1; ii <= n_a[12]; ii++) {
															error_approx_coarse11[ii] = 0.0;
														}
														// pre smothing
														for (integer iter = 0; iter < nu1; iter++) {
															integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10]+2*nnz_a[11];
															integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2*nnz_a[10] +2*nnz_a[11]+ nnz_a[12];
															seidel(A, n5, n6, error_approx_coarse11, residual_coarse11, flag, n_a[12]);
														}

														if (ilevel > 13) {
															// 11 сентября 2015 РИМИНИ пляж 

															// residual
															Real *residual_fine12 = new Real[n_a[12] + 1];
															integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11];
															integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + nnz_a[12];
															residual(A, n1, n2, error_approx_coarse11, residual_coarse11, flag, n_a[12], residual_fine12);


															Real *residual_coarse12 = new Real[n_a[13] + 1];

															// restriction
															integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11];
															integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12];
															restriction(R, n3, n4, flag, residual_fine12, residual_coarse12, n_a[12], n_a[13]);

															// A*e=r;
															Real* error_approx_coarse12 = new Real[n_a[13] + 1];
															for (integer ii = 1; ii <= n_a[13]; ii++) {
																error_approx_coarse12[ii] = 0.0;
															}
															// pre smothing
															for (integer iter = 0; iter < nu1; iter++) {
																integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12];
																integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + nnz_a[13];
																seidel(A, n5, n6, error_approx_coarse12, residual_coarse12, flag, n_a[13]);
															}


															if (ilevel > 14) {
																// 11 сентября 2015 РИМИНИ пляж 

																// residual
																Real *residual_fine13 = new Real[n_a[13] + 1];
																integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12];
																integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + nnz_a[13];
																residual(A, n1, n2, error_approx_coarse12, residual_coarse12, flag, n_a[13], residual_fine13);


																Real *residual_coarse13 = new Real[n_a[14] + 1];

																// restriction
																integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11]+nnz_aRP[12];
																integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12]+nnz_aRP[13];
																restriction(R, n3, n4, flag, residual_fine13, residual_coarse13, n_a[13], n_a[14]);

																// A*e=r;
																Real* error_approx_coarse13 = new Real[n_a[14] + 1];
																for (integer ii = 1; ii <= n_a[14]; ii++) {
																	error_approx_coarse13[ii] = 0.0;
																}
																// pre smothing
																for (integer iter = 0; iter < nu1; iter++) {
																	integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13];
																	integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13] + nnz_a[14];
																	seidel(A, n5, n6, error_approx_coarse13, residual_coarse13, flag, n_a[14]);
																}

																if (ilevel > 15) {
																	// 14 сентября 2015 Москва на работе в пн. 

																	// residual
																	Real *residual_fine14 = new Real[n_a[14] + 1];
																	integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13];
																	integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13] + nnz_a[14];
																	residual(A, n1, n2, error_approx_coarse13, residual_coarse13, flag, n_a[14], residual_fine14);


																	Real *residual_coarse14 = new Real[n_a[15] + 1];

																	// restriction
																	integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13];
																	integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14];
																	restriction(R, n3, n4, flag, residual_fine14, residual_coarse14, n_a[14], n_a[15]);

																	// A*e=r;
																	Real* error_approx_coarse14 = new Real[n_a[15] + 1];
																	for (integer ii = 1; ii <= n_a[15]; ii++) {
																		error_approx_coarse14[ii] = 0.0;
																	}
																	// pre smothing
																	for (integer iter = 0; iter < nu1; iter++) {
																		integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13] + 2 * nnz_a[14];
																		integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13] + 2 * nnz_a[14] + nnz_a[15];
																		seidel(A, n5, n6, error_approx_coarse14, residual_coarse14, flag, n_a[15]);
																	}

																	// post smothing
																	for (integer iter = 0; iter < nu2; iter++) {
																		integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13] + 2 * nnz_a[14];
																		integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13] + 2 * nnz_a[14] + nnz_a[15];
																		seidel(A, n5, n6, error_approx_coarse14, residual_coarse14, flag, n_a[15]);
																	}


																	// prolongation
																	// residual_r
																	Real *error_approx_fine14 = new Real[n_a[14] + 1];
																	for (integer ii = 1; ii <= n_a[14]; ii++) {
																		error_approx_fine14[ii] = 0.0;
																	}

																	//for (integer ii = 1; ii <= n_a[15]; ii++) {// debug
																	//printf("error_approx_coarse14[%d]=%e\n",ii, error_approx_coarse14[ii]);

																	//printf("residual_coarse14[%d]=%e\n", ii, residual_coarse14[ii]);
																	//getchar();
																	//}
																	//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12]+2*nnz_a[13]+2*nnz_a[14]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12] +2*nnz_a[13]+2*nnz_a[14]+ nnz_a[15]; ii++) {// debug
																	//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
																	//if (ii % 20 == 0) getchar();
																	//}

																	integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13];
																	integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14];
																	prolongation(P, n7, n8, flag, error_approx_fine14, error_approx_coarse14, n_a[14], n_a[15]);

																	// correction
																	for (integer ii = 1; ii <= n_a[14]; ii++) {
																		error_approx_coarse13[ii] += error_approx_fine14[ii];
																	}

																	// free
																	delete[] error_approx_fine14;
																	delete[] error_approx_coarse14;
																	delete[] residual_coarse14;
																	delete[] residual_fine14;

																}


																// post smothing
																for (integer iter = 0; iter < nu2; iter++) {
																	integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13];
																	integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13] + nnz_a[14];
																	seidel(A, n5, n6, error_approx_coarse13, residual_coarse13, flag, n_a[14]);
																}


																// prolongation
																// residual_r
																Real *error_approx_fine13 = new Real[n_a[13] + 1];
																for (integer ii = 1; ii <= n_a[13]; ii++) {
																	error_approx_fine13[ii] = 0.0;
																}

																//for (integer ii = 1; ii <= n_a[14]; ii++) {// debug
																//printf("error_approx_coarse13[%d]=%e\n",ii, error_approx_coarse13[ii]);

																//printf("residual_coarse13[%d]=%e\n", ii, residual_coarse13[ii]);
																//getchar();
																//}
																//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12]+2*nnz_a[13]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12] +2*nnz_a[13]+ nnz_a[14]; ii++) {// debug
																//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
																//if (ii % 20 == 0) getchar();
																//}

																integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12];
																integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13];
																prolongation(P, n7, n8, flag, error_approx_fine13, error_approx_coarse13, n_a[13], n_a[14]);

																// correction
																for (integer ii = 1; ii <= n_a[13]; ii++) {
																	error_approx_coarse12[ii] += error_approx_fine13[ii];
																}

																// free
																delete[] error_approx_fine13;
																delete[] error_approx_coarse13;
																delete[] residual_coarse13;
																delete[] residual_fine13;

															}


															// post smothing
															for (integer iter = 0; iter < nu2; iter++) {
																integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12];
																integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + nnz_a[13];
																seidel(A, n5, n6, error_approx_coarse12, residual_coarse12, flag, n_a[13]);
															}


															// prolongation
															// residual_r
															Real *error_approx_fine12 = new Real[n_a[12] + 1];
															for (integer ii = 1; ii <= n_a[12]; ii++) {
																error_approx_fine12[ii] = 0.0;
															}

															//for (integer ii = 1; ii <= n_a[13]; ii++) {// debug
															//printf("error_approx_coarse12[%d]=%e\n",ii, error_approx_coarse12[ii]);

															//printf("residual_coarse12[%d]=%e\n", ii, residual_coarse12[ii]);
															//getchar();
															//}
															//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12]+ nnz_a[13]; ii++) {// debug
															//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
															//if (ii % 20 == 0) getchar();
															//}

															integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11];
															integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12];
															prolongation(P, n7, n8, flag, error_approx_fine12, error_approx_coarse12, n_a[12], n_a[13]);

															// correction
															for (integer ii = 1; ii <= n_a[12]; ii++) {
																error_approx_coarse11[ii] += error_approx_fine12[ii];
															}

															// free
															delete[] error_approx_fine12;
															delete[] error_approx_coarse12;
															delete[] residual_coarse12;
															delete[] residual_fine12;

														}



														// post smothing
														for (integer iter = 0; iter < nu2; iter++) {
															integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11];
															integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + nnz_a[12];
															seidel(A, n5, n6, error_approx_coarse11, residual_coarse11, flag, n_a[12]);
														}


														// prolongation
														// residual_r
														Real *error_approx_fine11 = new Real[n_a[11] + 1];
														for (integer ii = 1; ii <= n_a[11]; ii++) {
															error_approx_fine11[ii] = 0.0;
														}

														//for (integer ii = 1; ii <= n_a[12]; ii++) {// debug
														//printf("error_approx_coarse11[%d]=%e\n",ii, error_approx_coarse11[ii]);

														//printf("residual_coarse11[%d]=%e\n", ii, residual_coarse11[ii]);
														//getchar();
														//}
														//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+ nnz_a[12]; ii++) {// debug
														//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
														//if (ii % 20 == 0) getchar();
														//}

														integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10];
														integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11];
														prolongation(P, n7, n8, flag, error_approx_fine11, error_approx_coarse11, n_a[11], n_a[12]);

														// correction
														for (integer ii = 1; ii <= n_a[11]; ii++) {
															error_approx_coarse10[ii] += error_approx_fine11[ii];
														}

														// free
														delete[] error_approx_fine11;
														delete[] error_approx_coarse11;
														delete[] residual_coarse11;
														delete[] residual_fine11;

													}


													// post smothing
													for (integer iter = 0; iter < nu2; iter++) {
														integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9]+2*nnz_a[10];
														integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] +2*nnz_a[10]+ nnz_a[11];
														seidel(A, n5, n6, error_approx_coarse10, residual_coarse10, flag, n_a[11]);
													}


													// prolongation
													// residual_r
													Real *error_approx_fine10 = new Real[n_a[10] + 1];
													for (integer ii = 1; ii <= n_a[10]; ii++) {
														error_approx_fine10[ii] = 0.0;
													}

													//for (integer ii = 1; ii <= n_a[11]; ii++) {// debug
													//printf("error_approx_coarse10[%d]=%e\n",ii, error_approx_coarse10[ii]);

													//printf("residual_coarse10[%d]=%e\n", ii, residual_coarse10[ii]);
													//getchar();
													//}
													//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+ nnz_a[11]; ii++) {// debug
													//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
													//if (ii % 20 == 0) getchar();
													//}

													integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8]+nnz_aRP[9];
													integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9]+nnz_aRP[10];
													prolongation(P, n7, n8, flag, error_approx_fine10, error_approx_coarse10, n_a[10], n_a[11]);

													// correction
													for (integer ii = 1; ii <= n_a[10]; ii++) {
														error_approx_coarse9[ii] += error_approx_fine10[ii];
													}

													// free
													delete[] error_approx_fine10;
													delete[] error_approx_coarse10;
													delete[] residual_coarse10;
													delete[] residual_fine10;

												}



												// post smothing
												for (integer iter = 0; iter < nu2; iter++) {
													integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8]+2*nnz_a[9];
													integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] +2*nnz_a[9]+ nnz_a[10];
													seidel(A, n5, n6, error_approx_coarse9, residual_coarse9, flag, n_a[10]);
												}


												// prolongation
												// residual_r
												Real *error_approx_fine9 = new Real[n_a[9] + 1];
												for (integer ii = 1; ii <= n_a[9]; ii++) {
													error_approx_fine9[ii] = 0.0;
												}

												//for (integer ii = 1; ii <= n_a[10]; ii++) {// debug
												//printf("error_approx_coarse9[%d]=%e\n",ii, error_approx_coarse9[ii]);

												//printf("residual_coarse9[%d]=%e\n", ii, residual_coarse9[ii]);
												//getchar();
												//}
												//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+ nnz_a[10]; ii++) {// debug
												//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
												//if (ii % 20 == 0) getchar();
												//}

												integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8];
												integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8]+nnz_aRP[9];
												prolongation(P, n7, n8, flag, error_approx_fine9, error_approx_coarse9, n_a[9], n_a[10]);

												// correction
												for (integer ii = 1; ii <= n_a[9]; ii++) {
													error_approx_coarse8[ii] += error_approx_fine9[ii];
												}

												// free
												delete[] error_approx_fine9;
												delete[] error_approx_coarse9;
												delete[] residual_coarse9;
												delete[] residual_fine9;

											}

											// post smothing
											for (integer iter = 0; iter < nu2; iter++) {
												integer n5 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8];
												integer n6 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + nnz_a[9];
												seidel(A, n5, n6, error_approx_coarse8, residual_coarse8, flag, n_a[9]);
											}


											// prolongation
											// residual_r
											Real *error_approx_fine8 = new Real[n_a[8] + 1];
											for (integer ii = 1; ii <= n_a[8]; ii++) {
												error_approx_fine8[ii] = 0.0;
											}

											//for (integer ii = 1; ii <= n_a[9]; ii++) {// debug
											//printf("error_approx_coarse8[%d]=%e\n",ii, error_approx_coarse8[ii]);

											//printf("residual_coarse8[%d]=%e\n", ii, residual_coarse8[ii]);
											//getchar();
											//}
											//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+ nnz_a[9]; ii++) {// debug
											//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
											//if (ii % 20 == 0) getchar();
											//}

											integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7];
											integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8];
											prolongation(P,n7, n8, flag, error_approx_fine8, error_approx_coarse8, n_a[8], n_a[9]);

											// correction
											for (integer ii = 1; ii <= n_a[8]; ii++) {
												error_approx_coarse7[ii] += error_approx_fine8[ii];
											}

											// free
											delete[] error_approx_fine8;
											delete[] error_approx_coarse8;
											delete[] residual_coarse8;
											delete[] residual_fine8;

										}

										// post smothing
										for (integer iter = 0; iter < nu2; iter++) {
											seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + nnz_a[8], error_approx_coarse7, residual_coarse7, flag, n_a[8]);
										}


										// prolongation
										// residual_r
										Real *error_approx_fine7 = new Real[n_a[7] + 1];
										for (integer ii = 1; ii <= n_a[7]; ii++) {
											error_approx_fine7[ii] = 0.0;
										}

										//for (integer ii = 1; ii <= n_a[8]; ii++) {// debug
										//printf("error_approx_coarse7[%d]=%e\n",ii, error_approx_coarse7[ii]);

										//printf("residual_coarse7[%d]=%e\n", ii, residual_coarse7[ii]);
										//getchar();
										//}
										//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+ nnz_a[8]; ii++) {// debug
										//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
										//if (ii % 20 == 0) getchar();
										//}

										prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6]+nnz_aRP[7], flag, error_approx_fine7, error_approx_coarse7, n_a[7], n_a[8]);

										// correction
										for (integer ii = 1; ii <= n_a[7]; ii++) {
											error_approx_coarse6[ii] += error_approx_fine7[ii];
										}

										// free
										delete[] error_approx_fine7;
										delete[] error_approx_coarse7;
										delete[] residual_coarse7;
										delete[] residual_fine7;

									}


									// post smothing
									for (integer iter = 0; iter < nu2; iter++) {
										seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + nnz_a[7], error_approx_coarse6, residual_coarse6, flag, n_a[7]);
									}


									// prolongation
									// residual_r
									//Real *error_approx_fine6 = new Real[n_a[6] + 1];
									for (integer ii = 1; ii <= n_a[6]; ii++) {
										error_approx_fine6[ii] = 0.0;
									}

									//for (integer ii = 1; ii <= n_a[7]; ii++) {// debug
									//printf("error_approx_coarse6[%d]=%e\n",ii, error_approx_coarse6[ii]);

									//printf("residual_coarse6[%d]=%e\n", ii, residual_coarse6[ii]);
									//getchar();
									//}
									//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+ nnz_a[7]; ii++) {// debug
									//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
									//if (ii % 20 == 0) getchar();
									//}

									prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5]+nnz_aRP[6], flag, error_approx_fine6, error_approx_coarse6, n_a[6], n_a[7]);

									// correction
									for (integer ii = 1; ii <= n_a[6]; ii++) {
										error_approx_coarse5[ii] += error_approx_fine6[ii];
									}

									// free
									//delete[] error_approx_fine6;
									//delete[] error_approx_coarse6;
									//delete[] residual_coarse6;
									//delete[] residual_fine6;

								}

								// post smothing
								for (integer iter = 0; iter < nu2; iter++) {
									//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + nnz_a[6], error_approx_coarse5, residual_coarse5, flag, n_a[6]);
									seidelq(A, 1, n_a[6], error_approx_coarse5, residual_coarse5, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5]);
								}

								// prolongation
								// residual_r
								//Real *error_approx_fine5 = new Real[n_a[5] + 1];
								for (integer ii = 1; ii <= n_a[5]; ii++) {
									error_approx_fine5[ii] = 0.0;
								}

								//for (integer ii = 1; ii <= n_a[6]; ii++) {// debug
								//printf("error_approx_coarse5[%d]=%e\n",ii, error_approx_coarse5[ii]);

								//printf("residual_coarse5[%d]=%e\n", ii, residual_coarse5[ii]);
								//getchar();
								//}
								//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+ nnz_a[6]; ii++) {// debug
								//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
								//if (ii % 20 == 0) getchar();
								//}

								prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], flag, error_approx_fine5, error_approx_coarse5, n_a[5], n_a[6]);

								// correction
								for (integer ii = 1; ii <= n_a[5]; ii++) {
									error_approx_coarse4[ii] += error_approx_fine5[ii];
								}

								// free
								//delete[] error_approx_fine5;
								//delete[] error_approx_coarse5;
								//delete[] residual_coarse5;
								//delete[] residual_fine5;

							}
							// post smothing
							for (integer iter = 0; iter < nu2; iter++) {
								//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + nnz_a[5], error_approx_coarse4, residual_coarse4, flag, n_a[5]);
								seidelq(A, 1, n_a[5], error_approx_coarse4, residual_coarse4, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4]);
							}

							// prolongation
							// residual_r
							//Real *error_approx_fine4 = new Real[n_a[4] + 1];
							for (integer ii = 1; ii <= n_a[4]; ii++) {
								error_approx_fine4[ii] = 0.0;
							}

							//for (integer ii = 1; ii <= n_a[5]; ii++) {// debug
							//printf("error_approx_coarse4[%d]=%e\n",ii, error_approx_coarse4[ii]);

							//printf("residual_coarse4[%d]=%e\n", ii, residual_coarse4[ii]);
							//getchar();
							//}
							//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ nnz_a[5]; ii++) {// debug
							//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
							//if (ii % 20 == 0) getchar();
							//}

							prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4], flag, error_approx_fine4, error_approx_coarse4, n_a[4], n_a[5]);

							// correction
							for (integer ii = 1; ii <= n_a[4]; ii++) {
								error_approx_coarse3[ii] += error_approx_fine4[ii];
							}

							// free
							//delete[] error_approx_fine4;
							//delete[] error_approx_coarse4;
							//delete[] residual_coarse4;
							//delete[] residual_fine4;

						}
						// post smothing
						for (integer iter = 0; iter < nu2; iter++) {
							//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + nnz_a[4], error_approx_coarse3, residual_coarse3, flag, n_a[4]);
							seidelq(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3]);
						}

						// prolongation
						// residual_r
						//Real *error_approx_fine3 = new Real[n_a[3] + 1];
						for (integer ii = 1; ii <= n_a[3]; ii++) {
							error_approx_fine3[ii] = 0.0;
						}

						//for (integer ii = 1; ii <= n_a[4]; ii++) {// debug
						//printf("error_approx_coarse3[%d]=%e\n",ii, error_approx_coarse3[ii]);

						//printf("residual_coarse3[%d]=%e\n", ii, residual_coarse3[ii]);
						//getchar();
						//}
						//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ nnz_a[4]; ii++) {// deug
						//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
						//if (ii % 20 == 0) getchar();
						//}

						prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], flag, error_approx_fine3, error_approx_coarse3, n_a[3], n_a[4]);

						// correction
						for (integer ii = 1; ii <= n_a[3]; ii++) {
							error_approx_coarse2[ii] += error_approx_fine3[ii];
						}

						// free
						//delete[] error_approx_fine3;
						//delete[] error_approx_coarse3;
						//delete[] residual_coarse3;
						//delete[] residual_fine3;

					}
					// post smothing
					for (integer iter = 0; iter < nu2; iter++) {
						//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + nnz_a[3], error_approx_coarse2, residual_coarse2, flag, n_a[3]);
						seidelq(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2]);

					}

					// prolongation
					// residual_r
					//Real *error_approx_fine2 = new Real[n_a[2] + 1];
					for (integer ii = 1; ii <= n_a[2]; ii++) {
						error_approx_fine2[ii] = 0.0;
					}

					//for (integer ii = 1; ii <= n_a[3]; ii++) {// deug
					//printf("error_approx_coarse2[%d]=%e\n",ii, error_approx_coarse2[ii]);

					//printf("residual_coarse2[%d]=%e\n", ii, residual_coarse2[ii]);
					//getchar();
					//}
					//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ nnz_a[3]; ii++) {// deug
					//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
					//if (ii % 20 == 0) getchar();
					//}

					prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], flag, error_approx_fine2, error_approx_coarse2, n_a[2], n_a[3]);

					// correction
					for (integer ii = 1; ii <= n_a[2]; ii++) {
						error_approx_coarse1[ii] += error_approx_fine2[ii];
					}

					// free
					//delete[] error_approx_fine2;
					//delete[] error_approx_coarse2;
					//delete[] residual_coarse2;
					//delete[] residual_fine2;

				}
				// post smothing
				for (integer iter = 0; iter < nu2; iter++) {
					//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1], 2 * nnz_a[0] + 2 * nnz_a[1] + nnz_a[2], error_approx_coarse1, residual_coarse1, flag, n_a[2]);
					seidelq(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0] + n_a[1]);
				}

				// prolongation
				// residual_r
				//Real *error_approx_fine1 = new Real[n_a[1] + 1];
				for (integer ii = 1; ii <= n_a[1]; ii++) {
					error_approx_fine1[ii] = 0.0;
				}

				//for (integer ii = 1; ii <= n_a[2]; ii++) {// deug
					//printf("error_approx_coarse1[%d]=%e\n",ii, error_approx_coarse1[ii]);

					//printf("residual_coarse1[%d]=%e\n", ii, residual_coarse1[ii]);
					//getchar();
				//}
				//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + nnz_a[2]; ii++) {// deug
					//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
					//if (ii % 20 == 0) getchar();
				//}

				prolongation(P, 1+nnz_aRP[0], nnz_aRP[0]+nnz_aRP[1], flag, error_approx_fine1, error_approx_coarse1, n_a[1], n_a[2]);

				// correction
				for (integer ii = 1; ii <= n_a[1]; ii++) {
					error_approx_coarse[ii] += error_approx_fine1[ii];
				}

				// free
				//delete[] error_approx_fine1;
				//delete[] error_approx_coarse1;
				//delete[] residual_coarse1;
				//delete[] residual_fine1;

			}
			
			// post smothing
			for (integer iter = 0; iter < nu2; iter++) {
				//seidel(A, 1 + 2 * nnz_a[0], 2 * nnz_a[0] + nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1]);
				seidelq(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0]);
			}

			// prolongation
			// residual_r
			//Real *error_approx_fine = new Real[n_a[0] + 1];
			for (integer ii = 1; ii <= n_a[0]; ii++) {
				error_approx_fine[ii] = 0.0;
			}

			prolongation(P, 1, nnz_aRP[0], flag, error_approx_fine, error_approx_coarse, n_a[0], n_a[1]);

			// correction
			for (integer ii = 1; ii <= n_a[0]; ii++) {
				x[ii] += error_approx_fine[ii];
			}

			// free
			//delete[] error_approx_fine;
			//delete[] error_approx_coarse;
			//delete[] residual_coarse;
			//delete[] residual_fine;
		}
		// post smother
		for (integer iter = 0; iter < nu2; iter++) {
			//seidel(A, 1, nnz_a[0], x, b, flag, n_a[0]);
			//quick seidel
			seidelq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);
		}
	}

	system("pause");

	// free
	delete[] error_approx_fine;
	if (ilevel > 1) {
		delete[] error_approx_coarse;
		delete[] residual_coarse;
		if (ilevel > 2) {
			// free
			delete[] error_approx_fine1;
			delete[] error_approx_coarse1;
			delete[] residual_coarse1;
			delete[] residual_fine1;
			if (ilevel > 3) {
				// free
				delete[] error_approx_fine2;
				delete[] error_approx_coarse2;
				delete[] residual_coarse2;
				delete[] residual_fine2;
				if (ilevel > 4) {
					// free
					delete[] error_approx_fine3;
					delete[] error_approx_coarse3;
					delete[] residual_coarse3;
					delete[] residual_fine3;
					if (ilevel > 5) {
						// free
						delete[] error_approx_fine4;
						delete[] error_approx_coarse4;
						delete[] residual_coarse4;
						delete[] residual_fine4;
						if (ilevel > 6) {
							// free
							delete[] error_approx_fine5;
							delete[] error_approx_coarse5;
							delete[] residual_coarse5;
							delete[] residual_fine5;
							if (ilevel > 7) {
								// free
								delete[] error_approx_fine6;
								delete[] error_approx_coarse6;
								delete[] residual_coarse6;
								delete[] residual_fine6;
							}
						}
					}
				}
			}
		}
	}
	delete[] residual_fine;

	delete[] row_ptr_start;
	delete[] row_ptr_end;
	return 0;

} // aggregative_amg

// возвращает максимум из двух целых чисел.
integer i_my_max(integer ia, integer ib) {
	if (ia > ib) {
		return ia;
	}
	else {
		return ib;
	}
}

// возвращает минимум из двух целых чисел.
integer i_my_min(integer ia, integer ib) {
	if (ia > ib) {
		return ib;
	}
	else {
		return ia;
	}
}

// 18 октября 2015. Полностью работоспособный мультигрид.
// Тестировалось на условиях Дирихле но должно работать на любых 
// краевых задачах. 18 октября датируется версия 0.04. Версия 0.04 на треть
// быстрее версии 0_03. Были ускорены как операции построения C-F разбиения, 
// так и нахождение оператора Галёркина. При нахождении С-F разбиения 
// учитывается уже построеннная его часть и поэтому число сканирований на
// на поздних циклах сокращается охватывая только не построенную часть.
// При нахождении произведеия Галёркина получена самая оптимальная по быстродествию версия,
// Основанная на алгоритме слияния отсортированных списков.
// 4 октября правильное построение последовательности вложенных графов.
// 30 сентября продолжаем исправление метода. Делаем классический 
// алгебраический многосеточный метод на основе  C-F разбиения.
// 16 сентября 2015 года обнаружено что операции 
// сгрубления и интерполяции сделаны совершенно неверно,
// и если сгрубление еще в какой-то мере проецирует то интерполяция просто никакая.
// Операции сгрубления и интерполляции будут сделаны заново на основе статьи К.Н. Волкова в новой версии солвера.
// 3 september 2015 Villa Borgese.
integer classic_aglomerative_amg1(Ak* &A, integer nnz,
	integer n, // dimension of vectors x and b.
	Ak* &R, // restriction
	Ak* &P, // prolongation
	Ak* &Atemp,
	Ak* &Atemp2,
	Real* &x, Real* &b,
	Real theta
	) {

	

	// контроль числа сильных связей между переменными.
	// Real theta = 0.25;  // 0.25 for 2D and 0.5 for 3D 

	bool blite_and_quick_version = false;
	// false даёт более качественное построение иерархии матриц, 
	// с ним количество итераций сокращается.
	const integer heapsortsizelimit = 10000000; // 10M

	const Real RealZERO = 1.0e-300;// 1.0e-10;
	const Real divisionZERO = 1.0e-300;
	const Real RealMAXIMUM = 1.0e300;
	// включение/отключение отладочного режима.
	bool debug_reshime = false;

	const integer max_sosed = 27850;
	// Мы будем поддерживать информацию о максимальном количестве соседей.
	integer Maximumsosedcount = -1;
	// мы получим порядковый номер элемента с максимальным число соседей в матрице А,
	// прямо в результате построения C-F разбиения.
	integer icandidate_shadow = -1;
	bool bmaxsosedinfoactive = false;

	// количество рекурсивных вызовов ограничено, поэтому QuickSort не подходит.
	bool bquicktypesort = false;


	// x_coarse[1..n_coarse]=R n_coarse x n_fine*x_fine[1..n_fine];
	// x_fine[1..n_fine]=P n_fine x n_coarse*x_coarse[1..n_coarse];

	// нумерация начинается с единицы.

	const integer maxlevel = 30;
	integer ilevel = 1;
	integer nnz_a[maxlevel];
	integer n_a[maxlevel];
	nnz_a[0] = nnz;
	n_a[0] = n;
	bool* flag = new bool[n + 1];
	//bool* flag_ = new bool[n + 1];
	bool* flag_shadow = new bool[n + 1];
	integer iadd = 0;
	integer nnzR = 1;
	integer iaddR = 0;
	integer nnz_aRP[maxlevel];
	bool bcontinue = true;
	bool* this_is_C_node = new bool[n + 1];
	bool* this_is_F_node = new bool[n + 1];

	// храним соседей новых F узлов
	const integer imaxpoolsize = 2000;
	integer pool_sosed[imaxpoolsize];
	integer imax_pool_ind;

	while ((ilevel < maxlevel - 1) && (n_a[ilevel - 1] > 50) && (bcontinue)) {

		Maximumsosedcount = -1;
		bmaxsosedinfoactive = false;

		if (ilevel > 1) {
			if (n_a[ilevel - 2] == n_a[ilevel - 1]) break;
		}

		if (((ilevel > 1) && (nnz_a[ilevel - 1] > nnz_a[ilevel - 2]))) {
			//break;
		}

		// level info.
		if (ilevel == 2) {
			printf("ilevel=%d",ilevel);
			printf("n_a[0]=%d n_a[1]=%d nnz_a[0]=%d nnz_a[1]=%d iadd=%d\n", n_a[0], n_a[1],  nnz_a[0], nnz_a[1],  iadd);
			if (debug_reshime) system("pause");
		}
		if (ilevel == 3) {
			printf("ilevel=%d", ilevel);
			printf("n_a[0]=%d n_a[1]=%d n_a[2]=%d nnz_a[0]=%d nnz_a[1]=%d nnz_a[2]=%d iadd=%d\n", n_a[0], n_a[1], n_a[2], nnz_a[0], nnz_a[1], nnz_a[2], iadd);
			if (debug_reshime) system("pause");
		}
		if (ilevel == 4) {
			printf("ilevel=%d", ilevel);
			printf("n_a[0]=%d n_a[1]=%d n_a[2]=%d n_a[3]=%d \n", n_a[0], n_a[1], n_a[2], n_a[3]);
			printf("nnz_a[0]=%d nnz_a[1]=%d nnz_a[2]=%d nnz_a[3]=%d iadd=%d\n",  nnz_a[0], nnz_a[1], nnz_a[2], nnz_a[3], iadd);
			if (debug_reshime) system("pause");
		}
		
		if ((ilevel == 5) || (ilevel == 6) || (ilevel == 7) || (ilevel == 8) || (ilevel == 9) || (ilevel == 10)) {
			printf("ilevel=%d\n", ilevel);
			for (integer i_1 = 0; i_1 < ilevel; i_1++) {
				printf("n_a[%d]=%d nnz_a[%d]=%d\n", i_1, n_a[i_1], i_1, nnz_a[i_1]);
			}
			if (debug_reshime) system("pause");
		}

		if ((ilevel == 11) || (ilevel == 12) || (ilevel == 13) || (ilevel == 14) || (ilevel == 15) ) {
			printf("ilevel=%d\n", ilevel);
			for (integer i_1 = 0; i_1 < ilevel; i_1++) {
				printf("n_a[%d]=%d nnz_a[%d]=%d\n", i_1, n_a[i_1], i_1, nnz_a[i_1]);
			}
			if (debug_reshime) system("pause");
		}

		if ((ilevel == 16) || (ilevel == 17) || (ilevel == 18) || (ilevel == 19) || (ilevel == 20)) {
			printf("ilevel=%d\n", ilevel);
			for (integer i_1 = 0; i_1 < ilevel; i_1++) {
				printf("n_a[%d]=%d nnz_a[%d]=%d\n", i_1, n_a[i_1], i_1, nnz_a[i_1]);
			}
			if (debug_reshime) system("pause");
		}

		if (ilevel > 17) {
			printf("very big matrix (mesh). no programming.\n");
			system("pause");
			exit(1);
		}

		//nnzR = 1;

		for (int ii = 1; ii <= n_a[ilevel-1]; ii++) {
			this_is_C_node[ii] = false;
			this_is_F_node[ii] = false;
		}

		for (int ii = n_a[ilevel - 1]; ii <= n; ii++) {
			this_is_C_node[ii] = false;
			this_is_F_node[ii] = false;
		}

		// сортировка исходной  А  по i.
		//heapsort(A, key=i*n_a[ilevel - 1] + j, 1, nnz_a[ilevel - 1]);
		if (bquicktypesort) {
			QuickSort(A, /*n_a[ilevel - 1],*/ 1 + iadd, nnz_a[ilevel - 1] + iadd);
		}
		else {
			if (nnz_a[ilevel - 1] < heapsortsizelimit) {
				HeapSort(A, n_a[ilevel - 1], 1 + iadd, nnz_a[ilevel - 1] + iadd);
			}
			else {
				Ak* Aorig = &A[1 + iadd];
				MergeSort(Aorig, nnz_a[ilevel - 1]);
				Aorig = NULL;
			}
		}

		// Экономично вычисляет полуширину ленты матрицы.
		integer band_size = -1;
		integer band_size_i = -1;
		for (integer i_1 = 1; i_1 <= n; i_1++) {
			flag[i_1] = false;
		}
		for (integer ii = 1 + iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
			if (!flag[A[ii].i]) {
				flag[A[ii].i] = true;
				integer istart = ii;
				while ((istart <= nnz_a[ilevel -1]  +iadd) && (A[istart].i == A[ii].i) && (A[istart].j != A[ii].i)) istart++;
				if ((istart <= nnz_a[ilevel - 1] + iadd) && (A[istart].i == A[ii].i)) {
					integer ifound = istart;
					istart = ii;
					while ((istart <= nnz_a[ilevel - 1]  + iadd) && (A[istart].i == A[ii].i)) {
						
						if (A[istart].j != A[ii].i) {
							integer ij = BinarySearchAi(A, A[istart].j, 1 + iadd, nnz_a[ilevel - 1] + iadd);
							while ((ij <= nnz_a[ilevel - 1] + iadd) &&(A[ij].i == A[istart].j) && (A[ij].j != A[istart].j)) ij++;
							if ((ij <= nnz_a[ilevel -1] + iadd) && (A[ij].i == A[istart].j)&&(A[ij].j == A[istart].j)) {
								if (abs(ij - ifound) > band_size) band_size = abs(ij - ifound);
								if (abs(A[ij].j - A[ifound].i) > band_size_i) band_size_i = abs(A[ij].j - A[ifound].i);
							}
						}
						istart++;
					}
				}
			}
		}
		printf("bandsize=%d\n",band_size);
		//band_size = -1; // OFF band_size acselerator.
		
		// Создаём копию в Atemp, копия будет отсортирована по j.
		for (int i_1 = 1 + iadd; i_1 <= nnz_a[ilevel - 1] + iadd; i_1++) {
			Atemp[i_1 - iadd] = A[i_1];
			Atemp[i_1 - iadd].ind = i_1; // запоминаем первоначальную позицию в А.
		}
		// Сортируем копию по j.
		// Мы сортируем по j, чтобы потом быстро искать по j.
		HeapSort_j(Atemp, 1, nnz_a[ilevel - 1]);

		// Экономично вычисляет полуширину ленты матрицы.
		integer band_sizej = -1;
		for (integer i_1 = 1; i_1 <= n; i_1++) {
			flag[i_1] = false;
		}
		for (integer ii = 1 ; ii <= nnz_a[ilevel - 1] ; ii++) {
			if (!flag[Atemp[ii].j]) {
				flag[Atemp[ii].j] = true;
				integer istart = ii;
				while ((istart <= nnz_a[ilevel - 1] ) && (Atemp[istart].j == Atemp[ii].j) && (Atemp[istart].i != Atemp[ii].j)) istart++;
				if ((istart <= nnz_a[ilevel - 1] ) && (Atemp[istart].j == Atemp[ii].j)) {
					integer ifound = istart;
					istart = ii;
					while ((istart <= nnz_a[ilevel - 1] ) && (Atemp[istart].j == Atemp[ii].j)) {

						if (Atemp[istart].i != Atemp[ii].j) {
							integer ij = BinarySearchAj(Atemp, Atemp[istart].i, 1 , nnz_a[ilevel - 1] );
							while ((ij <= nnz_a[ilevel - 1] ) && (Atemp[ij].j == Atemp[istart].i) && (Atemp[ij].i != Atemp[istart].i)) ij++;
							if ((ij <= nnz_a[ilevel - 1] ) && (Atemp[ij].j == Atemp[istart].i) && (Atemp[ij].i == Atemp[istart].i)) {
								if (abs(ij - ifound) > band_sizej) band_sizej = abs(ij - ifound);
							}
						}
						istart++;
					}
				}
			}
		}
		printf("bandsizej=%d\n", band_sizej);
		if ((band_size == -1) && (band_sizej == -1)) {
			bcontinue = false;
			// Судя по всему иерархия матриц уже построена.
			// Досрочный выход из этпа конструирования по причине
			// отсутствия связей в матрице... ?
			// Возможно эта неопределённость была причиной сбоя на реальном тестировании.
			break;
		}
		//getchar();
		//band_size = -1; // OFF band_size acselerator.
		

		//if (ilevel == 2) {
			// Матрица второго уровня составлена совершенно неверно.
			// Возможно неверно было составлено произведение Галёркина RAP.
			//for (integer ii7 = 1 + iadd; ii7 <= iadd + nnz_a[ilevel - 1]; ii7++) {
				//if (A[ii7].i > n_a[ilevel - 1]) {
					//		printf("matrix incorrect i\n");
				//}
				//if (A[ii7].j > n_a[ilevel - 1]) {
					//	printf("%d ",A[ii7].j);
				//}
			//	printf("A[%d].i=%d A[%d].j=%d A[%d].aij=%e\n", ii7, A[ii7].i, ii7, A[ii7].j, ii7, A[ii7].aij);
				//system("pause");
			//}
		//}

		//for (integer ii = 1 + iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
			//	printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n",ii,A[ii].aij,ii,A[ii].i,ii,A[ii].j);
			//if (ii % 20 == 0) getchar();
		//}

		for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
			flag[ii] = false;
		}

		// позиция начала каждой строки в матрице.
		integer* row_startA = new integer[n_a[ilevel - 1] + 1];
		for (integer ii = 1 + iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
			if (flag[A[ii].i] == false) {
				flag[A[ii].i] = true;
				row_startA[A[ii].i] = ii;
			}
		}

		for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
			flag[ii] = false;
		}
		

		// вычисляем для кадого узла число его соседей.
		integer* count_sosed = new integer[n_a[ilevel - 1]+1];
		for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
			count_sosed[ii] = 0; // нет соседей.
		}

		if (blite_and_quick_version)
		{
			// А при таком определении узел Дирихле имеет ноль соседей.
			// Соседей вычисляем на самой первой матрице А (самой левой).
			for (integer ii = 1 + iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
				if (flag[A[ii].i] == false) {
					integer ic = -1;
					for (integer is0 = ii; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == A[ii].i); is0++) {
						ic++; //i,j
					}
					count_sosed[A[ii].i] = ic;
					if (ic > Maximumsosedcount) {
						Maximumsosedcount = ic;
						icandidate_shadow = ii;
						bmaxsosedinfoactive = true;
					}
					flag[A[ii].i] = true;
				}
			}
		}
		else {


			// При таком коде узел Дирихле тоже имеет соседа, сосед это 
			// внутренний узел который связан с этим узлом Дирихле.
			// Соседей вычисляем на самой первой матрице А (самой левой).
			for (integer ii = 1 + iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
				if (flag[A[ii].i] == false) {
					integer ic = -1;
					integer cand[max_sosed];
					for (integer is0 = ii; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == A[ii].i); is0++) {
						ic++; //i,j
						cand[ic] = A[is0].j;
					}
					integer len_sosed = ic;
					// Найти столбец j который равен индексу A[ii].i
					//for (integer ii1 = 1 + iadd; ii1 <= nnz_a[ilevel - 1] + iadd; ii1++) {
						//if (A[ii1].i != A[ii].i) {
						//	if (A[ii1].j == A[ii].i) {
								// j,i
						//		bool foundsosed = false;
						//		for (integer i_1 = 0; i_1 <= len_sosed; i_1++) {
						//			if (A[ii1].j == cand[i_1]) foundsosed = true;
						//		}
						//		if (!foundsosed) {
						//			ic++;
						//			cand[ic] = A[ii1].j;
						//			len_sosed++;
						//		}
						//	}
						//}
					//}
					// Ускоренная версия с бинарным поиском по j.
					integer ii2 = BinarySearchAj(Atemp, A[ii].i, 1, nnz_a[ilevel - 1]);
					for (integer ii1 = ii2; (ii1 <= nnz_a[ilevel - 1])&&(Atemp[ii1].j==A[ii].i); ii1++) {
						if (Atemp[ii1].i != A[ii].i) {
							// j,i
							bool foundsosed = false;
							for (integer i_1 = 0; i_1 <= len_sosed; i_1++) {
								if (Atemp[ii1].j == cand[i_1]) foundsosed = true;
							}
							if (!foundsosed) {
								ic++;
								cand[ic] = Atemp[ii1].j;
								len_sosed++;
							}
						}
					}
					

					count_sosed[A[ii].i] = ic;
					if (ic > Maximumsosedcount) {
						Maximumsosedcount = ic;
						icandidate_shadow = ii;
						bmaxsosedinfoactive = true;
					}
					flag[A[ii].i] = true;
				}
			}
		}
	

		integer maxsosed = 0;
		integer icandidate = 0;
		for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
			flag[ii] = false; // init flag
		}
		// Находим узел с наибольшим числом соседей и запоминаем его.
		// Это первый встретившийся узел с наибольшим числом соседей.
		for (integer ii = 1 + iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
			if (flag[A[ii].i] == false) {
				if (count_sosed[A[ii].i] > maxsosed) {
					maxsosed = count_sosed[A[ii].i];
					icandidate = ii;
					if (bmaxsosedinfoactive) {
						// организуем досрочный выход из цикла for.
						// Это должно сильно сокращать количество сканирований.
						if (maxsosed == Maximumsosedcount) break;
					}
				}
				flag[A[ii].i] = true;
			}
		}

		for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
			flag[ii] = false; // init flag
		}


		// нужно выделить кто попал в coarse, а кто в этот раз попал в Fine Выделить всех кто соседствует
		// с новыми Fine увеличить им счётчик соседей.

		integer n_coarce = 1; // начальный номер C узла.
		nnzR = 1;
		
		const integer NULL_SOSED = -1;
		integer vacant = NULL_SOSED;
		bool bcontinue = true;

		// Построение C-F разбиения.
		//while (icandidate != 0)
		integer icountprohod = 0;
		// Мы будем заоминать с какой позиции начинаются ещё не помеченные узлы,
		// это сократит количество перебираемых элементов в поиске узла с максимальным 
		// количеством соседей.
		integer ibegining_start_index_found_maximum = 1 + iadd;
		// храним те узлы которые уже были пройдены при конструировании.
		bool *bmarkervisit = new bool[n + 1];
		// поначалу все узлы помечены как непосещенные.
		for (integer i_1 = 1; i_1 <= n; i_1++) bmarkervisit[i_1] = false;

		for (integer i_1 = 1; i_1 <= nnz_a[ilevel - 1]; i_1++) {
			Atemp2[i_1] = A[i_1 + iadd]; // copy
			Atemp2[i_1].ind = i_1 + iadd;
			A[i_1 + iadd].ind = i_1; // запоминаем обратную связь.
			//if (ilevel == 2) {
				//printf("%e %d %d %d\n", Atemp2[i_1].aij, Atemp2[i_1].i, Atemp2[i_1].j, Atemp2[i_1].ind);
				//getchar();
			//}
		}

		// увеличение быстродействия достигается 
		// сокращением пределов сканирования
		// здесь хранится индекс начала сканирования flag.
		integer istartflag_scan = 1;

		//if (bmaxsosedinfoactive == true) {
			//printf("bmaxsosedinfoactive==true\n");
			//getchar();
		//}
		List* plist = NULL;
		List* plist_current = NULL;

while (bcontinue)
{



	imax_pool_ind = 1;

	//printf("prohod count number %d\n", icountprohod);


	integer set[max_sosed]; // не более 20 узлов в одном агрегате.
	// инициализация убрана потому что она не нужна и она сильно тормозит быстродействие.
	//for (integer js = 0; js < max_sosed; js++) {
	//set[js] = NULL_SOSED;
	//}
	integer ic = 0;

	integer ii = icandidate;
	if (flag[A[ii].i] == false) {
		// Вычисляем по немодифицированной матрице А (хранящейся слева).
		integer nnzRl = nnzR + iaddR;


		ic=0; // Обязательная инициализация.
		set[ic] = A[ii].i;
		
		
		Real max_vnediagonal = -1.0; // максимальное значение модуля вне диагонального элемента. 
		// добавляем диагональный элемент.
		// узел set[0]==A[is0].i.
		// Нахождение значения максимального внедиагольного элемента, с 
		// учётом того что даже узел Дирихле связан с одним внутренним узлом расчётной области.
		this_is_C_node[set[0]] = true;
		bmarkervisit[set[0]] = true;
		pool_sosed[0] = ii;
		for (integer is0 = ii; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == set[0]); is0++) {
			if (A[is0].j == set[0]) {

				nnzRl++;
				break;
			}
			else {
				if (fabs(A[is0].aij) > max_vnediagonal) {
					max_vnediagonal = fabs(A[is0].aij); //i,j
				}
			}
			if (!blite_and_quick_version) {
				// Этот цикл является добавочным.
				// Найти связь j,i где кандидат i еще не был включен в строящееся С-F разбиение.
				// Медленный линейный поиск.
				//for (integer ii1 = 1 + iadd; ii1 <= nnz_a[ilevel - 1] + iadd; ii1++) {
					//if (A[ii1].i != set[0]) {
						//if (!flag[A[ii1].i]) {
						//	if (A[ii1].j == set[0]) {
						//		if (fabs(A[ii1].aij) > max_vnediagonal) {
						//			max_vnediagonal = fabs(A[ii1].aij); //j,i
						//		}
						//	}
						//}
					//}
				//}

				// Этот цикл является добавочным.
				// Найти связь j,i где кандидат i еще не был включен в строящееся С-F разбиение.
				// Ускоренная версия на основе двоичного поиска.
				integer ii2 = BinarySearchAj(Atemp, set[0], 1, nnz_a[ilevel - 1]);
				for (integer ii1 = ii2; (ii1 <= nnz_a[ilevel - 1]) && (Atemp[ii1].j == set[0]); ii1++) 
				{
					if (Atemp[ii1].i != set[0]) {
						if (!flag[Atemp[ii1].i]) {
								if (fabs(Atemp[ii1].aij) > max_vnediagonal) {
									max_vnediagonal = fabs(Atemp[ii1].aij); //j,i
								}
						}
					}
				}
				
			}
		}

		ic++;



		// если узел j ещё не был добавлен в агрегат.
		if (flag[A[ii].j] == false) {
			if ((A[ii].j != set[0]) && (fabs(A[ii].aij) >= theta*max_vnediagonal)) {
				vacant = A[ii].j;
				pool_sosed[imax_pool_ind] = ii;
				for (integer js = 0; js < ic; js++) {
					if (vacant == set[js]) {
						vacant = NULL_SOSED;
					}
				}
				if (vacant != NULL_SOSED) {
					set[ic] = vacant;
					imax_pool_ind++;
					if (imax_pool_ind > imaxpoolsize) {
						printf("pool index incorrupt: pool ind > pool_size\n");
						getchar();
					}
					nnzRl++;
					ic++;
				}
			}
		}
		integer iscan = ii + 1;
		while ((iscan <= nnz_a[ilevel - 1] + iadd) && (A[iscan].i == set[0])) {
			// если узел j ещё не был добавлен в агрегат.
			if (flag[A[iscan].j] == false) {
				if ((A[iscan].j != set[0]) && (fabs(A[iscan].aij) >= theta*max_vnediagonal)) {
					vacant = A[iscan].j;
					pool_sosed[imax_pool_ind] = iscan;
					for (integer js = 0; js < ic; js++) {
						if (vacant == set[js]) {
							vacant = NULL_SOSED;
						}
					}
					if (vacant != NULL_SOSED) {
						set[ic] = vacant;
						imax_pool_ind++;
						if (imax_pool_ind > imaxpoolsize) {
							printf("pool index incorrupt: pool ind > pool_size\n");
							getchar();
						}
						nnzRl++;
						ic++;

					}
				}
			}

			iscan++;

		} // while

		// Это была учтена только связь i,j

		if (!blite_and_quick_version) {

			// Учёт свяи j,i
			//for (integer ii1 = 1 + iadd; ii1 <= nnz_a[ilevel - 1]; ii1++) {
				//if ((A[ii1].i != set[0]) && (A[ii1].j == set[0])) {
					//if (!flag[A[ii1].i]) {
						//if (fabs(A[ii1].aij) >= theta*max_vnediagonal) {
						//	vacant = A[ii1].i;
						//	for (integer js = 0; js < ic; js++) {
						//		if (vacant == set[js]) {
						//			vacant = NULL_SOSED;
						//		}
						//	}
						//	if (vacant != NULL_SOSED) {
						//		set[ic] = vacant; // j,i связь.

						//		nnzRl++;
						//		ic++;
						//	}
					//	}
					//}
				//}
			//}

			// Учёт связи j,i
			// Медленная версия на основе линейного поиска.
			//for (integer ii1 = 1 + iadd; ii1 <= nnz_a[ilevel - 1]; ii1++) {
				//if ((A[ii1].i != set[0]) && (A[ii1].j == set[0])) {
					//if (!flag[A[ii1].i]) {
						//if (fabs(A[ii1].aij) >= theta*max_vnediagonal) {
						//	vacant = A[ii1].i;
						//	for (integer js = 0; js < ic; js++) {
						//		if (vacant == set[js]) {
						//			vacant = NULL_SOSED;
						//		}
						//	}
						//	if (vacant != NULL_SOSED) {
						//		set[ic] = vacant; // j,i связь.

						//		nnzRl++;
						//		ic++;
						//	}
						//}
					//}
				//}
			//}

			// Учёт свяи j,i
			// Более быстрая версия на основе двоичного поиска.
			integer ii1 = BinarySearchAj(Atemp, set[0], 1, nnz_a[ilevel - 1]);
			for (integer ii2 = ii1; (ii2 <= nnz_a[ilevel - 1]) && (Atemp[ii2].j==set[0]); ii2++) {
				if ((Atemp[ii2].i != set[0]) && (!flag[Atemp[ii2].i])) {
					if (fabs(Atemp[ii2].aij) >= theta*max_vnediagonal) {
						vacant = Atemp[ii2].i;
						pool_sosed[imax_pool_ind] = Atemp[ii2].ind; // Позиция в А.
						for (integer js = 0; js < ic; js++) {
							if (vacant == set[js]) {
								vacant = NULL_SOSED;
							}
						}
						if (vacant != NULL_SOSED) {
							set[ic] = vacant; // j,i связь.
							imax_pool_ind++;
							if (imax_pool_ind > imaxpoolsize) {
								printf("pool index incorrupt: pool ind > pool_size\n");
								getchar();
							}
							nnzRl++;
							ic++;
						}
					}
				}
			}
			

		}


				for (integer isc = 1; isc < ic; isc++) {
					this_is_F_node[set[isc]] = true; // это только новые F узлы.
					bmarkervisit[set[isc]] = true;
				}




				// Помечаем узлы как включённые в агрегат.
				for (integer js = 0; js < ic; js++) {
					flag[set[js]] = true;
				}
				/*
				for (integer isc = 1; isc < n_a[ilevel - 1]; isc++) {
					if (flag[isc] == false) {
						// найти соседей узла isc 
						integer ii1 = BinarySearchAi(A, isc, 1 + iadd, nnz_a[ilevel - 1] + iadd);
						//for (integer ii1 = 1 + iadd; ii1 <= nnz_a[ilevel - 1] + iadd; ii1++) {
						//	if (A[ii1].i == isc) {
								integer ic2 = 0;
								for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == A[ii1].i); is0++) {
									for (int js = 1; js < ic; js++) {
										if (A[is0].j == set[js]) {
											ic2++;
										}
									}
								}
								count_sosed[isc] += ic2;
						//	}
						//}
						// если среди них окажется из set[1..ic-1] 
						// значит увеличиваем счётчик count_sosed[isc] на единицу.
					}
				}
				*/

				//if (1) {
				if (band_size == -1) {

					// Нет никакой информации о ширине ленты.

					for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
						flag_shadow[i_1] = flag[i_1];
					}
					// Модификация счётчиков для соседей новых F узлов.
					// j есть новый F тогда i счетчик меняется. // ij связь.
					//for (integer ii1 = 1 + iadd; ii1 <= nnz_a[ilevel - 1] + iadd; ii1++) {
					// Более быстрый вариант кода.
					for (integer ii_2 = ibegining_start_index_found_maximum; ii_2 <= nnz_a[ilevel - 1] + iadd; ii_2++) {
						integer ii1 = Atemp2[ii_2 - iadd].ind;
						integer isc = A[ii1].i;
						if ((flag_shadow[isc] == false) && ((!bmarkervisit[isc]))) {
							flag_shadow[isc] = true;
							integer ic2 = 0;
							for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == A[ii1].i); is0++) {
								for (int js = 1; js < ic; js++) {
									if (A[is0].j == set[js]) {
										ic2++;
									}
								}
							}
							count_sosed[isc] += ic2;
							if (bmaxsosedinfoactive) {
								// Обновляем информацию о максимальном количестве соседей.
								if (count_sosed[isc]>=Maximumsosedcount) {
									Maximumsosedcount = count_sosed[isc];
									icandidate_shadow = ii1;

									if (plist == NULL) {
										plist = new List;
										plist->next = NULL;
										plist->prev = NULL;
										plist->ii = ii1;
										plist->countsosed = count_sosed[isc];
										plist->i = isc;
										plist_current = plist;
									}
									else {
										List *ptemp=new List;
										ptemp->ii = ii1;
										ptemp->i = isc;
										ptemp->countsosed = count_sosed[isc];
										ptemp->next = NULL;
										ptemp->prev = plist_current;
										plist_current = ptemp;
										ptemp = NULL;
									}
								}
							}
						}
					}

				}
				else {

					// Здесь в полной мере используется гипотеза локальности.
					// Она состоит в следующем : узел с порядковым номером в матрице ii 
					// имеет соседей лишь в окрестности +- band_size.


					// Ширина ленты равна band_size
					//for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
						//flag_shadow[i_1] = flag[i_1];
					//}
					if (ic <= 0) {
						// Проверено ic не может быть нулём, а случай ic==1 
						// ничем не отличается от общего случая где прекрасно работает гипотеза локальности.
						// Выключим данную ветвь кода перенеся её в общий случай, т.к. данная 
						// ветвь очень медленная т.к. в ней не используется гипотеза локальности.
						// 7 ноября 2015.


						if (ic == 0) {
							printf("ic==0 iformation\n");
							getchar();
						}
						// На завершающей стадии соседей нет поэтому должен работать
						// код которому ничего неизвестно о ширине ленты.


						// Нет никакой информации о ширине ленты.

						for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
							flag_shadow[i_1] = flag[i_1];
						}
						// Модификация счётчиков для соседей новых F узлов.
						// j есть новый F тогда i счетчик меняется. // ij связь.
						//for (integer ii1 = 1 + iadd; ii1 <= nnz_a[ilevel - 1] + iadd; ii1++) {
						// Более быстрый вариант кода.
						for (integer ii_2 = ibegining_start_index_found_maximum; ii_2 <= nnz_a[ilevel - 1] + iadd; ii_2++) {
							integer ii1 = Atemp2[ii_2 - iadd].ind;
							integer isc = A[ii1].i;
							if ((flag_shadow[isc] == false) && ((!bmarkervisit[isc]))) {
								flag_shadow[isc] = true;
								integer ic2 = 0;
								for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == A[ii1].i); is0++) {
									for (int js = 1; js < ic; js++) {
										if (A[is0].j == set[js]) {
											ic2++;
										}
									}
								}
								count_sosed[isc] += ic2;
								if (bmaxsosedinfoactive) {
									// Обновляем информацию о максимальном количестве соседей.
									if (count_sosed[isc]>=Maximumsosedcount) {
										Maximumsosedcount = count_sosed[isc];
										icandidate_shadow = ii1;

										if (plist == NULL) {
											plist = new List;
											plist->next = NULL;
											plist->prev = NULL;
											plist->ii = ii1;
											plist->countsosed = count_sosed[isc];
											plist->i = isc;
											plist_current = plist;
										}
										else {
											List *ptemp = new List;
											ptemp->ii = ii1;
											ptemp->i = isc;
											ptemp->countsosed = count_sosed[isc];
											ptemp->next = NULL;
											ptemp->prev = plist_current;
											plist_current = ptemp;
											ptemp = NULL;
										}
									}
								}
							}
						}

					}
					else {

						if (0) {
							// медленный вариант кода.


							for (int js = 1; js < ic; js++) {
								for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
									flag_shadow[i_1] = flag[i_1];
								}

								integer istart = i_my_max(1 + iadd, pool_sosed[js] - band_size - 1);
								integer iend = i_my_min(nnz_a[ilevel - 1] + iadd, pool_sosed[js] + band_size + 1);
								// Ищем только среди ближайшего окружения вновь добавленного F узла.
								for (integer ii_2 = istart; ii_2 <= iend; ii_2++) {
									integer isc = A[ii_2].i;
									if ((flag_shadow[isc] == false)/* && ((!bmarkervisit[isc]))*/) {
										flag_shadow[isc] = true;
										integer ic2 = 0;
										integer iend2 = i_my_min(nnz_a[ilevel - 1] + iadd, ii_2 + band_size + 1);
										integer istart2 = ii_2;
										while ((istart2 >= 1 + iadd) && (A[istart2].i == A[ii_2].i)) istart2--;
										istart2++;
										for (integer is0 = istart2; (is0 <= iend2) && (A[is0].i == A[ii_2].i); is0++) {
											if (A[is0].j == set[js]) {
												ic2++;
											}
										}

										count_sosed[isc] += ic2;
										if (bmaxsosedinfoactive) {
											// Обновляем информацию о максимальном количестве соседей.
											if (count_sosed[isc] >= Maximumsosedcount) {
												Maximumsosedcount = count_sosed[isc];
												icandidate_shadow = ii_2;

												if (plist == NULL) {
													plist = new List;
													plist->next = NULL;
													plist->prev = NULL;
													plist->ii = ii_2;
													plist->countsosed = count_sosed[isc];
													plist->i = isc;
													plist_current = plist;
												}
												else {
													List *ptemp = new List;
													ptemp->ii = ii_2;
													ptemp->i = isc;
													ptemp->countsosed = count_sosed[isc];
													ptemp->next = NULL;
													ptemp->prev = plist_current;
													plist_current = ptemp;
													ptemp = NULL;
												}
											}
										}
									}
								}
							}
						}
						else {

							for (int js = 1; js < ic; js++) {
								//for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
									//flag_shadow[i_1] = flag[i_1];
								//}

								//integer istart = i_my_max(1 + iadd, pool_sosed[js] - band_size - 1);
								//integer iend = i_my_min(nnz_a[ilevel - 1] + iadd, pool_sosed[js] + band_size + 1);

								//integer i3 = 1;
								//while ((i3 <= n_a[ilevel - 1]) && (row_startA[i3] < istart)) i3++;
								//integer i4 = n_a[ilevel - 1];
								//while ((i4 >= 1) && (row_startA[i4]>iend)) i4--;

								integer i3 = i_my_max(1, set[js] - band_size_i - 1);
								integer i4 = i_my_min(n_a[ilevel - 1], set[js] + band_size_i + 1);

								// Ищем только среди ближайшего окружения вновь добавленного F узла.
								//for (integer ii_2 = istart; ii_2 <= iend; ii_2++) {

								// ищем соседа узла set[js].
								for (integer i5 = i3; i5 <= i4; i5++) {
									//integer isc = A[ii_2].i;
									integer ii_2 = row_startA[i5];
									//integer ii_2 = BinarySearchAi(A, i5, 1 + iadd, nnz_a[ilevel - 1] + iadd);
									integer isc = i5;
									//if ((flag_shadow[isc] == false)/* && ((!bmarkervisit[isc]))*/) {
									if (flag[isc]==false) {
										//flag_shadow[isc] = true;
										integer ic2 = 0;
										integer iend2 = i_my_min(nnz_a[ilevel - 1] + iadd, ii_2 + band_size + 1);
										integer istart2 = ii_2;
										while ((istart2 >= 1 + iadd) && (A[istart2].i == A[ii_2].i)) istart2--;
										istart2++;
										for (integer is0 = istart2; (is0 <= iend2) && (A[is0].i == A[ii_2].i); is0++) {
											if (A[is0].j == set[js]) {
												ic2++;
											}
										}

										count_sosed[isc] += ic2;
										if (bmaxsosedinfoactive) {
											// Обновляем информацию о максимальном количестве соседей.
											if (count_sosed[isc] >= Maximumsosedcount) {
												Maximumsosedcount = count_sosed[isc];
												icandidate_shadow = ii_2;

												if (plist == NULL) {
													plist = new List;
													plist->next = NULL;
													plist->prev = NULL;
													plist->ii = ii_2;
													plist->countsosed = count_sosed[isc];
													plist->i = isc;
													plist_current = plist;
												}
												else {
													List *ptemp = new List;
													ptemp->ii = ii_2;
													ptemp->i = isc;
													ptemp->countsosed = count_sosed[isc];
													ptemp->next = NULL;
													ptemp->prev = plist_current;
													plist_current = ptemp;
													ptemp = NULL;
												}
											}
										}
									}
								}
							}
						}
					}

				}


				

				n_coarce++; // Увеличено количество С узлов.

				// Один агрегат создан.

			} // узел не был ещё включён в агрегат.


		
        bcontinue = false;
		for (integer i_1 = istartflag_scan; i_1 <= n_a[ilevel - 1]; i_1++) {
			if (flag[i_1] == false) {
				bcontinue = true;
				istartflag_scan = i_1; // сокращаем пределы сканирования.
				break; // досрочный выход из цикла for.
				if (maxsosed == -1) {
					printf("ERROR!!!!  i_1=%d\n", i_1);
					system("pause");
				}
			}
        }

			// Вычисление узла с максимальным количеством соседей.
			maxsosed = -1;
			icandidate = 0;


			// TODO 6 november
			
			const integer ipool_size_limit = 64000;
			integer ipool[ipool_size_limit];
			integer isize_p = -1;
			for (integer isc = 0; isc < ic; isc++) {
				//integer ii_s = pool_sosed[isc];// позиция в А.
				integer ii_s = BinarySearchAi(A, set[isc], 1 + iadd, iadd+nnz_a[ilevel - 1]);
				integer ii_c = ii_s;
				//while ((ii_c >= 1 + iadd) && (A[ii_c].i == A[ii_s].i)) ii_c--;
				//ii_c++;
				// Вся строка матрицы начиная с позиции ii_c должна быть помещена в Atemp2;
				while ((ii_c <= iadd + nnz_a[ilevel - 1]) && (A[ii_c].i == A[ii_s].i)) {
					integer icandidateq = ii_c;
					bool found1 = false;
					for (integer i_7 = 0; i_7 <= isize_p; i_7++) {
						if (i_7 < ipool_size_limit) {
							if (ipool[i_7] == icandidateq) {
								found1 = true;
							}
						}
						else {
							printf("nado uvelichiti ipool bolee %d, string 4118\n", ipool_size_limit);
							system("pause");
							exit(1);
						}
					}
					// элемент в списке не обнаружен, поэтому запомним индекс.
					if (found1 == false) {
						isize_p++;
						if (isize_p < ipool_size_limit) {
							ipool[isize_p] = icandidateq;
						}
						else {
							printf("nado uvelichiti ipool bolee %d, string 4118\n", ipool_size_limit);
							system("pause");
							exit(1);
						}
					}
					ii_c++;
				}
			}
			// ipool хранит кандидатов для добавления в Atemp2.
			for (integer i_7 = 0; i_7 <= isize_p; i_7++) {
				// Позиции элементов в А остаются на месте, поэтому коллизий быть недолжно.

				Ak temp = Atemp2[ibegining_start_index_found_maximum - iadd + i_7];
				Atemp2[ibegining_start_index_found_maximum - iadd + i_7] = A[ipool[i_7]];
				Atemp2[ibegining_start_index_found_maximum - iadd + i_7].ind = ipool[i_7];
				
				Atemp2[A[ipool[i_7]].ind] = temp;
				A[temp.ind].ind = A[ipool[i_7]].ind;
				A[ipool[i_7]].ind = ibegining_start_index_found_maximum - iadd + i_7;
			}
			ibegining_start_index_found_maximum += isize_p + 1;
			


			// Мы будем делать досрочный выход из цикла for 19 раз из 20,
			// и лишь один раз апдейтить информацию о максимальном количестве соседей.
			// Данная модификация придумана 17 октября 2015 года.
			// Эта эвристика обеспечила ускорение на 8% от суммарного времени исполнения 
			// всего алгоритма.
			if (icountprohod % 200 == 0) {

				

				// TODO 6november start delete code
				/*
				// Сначала пишем те узлы которые были посещены.
				integer i_2 = 1;
				for (integer i_1 = 1 + iadd; i_1 <= nnz_a[ilevel - 1] + iadd; i_1++) {
					if (bmarkervisit[A[i_1].i]) {
						Atemp2[i_2] = A[i_1];
						Atemp2[i_2].ind = i_1; // обязательно запоминаем первоначальный индекс.
						i_2++;
					}
				}
				// Потом пишем все другие узлы.
				for (integer i_1 = 1 + iadd; i_1 <= nnz_a[ilevel - 1] + iadd; i_1++) {
					if (!bmarkervisit[A[i_1].i]) {
						Atemp2[i_2] = A[i_1];
						Atemp2[i_2].ind = i_1; // обязательно запоминаем первоначальный индекс.
						i_2++;
					}
				}
				// TODO 6 november end delete code.
				*/

				/* // 7 novemver 2015.
				bool bfirst_loc = true; // нам нужно именно первое значение позиции.

				for (integer i_2 = ibegining_start_index_found_maximum; i_2 <= nnz_a[ilevel - 1] + iadd; i_2++) {
					integer  i_1 = i_2 - iadd;


					if (!bmarkervisit[Atemp2[i_1].i]) {
						// Мы запоминаем стартовую позицию начиная с которой начинаются ещё
						// непомеченные узлы.
						if (bfirst_loc) {
							ibegining_start_index_found_maximum = i_1 + iadd;
							bfirst_loc = false;
							break;
							//printf("diagnostic =%d\n",i_1);
							//getchar();
						}
					}
				}
				*/

				
				if ((flag[A[icandidate_shadow].i] == false) && (count_sosed[A[icandidate_shadow].i] == Maximumsosedcount)) {
					// Дело в том что мы уже нашли узел с наибольшим числом соседей прямо в ходе
					// модификации счётчиков соседей новых F узлов.
					// просто сразу это было неочевидно поэтому мы организовали поиски.
					icandidate = icandidate_shadow;
				}
				else {

					for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
						flag_shadow[i_1] = flag[i_1];
					}


					for (integer i_2 = ibegining_start_index_found_maximum; i_2 <= nnz_a[ilevel - 1] + iadd; i_2++) {
						integer  i_1 = i_2 - iadd;



						if (flag_shadow[Atemp2[i_1].i] == false) {
							if (count_sosed[Atemp2[i_1].i] > maxsosed) {
								maxsosed = count_sosed[Atemp2[i_1].i];
								//icandidate = i_1;
								icandidate = Atemp2[i_1].ind;
							}
							flag_shadow[Atemp2[i_1].i] = true;
						}
					}

					Maximumsosedcount = maxsosed;
					bmaxsosedinfoactive = true;
				}

			}
			else {

				if ((flag[A[icandidate_shadow].i] == false) && 
					  (count_sosed[A[icandidate_shadow].i] == Maximumsosedcount)) {
					// Дело в том что мы уже нашли узел с наибольшим числом соседей прямо в ходе
					// модификации счётчиков соседей новых F узлов.
					// просто сразу это было неочевидно поэтому мы организовали поиски.
					icandidate = icandidate_shadow;
				}
				else {

					
					bool found_candidate = false;
					if (plist!=NULL) {
						// Сканируем короткий пулл кандидатов.
						while ((plist_current != NULL) && (!(((flag[A[plist_current->ii].i] == false) &&
							(count_sosed[A[plist_current->ii].i] == Maximumsosedcount))))) {
							List* temp;
							temp = plist_current;
							plist_current = plist_current->prev;
							if (plist_current != NULL) {
								plist_current->next = NULL;
							}
							temp->prev = NULL;
							if (plist != temp) {
								delete temp;
							}
							else {
								plist = NULL;
								delete temp;
							}
						}
						if (plist_current != NULL) {
							icandidate_shadow = plist_current->ii;
							icandidate = icandidate_shadow;
							found_candidate = true;
						}
					}

					if (found_candidate==false) 
					{

						for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) {
							flag_shadow[i_1] = flag[i_1];
						}




						// если это не так мы всегда можем запустить старый проверенный код.

						//bool bfirst_loc = true; // нам нужно именно первое значение позиции.

						for (integer i_2 = ibegining_start_index_found_maximum; i_2 <= nnz_a[ilevel - 1] + iadd; i_2++) {


							integer i_1 = i_2 - iadd;
							// слишком частые коректировки начала списка замедляют код вдвое.
							// мы будем коректировать начальную позицию раз в 20 раз.
							//if (!bmarkervisit[Atemp2[i_1].i]) {
							// Мы запоминаем стартовую позицию начиная с которой начинаются ещё
							// непомеченные узлы.
							//if (bfirst_loc) {
							//	ibegining_start_index_found_maximum = i_1+iadd;
							//	bfirst_loc = false;
							//}
							//}

							if (flag_shadow[Atemp2[i_1].i] == false) {

								if (count_sosed[Atemp2[i_1].i] > maxsosed) {
									maxsosed = count_sosed[Atemp2[i_1].i];
									//icandidate = i_1;
									// icandidate номер элемента в матрице А отвечающий за найбольшее число соседей.
									icandidate = Atemp2[i_1].ind;
									if (bmaxsosedinfoactive) {
										// Организуем досрочный выход из цикла,
										// это длжно ускорить построение C-F разбиения.
										if (maxsosed == Maximumsosedcount) break;
									}

								}
								flag_shadow[Atemp2[i_1].i] = true;
							}
						}

						// statistics
						//integer i_11 = 0;
						//for (integer i_10 = 1; i_10 <= n_a[ilevel - 1]; i_10++) {
						//if (flag[i_10] == false) i_11++;
						//}
						//printf("%d %d\n", A[icandidate_shadow].i, icandidate); // 97%
						//printf("procent %f\n",(float)(100*i_11/n_a[ilevel-1]));
						//getchar();
					}

				}
			}

			//printf("maximum number of sosed=%d\n",maxsosed);
			if (maxsosed == -1) if (debug_reshime) system("pause");
			//getchar();

			if ((icandidate == 0) && (maxsosed == -1)) {
				bcontinue = false;
			}

			icountprohod++;

		} // Построение C-F разбиения. создано.

		delete[] bmarkervisit;
		// Уничтожение памяти из по пула кандидатов.
		if (plist != NULL) {
			
			while (plist != NULL) {
				plist_current = plist;
				plist = plist->next;
				plist_current->next = NULL;
				plist->prev = NULL;
				delete plist_current;
			}
		}

		// В методе стандартной интерполяции присутствует шаг уменьшения разреженности,
		// для того чтобы правильно аппроксимировать все F переменные C переменными надо
		// увеличить количество С переменных.
		int ipromah = 0;
		int ipromah_one = 0;
		bool bweSholdbeContinue = true;
		while (bweSholdbeContinue) {
			bweSholdbeContinue = false;

			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) if (this_is_F_node[i_1] == true) {
				// Найти соседей данного F-node которые C-node.
				integer icsos = 0;
				integer i_2 = BinarySearchAi(A, i_1, 1 + iadd, nnz_a[ilevel - 1] + iadd);
			
				bool bvisit = false;
				for (integer is0 = i_2; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == A[i_2].i); is0++) {
					if (A[is0].j != A[i_2].i) {
						bvisit = true;
						if (this_is_C_node[A[is0].j] == true) {
							icsos++;
						}
						else {
							ipromah++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) ipromah_one++; // количество F узлов с одним единственным С соседом.
				// Если bvisit то внедиагональные элементы есть но они все Fnodes. Иначе там обособленное условие Дирихле.
				if ((icsos == 0) && (bvisit)) {

					// А если он F узел дирихле без соседей, то сумма тоже может быть нулевой и это вызовет деление на ноль.
					// Узлы Дирихле могли быть без соседей на начальных уровнях, они располагались в конце списка и были
					// поглощены агломератами внутренних узлов и всё было впорядке.
					// Чтобы преодолеть это затруднение нужен алгоритм с обратной связью.

					// Нет С соседей, этот узел станет С узлом.
					this_is_F_node[i_1] = false;
					this_is_C_node[i_1] = true;
					// F node стал C_node!!! Идея стандартной интерполяции 
					// приводит к уменьшению разреженности оператора Галёркина.
					bweSholdbeContinue = true;
				}

			}

			if (bweSholdbeContinue) {
				printf(" prohod succseful\n");
			}
			else {
				printf("prohod empty\n");
			}

		}


		// Нужно корректно обработать узлы Дирихле,
		// Если F узел окажется узлом Дирихле без соседей то его надо сделать С узлом,
		// Но узнать такой узел можно лишь в процессе выполнения алгоритма дальше по ходу исполнения.
		// Поэтому может потребоваться вернуться и начать заново (обратная связь).


		integer* C_numerate =  new integer[n_a[ilevel - 1] + 1];
		integer icounter=1;
		integer icount1;
		integer numberofcoarcenodes;
		Real* ap_coarse = NULL;

		bweSholdbeContinue = true;
		while (bweSholdbeContinue) {
			bweSholdbeContinue = false;

			for (int i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) if (this_is_C_node[i_1]) { n_coarce++; }
			n_coarce--;


			// debug
			// проверка качества C-F разбиения.
			//Real* exp1 = new Real[n_a[ilevel - 1] + 1];
			//for (int i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) exp1[i_1] = 0.0;
			//for (int i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) if (this_is_C_node[i_1]) exp1[i_1] = 2.0;
			//for (int i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) if (this_is_F_node[i_1]) exp1[i_1] = 1.0;
			//exporttecplot(exp1,n);
			//delete[] exp1;

			//printf("export ready");

			//system("pause");


			// C-F разбиение построено, самое время построить оператор интерполяции.
			// потом найти оператор проекции, как транспонированный оператор интерполяции.
			// Всё завершает построение матрицы нового сеточного уровня и можно запускать новый уровень.

			// Построение оператора интерполляции : 
			// coarse 2 fine.
			//P*coarse==fine
			
			for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) C_numerate[i_1] = 0;
			icounter = 1;
			for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) if (this_is_C_node[i_1] == true) {
				//printf("C ind= %d", i_1); getchar();
				C_numerate[i_1] = icounter;
				icounter++;
			}




			// C_numerate - перенумерация на множестве Coarse узлов.
			// Построение пролонгации для узлов которые составляют грубую сетку.
			icount1 = 1 + iaddR; // nnz_R
			for (integer i_1 = 1; i_1 <= n_a[ilevel - 1]; i_1++) if (this_is_C_node[i_1] == true) {
				P[icount1].aij = 1.0;
				P[icount1].i = C_numerate[i_1]; // coarse number
				P[icount1].j = i_1; // fine number.
				icount1++;
			}

			// значение icount1 нужно далее.НЕ трогать !!!.
			numberofcoarcenodes = icount1 - 1 - iaddR;

			// Для модификации R  надо transpose(P)/ap.
			printf("countloc=%d\n", numberofcoarcenodes);
			if (debug_reshime) system("pause");
			
			ap_coarse = new Real[numberofcoarcenodes + 1];
			if (ap_coarse == NULL) {
				printf("error cannot memory allocate.");
				system("pause");
			}
			ap_coarse[0] = 0.0;
			//integer countloc = 1;
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_C_node[i8] == true) {
				integer ii1 = BinarySearchAi(A, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				//integer ii2 = ii1 - 1;
				//if ((ii2 >= 1 + iadd) && (A[ii2].i == A[ii1].i)) {
				//printf("koren zla\n"); // бинарный поиск должен гарантированно находить самого левого представителя.
				//getchar();
				//}
				for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == A[ii1].i); is0++) {
					//printf("i=%d j=%d A[is0].aij=%e ", A[is0].i, A[is0].j, A[is0].aij);
					if (A[is0].j == A[ii1].i) {
						//if (countloc > icount1 - 1) { printf("system error\n"); getchar(); }
						if (fabs(A[is0].aij) > RealMAXIMUM) {
							printf("perepolnenie error!");
							getchar();
						}
						ap_coarse[C_numerate[i8]] = fabs(A[is0].aij);
						//printf("find = %e", fabs(A[is0].aij));
					}
				}
				//printf("\n");
				//getchar();

				//countloc++;
			}

			// верно 2 октября.
			//for (integer i25 = 1; i25 < icount1; i25++) {
			//if (ap_coarse[i25]>1) {
			//printf("ap_coarse[%d]=%e\n", i25, ap_coarse[i25]);
			//getchar();
			//}
			//}

			ipromah = 0;
			ipromah_one = 0;
			// Построение пролонгации для узлов которые составляют F nodes.
			// Каждый F-nodes окружён C-nodes.
			for (integer i8 = 1; i8 <= n_a[ilevel - 1]; i8++) if (this_is_F_node[i8] == true) {
				// Найти соседей данного F-node которые C-node.
				integer icsos = 0;
				integer ii1 = BinarySearchAi(A, i8, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				Real sumP = 0.0;
				for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == A[ii1].i); is0++) {
					if (A[is0].j != A[ii1].i) {
						if (this_is_C_node[A[is0].j] == true) {
							sumP += fabs(A[is0].aij); // сумма модулей внедиагональных элементов которые принадлежат С узлам.
							icsos++;
						}
						else {
							ipromah++; // подсчитываем проблемы интерполяции 
						}
					}
				}
				if (icsos == 1) ipromah_one++; // количество F узлов с одним единственным С соседом.

				for (integer is0 = ii1; (is0 <= nnz_a[ilevel - 1] + iadd) && (A[is0].i == A[ii1].i); is0++) {
					if (A[is0].j != A[ii1].i) {
						if (this_is_C_node[A[is0].j] == true) {

							// Внедиагональный элемент из множества С узлов.

							if (fabs(sumP) < RealZERO) {
								printf("error interpolation zero diagonal sumP.\n");
								printf("Fnode all sosed is F");
								//system("pause");
								printf("i8 is Dirichlet node\n");
								this_is_F_node[i8] = false; // Этот узел Дирихле станет С нодом.
								this_is_C_node[i8] = true;
								bweSholdbeContinue = true;
								//exit(1);
								// здесь нужна непрямая интерполляция.
							}
							else {

								P[icount1].j = i8;
								P[icount1].i = C_numerate[A[is0].j];
								P[icount1].aij = fabs(A[is0].aij) / sumP;
								icount1++;
							}
						}

					}
				}


			}

			if (bweSholdbeContinue) {
				delete[] ap_coarse;
				ap_coarse = NULL;
				printf("obratnaq svqz restart...\n");
			}

		}

		nnzR = icount1-iaddR;

		

		// нужно определить nnzR количество ненулевых элементов в матрице R и P.

		// оператор restriction построен и он упорядочен по i.
		// число ненулевых элементов nnzR-1.
		// P=Copy(R);
		for (integer ii = 1 + iaddR; ii <= iaddR + nnzR - 1; ii++) {
			R[ii] = P[ii];
			//if ((ilevel == 2) && (ii >= iaddR + nnzR - 1-30)) {
				//printf("ii=%d aij=%e, i=%d j=%d\n", ii, P[ii].aij, P[ii].i, P[ii].j);
				//getchar();
			//}
		}

		// heapsort(P,key==j,iaddR+1,iaddR+nnzR - 1);
		if (bquicktypesort) {
			QuickSort_j(P, 1 + iaddR, iaddR + nnzR - 1);
		}
		else {
			HeapSort_j(P, /*n_a[ilevel - 1],*/ 1 + iaddR, iaddR + nnzR - 1);
		}

		/*
		if (bquicktypesort) {
			QuickSort_j(R,  1 + iaddR, iaddR + nnzR - 1);
		}
		else {
			HeapSort_j(R, n_a[ilevel - 1], 1 + iaddR, iaddR + nnzR - 1);
		}

		printf("start now\n");
		for (integer ii5 = 1 + iaddR; ii5 <= iaddR + nnzR - 1; ii5++) {
			printf("R[%d].i=%d R[%d].j=%d R[%d].aij=%e\n", ii5, R[ii5].i, ii5, R[ii5].j, ii5, R[ii5].aij);
			system("pause");
		}
		*/
		// где то надо разделить на ap, т.к. 
		// R=P/ap.

		
		// heapsort(R,key==i,iaddR+1,iaddR+nnzR - 1);
		if (bquicktypesort) {
			QuickSort(R,  1 + iaddR, iaddR + nnzR - 1);
		}
		else {
			HeapSort(R, n_a[ilevel - 1], 1 + iaddR, iaddR + nnzR - 1);
		}

		
		printf("first level size n=%d numberofcoarcenodes=%d\n", n, numberofcoarcenodes);
		if (debug_reshime) system("pause");
		for (integer i_1 = 1; i_1 <= numberofcoarcenodes; i_1++) {
			flag[i_1] = false; // init flag.
		}
		
		for (integer i_1 = 1 + iaddR; i_1 <= iaddR + nnzR - 1; i_1++) {
			if (flag[R[i_1].i] == false)
			{
				for (integer i_2 = i_1; (i_2 <= iaddR + nnzR - 1) && (R[i_2].i == R[i_1].i); i_2++) {
					if ((R[i_1].i<1) || (R[i_1].i>numberofcoarcenodes)) {
						printf("error R[%d].i=%d\n",i_1,R[i_1].i);
						system("pause");
					}
					if (fabs(ap_coarse[R[i_1].i]) < divisionZERO) {
						printf("error division by zero\n");
						system("pause");
					}
					Real delitel;
					if (ap_coarse[R[i_1].i] > 1.0) {
						// internal node
						delitel = 0.5*ap_coarse[R[i_1].i];
					}
					else {
						// Dirichlet
						delitel = ap_coarse[R[i_1].i];
					}
					R[i_2].aij = R[i_2].aij / (delitel);
				}
				flag[R[i_1].i] = true;
			}
		}
	

		delete[] ap_coarse;
		

		// debug.
		//for (integer i_1 = 1 + iaddR; i_1 <= iaddR + nnzR - 1; i_1++) {
			//printf("R[%d].i=%d R[%d].j=%d R[%d].aij=%e\n", i_1, R[i_1].i, i_1, R[i_1].j, i_1, R[i_1].aij);
			//system("pause");
		//}

		
		//for (integer i_1 = iaddR + nnzR - 1-20; i_1 <= iaddR + nnzR - 1; i_1++) {
			//printf("R[%d].i=%d R[%d].j=%d R[%d].aij=%e\n", i_1, R[i_1].i, i_1, R[i_1].j, i_1, R[i_1].aij);
			//system("pause");
		//}
		
		
		
		//for (integer i_1 = 1 + iaddR; i_1 <= iaddR + nnzR - 1; i_1++) {
			//printf("P[%d].i=%d P[%d].j=%d P[%d].aij=%e\n", i_1, P[i_1].i, i_1, P[i_1].j, i_1, P[i_1].aij);
			//system("pause");
		//}
		

		// Сортировка А по j.
		//heapsort(A, key=j*n_a[ilevel - 1] + i, 1, nnz_a[ilevel - 1]);
		if (bquicktypesort) {
			QuickSort_j(A, 1 + iadd, nnz_a[ilevel - 1] + iadd);
		}
		else {
			if (nnz_a[ilevel - 1] < heapsortsizelimit) {
				HeapSort_j(A, /*n_a[ilevel - 1]*/ 1 + iadd, nnz_a[ilevel - 1] + iadd);
			}
			else {
				Ak* Aorig = &A[1 + iadd];
				MergeSort_j(Aorig, nnz_a[ilevel - 1]);
				Aorig = NULL;
			}
		}

		// Нахождение матрицы грубосеточного уровня :
		// Acorse=R*Afine*P;
		// часть 1 : R*Afine.
		//         xxxxxx
        //         xxxxxx
		//  xxxxxx xxxxxx xxxxxx
		//  xxxxxx xxxxxx xxxxxx
		//         xxxxxx
		//         xxxxxx
		//    R       A     [RA]
		/*
		integer istartAnew = nnz_a[ilevel - 1] + 1+iadd;
		for (integer jstr = 1; jstr <= n_a[ilevel - 1]; jstr++) {
			// jstr - столбец матрицы А.
			// icounter-1 - число узлов на грубом уровне.
			for (integer i = 1; i <= icounter - 1; i++) {
				flag[i] = false;
			}
			// А (на первой позиции) должна быть отсортирована по j.
			integer ii1 = BinarySearchAj(A, jstr, 1 + iadd, nnz_a[ilevel - 1] + iadd);
			// это справедливо только для первого уровня.
			integer istart = 1 + iaddR;
			integer iend = nnzR - 1 + iaddR;
			for (integer ii = istart; ii <= iend; ii++) {
				if (flag[R[ii].i] == false) {
					integer istr = R[ii].i;
					integer ic = ii;
					// i-coarse, j-fine
					Real sum1 = 0.0;
					while ((ic <= iend) && (R[ic].i == istr)) {
						// R[R[ii].i][R[ic].j]*A[R[ic].j][jstr]

						for (integer ii2 = ii1; (ii2 <= nnz_a[ilevel - 1] + iadd) && (A[ii2].j == jstr); ii2++) {
							// (A[ii2].j == jstr) {
							if (A[ii2].i == R[ic].j) {
								sum1 += R[ic].aij*A[ii2].aij;
							}
							//}
						}
						ic++;
					}
					if (fabs(sum1) > RealZERO) {
						A[istartAnew].aij = sum1;
						A[istartAnew].i = istr;
						A[istartAnew].j = jstr;
						if (jstr < 0) {
							printf("fatal error");
							printf("i=%d j=%d aij=%e\n", istr, jstr, sum1);
							system("pause");
						}
						istartAnew++;
					}
					flag[R[ii].i] = true;
				}
			}
		} 
		*/

		// Идея droptolerance состоит в отсечении в матрице элементов 
		// которые меньше чем одна тысячная от диагонального
		// А если диагонального элемента нет то берётся норма Чебышева строки.
		//
		//
		//

		// часть 1 : R*Afine.
		//         xxxxxx
		//         xxxxxx
		//  xxxxxx xxxxxx xxxxxx
		//  xxxxxx xxxxxx xxxxxx
		//         xxxxxx
		//         xxxxxx
		//    R       A     [RA]
		//integer istartAnew = nnz_a[ilevel - 1] + 1 + iadd;
		//for (integer jstr = 1; jstr <= n_a[ilevel - 1]; jstr++) {
			// jstr - столбец матрицы А.
			// icounter-1 - число узлов на грубом уровне.
			//if ((icounter - 1 > n) || (icounter - 1 < 0)) {
				//printf("flag incorrupt 4...\n");
				//system("pause");
				//exit(1);
			//}
			//for (integer i = 1; i <= icounter - 1; i++) {
				//flag[i] = false;
			//}
			// А (на первой позиции) должна быть отсортирована по j.
			//integer ii1 = BinarySearchAj(A, jstr, 1 + iadd, nnz_a[ilevel - 1] + iadd);
			// это справедливо только для первого уровня.
			//integer istart = 1 + iaddR;
			//integer iend = nnzR - 1 + iaddR;
			//for (integer ii = istart; ii <= iend; ii++) {
				//if ((R[ii].i > n) || (R[ii].i < 0)) {
					//printf("flag incorrupt 3...\n");
					//system("pause");
					//exit(1);
				//}
				//if (flag[R[ii].i] == false) {
					//integer istr = R[ii].i;
					//integer ic = ii;
					// i-coarse, j-fine
					//Real sum1 = 0.0;
					
					// это медленная версия кода.
					//while ((ic <= iend) && (R[ic].i == istr)) {
						// R[R[ii].i][R[ic].j]*A[R[ic].j][jstr]

						//for (integer ii2 = ii1; (ii2 <= nnz_a[ilevel - 1] + iadd) && (A[ii2].j == jstr); ii2++) {
							
							//if (A[ii2].i == R[ic].j) {
								//sum1 += R[ic].aij*A[ii2].aij;
							//}
							
						//}
						//ic++;
					//}
					
					// Это более быстрый код.
					//integer ks = ic;
					//integer ls = ii1;
					//integer kf = ic;
					//bool bvis = false;
					//Real retalon = 0.0;
					//while ((kf <= iend) && (R[kf].i == istr)) {
						//if (R[kf].j == istr) {
						//	retalon = fabs(R[kf].aij);
						//	bvis = true;
						//}
						//kf++;
					//}
					//kf--;
					
					//if (bvis == false) {
						//kf = ic;
						//while ((kf <= iend) && (R[kf].i == istr)) {
						//	if (fabs(R[kf].aij) > retalon) retalon = fabs(R[kf].aij);
						//	kf++;
						//}
						//kf--;
					//}

					//integer lf = ii1;
					//for (integer ii2 = ii1; (ii2 <= nnz_a[ilevel - 1] + iadd) && (A[ii2].j == jstr); ii2++) {
						//if (fabs(A[ii2].aij) > retalon) retalon = fabs(A[ii2].aij);
						//if (A[ii2].i == istr) retalon *= fabs(A[ii2].aij);
						//lf++;
					//}
					//lf--;
					
					//while ((ks <= kf) && (ls <= lf)) {
						//if (A[ls].i < R[ks].j) {
						//	ls++;
						//}
						//else if (A[ls].i > R[ks].j) {
						//	ks++;
						//}
						//else {
						//	sum1 += R[ks].aij*A[ls].aij;
						//	ks++;
						//	ls++;
						//}
					//}




					//if (fabs(sum1) > 0.001*retalon) {
					//if (fabs(sum1)>1.0e-30) {
						//A[istartAnew].aij = sum1;
						//A[istartAnew].i = istr;
						//A[istartAnew].j = jstr;
						//if (jstr < 0) {
						//	printf("fatal error");
						//	printf("i=%d j=%d aij=%e\n", istr, jstr, sum1);
						//	system("pause");
						//}
						//istartAnew++;
					//}
					//flag[R[ii].i] = true;
				//}
			//}
		//}


        printf("nnz left operand=%d, nnz right operand=%d\n",nnzR,nnz_a[ilevel-1]);

		integer istartAnew;
		integer* kf_array;

	

		if (0) {

			// Самая быстрая версия на основе слияния упорядоченных массивов.

			//Более быстрый вариант алгоритма.
			// Быстрее этого кода на основе идеи слияния списков уже не будет.
			// 17 октября 2015. Нужно двигаться в сторону Писсанецки.
			// часть 1 : R*Afine.
			//         xxxxxx
			//         xxxxxx
			//  xxxxxx xxxxxx xxxxxx
			//  xxxxxx xxxxxx xxxxxx
			//         xxxxxx
			//         xxxxxx
			//    R       A     [RA]
			 kf_array = new integer[numberofcoarcenodes + 1];
			integer istart1 = 1 + iaddR;
			integer iend1 = nnzR - 1 + iaddR;
			for (integer i = 1; i <= icounter - 1; i++) {
				flag[i] = false;
			}
			for (integer ii = istart1; ii <= iend1; ii++) {
				if (flag[R[ii].i] == false) {
					integer istr = R[ii].i;
					integer ic = ii;
					integer kf = ic;

					while ((kf <= iend1) && (R[kf].i == istr)) {
						kf++;
					}
					kf--;
					kf_array[istr] = kf;
					flag[R[ii].i] = true;
				}
			}

			integer *start_position_i_string_in_R = new integer[numberofcoarcenodes + 1];
			for (integer i = 1; i <= icounter - 1; i++) {
				flag[i] = false;
			}
			integer istart3 = 1 + iaddR;
			integer iend3 = nnzR - 1 + iaddR;
			for (integer ii = istart3; ii <= iend3; ii++) {
				if (flag[R[ii].i] == false) {
					start_position_i_string_in_R[R[ii].i] = ii;
					flag[R[ii].i] = true;
				}
			}


			 istartAnew = nnz_a[ilevel - 1] + 1 + iadd;
			for (integer jstr = 1; jstr <= n_a[ilevel - 1]; jstr++) {
				// jstr - столбец матрицы А.
				// icounter-1 - число узлов на грубом уровне.
				//if ((icounter - 1 > n) || (icounter - 1 < 0)) {
				//printf("flag incorrupt 4...\n");
				//system("pause");
				//exit(1);
				//}
				for (integer i = 1; i <= icounter - 1; i++) {
					flag[i] = false;
				}
				// А (на первой позиции) должна быть отсортирована по j.
				integer ii1 = BinarySearchAj(A, jstr, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				integer lf = ii1;
				for (integer ii2 = ii1; (ii2 <= nnz_a[ilevel - 1] + iadd) && (A[ii2].j == jstr); ii2++) {
					lf++;
				}
				lf--;

				// это справедливо  для любого  уровня.
				//integer istart = 1 + iaddR;
				//integer iend = nnzR - 1 + iaddR;
				//for (integer ii = istart; ii <= iend; ii++) {

				//if (flag[R[ii].i] == false) {
				// Внимание не забыть или открыть флаг.
				//flag[R[ii].i] = true;

				for (integer i_2 = 1; i_2 <= numberofcoarcenodes; i_2++) {
					//integer istr = R[ii].i;
					//integer ic = ii;

					integer istr = i_2;
					integer ic = start_position_i_string_in_R[i_2];

					// i-coarse, j-fine
					Real sum1 = 0.0;

					// это медленная версия кода.
					//while ((ic <= iend) && (R[ic].i == istr)) {
					// R[R[ii].i][R[ic].j]*A[R[ic].j][jstr]

					//for (integer ii2 = ii1; (ii2 <= nnz_a[ilevel - 1] + iadd) && (A[ii2].j == jstr); ii2++) {

					//if (A[ii2].i == R[ic].j) {
					//sum1 += R[ic].aij*A[ii2].aij;
					//}

					//}
					//ic++;
					//}

					// Это более быстрый код.
					integer ks = ic;

					//integer kf = ic;

					//while ((kf <= iend) && (R[kf].i == istr)) {
					//kf++;
					//}
					//kf--;

					integer kf = kf_array[istr];

					integer ls = ii1;



					while ((ks <= kf) && (ls <= lf)) {

						if (A[ls].i < R[ks].j) {
							ls++;
						}
						else if (A[ls].i > R[ks].j) {
							ks++;
						}
						else /*if (A[ls].i == R[ks].j)*/ {
							sum1 += R[ks].aij*A[ls].aij;
							ks++;
							ls++;
						}

					}




					if (fabs(sum1) > 1.0e-30) {
						A[istartAnew].aij = sum1;
						A[istartAnew].i = istr;
						A[istartAnew].j = jstr;
						if (jstr < 0) {
							printf("fatal error");
							printf("i=%d j=%d aij=%e\n", istr, jstr, sum1);
							system("pause");
						}
						istartAnew++;
					}

					//}
				}
			}


			delete[] kf_array;
			delete[] start_position_i_string_in_R;

		}
		else if (0) {

			// Идея Писсанецки.

			//Более быстрый вариант алгоритма.
			// Быстрее этого кода на основе идеи слияния списков уже не будет.
			// 22 октября 2015. Нужно двигаться в сторону Писсанецки.
			// часть 1 : R*Afine.
			//         xxxxxx
			//         xxxxxx
			//  xxxxxx xxxxxx xxxxxx
			//  xxxxxx xxxxxx xxxxxx
			//         xxxxxx
			//         xxxxxx
			//    R       A     [RA]
			kf_array = new integer[numberofcoarcenodes + 1];
			integer istart1 = 1 + iaddR;
			integer iend1 = nnzR - 1 + iaddR;
			for (integer i = 1; i <= icounter - 1; i++) {
				flag[i] = false;
			}
			for (integer ii = istart1; ii <= iend1; ii++) {
				if (flag[R[ii].i] == false) {
					integer istr = R[ii].i;
					integer ic = ii;
					integer kf = ic;

					while ((kf <= iend1) && (R[kf].i == istr)) {
						kf++;
					}
					kf--;
					kf_array[istr] = kf;
					flag[R[ii].i] = true;
				}
			}

			integer *start_position_i_string_in_R = new integer[numberofcoarcenodes + 1];
			for (integer i = 1; i <= icounter - 1; i++) {
				flag[i] = false;
			}
			integer istart3 = 1 + iaddR;
			integer iend3 = nnzR - 1 + iaddR;
			for (integer ii = istart3; ii <= iend3; ii++) {
				if (flag[R[ii].i] == false) {
					start_position_i_string_in_R[R[ii].i] = ii;
					flag[R[ii].i] = true;
				}
			}

			integer* ind = new integer[n_a[ilevel - 1]+1];


			 istartAnew = nnz_a[ilevel - 1] + 1 + iadd;
			for (integer jstr = 1; jstr <= n_a[ilevel - 1]; jstr++) {
				// jstr - столбец матрицы А.
				// icounter-1 - число узлов на грубом уровне.
				//if ((icounter - 1 > n) || (icounter - 1 < 0)) {
				//printf("flag incorrupt 4...\n");
				//system("pause");
				//exit(1);
				//}

				for (integer i = 1; i <= n_a[ilevel - 1]; i++) {
					ind[i] = -1; // инициализация.
				}

				for (integer i = 1; i <= icounter - 1; i++) {
					flag[i] = false;
				}
				// А (на первой позиции) должна быть отсортирована по j.
				integer ii1 = BinarySearchAj(A, jstr, 1 + iadd, nnz_a[ilevel - 1] + iadd);
				//integer lf = ii1;
				for (integer ii2 = ii1; (ii2 <= nnz_a[ilevel - 1] + iadd) && (A[ii2].j == jstr); ii2++) {
					ind[A[ii2].i] = ii2; // запоминание индекса.
					//lf++;
				}
				//lf--;

				// это справедливо  для любого  уровня.
				//integer istart = 1 + iaddR;
				//integer iend = nnzR - 1 + iaddR;
				//for (integer ii = istart; ii <= iend; ii++) {

				//if (flag[R[ii].i] == false) {
				// Внимание не забыть или открыть флаг.
				//flag[R[ii].i] = true;

				for (integer i_2 = 1; i_2 <= numberofcoarcenodes; i_2++) {
					//integer istr = R[ii].i;
					//integer ic = ii;

					integer istr = i_2;
					integer ic = start_position_i_string_in_R[i_2];

					// i-coarse, j-fine
					Real sum1 = 0.0;

					

					// Это более быстрый код.
					integer ks = ic;

					
					integer kf = kf_array[istr];

					integer ls = ii1;


					while (ks <= kf) {
						if (ind[R[ks].j] != -1) {
							sum1 += R[ks].aij*A[ind[R[ks].j]].aij;
						}
						ks++;
					}

					/*
					while ((ks <= kf) && (ls <= lf)) {

						if (A[ls].i < R[ks].j) {
							ls++;
						}
						else if (A[ls].i > R[ks].j) {
							ks++;
						}
						else {
							sum1 += R[ks].aij*A[ls].aij;
							ks++;
							ls++;
						}

					}
				    */



					if (fabs(sum1) > 1.0e-30) {
						A[istartAnew].aij = sum1;
						A[istartAnew].i = istr;
						A[istartAnew].j = jstr;
						if (jstr < 0) {
							printf("fatal error");
							printf("i=%d j=%d aij=%e\n", istr, jstr, sum1);
							system("pause");
						}
						istartAnew++;
					}

					//}
				}
			}


			delete[] kf_array;
			delete[] start_position_i_string_in_R;
			delete[] ind;
		}
		else {

			if (0) {

				// Фред Густавсон IBM 1978.
				// 23 октября 2015 года.

				// часть 1 : R*Afine.
				//         xxxxxx
				//         xxxxxx
				//  xxxxxx xxxxxx xxxxxx
				//  xxxxxx xxxxxx xxxxxx
				//         xxxxxx
				//         xxxxxx
				//    R       A     [RA]
				// Сортировка А по строкам.
				HeapSort(A, n_a[ilevel - 1], 1 + iadd, nnz_a[ilevel - 1] + iadd);

				// Преобразование к формату CRS.

				integer* row_ind_SR = new integer[numberofcoarcenodes + 1];
				integer* row_ind_ER = new integer[numberofcoarcenodes + 1];
				for (integer i_1 = 1; i_1 <= numberofcoarcenodes; i_1++) {
					row_ind_SR[i_1] = -1;
					row_ind_ER[i_1] = -2;
				}
				integer istart1 = 1 + iaddR;
				integer iend1 = nnzR - 1 + iaddR;
				for (integer i = 1; i <= icounter - 1; i++) {
					flag[i] = false;
				}
				for (integer ii = istart1; ii <= iend1; ii++) {
					if (flag[R[ii].i] == false) {
						integer istr = R[ii].i;
						integer ic = ii;
						integer kf = ic;

						while ((kf <= iend1) && (R[kf].i == istr)) {
							kf++;
						}
						kf--;
						row_ind_SR[istr] = ic;
						row_ind_ER[istr] = kf;
						flag[R[ii].i] = true;
					}
				}

				// Пустые строки просто отсутствуют.
				//for (integer i_1 = 1; i_1 <= numberofcoarcenodes; i_1++) {
				//if (row_ind_SR[i_1] == -1) {
				//printf("empty string %d\n", row_ind_ER[i_1]);
				//}
				//}

				integer* row_ind_SA = new integer[n_a[ilevel - 1] + 1];
				integer* row_ind_EA = new integer[n_a[ilevel - 1] + 1];
				integer istart3 = 1 + iadd;
				integer iend3 = nnz_a[ilevel - 1] + iadd;
				for (integer i = 1; i <= n_a[ilevel - 1]; i++) {
					flag[i] = false;
				}
				for (integer ii = istart3; ii <= iend3; ii++) {
					if (flag[A[ii].i] == false) {
						integer istr = A[ii].i;
						integer ic = ii;
						integer kf = ic;

						while ((kf <= iend3) && (A[kf].i == istr)) {
							kf++;
						}
						kf--;
						row_ind_SA[istr] = ic;
						row_ind_EA[istr] = kf;
						flag[A[ii].i] = true;
					}
				}

				istartAnew = nnz_a[ilevel - 1] + 1 + iadd;

				// Данные используемые для частичного формирователя суммы.
				Real* vector_sum = new Real[n_a[ilevel - 1] + 1];
				//bool* b_visit_vec_sum = new bool[n_a[ilevel - 1] + 1];
				integer size_v = sizeof(Real)*(1 + n_a[ilevel - 1]);

				// Сканируем первый операнд построчно.
				for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

					// Начинаем обрабатывать новую строку.
					// Сброс формирователя суммы в ноль.
					//#pragma omp parallel for
					for (integer i_2 = 1; i_2 <= n_a[ilevel - 1]; i_2++) {
						vector_sum[i_2] = 0.0; // инициализация // 18.5
						//b_visit_vec_sum[i_2] = false;
					}
					// Более быстрое обнуление
					// По времени в лучшем случае также. Но может сильно ухудшиться из-за занятости системной dll-ки.
					//memset(vector_sum, 0,size_v );

					// Сканируем текущую i-ую строку поэлементно
					for (integer ii = row_ind_SR[istr]; ii <= row_ind_ER[istr]; ii++) {
						integer col_ind = R[ii].j;
						// Сканируем col_ind строку второго операнда

						for (integer i_1 = row_ind_SA[col_ind]; i_1 <= row_ind_EA[col_ind]; i_1++) {
							Real left_operand = R[ii].aij;
							Real right_operand = A[i_1].aij;
							vector_sum[A[i_1].j] += left_operand*right_operand;
							//b_visit_vec_sum[A[i_1].j] = true;
						}
					}

					for (integer jstr = 1; jstr <= n_a[ilevel - 1]; jstr++) {
						// сначала быстрая проверка 
						//if (b_visit_vec_sum[jstr]) {
						// а потом медленная.
						if (fabs(vector_sum[jstr]) > 1.0e-30) { // 36.3
							A[istartAnew].aij = vector_sum[jstr];
							A[istartAnew].i = istr;
							A[istartAnew].j = jstr;
							istartAnew++;
						}
						//}
					}
				}

				delete[] row_ind_SR;
				delete[] row_ind_ER;
				delete[] row_ind_SA;
				delete[] row_ind_EA;
				delete[] vector_sum;
				//delete[] b_visit_vec_sum;

			}
			else {
				// Фред Густавсон IBM 1978.
				// 23 октября 2015 года.

				// часть 1 : R*Afine.
				//         xxxxxx
				//         xxxxxx
				//  xxxxxx xxxxxx xxxxxx
				//  xxxxxx xxxxxx xxxxxx
				//         xxxxxx
				//         xxxxxx
				//    R       A     [RA]
				// Сортировка А по строкам.
				HeapSort(A, n_a[ilevel - 1], 1 + iadd, nnz_a[ilevel - 1] + iadd);

				// Преобразование к формату CRS.

				integer* row_ind_SR = new integer[numberofcoarcenodes + 1];
				integer* row_ind_ER = new integer[numberofcoarcenodes + 1];
				for (integer i_1 = 1; i_1 <= numberofcoarcenodes; i_1++) {
					row_ind_SR[i_1] = -1;
					row_ind_ER[i_1] = -2;
				}
				integer istart1 = 1 + iaddR;
				integer iend1 = nnzR - 1 + iaddR;
				for (integer i = 1; i <= icounter - 1; i++) {
					flag[i] = false;
				}
				for (integer ii = istart1; ii <= iend1; ii++) {
					if (flag[R[ii].i] == false) {
						integer istr = R[ii].i;
						integer ic = ii;
						integer kf = ic;

						while ((kf <= iend1) && (R[kf].i == istr)) {
							kf++;
						}
						kf--;
						row_ind_SR[istr] = ic;
						row_ind_ER[istr] = kf;
						flag[R[ii].i] = true;
					}
				}

				// Пустые строки просто отсутствуют.
				//for (integer i_1 = 1; i_1 <= numberofcoarcenodes; i_1++) {
				//if (row_ind_SR[i_1] == -1) {
				//printf("empty string %d\n", row_ind_ER[i_1]);
				//}
				//}

				integer* row_ind_SA = new integer[n_a[ilevel - 1] + 1];
				integer* row_ind_EA = new integer[n_a[ilevel - 1] + 1];
				integer istart3 = 1 + iadd;
				integer iend3 = nnz_a[ilevel - 1] + iadd;
				for (integer i = 1; i <= n_a[ilevel - 1]; i++) {
					flag[i] = false;
				}
				for (integer ii = istart3; ii <= iend3; ii++) {
					if (flag[A[ii].i] == false) {
						integer istr = A[ii].i;
						integer ic = ii;
						integer kf = ic;

						while ((kf <= iend3) && (A[kf].i == istr)) {
							kf++;
						}
						kf--;
						row_ind_SA[istr] = ic;
						row_ind_EA[istr] = kf;
						flag[A[ii].i] = true;
					}
				}

				istartAnew = nnz_a[ilevel - 1] + 1 + iadd;

				// Данные используемые для частичного формирователя суммы.
				Real* vector_sum = new Real[n_a[ilevel - 1] + 1];
				//bool* b_visit_vec_sum = new bool[n_a[ilevel - 1] + 1];
				integer size_v = sizeof(Real)*(1 + n_a[ilevel - 1]);
				// Храним индексы ненулевых элементов в отсортированном порядке.
				integer* index_visit = new integer[n_a[ilevel - 1] + 1];
				integer index_size = 0;

				// Сканируем первый операнд построчно.
				for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

					// Начинаем обрабатывать новую строку.
					// Сброс формирователя суммы в ноль.
					//#pragma omp parallel for
					//for (integer i_2 = 1; i_2 <= n_a[ilevel - 1]; i_2++) {
					//vector_sum[i_2] = 0.0; // инициализация // 18.5
					//b_visit_vec_sum[i_2] = false;
					//}
					// Более быстрое обнуление
					// По времени в лучшем случае также. Но может сильно ухудшиться из-за занятости системной dll-ки.
					//memset(vector_sum, 0,size_v );

					// Сканируем текущую i-ую строку поэлементно
					for (integer ii = row_ind_SR[istr]; ii <= row_ind_ER[istr]; ii++) {
						integer col_ind = R[ii].j;
						// Сканируем col_ind строку второго операнда

						for (integer i_1 = row_ind_SA[col_ind]; i_1 <= row_ind_EA[col_ind]; i_1++) {
							Real left_operand = R[ii].aij;
							Real right_operand = A[i_1].aij;
							integer iaddind = A[i_1].j;
							bool foundnow = false;
							integer ifoundind = -1;
							// линейный поиск позиции в массиве на добавление.
							for (integer i_6 = 1; i_6 <= index_size; i_6++) {
								if (index_visit[i_6] == iaddind) {
									foundnow = true;
									ifoundind = i_6;
									break;
								}
							}
							if (foundnow) {
								vector_sum[index_visit[ifoundind]] += left_operand*right_operand;
							}
							else {
								// Первое добавление.
								index_size++;
								index_visit[index_size] = iaddind;
								ifoundind = index_size;
								vector_sum[index_visit[ifoundind]] = left_operand*right_operand;
							}

							//vector_sum[A[i_1].j] += left_operand*right_operand;
							//b_visit_vec_sum[A[i_1].j] = true;
						}
					}

					for (integer i_6 = 1; i_6 <= index_size; i_6++) {
						integer jstr = index_visit[i_6];
						if (fabs(vector_sum[jstr]) > 1.0e-30) {
							A[istartAnew].aij = vector_sum[jstr];
							A[istartAnew].i = istr;
							A[istartAnew].j = jstr;
							istartAnew++;
						}
					}
					index_size = 0;


				}
				delete[] index_visit;
				delete[] row_ind_SR;
				delete[] row_ind_ER;
				delete[] row_ind_SA;
				delete[] row_ind_EA;
				delete[] vector_sum;
			}

		}

		// Сортировка А по i (оригинала.).
		//heapsort(A, key=i*n_a[ilevel - 1] + j, 1, nnz_a[ilevel - 1]);
		if (bquicktypesort) {
			QuickSort(A, 1 + iadd, nnz_a[ilevel - 1] + iadd);
		}
		else {
			if (nnz_a[ilevel - 1] < heapsortsizelimit) {
				HeapSort(A, n_a[ilevel - 1], 1 + iadd, nnz_a[ilevel - 1] + iadd);
			}
			else {
				Ak* Aorig = &A[1 + iadd];
				MergeSort(Aorig, nnz_a[ilevel - 1]);
				Aorig = NULL;
			}
		}

        // Проверка пройдена успешно.
		for (integer i_1 = nnz_a[ilevel - 1] + 1 + iadd; i_1 <= istartAnew - 1; i_1++) {
			if (A[i_1].j < 0) {
				printf("error : negativ j index\n");
				system("pause");
			}
			//printf("A[%d].i=%d A[%d].j=%d A[%d].aij=%e\n", i_1, A[i_1].i, i_1, A[i_1].j, i_1, A[i_1].aij);
			//system("pause");
		}


		
		// Часть 2. [R*Afine]*P=Abuf*P.
		// Сортировка [R*А] по i.
		//heapsort(A, key=i*n_coarce + j, 1, nnz_a[ilevel - 1]);
		if (bquicktypesort) {
			QuickSort(A, nnz_a[ilevel - 1] + 1 + iadd, istartAnew - 1);
		}
		else {
			if (nnz_a[ilevel - 1] < heapsortsizelimit) {
		        HeapSort(A, icounter - 1, nnz_a[ilevel - 1] + 1 + iadd, istartAnew - 1);
			}
			else {
				Ak* Aorig = &A[nnz_a[ilevel - 1] + 1 + iadd];
				MergeSort(Aorig, istartAnew -1-( nnz_a[ilevel - 1]+1+iadd)+1);
				Aorig = NULL;
			}
		}


		// Проверка пройдена успешно.
		for (integer i_1 = nnz_a[ilevel - 1] + 1 + iadd; i_1 <= istartAnew - 1; i_1++) {
			if (A[i_1].j < 0) {
				printf("error : negativ j index\n");
				system("pause");
			}
		//-->	//printf("A[%d].i=%d A[%d].j=%d A[%d].aij=%e\n", i_1, A[i_1].i, i_1, A[i_1].j, i_1, A[i_1].aij);
			//system("pause");
		}

		// Самый медленный в мире BubbleSort.
		// нужен для проверки более быстрых алгоритмов сортировки.
		//integer k3 = 0;
		//for (integer k1 = nnz_a[ilevel - 1] + 1 + iadd; k1 < istartAnew - 1; k1++,k3++) {
			//for (integer k2 = nnz_a[ilevel - 1] + 1 + iadd; k2 < istartAnew - 1 - k3; k2++) {
				//if (A[k2].i>A[k2 + 1].i) {
					// change
					//Ak Temp = A[k2];
				//	A[k2] = A[k2 + 1];
					//A[k2 + 1] = Temp;
				//}
				//else if (A[k2].i == A[k2 + 1].i) {
					//if (A[k2].j>A[k2 + 1].j) {
						// change
						//Ak Temp = A[k2];
						//A[k2] = A[k2 + 1];
						//A[k2 + 1] = Temp;
					//}
				//}
			//}
		//}

		// Контроль [R*A] debug
		//for (int i_1 = nnz_a[ilevel - 1] + 1 + iadd; i_1 <= istartAnew - 1; i_1++) {
			//printf("A[%d].i=%d A[%d].j=%d A[%d].aij=%e\n", i_1, A[i_1].i, i_1, A[i_1].j, i_1, A[i_1].aij);
			//system("pause");
		//}

		if (bquicktypesort) {
			QuickSort(P, 1 + iaddR, iaddR + nnzR - 1);
		}
		else {
			HeapSort(P, n_a[ilevel - 1], 1 + iaddR, iaddR + nnzR - 1);
		}

		// Prolongation должна быть упорядочена по j.
		// Начальная позиция элементов матрицы грубосеточного уровня.
		integer istartAnew2 = istartAnew; 

		/*
		// icounter - 1; // Количество узлов на грубосеточном уровне.
		for (integer i = 1; i <= n_coarce; i++) {
			flag[i] = false;
		}
		// P должен быть упорядочен по строкам.
		// проверка диагонали в P
		for (integer ii77 = 1; ii77 <= nnzR; ii77++) {
			printf("i=%d j=%d %e\n",P[ii77].i,P[ii77].j,P[ii77].aij);
			getchar();
		}
		*/

		//begin medlennji
		/*
		// сканируем по столбцам
		for (integer jstr = 1; jstr <= icounter-1; jstr++) {
			for (integer i = 1; i <= n_coarce; i++) {
				flag[i] = false;
			}
			//integer ii1 = BinarySearchAj(P, jstr, 1+iaddR, nnzR - 1+iaddR);
			// Дело в том что у матрицы P всё хранится в перепутанном (с точностью до транспонирования виде).
			integer ii1 = BinarySearchAi(P, jstr, 1 + iaddR, nnzR - 1 + iaddR);
			// это справедливо только для первого уровня.
			integer istart = nnz_a[ilevel - 1] + 1 + iadd;
			integer iend = istartAnew - 1;
			for (integer ii = istart; ii <= iend; ii++) {
				if (flag[A[ii].i] == false) {
					// сканируем построчно.
					integer istr = A[ii].i;
					integer ic = ii;
					// i-coarse, j-fine
					Real sum1 = 0.0;
					while ((ic <= iend) && (A[ic].i == istr)) {
						// [R*A][A[ii].i][A[ic].j]*P[A[ic].j][jstr]
						// у матрицы P всё перепутано с точностью до транспонирования.
						// [R*A][A[ii].i][A[ic].j]*P[jstr][A[ic].j]
						for (integer ii2 = ii1; (ii2 <= nnzR- 1+iaddR) && (P[ii2].i == jstr); ii2++) {
							
							if (P[ii2].j == A[ic].j) {
								sum1 += A[ic].aij*P[ii2].aij;
								//printf("%e i=%d j=%d k=%d\n", A[ic].aij*P[ii2].aij,istr,P[ii2].j,jstr);
								//getchar();
							}
						}

						ic++;
					}

					if (fabs(sum1) > RealZERO) {
						A[istartAnew2].aij = sum1;
						A[istartAnew2].i = istr;
						A[istartAnew2].j = jstr;
						istartAnew2++;
				    }
					flag[A[ii].i] = true;
				}
			}
		}
		*/



		
		// сканируем по столбцам
		//for (integer jstr = 1; jstr <= icounter - 1; jstr++) {
			
			//for (integer i = 1; i <= n; i++) {
				//flag[i] = false;
			//}
			//integer ii1 = BinarySearchAj(P, jstr, 1+iaddR, nnzR - 1+iaddR);
			// Дело в том что у матрицы P всё хранится в перепутанном (с точностью до транспонирования виде).
			//integer ii1 = BinarySearchAi(P, jstr, 1 + iaddR, nnzR - 1 + iaddR);
			// это справедливо только для первого уровня.
			//integer istart = nnz_a[ilevel - 1] + 1 + iadd;
			//integer iend = istartAnew - 1;
			//for (integer ii = istart; ii <= iend; ii++) {
				//if ((A[ii].i > n)||(A[ii].i<0)) {
					//printf("flag incorrupt 1...\n");
					//system("pause");
					//exit(1);
				//}
				//if (flag[A[ii].i] == false) {
					// сканируем построчно.
					//integer istr = A[ii].i;
					//integer ic = ii;
					// i-coarse, j-fine
				//	Real sum1 = 0.0;

					// Заменим медленный двойной цикл на слияние.
					//while ((ic <= iend) && (A[ic].i == istr)) {
						// [R*A][A[ii].i][A[ic].j]*P[A[ic].j][jstr]
						// у матрицы P всё перепутано с точностью до транспонирования.
						// [R*A][A[ii].i][A[ic].j]*P[jstr][A[ic].j]
						//for (integer ii2 = ii1; (ii2 <= nnzR - 1 + iaddR) && (P[ii2].i == jstr); ii2++) {

						//	if (P[ii2].j == A[ic].j) {
								//sum1 += A[ic].aij*P[ii2].aij;
						//	}
						//}

						//ic++;
					//}

					// Более быстрый код на основе слияния.
					//integer ks = ic;
					//integer ls = ii1;
					//integer kf = ic;
					//Real retalon = 0.0;
					//bool bvis = false;
					//while ((kf <= iend) && (A[kf].i == istr)) {
						//if (fabs(A[kf].aij) > retalon) retalon = fabs(A[kf].aij);
						//if (A[kf].j == istr) {
						//	retalon = fabs(A[kf].aij);
						//	bvis = true;
						//}
						///kf++;
					//}
					//kf--;
					//if (bvis == false) {
						//kf = ic;
						//while ((kf <= iend) && (A[kf].i == istr)) {
						//	if (fabs(A[kf].aij) > retalon) retalon = fabs(A[kf].aij);
						//	kf++;
						//}
						//kf--;
					//}
					//integer lf = ii1;
					//for (integer ii2 = ii1; (ii2 <= nnzR - 1 + iaddR) && (P[ii2].i == jstr); ii2++) {
						//if (P[ii2].j == istr) retalon *= fabs(P[ii2].aij); // это диагональный элемент.
						//lf++;
					//}
					//lf--;
					//while ((ks <= kf) && (ls <= lf)) {
						//if (P[ls].j<A[ks].j) {
						//	ls++;
						//}
						//else if (P[ls].j > A[ks].j) {
						//	ks++;
						//}
						//else {
						//	sum1 += A[ks].aij*P[ls].aij;
						//	ks++;
						//	ls++;
						//}
					//}
					//if (fabs(retalon) < 1.0e-30) {
						//printf("RAP retalon=%e string %d is zero\n",retalon,istr);
						//getchar();
					//}

					//---//if (fabs(sum1) > 0.001*retalon) {
					//if (fabs(sum1)>1.0e-30) {
						//A[istartAnew2].aij = sum1;
						//A[istartAnew2].i = istr;
						//A[istartAnew2].j = jstr;
						//istartAnew2++;
					//}
					//flag[A[ii].i] = true;
				//}
			//}
		//}

// Быстрее этого кода на основе идеи слияния списков уже не будет.
// 17 октября 2015. Нужно двигаться в сторону Писсанецки.
printf("nnz left operand=%d, nnz right operand=%d\n", istartAnew - (nnz_a[ilevel - 1] + 1 + iadd), nnzR);

if (0) {

	// Код на основе слияния упорядоченных массивов.

	kf_array = new integer[numberofcoarcenodes + 1];
	integer istart2 = nnz_a[ilevel - 1] + 1 + iadd;
	integer iend2 = istartAnew - 1;
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	for (integer ii = istart2; ii <= iend2; ii++) {
		if (flag[A[ii].i] == false) {
			// сканируем построчно.
			integer istr = A[ii].i;
			integer ic = ii;

			integer kf = ic;

			while ((kf <= iend2) && (A[kf].i == istr)) {
				kf++;
			}
			kf--;
			kf_array[istr] = kf;
			flag[A[ii].i] = true;

		}
	}

	// Количество строк в матрице А есть точно numberofcoarcenodes
	integer *start_position_i_string_in_RA = new integer[numberofcoarcenodes + 1];
	for (integer i = 1; i <= icounter - 1; i++) {
		flag[i] = false;
	}
	integer istart4 = nnz_a[ilevel - 1] + 1 + iadd;
	integer iend4 = istartAnew - 1;
	for (integer ii = istart4; ii <= iend4; ii++) {
		if (flag[A[ii].i] == false) {
			start_position_i_string_in_RA[A[ii].i] = ii;
			flag[A[ii].i] = true;
		}
	}


	// Более быстрая версия кода : 15 октября 2015
	// сканируем по столбцам
	for (integer jstr = 1; jstr <= icounter - 1; jstr++) {

		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		//integer ii1 = BinarySearchAj(P, jstr, 1+iaddR, nnzR - 1+iaddR);
		// Дело в том что у матрицы P всё хранится в перепутанном (с точностью до транспонирования виде).
		integer ii1 = BinarySearchAi(P, jstr, 1 + iaddR, nnzR - 1 + iaddR);

		integer lf = ii1;
		for (integer ii2 = ii1; (ii2 <= nnzR - 1 + iaddR) && (P[ii2].i == jstr); ii2++) {
			lf++;
		}
		lf--;

		// это справедливо только для первого уровня.
		//integer istart = nnz_a[ilevel - 1] + 1 + iadd;
		//integer iend = istartAnew - 1;
		//for (integer ii = istart; ii <= iend; ii++) {

		//if (flag[A[ii].i] == false) {
		// Ни в коем случае не забываем модифицировать флаг.
		//flag[A[ii].i] = true;
		// сканируем построчно.
		//integer istr = A[ii].i;
		//integer ic = ii;

		for (integer i_2 = 1; i_2 <= numberofcoarcenodes; i_2++) {

			integer istr = i_2;
			integer ic = start_position_i_string_in_RA[i_2];

			// i-coarse, j-fine
			Real sum1 = 0.0;

			// Заменим медленный двойной цикл на слияние.
			//while ((ic <= iend) && (A[ic].i == istr)) {
			// [R*A][A[ii].i][A[ic].j]*P[A[ic].j][jstr]
			// у матрицы P всё перепутано с точностью до транспонирования.
			// [R*A][A[ii].i][A[ic].j]*P[jstr][A[ic].j]
			//for (integer ii2 = ii1; (ii2 <= nnzR - 1 + iaddR) && (P[ii2].i == jstr); ii2++) {

			//	if (P[ii2].j == A[ic].j) {
			//sum1 += A[ic].aij*P[ii2].aij;
			//	}
			//}

			//ic++;
			//}

			// Более быстрый код на основе слияния.
			integer ks = ic;

			//integer kf = ic;

			//while ((kf <= iend) && (A[kf].i == istr)) {
			//kf++;
			//}
			//kf--;
			integer kf = kf_array[istr];

			integer ls = ii1;

			while ((ks <= kf) && (ls <= lf)) {

				if (P[ls].j < A[ks].j) {
					ls++;
				}
				else if (P[ls].j > A[ks].j) {
					ks++;
				}
				else /*if (P[ls].j==A[ks].j)*/ {
					sum1 += A[ks].aij*P[ls].aij;
					ks++;
					ls++;
				}

			}

			if (fabs(sum1) > 1.0e-30) {
				A[istartAnew2].aij = sum1;
				A[istartAnew2].i = istr;
				A[istartAnew2].j = jstr;
				istartAnew2++;
			}

			//}
		}
	}

	delete[] kf_array;
	delete[] start_position_i_string_in_RA;

}
else if (0) {
	// Идея Писсанецки 22 октября 2015.

	kf_array = new integer[numberofcoarcenodes + 1];
	integer istart2 = nnz_a[ilevel - 1] + 1 + iadd;
	integer iend2 = istartAnew - 1;
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	for (integer ii = istart2; ii <= iend2; ii++) {
		if (flag[A[ii].i] == false) {
			// сканируем построчно.
			integer istr = A[ii].i;
			integer ic = ii;

			integer kf = ic;

			while ((kf <= iend2) && (A[kf].i == istr)) {
				kf++;
			}
			kf--;
			kf_array[istr] = kf;
			flag[A[ii].i] = true;

		}
	}

	// Количество строк в матрице А есть точно numberofcoarcenodes
	integer *start_position_i_string_in_RA = new integer[numberofcoarcenodes + 1];
	for (integer i = 1; i <= icounter - 1; i++) {
		flag[i] = false;
	}
	integer istart4 = nnz_a[ilevel - 1] + 1 + iadd;
	integer iend4 = istartAnew - 1;
	for (integer ii = istart4; ii <= iend4; ii++) {
		if (flag[A[ii].i] == false) {
			start_position_i_string_in_RA[A[ii].i] = ii;
			flag[A[ii].i] = true;
		}
	}

	integer *ind = new integer[n_a[ilevel-1] + 1];


	// Более быстрая версия кода : 15 октября 2015
	// сканируем по столбцам
	for (integer jstr = 1; jstr <= icounter - 1; jstr++) {

		for (integer i = 1; i <= n_a[ilevel-1]; i++) {
			ind[i] = -1; // инициализация.
		}

		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		//integer ii1 = BinarySearchAj(P, jstr, 1+iaddR, nnzR - 1+iaddR);
		// Дело в том что у матрицы P всё хранится в перепутанном (с точностью до транспонирования виде).
		integer ii1 = BinarySearchAi(P, jstr, 1 + iaddR, nnzR - 1 + iaddR);

		//integer lf = ii1;
		for (integer ii2 = ii1; (ii2 <= nnzR - 1 + iaddR) && (P[ii2].i == jstr); ii2++) {
			ind[P[ii2].j] = ii2;
			//lf++;
		}
		//lf--;

		// это справедливо только для первого уровня.
		//integer istart = nnz_a[ilevel - 1] + 1 + iadd;
		//integer iend = istartAnew - 1;
		//for (integer ii = istart; ii <= iend; ii++) {

		//if (flag[A[ii].i] == false) {
		// Ни в коем случае не забываем модифицировать флаг.
		//flag[A[ii].i] = true;
		// сканируем построчно.
		//integer istr = A[ii].i;
		//integer ic = ii;

		for (integer i_2 = 1; i_2 <= numberofcoarcenodes; i_2++) {

			integer istr = i_2;
			integer ic = start_position_i_string_in_RA[i_2];

			// i-coarse, j-fine
			Real sum1 = 0.0;

			

			// Более быстрый код на основе слияния.
			integer ks = ic;

			
			integer kf = kf_array[istr];

			//integer ls = ii1;

			while (ks <= kf) {
				if (ind[A[ks].j] != -1) {
					sum1 += A[ks].aij*P[ind[A[ks].j]].aij;
				}
				ks++;
			}

			//while ((ks <= kf) && (ls <= lf)) {

				//if (P[ls].j < A[ks].j) {
					//ls++;
				//}
				//else if (P[ls].j > A[ks].j) {
					//ks++;
				//}
				//else /*if (P[ls].j==A[ks].j)*/ {
					//sum1 += A[ks].aij*P[ls].aij;
					//ks++;
					//ls++;
				//}

			//}

			if (fabs(sum1) > 1.0e-30) {
				A[istartAnew2].aij = sum1;
				A[istartAnew2].i = istr;
				A[istartAnew2].j = jstr;
				istartAnew2++;
			}

			//}
		}
	}

	delete[] kf_array;
	delete[] start_position_i_string_in_RA;
	delete[] ind;

}
else {
	// Фред Густавсон IBM 1978
	// В ядре кода Густавсона нету ни одного ветвления,
	// а мы знаем что в результате профайлинга предыдущих версий кода :
	// (наивный, слияние, Писсанецки) львиная доля вычислительной работы уходила
	// на сравнения (ветвления) в отношении примерно 30 к 1. 30 сравнений на одно сумирование.
	// 23 октября 2015 года.

	if (0) {
		// Преобразование обоих матриц в формат CRS.
		HeapSort_j(P, 1 + iaddR, iaddR + nnzR - 1);

		integer* row_ind_AS = new integer[numberofcoarcenodes + 1];
		integer* row_ind_AE = new integer[numberofcoarcenodes + 1];
		integer istart2 = nnz_a[ilevel - 1] + 1 + iadd;
		integer iend2 = istartAnew - 1;
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		for (integer ii = istart2; ii <= iend2; ii++) {
			if (flag[A[ii].i] == false) {
				// сканируем построчно.
				integer istr = A[ii].i;
				integer ic = ii;

				integer kf = ic;

				while ((kf <= iend2) && (A[kf].i == istr)) {
					kf++;
				}
				kf--;
				row_ind_AS[istr] = ic;
				row_ind_AE[istr] = kf;
				flag[A[ii].i] = true;

			}
		}

		integer* row_ind_PS = new integer[n_a[ilevel - 1] + 1];
		integer* row_ind_PE = new integer[n_a[ilevel - 1] + 1];
		// Инициализация чрезвычайно важна, т.к. 
		// обязательно присутствуют пустые строки которые
		// надо корректно обрабатывать.
		for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
			row_ind_PS[ii] = -1; // инициализация.
			row_ind_PE[ii] = -2;
		}
		integer istart4 = 1 + iaddR;
		integer iend4 = nnzR - 1 + iaddR;
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		for (integer ii = istart4; ii <= iend4; ii++) {
			if (flag[P[ii].j] == false) {
				// сканируем построчно.
				integer istr = P[ii].j;
				integer ic = ii;

				integer kf = ic;

				while ((kf <= iend4) && (P[kf].j == istr)) {
					kf++;
				}
				kf--;
				row_ind_PS[istr] = ic;
				row_ind_PE[istr] = kf;
				flag[P[ii].j] = true;

			}
		}

		// Данный код подтверждает что обязательно присутствуют
		// пустые строки.
		//for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
		//if (row_ind_PS[ii] == -1) {
		//printf("propusk string %d\n", row_ind_PE[ii]);
		//getchar();
		//}
		//}

		// Накопитель результата.
		Real* vector_sum = new Real[numberofcoarcenodes + 1];
		integer size_v = sizeof(Real)*(1 + numberofcoarcenodes);

		// Мы будем сканировать левый операнд построчно, а
		// после окончания обработки одной строки левого операнда
		// получать готовую строку результата.

		for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

			// Переход к новой строке, сброс накопителя результата.
			//#pragma omp parallel for
			for (integer i_1 = 1; i_1 <= numberofcoarcenodes; i_1++) {
				vector_sum[i_1] = 0.0; // обнуление.
			}
			// хоть профайлер и показывает что на memset тратиться меньше процессорного времени,
			// но там возникает системный вызов на который тратиться столько же времени что и раньше.
			// работа перекочевала из одного места в другое.
			// Более быстрое обнуление.
			//memset(vector_sum, 0, size_v);

			// сканируем все элементы строки левого операнда.
			for (integer ii1 = row_ind_AS[istr]; ii1 <= row_ind_AE[istr]; ii1++) {
				integer col_ind = A[ii1].j;

				// Сканируем col_ind строку правого операнда накапливая сумму.
				for (integer ii2 = row_ind_PS[col_ind]; ii2 <= row_ind_PE[col_ind]; ii2++) {
					Real rleft = A[ii1].aij;
					Real rright = P[ii2].aij;

					vector_sum[P[ii2].i] += rleft*rright;
				}
			}

			for (integer jstr = 1; jstr <= numberofcoarcenodes; jstr++) {
				if (fabs(vector_sum[jstr]) > 1.0e-30) {
					A[istartAnew2].aij = vector_sum[jstr];
					A[istartAnew2].i = istr;
					A[istartAnew2].j = jstr;
					istartAnew2++;
				}
			}

		}

		delete[] vector_sum;



		delete[] row_ind_AS;
		delete[] row_ind_AE;
		delete[] row_ind_PS;
		delete[] row_ind_PE;
	}
	else {

		// Рабочая версия алгоритма Фреда Густавсона.
		// IBM 1978 Sparse Matrix multiplication.

		// Преобразование обоих матриц в формат CRS.
		HeapSort_j(P, 1 + iaddR, iaddR + nnzR - 1);

		integer* row_ind_AS = new integer[numberofcoarcenodes + 1];
		integer* row_ind_AE = new integer[numberofcoarcenodes + 1];
		integer istart2 = nnz_a[ilevel - 1] + 1 + iadd;
		integer iend2 = istartAnew - 1;
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		for (integer ii = istart2; ii <= iend2; ii++) {
			if (flag[A[ii].i] == false) {
				// сканируем построчно.
				integer istr = A[ii].i;
				integer ic = ii;

				integer kf = ic;

				while ((kf <= iend2) && (A[kf].i == istr)) {
					kf++;
				}
				kf--;
				row_ind_AS[istr] = ic;
				row_ind_AE[istr] = kf;
				flag[A[ii].i] = true;

			}
		}

		integer* row_ind_PS = new integer[n_a[ilevel - 1] + 1];
		integer* row_ind_PE = new integer[n_a[ilevel - 1] + 1];
		// Инициализация чрезвычайно важна, т.к. 
		// обязательно присутствуют пустые строки которые
		// надо корректно обрабатывать.
		for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
			row_ind_PS[ii] = -1; // инициализация.
			row_ind_PE[ii] = -2;
		}
		integer istart4 = 1 + iaddR;
		integer iend4 = nnzR - 1 + iaddR;
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		for (integer ii = istart4; ii <= iend4; ii++) {
			if (flag[P[ii].j] == false) {
				// сканируем построчно.
				integer istr = P[ii].j;
				integer ic = ii;

				integer kf = ic;

				while ((kf <= iend4) && (P[kf].j == istr)) {
					kf++;
				}
				kf--;
				row_ind_PS[istr] = ic;
				row_ind_PE[istr] = kf;
				flag[P[ii].j] = true;

			}
		}

		// Данный код подтверждает что обязательно присутствуют
		// пустые строки.
		//for (integer ii = 1; ii <= n_a[ilevel - 1]; ii++) {
		//if (row_ind_PS[ii] == -1) {
		//printf("propusk string %d\n", row_ind_PE[ii]);
		//getchar();
		//}
		//}

		// Накопитель результата.
		Real* vector_sum = new Real[numberofcoarcenodes + 1];
		integer size_v = sizeof(Real)*(1 + numberofcoarcenodes);
		// Храним индексы ненулевых элементов в отсортированном порядке.
		integer* index_visit = new integer[n_a[ilevel - 1] - 1];
		integer index_size = 0;

		// Мы будем сканировать левый операнд построчно, а
		// после окончания обработки одной строки левого операнда
		// получать готовую строку результата.

		for (integer istr = 1; istr <= numberofcoarcenodes; istr++) {

			// Переход к новой строке, сброс накопителя результата.
			//#pragma omp parallel for
			//for (integer i_1 = 1; i_1 <= numberofcoarcenodes; i_1++) {
				//vector_sum[i_1] = 0.0; // обнуление.
			//}
			// хоть профайлер и показывает что на memset тратиться меньше процессорного времени,
			// но там возникает системный вызов на который тратиться столько же времени что и раньше.
			// работа перекочевала из одного места в другое.
			// Более быстрое обнуление.
			//memset(vector_sum, 0, size_v);

			// сканируем все элементы строки левого операнда.
			for (integer ii1 = row_ind_AS[istr]; ii1 <= row_ind_AE[istr]; ii1++) {
				integer col_ind = A[ii1].j;

				// Сканируем col_ind строку правого операнда накапливая сумму.
				for (integer ii2 = row_ind_PS[col_ind]; ii2 <= row_ind_PE[col_ind]; ii2++) {
					Real left_operand = A[ii1].aij;
					Real right_operand = P[ii2].aij;

					integer iaddind = P[ii2].i;
					bool foundnow = false;
					integer ifoundind = -1;
					// линейный поиск позиции в массиве на добавление.
					for (integer i_6 = 1; i_6 <= index_size; i_6++) {
						if (index_visit[i_6] == iaddind) {
							foundnow = true;
							ifoundind = i_6;
							break;
						}
					}
					if (foundnow) {
						vector_sum[index_visit[ifoundind]] += left_operand*right_operand;
					}
					else {
						// Первое добавление.
						index_size++;
						index_visit[index_size] = iaddind;
						ifoundind = index_size;
						vector_sum[index_visit[ifoundind]] = left_operand*right_operand;
					}



					//vector_sum[P[ii2].i] += rleft*rright;
				}
			}

			
			for (integer i_6 = 1; i_6 <= index_size; i_6++) {
				integer jstr = index_visit[i_6];
				if (fabs(vector_sum[jstr]) > 1.0e-30) {
					A[istartAnew2].aij = vector_sum[jstr];
					A[istartAnew2].i = istr;
					A[istartAnew2].j = jstr;
					istartAnew2++;
				}
			}
			index_size = 0;

		}

		delete[] vector_sum;
		delete[] index_visit;


		delete[] row_ind_AS;
		delete[] row_ind_AE;
		delete[] row_ind_PS;
		delete[] row_ind_PE;
	}

}

		if (bquicktypesort) {
			QuickSort_j(P, 1 + iaddR, iaddR + nnzR - 1);
		}
		else {
			HeapSort_j(P, /*n_a[ilevel - 1],*/ 1 + iaddR, iaddR + nnzR - 1);
		}

		// Копируем матрицу А следующего уровня влево вплотную к матице первоначального уровня.
		//integer icounter3 = 1;
		integer nsize = istartAnew2 - (istartAnew);
		for (integer i_1 = nnz_a[ilevel - 1] + 1+iadd, i_2=1; i_2 <= nsize; i_1++, i_2++) {
			integer i_right_position = istartAnew - 1 + i_2;
			A[i_1] = A[i_right_position];
		}


		printf("Prolongation is construct.\n");
		printf("Error interpolation is count %d\n", ipromah);
		printf("diagnostic ipromah_one=%d\n", ipromah_one);
		if (debug_reshime) system("pause");
	
		delete[] C_numerate;

		nnz_aRP[ilevel - 1] = nnzR - 1;
		iaddR += nnzR - 1;
		n_a[ilevel] = icounter-1;
		nnz_a[ilevel] = nsize;
		iadd += nnz_a[ilevel - 1];


		ilevel++;

		// сортировка А по i.
		//heapsort(A, key=i*n_a[ilevel - 1] + j, 1, nnz_a[ilevel - 1]);
		if (bquicktypesort) {
			QuickSort(A, /*n_a[ilevel - 1],*/ 1 + iadd, nnz_a[ilevel - 1] + iadd);
		}
		else {
			if (nnz_a[ilevel - 1] < heapsortsizelimit) {
				HeapSort(A, n_a[ilevel - 1], 1 + iadd, nnz_a[ilevel - 1] + iadd);
			}
			else {
				Ak* Aorig = &A[1 + iadd];
				MergeSort(Aorig, nnz_a[ilevel - 1]);
				Aorig = NULL;
			}
		}


		// debug
		//for (integer i_1 = 1 + iadd; i_1 <= iadd + nnz_a[ilevel - 1]; i_1++) {
			//if (A[i_1].i > n_a[ilevel - 1]) {
		//		printf("matrix incorrect i\n");
			//}
			//if (A[i_1].j > n_a[ilevel - 1]) {
			//	printf("%d ",A[i_1].j);
			//}
			//printf("A[%d].i=%d A[%d].j=%d A[%d].aij=%e\n",i_1,A[i_1].i,i_1,A[i_1].j,i_1,A[i_1].aij);
			//system("pause");
		//}
		

		//exit(1);
		printf("one level construct OK.\n");
		if (debug_reshime) system("pause");


		//printf("export b\n");
		//exporttecplot(b,n);

		//Real *test_coarse = new Real[n_a[1] + 1];

		// restriction
		//restriction(R, 1, nnz_aRP[0], flag, b, test_coarse, n_a[0], n_a[1]);
		//for (integer ii = 1; ii <= n_a[0]; ii++) {
			//b[ii] = 0.0;
		//}

	/*{
		Real *test_coarse1 = new Real[n_a[2] + 1];

		// restriction
		restriction(R, 1 + nnz_aRP[0], nnz_aRP[0] + nnz_aRP[1], flag, test_coarse, test_coarse1, n_a[1], n_a[2]);
		for (integer ii = 1; ii <= n_a[1]; ii++) {
			test_coarse[ii] = 0.0;
		}

	{
		Real *test_coarse2 = new Real[n_a[3] + 1];

		// restriction
		restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], flag, test_coarse1, test_coarse2, n_a[2], n_a[3]);
		for (integer ii = 1; ii <= n_a[2]; ii++) {
			test_coarse1[ii] = 0.0;
		}

		prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], flag, test_coarse1, test_coarse2, n_a[2], n_a[3]);
	}

	prolongation(P, 1 + nnz_aRP[0], nnz_aRP[0] + nnz_aRP[1], flag, test_coarse, test_coarse1, n_a[1], n_a[2]);
	}
	*/
	//prolongation(P, 1, nnz_aRP[0], flag, b, test_coarse, n_a[0], n_a[1]);

	//exporttecplot(b, n);
		// proverka start
		// на сетке 81х81 проверка успешно проходится.


		printf("proverka start\n");
		// первый уровень вложенности.
		if (ilevel > 1) {
			for (integer i = 1; i <= n; i++) {
				flag[i] = false;
			}

			for (integer ii77 = nnz_a[0] + 1; ii77 <= nnz_a[0] + nnz_a[1]; ii77++) {
				if (flag[A[ii77].i] == false) {
					integer istr77 = A[ii77].i;
					integer ic77 = ii77;
					integer icdiag = ii77;
					Real ap = 0.0;
					//x[istr] = b[istr];
					while ((ic77 <= nnz_a[0] + nnz_a[1]) && (A[ic77].i == istr77)) {
						if (A[ic77].j != istr77) {
							//x[istr] += -A[ic].aij*x[A[ic].j];
						}
						else {
							ap = A[ic77].aij;
							icdiag = ic77;
						}
						ic77++;
					}
					if (fabs(ap) < RealZERO) {
						printf("zero diagonal element %e in string %d in level 1 matrix", ap, istr77);
						system("PAUSE");
						//exit(1);
					}
					
					flag[A[ii77].i] = true;

				}
			}
		}
		

	//проверка конец

		delete[] count_sosed;
		delete[] row_startA;

	}// иерархия сеток построена.

	ilevel--;

	// 31.224s [50.986] 2D m=81 debug x64 acumulqtor
	// 13.792 [18.156] 2D m=81 realese x64 acumulqtor
	// 8.028s 2D m=81 debug x64 rozetka
	// 3.827 2D m=81 realese x64 rozetka

	printf("ilevel=%d\n", ilevel);
	for (integer i_1 = 0; i_1 <= ilevel; i_1++) {
		printf("n_a[%d]=%d nnz_a[%d]=%d\n",i_1, n_a[i_1], i_1, nnz_a[i_1]);
	}
	printf("Graph(Mesh) ierarhion is construct sucsseful...\n");
	if (debug_reshime) system("pause");
	//system("pause");
	//exit(1);

		/*
		// отладочная печать в рабочей версии требуется закоментировать.
		for (integer ii = 1; ii <= n; ii++) {
			flag_[ii] = false;
		}
		for (integer ii = 1 + iaddR; ii <= iaddR + nnzR - 1; ii++) {
			if (flag_[R[ii].i] == false) {
				integer istr = R[ii].i;
				integer ic7 = ii;
				while ((ic7 <= iaddR + nnzR - 1) && (R[ic7].i == istr)) {
					printf("%e ", R[ic7].aij);
					ic7++;
				}
				printf("\n");
				system("pause");
				flag_[R[ii].i] = true;
			}
		}
		*/

		

	/*
	//exporttecplot(b,n);

	Real *test_coarse = new Real[n_a[1] + 1];

	// restriction
	restriction(R, 1, nnz_aRP[0], flag, b, test_coarse, n_a[0], n_a[1]);
	for (integer ii = 1; ii <= n_a[0]; ii++) {
	b[ii] = 0.0;
	}

	{
	Real *test_coarse1 = new Real[n_a[2] + 1];

	// restriction
	restriction(R, 1+nnz_aRP[0], nnz_aRP[0]+nnz_aRP[1], flag, test_coarse, test_coarse1, n_a[1], n_a[2]);
	for (integer ii = 1; ii <= n_a[1]; ii++) {
	test_coarse[ii] = 0.0;
	}

	{
	Real *test_coarse2 = new Real[n_a[3] + 1];

	// restriction
	restriction(R, 1 + nnz_aRP[0]+nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1]+nnz_aRP[2], flag, test_coarse1, test_coarse2, n_a[2], n_a[3]);
	for (integer ii = 1; ii <= n_a[2]; ii++) {
	test_coarse1[ii] = 0.0;
	}

	prolongation(P, 1 + nnz_aRP[0]+ nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], flag, test_coarse1, test_coarse2, n_a[2], n_a[3]);
	}

	prolongation(P, 1 + nnz_aRP[0], nnz_aRP[0] + nnz_aRP[1], flag, test_coarse, test_coarse1, n_a[1], n_a[2]);
	}

	prolongation(P, 1, nnz_aRP[0], flag, b, test_coarse, n_a[0], n_a[1]);

	exporttecplot(b, n);
	*/

	// подготовка матрицы к cycling:
	/*
	// smoother.
	// 1 september 2015.
	void seidel(Ak* &A, integer istart, integer iend, Real* &x, Real* &b, bool* &flag, integer n)
	{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	for (integer i = 1; i <= n; i++) {
	flag[i] = false;
	}
	for (integer ii = istart; ii <= iend; ii++) {
	if (flag[A[ii].i] == false) {
	integer istr = A[ii].i;
	integer ic = ii;
	Real ap = 0.0;
	x[istr] = b[istr];
	while ((ic<=iend)&&(A[ic].i == istr)) {
	if (A[ic].j != istr) {
	x[istr] += -A[ic].aij*x[A[ic].j];
	}
	else ap = A[ic].aij;
	ic++;
	}
	if (fabs(ap) < RealZERO) {
	printf("zero diagonal elements in string %d",istr);
	getchar();
	exit(1);
	}
	else {
	x[istr] /= ap;
	}
	flag[A[ii].i] = true;
	}
	}


	} // seidel

	*/
	integer *row_ptr_start = new integer[4 * n_a[0] + 1];
	integer *row_ptr_end = new integer[4 * n_a[0] + 1];
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	for (integer ii = 1; ii <= nnz_a[0]; ii++) {
		if (flag[A[ii].i] == false) {
			integer istr = A[ii].i;
			integer ic = ii;
			integer icdiag = ii;
			row_ptr_start[istr] = ii;
			Real ap = 0.0;
			//x[istr] = b[istr];
			while ((ic <= nnz_a[0]) && (A[ic].i == istr)) {
				if (A[ic].j != istr) {
					//x[istr] += -A[ic].aij*x[A[ic].j];
				}
				else {
					ap = A[ic].aij;
					icdiag = ic;
				}
				ic++;
			}
			row_ptr_end[istr] = ic - 1;
			if (fabs(ap) < RealZERO) {
				printf("zero diagonal elements in string %d in basic matrix", istr);
				system("PAUSE");
				exit(1);
			}
			else {
				//x[istr] /= ap;
			}

			flag[A[ii].i] = true;
			Ak temp = A[ii];
			A[ii] = A[icdiag];
			A[icdiag] = temp;
			A[ii].aij = 1.0 / ap; // умножение быстрей деления.
		}
	}

	bool bstop = false;
	// первый уровень вложенности.
	if (ilevel > 1) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		
		for (integer ii =  nnz_a[0] + 1; ii <=  nnz_a[0] + nnz_a[1]; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <=  nnz_a[0] + nnz_a[1]) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0]] = ic - 1;
				if (fabs(ap) < RealZERO) {
					printf("zero diagonal element %e in string %d in level 1 matrix", ap, istr);
					system("PAUSE");
					//exit(1);
					bstop = true;
				}
				else {
					
					flag[A[ii].i] = true;
					Ak temp = A[ii];
					A[ii] = A[icdiag];
					A[icdiag] = temp;
					A[ii].aij = 1.0 / ap; // умножение быстрей деления.

				}

				
			}
		}
	}

	if (bstop) exit(1);

	// второй уровень вложенности.

	if (ilevel > 2) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		for (integer ii =  nnz_a[0] + nnz_a[1] + 1; ii <=  nnz_a[0] +  nnz_a[1] + nnz_a[2]; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0] + n_a[1]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <=  nnz_a[0] +  nnz_a[1] + nnz_a[2]) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0] + n_a[1]] = ic - 1;
				if (fabs(ap) < RealZERO) {
					printf("zero diagonal elements in string %d in level 2 matrix", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}


	// третий уровень вложенности.

	if (ilevel > 3) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		for (integer ii =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] + 1; ii <=  nnz_a[0] +  nnz_a[1] +  nnz_a[2] + nnz_a[3]; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0] + n_a[1] + n_a[2]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <=  nnz_a[0] +  nnz_a[1] +  nnz_a[2] + nnz_a[3]) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0] + n_a[1] + n_a[2]] = ic - 1;
				if (fabs(ap) < RealZERO) {
					printf("zero diagonal elements in string %d in level 3 matrix", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}


	// 14 сентября 2015 понедельник
	// четвёртый уровень вложенности.

	if (ilevel > 4) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		integer ist = nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] + 1;
		integer iend =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] + nnz_a[4];
		for (integer ii = ist; ii <= iend; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <= iend) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3]] = ic - 1;
				if (fabs(ap) < RealZERO) {
					printf("zero diagonal elements in string %d in level 4 matrix", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}


	// пятый уровень вложенности.

	if (ilevel > 5) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		integer ist = nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] + 1;
		integer iend =  nnz_a[0] +  nnz_a[1] + nnz_a[2] +  nnz_a[3] +  nnz_a[4] + nnz_a[5];
		for (integer ii = ist; ii <= iend; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <= iend) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4]] = ic - 1;
				if (fabs(ap) < RealZERO) {
					printf("zero diagonal elements in string %d in level 5 matrix", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}

	// шестой уровень вложенности.

	if (ilevel > 6) {
		for (integer i = 1; i <= n; i++) {
			flag[i] = false;
		}
		integer ist =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] + 1;
		integer iend =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] + nnz_a[6];
		for (integer ii = ist; ii <= iend; ii++) {
			if (flag[A[ii].i] == false) {
				integer istr = A[ii].i;
				integer ic = ii;
				integer icdiag = ii;
				row_ptr_start[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5]] = ii;
				Real ap = 0.0;
				//x[istr] = b[istr];
				while ((ic <= iend) && (A[ic].i == istr)) {
					if (A[ic].j != istr) {
						//x[istr] += -A[ic].aij*x[A[ic].j];
					}
					else {
						ap = A[ic].aij;
						icdiag = ic;
					}
					ic++;
				}
				row_ptr_end[istr + n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5]] = ic - 1;
				if (fabs(ap) < RealZERO) {
					printf("zero diagonal elements in string %d in level 6 matrix", istr);
					system("PAUSE");
					exit(1);
				}
				else {
					//x[istr] /= ap;
				}

				flag[A[ii].i] = true;
				Ak temp = A[ii];
				A[ii] = A[icdiag];
				A[icdiag] = temp;
				A[ii].aij = 1.0 / ap; // умножение быстрей деления.
			}
		}
	}


	

	// smoother.
	// 9 september 2015.
	// q - quick.
	// seidelq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);
	//void seidelq(Ak* &A, integer istartq, integer iendq, Real* &x, Real* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
	//{
	    // istart - начальная позиция ненулевых элементов в матрице А.
	    // iend - конечная позиция ненулевых элементов в матрице А.
	    //integer startpos = istartq + iadd;
	    //integer endpos = iendq+iadd;
	    //for (integer ii = startpos; ii <= endpos; ii++) {
	         //integer istr = ii - iadd;
	         //x[istr] = b[istr];
	         //for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++)
	         //{
	              //x[istr] += -A[ii1].aij*x[A[ii1].j];
	         //}
	         //x[istr] *= A[row_ptr_start[ii]].aij;
	    //}


	//} // seidelq
	

	delete[] this_is_C_node;
	delete[] this_is_F_node;


	printf("cycling: V cycle.\n");
	printf("level=%d\n", ilevel);
	printf("multigrid R.P.Fedorenko 1961.\n");
	printf("standart aglomerative algebraic multigrid method.\n");
	if (debug_reshime) system("pause");
	//exit(1);

	// 10 11 21 multigrid tutorial Вильм Бригг.

	integer nu1 = 4;
	integer nu2 = 3;

	//ilevel = 1; //debug
	Real rho = 1.0;
	Real dres = 1.0;
	int iiter = 1;
	const Real tolerance = 1.0e-12;


	Real *residual_fine = new Real[n_a[0] + 1];
	Real *residual_coarse = NULL;
	Real* error_approx_coarse = NULL;
	Real *residual_fine1 = NULL;
	Real *residual_coarse1 = NULL;
	Real* error_approx_coarse1 = NULL;
	Real *error_approx_fine1 = NULL;
	Real *residual_fine2 = NULL;
	Real *residual_coarse2 = NULL;
	Real* error_approx_coarse2 = NULL;
	Real *error_approx_fine2 = NULL;
	Real *residual_fine3 = NULL;
	Real *residual_coarse3 = NULL;
	Real* error_approx_coarse3 = NULL;
	Real *error_approx_fine3 = NULL;
	Real *residual_fine4 = NULL;
	Real *residual_coarse4 = NULL;
	Real *error_approx_coarse4 = NULL;
	Real *error_approx_fine4 = NULL;
	Real *residual_fine5 = NULL;
	Real *residual_coarse5 = NULL;
	Real* error_approx_coarse5 = NULL;
	Real *error_approx_fine5 = NULL;
	Real *residual_fine6 = NULL;
	Real *residual_coarse6 = NULL;
	Real* error_approx_coarse6 = NULL;
	Real *error_approx_fine6 = NULL;

	if (ilevel > 1) {
		residual_coarse = new Real[n_a[1] + 1];
		error_approx_coarse = new Real[n_a[1] + 1];
		if (ilevel > 2) {
			// residual
			residual_fine1 = new Real[n_a[1] + 1];
			residual_coarse1 = new Real[n_a[2] + 1];
			error_approx_coarse1 = new Real[n_a[2] + 1];
			error_approx_fine1 = new Real[n_a[1] + 1];
			if (ilevel > 3) {
				// residual
				residual_fine2 = new Real[n_a[2] + 1];
				residual_coarse2 = new Real[n_a[3] + 1];
				error_approx_coarse2 = new Real[n_a[3] + 1];
				error_approx_fine2 = new Real[n_a[2] + 1];
				if (ilevel > 4) {
					// residual
					residual_fine3 = new Real[n_a[3] + 1];
					residual_coarse3 = new Real[n_a[4] + 1];
					error_approx_coarse3 = new Real[n_a[4] + 1];
					error_approx_fine3 = new Real[n_a[3] + 1];
					if (ilevel > 5) {
						// residual
						residual_fine4 = new Real[n_a[4] + 1];
						residual_coarse4 = new Real[n_a[5] + 1];
						error_approx_coarse4 = new Real[n_a[5] + 1];
						error_approx_fine4 = new Real[n_a[4] + 1];
						if (ilevel > 6) {
							// residual
							residual_fine5 = new Real[n_a[5] + 1];
							residual_coarse5 = new Real[n_a[6] + 1];
							error_approx_coarse5 = new Real[n_a[6] + 1];
							error_approx_fine5 = new Real[n_a[5] + 1];
							if (ilevel > 7) {
								// residual
								residual_fine6 = new Real[n_a[6] + 1];
								residual_coarse6 = new Real[n_a[7] + 1];
								error_approx_coarse6 = new Real[n_a[7] + 1];
								error_approx_fine6 = new Real[n_a[6] + 1];
							}
						}
					}
				}
			}
		}
	}
	Real *error_approx_fine = new Real[n_a[0] + 1];


	//for (integer iprohod = 0; iprohod < 20; iprohod++) {
	while (dres>tolerance) {

		// smother
		for (integer iter = 0; iter < nu1; iter++) {
			//seidel(A, 1, nnz_a[0], x, b, flag, n_a[0]);
			//quick seidel
			seidelq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);

		}

		//exporttecplot(x, n);

		// residual_r
		//Real *residual_fine = new Real[n_a[0] + 1];
		//residual(A, 1, nnz_a[0], x, b, flag, n_a[0], residual_fine);
		residualq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine);
		dres = norma(residual_fine, n_a[0]);
		printf("%d %e rho=%e\n", iiter, dres, dres / rho);
		iiter++;
		//rho=norma(residual_fine, n_a[0]);
		rho = dres;
		//if (iprohod%5==0) getchar();
		if (ilevel > 1) {

			//Real *residual_coarse = new Real[n_a[1] + 1];

			// restriction
			restriction(R, 1, nnz_aRP[0], flag, residual_fine, residual_coarse, n_a[0], n_a[1]);

			// A*e=r;
			//Real* error_approx_coarse = new Real[n_a[1] + 1];
			for (integer ii = 1; ii <= n_a[1]; ii++) {
				error_approx_coarse[ii] = 0.0;
			}
			// pre smothing
			for (integer iter = 0; iter < nu1; iter++) {
				//seidel(A, 1 + 2 * nnz_a[0], 2 * nnz_a[0] + nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1]);
				seidelq(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0]);
			}

			if (ilevel > 2) {
				// residual
				//Real *residual_fine1 = new Real[n_a[1] + 1];
				//residual(A, 1+2*nnz_a[0], 2*nnz_a[0]+nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1], residual_fine1);
				residualq(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0], residual_fine1);
				//residualqspeshial(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0], residual_fine1);


				//Real *residual_coarse1 = new Real[n_a[2] + 1];

				// restriction
				restriction(R, 1 + nnz_aRP[0], nnz_aRP[0] + nnz_aRP[1], flag, residual_fine1, residual_coarse1, n_a[1], n_a[2]);

				// A*e=r;
				//Real* error_approx_coarse1 = new Real[n_a[2] + 1];
				for (integer ii = 1; ii <= n_a[2]; ii++) {
					error_approx_coarse1[ii] = 0.0;
				}
				// pre smothing
				for (integer iter = 0; iter < nu1; iter++) {
					//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1], 2 * nnz_a[0] + 2 * nnz_a[1] + nnz_a[2], error_approx_coarse1, residual_coarse1, flag, n_a[2]);
					seidelq(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0] + n_a[1]);
				}
				if (ilevel > 3) {
					// residual
					//Real *residual_fine2 = new Real[n_a[2] + 1];
					//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1], 2 * nnz_a[0] + 2 * nnz_a[1] + nnz_a[2], error_approx_coarse1, residual_coarse1, flag, n_a[2], residual_fine2);
					residualq(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0]+n_a[1], residual_fine2);
					//residualqspeshial(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0] + n_a[1], residual_fine2);

					//Real *residual_coarse2 = new Real[n_a[3] + 1];

					// restriction
					restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], flag, residual_fine2, residual_coarse2, n_a[2], n_a[3]);

					// A*e=r;
					//Real* error_approx_coarse2 = new Real[n_a[3] + 1];
					for (integer ii = 1; ii <= n_a[3]; ii++) {
						error_approx_coarse2[ii] = 0.0;
					}
					// pre smothing
					for (integer iter = 0; iter < nu1; iter++) {
						//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + nnz_a[3], error_approx_coarse2, residual_coarse2, flag, n_a[3]);
						seidelq(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2]);

					}
					if (ilevel > 4) {
						// residual
						//Real *residual_fine3 = new Real[n_a[3] + 1];
						//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + nnz_a[3], error_approx_coarse2, residual_coarse2, flag, n_a[3], residual_fine3);
						residualq(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2], residual_fine3);
						//speshial
						//residualqspeshial(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2], residual_fine3);



						//Real *residual_coarse3 = new Real[n_a[4] + 1];

						// restriction
						restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], flag, residual_fine3, residual_coarse3, n_a[3], n_a[4]);

						// A*e=r;
						//Real* error_approx_coarse3 = new Real[n_a[4] + 1];
						for (integer ii = 1; ii <= n_a[4]; ii++) {
							error_approx_coarse3[ii] = 0.0;
						}
						// pre smothing
						for (integer iter = 0; iter < nu1; iter++) {
							//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + nnz_a[4], error_approx_coarse3, residual_coarse3, flag, n_a[4]);
							seidelq(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3]);
						}
						if (ilevel > 5) {
							// residual
							//Real *residual_fine4 = new Real[n_a[4] + 1];
							//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + nnz_a[4], error_approx_coarse3, residual_coarse3, flag, n_a[4], residual_fine4);
							residualq(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3], residual_fine4);
							//speshial 14 september 2015.
							//residualqspeshial(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3], residual_fine4);


							//Real *residual_coarse4 = new Real[n_a[5] + 1];

							// restriction
							restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4], flag, residual_fine4, residual_coarse4, n_a[4], n_a[5]);

							// A*e=r;
							//Real* error_approx_coarse4 = new Real[n_a[5] + 1];
							for (integer ii = 1; ii <= n_a[5]; ii++) {
								error_approx_coarse4[ii] = 0.0;
							}
							// pre smothing
							for (integer iter = 0; iter < nu1; iter++) {
								//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + nnz_a[5], error_approx_coarse4, residual_coarse4, flag, n_a[5]);
								seidelq(A, 1, n_a[5], error_approx_coarse4, residual_coarse4, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4]);
							}
							if (ilevel > 6) {
								// residual
								//Real *residual_fine5 = new Real[n_a[5] + 1];
								//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + nnz_a[5], error_approx_coarse4, residual_coarse4, flag, n_a[5], residual_fine5);
								//if (ilevel <= 15) {
								residualq(A, 1, n_a[5], error_approx_coarse4, residual_coarse4, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4], residual_fine5);
								//}
								//else {
								// приводит к расходимости.
								//speshial 14 september 2015.
								// это уже приводит к увеличению числа итераций на примере сетки в 1млн узлов. остановимся.
								//residualqspeshial(A, 1, n_a[5], error_approx_coarse4, residual_coarse4, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4], residual_fine5);
								//}

								//Real *residual_coarse5 = new Real[n_a[6] + 1];

								// restriction
								restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], flag, residual_fine5, residual_coarse5, n_a[5], n_a[6]);

								// A*e=r;
								//Real* error_approx_coarse5 = new Real[n_a[6] + 1];
								for (integer ii = 1; ii <= n_a[6]; ii++) {
									error_approx_coarse5[ii] = 0.0;
								}
								// pre smothing
								for (integer iter = 0; iter < nu1; iter++) {
									//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + nnz_a[6], error_approx_coarse5, residual_coarse5, flag, n_a[6]);
									seidelq(A, 1, n_a[6], error_approx_coarse5, residual_coarse5, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5]);
								}

								if (ilevel > 7) {
									// residual
									//Real *residual_fine6 = new Real[n_a[6] + 1];
									//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] +2*nnz_a[5], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + nnz_a[6], error_approx_coarse5, residual_coarse5, flag, n_a[6], residual_fine6);
									residualq(A, 1, n_a[6], error_approx_coarse5, residual_coarse5, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5], residual_fine6);

									//Real *residual_coarse6 = new Real[n_a[7] + 1];

									// restriction
									restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6], flag, residual_fine6, residual_coarse6, n_a[6], n_a[7]);

									// A*e=r;
									//Real* error_approx_coarse6 = new Real[n_a[7] + 1];
									for (integer ii = 1; ii <= n_a[7]; ii++) {
										error_approx_coarse6[ii] = 0.0;
									}
									// pre smothing
									for (integer iter = 0; iter < nu1; iter++) {
										seidel(A, 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] + nnz_a[4] +  nnz_a[5] +  nnz_a[6],  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] + nnz_a[4] +  nnz_a[5] +  nnz_a[6] + nnz_a[7], error_approx_coarse6, residual_coarse6, flag, n_a[7]);
									}

									if (ilevel > 8) {
										// residual
										Real *residual_fine7 = new Real[n_a[7] + 1];
										residual(A, 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6],  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] + nnz_a[7], error_approx_coarse6, residual_coarse6, flag, n_a[7], residual_fine7);


										Real *residual_coarse7 = new Real[n_a[8] + 1];

										// restriction
										restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7], flag, residual_fine7, residual_coarse7, n_a[7], n_a[8]);

										// A*e=r;
										Real* error_approx_coarse7 = new Real[n_a[8] + 1];
										for (integer ii = 1; ii <= n_a[8]; ii++) {
											error_approx_coarse7[ii] = 0.0;
										}
										// pre smothing
										for (integer iter = 0; iter < nu1; iter++) {
											seidel(A, 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7],  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] + nnz_a[8], error_approx_coarse7, residual_coarse7, flag, n_a[8]);
										}


										if (ilevel > 9) {
											// residual
											Real *residual_fine8 = new Real[n_a[8] + 1];
											integer n1 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7];
											integer n2 =  nnz_a[0] +  nnz_a[1] + nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] + nnz_a[8];
											residual(A, n1, n2, error_approx_coarse7, residual_coarse7, flag, n_a[8], residual_fine8);


											Real *residual_coarse8 = new Real[n_a[9] + 1];

											// restriction
											integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7];
											integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8];
											restriction(R, n3, n4, flag, residual_fine8, residual_coarse8, n_a[8], n_a[9]);

											// A*e=r;
											Real* error_approx_coarse8 = new Real[n_a[9] + 1];
											for (integer ii = 1; ii <= n_a[9]; ii++) {
												error_approx_coarse8[ii] = 0.0;
											}
											// pre smothing
											for (integer iter = 0; iter < nu1; iter++) {
												integer n5 = 1 + nnz_a[0] + nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8];
												integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] + nnz_a[9];
												seidel(A, n5, n6, error_approx_coarse8, residual_coarse8, flag, n_a[9]);
											}

											if (ilevel > 10) {
												// 8 сентября 2015 РИМИНИ пляж 

												// residual
												Real *residual_fine9 = new Real[n_a[9] + 1];
												integer n1 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8];
												integer n2 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] + nnz_a[9];
												residual(A, n1, n2, error_approx_coarse8, residual_coarse8, flag, n_a[9], residual_fine9);


												Real *residual_coarse9 = new Real[n_a[10] + 1];

												// restriction
												integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8];
												integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9];
												restriction(R, n3, n4, flag, residual_fine9, residual_coarse9, n_a[9], n_a[10]);

												// A*e=r;
												Real* error_approx_coarse9 = new Real[n_a[10] + 1];
												for (integer ii = 1; ii <= n_a[10]; ii++) {
													error_approx_coarse9[ii] = 0.0;
												}
												// pre smothing
												for (integer iter = 0; iter < nu1; iter++) {
													integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9];
													integer n6 = nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] + nnz_a[10];
													seidel(A, n5, n6, error_approx_coarse9, residual_coarse9, flag, n_a[10]);
												}

												if (ilevel > 11) {
													// 8 сентября 2015 РИМИНИ пляж 

													// residual
													Real *residual_fine10 = new Real[n_a[10] + 1];
													integer n1 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9];
													integer n2 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] + nnz_a[10];
													residual(A, n1, n2, error_approx_coarse9, residual_coarse9, flag, n_a[10], residual_fine10);


													Real *residual_coarse10 = new Real[n_a[11] + 1];

													// restriction
													integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9];
													integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10];
													restriction(R, n3, n4, flag, residual_fine10, residual_coarse10, n_a[10], n_a[11]);

													// A*e=r;
													Real* error_approx_coarse10 = new Real[n_a[11] + 1];
													for (integer ii = 1; ii <= n_a[11]; ii++) {
														error_approx_coarse10[ii] = 0.0;
													}
													// pre smothing
													for (integer iter = 0; iter < nu1; iter++) {
														integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10];
														integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] + nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] + nnz_a[11];
														seidel(A, n5, n6, error_approx_coarse10, residual_coarse10, flag, n_a[11]);
													}

													if (ilevel > 12) {
														// 11 сентября 2015 РИМИНИ пляж 

														// residual
														Real *residual_fine11 = new Real[n_a[11] + 1];
														integer n1 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] + nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10];
														integer n2 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] + nnz_a[11];
														residual(A, n1, n2, error_approx_coarse10, residual_coarse10, flag, n_a[11], residual_fine11);


														Real *residual_coarse11 = new Real[n_a[12] + 1];

														// restriction
														integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10];
														integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11];
														restriction(R, n3, n4, flag, residual_fine11, residual_coarse11, n_a[11], n_a[12]);

														// A*e=r;
														Real* error_approx_coarse11 = new Real[n_a[12] + 1];
														for (integer ii = 1; ii <= n_a[12]; ii++) {
															error_approx_coarse11[ii] = 0.0;
														}
														// pre smothing
														for (integer iter = 0; iter < nu1; iter++) {
															integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11];
															integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] + nnz_a[12];
															seidel(A, n5, n6, error_approx_coarse11, residual_coarse11, flag, n_a[12]);
														}

														if (ilevel > 13) {
															// 11 сентября 2015 РИМИНИ пляж 

															// residual
															Real *residual_fine12 = new Real[n_a[12] + 1];
															integer n1 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11];
															integer n2 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] + nnz_a[12];
															residual(A, n1, n2, error_approx_coarse11, residual_coarse11, flag, n_a[12], residual_fine12);


															Real *residual_coarse12 = new Real[n_a[13] + 1];

															// restriction
															integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11];
															integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12];
															restriction(R, n3, n4, flag, residual_fine12, residual_coarse12, n_a[12], n_a[13]);

															// A*e=r;
															Real* error_approx_coarse12 = new Real[n_a[13] + 1];
															for (integer ii = 1; ii <= n_a[13]; ii++) {
																error_approx_coarse12[ii] = 0.0;
															}
															// pre smothing
															for (integer iter = 0; iter < nu1; iter++) {
																integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12];
																integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] + nnz_a[13];
																seidel(A, n5, n6, error_approx_coarse12, residual_coarse12, flag, n_a[13]);
															}


															if (ilevel > 14) {
																// 11 сентября 2015 РИМИНИ пляж 

																// residual
																Real *residual_fine13 = new Real[n_a[13] + 1];
																integer n1 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] + nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12];
																integer n2 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] + nnz_a[13];
																residual(A, n1, n2, error_approx_coarse12, residual_coarse12, flag, n_a[13], residual_fine13);


																Real *residual_coarse13 = new Real[n_a[14] + 1];

																// restriction
																integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12];
																integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13];
																restriction(R, n3, n4, flag, residual_fine13, residual_coarse13, n_a[13], n_a[14]);

																// A*e=r;
																Real* error_approx_coarse13 = new Real[n_a[14] + 1];
																for (integer ii = 1; ii <= n_a[14]; ii++) {
																	error_approx_coarse13[ii] = 0.0;
																}
																// pre smothing
																for (integer iter = 0; iter < nu1; iter++) {
																	integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] + nnz_a[9] +  nnz_a[10] + nnz_a[11] +  nnz_a[12] +  nnz_a[13];
																	integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] +  nnz_a[13] + nnz_a[14];
																	seidel(A, n5, n6, error_approx_coarse13, residual_coarse13, flag, n_a[14]);
																}

																if (ilevel > 15) {
																	// 14 сентября 2015 Москва на работе в пн. 

																	// residual
																	Real *residual_fine14 = new Real[n_a[14] + 1];
																	integer n1 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] +  nnz_a[13];
																	integer n2 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] + nnz_a[10] +  nnz_a[11] +  nnz_a[12] +  nnz_a[13] + nnz_a[14];
																	residual(A, n1, n2, error_approx_coarse13, residual_coarse13, flag, n_a[14], residual_fine14);


																	Real *residual_coarse14 = new Real[n_a[15] + 1];

																	// restriction
																	integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13];
																	integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14];
																	restriction(R, n3, n4, flag, residual_fine14, residual_coarse14, n_a[14], n_a[15]);

																	// A*e=r;
																	Real* error_approx_coarse14 = new Real[n_a[15] + 1];
																	for (integer ii = 1; ii <= n_a[15]; ii++) {
																		error_approx_coarse14[ii] = 0.0;
																	}
																	// pre smothing
																	for (integer iter = 0; iter < nu1; iter++) {
																		integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] +  nnz_a[13] +  nnz_a[14];
																		integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] +  nnz_a[13] +  nnz_a[14] + nnz_a[15];
																		seidel(A, n5, n6, error_approx_coarse14, residual_coarse14, flag, n_a[15]);
																	}

																	if (ilevel > 16) {
																		// 10 октября 2015. 

																		// residual
																		Real *residual_fine15 = new Real[n_a[15] + 1];
																		integer n1 = 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13]+nnz_a[14];
																		integer n2 = nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14]+nnz_a[15];
																		residual(A, n1, n2, error_approx_coarse14, residual_coarse14, flag, n_a[15], residual_fine15);


																		Real *residual_coarse15 = new Real[n_a[16] + 1];

																		// restriction
																		integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13]+nnz_aRP[14];
																		integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14]+nnz_aRP[15];
																		restriction(R, n3, n4, flag, residual_fine15, residual_coarse15, n_a[15], n_a[16]);

																		// A*e=r;
																		Real* error_approx_coarse15 = new Real[n_a[16] + 1];
																		for (integer ii = 1; ii <= n_a[16]; ii++) {
																			error_approx_coarse15[ii] = 0.0;
																		}
																		// pre smothing
																		for (integer iter = 0; iter < nu1; iter++) {
																			integer n5 = 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15];
																			integer n6 = nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15] + nnz_a[16];
																			seidel(A, n5, n6, error_approx_coarse15, residual_coarse15, flag, n_a[16]);
																		}

																		if (ilevel > 17) {
																			// 10 октября 2015. 

																			// residual
																			Real *residual_fine16 = new Real[n_a[16] + 1];
																			integer n1 = 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15];
																			integer n2 = nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15] + nnz_a[16];
																			residual(A, n1, n2, error_approx_coarse15, residual_coarse15, flag, n_a[16], residual_fine16);


																			Real *residual_coarse16 = new Real[n_a[17] + 1];

																			// restriction
																			integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15];
																			integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16];
																			restriction(R, n3, n4, flag, residual_fine16, residual_coarse16, n_a[16], n_a[17]);

																			// A*e=r;
																			Real* error_approx_coarse16 = new Real[n_a[17] + 1];
																			for (integer ii = 1; ii <= n_a[17]; ii++) {
																				error_approx_coarse16[ii] = 0.0;
																			}
																			// pre smothing
																			for (integer iter = 0; iter < nu1; iter++) {
																				integer n5 = 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15] + nnz_a[16];
																				integer n6 = nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15] + nnz_a[16] + nnz_a[17];
																				seidel(A, n5, n6, error_approx_coarse16, residual_coarse16, flag, n_a[17]);
																			}

																			// post smothing
																			for (integer iter = 0; iter < nu2; iter++) {
																				integer n5 = 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15] + nnz_a[16];
																				integer n6 = nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15] + nnz_a[16] + nnz_a[17];
																				seidel(A, n5, n6, error_approx_coarse16, residual_coarse16, flag, n_a[17]);
																			}


																			// prolongation
																			// residual_r
																			Real *error_approx_fine16 = new Real[n_a[16] + 1];
																			for (integer ii = 1; ii <= n_a[16]; ii++) {
																				error_approx_fine16[ii] = 0.0;
																			}



																			integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15];
																			integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16];
																			prolongation(P, n7, n8, flag, error_approx_fine16, error_approx_coarse16, n_a[16], n_a[17]);

																			// correction
																			for (integer ii = 1; ii <= n_a[16]; ii++) {
																				error_approx_coarse15[ii] += error_approx_fine16[ii];
																			}

																			// free
																			delete[] error_approx_fine16;
																			delete[] error_approx_coarse16;
																			delete[] residual_coarse16;
																			delete[] residual_fine16;

																		}


																		// post smothing
																		for (integer iter = 0; iter < nu2; iter++) {
																			integer n5 = 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15];
																			integer n6 = nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7] + nnz_a[8] + nnz_a[9] + nnz_a[10] + nnz_a[11] + nnz_a[12] + nnz_a[13] + nnz_a[14] + nnz_a[15] + nnz_a[16];
																			seidel(A, n5, n6, error_approx_coarse15, residual_coarse15, flag, n_a[16]);
																		}


																		// prolongation
																		// residual_r
																		Real *error_approx_fine15 = new Real[n_a[15] + 1];
																		for (integer ii = 1; ii <= n_a[15]; ii++) {
																			error_approx_fine15[ii] = 0.0;
																		}



																		integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13]+nnz_aRP[14];
																		integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15];
																		prolongation(P, n7, n8, flag, error_approx_fine15, error_approx_coarse15, n_a[15], n_a[16]);

																		// correction
																		for (integer ii = 1; ii <= n_a[15]; ii++) {
																			error_approx_coarse14[ii] += error_approx_fine15[ii];
																		}

																		// free
																		delete[] error_approx_fine15;
																		delete[] error_approx_coarse15;
																		delete[] residual_coarse15;
																		delete[] residual_fine15;

																	}

																	// post smothing
																	for (integer iter = 0; iter < nu2; iter++) {
																		integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] + nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] +  nnz_a[13] +  nnz_a[14];
																		integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] +  nnz_a[13] +  nnz_a[14] + nnz_a[15];
																		seidel(A, n5, n6, error_approx_coarse14, residual_coarse14, flag, n_a[15]);
																	}


																	// prolongation
																	// residual_r
																	Real *error_approx_fine14 = new Real[n_a[14] + 1];
																	for (integer ii = 1; ii <= n_a[14]; ii++) {
																		error_approx_fine14[ii] = 0.0;
																	}

																	

																	integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13];
																	integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14];
																	prolongation(P, n7, n8, flag, error_approx_fine14, error_approx_coarse14, n_a[14], n_a[15]);

																	// correction
																	for (integer ii = 1; ii <= n_a[14]; ii++) {
																		error_approx_coarse13[ii] += error_approx_fine14[ii];
																	}

																	// free
																	delete[] error_approx_fine14;
																	delete[] error_approx_coarse14;
																	delete[] residual_coarse14;
																	delete[] residual_fine14;

																}


																// post smothing
																for (integer iter = 0; iter < nu2; iter++) {
																	integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] +  nnz_a[13];
																	integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] +  nnz_a[13] + nnz_a[14];
																	seidel(A, n5, n6, error_approx_coarse13, residual_coarse13, flag, n_a[14]);
																}


																// prolongation
																// residual_r
																Real *error_approx_fine13 = new Real[n_a[13] + 1];
																for (integer ii = 1; ii <= n_a[13]; ii++) {
																	error_approx_fine13[ii] = 0.0;
																}

																//for (integer ii = 1; ii <= n_a[14]; ii++) {// debug
																//printf("error_approx_coarse13[%d]=%e\n",ii, error_approx_coarse13[ii]);

																//printf("residual_coarse13[%d]=%e\n", ii, residual_coarse13[ii]);
																//getchar();
																//}
																//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12]+2*nnz_a[13]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12] +2*nnz_a[13]+ nnz_a[14]; ii++) {// debug
																//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
																//if (ii % 20 == 0) getchar();
																//}

																integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12];
																integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13];
																prolongation(P, n7, n8, flag, error_approx_fine13, error_approx_coarse13, n_a[13], n_a[14]);

																// correction
																for (integer ii = 1; ii <= n_a[13]; ii++) {
																	error_approx_coarse12[ii] += error_approx_fine13[ii];
																}

																// free
																delete[] error_approx_fine13;
																delete[] error_approx_coarse13;
																delete[] residual_coarse13;
																delete[] residual_fine13;

															}


															// post smothing
															for (integer iter = 0; iter < nu2; iter++) {
																integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] + nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12];
																integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] +  nnz_a[12] + nnz_a[13];
																seidel(A, n5, n6, error_approx_coarse12, residual_coarse12, flag, n_a[13]);
															}


															// prolongation
															// residual_r
															Real *error_approx_fine12 = new Real[n_a[12] + 1];
															for (integer ii = 1; ii <= n_a[12]; ii++) {
																error_approx_fine12[ii] = 0.0;
															}

															//for (integer ii = 1; ii <= n_a[13]; ii++) {// debug
															//printf("error_approx_coarse12[%d]=%e\n",ii, error_approx_coarse12[ii]);

															//printf("residual_coarse12[%d]=%e\n", ii, residual_coarse12[ii]);
															//getchar();
															//}
															//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+2*nnz_a[12]+ nnz_a[13]; ii++) {// debug
															//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
															//if (ii % 20 == 0) getchar();
															//}

															integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11];
															integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12];
															prolongation(P, n7, n8, flag, error_approx_fine12, error_approx_coarse12, n_a[12], n_a[13]);

															// correction
															for (integer ii = 1; ii <= n_a[12]; ii++) {
																error_approx_coarse11[ii] += error_approx_fine12[ii];
															}

															// free
															delete[] error_approx_fine12;
															delete[] error_approx_coarse12;
															delete[] residual_coarse12;
															delete[] residual_fine12;

														}



														// post smothing
														for (integer iter = 0; iter < nu2; iter++) {
															integer n5 = 1 +  nnz_a[0] + nnz_a[1] +  nnz_a[2] +  nnz_a[3] + nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11];
															integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] + nnz_a[4] +  nnz_a[5] + nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] +  nnz_a[11] + nnz_a[12];
															seidel(A, n5, n6, error_approx_coarse11, residual_coarse11, flag, n_a[12]);
														}


														// prolongation
														// residual_r
														Real *error_approx_fine11 = new Real[n_a[11] + 1];
														for (integer ii = 1; ii <= n_a[11]; ii++) {
															error_approx_fine11[ii] = 0.0;
														}

														//for (integer ii = 1; ii <= n_a[12]; ii++) {// debug
														//printf("error_approx_coarse11[%d]=%e\n",ii, error_approx_coarse11[ii]);

														//printf("residual_coarse11[%d]=%e\n", ii, residual_coarse11[ii]);
														//getchar();
														//}
														//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+2*nnz_a[11]+ nnz_a[12]; ii++) {// debug
														//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
														//if (ii % 20 == 0) getchar();
														//}

														integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10];
														integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11];
														prolongation(P, n7, n8, flag, error_approx_fine11, error_approx_coarse11, n_a[11], n_a[12]);

														// correction
														for (integer ii = 1; ii <= n_a[11]; ii++) {
															error_approx_coarse10[ii] += error_approx_fine11[ii];
														}

														// free
														delete[] error_approx_fine11;
														delete[] error_approx_coarse11;
														delete[] residual_coarse11;
														delete[] residual_fine11;

													}


													// post smothing
													for (integer iter = 0; iter < nu2; iter++) {
														integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10];
														integer n6 = nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] + nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] +  nnz_a[10] + nnz_a[11];
														seidel(A, n5, n6, error_approx_coarse10, residual_coarse10, flag, n_a[11]);
													}


													// prolongation
													// residual_r
													Real *error_approx_fine10 = new Real[n_a[10] + 1];
													for (integer ii = 1; ii <= n_a[10]; ii++) {
														error_approx_fine10[ii] = 0.0;
													}

													//for (integer ii = 1; ii <= n_a[11]; ii++) {// debug
													//printf("error_approx_coarse10[%d]=%e\n",ii, error_approx_coarse10[ii]);

													//printf("residual_coarse10[%d]=%e\n", ii, residual_coarse10[ii]);
													//getchar();
													//}
													//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+2*nnz_a[10]+ nnz_a[11]; ii++) {// debug
													//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
													//if (ii % 20 == 0) getchar();
													//}

													integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9];
													integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10];
													prolongation(P, n7, n8, flag, error_approx_fine10, error_approx_coarse10, n_a[10], n_a[11]);

													// correction
													for (integer ii = 1; ii <= n_a[10]; ii++) {
														error_approx_coarse9[ii] += error_approx_fine10[ii];
													}

													// free
													delete[] error_approx_fine10;
													delete[] error_approx_coarse10;
													delete[] residual_coarse10;
													delete[] residual_fine10;

												}



												// post smothing
												for (integer iter = 0; iter < nu2; iter++) {
													integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9];
													integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] +  nnz_a[9] + nnz_a[10];
													seidel(A, n5, n6, error_approx_coarse9, residual_coarse9, flag, n_a[10]);
												}


												// prolongation
												// residual_r
												Real *error_approx_fine9 = new Real[n_a[9] + 1];
												for (integer ii = 1; ii <= n_a[9]; ii++) {
													error_approx_fine9[ii] = 0.0;
												}

												//for (integer ii = 1; ii <= n_a[10]; ii++) {// debug
												//printf("error_approx_coarse9[%d]=%e\n",ii, error_approx_coarse9[ii]);

												//printf("residual_coarse9[%d]=%e\n", ii, residual_coarse9[ii]);
												//getchar();
												//}
												//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+2*nnz_a[9]+ nnz_a[10]; ii++) {// debug
												//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
												//if (ii % 20 == 0) getchar();
												//}

												integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8];
												integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9];
												prolongation(P, n7, n8, flag, error_approx_fine9, error_approx_coarse9, n_a[9], n_a[10]);

												// correction
												for (integer ii = 1; ii <= n_a[9]; ii++) {
													error_approx_coarse8[ii] += error_approx_fine9[ii];
												}

												// free
												delete[] error_approx_fine9;
												delete[] error_approx_coarse9;
												delete[] residual_coarse9;
												delete[] residual_fine9;

											}

											// post smothing
											for (integer iter = 0; iter < nu2; iter++) {
												integer n5 = 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8];
												integer n6 =  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] +  nnz_a[8] + nnz_a[9];
												seidel(A, n5, n6, error_approx_coarse8, residual_coarse8, flag, n_a[9]);
											}


											// prolongation
											// residual_r
											Real *error_approx_fine8 = new Real[n_a[8] + 1];
											for (integer ii = 1; ii <= n_a[8]; ii++) {
												error_approx_fine8[ii] = 0.0;
											}

											//for (integer ii = 1; ii <= n_a[9]; ii++) {// debug
											//printf("error_approx_coarse8[%d]=%e\n",ii, error_approx_coarse8[ii]);

											//printf("residual_coarse8[%d]=%e\n", ii, residual_coarse8[ii]);
											//getchar();
											//}
											//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+2*nnz_a[8]+ nnz_a[9]; ii++) {// debug
											//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
											//if (ii % 20 == 0) getchar();
											//}

											integer n7 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7];
											integer n8 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8];
											prolongation(P, n7, n8, flag, error_approx_fine8, error_approx_coarse8, n_a[8], n_a[9]);

											// correction
											for (integer ii = 1; ii <= n_a[8]; ii++) {
												error_approx_coarse7[ii] += error_approx_fine8[ii];
											}

											// free
											delete[] error_approx_fine8;
											delete[] error_approx_coarse8;
											delete[] residual_coarse8;
											delete[] residual_fine8;

										}

										// post smothing
										for (integer iter = 0; iter < nu2; iter++) {
											seidel(A, 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] + nnz_a[6] +  nnz_a[7],  nnz_a[0] + nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] +  nnz_a[7] + nnz_a[8], error_approx_coarse7, residual_coarse7, flag, n_a[8]);
										}


										// prolongation
										// residual_r
										Real *error_approx_fine7 = new Real[n_a[7] + 1];
										for (integer ii = 1; ii <= n_a[7]; ii++) {
											error_approx_fine7[ii] = 0.0;
										}

										//for (integer ii = 1; ii <= n_a[8]; ii++) {// debug
										//printf("error_approx_coarse7[%d]=%e\n",ii, error_approx_coarse7[ii]);

										//printf("residual_coarse7[%d]=%e\n", ii, residual_coarse7[ii]);
										//getchar();
										//}
										//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+2*nnz_a[7]+ nnz_a[8]; ii++) {// debug
										//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
										//if (ii % 20 == 0) getchar();
										//}

										prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7], flag, error_approx_fine7, error_approx_coarse7, n_a[7], n_a[8]);

										// correction
										for (integer ii = 1; ii <= n_a[7]; ii++) {
											error_approx_coarse6[ii] += error_approx_fine7[ii];
										}

										// free
										delete[] error_approx_fine7;
										delete[] error_approx_coarse7;
										delete[] residual_coarse7;
										delete[] residual_fine7;

									}


									// post smothing
									for (integer iter = 0; iter < nu2; iter++) {
										seidel(A, 1 +  nnz_a[0] +  nnz_a[1] +  nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6],  nnz_a[0] +  nnz_a[1] + nnz_a[2] +  nnz_a[3] +  nnz_a[4] +  nnz_a[5] +  nnz_a[6] + nnz_a[7], error_approx_coarse6, residual_coarse6, flag, n_a[7]);
									}


									// prolongation
									// residual_r
									//Real *error_approx_fine6 = new Real[n_a[6] + 1];
									for (integer ii = 1; ii <= n_a[6]; ii++) {
										error_approx_fine6[ii] = 0.0;
									}

									//for (integer ii = 1; ii <= n_a[7]; ii++) {// debug
									//printf("error_approx_coarse6[%d]=%e\n",ii, error_approx_coarse6[ii]);

									//printf("residual_coarse6[%d]=%e\n", ii, residual_coarse6[ii]);
									//getchar();
									//}
									//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+2*nnz_a[6]+ nnz_a[7]; ii++) {// debug
									//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
									//if (ii % 20 == 0) getchar();
									//}

									prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6], flag, error_approx_fine6, error_approx_coarse6, n_a[6], n_a[7]);

									// correction
									for (integer ii = 1; ii <= n_a[6]; ii++) {
										error_approx_coarse5[ii] += error_approx_fine6[ii];
									}

									// free
									//delete[] error_approx_fine6;
									//delete[] error_approx_coarse6;
									//delete[] residual_coarse6;
									//delete[] residual_fine6;

								}

								// post smothing
								for (integer iter = 0; iter < nu2; iter++) {
									//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + nnz_a[6], error_approx_coarse5, residual_coarse5, flag, n_a[6]);
									seidelq(A, 1, n_a[6], error_approx_coarse5, residual_coarse5, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5]);
								}

								// prolongation
								// residual_r
								//Real *error_approx_fine5 = new Real[n_a[5] + 1];
								for (integer ii = 1; ii <= n_a[5]; ii++) {
									error_approx_fine5[ii] = 0.0;
								}

								//for (integer ii = 1; ii <= n_a[6]; ii++) {// debug
								//printf("error_approx_coarse5[%d]=%e\n",ii, error_approx_coarse5[ii]);

								//printf("residual_coarse5[%d]=%e\n", ii, residual_coarse5[ii]);
								//getchar();
								//}
								//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ 2 * nnz_a[5]+ nnz_a[6]; ii++) {// debug
								//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
								//if (ii % 20 == 0) getchar();
								//}

								prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], flag, error_approx_fine5, error_approx_coarse5, n_a[5], n_a[6]);

								// correction
								for (integer ii = 1; ii <= n_a[5]; ii++) {
									error_approx_coarse4[ii] += error_approx_fine5[ii];
								}

								// free
								//delete[] error_approx_fine5;
								//delete[] error_approx_coarse5;
								//delete[] residual_coarse5;
								//delete[] residual_fine5;

							}
							// post smothing
							for (integer iter = 0; iter < nu2; iter++) {
								//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + nnz_a[5], error_approx_coarse4, residual_coarse4, flag, n_a[5]);
								seidelq(A, 1, n_a[5], error_approx_coarse4, residual_coarse4, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4]);
							}

							// prolongation
							// residual_r
							//Real *error_approx_fine4 = new Real[n_a[4] + 1];
							for (integer ii = 1; ii <= n_a[4]; ii++) {
								error_approx_fine4[ii] = 0.0;
							}

							//for (integer ii = 1; ii <= n_a[5]; ii++) {// debug
							//printf("error_approx_coarse4[%d]=%e\n",ii, error_approx_coarse4[ii]);

							//printf("residual_coarse4[%d]=%e\n", ii, residual_coarse4[ii]);
							//getchar();
							//}
							//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ 2 * nnz_a[4]+ nnz_a[5]; ii++) {// debug
							//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
							//if (ii % 20 == 0) getchar();
							//}

							prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4], flag, error_approx_fine4, error_approx_coarse4, n_a[4], n_a[5]);

							// correction
							for (integer ii = 1; ii <= n_a[4]; ii++) {
								error_approx_coarse3[ii] += error_approx_fine4[ii];
							}

							// free
							//delete[] error_approx_fine4;
							//delete[] error_approx_coarse4;
							//delete[] residual_coarse4;
							//delete[] residual_fine4;

						}
						// post smothing
						for (integer iter = 0; iter < nu2; iter++) {
							//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + nnz_a[4], error_approx_coarse3, residual_coarse3, flag, n_a[4]);
							seidelq(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3]);
						}

						// prolongation
						// residual_r
						//Real *error_approx_fine3 = new Real[n_a[3] + 1];
						for (integer ii = 1; ii <= n_a[3]; ii++) {
							error_approx_fine3[ii] = 0.0;
						}

						//for (integer ii = 1; ii <= n_a[4]; ii++) {// debug
						//printf("error_approx_coarse3[%d]=%e\n",ii, error_approx_coarse3[ii]);

						//printf("residual_coarse3[%d]=%e\n", ii, residual_coarse3[ii]);
						//getchar();
						//}
						//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]+ 2 * nnz_a[3]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ 2 * nnz_a[3]+ nnz_a[4]; ii++) {// deug
						//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
						//if (ii % 20 == 0) getchar();
						//}

						prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], flag, error_approx_fine3, error_approx_coarse3, n_a[3], n_a[4]);

						// correction
						for (integer ii = 1; ii <= n_a[3]; ii++) {
							error_approx_coarse2[ii] += error_approx_fine3[ii];
						}

						// free
						//delete[] error_approx_fine3;
						//delete[] error_approx_coarse3;
						//delete[] residual_coarse3;
						//delete[] residual_fine3;

					}
					// post smothing
					for (integer iter = 0; iter < nu2; iter++) {
						//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + nnz_a[3], error_approx_coarse2, residual_coarse2, flag, n_a[3]);
						seidelq(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2]);

					}

					// prolongation
					// residual_r
					//Real *error_approx_fine2 = new Real[n_a[2] + 1];
					for (integer ii = 1; ii <= n_a[2]; ii++) {
						error_approx_fine2[ii] = 0.0;
					}

					//for (integer ii = 1; ii <= n_a[3]; ii++) {// deug
					//printf("error_approx_coarse2[%d]=%e\n",ii, error_approx_coarse2[ii]);

					//printf("residual_coarse2[%d]=%e\n", ii, residual_coarse2[ii]);
					//getchar();
					//}
					//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]+ 2 * nnz_a[2]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2]+ nnz_a[3]; ii++) {// deug
					//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
					//if (ii % 20 == 0) getchar();
					//}

					prolongation(P, 1 + nnz_aRP[0] + nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], flag, error_approx_fine2, error_approx_coarse2, n_a[2], n_a[3]);

					// correction
					for (integer ii = 1; ii <= n_a[2]; ii++) {
						error_approx_coarse1[ii] += error_approx_fine2[ii];
					}

					// free
					//delete[] error_approx_fine2;
					//delete[] error_approx_coarse2;
					//delete[] residual_coarse2;
					//delete[] residual_fine2;

				}
				// post smothing
				for (integer iter = 0; iter < nu2; iter++) {
					//seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1], 2 * nnz_a[0] + 2 * nnz_a[1] + nnz_a[2], error_approx_coarse1, residual_coarse1, flag, n_a[2]);
					seidelq(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0] + n_a[1]);
				}

				// prolongation
				// residual_r
				//Real *error_approx_fine1 = new Real[n_a[1] + 1];
				for (integer ii = 1; ii <= n_a[1]; ii++) {
					error_approx_fine1[ii] = 0.0;
				}

				//for (integer ii = 1; ii <= n_a[2]; ii++) {// deug
				//printf("error_approx_coarse1[%d]=%e\n",ii, error_approx_coarse1[ii]);

				//printf("residual_coarse1[%d]=%e\n", ii, residual_coarse1[ii]);
				//getchar();
				//}
				//for (integer ii = 1 + 2 * nnz_a[0] + 2 * nnz_a[1]; ii <= 2 * nnz_a[0] + 2 * nnz_a[1] + nnz_a[2]; ii++) {// deug
				//printf("A[%d].aij=%e, A[%d].i=%d, A[%d].j=%d\n", ii, A[ii].aij, ii, A[ii].i, ii, A[ii].j);
				//if (ii % 20 == 0) getchar();
				//}

				prolongation(P, 1 + nnz_aRP[0], nnz_aRP[0] + nnz_aRP[1], flag, error_approx_fine1, error_approx_coarse1, n_a[1], n_a[2]);

				// correction
				for (integer ii = 1; ii <= n_a[1]; ii++) {
					error_approx_coarse[ii] += error_approx_fine1[ii];
				}

				// free
				//delete[] error_approx_fine1;
				//delete[] error_approx_coarse1;
				//delete[] residual_coarse1;
				//delete[] residual_fine1;

			}

			// post smothing
			for (integer iter = 0; iter < nu2; iter++) {
				//seidel(A, 1 + 2 * nnz_a[0], 2 * nnz_a[0] + nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1]);
				seidelq(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0]);
			}

			// prolongation
			// residual_r
			//Real *error_approx_fine = new Real[n_a[0] + 1];
			for (integer ii = 1; ii <= n_a[0]; ii++) {
				error_approx_fine[ii] = 0.0;
			}

			prolongation(P, 1, nnz_aRP[0], flag, error_approx_fine, error_approx_coarse, n_a[0], n_a[1]);

			// correction
			for (integer ii = 1; ii <= n_a[0]; ii++) {
				x[ii] += error_approx_fine[ii];
			}

			// free
			//delete[] error_approx_fine;
			//delete[] error_approx_coarse;
			//delete[] residual_coarse;
			//delete[] residual_fine;
		}
		// post smother
		for (integer iter = 0; iter < nu2; iter++) {
			//seidel(A, 1, nnz_a[0], x, b, flag, n_a[0]);
			//quick seidel
			seidelq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);
		}

		//system("pause");
	}

	if (debug_reshime) system("pause");
	//system("pause");

	// free
	delete[] error_approx_fine;
	if (ilevel > 1) {
		delete[] error_approx_coarse;
		delete[] residual_coarse;
		if (ilevel > 2) {
			// free
			delete[] error_approx_fine1;
			delete[] error_approx_coarse1;
			delete[] residual_coarse1;
			delete[] residual_fine1;
			if (ilevel > 3) {
				// free
				delete[] error_approx_fine2;
				delete[] error_approx_coarse2;
				delete[] residual_coarse2;
				delete[] residual_fine2;
				if (ilevel > 4) {
					// free
					delete[] error_approx_fine3;
					delete[] error_approx_coarse3;
					delete[] residual_coarse3;
					delete[] residual_fine3;
					if (ilevel > 5) {
						// free
						delete[] error_approx_fine4;
						delete[] error_approx_coarse4;
						delete[] residual_coarse4;
						delete[] residual_fine4;
						if (ilevel > 6) {
							// free
							delete[] error_approx_fine5;
							delete[] error_approx_coarse5;
							delete[] residual_coarse5;
							delete[] residual_fine5;
							if (ilevel > 7) {
								// free
								delete[] error_approx_fine6;
								delete[] error_approx_coarse6;
								delete[] residual_coarse6;
								delete[] residual_fine6;
							}
						}
					}
				}
			}
		}
	}
	

	delete[] residual_fine;

	delete[] row_ptr_start;
	delete[] row_ptr_end;

	
	delete[] flag_shadow;
	delete[] flag;
	//delete[] flag_;
	return 0;

} // classic_aglomerative_amg1

// 21 окт 2015. Нужно реализовать Писсанецки и Густавсона, поиск можно сделать через хеш.
// 6 september 2015 кажется заработало.
integer _tmain(integer argc, _TCHAR* argv[])
{

	// Теорема : band_size в результате построения вложенных операторов может только уменьшиться.
	// Доказательство требуется.

	// Тестовая обвязка.
	// недостаток данного тестирования отсутствие условий Неймана, условий задания ненулевого теплового потока.
	// Отсутствие сборки на примере матрицы полученной аппроксимацией конвективной задачи с одной из противопоточных схем.
	// 14 сентября 2015 года реальное тестирование на примере расчёта теплового сопротивления транзисторной структуры, 
	// вставляем код в AliceFlow.

	// Условия Неймана действительно сильно замедляют сходимость. Пока можно сказать что константа 1/16  не подходит и с ней 
	// алгоритм расходится.


	// 30min 25s 40x40x40. (7min 50s Писанецки) (1min 18s , 3min 32s Густавсон) (1min 04s Густавсон от сети)
	// 6min 201x201 5.40 lite and quick  15 it
	// 201x201  37.132, 35.761
	// Густавсон ускорение 50.0% по сравнению с версией 0.04.
	// 1.01min 121x121 v0_03time=29.983, v0_04time=20.116; PGO=19.695. (v0_05=12.136 Густавсон).

    // m=81; 9.3Mb.
	// m=101; 14.3Mb.
	// speed up 38%
	// Идея Писсанецки ускорила код на 15% от версии 0.04.
	// time cube (40) 33.54s band_size ON.
	// Сравнение быстродействия cube(40) от сети:
	// без гипотезы локальности : 1min 12s. (11 iteracii)
	// с гипотезой локальности (band_size ON) : 
	integer m = 1000; // 37.6с     (Слияние 2.24-2.29min) (2.02-2.07 Писсанецки).
	Real h = 1.0/ (m - 1);
	Real h2 = h*h;

	// 4040

	Real theta = 0.25; // Контроль числа сильных связей между переменными.

	const integer dim_2D = 1;
	integer band_size = -1; // нет информации о ленте.

	switch (dim_2D) {
	case 1: theta = 0.25; break; // 2D
	case 0: theta = 0.5; break; // 3D
	}

	Ak* A=NULL;
	Ak* Atemp = NULL;
	Ak* Atemp2 = NULL;
	Ak* R=NULL;
	Ak* P=NULL;
	integer n, nnz;
	Real* x = NULL;
	Real* b = NULL;

	if (dim_2D) {
		// Finite Volume Method matrix assemble. 6 september 2015

		printf("%dx%d\n", m, m);


		band_size = 5*m; // шаблон пятиточечный.
		 n = m*m;
		 nnz = 5 * ((m - 2)*(m - 2)) + 4 + 4 * (m - 2);
		 printf("nnz=%d\n",nnz);
		 //getchar();

		// 31 aug 2015
		A = new Ak[8 * nnz + 1];
		Atemp = new Ak[3* nnz + 1];
		Atemp2 = new Ak[3*nnz + 1];
		R = new Ak[10 * n + 1];
		P = new Ak[10 * n + 1];
		// 19684 3*n+1.
		//printf("%d",3*n+1);
		//getchar();
		integer ic = 1;
		for (integer i = 1; i <= m; i++) {
			for (integer j = 1; j <= m; j++) {
				if ((i > 1) && (j > 1) && (i < m) && (j < m)) {
					//A[ic].aij = 4.0 / h2;
					A[ic].aij = 4.0; // h/h
					A[ic].i = (i - 1)*m + j;
					A[ic].j = (i - 1)*m + j;
					ic++;
					//A[ic].aij = -1.0 / h2;
					A[ic].aij = -1.0;
					A[ic].i = (i - 1)*m + j;
					A[ic].j = (i - 1)*m + j + 1;
					ic++;
					//A[ic].aij = -1.0 / h2;
					A[ic].aij = -1.0;
					A[ic].i = (i - 1)*m + j;
					A[ic].j = (i - 1)*m + j - 1;
					ic++;
					//A[ic].aij = -1.0 / h2;
					A[ic].aij = -1.0;
					A[ic].i = (i - 1)*m + j;
					A[ic].j = (i - 1)*m + j + m;
					ic++;
					//A[ic].aij = -1.0 / h2;
					A[ic].aij = -1.0;
					A[ic].i = (i - 1)*m + j;
					A[ic].j = (i - 1)*m + j - m;
					ic++;
				}
				if (((i == m) && (j == m)) || ((i == m) && (j == 1)) || ((i == 1) && (j == m)) || ((i == 1) && (j == 1)) || ((i == 1) && (j > 1) && (j < m)) || ((i == m) && (j>1) && (j < m)) || ((j == 1) && (i>1) && (i < m)) || ((j == m) && (i>1) && (i < m))) {
					A[ic].aij = 1.0;
					A[ic].i = (i - 1)*m + j;
					A[ic].j = (i - 1)*m + j;
					ic++;

				}
			}
		}

		printf("nnz=%d ic-1=%d\n", nnz, ic - 1);


		x = new Real[n+1]; // решение.
		for (integer i = 0; i <= n; i++) x[i] = 0.0;
		b = new Real[n+1]; // Правая часть.

		ic = 1;
		for (integer i = 1; i <= m; i++) {
			for (integer j = 1; j <= m; j++) {
				if ((i > 1) && (j > 1) && (i < m) && (j < m)) {
					b[ic++] = 8.0*3.141*3.141*sin(2 * 3.141*(i - 1)*h)*sin(2 * 3.141*(j - 1)*h)*h2;
				}
				if (((i == m) && (j == m)) || ((i == m) && (j == 1)) || ((i == 1) && (j == m)) || ((i == 1) && (j == 1)) || ((i == 1) && (j>1) && (j < m)) || ((i == m) && (j>1) && (j < m)) || ((j == 1) && (i>1) && (i < m)) || ((j == m) && (i>1) && (i < m))) {
					b[ic++] = sin(2 * 3.141*(i - 1)*h)*sin(2 * 3.141*(j - 1)*h);
				}
			}
		}

	}
	else {
		// 6 сентября пляж РИМИНИ

		// Finite volume method matrix assemble. 
		// volume = h2*h; Square=h2; delta_x=delta_y=delta_z=h;
		// диагональный член положителен остальные отрицательны.

		// 3D
		band_size = 7*m*m; // шаблон семиточечный.
		 n = m*m*m;
		 nnz = 7 * ((m - 2)*(m - 2)*(m - 2)) + 8 + 6 * (m - 2)*(m-2) + 12*(m-2);
		 // куб : 8 вершин, 12 рёбер, 6 граней.

		 printf("%dx%dx%d\n",m,m,m);
		 printf("nnz=%d\n",nnz);

		// 31 aug 2015
		 // 6 3 3 с запасом. реально 4.52 2.26 2.26
		 // real size 22.6 on thermal resistance.
		A = new Ak[(integer)(13 * nnz) + 1]; // 6
		// Иногда на следующем уровне вложенности число ненулевых элементов больше чем 
		// в начальной матрице.
		// real size 9.4 for resistor thermal resistor.
		// 
		Atemp = new Ak[3* nnz + 1];
		// real size 9.4 for resistor thermal resistor.
		Atemp2 = new Ak[3*nnz + 1];
		R = new Ak[(integer)(10 * n) + 1]; // 3*nnz 2.4
		P = new Ak[(integer)(10 * n) + 1]; // 3*nnz 2.4
		integer ic = 1;
		for (integer i = 1; i <= m; i++) {
			for (integer j = 1; j <= m; j++) {
				for (integer k = 1; k <= m; k++) {
					if ((i > 1) && (j > 1) && (k > 1) && (i < m) && (j < m) &&( k < m)) {
						A[ic].aij = 6.0 *h;
						A[ic].i = (k-1)*m*m + (i - 1)*m + j;
						A[ic].j = (k - 1)*m*m+ (i - 1)*m + j;
						ic++;
						A[ic].aij = -1.0 *h;
						A[ic].i = (k - 1)*m*m + (i - 1)*m + j;
						A[ic].j = (k - 1)*m*m + (i - 1)*m + j + 1;
						ic++;
						A[ic].aij = -1.0 *h;
						A[ic].i = (k - 1)*m*m + (i - 1)*m + j;
						A[ic].j = (k - 1)*m*m + (i - 1)*m + j - 1;
						ic++;
						A[ic].aij = -1.0 *h;
						A[ic].i = (k - 1)*m*m + (i - 1)*m + j;
						A[ic].j = (k - 1)*m*m + (i - 1)*m + j + m;
						ic++;
						A[ic].aij = -1.0 *h;
						A[ic].i = (k - 1)*m*m + (i - 1)*m + j;
						A[ic].j = (k - 1)*m*m + (i - 1)*m + j - m;
						ic++;
						A[ic].aij = -1.0 *h;
						A[ic].i = (k - 1)*m*m + (i - 1)*m + j;
						A[ic].j = (k - 1)*m*m + (i - 1)*m + j + m*m;
						ic++;
						A[ic].aij = -1.0 *h;
						A[ic].i = (k - 1)*m*m + (i - 1)*m + j;
						A[ic].j = (k - 1)*m*m + (i - 1)*m + j - m*m;
						ic++;
					}
					else {
					//if (((i == m) && (j == m)) || ((i == m) && (j == 1)) || ((i == 1) && (j == m)) || ((i == 1) && (j == 1)) || ((i == 1) && (j > 1) && (j < m)) || ((i == m) && (j>1) && (j < m)) || ((j == 1) && (i>1) && (i < m)) || ((j == m) && (i>1) && (i < m))) {
						A[ic].aij = 1.0;
						A[ic].i = (k - 1)*m*m + (i - 1)*m + j;
						A[ic].j = (k - 1)*m*m + (i - 1)*m + j;
						ic++;

					}
				}
			}
		}

		printf("nnz=%d ic-1=%d\n",nnz,ic-1);


		x = new Real[n+1]; // решение.
		for (integer i = 0; i <= n; i++) x[i] = 0.0; // инициалзация.
		b = new Real[n+1]; // Правая часть.

		ic = 1;
		for (integer i = 1; i <= m; i++) {
			for (integer j = 1; j <= m; j++) {
				for (integer k = 1; k <= m; k++) {
					if ((i > 1) && (j > 1)&& (k > 1) && (i < m) && (j < m) && (k < m)) {
						b[ic++] = 16.0*3.141*3.141*3.141*sin(2 * 3.141*(i - 1)*h)*sin(2 * 3.141*(j - 1)*h)*sin(2 * 3.141*(k - 1)*h)*h2*h;
					}
					else {
					//if (((i == m) && (j == m)) || ((i == m) && (j == 1)) || ((i == 1) && (j == m)) || ((i == 1) && (j == 1)) || ((i == 1) && (j>1) && (j < m)) || ((i == m) && (j>1) && (j < m)) || ((j == 1) && (i>1) && (i < m)) || ((j == m) && (i>1) && (i < m))) {
						b[ic++] = sin(2 * 3.141*(k - 1)*h)*sin(2 * 3.141*(i - 1)*h)*sin(2 * 3.141*(j - 1)*h);
					}
				}
			}
		}
	}

	classic_aglomerative_amg1(A, nnz, n, R, P, Atemp, Atemp2, x, b, theta);

	// m in 3D memorysize mb
	// 21  14.8Mb
	// 31  49.7Mb
	// 41  119.4Mb

	//getchar();
	delete[] A;
	delete[] R;
	delete[] P;
	delete[] Atemp;
	delete[] Atemp2;
	delete[] x;
	delete[] b;
	
	// 2D m=81 debug от акумулятора. 
	// nu1=1 nu2=0 47.904
	// nu1=1 nu2=1 34.548
	// nu1=2 nu2=1 29.465
	// nu1=2 nu2=2 27.76
	// 2D m=81 realese от акумулятора.
	// nu1=2 nu2=1 8.898
	// nu1=2 nu2=2 8.109
	// nu1=3 nu2=2 7.247
	// nu1=3 nu2=3 6.962
	// nu1=4 nu2=3 6.555 6.799
	// seidelq 1 level
	// nu1=4 nu2=3 5.691 5.605 
	  // от сети.
	  // seidelq 2 level
	  // nu1=4 nu2=3 1.121 1.9 1.672
	  // 1.9 792 0.982
	  // 1.174 273 0.94248 // оставлен residual на втором левеле а не residualq.
	// от аккумулятора.
    // nu1=4 nu2=3 3.412  оставлен residual на втором левеле а не residualq.
	// nu1=4 nu2=3 5.503 residualq.
	// оставлен residual на втором левеле а не residualq.
	// m=81 12 level
	// nu1=4 nu2=3 2.916 3.081 2.899 2.89
	// nu1=4 nu2=3 2.487  (третий патч 1.715)


	// от акумулятора 1.091
	//(msvcr120.dll 3.94 2.06 2.33 3.1) 1.234 1.072
	// 2.893 2.788 2.892 2.8 2.814
	// m=121 от акумулятора. 13 level.
	// 2.556 2.546 2.627 2.533 2.602 2.501 2.542 2.487 2.589 5.131
	// 1.686 1.742 от сети. 1.763
	// m=221 2D 14 level
	// 33.857 32.649 (6.809s 3 патча. 6.984).


	// m=24 3D от сети
	// без патча 957мс.
	// без патча m=48 8.713c. с патчем 5.7s. 110тыс узлов.

	// etalon 27s. или 36с. от сети. etalon 47s от аккумулятора.
	// m=100 3D 1 million nodes
	// 15 level 117 iteration.
	// 2 min 42s. в 6 раз медленнее.
	// 47.16s от сети. 47.276
	//system("pause");
	return 0;
}


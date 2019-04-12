// aggregativeAMG.cpp: определяет точку входа для консольного приложения.
//

//#include "stdafx.h"
//#include <math.h>
//#include <stdlib.h>
//#include <windows.h>
//#include <omp.h>

//#define doublereal double
//#define integer int

typedef struct TAk {
	integer i, j;
	doublereal aij;
	integer ind; // позиция в первоначальной сортировке.
} Ak;

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

  // Вот алгоритм PivotList
integer PivotList(Ak * &A, integer n, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue = A[first].i*n+A[first].j;
	integer PivotPointeger = first;

	for (integer index = (first + 1); index <= last; index++) {
		if (A[index].i*n+A[index].j<PivotValue) {
			PivotPoint++;
			swap(A, PivotPoint, index);
		}
	}

	swap(A, first, PivotPoint);

	return PivotPoint;
} // PivotListamg

  // Вот алгоритм PivotList
integer PivotList_j(Ak * &A, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue_j = A[first].j;
	integer PivotValue_i = A[first].i;
	integer PivotPointeger = first;

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
void HeapSort_j(Ak * &A, integer n, integer first, integer last)
{

	Ak maxelm; // элемент с наибольшим значением ключа

	integer iadd = first - 1;
			   // конструирование пирамиды
	for (integer i =  ((last - first + 1 )/ 2); i >= 1; i--)
	{
		FixHeap_j(A, i, A[i+iadd], last - first + 1, n, iadd);
	}
	for (integer i = last-first+1; i >= 2; i--)
	{
		// скопировать корень пирамиды в список
		// переформировать пирамиду
		maxelm = A[1+iadd];
		FixHeap_j(A, 1, A[i+iadd], i - 1, n, iadd);
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
	Ak* Bm = Aorig + size;
	Ak* C;

	// A - сортируемый массив, B - вспомогательная память (справа от основных данных в А такого-же размера что и основные данные в А).
	// C - указатель для обмена.
	for (integer i = 1; i < size; i = i * 2) // размер объединяемых фрагментов
	{
		for (integer j = 0; j < size; j = j + 2 * i) // начало первого из объединяемых 
			// фрагментов
		{
			integer r = j + i; // начало второго из объединяемых фрагментов
			integer n1=0, n2=0;
			if (i < size - j) { n1 = i; }
			else { n1 = size - j; };
			if (i < size - r) { n2 = i; }
			else { n2 = size - r; };

			if (n1 < 0) n1 = 0;
			if (n2 < 0) n2 = 0;

			// слияние упорядоченных фрагментов
			for (integer ia = 0, ib = 0, k = 0; k < n1 + n2; k++)
			{
				if (ia >= n1) Bm[j + k] = A[r + ib++];
				else
					if (ib >= n2) Bm[j + k] = A[j + ia++];
					else {
						bool compare_result = false;
						integer lCadd = j + ia;
						integer lCadd1 = r + ib;
						if (A[lCadd1].i > A[lCadd].i) {
							compare_result = true;
						}
						else if (A[lCadd1].i == A[lCadd].i) {
							if (A[lCadd1].j > A[lCadd].j) {
								compare_result = true;
							}
						}
						if (compare_result) {
							Bm[j + k] = A[j + ia++];
						}
						else {
							Bm[j + k] = A[r+ib++];
						}
					}
			}



		}
		C = A; A = Bm; Bm = C;

	}

	C = A; A = Bm; Bm = C;

	// Копирование, если результат размещен не в основном а в вспомогательном массиве
	if (Bm != Aorig)
		memcpy(Aorig, Bm, size*sizeof(Ak));

	A = NULL; Bm = NULL; C = NULL;
} // MergeSort


// Сортировка слиянием. 
// Требует дополнительной пямяти.
void MergeSort_j(Ak * &Aorig, integer size) {
	// предполагается индексация от нуля до size-1.
	// Массив А предполагается не менее двойного размера.
	// Двойная память это недостаток данного алгоритма.
	Ak* A = Aorig;
	Ak* Bm = A + size;
	Ak* C;

	// A - сортируемый массив, B - вспомогательная память (справа от основных данных в А такого-же размера что и основные данные в А).
	// C - указатель для обмена.
	for (integer i = 1; i < size; i = i * 2) // размер объединяемых фрагментов
	{
		for (integer j = 0; j < size; j = j + 2 * i) // начало первого из объединяемых 
			// фрагментов
		{
			integer r = j + i; // начало второго из объединяемых фрагментов
			integer n1 = 0, n2 = 0;
			if (i < size - j) { n1 = i; }
			else { n1 = size - j; };
			if (i < size - r) { n2 = i; }
			else { n2 = size - r; };

			if (n1 < 0) n1 = 0;
			if (n2 < 0) n2 = 0;

			// слияние упорядоченных фрагментов
			for (integer ia = 0, ib = 0, k = 0; k < n1 + n2; k++)
			{
				if (ia >= n1) Bm[j + k] = A[r + ib++];
				else
					if (ib >= n2) Bm[j + k] = A[j + ia++];
					else {
						bool compare_result = false;
						integer lCadd = j + ia;
						integer lCadd1 = r + ib;
						if (A[lCadd1].j > A[lCadd].j) {
							compare_result = true;
						}
						else if (A[lCadd1].j == A[lCadd].j) {
							if (A[lCadd1].i > A[lCadd].i) {
								compare_result = true;
							}
						}
						if (compare_result) {
							Bm[j + k] = A[j + ia++];
						}
						else {
							Bm[j + k] = A[r + ib++];
						}
					}
			}



		}
		C = A; A = Bm; Bm = C;

	}

	C = A; A = Bm; Bm = C;

	// Копирование, если результат размещен не в основном а в вспомогательном массиве
	if (Bm != Aorig)
		memcpy(Aorig, Bm, size*sizeof(Ak));

	A = NULL; Bm = NULL; C = NULL;
} // MergeSort_j



  // Быстрая сортировка Хоара.
  // Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
  // стр. 106.
void QuickSort(Ak * &A, integer n, integer first, integer last) {
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
				if (A[i].i*n+A[i].j>A[i + 1].i*n+A[i+1].j) {
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
			pivot = PivotList(A, n, first, last);
			QuickSort(A, n, first, pivot - 1);
			QuickSort(A, n, pivot + 1, last);
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
			pivot = PivotList_j(A, first, last);
			QuickSort_j(A, first, pivot - 1);
			QuickSort_j(A, pivot + 1, last);
		}
	}
} // QuickSort_j

/*
// Это исторически первоначальный код содержащий лишь построение последовательности 
// вложенных матриц и не иодержащий построения операторов restriction и prolongation. 
// создание этого кода завершено 1 сентября 2015 года.

integer aggregative_amg(Ak* &A, integer nnz, integer n, doublereal* &x, doublereal* &b) {


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
void seidel(Ak* &A, integer istart, integer iend, doublereal* &x, doublereal* &b, bool* &flag, integer n)
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
			doublereal ap = 0.0;
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
		doublereal ap = 0.0;
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
void seidelqSOR(Ak* &A, integer istartq, integer iendq, doublereal* &x, doublereal* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// Сделаем метод последовательной верхней релаксации.
	doublereal omega = 1.855;
	// Верхняя итерация на каждом из уровней.

	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		doublereal xnew = b[istr];
		//x[istr] = b[istr];

		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++)
		{
			xnew += -A[ii1].aij*x[A[ii1].j];
			//x[istr] += -A[ii1].aij*x[A[ii1].j];
		}
		//x[istr] *= A[row_ptr_start[ii]].aij;
		x[istr] = (1.0 - omega)*x[istr] + omega*(A[row_ptr_start[ii]].aij*xnew);
	}


} // seidelq


// smoother.
// 9 september 2015.
// q - quick.
void seidelq(Ak* &A, integer istartq, integer iendq, doublereal* &x, doublereal* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
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
void seidelq(Ak* &A, integer istartq, integer iendq, doublereal* &x, doublereal* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	//doublereal *sum=new doublereal[endpos-startpos+1];

#pragma omp parallel for
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		x[istr] = b[istr];
		//	doublereal sum = 0.0;
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
void residual(Ak* &A, integer istart, integer iend, doublereal* &x, doublereal* &b, bool* &flag, integer n, doublereal* &residual)
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
void residualq(Ak* &A, integer istartq, integer iendq, doublereal* &x, doublereal* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, doublereal* &residual)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		// r=b+sum(anb*Fnb)-ap*Fp; ap, anb>0.
		residual[istr] = b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++) { residual[istr] += -A[ii1].aij*x[A[ii1].j]; }
		residual[istr] += (-1.0/A[row_ptr_start[ii]].aij)*x[istr];
	}


} // residualq

// smoother.
// 9 september 2015.
// q - quick.
void residualqspeshial(Ak* &A, integer istartq, integer iendq, doublereal* &x, doublereal* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, doublereal* &residual)
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
		doublereal omega = 1.855; // SOR
		//residual[istr] += ((-1.0 / A[row_ptr_start[ii]].aij)*x[istr]); // верный вариант.
		residual[istr] += ((-1.0 / A[row_ptr_start[ii]].aij)*x[istr])/7.0;
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
void restriction(Ak* &R, integer istart, integer iend, bool* &flag, doublereal * &x_fine, doublereal * &x_coarse, integer n_fine, integer n_coarse) {
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
void prolongation(Ak* &P, integer istart, integer iend, bool* &flag, doublereal * &x_fine, doublereal * &x_coarse, integer n_fine, integer n_coarse) {
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
doublereal norma(doublereal * &r, integer n) {
	doublereal ret = 0.0;
	for (integer ii = 1; ii <= n; ii++) {
		ret += r[ii] * r[ii] / n;
	}
	ret = sqrt(ret);
	return ret;
}

// экспорт полевой величины u в программу tecplot 360.
void exporttecplot(doublereal* u, integer n_size) {
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
			doublereal h = 1.0 / (m - 1);
			fprintf_s(fp, "I=%d, J=%d, K=1, F=TOCHKA\n", m, n);
			for (integer j = 0; j < n; j++) for (integer i = 0; i < m; i++)   fprintf_s(fp, "%e %e %e\n", i*h, j*h, u[i*m + j + 1]);
			fclose(fp);
			//-->WinExec("C:\\Program Files\\Tecplot\\Tecplot 360 EX 2014 R1\\bin\\tec360.exe fedorenko1.PLT", SW_NORMAL);
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
		else return mid;
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
		else return mid;
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
		else return mid;
	}
	return -1;
}

// 3 september 2015 Villa Borgese.
integer aggregative_amg(Ak* &A, integer nnz, integer n,
	                Ak* &R, // restriction
	                Ak* &P, // prolongation
	                doublereal* &x, doublereal* &b)
{

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
			QuickSort(A, n_a[ilevel - 1], 1 + iadd, nnz_a[ilevel - 1] + iadd);
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
			HeapSort_j(A, n_a[ilevel - 1], 1 + iadd + 2 * nnz_a[ilevel - 1], nnz_a[ilevel - 1] + iadd + 2 * nnz_a[ilevel - 1]);
		}
		else {
			Ak* Aorig = &A[1 + iadd + 2 * nnz_a[ilevel - 1]];
			MergeSort_j(Aorig, nnz_a[ilevel - 1]);
			Aorig = NULL;
		}


		integer n_coarce = 1; // номер агрегата.
		nnzR = 1;
		const integer max_sosed = 37850;
		const integer NULL_SOSED = -1;
		integer vacant = NULL_SOSED;
		for (integer ii = 1 + iadd; ii <= nnz_a[ilevel - 1] + iadd; ii++) {
			if (flag[A[ii].i] == false) {
				// Вычисляем по немодифиуцированной матрице А (хранящейся слева).

				doublereal sum = 0.0;
				integer nnzRl = nnzR + iaddR;

				integer set[max_sosed]; // не более 20 узлов в одном агрегате.
				// инициализация убрана потомучто она не нужна и она сильно тормозит быстродействие.
				//for (integer js = 0; js < max_sosed; js++) {
					//set[js] = NULL_SOSED;
				//}
				integer ic = 0;
				set[ic] = A[ii].i;
				doublereal theta = 0.25; // контроль числа сильных связей между переменными.
				doublereal max_vnediagonal = -1.0; // максимальное значение модуля вне диагонального элемента. 
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
				
				/*
				{
					integer* Aset = new integer[ic];
					for (integer ii3 = 0; ii3 < ic; ii3++) Aset[ii3] = set[ii3]; // copy
					// HeapSort(Aset,0,ic-1);
					quickSort_set(Aset, 0, ic - 1);

					// bynarySearh 88.56%
					// aggregativeamg 8.89%
					// seidel 0.87%
					// (i,j) -> (I,J)
					// модифицируем копию A находящуюся справа.
					for (integer k = nnz_a[ilevel - 1] + 1 + iadd; k <= 2 * nnz_a[ilevel - 1] + iadd; k++) {
						bool found = false;

						if (BinarySearch(Aset, A[k - nnz_a[ilevel - 1]].i, ic - 1) > -1) found = true;
						//for (integer k1 = 0; k1 < ic; k1++) {
						//	if (A[k - nnz_a[ilevel - 1]].i == set[k1]) found = true;
						//	}
						if (found) A[k].i = n_coarce;
						found = false;
						if (BinarySearch(Aset, A[k - nnz_a[ilevel - 1]].j, ic - 1) > -1) found = true;
						
						//for (integer k1 = 0; k1 < ic; k1++) {
						//if (A[k - nnz_a[ilevel - 1]].j == set[k1]) found = true;
						//}
						
						if (found) A[k].j = n_coarce;

					}

					delete[] Aset;
				}
				*/
				
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
		    HeapSort_j(P, n_a[ilevel-1], 1 + iaddR, iaddR + nnzR - 1);
		}

		// оператор интерполяции это не просто транспонированный оператор проекции а
		// а транспонированный оператор проекции умноженный на константу. Константа 
		// определяется из следующего соображения : если сумма элементов оператора рестрикции в стоке единица,
		// то соответственно в столбце у оператора интерполяции максимальный элемент равен единица.
		// этот код обязательно должен быть включён чтобы пара рестрикция-интерполяция была верна.
		// действительо рабочий 15 сентября 2015.
		for (integer ii = 1; ii <= n; ii++) {
			flag_[ii] = false;
		}
		doublereal mul = -1.e30;
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
		
		
		nnz_aRP[ilevel - 1] = nnzR - 1;
		iaddR += nnzR - 1;

		//printf("%d %d\n",n,n_coarce-1);
		//getchar();
		n_a[ilevel] = n_coarce - 1;

		// сортировка по новому ключу key=i*(iglcoarce_number-1)+j;
		// в позиции ind сохранён индекс предыдущей позиции.
		//heapsort(A, key = i*(n_coarce - 1) + j, nnz + 1, 2 * nnz);
		if (bquicktypesort) {
			QuickSort(A, n_coarce - 1, nnz_a[ilevel - 1] + 1 + iadd, 2 * nnz_a[ilevel - 1] + iadd);
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

/*
	//exporttecplot(b,n);

	doublereal *test_coarse = new doublereal[n_a[1] + 1];

	// restriction
	restriction(R, 1, nnz_aRP[0], flag, b, test_coarse, n_a[0], n_a[1]);
	for (integer ii = 1; ii <= n_a[0]; ii++) {
		b[ii] = 0.0;
	}

	{
		doublereal *test_coarse1 = new doublereal[n_a[2] + 1];

		// restriction
		restriction(R, 1+nnz_aRP[0], nnz_aRP[0]+nnz_aRP[1], flag, test_coarse, test_coarse1, n_a[1], n_a[2]);
		for (integer ii = 1; ii <= n_a[1]; ii++) {
			test_coarse[ii] = 0.0;
		}

		{
			doublereal *test_coarse2 = new doublereal[n_a[3] + 1];

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
void seidel(Ak* &A, integer istart, integer iend, doublereal* &x, doublereal* &b, bool* &flag, integer n)
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
			doublereal ap = 0.0;
			x[istr] = b[istr];
			while ((ic<=iend)&&(A[ic].i == istr)) {
				if (A[ic].j != istr) {
					x[istr] += -A[ic].aij*x[A[ic].j];
				}
				else ap = A[ic].aij;
				ic++;
			}
			if (fabs(ap) < 1.0e-30) {
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
	integer *row_ptr_start = new integer[4 * n_a[0]+1];
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
			doublereal ap = 0.0;
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
				doublereal ap = 0.0;
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
				doublereal ap = 0.0;
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
				doublereal ap = 0.0;
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
				doublereal ap = 0.0;
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
				doublereal ap = 0.0;
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
				doublereal ap = 0.0;
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


	/*
	
	// smoother.
// 9 september 2015.
// q - quick.
// seidelq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);
void seidelq(Ak* &A, integer istartq, integer iendq, doublereal* &x, doublereal* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
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
	*/

	


	printf("cycling: V cycle.\n");
	printf("level=%d\n", ilevel);
	printf("multigrid R.P.Fedorenko 1961.\n");
	printf("aggregative algebraic multigrid method.\n");
	//system("pause");

	// 10 11 21 multigrid tutorial Вильм Бригг.

	integer nu1 = 4; // 4
	integer nu2 = 3; // 3

//	ilevel = 1; //debug
	doublereal rho = 1.0;
	doublereal dres = 1.0;
	integer iiter = 1;
	const doublereal tolerance = 1.0e-12;// 1.0e-12;


	doublereal *residual_fine = new doublereal[n_a[0] + 1];
	doublereal *residual_coarse = NULL;
	doublereal* error_approx_coarse = NULL;
	doublereal *residual_fine1 = NULL;
	doublereal *residual_coarse1 = NULL;
	doublereal* error_approx_coarse1 = NULL;
	doublereal *error_approx_fine1 = NULL;
	doublereal *residual_fine2 = NULL;
	doublereal *residual_coarse2 = NULL;
	doublereal* error_approx_coarse2 = NULL;
	doublereal *error_approx_fine2 = NULL;
	doublereal *residual_fine3 = NULL;
	doublereal *residual_coarse3 = NULL;
	doublereal* error_approx_coarse3 = NULL;
	doublereal *error_approx_fine3 = NULL;
	doublereal *residual_fine4 = NULL;
	doublereal *residual_coarse4 = NULL;
	doublereal *error_approx_coarse4 = NULL;
	doublereal *error_approx_fine4 = NULL;
	doublereal *residual_fine5 = NULL;
	doublereal *residual_coarse5 = NULL;
	doublereal* error_approx_coarse5 = NULL;
	doublereal *error_approx_fine5 = NULL;
	doublereal *residual_fine6 = NULL;
	doublereal *residual_coarse6 = NULL;
	doublereal* error_approx_coarse6 = NULL;
	doublereal *error_approx_fine6 = NULL;

	if (ilevel > 1) {
		residual_coarse = new doublereal[n_a[1] + 1];
		error_approx_coarse = new doublereal[n_a[1] + 1];
		if (ilevel > 2) {
			// residual
			residual_fine1 = new doublereal[n_a[1] + 1];
			residual_coarse1 = new doublereal[n_a[2] + 1];
			error_approx_coarse1 = new doublereal[n_a[2] + 1];
			error_approx_fine1 = new doublereal[n_a[1] + 1];
			if (ilevel > 3) {
				// residual
				residual_fine2 = new doublereal[n_a[2] + 1];
				residual_coarse2 = new doublereal[n_a[3] + 1];
				error_approx_coarse2 = new doublereal[n_a[3] + 1];
				error_approx_fine2 = new doublereal[n_a[2] + 1];
				if (ilevel > 4) {
					// residual
					residual_fine3 = new doublereal[n_a[3] + 1];
					residual_coarse3 = new doublereal[n_a[4] + 1];
					error_approx_coarse3 = new doublereal[n_a[4] + 1];
					error_approx_fine3 = new doublereal[n_a[3] + 1];
					if (ilevel > 5) {
						// residual
						residual_fine4 = new doublereal[n_a[4] + 1];
						residual_coarse4 = new doublereal[n_a[5] + 1];
						error_approx_coarse4 = new doublereal[n_a[5] + 1];
						error_approx_fine4 = new doublereal[n_a[4] + 1];
						if (ilevel > 6) {
							// residual
							residual_fine5 = new doublereal[n_a[5] + 1];
							residual_coarse5 = new doublereal[n_a[6] + 1];
							error_approx_coarse5 = new doublereal[n_a[6] + 1];
							error_approx_fine5 = new doublereal[n_a[5] + 1];
							if (ilevel > 7) {
								// residual
								residual_fine6 = new doublereal[n_a[6] + 1];
								residual_coarse6 = new doublereal[n_a[7] + 1];
								error_approx_coarse6 = new doublereal[n_a[7] + 1];
								error_approx_fine6 = new doublereal[n_a[6] + 1];
							}
						}
					}
				}
			}
		}
	}
	doublereal *error_approx_fine = new doublereal[n_a[0] + 1];

	// если это здесь раскоментировать то метод станет расходящимся.
	//for (integer ii = 1; ii <= n_a[1]; ii++) {
		//error_approx_coarse[ii] = 0.0;
	//}
	//for (integer ii = 1; ii <= n_a[2]; ii++) {
		//error_approx_coarse1[ii] = 0.0;
	//}


	//ilevel = 2;

	

	//for (integer iprohod = 0; iprohod < 20; iprohod++) {
	while (dres>tolerance) {

		

		// smother
		for (integer iter = 0; iter < nu1; iter++) {
			//seidel(A, 1, nnz_a[0], x, b, flag, n_a[0]);
			//quick seidel
			seidelq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);
			//seidelqSOR(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);
		}
		
		//exporttecplot(x, n);

		// residual_r
		//doublereal *residual_fine = new doublereal[n_a[0] + 1];
		//residual(A, 1, nnz_a[0], x, b, flag, n_a[0], residual_fine);
		residualq(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0, residual_fine);
		dres = norma(residual_fine, n_a[0]);
		printf("%d %e rho=%e\n",iiter, dres, dres/rho);
		iiter++;

		//if (iiter > 6) break;


		//rho=norma(residual_fine, n_a[0]);
		rho = dres;
		//if (iprohod%5==0) getchar();
		if (ilevel > 1) {

			//doublereal *residual_coarse = new doublereal[n_a[1] + 1];

			// restriction
			restriction(R, 1, nnz_aRP[0], flag, residual_fine, residual_coarse, n_a[0], n_a[1]);

			// A*e=r;
			//doublereal* error_approx_coarse = new doublereal[n_a[1] + 1];
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
				//doublereal *residual_fine1 = new doublereal[n_a[1] + 1];
				//residual(A, 1+2*nnz_a[0], 2*nnz_a[0]+nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1], residual_fine1);
				//residualq(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0], residual_fine1);
				residualqspeshial(A, 1, n_a[1], error_approx_coarse, residual_coarse, row_ptr_start, row_ptr_end, n_a[0], residual_fine1);


				//doublereal *residual_coarse1 = new doublereal[n_a[2] + 1];

				// restriction
				restriction(R, 1+nnz_aRP[0], nnz_aRP[0]+nnz_aRP[1], flag, residual_fine1, residual_coarse1, n_a[1], n_a[2]);
			
				// A*e=r;
				//doublereal* error_approx_coarse1 = new doublereal[n_a[2] + 1];
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
					//doublereal *residual_fine2 = new doublereal[n_a[2] + 1];
					//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1], 2 * nnz_a[0] + 2 * nnz_a[1] + nnz_a[2], error_approx_coarse1, residual_coarse1, flag, n_a[2], residual_fine2);
					//residualq(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0]+n_a[1], residual_fine2);
					residualqspeshial(A, 1, n_a[2], error_approx_coarse1, residual_coarse1, row_ptr_start, row_ptr_end, n_a[0] + n_a[1], residual_fine2);

					//doublereal *residual_coarse2 = new doublereal[n_a[3] + 1];

					// restriction
					restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], flag, residual_fine2, residual_coarse2, n_a[2], n_a[3]);

					// A*e=r;
					//doublereal* error_approx_coarse2 = new doublereal[n_a[3] + 1];
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
						//doublereal *residual_fine3 = new doublereal[n_a[3] + 1];
						//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + nnz_a[3], error_approx_coarse2, residual_coarse2, flag, n_a[3], residual_fine3);
						//residualq(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2], residual_fine3);
						//speshial
						residualqspeshial(A, 1, n_a[3], error_approx_coarse2, residual_coarse2, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2], residual_fine3);



						//doublereal *residual_coarse3 = new doublereal[n_a[4] + 1];

						// restriction
						restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], flag, residual_fine3, residual_coarse3, n_a[3], n_a[4]);

						// A*e=r;
						//doublereal* error_approx_coarse3 = new doublereal[n_a[4] + 1];
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
							//doublereal *residual_fine4 = new doublereal[n_a[4] + 1];
							//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + nnz_a[4], error_approx_coarse3, residual_coarse3, flag, n_a[4], residual_fine4);
							residualq(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3], residual_fine4);
							//speshial 14 september 2015.
							//residualqspeshial(A, 1, n_a[4], error_approx_coarse3, residual_coarse3, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3], residual_fine4);


							//doublereal *residual_coarse4 = new doublereal[n_a[5] + 1];

							// restriction
							restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3]+ nnz_aRP[4], flag, residual_fine4, residual_coarse4, n_a[4], n_a[5]);

							// A*e=r;
							//doublereal* error_approx_coarse4 = new doublereal[n_a[5] + 1];
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
								//doublereal *residual_fine5 = new doublereal[n_a[5] + 1];
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

								//doublereal *residual_coarse5 = new doublereal[n_a[6] + 1];

								// restriction
								restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], flag, residual_fine5, residual_coarse5, n_a[5], n_a[6]);

								// A*e=r;
								//doublereal* error_approx_coarse5 = new doublereal[n_a[6] + 1];
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
									//doublereal *residual_fine6 = new doublereal[n_a[6] + 1];
									//residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] +2*nnz_a[5], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + nnz_a[6], error_approx_coarse5, residual_coarse5, flag, n_a[6], residual_fine6);
									residualq(A, 1, n_a[6], error_approx_coarse5, residual_coarse5, row_ptr_start, row_ptr_end, n_a[0] + n_a[1] + n_a[2] + n_a[3] + n_a[4] + n_a[5], residual_fine6);

									//doublereal *residual_coarse6 = new doublereal[n_a[7] + 1];

									// restriction
									restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5]+nnz_aRP[6], flag, residual_fine6, residual_coarse6, n_a[6], n_a[7]);

									// A*e=r;
									//doublereal* error_approx_coarse6 = new doublereal[n_a[7] + 1];
									for (integer ii = 1; ii <= n_a[7]; ii++) {
										error_approx_coarse6[ii] = 0.0;
									}
									// pre smothing
									for (integer iter = 0; iter < nu1; iter++) {
										seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + nnz_a[7], error_approx_coarse6, residual_coarse6, flag, n_a[7]);
									}

									if (ilevel > 8) {
										// residual
										doublereal *residual_fine7 = new doublereal[n_a[7] + 1];
										residual(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5]+2*nnz_a[6], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] +2*nnz_a[6]+ nnz_a[7], error_approx_coarse6, residual_coarse6, flag, n_a[7], residual_fine7);


										doublereal *residual_coarse7 = new doublereal[n_a[8] + 1];

										// restriction
										restriction(R, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6]+nnz_aRP[7], flag, residual_fine7, residual_coarse7, n_a[7], n_a[8]);

										// A*e=r;
										doublereal* error_approx_coarse7 = new doublereal[n_a[8] + 1];
										for (integer ii = 1; ii <= n_a[8]; ii++) {
											error_approx_coarse7[ii] = 0.0;
										}
										// pre smothing
										for (integer iter = 0; iter < nu1; iter++) {
											seidel(A, 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6]+2*nnz_a[7], 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] +2*nnz_a[7]+ nnz_a[8], error_approx_coarse7, residual_coarse7, flag, n_a[8]);
										}


										if (ilevel > 9) {
											// residual
											doublereal *residual_fine8 = new doublereal[n_a[8] + 1];
											integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7];
											integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + nnz_a[8];
											residual(A, n1, n2, error_approx_coarse7, residual_coarse7, flag, n_a[8], residual_fine8);


											doublereal *residual_coarse8 = new doublereal[n_a[9] + 1];

											// restriction
											integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6]+nnz_aRP[7];
											integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7]+nnz_aRP[8];
											restriction(R,n3 ,n4 , flag, residual_fine8, residual_coarse8, n_a[8], n_a[9]);

											// A*e=r;
											doublereal* error_approx_coarse8 = new doublereal[n_a[9] + 1];
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
												doublereal *residual_fine9 = new doublereal[n_a[9] + 1];
												integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7]+2*nnz_a[8];
												integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] +2*nnz_a[8]+ nnz_a[9];
												residual(A, n1, n2, error_approx_coarse8, residual_coarse8, flag, n_a[9], residual_fine9);


												doublereal *residual_coarse9 = new doublereal[n_a[10] + 1];

												// restriction
												integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7]+nnz_aRP[8];
												integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8]+nnz_aRP[9];
												restriction(R, n3, n4, flag, residual_fine9, residual_coarse9, n_a[9], n_a[10]);

												// A*e=r;
												doublereal* error_approx_coarse9 = new doublereal[n_a[10] + 1];
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
													doublereal *residual_fine10 = new doublereal[n_a[10] + 1];
													integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8]+2*nnz_a[9];
													integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] +2*nnz_a[9]+ nnz_a[10];
													residual(A, n1, n2, error_approx_coarse9, residual_coarse9, flag, n_a[10], residual_fine10);


													doublereal *residual_coarse10 = new doublereal[n_a[11] + 1];

													// restriction
													integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8]+nnz_aRP[9];
													integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9]+nnz_aRP[10];
													restriction(R, n3, n4, flag, residual_fine10, residual_coarse10, n_a[10], n_a[11]);

													// A*e=r;
													doublereal* error_approx_coarse10 = new doublereal[n_a[11] + 1];
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
														doublereal *residual_fine11 = new doublereal[n_a[11] + 1];
														integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9]+2*nnz_a[10];
														integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] +2*nnz_a[10]+ nnz_a[11];
														residual(A, n1, n2, error_approx_coarse10, residual_coarse10, flag, n_a[11], residual_fine11);


														doublereal *residual_coarse11 = new doublereal[n_a[12] + 1];

														// restriction
														integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9]+nnz_aRP[10];
														integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] +nnz_aRP[10]+ nnz_aRP[11];
														restriction(R, n3, n4, flag, residual_fine11, residual_coarse11, n_a[11], n_a[12]);

														// A*e=r;
														doublereal* error_approx_coarse11 = new doublereal[n_a[12] + 1];
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
															doublereal *residual_fine12 = new doublereal[n_a[12] + 1];
															integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11];
															integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + nnz_a[12];
															residual(A, n1, n2, error_approx_coarse11, residual_coarse11, flag, n_a[12], residual_fine12);


															doublereal *residual_coarse12 = new doublereal[n_a[13] + 1];

															// restriction
															integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11];
															integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12];
															restriction(R, n3, n4, flag, residual_fine12, residual_coarse12, n_a[12], n_a[13]);

															// A*e=r;
															doublereal* error_approx_coarse12 = new doublereal[n_a[13] + 1];
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
																doublereal *residual_fine13 = new doublereal[n_a[13] + 1];
																integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12];
																integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + nnz_a[13];
																residual(A, n1, n2, error_approx_coarse12, residual_coarse12, flag, n_a[13], residual_fine13);


																doublereal *residual_coarse13 = new doublereal[n_a[14] + 1];

																// restriction
																integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11]+nnz_aRP[12];
																integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12]+nnz_aRP[13];
																restriction(R, n3, n4, flag, residual_fine13, residual_coarse13, n_a[13], n_a[14]);

																// A*e=r;
																doublereal* error_approx_coarse13 = new doublereal[n_a[14] + 1];
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
																	doublereal *residual_fine14 = new doublereal[n_a[14] + 1];
																	integer n1 = 1 + 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13];
																	integer n2 = 2 * nnz_a[0] + 2 * nnz_a[1] + 2 * nnz_a[2] + 2 * nnz_a[3] + 2 * nnz_a[4] + 2 * nnz_a[5] + 2 * nnz_a[6] + 2 * nnz_a[7] + 2 * nnz_a[8] + 2 * nnz_a[9] + 2 * nnz_a[10] + 2 * nnz_a[11] + 2 * nnz_a[12] + 2 * nnz_a[13] + nnz_a[14];
																	residual(A, n1, n2, error_approx_coarse13, residual_coarse13, flag, n_a[14], residual_fine14);


																	doublereal *residual_coarse14 = new doublereal[n_a[15] + 1];

																	// restriction
																	integer n3 = 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13];
																	integer n4 = nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14];
																	restriction(R, n3, n4, flag, residual_fine14, residual_coarse14, n_a[14], n_a[15]);

																	// A*e=r;
																	doublereal* error_approx_coarse14 = new doublereal[n_a[15] + 1];
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
																	doublereal *error_approx_fine14 = new doublereal[n_a[14] + 1];
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
																doublereal *error_approx_fine13 = new doublereal[n_a[13] + 1];
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
															doublereal *error_approx_fine12 = new doublereal[n_a[12] + 1];
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
														doublereal *error_approx_fine11 = new doublereal[n_a[11] + 1];
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
													doublereal *error_approx_fine10 = new doublereal[n_a[10] + 1];
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
												doublereal *error_approx_fine9 = new doublereal[n_a[9] + 1];
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
											doublereal *error_approx_fine8 = new doublereal[n_a[8] + 1];
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
										doublereal *error_approx_fine7 = new doublereal[n_a[7] + 1];
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
									//doublereal *error_approx_fine6 = new doublereal[n_a[6] + 1];
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
								//doublereal *error_approx_fine5 = new doublereal[n_a[5] + 1];
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
							//doublereal *error_approx_fine4 = new doublereal[n_a[4] + 1];
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
						//doublereal *error_approx_fine3 = new doublereal[n_a[3] + 1];
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
					//doublereal *error_approx_fine2 = new doublereal[n_a[2] + 1];
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
				//doublereal *error_approx_fine1 = new doublereal[n_a[1] + 1];
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
			//doublereal *error_approx_fine = new doublereal[n_a[0] + 1];
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
			//seidelqSOR(A, 1, n_a[0], x, b, row_ptr_start, row_ptr_end, 0);
		}

		system("pause");
		
	}

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
	return 0;

} // aggregative_amg


/*
// 6 september 2015 кажется заработало.
integer _tmain(integer argc, _TCHAR* argv[])
{
	// Тестовая обвязка.
	// недостаток данного тестирования отсутствие условий Неймана, условий задания ненулевого теплового потока.
	// Отсутствие сборки на примере матрицы полученной аппроксимацией конвективной задачи с одной из противопоточных схем.
	// 14 сентября 2015 года реальное тестирование на примере расчёта теплового сопротивления транзисторной структуры, 
	// вставляем код в AliceFlow.


    // m=81; 9.3Mb.
	// m=101; 14.3Mb.
	integer m = 21; // 48
	doublereal h = 1.0/ (m - 1);
	doublereal h2 = h*h;

	const integer dim_2D = 1;
	Ak* A=NULL;
	Ak* R=NULL;
	Ak* P=NULL;
	integer n, nnz;
	doublereal* x = NULL;
	doublereal* b = NULL;

	if (dim_2D) {
		// Finite Volume Method matrix assemble. 6 september 2015

		 n = m*m;
		 nnz = 5 * ((m - 2)*(m - 2)) + 4 + 4 * (m - 2);

		// 31 aug 2015
		A = new Ak[6 * nnz + 1];
		R = new Ak[3 * n + 1];
		P = new Ak[3 * n + 1];
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


		x = new doublereal[n+1]; // решение.
		for (integer i = 0; i <= n; i++) x[i] = 0.0;
		b = new doublereal[n+1]; // Правая часть.

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
		 n = m*m*m;
		 nnz = 7 * ((m - 2)*(m - 2)*(m - 2)) + 8 + 6 * (m - 2)*(m-2) + 12*(m-2);
		 // куб : 8 вершин, 12 рёбер, 6 граней.

		// 31 aug 2015
		 // 6 3 3 с запасом. реально 4.52 2.26 2.26
		A = new Ak[(integer)(6 * nnz) + 1]; // 6
		R = new Ak[(integer)(3 * n) + 1]; // 3*nnz 2.4
		P = new Ak[(integer)(3 * n) + 1]; // 3*nnz 2.4
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


		x = new doublereal[n+1]; // решение.
		for (integer i = 0; i <= n; i++) x[i] = 0.0; // инициалзация.
		b = new doublereal[n+1]; // Правая часть.

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

	aggregative_amg(A, nnz, n, R, P, x,b);

	// m in 3D memorysize mb
	// 21  14.8Mb
	// 31  49.7Mb
	// 41  119.4Mb

	//getchar();
	delete[] A;
	delete[] R;
	delete[] P;
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

	return 0;
}
*/


// Здесь содержится обвязка вызывающая amg1r5.
// локальное выдление памяти :всё внутри, многократные alloc и free.
void my_agr_amg_loc_memory(equation3D* &sl, equation3D_bon* &slb,
	integer maxelm, integer maxbound,
	doublereal *dV, doublereal* &dX0,
	doublereal alpharelax, integer iVar, bool bLRfree, QuickMemVorst& m)
{

	// Замер времени.
	unsigned integer calculation_main_start_time; // начало счёта мс.
	unsigned integer calculation_main_end_time; // окончание счёта мс.

	calculation_main_start_time = clock(); // момент начала счёта.



	// На случай если память не была выделена.
	if (dX0 == NULL) {
		dX0 = new doublereal[maxelm + maxbound];
		for (integer i = 0; i<maxelm + maxbound; i++) {
			dX0[i] = 0.0;
		}
	}


	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля
	doublereal res_sum = 0.0;
	res_sum = 0.0;
	for (integer i = 0; i<maxelm; i++) {
		// внутренность матрицы.
		doublereal buf = 0.0;
		buf = (sl[i].ap*dX0[sl[i].iP] - dV[sl[i].iP]);
		if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) buf -= sl[i].ab*dX0[sl[i].iB];
		if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) buf -= sl[i].ae*dX0[sl[i].iE];
		if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) buf -= sl[i].an*dX0[sl[i].iN];
		if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) buf -= sl[i].as*dX0[sl[i].iS];
		if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) buf -= sl[i].at*dX0[sl[i].iT];
		if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) buf -= sl[i].aw*dX0[sl[i].iW];
		buf *= buf;
		res_sum += buf;
	}
	for (integer i = 0; i<maxbound; i++) {
		// граничные узлы.
		doublereal buf = 0.0;
		buf = slb[i].aw*dX0[slb[i].iW] - dV[slb[i].iW];
		if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) buf -= slb[i].ai*dX0[slb[i].iI];
		buf *= buf;
		res_sum += buf;
	}
	res_sum = sqrt(res_sum);
	//printf("residual start=%1.4e\n",res_sum);
	//getchar();

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


	//if (res_sum>1.0E-10) 
	if (res_sum>1.05*finish_residual) // защита от повторного холостого запуска экономит время конечного пользователя.
	{

		//yes_print_amg=false;
		//yes_print_amg = false;



		

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


		integer nna = 0; // количество ненулевых элементов в матрице СЛАУ.


						 // подсчёт числа ненулевых элементов в матрице.
		nna = 0;
		for (integer i = 0; i<maxelm; i++) {
			// внутренность матрицы.
			if ((sl[i].iB>-1) && (fabs(sl[i].ab) > nonzeroEPS)) (nna)++;
			if ((sl[i].iE>-1) && (fabs(sl[i].ae) > nonzeroEPS)) (nna)++;
			if ((sl[i].iN>-1) && (fabs(sl[i].an) > nonzeroEPS)) (nna)++;
			if ((sl[i].iS>-1) && (fabs(sl[i].as) > nonzeroEPS)) (nna)++;
			if ((sl[i].iT>-1) && (fabs(sl[i].at) > nonzeroEPS)) (nna)++;
			if ((sl[i].iW>-1) && (fabs(sl[i].aw) > nonzeroEPS)) (nna)++;
			if ((sl[i].iP>-1) && (fabs(sl[i].ap) > nonzeroEPS)) (nna)++;
		}
		for (integer i = 0; i<maxbound; i++) {
			// граничные узлы.
			if ((slb[i].iW>-1) && (fabs(slb[i].aw) > nonzeroEPS)) (nna)++;
			if ((slb[i].iI>-1) && (fabs(slb[i].ai) > nonzeroEPS)) (nna)++;
		}

		integer nnu = 0; // число неизвестных.
		nnu = maxelm + maxbound;

		/*
		// Рекомендуемые по умолчанию параметры.
		integer nda=0; // память под вектор значений матрицы слау.
		nda=3*(nna)+5*(nnu);
		integer ndia=0;
		ndia=(integer)(2.2*(nnu));
		integer ndja=0;
		ndja=3*(nna)+5*(nnu);
		integer ndu=0;
		ndu=(integer)(2.2*(nnu));
		integer ndf=0;
		ndf=(integer)(2.2*(nnu));
		integer ndig=0;
		ndig=(integer)(5.4*(nnu));
		*/

		/*
		// в двое больше памяти чем рекомендовано.
		integer nda=0; // память под вектор значений матрицы слау.
		nda=6*(nna)+10*(nnu);
		integer ndia=0;
		ndia=(integer)(4.4*(nnu));
		integer ndja=0;
		ndja=6*(nna)+10*(nnu);
		integer ndu=0;
		ndu=(integer)(4.4*(nnu));
		integer ndf=0;
		ndf=(integer)(4.4*(nnu));
		integer ndig=0;
		ndig=(integer)(10.8*(nnu));
		*/

		// данная константа работоспособна вплоть до размерностей сетки равных 34млн 463тысячи 250узлов.
		//doublereal rsize=1.51; // 1048416
		// Вынужденные течения достаточно 2.5. 
		// на задаче Концевого Ю.А. Электростатика со столбиком в случае сетки со сгущением достаточно 2.0.

		
		Ak* A = NULL;
		Ak* R = NULL;
		Ak* P = NULL;
		doublereal* rthdsd_amg = NULL;
		doublereal* result_amg = NULL;

		
		
		

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
		/*          IPRinteger = 10606 */

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
        /*
		integer iswtch = 0;
		iswtch = 4;
		integer iout = 0;
		iout = 12;
		integer iprinteger = 0;
		iprinteger = 10606;
		integer levelx = 0;
		levelx = 25;
		integer ifirst = 0;
		// начальное приближение :
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
		//matrix=11; // symmetric SPD.
		matrix = 22;

		if ((iVar == PAM) && (bLRfree)) {
			//printf("work amg1r5\n");
			//getchar();
			// Симметричная положительно определённая матрица это такая матрица
			// которая возникает для поправки давления при решении вязких несжимаемых уравнений Навье-Стокса в 
			// случае задач : каверна, тест Валь-Девиса. Для задач промышленного масштаба это всякие естественные
			// конвекции охлаждающие висящие в воздухе без контакта с теплоотводом греющиеся изделия.
			// Это особый специфический класс задач.
			matrix = 11;
		}
		*/
		// allocate memory.
		
		A = new Ak[(integer)(6 * nna) + 1]; // 6
		if (A == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for A matrix in my_agregat_amg.cpp algorithm...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
		R = new Ak[(integer)(3 * nnu) + 1]; // 3*nnz 2.4
		if (R == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for R matrix in my_agregat_amg.cpp algorithm...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
		P = new Ak[(integer)(3 * nnu) + 1]; // 3*nnz 2.4
		if (P == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for P matrix in my_agregat_amg.cpp algorithm...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}

		result_amg = new doublereal[nnu + 1];
		if (result_amg == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for result_amg vector in my_agregat_amg.cpp algorithm...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}

		rthdsd_amg = new doublereal[nnu + 1];
		if (result_amg == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for rthdsd_amg vector in my_agregat_amg.cpp algorithm...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}


		// правая часть.
		/*
		for (integer i = 0; i < nnu; i++) {
			rthdsd_amg[i] = 0.0;
			if (i<maxelm + maxbound) {
				// обязательно нужно проверить была ли выделена оперативная память. 
				rthdsd_amg[i + 1] = dV[i];
			}
		}
		rthdsd_amg[0] = 0.0;
		*/
		/*
		// вектор решения. 
		for (integer i = 0; i < nnu; i++) {
			result_amg[i] = 0.0;
			if (i<maxelm + maxbound) {
				// обязательно нужно проверить была ли выделена оперативная память. 
				result_amg[i + 1] = dX0[i];
			}
		}
		result_amg[0] = 0.0;
		*/

		// см. equation3DtoCRS.

		integer ik = 0; // счётчик ненулевых элементов СЛАУ
		integer id = 1;

		// для внутренних узлов расчётной области:
		for (integer k = 0; k<maxelm; k++) {
			//printf("%e %e %e %e %e %e %e\n", sl[k].ap / alpharelax, -sl[k].ae, -sl[k].aw, -sl[k].an, -sl[k].as, -sl[k].at, -sl[k].ab);
			//system("pause");

			if (fabs(sl[k].ap) > nonzeroEPS) {
				A[ik + id].aij = sl[k].ap / alpharelax;
				A[ik + id].j = sl[k].iP + 1;
				A[ik + id].i = sl[k].iP + 1;
				rthdsd_amg[sl[k].iP + 1] = dV[sl[k].iP];
				result_amg[sl[k].iP + 1] = dX0[sl[k].iP];
				ik++;
			}
			if ((sl[k].iE>-1) && (fabs(sl[k].ae) > nonzeroEPS)) {
				A[ik + id].aij = -sl[k].ae;
				A[ik + id].j = sl[k].iE + 1;
				A[ik + id].i = sl[k].iP + 1;
				ik++;
			}
			if ((sl[k].iN>-1) && (fabs(sl[k].an) > nonzeroEPS)) {
				A[ik + id].aij = -sl[k].an;
				A[ik + id].j = sl[k].iN + 1;
				A[ik + id].i = sl[k].iP + 1;
				ik++;
			}
			if ((sl[k].iT>-1) && (fabs(sl[k].at) > nonzeroEPS)) {
				A[ik + id].aij = -sl[k].at;
				A[ik + id].j = sl[k].iT + 1;
				A[ik + id].i = sl[k].iP + 1;
				ik++;
			}
			if ((sl[k].iS>-1) && (fabs(sl[k].as) > nonzeroEPS)) {
				A[ik + id].aij = -sl[k].as;
				A[ik + id].j = sl[k].iS + 1;
				A[ik + id].i = sl[k].iP + 1;
				ik++;
			}
			if ((sl[k].iW>-1) && (fabs(sl[k].aw) > nonzeroEPS)) {
				A[ik + id].aij = -sl[k].aw;
				A[ik + id].j = sl[k].iW + 1;
				A[ik + id].i = sl[k].iP + 1;
				ik++;
			}
			if ((sl[k].iB>-1) && (fabs(sl[k].ab) > nonzeroEPS)) {
				A[ik + id].aij = -sl[k].ab;
				A[ik + id].j = sl[k].iB + 1;
				A[ik + id].i = sl[k].iP + 1;
				ik++;
			}


		}


		// для внутренних узлов расчётной области:
		for (integer k = 0; k<maxbound; k++) {
			if (fabs(slb[k].aw) > nonzeroEPS) {
				// val[ik]=slb[k].aw/alpharelax;
				A[ik + id].aij = slb[k].aw; // релаксация для граничных узлов не применяется.
										/*if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
										// Внимание !!! было произведено тестирование : один вариант был с нижней релаксацией для граничных узлов,
										// а второй вариант был без нижней релаксации на граничных узлах. Было выяснено, что для сходимости
										// более благоприятен вариант без нижней релаксации на граничных узлах.
										// Данное изменение согласовано с функцией solve.

										val[ik]/=alpharelax; // Если условия Неймана то нижняя релаксация.
										}*/
				A[ik + id].j = slb[k].iW + 1;
				A[ik + id].i = slb[k].iW + 1; // dirichlet
				rthdsd_amg[slb[k].iW + 1] = dV[slb[k].iW];
				result_amg[slb[k].iW + 1] = dX0[slb[k].iW];
				ik++;
			}
			if ((slb[k].iI>-1) && (fabs(slb[k].ai) > nonzeroEPS)) {
				A[ik + id].aij = -slb[k].ai;
				A[ik + id].j = slb[k].iI + 1;
				A[ik + id].i = slb[k].iW + 1;
				rthdsd_amg[slb[k].iW + 1] = dV[slb[k].iW];
				result_amg[slb[k].iW + 1] = dX0[slb[k].iW];
				// Это очень важный вопрос и он требует проверки !

				ik++;
			}

		}


		// TODO : 
		// нужно акуратно прописать выделения и уничтожения памяти с учётом того что было сделано в BiCGStabP.

		



		//printf("getready ...");
		//getchar();
		// amg - особенно хорош для поправки давления в SIMPLE алгоритме.
		// свой алгоритм 2015 года.
	

		aggregative_amg(A, nna, nnu, R, P, result_amg, rthdsd_amg);

		system("pause");
		
		// возвращаем решение СЛАУ.
		// вектор решения. 
		for (integer i = 0; i < nnu; i++) {
			if (i<maxelm + maxbound) {
				// обязательно нужно проверить была ли выделена оперативная память. 
				dX0[i]=result_amg[i + 1];
			}
		}
		


		// освобождение памяти.
		if (A != NULL) {
			delete[] A;
		}
		if (R != NULL) {
			delete[] R;
		}
		if (P != NULL) {
			delete[] P;
		}
		if (result_amg != NULL) {
			delete[] result_amg;
		}
		if (rthdsd_amg != NULL) {
			delete[] rthdsd_amg;
		}
		


		res_sum = 0.0;
		for (integer i1 = 0; i1<maxelm; i1++) {
			// внутренность матрицы.
			doublereal buf = 0.0;
			buf = (sl[i1].ap*dX0[sl[i1].iP] - dV[sl[i1].iP]);
			if ((sl[i1].iB>-1) && (fabs(sl[i1].ab) > nonzeroEPS)) buf -= sl[i1].ab*dX0[sl[i1].iB];
			if ((sl[i1].iE>-1) && (fabs(sl[i1].ae) > nonzeroEPS)) buf -= sl[i1].ae*dX0[sl[i1].iE];
			if ((sl[i1].iN>-1) && (fabs(sl[i1].an) > nonzeroEPS)) buf -= sl[i1].an*dX0[sl[i1].iN];
			if ((sl[i1].iS>-1) && (fabs(sl[i1].as) > nonzeroEPS)) buf -= sl[i1].as*dX0[sl[i1].iS];
			if ((sl[i1].iT>-1) && (fabs(sl[i1].at) > nonzeroEPS)) buf -= sl[i1].at*dX0[sl[i1].iT];
			if ((sl[i1].iW>-1) && (fabs(sl[i1].aw) > nonzeroEPS)) buf -= sl[i1].aw*dX0[sl[i1].iW];
			buf *= buf;
			res_sum += buf;
		}
		for (integer i1 = 0; i1<maxbound; i1++) {
			// граничные узлы.
			doublereal buf = 0.0;
			buf = slb[i1].aw*dX0[slb[i1].iW] - dV[slb[i1].iW];
			if ((slb[i1].iI>-1) && (fabs(slb[i1].ai) > nonzeroEPS)) buf -= slb[i1].ai*dX0[slb[i1].iI];
			buf *= buf;
			res_sum += buf;
		}
		res_sum = sqrt(res_sum);
		//printf("residual finish=%1.4e\n",res_sum);
		//getchar();
		if (bsolid_static_only) {
			// используется только для теплопередачи в твёрдом теле для ускорения
			// решения задачи - защита от рестарта.
			finish_residual = res_sum; // значение невязки решённой задачи.
		}

	}

	calculation_main_end_time = clock();
	calculation_vorst_seach_time += calculation_main_end_time - calculation_main_start_time;

} // my_agr_amg_loc_memory
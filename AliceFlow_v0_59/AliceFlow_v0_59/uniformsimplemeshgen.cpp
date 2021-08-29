// Файл uniformsimplemeshgen.cpp
// содержит генератор равномерной 3D расчётной сетки.
// 19 декабря 2011 добавлена возможность использовать неравномерную сетку
// в соответствии с законом геометрической прогрессии (сгущение к двум границам).

#pragma once
#ifndef _UNIFORMSIMPLEMESHGEN_CPP_
#define _UNIFORMSIMPLEMESHGEN_CPP_ 1


void my_solid_properties(doublereal TiP, float &rho, float &cp, float &lam, integer ilibident);

// минимальное количество клеток по 
// каждому из направлений на минимальном
// элементе разбиения
int min_elem_in_x_element=4; // обязательно не меньше 4!
int min_elem_in_y_element=4; // обязательно не меньше 4!
int min_elem_in_z_element=4; // обязательно не меньше 4!
// ручная адаптация сетки
int adapt_x=0;
int adapt_y=0;
int adapt_z=0;
// разбиение диаметра цилиндра.
const doublereal dcount_diametr_cylinder = 8.0;
const bool bcylinder_meshing = true;

/*const*/ /*doublereal etalon_max_size_ratio=2.0;*/ // см inputlaplas.cpp



// реализация функции SetLength из Delphi 
// для изменения размеров динамического массива ra. 
// Возвращает модифицированный массив ra необходимого 
//  размера. 
// Во время работы функции утечки памяти не происходит.
void SetLength(doublereal* &ra, integer isizeold, integer isize) 
{

	// isize - новая длина динамического массива.
	doublereal *temp=nullptr;
	temp = new doublereal[isize];
    integer i;
	for (i=0; i<isize; i++) temp[i]=0.0; // инициализация
    /*
	integer isizeold;
	if (ra != nullptr) {
       isizeold = sizeof(ra)/sizeof(ra[0]); // длина старого массива
	   #if doubleintprecision == 1
			printf("%lld  ",isizeold); // не хочет правильно определять размер массива
	   #else
			printf("%d  ",isizeold); // не хочет правильно определять размер массива
	   #endif

       
	} 
	  else
	{
       isizeold=0;
	}
	*/
	if (isize < isizeold) isizeold=isize;
	for (i=0; i<isizeold; i++) temp[i]=ra[i];
	
	//if (ra != nullptr) {
		delete[]  ra; // уничтожение объекта
		//ra = nullptr;
	//}
	ra = new doublereal[isize]; // выделение памяти
	for (i=0; i<isize; i++) ra[i]=temp[i]; // копирование
	delete[] temp; // освобождение памяти
	temp = nullptr;
	
} // SetLength

// добавляет несуществующую границу к массиву
// Теперь с учётом координатного направления 13.08.2019
//template <class ITYPE>
void addboundary(doublereal* &rb, int &in, doublereal g, integer iDir,
	BLOCK* &b, int &lb, WALL* &w, int &lw, SOURCE* &s, int &ls) {
	// rb - модифицируемый массив границ,
	// in - номер последней границы в массиве, нумерация начинается с нуля.
	// g - граница, кандидат на добавление.
	// iDir - координатное направление XY(Z directional), XZ(Y directional), YZ(X directional).
	
	// Зачем расстояния менее одной десятой мкм. 
	// Гипотеза в том что это не физично для обычных моделей.
	// const doublereal eps = 0.1e-6;//1e-30; // для отделения вещественного нуля.
	// Это наиважнейший параметр влияющий на быстродействие и точность. 
	// Разумно подобрав его быстродействие возрастает более чем на порядок,
	// а те задачи которые раньше невозможно было посчитать начинают сходиться.
	// 13.08.2019
	doublereal eps = 1.0e-10;// shorter_length_for_simplificationX(g);
	if (b_adhesion_Mesh) {
		// Долгое время успешно работало для значения eps=1.0e-10 для всех направлений.
		switch (iDir) {
		case XY_PLANE: eps = shorter_length_for_simplificationZ(g, b, lb, w, lw, s, ls);
			break;
		case XZ_PLANE: eps = shorter_length_for_simplificationY(g, b, lb, w, lw, s, ls);
			break;
		case YZ_PLANE: eps = shorter_length_for_simplificationX(g, b, lb, w, lw, s, ls);
			break;
		default:
			printf("fatal error!!! unknown directional on function addboundary(...) in module uniformsimplemeshgen.cpp\n");
			system("pause");
			exit(1);
			break;
		}
	}
	else {
		eps= 1.0e-10;
	}

	bool bfind=false;
	for (integer i=0; i<=in; i++) if (fabs(rb[i]-g) < eps/*admission*/) bfind=true;
	if (!bfind) {
        SetLength(rb, in+1, in+2);
		in++;
		rb[in]=g; // запись добавляемой границы в конец динамического массива.
	}
} // addboundary



// добавляет несуществующую границу к массиву
// Теперь с учётом координатного направления 13.08.2019
void addboundary_rudiment(doublereal*& rb, int& in, doublereal g, integer iDir,
	BLOCK*& b, int& lb, WALL*& w, int& lw, SOURCE*& s, int& ls) {
	// rb - модифицируемый массив границ,
	// in - номер последней границы в массиве, нумерация начинается с нуля.
	// g - граница, кандидат на добавление.
	// iDir - координатное направление XY(Z directional), XZ(Y directional), YZ(X directional).

	// Зачем расстояния менее одной десятой мкм. 
	// Гипотеза в том что это не физично для обычных моделей.
	// const doublereal eps = 0.1e-6;//1e-30; // для отделения вещественного нуля.
	// Это наиважнейший параметр влияющий на быстродействие и точность. 
	// Разумно подобрав его быстродействие возрастает более чем на порядок,
	// а те задачи которые раньше невозможно было посчитать начинают сходиться.
	// 13.08.2019
	doublereal eps = 1.0e-12;// shorter_length_for_simplificationX(g);
	// Долгое время успешно работало для значения eps=1.0e-10 для всех направлений.
	/*switch (iDir) {
	case XY: eps = shorter_length_for_simplificationZ(g, b, lb, w, lw, s, ls);
		break;
	case XZ: eps = shorter_length_for_simplificationY(g, b, lb, w, lw, s, ls);
		break;
	case YZ: eps = shorter_length_for_simplificationX(g, b, lb, w, lw, s, ls);
		break;
	default:
		printf("fatal error!!! unknown directional on function addboundary(...) in module uniformsimplemeshgen.cpp\n");
		system("pause");
		exit(1);
		break;
	}
	*/
	// Отличие в том что данный метод добавляет почти всегда, а не добавляет лишь если расстояние между линиями менее 1.0Е-10.
	bool bfind = false;
	for (int i = 0; i <= in; i++) if (fabs(rb[i] - g) < eps/*admission*/) bfind = true;
	if (!bfind) {
		SetLength(rb, in + 1, in + 2);
		in++;
		rb[in] = g; // запись добавляемой границы в конец динамического массива.
	}
	
} // addboundary_rudiment


/*
  // добавляет несуществующую границу к массиву
  // функция устарела 13.08.2019.
void addboundary(doublereal* &rb, integer &in, doublereal g) {
	// rb - модифицируемый массив границ,
	// in - номер последней границы в массиве, нумерация начинается с нуля.
	// g - граница, кандидат на добавление.

	// Зачем расстояния менее одной десятой мкм. 
	// Гипотеза в том что это не физично для обычных моделей.
	//const doublereal eps = 0.1e-6;//1e-30; // для отделения вещественного нуля.
	// Это наиважнейший параметр влияющий на быстродействие и точность. 
	// Разумно подобрав его быстродействие возрастает более чем на порядок,
	// а те задачи которые раньше невозможно было посчитать начинают сходиться.
	// 13.08.2019
	const doublereal eps = shorter_length_for_simplification;
	// Долгое время успешно работало для значения eps=1.0e-10.
	bool bfind = false;
	for (integer i = 0; i <= in; i++) if (fabs(rb[i] - g)<eps) bfind = true;
	if (!bfind) {
		SetLength(rb, in + 1, in + 2);
		in++;
		rb[in] = g; // запись добавляемой границы в конец динамического массива.
	}
} // addboundary
*/

//С помощью неполной НЕУБЫВАЮЩЕЙ кучи 
//крупные элементы закидываем поближе к концу массива 
template <typename myARRT>
void reheap(myARRT a[], int length, int i)  {
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	myARRT Temp = a[i];
	int parent = i;
	int child = 2 * (i + 1) - 1;
	//Просматриваем потомков, а также потомков потомков 
	//и сравниваем их с родителем (если что - передвигаем потомков влево) //Цикл продолжается пока не выпадем за пределы массива
	 //или пока не обменяем какого-нибудь потомка на родителя. 
	while ((child < length) && (!done)) {
		//Если правый потомок в пределах массива 
		if (child < length - 1) {
			//То из левого и правого потомка выбираем наименьшего 
			if (a[child] >= a[child + 1]) { child += 1; }
		}
		//Родитель меньше потомков? 
		if (Temp < a[child]) {
			//Тогда с этим родителем и его потомками разобрались
			done = true;
			//Родитель НЕ меньше чем наименьший из его потомков. 
			//Перемещаем потомка на место родителя 
			//и с родителем в цикле сравниваем уже потомков этого потомка 
		}
		else
		{
			a[parent] = a[child];
			parent = child;
			child = 2 * (parent + 1) - 1;
		}
	}
	//Родитель, с которого всё начиналось 
	//передвигается ближе к концу массива 
	//(или остаётся на месте если не повезло) 
	a[parent] = Temp;
}

//С помощью неполной НЕВОЗРАСТАЮЩЕЙ кучи 
//мелкие элементы закидываем поближе к началу массива 
template <typename myARRT>
void invreheap(myARRT a[], int length, int i)
{
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	myARRT Temp = a[length - 1 - i];
	int parent = i;
	int child = 2 * (i + 1) - 1;
	//Просматриваем потомков, а также потомков потомков 
	//и сравниваем их с родителем (если что - передвигаем потомков) 
	//Цикл продолжается пока не выпадем за пределы массива 
	//или пока не обменяем какого-нибудь потомка на родителя. 
	while ((child < length) && (!done))
	{
		//Если левый потомок в пределах массива 
		if (child < length - 1)
		{
			//То из левого и правого потомка выбираем наибольшего 
			if (a[length - 1 - child] <= a[length - 1 - (child + 1)])
			{
				child += 1;
			}
		}
		//Родитель больше потомков? 
		if (Temp > a[length - 1 - child])
		{
			//Тогда с этим родителем и его потомками разобрались 
			done = true;
		}
		else
		{
			//Родитель НЕ больше чем наибольший из его потомков. 
			//Перемещаем потомка на место родителя 
			//и с родителем в цикле сравниваем уже потомков этого потомка 
			a[length - 1 - parent] = a[length - 1 - child];
			parent = child;
			child = 2 * (parent + 1) - 1;
		}
	}
	//Родитель, с которого всё начиналось 
	//передвигается ближе к началу массива 
	//(или остаётся на месте если не повезло) 
	a[length - 1 - parent] = Temp;
}


/** * Демонстраицонный алгоритм для J-сортировки (JSort).
* Автор алгоритма - Джейсон Моррисон (Jason Morrison)
* <http://www.scs.carleton.ca/~morrison>
* * JSortAlgorithm.java
* * Автор реализации - Патрик Морин 
* @author Patrick Morin */

//Основная процедура сортировки 
template <typename myARRT>
void sortJ(myARRT a[], integer in)
{
	//Строим неубывающую кучу 
	//Большие элементы из начала массива 
	//закидываем поближе к концу 
	for (integer i = in; i >= 0; i--) reheap(a, in+1, i);
	//Строим невозрастающую кучу
	//Меньшие элементы из конца массива 
	//закидываем поближе к началу 
	for (integer i = in; i >= 0; i--) invreheap(a, in+1, i);
	//Массив ПОЧТИ упорядочен 
	//Досортировываем вставками 
	for (integer j = 1; j <= in; j++)
	{
		myARRT Temp = a[j];
		integer i = j - 1;
		while (i >= 0 && a[i] > Temp)
		{
			a[i + 1] = a[i];
			i -= 1;
		}
		a[i + 1] = Temp;
	}
} // sortJ


// Сортировка по возрастанию прямым обменом -
// пузырьковая улучшенная.
// Основной принцип: если элементы уже упорядочены, 
// сортировка прекращается.
// in - предполагается достаточно малым меньше 500,
// именно поэтому применима пузырьковая сортировка.
// rb сортируемый массив от нуля до in включительно.
// rb [left..right]
template <typename myARRT>
void BubbleEnhSort(myARRT* &rb, integer first, integer last) {
	integer i,j,k;
	myARRT x;
	bool swapped; // факт обмена

	for (i= first+1; i<= last; i++) {
		k=0; // количество обменов данного прохода
		swapped=false; // обменов не было
		// i пузырьков уже всплыло, их просматривать ненужно
		for (j= last; j>=i; j--) {
			if (rb[j-1]>rb[j]) {
				// SWAP
				x=rb[j-1];
				rb[j-1]=rb[j];
				rb[j]=x;
				k++;
				swapped=true;
			}
		}
		if (!swapped) break; // выход из цикла
	}
} // BubbleEnhSort

// Сортировка Дональда Шелла из Седжвика 1959год.
// rb[0..in]
template <typename myARRT>
void ShellSort(myARRT* &rb, integer in) {
	integer i, j;
	myARRT x;
	integer h;

	//for (h = 1; h <= in / 9; h = 3 * h + 1);
	h = 1;
	while (h <= in / 9) {
		h = 3 * h + 1;
	}
	for (; h > 0; h /= 3) {
		for (i = h; i <= in; i++) {
			j = i; 
			x = rb[i];

			while (j >= h && x < rb[j-h]) {
				rb[j] = rb[j - h]; j -= h;
			}
			rb[j] = x;
		}
	}
} // ShellSort 1959


// Сортировка Тима (Тим Петерсон).
const integer RUN = 32;
const bool INS_Sort = true;

// this function sorts array from left index to
// to right index which is of size atmost RUN
template <typename myARRT>
void insertionSortTim(myARRT arr[], integer left, integer right)
{
	for (integer i = left + 1; i <= right; i++)
	{
		myARRT temp = arr[i];
		integer j = i - 1;
		while (arr[j] > temp && j >= left)
		{
			arr[j + 1] = arr[j];
			j--;
		}
		arr[j + 1] = temp;
	}
}

// merge function merges the sorted runs
template <typename myARRT>
void mergeTim(myARRT arr[], integer l, integer m, integer r)
{
	// original array is broken in two parts
	// left and right array
	integer len1 = m - l + 1, len2 = r - m;

	//myARRT left[len1], right[len2];
	myARRT *left = nullptr; 
	if (len1 >= 0) {
		left = new myARRT[len1 + 1];
	}
	myARRT *right = nullptr;
	if (len2 >= 0) {
		right = new myARRT[len2 + 1];
	}

	for (integer i = 0; i < len1; i++)
		left[i] = arr[l + i];
	for (integer i = 0; i < len2; i++)
		right[i] = arr[m + 1 + i];

	integer i = 0;
	integer j = 0;
	integer k = l;

	// after comparing, we merge those two array
	// in larger sub array
	while (i < len1 && j < len2)
	{
		if (left[i] <= right[j])
		{
			arr[k] = left[i];
			i++;
		}
		else
		{
			arr[k] = right[j];
			j++;
		}
		k++;
	}

	// copy remaining elements of left, if any
	while (i < len1)
	{
		arr[k] = left[i];
		k++;
		i++;
	}

	// copy remaining element of right, if any
	while (j < len2)
	{
		arr[k] = right[j];
		k++;
		j++;
	}

	delete[] left;
	delete[] right;
}

// Возвращает минимум из двух целых чисел.
// 8.05.2018
//integer min(integer a, integer b) {
	//if (a <= b) return a;
	//else return b;
//}


// Требует дополнительной памяти O(n).
// iterative Timsort function to sort the
// array[0...n-1] (similar to merge sort)
template <typename myARRT>
void timSort(myARRT arr[], integer n)
{
	// Sort individual subarrays of size RUN
	for (integer i = 0; i < n; i += RUN)
		insertionSortTim(arr, i, ((i + RUN - 1) < (n - 1) ? (i + RUN - 1) : (n - 1)));

	// start merging from size RUN (or 32). It will merge
	// to form size 64, then 128, 256 and so on ....
	for (integer size = RUN; size < n; size = 2 * size)
	{
		// pick starting point of left sub array. We
		// are going to merge arr[left..left+size-1]
		// and arr[left+size, left+2*size-1]
		// After every merge, we increase left by 2*size
		for (integer left = 0; left < n; left += 2 * size)
		{
			// find ending point of left sub array
			// mid+1 is starting point of right sub array
			integer mid = left + size - 1;
			integer right = ((left + 2 * size - 1) < (n - 1) ? (left + 2 * size - 1) : (n - 1));

			// merge sub array arr[left.....mid] &
			// arr[mid+1....right]
			if ((mid < right) && (mid >= left)) {
				mergeTim(arr, left, mid, right);
			}
		}
	}
}

// utility function to print the Array
template <typename myARRT>
void printArray(myARRT arr[], integer n)
{
	for (integer i = 0; i < n; i++)
		printf("%lld  ", arr[i]);
	printf("\n");
}

// Driver program to test above function
integer test_Tim_Sort()
{
	integer arr[] = { 5, 21, 7, 23, 19 };
	integer n = sizeof(arr) / sizeof(arr[0]);
	printf("Given Array is\n");
	printArray<integer>(arr, n);

	timSort<integer>(arr, n);

	printf("After Sorting Array is\n");
	printArray(arr, n);
	return 0;
}


  // Для генерации матрицы СЛАУ требуется в случае реализации
  // на динамических массивах переупорядочивание элементов:
  // сортировка. Здесь будет реализована быстрая сортировка.
  // Брайан Керниган и Денис Ритчи "The C programming language".
  // swap: Обмен местами v[i] и v[j]
template <typename myARRT>
void swap(myARRT* &v, integer i, integer j)
{
	myARRT temp;

	// change v[i] <-> v[j]
	temp = v[i];
	v[i] = v[j];
	v[j] = temp;

} // swap

  // Вот алгоритм PivotList
template <typename myARRT>
integer PivotList(myARRT* &list, integer first, integer last) {
	// list обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	myARRT PivotValue = list[first];
	integer PivotPoint = first;

	for (integer index = (first + 1); index <= last; index++) {
		if (list[index]<PivotValue) {
			PivotPoint++;
			swap<myARRT>(list, PivotPoint, index);
		}
	}

	swap<myARRT>(list, first, PivotPoint);

	return PivotPoint;
} // PivotList


  // Быстрая сортировка Хоара.
  // Запрограммировано с использованием ДЖ. Макконелл Анализ алгоритмов
  // стр. 106.
template <typename myARRT>
void QuickSort(myARRT* &list, integer first, integer last) {
	// list упорядочиваемый список элементов
	// first номер первого элемента в сортируемой части списка
	// last номер последнего элемента в сортируемой части списка

	integer pivot;

	if (last - first > 15) {
		if (first < last) {
			pivot = PivotList<myARRT>(list, first, last);
			QuickSort<myARRT>(list, first, pivot - 1);
			QuickSort<myARRT>(list, pivot + 1, last);
		}
	}
	else {
		// менее 16 элементов
		BubbleEnhSort<myARRT>(list, first, last);
	}
} // QuickSort

// Обётка для сортировки (выбор метода сортировки).
// не использует дополнительной оперативной памяти.
template <typename myARRT>
void Sort_method(myARRT* &rb, integer in) {
	//BubbleEnhSort<myARRT>(rb, 0, in);
	// 24 03 2017
	ShellSort<myARRT>(rb, in);
	// 15.03.2019 
	// Тим Петерсон. 
	//timSort<myARRT>(rb, in+1);
	// 15.03.2019
	//QuickSort<myARRT>(rb, 0, in);
	// 15.03.2019
	//sortJ<myARRT>(rb, in);
}

// иногда сетку вблизи источников тепла
// желательно сделать более подробной.
// Данная адаптация подходит для следующего расположения 
// источников тепла: в группе два источника с расстоянием между центрами dx1,
// Объединение состоит из нескольких групп с расстоянием между группами dx2.
// Крайние источники объединения граничат с грубой сеткой, для этого 
// с краёв объединения добавляется три клетки величиной dx2/n, где n - количество
// клеток в зазоре между группами.
// Плюс к этому для каждого источника каждая внешняя клетка разбивается пополам.
// Ограничение для такой адаптации: источники должны быть расположены только вдоль 
// одного направления и только в плоскости XY.
void simplecorrect_meshgen_x(doublereal* &xpos,  int &inx, 
				   int lb, int ls, int lw, BLOCK* b, SOURCE* s, WALL* w) {
	// вызывается сразу после simplemeshgen()
	int inxcopy=inx;
	doublereal *xposcopy = new doublereal[inx+1]; // копия xpos
	int i=0;
	for (i=0; i<=inx; i++) xposcopy[i]=xpos[i]; // запоминание массива координат

    int ipol=0;
	integer j,i1,i2; // счётчики цикла for

    // источники
	for (i=0; i<ls; i++) {
		if (s[i].iPlane==XY_PLANE) {
		   // в плоскости XY
           ipol++;
		}
	}

	integer *ileft=new integer[ipol+1]; ileft[ipol]=-1; // заглушка
	integer *iright=new integer[ipol+1]; iright[ipol]=-1; // заглушка

	i1=0; i2=0;
    // источники
	for (i=0; i<ls; i++) {
		if (s[i].iPlane==XY_PLANE) {
		   // в плоскости XY
           for (j=0; j<=inx; j++) {
			   // левые границы источников
			   if (fabs(xpos[j]-s[i].g.xS)<admission) { 
				   ileft[i1]=j;
				   i1++;
			   }
               // правые границы источников
			   if (fabs(xpos[j]-s[i].g.xE)<admission) { 
				   iright[i2]=j;
				   i2++;
			   }
		   }
		}
	}

	// массивы ileft, iright должны быть отсортированы
    //BubbleEnhSort<integer>(ileft, 0, ipol-1);
    //BubbleEnhSort<integer>(iright, 0, ipol-1);
	QuickSort<integer>(ileft, 0, ipol - 1);
	QuickSort<integer>(iright, 0, ipol - 1);

	// distmin - размер минимальной клетки
	// в зазоре между группами источников. 
	// В группе 2 источника с расстоянием между ними как правило 26 мкм. 
	doublereal distmin=1e20;
    integer k=0;
    for (j=0; j<=inxcopy; j++) {
		if (j==ileft[k]) {
			// найдена левая граница
			if ((ileft[k]!=0)&&(distmin>(xposcopy[j]-xposcopy[j-1]))) distmin=(xposcopy[j]-xposcopy[j-1]);
		}
        if (j==iright[k]) {
			// найдена левая граница
            if ((iright[k]!=inx)&&(distmin>(xposcopy[j+1]-xposcopy[j]))) distmin=(xposcopy[j+1]-xposcopy[j]);
			k++;
		}
	}

    // вставка дополнительного разбиения так чтобы 
	// крайние источники оказались в равноправном положении
	// по мелкости сетки как и центральные источники.
	bool bflag=false;
    k=0;
    for (j=0; j<=inxcopy; j++) {
		if ((j==ileft[k])&&(ileft[k]!=0)) {
			// найдена левая граница
			if (3.0*distmin<(xposcopy[j]-xposcopy[j-1])) {
               // Вставка трёх интервалов
               SetLength(xpos,inx+1,inx+4);
			   xpos[inx+1]=xposcopy[j]-3.0*distmin;
               xpos[inx+2]=xposcopy[j]-2.0*distmin;
               xpos[inx+3]=xposcopy[j]-distmin;
			   inx=inx+3;
               bflag=true;
			}
		}
		// 2 november 2016.
		//if ((j == iright[k]) && (iright[k]!=inxcopy)) {
        if ((j==iright[k])&&(iright[k]<inxcopy)) {
			// найдена левая граница
			if (3.0*distmin<(xposcopy[j+1]-xposcopy[j])) {
				// Вставка трёх интервалов
                SetLength(xpos,inx+1,inx+4);
				xpos[inx+3]=xposcopy[j]+distmin;
                xpos[inx+2]=xposcopy[j]+2.0*distmin;
                xpos[inx+1]=xposcopy[j]+3.0*distmin;
			    inx=inx+3;
                bflag=true;
			}
			k++;
		}
	}

	if (bflag) {
		// сортировка по возрастанию
        //BubbleEnhSort<doublereal>(xpos, 0, inx);
		Sort_method<doublereal>(xpos,inx);
	    SetLength(xposcopy,inxcopy+1,inx+1);
	    inxcopy=inx;
	}
	
    for (i=0; i<=inx; i++) xposcopy[i]=xpos[i]; // запоминание массива координат
	inx=inx+2*ipol;
	delete[] xpos;
	xpos = new doublereal[inx+1];
	i=0; k=0; 
	for (j=0; j<=inxcopy; j++) {
		
			if ((j == ileft[k]) && (ileft[k] != 0)) {
				// найдена левая граница
				if (i <= inx) {
					xpos[i++] = 0.5*(xposcopy[j - 1] + xposcopy[j]);
				}
			}
			if (i <= inx) {
				xpos[i++] = xposcopy[j];
			}
			if ((j == iright[k]) && (iright[k] != inxcopy)) {
				// найдена левая граница
				if (i <= inx) {
					xpos[i++] = 0.5*(xposcopy[j] + xposcopy[j + 1]);
				}
				k++;
			}
		
	}

} // simplecorrect_meshgen_x

/*
doublereal fmin(doublereal ra, doublereal rb) {
	doublereal r=ra;
	if (rb<ra) r=rb;
	return rb;
} // fmin
*/

// Источники расположены только в плоскости XY и только на одном уровне по Z.
// можно запускать подряд дважды и более.
void simplecorrect_meshgen_z(doublereal* &zpos,  int &inz, 
				   int lb, int ls, int lw, BLOCK* b, SOURCE* s, WALL* w) {

	// вызывается сразу после simplemeshgen()
	int i=0;

	bool bfindsourseXY=false;
	// источники
	for (i=0; i<ls; i++) {
		if (s[i].iPlane==XY_PLANE) {
			bfindsourseXY=true;
		}
	}

	if (bfindsourseXY) {
	
		// Только если обнаружены источники тепла лежащие в плоскости XOY.

	    doublereal rzlevel;
        // источники
	    for (i=0; i<ls; i++) {
		    if (s[i].iPlane==XY_PLANE) {
		       // в плоскости XY
			   rzlevel=s[i].g.zS;
			   break; // выход из цикла for
		    }
	    }

	    integer izlevel=0;
	    for (i=0; i<=inz; i++) if (fabs(zpos[i]-rzlevel)<admission) { 
		    izlevel=i;
	    	break; // досрочный выход из цикла for
	    }

#if doubleintprecision == 1
		//printf("izlevel==%lld inz==%lld\n",izlevel, inz);
#else
		//printf("izlevel==%d inz==%d\n",izlevel, inz);
#endif
       
        //getchar(); // debug

	    doublereal dzminus, dzplus;
	    dzminus=zpos[izlevel]-zpos[izlevel-1];
        if (izlevel==inz) {
            // Вверх увеличивать нельзя там нету расчётной области.

            // увеличиваем подробность снизу на пол клетки.
            SetLength(zpos,inz+1, inz+2);
            zpos[inz+1]=zpos[inz];
            zpos[inz]=zpos[inz-1]+0.5*dzminus;
            inz++;
	    }
	    else {
            // Мельчим сетку сверху над источником в воздухе.
	        dzplus=zpos[izlevel+1]-zpos[izlevel];
	        if (dzplus>3*dzminus) {
		       SetLength(zpos,inz+1, inz+4);
		       zpos[inz+1]=rzlevel+dzminus;
		       zpos[inz+2]=rzlevel+2.0*dzminus;
		       zpos[inz+3]=rzlevel+3.0*dzminus;
		       inz=inz+3;
               //BubbleEnhSort<doublereal>(zpos,0, inz);
			   Sort_method<doublereal>(zpos, inz);
	        }
	    }
	    /*else if (dzminus>3*dzplus) {
          SetLength(zpos,inz+1, inz+4);
		  zpos[inz+1]=rzlevel-dzplus;
		  zpos[inz+2]=rzlevel-2.0*dzplus;
		  zpos[inz+3]=rzlevel-3.0*dzplus;
		  inz=inz+3;
          //BubbleEnhSort<doublereal>(zpos,0, inz);
		  Sort_method<doublereal>(zpos, inz);
	    }*/

	    for (i=0; i<=inz; i++) if (fabs(zpos[i]-rzlevel)<admission) { 
	    	izlevel=i;
		    break; // досрочный выход из цикла for
	    }

        if (izlevel!=inz) {
           // Если сверху после источника есть расчётная область !
	       dzminus=fmin(0.5*(zpos[izlevel]-zpos[izlevel-1]),0.5*(zpos[izlevel+1]-zpos[izlevel]));
           SetLength(zpos,inz+1, inz+2);
	       zpos[inz+1]=rzlevel+dzminus;
	       //zpos[inz+2]=rzlevel-dzminus;
           inz=inz+1;
	    }
        //BubbleEnhSort<doublereal>(zpos, 0, inz);
		Sort_method<doublereal>(zpos, inz);


		// В плоскости XOY были обнаружены источники тепла.
	}

} // simplecorrect_meshgen_z

// Предпологается что плоские источники тепла расположеннные
// группами вдоль выделенного направления расположены в плоскости XY.
// запускать можно повторно и более раз.
void simplecorrect_meshgen_y(doublereal* &ypos,  int &iny, 
				   int lb, int ls, int lw, BLOCK* b, SOURCE* s, WALL* w) {

    // Адаптация  сетки только для первого источника и только с боков по оси OY 
    // источник при этом лежит в плоскости OXY и его нормаль по оси Oz.
    // Причём источник должен быть строго внутренний иначе ошибка при выходе за границы расчётной области.

   // вызывается сразу после simplemeshgen()
	int i=0;

	doublereal rymin=1.0e30, rymax=-1.0e30;
    // источники
	// Правил 14.08.2019
	for (i=0; i<ls; i++) {
		if (s[i].iPlane==XY_PLANE) {
		   // в плоскости XY
			if (s[i].g.yS<rymin) rymin=s[i].g.yS;
			if (s[i].g.yE>rymax) rymax=s[i].g.yE;
			break; // выход из цикла for
		}
	}

    doublereal dyminus, dyplus;

	integer iup=0,idouwn=0;
	for (i=0; i<=iny; i++) if (fabs(ypos[i]-rymin)<admission) { 
		idouwn=i;
		break; // досрочный выход из цикла for
	}
    
    if (idouwn!=0) {
	   
	   dyminus=ypos[idouwn]-ypos[idouwn-1];
	   dyplus=ypos[idouwn+1]-ypos[idouwn];
       if (dyminus>3*dyplus) {
          SetLength(ypos,iny+1, iny+4);
		  ypos[iny+1]=rymin-dyplus;
		  ypos[iny+2]=rymin-2.0*dyplus;
		  ypos[iny+3]=rymin-3.0*dyplus;
		  iny=iny+3;
          //BubbleEnhSort<doublereal>(ypos,0, iny);
		  Sort_method<doublereal>(ypos, iny);
	   }
	}

	for (i=0; i<=iny; i++) if (fabs(ypos[i]-rymax)<admission) { 
		iup=i;
		break; // досрочный выход из цикла for
	}

    if (iup!=iny) {

       dyminus=ypos[iup]-ypos[iup-1];
	   dyplus=ypos[iup+1]-ypos[iup];
       if (dyplus>3*dyminus) {
          SetLength(ypos,iny+1, iny+4);
		  ypos[iny+1]=rymax+dyminus;
		  ypos[iny+2]=rymax+2.0*dyminus;
		  ypos[iny+3]=rymax+3.0*dyminus;
		  iny=iny+3;
          //BubbleEnhSort<doublereal>(ypos, 0, iny);
		  Sort_method<doublereal>(ypos, iny);
	   }
	}

    for (i=0; i<=iny; i++) if (fabs(ypos[i]-rymin)<admission) { 
		idouwn=i;
		break; // досрочный выход из цикла for
	}

	for (i=0; i<=iny; i++) if (fabs(ypos[i]-rymax)<admission) { 
		iup=i;
		break; // досрочный выход из цикла for
	}

    if ((idouwn!=0)&&(iup!=iny)) {
       SetLength(ypos, iny+1, iny+5);
	   ypos[iny+1]=0.5*(ypos[idouwn-1]+ypos[idouwn]);
	   ypos[iny+2]=0.5*(ypos[idouwn]+ypos[idouwn+1]);
       ypos[iny+3]=0.5*(ypos[iup-1]+ypos[iup]);
	   ypos[iny+4]=0.5*(ypos[iup]+ypos[iup+1]);
	   iny=iny+4;
	}
    //BubbleEnhSort<doublereal>(ypos,0, iny);
	Sort_method<doublereal>(ypos, iny);

} // simplecorrect_meshgen_y

// Выдаёт true если теплопроводность блока ib1 больше 
// теплопроводности блока ib2.
bool comparison_lam(TPROP* matlist, BLOCK* b, integer ib1, integer ib2, doublereal t_C) {
	if ((matlist[b[ib1].imatid].blibmat == 0) && (matlist[b[ib2].imatid].blibmat == 0)) {
		//doublereal lam1 = matlist[b[ib1].imatid].lam;
		//doublereal lam2 = matlist[b[ib2].imatid].lam;
		float lam1 = get_lam(matlist[b[ib1].imatid].n_lam, matlist[b[ib1].imatid].temp_lam, matlist[b[ib1].imatid].arr_lam, t_C);
		float lam2 = get_lam(matlist[b[ib2].imatid].n_lam, matlist[b[ib2].imatid].temp_lam, matlist[b[ib2].imatid].arr_lam, t_C);
		if (lam1 > lam2) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else if ((matlist[b[ib1].imatid].blibmat == 1) && (matlist[b[ib2].imatid].blibmat == 0)) {
		float lam1=0.025f, rho=1.1614f, cp=1005.0f, Tnode=20.0;
		// библиотечный, внутрипрограммный материал.
		my_solid_properties(Tnode, rho, cp, lam1, matlist[b[ib1].imatid].ilibident);
		
		//doublereal lam2 = matlist[b[ib2].imatid].lam;
		float lam2= get_lam(matlist[b[ib2].imatid].n_lam, matlist[b[ib2].imatid].temp_lam, matlist[b[ib2].imatid].arr_lam, t_C);
		if (lam1 > lam2) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else if ((matlist[b[ib1].imatid].blibmat == 0) && (matlist[b[ib2].imatid].blibmat == 1)) {
		float lam2=0.025f, rho=1.1614f, cp=1005.0f, Tnode=20.0;
		// библиотечный, внутрипрограммный материал.
		my_solid_properties(Tnode, rho, cp, lam2, matlist[b[ib2].imatid].ilibident);

		//doublereal lam1 = matlist[b[ib1].imatid].lam;
		float lam1 = get_lam(matlist[b[ib1].imatid].n_lam, matlist[b[ib1].imatid].temp_lam, matlist[b[ib1].imatid].arr_lam, t_C);
		if (lam1 > lam2) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		float lam1=0.025f, rho1=1.1614f, cp1=1005.0f,Tnode=20.0;
		// библиотечный, внутрипрограммный материал.
		my_solid_properties(Tnode, rho1, cp1, lam1, matlist[b[ib1].imatid].ilibident);
		float lam2=0.025f, rho2=1.1614f, cp2=1005.0f;
		// библиотечный, внутрипрограммный материал.
		my_solid_properties(Tnode, rho2, cp2, lam2, matlist[b[ib2].imatid].ilibident);
		if (lam1 > lam2) {
			return 1;
		}
		else {
			return 0;
		}
	}
} // comparison_lam


void additional_line(doublereal ratio_quality, int inx, int iny, int inz,
	bool &bcont, integer *& i_1, integer *& j_1, integer *& k_1, integer ic9, integer &icorrect_stat,
	doublereal *&xpos, doublereal *&ypos, doublereal *&zpos) {

	for (integer ic7 = 0; ic7<ic9; ic7++)
	{

		//doublereal dx = 1, dy = 1, dz = 1;
		doublereal dx = xpos[i_1[ic7] + 1] - xpos[i_1[ic7]];
		doublereal dy = ypos[j_1[ic7] + 1] - ypos[j_1[ic7]];
		doublereal dz = zpos[k_1[ic7] + 1] - zpos[k_1[ic7]];


		if ((dx >= dy) && (dy >= dz)) {
			if (dx / dz > ratio_quality) {

				bool badd = true;
				for (integer i97 = 0; i97 <= inx; i97++) if (fabs(xpos[i97] - (xpos[i_1[ic7]] + 0.5*dx)) < 1.0e-36) badd = false;

				if (badd) {
					bcont = true;
					SetLength(xpos, inx + 1, inx + 2);
					//for (integer l = inx; l >= i_1[ic7] + 1; l--) xpos[l + 1] = xpos[l];
					//xpos[inx + 1] = xpos[i_1[ic7]] + 0.5*dx;
					//xpos[i_1[ic7] + 1] = xpos[i_1[ic7]] + 0.5*dx;
					// добавляем в конец
					xpos[inx + 1] = xpos[i_1[ic7]] + 0.5*dx;
					inx = inx + 1;
					//BubbleEnhSort<doublereal>(xpos, 0, inx);
					//Sort_method<doublereal>(xpos, inx);
					icorrect_stat++;

					continue;
				}
			}
		}
		else if ((dx >= dz) && (dz >= dy)) {
			if (dx / dy > ratio_quality) {


				bool badd = true;
				for (integer i97 = 0; i97 <= inx; i97++) if (fabs(xpos[i97] - (xpos[i_1[ic7]] + 0.5*dx)) < 1.0e-36) badd = false;

				if (badd) {

					bcont = true;
					SetLength(xpos, inx + 1, inx + 2);
					//for (integer l = inx; l >= i_1[ic7] + 1; l--) xpos[l + 1] = xpos[l];
					//xpos[inx + 1] = xpos[i_1[ic7]] + 0.5*dx;
					//xpos[i_1[ic7] + 1] = xpos[i_1[ic7]] + 0.5*dx;
					// добавляем в конец
					xpos[inx + 1] = xpos[i_1[ic7]] + 0.5*dx;
					inx = inx + 1;
					//BubbleEnhSort<doublereal>(xpos,0, inx);
					//Sort_method<doublereal>(xpos, inx);
					icorrect_stat++;

					continue;
				}
			}
		}
		if ((dy >= dx) && (dx >= dz)) {
			if (dy / dz > ratio_quality) {


				bool badd = true;
				for (integer i97 = 0; i97 <= iny; i97++) if (fabs(ypos[i97] - (ypos[j_1[ic7]] + 0.5*dy)) < 1.0e-36) badd = false;

				if (badd) {

					bcont = true;
					SetLength(ypos, iny + 1, iny + 2);
					//for (integer l = iny; l >= j_1[ic7] + 1; l--) ypos[l + 1] = ypos[l];
					//ypos[iny + 1] = ypos[j_1[ic7]] + 0.5*dy;
					//ypos[j_1[ic7] + 1] = ypos[j_1[ic7]] + 0.5*dy;
					// добавляем в конец
					ypos[iny + 1] = ypos[j_1[ic7]] + 0.5*dy;
					iny = iny + 1;
					//BubbleEnhSort<doublereal>(ypos,0, iny);
					//Sort_method<doublereal>(ypos, iny);
					icorrect_stat++;

					continue;
				}
			}
		}
		else if ((dy >= dz) && (dz >= dx)) {
			if (dy / dx > ratio_quality) {

				bool badd = true;
				for (integer i97 = 0; i97 <= iny; i97++) if (fabs(ypos[i97] - (ypos[j_1[ic7]] + 0.5*dy)) < 1.0e-36) badd = false;

				if (badd) {

					bcont = true;
					SetLength(ypos, iny + 1, iny + 2);
					//for (integer l = iny; l >= j_1[ic7] + 1; l--) ypos[l + 1] = ypos[l];
					//ypos[iny + 1] = ypos[j_1[ic7]] + 0.5*dy;
					//ypos[j_1[ic7] + 1] = ypos[j_1[ic7]] + 0.5*dy;
					// добавляем в конец
					ypos[iny + 1] = ypos[j_1[ic7]] + 0.5*dy;
					iny = iny + 1;
					//BubbleEnhSort<doublereal>(ypos, 0, iny);
					//Sort_method<doublereal>(ypos, iny);
					icorrect_stat++;

					continue;
				}
			}
		}
		if ((dz >= dy) && (dy >= dx)) {
			if (dz / dx > ratio_quality) {

				bool badd = true;
				for (integer i97 = 0; i97 <= inz; i97++) if (fabs(zpos[i97] - (zpos[k_1[ic7]] + 0.5*dz)) < 1.0e-36) badd = false;

				if (badd) {

					bcont = true;
					SetLength(zpos, inz + 1, inz + 2);
					//for (integer l = inz; l >= k_1[ic7] + 1; l--) zpos[l + 1] = zpos[l];
					//zpos[inz + 1] = zpos[k_1[ic7]] + 0.5*dz;
					//zpos[k_1[ic7] + 1] = zpos[k_1[ic7]] + 0.5*dz;
					// добавляем в конец
					zpos[inz + 1] = zpos[k_1[ic7]] + 0.5*dz;
					inz = inz + 1;
					//BubbleEnhSort<doublereal>(zpos, 0, inz);
					//Sort_method<doublereal>(zpos, inz);
					icorrect_stat++;

					continue;
				}
			}
		}
		else if ((dz >= dx) && (dx >= dy)) {
			if (dz / dy > ratio_quality) {

				bool badd = true;
				for (integer i97 = 0; i97 <= inz; i97++) if (fabs(zpos[i97] - (zpos[k_1[ic7]] + 0.5*dz)) < 1.0e-36) badd = false;

				if (badd) {

					bcont = true;
					SetLength(zpos, inz + 1, inz + 2);
					//for (integer l = inz; l >= k_1[ic7] + 1; l--) zpos[l + 1] = zpos[l];
					//zpos[inz + 1] = zpos[k_1[ic7]] + 0.5*dz;
					//zpos[k_1[ic7] + 1] = zpos[k_1[ic7]] + 0.5*dz;
					// добавляем в конец
					zpos[inz + 1] = zpos[k_1[ic7]] + 0.5*dz;
					inz = inz + 1;
					//BubbleEnhSort<doublereal>(zpos, 0, inz);
					//Sort_method<doublereal>(zpos, inz);
					icorrect_stat++;

					continue;
				}
			}
		}


	}

}

void marker_for_additional_line(integer &ic9_shadow,
	integer &ic9, int inx, int iny, int inz,
	doublereal *&xpos, doublereal *&ypos, doublereal *&zpos,
	doublereal ratio_quality, const integer iSIZE_DIRECTIONAL,
	bool *& i_b, bool *& j_b, bool *& k_b,
	integer *& i_1, integer *& j_1, integer *& k_1,
	doublereal &ratio_start_check) {

	bool restart = false;
	

START_LABEL_FOR_ADAPT:

	restart = false;
	bool bfinish = false;

	int irandom = rand() % 3;
	ic9_shadow = ic9;

	if (irandom == 0) {
#pragma omp parallel for
		for (int i_1l = 0; i_1l < inx; i_1l++) {
			for (int j_1l = 0; j_1l < iny; j_1l++) {
				for (int k_1l = 0; k_1l < inz; k_1l++) {
					if ((!restart)&&(!bfinish)) {
					// сокращает время построения сетки.

					TOCHKA p;
					p.x = 0.5*(xpos[i_1l + 1] + xpos[i_1l]);
					p.y = 0.5*(ypos[j_1l + 1] + ypos[j_1l]);
					p.z = 0.5*(zpos[k_1l + 1] + zpos[k_1l]);

					doublereal dx = xpos[i_1l + 1] - xpos[i_1l];
					doublereal dy = ypos[j_1l + 1] - ypos[j_1l];
					doublereal dz = zpos[k_1l + 1] - zpos[k_1l];
					// Порядок узлов тоже важен. Модуль не нужен.
#if doubleintprecision == 1
					//if (dx < shorter_length_for_simplificationX(p.x)) {
					if (dx < 1.0e-36) {
						printf("Error: Slipanie geometrii...dX \n");
						//printf("x[%lld]=%e x[%lld]=%e x[%lld]=%e\n", i_1l - 1, xpos[i_1l - 1], i_1l, xpos[i_1l], i_1l + 1, xpos[i_1l + 1]);
						std::cout << "x[" << i_1l - 1 << "]=" << xpos[i_1l - 1] << " x[" << i_1l << "]=" << xpos[i_1l] << " x[" << i_1l + 1 << "]=" << xpos[i_1l + 1] << std::endl;
						system("pause");
					}
					//if (dy < shorter_length_for_simplificationY(p.y)) {
					if (dy < 1.0e-36) {
						printf("Error: Slipanie geometrii...dY \n");
						//printf("y[%lld]=%e y[%lld]=%e y[%lld]=%e\n", j_1l - 1, ypos[j_1l - 1], j_1l, ypos[j_1l], j_1l + 1, ypos[j_1l + 1]);
						std::cout << "y[" << j_1l - 1 << "]=" << ypos[j_1l - 1] << " y[" << j_1l << "]=" << ypos[j_1l] << " y[" << j_1l + 1 << "]=" << ypos[j_1l + 1] << std::endl;
						system("pause");
					}
					//if (dz < shorter_length_for_simplificationZ(p.z)) {
					if (dz < 1.0e-36) {
						printf("Error: Slipanie geometrii...dZ \n");
						//printf("z[%lld]=%e z[%lld]=%e z[%lld]=%e\n", k_1l - 1, zpos[k_1l - 1], k_1l, zpos[k_1l], k_1l + 1, zpos[k_1l + 1]);
						std::cout << "z[" << k_1l - 1 << "]=" << zpos[k_1l - 1] << " z[" << k_1l << "]=" << zpos[k_1l] << " z[" << k_1l + 1 << "]=" << zpos[k_1l + 1] << std::endl;
						system("pause");
					}
#else
					//if (dx < shorter_length_for_simplificationX(p.x)) {
					if (dx < 1.0e-36) {
						printf("Error: Slipanie geometrii...dX \n");
						printf("x[%d]=%e x[%d]=%e x[%d]=%e\n", i_1l - 1, xpos[i_1l - 1], i_1l, xpos[i_1l], i_1l + 1, xpos[i_1l + 1]);
						system("pause");
					}
					//if (dy < shorter_length_for_simplificationY(p.y)) {
					if (dy < 1.0e-36) {
						printf("Error: Slipanie geometrii...dY \n");
						printf("y[%d]=%e y[%d]=%e y[%d]=%e\n", j_1l - 1, ypos[j_1l - 1], j_1l, ypos[j_1l], j_1l + 1, ypos[j_1l + 1]);
						system("pause");
					}
					//if (dz < shorter_length_for_simplificationZ(p.z)) {
					if (dz < 1.0e-36) {
						printf("Error: Slipanie geometrii...dZ \n");
						printf("z[%d]=%e z[%d]=%e z[%d]=%e\n", k_1l - 1, zpos[k_1l - 1], k_1l, zpos[k_1l], k_1l + 1, zpos[k_1l + 1]);
						system("pause");
					}
#endif

					if ((dx >= dy) && (dy >= dz)) {
						if (dx / dz > ratio_quality) {
							if (i_b[i_1l] == false) {
								if (fabs(dz) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dx / dz;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										i_b[i_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
											
										}
									}
								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#endif
									system("pause");
								}
							}
						}
					}
					if ((dx >= dz) && (dz >= dy)) {
						if (dx / dy > ratio_quality) {
							if (i_b[i_1l] == false) {
								if (fabs(dy) > 1.0e-36) {

#pragma omp critical 
									{

										ratio_start_check = dx / dy;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										i_b[i_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
											
										}
									}
								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#endif
									system("pause");
								}
							}
						}
					}
					if ((dy >= dx) && (dx >= dz)) {
						if (dy / dz > ratio_quality) {
							if (j_b[j_1l] == false) {
								if (fabs(dz) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dy / dz;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										j_b[j_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
											
										}
									}
								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#endif
									system("pause");
								}
							}
						}
					}
					if ((dy >= dz) && (dz >= dx)) {
						if (dy / dx > ratio_quality) {
							if (j_b[j_1l] == false) {
								if (fabs(dx) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dy / dx;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										j_b[j_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
											
										}
									}

								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#endif
									system("pause");
								}
							}
						}
					}
					if ((dz >= dx) && (dx >= dy)) {
						if (dz / dy > ratio_quality) {
							if (k_b[k_1l] == false) {
								if (fabs(dy) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dz / dy;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										k_b[k_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
										
										}
									}

								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#endif
									system("pause");
								}
							}
						}
					}
					if ((dz >= dy) && (dy >= dx)) {
						if (dz / dx > ratio_quality) {
							if (k_b[k_1l] == false) {
								if (fabs(dx) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dz / dx;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										k_b[k_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
											
										}
									}

								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;


#endif
									system("pause");
								}
							}
						}
					}
					if (ic9 > 93) {
						// STOP process.
#pragma omp critical 
						{
							bfinish = true;
						}
					}
				}

				}
			}
		}
	}
	else if (irandom == 1) {

#pragma omp parallel for
		for (int j_1l = 0; j_1l < iny; j_1l++) {
			for (int k_1l = 0; k_1l < inz; k_1l++) {
				for (int i_1l = 0; i_1l < inx; i_1l++) {
					if ((!restart) && (!bfinish)) {
						// сокращает время построения сетки.

						TOCHKA p;
						p.x = 0.5 * (xpos[i_1l + 1] + xpos[i_1l]);
						p.y = 0.5 * (ypos[j_1l + 1] + ypos[j_1l]);
						p.z = 0.5 * (zpos[k_1l + 1] + zpos[k_1l]);

						doublereal dx = xpos[i_1l + 1] - xpos[i_1l];
						doublereal dy = ypos[j_1l + 1] - ypos[j_1l];
						doublereal dz = zpos[k_1l + 1] - zpos[k_1l];
						// Порядок узлов тоже важен. Модуль не нужен.
#if doubleintprecision == 1
					//if (dx < shorter_length_for_simplificationX(p.x)) {
						if (dx < 1.0e-36) {
							printf("Error: Slipanie geometrii...dX \n");
							//printf("x[%lld]=%e x[%lld]=%e x[%lld]=%e\n", i_1l - 1, xpos[i_1l - 1], i_1l, xpos[i_1l], i_1l + 1, xpos[i_1l + 1]);
							std::cout << "x[" << i_1l - 1 << "]=" << xpos[i_1l - 1] << " x[" << i_1l << "]=" << xpos[i_1l] << " x[" << i_1l + 1 << "]=" << xpos[i_1l + 1] << std::endl;
							system("pause");
						}
						//if (dy < shorter_length_for_simplificationY(p.y)) {
						if (dy < 1.0e-36) {
							printf("Error: Slipanie geometrii...dY \n");
							//printf("y[%lld]=%e y[%lld]=%e y[%lld]=%e\n", j_1l - 1, ypos[j_1l - 1], j_1l, ypos[j_1l], j_1l + 1, ypos[j_1l + 1]);
							std::cout << "y[" << j_1l - 1 << "]=" << ypos[j_1l - 1] << " y[" << j_1l << "]=" << ypos[j_1l] << " y[" << j_1l + 1 << "]=" << ypos[j_1l + 1] << std::endl;
							system("pause");
						}
						//if (dz < shorter_length_for_simplificationZ(p.z)) {
						if (dz < 1.0e-36) {
							printf("Error: Slipanie geometrii...dZ \n");
							//printf("z[%lld]=%e z[%lld]=%e z[%lld]=%e\n", k_1l - 1, zpos[k_1l - 1], k_1l, zpos[k_1l], k_1l + 1, zpos[k_1l + 1]);
							std::cout << "z[" << k_1l - 1 << "]=" << zpos[k_1l - 1] << " z[" << k_1l << "]=" << zpos[k_1l] << " z[" << k_1l + 1 << "]=" << zpos[k_1l + 1] << std::endl;
							system("pause");
						}
#else
					//if (dx < shorter_length_for_simplificationX(p.x)) {
						if (dx < 1.0e-36) {
							printf("Error: Slipanie geometrii...dX \n");
							printf("x[%d]=%e x[%d]=%e x[%d]=%e\n", i_1l - 1, xpos[i_1l - 1], i_1l, xpos[i_1l], i_1l + 1, xpos[i_1l + 1]);
							system("pause");
						}
						//if (dy < shorter_length_for_simplificationY(p.y)) {
						if (dy < 1.0e-36) {
							printf("Error: Slipanie geometrii...dY \n");
							printf("y[%d]=%e y[%d]=%e y[%d]=%e\n", j_1l - 1, ypos[j_1l - 1], j_1l, ypos[j_1l], j_1l + 1, ypos[j_1l + 1]);
							system("pause");
						}
						//if (dz < shorter_length_for_simplificationZ(p.z)) {
						if (dz < 1.0e-36) {
							printf("Error: Slipanie geometrii...dZ \n");
							printf("z[%d]=%e z[%d]=%e z[%d]=%e\n", k_1l - 1, zpos[k_1l - 1], k_1l, zpos[k_1l], k_1l + 1, zpos[k_1l + 1]);
							system("pause");
						}
#endif

						if ((dx >= dy) && (dy >= dz)) {
							if (dx / dz > ratio_quality) {
								if (i_b[i_1l] == false) {
									if (fabs(dz) > 1.0e-36) {

#pragma omp critical 
										{

											ratio_start_check = dx / dz;
											i_1[ic9] = i_1l;
											j_1[ic9] = j_1l;
											k_1[ic9] = k_1l;
											i_b[i_1l] = true;
											ic9++;

											if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
												//goto START_LABEL_FOR_ADAPT;

												restart = true;

											}
										}
									}
									else {
#if doubleintprecision == 1
										//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#else
										//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#endif
										system("pause");
									}
						}
							}
						}
						if ((dx >= dz) && (dz >= dy)) {
							if (dx / dy > ratio_quality) {
								if (i_b[i_1l] == false) {
									if (fabs(dy) > 1.0e-36) {
#pragma omp critical 
										{

											ratio_start_check = dx / dy;
											i_1[ic9] = i_1l;
											j_1[ic9] = j_1l;
											k_1[ic9] = k_1l;
											i_b[i_1l] = true;
											ic9++;

											if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
												//goto START_LABEL_FOR_ADAPT;

												restart = true;

											}
										}
									}
									else {
#if doubleintprecision == 1
										//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#else
										//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#endif
										system("pause");
									}
								}
							}
					}
						if ((dy >= dx) && (dx >= dz)) {
							if (dy / dz > ratio_quality) {
								if (j_b[j_1l] == false) {
									if (fabs(dz) > 1.0e-36) {
#pragma omp critical 
										{

											ratio_start_check = dy / dz;
											i_1[ic9] = i_1l;
											j_1[ic9] = j_1l;
											k_1[ic9] = k_1l;
											j_b[j_1l] = true;
											ic9++;

											if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
												//goto START_LABEL_FOR_ADAPT;

												restart = true;

											}
										}
									}
									else {
#if doubleintprecision == 1
										//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#else
										//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#endif
										system("pause");
									}
				}
			}
		}
						if ((dy >= dz) && (dz >= dx)) {
							if (dy / dx > ratio_quality) {
								if (j_b[j_1l] == false) {
									if (fabs(dx) > 1.0e-36) {
#pragma omp critical 
										{

											ratio_start_check = dy / dx;
											i_1[ic9] = i_1l;
											j_1[ic9] = j_1l;
											k_1[ic9] = k_1l;
											j_b[j_1l] = true;
											ic9++;

											if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
												//goto START_LABEL_FOR_ADAPT;

												restart = true;

											}
										}

							}
									else {
#if doubleintprecision == 1
										//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#else
										//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#endif
										system("pause");
									}
						}
	}
		}
						if ((dz >= dx) && (dx >= dy)) {
							if (dz / dy > ratio_quality) {
								if (k_b[k_1l] == false) {
									if (fabs(dy) > 1.0e-36) {

#pragma omp critical 
										{
											ratio_start_check = dz / dy;
											i_1[ic9] = i_1l;
											j_1[ic9] = j_1l;
											k_1[ic9] = k_1l;
											k_b[k_1l] = true;
											ic9++;

											if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
												//goto START_LABEL_FOR_ADAPT;

												restart = true;

											}
								}

							}
									else {
#if doubleintprecision == 1
										//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#else
										//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#endif
										system("pause");
									}
						}
	}
}
						if ((dz >= dy) && (dy >= dx)) {
							if (dz / dx > ratio_quality) {
								if (k_b[k_1l] == false) {
									if (fabs(dx) > 1.0e-36) {
#pragma omp critical 
										{

											ratio_start_check = dz / dx;
											i_1[ic9] = i_1l;
											j_1[ic9] = j_1l;
											k_1[ic9] = k_1l;
											k_b[k_1l] = true;
											ic9++;

											if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
												//goto START_LABEL_FOR_ADAPT;

												restart = true;

											}
									}

								}
									else {
#if doubleintprecision == 1
										//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#else
										//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
										std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;


#endif
										system("pause");
									}
							}
						}
					}
						if (ic9 > 93) {
							// STOP process.
#pragma omp critical 
							{
								bfinish = true;
							}
						}
					}

				}
			}
		}
	}
	else {
#pragma omp parallel for
		for (int k_1l = 0; k_1l < inz; k_1l++) {
			for (int i_1l = 0; i_1l < inx; i_1l++) {
				for (int j_1l = 0; j_1l < iny; j_1l++) {
					if ((!restart) && (!bfinish)) {
					// сокращает время построения сетки.

					TOCHKA p;
					p.x = 0.5*(xpos[i_1l + 1] + xpos[i_1l]);
					p.y = 0.5*(ypos[j_1l + 1] + ypos[j_1l]);
					p.z = 0.5*(zpos[k_1l + 1] + zpos[k_1l]);

					doublereal dx = xpos[i_1l + 1] - xpos[i_1l];
					doublereal dy = ypos[j_1l + 1] - ypos[j_1l];
					doublereal dz = zpos[k_1l + 1] - zpos[k_1l];
					// Порядок узлов тоже важен. Модуль не нужен.
#if doubleintprecision == 1
					//if (dx < shorter_length_for_simplificationX(p.x)) {
					if (dx < 1.0e-36) {
						printf("Error: Slipanie geometrii...dX \n");
						//printf("x[%lld]=%e x[%lld]=%e x[%lld]=%e\n", i_1l - 1, xpos[i_1l - 1], i_1l, xpos[i_1l], i_1l + 1, xpos[i_1l + 1]);
						std::cout << "x[" << i_1l - 1 << "]=" << xpos[i_1l - 1] << " x[" << i_1l << "]=" << xpos[i_1l] << " x[" << i_1l + 1 << "]=" << xpos[i_1l + 1] << std::endl;
						system("pause");
					}
					//if (dy < shorter_length_for_simplificationY(p.y)) {
					if (dy < 1.0e-36) {
						printf("Error: Slipanie geometrii...dY \n");
						//printf("y[%lld]=%e y[%lld]=%e y[%lld]=%e\n", j_1l - 1, ypos[j_1l - 1], j_1l, ypos[j_1l], j_1l + 1, ypos[j_1l + 1]);
						std::cout << "y[" << j_1l - 1 << "]=" << ypos[j_1l - 1] << " y[" << j_1l << "]=" << ypos[j_1l] << " y[" << j_1l + 1 << "]=" << ypos[j_1l + 1] << std::endl;
						system("pause");
					}
					//if (dz < shorter_length_for_simplificationZ(p.z)) {
					if (dz < 1.0e-36) {
						printf("Error: Slipanie geometrii...dZ \n");
						//printf("z[%lld]=%e z[%lld]=%e z[%lld]=%e\n", k_1l - 1, zpos[k_1l - 1], k_1l, zpos[k_1l], k_1l + 1, zpos[k_1l + 1]);
						std::cout << "z[" << k_1l - 1 << "]=" << zpos[k_1l - 1] << " z[" << k_1l << "]=" << zpos[k_1l] << " z[" << k_1l + 1 << "]=" << zpos[k_1l + 1] << std::endl;
						system("pause");
					}
#else
					//if (dx < shorter_length_for_simplificationX(p.x)) {
					if (dx < 1.0e-36) {
						printf("Error: Slipanie geometrii...dX \n");
						printf("x[%d]=%e x[%d]=%e x[%d]=%e\n", i_1l - 1, xpos[i_1l - 1], i_1l, xpos[i_1l], i_1l + 1, xpos[i_1l + 1]);
						system("pause");
					}
					//if (dy < shorter_length_for_simplificationY(p.y)) {
					if (dy < 1.0e-36) {
						printf("Error: Slipanie geometrii...dY \n");
						printf("y[%d]=%e y[%d]=%e y[%d]=%e\n", j_1l - 1, ypos[j_1l - 1], j_1l, ypos[j_1l], j_1l + 1, ypos[j_1l + 1]);
						system("pause");
					}
					//if (dz < shorter_length_for_simplificationZ(p.z)) {
					if (dz < 1.0e-36) {
						printf("Error: Slipanie geometrii...dZ \n");
						printf("z[%d]=%e z[%d]=%e z[%d]=%e\n", k_1l - 1, zpos[k_1l - 1], k_1l, zpos[k_1l], k_1l + 1, zpos[k_1l + 1]);
						system("pause");
					}
#endif

					if ((dx >= dy) && (dy >= dz)) {
						if (dx / dz > ratio_quality) {
							if (i_b[i_1l] == false) {
								if (fabs(dz) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dx / dz;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										i_b[i_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
										}
									}
								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#endif
									system("pause");
								}
							}
						}
					}
					if ((dx >= dz) && (dz >= dy)) {
						if (dx / dy > ratio_quality) {
							if (i_b[i_1l] == false) {
								if (fabs(dy) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dx / dy;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										i_b[i_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
											
										}
									}
								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#endif
									system("pause");
								}
							}
						}
					}
					if ((dy >= dx) && (dx >= dz)) {
						if (dy / dz > ratio_quality) {
							if (j_b[j_1l] == false) {
								if (fabs(dz) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dy / dz;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										j_b[j_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
											
										}
									}
								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
#endif
									system("pause");
								}
							}
						}
					}
					if ((dy >= dz) && (dz >= dx)) {
						if (dy / dx > ratio_quality) {
							if (j_b[j_1l] == false) {
								if (fabs(dx) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dy / dx;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										j_b[j_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
										
										}
									}

								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#endif
									system("pause");
								}
							}
						}
					}
					if ((dz >= dx) && (dx >= dy)) {
						if (dz / dy > ratio_quality) {
							if (k_b[k_1l] == false) {
								if (fabs(dy) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dz / dy;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										k_b[k_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;
											
										}
									}

								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#endif
									system("pause");
								}
							}
						}
					}
					if ((dz >= dy) && (dy >= dx)) {
						if (dz / dx > ratio_quality) {
							if (k_b[k_1l] == false) {
								if (fabs(dx) > 1.0e-36) {
#pragma omp critical 
									{

										ratio_start_check = dz / dx;
										i_1[ic9] = i_1l;
										j_1[ic9] = j_1l;
										k_1[ic9] = k_1l;
										k_b[k_1l] = true;
										ic9++;

										if (ic9 - ic9_shadow > iSIZE_DIRECTIONAL) {
											//goto START_LABEL_FOR_ADAPT;

											restart = true;

										}
									}

								}
								else {
#if doubleintprecision == 1
									//printf("i=%lld, j=%lld, k=%lld, inx=%lld, iny=%lld, inz=%lld, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;

#else
									//printf("i=%d, j=%d, k=%d, inx=%d, iny=%d, inz=%d, dx=%e, dy=%e, dz=%e\n", i_1l, j_1l, k_1l, inx, iny, inz, dx, dy, dz);
									std::cout << "i=" << i_1l << ", j=" << j_1l << ", k=" << k_1l << ", inx=" << inx << ", iny=" << iny << ", inz=" << inz << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;


#endif
									system("pause");
								}
										}
									}
								}
					if (ic9 > 93) {
						// STOP process.
#pragma omp critical 
						{

							bfinish = true;
						}
					}
					}

				}
			}
		}
	}

	if (restart) {
		goto START_LABEL_FOR_ADAPT;
	}

} // marker_for_additional_line

void cycle_additional_line(integer &iteration_number, integer &ic4,
	bool &bcont, doublereal &ratio_start_check,
	int inx, int iny, int inz,
	doublereal *& xpos, doublereal *& ypos, doublereal *& zpos,
	const doublereal ratio_quality, integer &icorrect_stat) {

START_LAB:

#if doubleintprecision == 1
	//printf("%lld ",0);
#else
	//printf("%d ",0);
#endif

	ic4++;
	ratio_start_check = 0.0;
	bcont = false;

	integer *i_1 = nullptr, *j_1 = nullptr, *k_1 = nullptr;
	i_1 = new integer[100];
	j_1 = new integer[100];
	k_1 = new integer[100];

	bool *i_b = nullptr, *j_b = nullptr, *k_b = nullptr;
	i_b = new bool[inx];
	j_b = new bool[iny];
	k_b = new bool[inz];

#pragma omp parallel for
	for (integer i_1l = 0; i_1l < inx; i_1l++)
	{
		i_b[i_1l] = false;
	}
#pragma omp parallel for
	for (integer j_1l = 0; j_1l < iny; j_1l++)
	{
		j_b[j_1l] = false;
	}
#pragma omp parallel for
	for (integer k_1l = 0; k_1l < inz; k_1l++)
	{
		k_b[k_1l] = false;
	}

	srand((unsigned int)(time(0)));

	const integer iSIZE_DIRECTIONAL = 10;
	integer ic9 = 0;
	integer ic9_shadow = ic9;


	marker_for_additional_line(ic9_shadow,
		ic9, inx, iny, inz, xpos, ypos, zpos,
		ratio_quality, iSIZE_DIRECTIONAL,
		i_b, j_b, k_b, i_1, j_1, k_1,
		ratio_start_check);


	delete[] i_b;
	delete[] j_b;
	delete[] k_b;

#if doubleintprecision == 1
	//printf("ic9=%lld %e ", ic9, ratio_start_check);
	std::cout << iteration_number << "  ic9=" << ic9 << " ratio_check=" << ratio_start_check << std::endl;
#else
	//printf("ic9=%d %e ", ic9, ratio_start_check);
	std::cout << iteration_number << "  ic9=" << ic9 << " ratio_check=" << ratio_start_check << std::endl;
#endif


	additional_line(ratio_quality, inx, iny, inz, bcont, i_1, j_1, k_1, ic9, icorrect_stat, xpos, ypos, zpos);

	delete[] i_1;
	delete[] j_1;
	delete[] k_1;

	//BubbleEnhSort<doublereal>(xpos, 0, inx);
	//BubbleEnhSort<doublereal>(ypos, 0, iny);
	//BubbleEnhSort<doublereal>(zpos, 0, inz);
#pragma omp parallel sections
	{
#pragma omp section
		{
			Sort_method<doublereal>(xpos, inx);
		}
#pragma omp section
		{
			Sort_method<doublereal>(ypos, iny);
		}
#pragma omp section
		{
			Sort_method<doublereal>(zpos, inz);
		}
	}
	if (bcont) {

		iteration_number++;

		// Защита от бесконечного цикла на полигональных моделях.
		if (iteration_number < 40) {
			goto START_LAB;
		}
	}
}


  // 11.02.2017
void quolite_refinement(int &inx, int &iny, int &inz, doublereal* &xpos, doublereal* &ypos, doublereal* &zpos) {

#if doubleintprecision == 1
	printf("apriory mesh size: inx=%d iny=%d inz=%d\n", inx, iny, inz);
#else
	printf("apriory mesh size: inx=%d iny=%d inz=%d\n", inx, iny, inz);
#endif

	// check_mesh
	doublereal ratio_start_check = 0.0;
#pragma omp parallel for
	for (int i_1 = 0; i_1 < inx; ++i_1) {
		for (int j_1 = 0; j_1 < iny; ++j_1) {
			for (int k_1 = 0; k_1 < inz; ++k_1) {
				
				doublereal dx = xpos[i_1 + 1] - xpos[i_1];
				doublereal dy = ypos[j_1 + 1] - ypos[j_1];
				doublereal dz = zpos[k_1 + 1] - zpos[k_1];
				if ((dx >= dy) && (dy >= dz)) {
					if (dx / dz > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dx / dz;
						}
					}
				}
				if ((dx >= dz) && (dz >= dy)) {
					if (dx / dy > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dx / dy;
						}
					}
				}
				if ((dy >= dx) && (dx >= dz)) {
					if (dy / dz > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dy / dz;
						}
					}
				}
				if ((dy >= dz) && (dz >= dx)) {
					if (dy / dx > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dy / dx;
						}
					}
				}
				if ((dz >= dx) && (dx >= dy)) {
					if (dz / dy > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dz / dy;
						}
					}
				}
				if ((dz >= dy) && (dy >= dx)) {
					if (dz / dx > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dz / dx;
						}
					}
				}
			}
		}
	}
	//printf("mesh apriority quolite=%e\n", ratio_start_check);
	std::cout << "mesh apriority quolite=" << ratio_start_check << std::endl;
	if (ratio_start_check > 5.0e6) {
		// Предупреждение о возможных проблемах сходимости вычислительного процесса.
		printf("WARNING!!! Your mesh apriority quolite is very big!!!\n");
		printf("May be your model is incorrect. Please, patch your model\n");
		system("pause");
	}


	// Автоматическая коррекция качества сетки.
	// 10_02_2017
	// (etalon_max_size_ratio2==30 optimum) это я вычитал во flowvision
	const doublereal ratio_quality = etalon_max_size_ratio2;
	integer icorrect_stat = 0;
	bool bcont = true;

	integer ic4 = 0;

	if (0) {
		while (bcont) {
			integer i_progress = 0;
		START_LAB1:
			bcont = false;
			for (int i_1 = 0; i_1 < inx; ++i_1) {
				for (int j_1 = 0; j_1 < iny; ++j_1) {
					for (int k_1 = 0; k_1 < inz; ++k_1) {
						// сокращает время построения сетки.
						if (i_1*iny*inz + j_1*inz + k_1 >= i_progress) {
							i_progress = i_1*iny*inz + j_1*inz + k_1;
							//doublereal dx = 1, dy = 1, dz = 1;
							doublereal dx = xpos[i_1 + 1] - xpos[i_1];
							doublereal dy = ypos[j_1 + 1] - ypos[j_1];
							doublereal dz = zpos[k_1 + 1] - zpos[k_1];
							if ((dx >= dy) && (dy >= dz)) {
								if (dx / dz > ratio_quality) {
									bcont = true;
									SetLength(xpos, inx + 1, inx + 2);
									for (integer l = inx; l >= i_1 + 1; l--) xpos[l + 1] = xpos[l];
									//xpos[inx + 1] = xpos[i_1] + 0.5*dx;
									xpos[i_1 + 1] = xpos[i_1] + 0.5*dx;
									inx = inx + 1;
									//BubbleEnhSort<doublereal>(xpos,0, inx);
									//Sort_method<doublereal>(xpos, inx);
									icorrect_stat++;
									goto START_LAB1;
								}
							}
							if ((dx >= dz) && (dz >= dy)) {
								if (dx / dy > ratio_quality) {
									bcont = true;
									SetLength(xpos, inx + 1, inx + 2);
									for (integer l = inx; l >= i_1 + 1; l--) xpos[l + 1] = xpos[l];
									//xpos[inx + 1] = xpos[i_1] + 0.5*dx;
									xpos[i_1 + 1] = xpos[i_1] + 0.5*dx;
									inx = inx + 1;
									//BubbleEnhSort<doublereal>(xpos,0, inx);
									//Sort_method<doublereal>(xpos,inx);
									icorrect_stat++;
									goto START_LAB1;
								}
							}
							if ((dy >= dx) && (dx >= dz)) {
								if (dy / dz > ratio_quality) {
									bcont = true;
									SetLength(ypos, iny + 1, iny + 2);
									for (integer l = iny; l >= j_1 + 1; l--) ypos[l + 1] = ypos[l];
									//ypos[iny + 1] = ypos[j_1] + 0.5*dy;
									ypos[j_1 + 1] = ypos[j_1] + 0.5*dy;
									iny = iny + 1;
									//BubbleEnhSort<doublereal>(ypos, 0,iny);
									//Sort_method<doublereal>(ypos,iny);
									icorrect_stat++;
									goto START_LAB1;
								}
							}
							if ((dy >= dz) && (dz >= dx)) {
								if (dy / dx > ratio_quality) {
									bcont = true;
									SetLength(ypos, iny + 1, iny + 2);
									for (integer l = iny; l >= j_1 + 1; l--) ypos[l + 1] = ypos[l];
									//ypos[iny + 1] = ypos[j_1] + 0.5*dy;
									ypos[j_1 + 1] = ypos[j_1] + 0.5*dy;
									iny = iny + 1;
									//BubbleEnhSort<doublereal>(ypos,0, iny);
									icorrect_stat++;
									goto START_LAB1;
								}
							}
							if ((dz >= dy) && (dy >= dx)) {
								if (dz / dx > ratio_quality) {
									bcont = true;
									SetLength(zpos, inz + 1, inz + 2);
									for (integer l = inz; l >= k_1 + 1; l--) zpos[l + 1] = zpos[l];
									//zpos[inz + 1] = zpos[k_1] + 0.5*dz;
									zpos[k_1 + 1] = zpos[k_1] + 0.5*dz;
									inz = inz + 1;
									//BubbleEnhSort<doublereal>(zpos, 0,inz);
									//Sort_method<doublereal>(zpos,inz);
									icorrect_stat++;
									goto START_LAB1;
								}
							}
							if ((dz >= dx) && (dx >= dy)) {
								if (dz / dy > ratio_quality) {
									bcont = true;
									SetLength(zpos, inz + 1, inz + 2);
									for (integer l = inz; l >= k_1 + 1; l--) zpos[l + 1] = zpos[l];
									//zpos[inz + 1] = zpos[k_1] + 0.5*dz;
									zpos[k_1 + 1] = zpos[k_1] + 0.5*dz;
									inz = inz + 1;
									//BubbleEnhSort<doublereal>(zpos, 0,inz);
									//Sort_method<doublereal>(zpos, inz);
									icorrect_stat++;
									goto START_LAB1;
								}
							}
						}
					}
				}
			}
		}
	}
	else if (1) {
		// Согласно этой стратегии мы сначала лечим наихудшие ячейки, а потом постепенно
		// улучшаем качество сетки и переходим к всё более и более лучшим.
		// Подходмаксимизации при котором в первую очередь лечатся наихудшие ячейки.
		// Даёт тоже самое разбиение.

		// Данная версия на момент 18 марта 2017 обладает наилучшим быстродействием, 
		// при сохранении корректности работы. В то время как версия при которой за проход добавляется лишь одна 
		// наилучшая линия сетки обладает очень плохим быстродействием, прога с ней работает неприемлемо долго.


		//while (bcont) {

		integer iteration_number = 0;

		cycle_additional_line(iteration_number, ic4,
			bcont, ratio_start_check, inx, iny, inz, xpos, ypos, zpos,
			ratio_quality, icorrect_stat);
	
		//}
	}
	else  {
		// Согласно этой стратегии мы сначала лечим наихудшие ячейки, а потом постепенно
		// улучшаем качество сетки и переходим к всё более и более лучшим.
		// Подходмаксимизации при котором в первую очередь лечатся наихудшие ячейки.
		// Даёт тоже самое разбиение.
		

		while (bcont) {

		START_LAB2:
#if doubleintprecision == 1
			printf("%d ", 0);
#else
			printf("%d ", 0);
#endif
			
			ic4++;
			ratio_start_check = 0.0;
			bcont = false;

			//integer i_1 = 0, j_1 = 0, k_1 = 0;
			// максиммально вытянутая ячейка.
			integer i_maximally_elongated_cell_1, j_maximally_elongated_cell_1, k_maximally_elongated_cell_1;

			for (int i_1l = 0; i_1l < inx; ++i_1l) {
				for (int j_1l = 0; j_1l < iny; ++j_1l) {
					for (int k_1l = 0; k_1l < inz; ++k_1l) {
						// сокращает время построения сетки.

						//doublereal dx = 1, dy = 1, dz = 1;
						doublereal dx = xpos[i_1l + 1] - xpos[i_1l];
						doublereal dy = ypos[j_1l + 1] - ypos[j_1l];
						doublereal dz = zpos[k_1l + 1] - zpos[k_1l];
						if ((dx >= dy) && (dy >= dz)) {
							if (dx / dz > ratio_start_check) {
								ratio_start_check = dx / dz;
								i_maximally_elongated_cell_1 = i_1l;
								j_maximally_elongated_cell_1 = j_1l;
								k_maximally_elongated_cell_1 = k_1l;
							}
						}
						if ((dx >= dz) && (dz >= dy)) {
							if (dx / dy > ratio_start_check) {
								ratio_start_check = dx / dy;
								i_maximally_elongated_cell_1 = i_1l;
								j_maximally_elongated_cell_1 = j_1l;
								k_maximally_elongated_cell_1 = k_1l;
							}
						}
						if ((dy >= dx) && (dx >= dz)) {
							if (dy / dz > ratio_start_check) {
								ratio_start_check = dy / dz;
								i_maximally_elongated_cell_1 = i_1l;
								j_maximally_elongated_cell_1 = j_1l;
								k_maximally_elongated_cell_1 = k_1l;
							}
						}
						if ((dy >= dz) && (dz >= dx)) {
							if (dy / dx > ratio_start_check) {
								ratio_start_check = dy / dx;
								i_maximally_elongated_cell_1 = i_1l;
								j_maximally_elongated_cell_1 = j_1l;
								k_maximally_elongated_cell_1 = k_1l;
							}
						}
						if ((dz >= dx) && (dx >= dy)) {
							if (dz / dy > ratio_start_check) {
								ratio_start_check = dz / dy;
								i_maximally_elongated_cell_1 = i_1l;
								j_maximally_elongated_cell_1 = j_1l;
								k_maximally_elongated_cell_1 = k_1l;
							}
						}
						if ((dz >= dy) && (dy >= dx)) {
							if (dz / dx > ratio_start_check) {
								ratio_start_check = dz / dx;
								i_maximally_elongated_cell_1 = i_1l;
								j_maximally_elongated_cell_1 = j_1l;
								k_maximally_elongated_cell_1 = k_1l;
							}
						}

					}
				}
			}

			
			{
				//doublereal dx = 1, dy = 1, dz = 1;
				doublereal dx = xpos[i_maximally_elongated_cell_1 + 1] - xpos[i_maximally_elongated_cell_1];
				doublereal dy = ypos[j_maximally_elongated_cell_1 + 1] - ypos[j_maximally_elongated_cell_1];
				doublereal dz = zpos[k_maximally_elongated_cell_1 + 1] - zpos[k_maximally_elongated_cell_1];
				if ((dx >= dy) && (dy >= dz)) {
					if (dx / dz > ratio_quality) {
						bcont = true;
						SetLength(xpos, inx + 1, inx + 2);
						for (integer l = inx; l >= i_maximally_elongated_cell_1 + 1; l--) xpos[l + 1] = xpos[l];
						//xpos[inx + 1] = xpos[i_maximally_elongated_cell_1] + 0.5*dx;
						xpos[i_maximally_elongated_cell_1 + 1] = xpos[i_maximally_elongated_cell_1] + 0.5*dx;
						
						inx = inx + 1;
						//BubbleEnhSort<doublereal>(xpos,0, inx);
						//Sort_method<doublereal>(xpos,inx);
						icorrect_stat++;
						goto START_LAB2;
					}
				}
				if ((dx >= dz) && (dz >= dy)) {
					if (dx / dy > ratio_quality) {
						bcont = true;
						SetLength(xpos, inx + 1, inx + 2);
						for (integer l = inx; l >= i_maximally_elongated_cell_1 + 1; l--) xpos[l + 1] = xpos[l];
						//xpos[inx + 1] = xpos[i_maximally_elongated_cell_1] + 0.5*dx;
						xpos[i_maximally_elongated_cell_1 + 1] = xpos[i_maximally_elongated_cell_1] + 0.5*dx;
						
						inx = inx + 1;
						//BubbleEnhSort<doublereal>(xpos, 0, inx);
						//Sort_method<doublereal>(xpos, inx);
						icorrect_stat++;
						goto START_LAB2;
					}
				}
				if ((dy >= dx) && (dx >= dz)) {
					if (dy / dz > ratio_quality) {
						bcont = true;
						SetLength(ypos, iny + 1, iny + 2);
						for (integer l = iny; l >= j_maximally_elongated_cell_1 + 1; l--) ypos[l + 1] = ypos[l];
						//ypos[iny + 1] = ypos[j_maximally_elongated_cell_1] + 0.5*dy;
						ypos[j_maximally_elongated_cell_1 + 1] = ypos[j_maximally_elongated_cell_1] + 0.5*dy;
						
						iny = iny + 1;
						//BubbleEnhSort<doublereal>(ypos, 0, iny);
						//Sort_method<doublereal>(ypos,iny);
						icorrect_stat++;
						goto START_LAB2;
					}
				}
				if ((dy >= dz) && (dz >= dx)) {
					if (dy / dx > ratio_quality) {
						bcont = true;
						SetLength(ypos, iny + 1, iny + 2);
						for (integer l = iny; l >= j_maximally_elongated_cell_1 + 1; l--) ypos[l + 1] = ypos[l];
						//ypos[iny + 1] = ypos[j_maximally_elongated_cell_1] + 0.5*dy;
						ypos[j_maximally_elongated_cell_1 + 1] = ypos[j_maximally_elongated_cell_1] + 0.5*dy;
						
						iny = iny + 1;
						//BubbleEnhSort<doublereal>(ypos, 0, iny);
						//Sort_method<doublereal>(ypos,iny);
						icorrect_stat++;
						goto START_LAB2;
					}
				}
				if ((dz >= dy) && (dy >= dx)) {
					if (dz / dx > ratio_quality) {
						bcont = true;
						SetLength(zpos, inz + 1, inz + 2);
						for (integer l = inz; l >= k_maximally_elongated_cell_1 + 1; l--) zpos[l + 1] = zpos[l];
						//zpos[inz + 1] = zpos[k_maximally_elongated_cell_1] + 0.5*dz;
						zpos[k_maximally_elongated_cell_1 + 1] = zpos[k_maximally_elongated_cell_1] + 0.5*dz;
						
						inz = inz + 1;
						//BubbleEnhSort<doublereal>(zpos, 0, inz);
						//Sort_method<doublereal>(zpos,inz);
						icorrect_stat++;
						goto START_LAB2;
					}
				}
				if ((dz >= dx) && (dx >= dy)) {
					if (dz / dy > ratio_quality) {
						bcont = true;
						SetLength(zpos, inz + 1, inz + 2);
						for (integer l = inz; l >= k_maximally_elongated_cell_1 + 1; l--) zpos[l + 1] = zpos[l];
						//zpos[inz + 1] = zpos[k_maximally_elongated_cell_1] + 0.5*dz;
						zpos[k_maximally_elongated_cell_1 + 1] = zpos[k_maximally_elongated_cell_1] + 0.5*dz;
						
						inz = inz + 1;
						//BubbleEnhSort<doublereal>(zpos, 0, inz);
						//Sort_method<doublereal>(zpos,inz);
						icorrect_stat++;
						goto START_LAB2;
					}
				}
			}
			//BubbleEnhSort<doublereal>(xpos, 0, inx);
			//BubbleEnhSort<doublereal>(ypos, 0, iny);
			//BubbleEnhSort<doublereal>(zpos, 0, inz);
			//Sort_method<doublereal>(xpos,inx);
			//Sort_method<doublereal>(ypos,iny);
			//Sort_method<doublereal>(zpos,inz);
			//goto START_LAB;
		}

#if doubleintprecision == 1
		printf("%lld \n", ic4);
#else
		printf("%d \n", ic4);
#endif
		
	}
#if doubleintprecision == 1
	//printf("automatic correct mesh ratio_quality=%e ipass=%lld\n", ratio_quality, icorrect_stat);
	std::cout << "automatic correct mesh ratio_quality=" << ratio_quality << " ipass=" << icorrect_stat << std::endl;
#else
	//printf("automatic correct mesh ratio_quality=%e ipass=%d\n", ratio_quality, icorrect_stat);
	std::cout << "automatic correct mesh ratio_quality=" << ratio_quality << " ipass=" << icorrect_stat << std::endl;
#endif
	

	ratio_start_check = 0.0;

#pragma omp parallel for
	for (int i_1 = 0; i_1 < inx; ++i_1) {
		for (int j_1 = 0; j_1 < iny; ++j_1) {
			for (int k_1 = 0; k_1 < inz; ++k_1) {
				//doublereal dx = 1, dy = 1, dz = 1;
				doublereal dx = xpos[i_1 + 1] - xpos[i_1];
				doublereal dy = ypos[j_1 + 1] - ypos[j_1];
				doublereal dz = zpos[k_1 + 1] - zpos[k_1];
				if ((dx >= dy) && (dy >= dz)) {
					if (dx / dz > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dx / dz;
						}
					}
				}
				if ((dx >= dz) && (dz >= dy)) {
					if (dx / dy > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dx / dy;
						}
					}
				}
				if ((dy >= dx) && (dx >= dz)) {
					if (dy / dz > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dy / dz;
						}
					}
				}
				if ((dy >= dz) && (dz >= dx)) {
					if (dy / dx > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dy / dx;
						}
					}
				}
				if ((dz >= dx) && (dx >= dy)) {
					if (dz / dy > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dz / dy;
						}
					}
				}
				if ((dz >= dy) && (dy >= dx)) {
					if (dz / dx > ratio_start_check) {
#pragma omp critical
						{
							ratio_start_check = dz / dx;
						}
					}
				}
			}
		}
	}
	//printf("mesh post refinement quolite=%e\n", ratio_start_check);
	std::cout << "mesh post refinement quolite=" << ratio_start_check << std::endl;
}

// Работает только для coarse grid генератора сетки.
// Улучшает качество сетки. Если блок разбивается только одной
// клеткой в заданном направлении то включается вставка 
// дополнительной сеточной линии в середину такого блока.
void patch_mesh_refinement_21_11_2019(int &inx,
	int &iny, int &inz, doublereal* &xpos,
	doublereal* &ypos, doublereal* &zpos, int lb,
	BLOCK* &b)
{
	// количество добавляемых клеток
	const int icell_insert = 2;// 1 or 2
	const doublereal epsilon = 1.0e-30;

	std::cout << "apriori refinement: inx=" << inx << " iny=" << iny << " inz=" << inz << "\n";

	// Дополнительное сгущение только в том случае если мы не используем snap_TO
	if (!(((bsnap_TO_global == 1) || (bsnap_TO_global == 3)))) {
		// Ox
		bool bweShouldbecontinue1 = true;
		while (bweShouldbecontinue1) {
			bweShouldbecontinue1 = false;

			for (integer i = 0; i < lb; i++) {
				if (b[i].g.itypegeom == PRISM) {
					doublereal xfound = b[i].g.xS;
					for (integer j = 0; j <= inx; j++) {
						//if (fabs(xpos[j] - xfound) < epsilon) {
						// Массив xpos должен быть отсортирован по возрастанию.
						if (((xfound - xpos[j]) < epsilon)&&((xpos[j] - xfound) < epsilon)) {
							// xfound  ==xS найден
							if ((j < inx) && (fabs(xpos[j + 1] - b[i].g.xE) < epsilon)) 
							{
								if (icell_insert == 1) {
									// Обнаружена всего одна клетка в блоке
									SetLength(xpos, inx + 1, inx + 2);
									// Добавляем в середине блока всего одну клетку.
									xpos[inx + 1] = 0.5*(b[i].g.xS + b[i].g.xE);
									inx = inx + 1;
									bweShouldbecontinue1 = true;
								}
								if (icell_insert == 2) {
									// Обнаружена всего одна клетка в блоке
									SetLength(xpos, inx + 1, inx + 3);
									// Добавляем в середине блока всего одну клетку.
									xpos[inx + 1] = b[i].g.xS + 0.3333*fabs(b[i].g.xE - b[i].g.xS);
									xpos[inx + 2] = b[i].g.xS + 0.6666*fabs(b[i].g.xE - b[i].g.xS);
									inx = inx + 2;
									bweShouldbecontinue1 = true;
								}
								//BubbleEnhSort<doublereal>(xpos, 0, inx);
								Sort_method<doublereal>(xpos, inx);
								
								goto END_WHILE_MARKER_001;
							}
						}
					}
				}
			}

		END_WHILE_MARKER_001:;
			//marker Ox
		}

		bweShouldbecontinue1 = true;
		while (bweShouldbecontinue1) {
			bweShouldbecontinue1 = false;

			for (integer i = 0; i < lb; i++) {
				if (b[i].g.itypegeom == PRISM) {
					doublereal yfound = b[i].g.yS;
					for (integer j = 0; j <= iny; j++) {
						if (fabs(ypos[j] - yfound) < epsilon) {
							// yfound  == yS найден
							if ((j < iny) && (fabs(ypos[j + 1] - b[i].g.yE) < epsilon)) {
								// Обнаружена всего одна клетка в блоке
								if (icell_insert == 1) {
									SetLength(ypos, iny + 1, iny + 2);
									// Добавляем в середине блока всего одну клетку.
									ypos[iny + 1] = 0.5*(b[i].g.yS + b[i].g.yE);
									iny = iny + 1;
									bweShouldbecontinue1 = true;
								}
								if (icell_insert == 2) {
									SetLength(ypos, iny + 1, iny + 3);
									// Добавляем в середине блока всего одну клетку.
									ypos[iny + 1] = b[i].g.yS + 0.3333*fabs(b[i].g.yE - b[i].g.yS);
									ypos[iny + 2] = b[i].g.yS + 0.6666*fabs(b[i].g.yE - b[i].g.yS);
									iny = iny + 2;
									bweShouldbecontinue1 = true;
								}
								//BubbleEnhSort<doublereal>(ypos, 0, iny);
								Sort_method<doublereal>(ypos, iny);								
								goto END_WHILE_MARKER_002;
							}
						}
					}
				}
			}

		END_WHILE_MARKER_002:;
			//marker Oy
		}


		bweShouldbecontinue1 = true;
		while (bweShouldbecontinue1) {
			bweShouldbecontinue1 = false;

			for (integer i = 0; i < lb; i++) {
				if (b[i].g.itypegeom == PRISM) {
					doublereal zfound = b[i].g.zS;
					for (integer j = 0; j <= inz; j++) {
						if (fabs(zpos[j] - zfound) < epsilon) {
							// zfound  == zS найден
							if ((j < inz) && (fabs(zpos[j + 1] - b[i].g.zE) < epsilon)) {
								// Обнаружена всего одна клетка в блоке
								if (icell_insert == 1) {
									SetLength(zpos, inz + 1, inz + 2);
									// Добавляем в середине блока всего одну клетку.
									zpos[inz + 1] = 0.5*(b[i].g.zS + b[i].g.zE);
									inz = inz + 1;
									bweShouldbecontinue1 = true;
								}
								if (icell_insert == 2) {
									SetLength(zpos, inz + 1, inz + 3);
									// Добавляем в середине блока всего одну клетку.
									zpos[inz + 1] = b[i].g.zS + 0.3333*(b[i].g.zE- b[i].g.zS);
									zpos[inz + 2] = b[i].g.zS + 0.6666*(b[i].g.zE - b[i].g.zS);
									inz = inz + 2;
									bweShouldbecontinue1 = true;
								}
								//BubbleEnhSort<doublereal>(zpos, 0, inz);
								Sort_method<doublereal>(zpos, inz);								
								goto END_WHILE_MARKER_003;
							}
						}
					}
				}
			}

		END_WHILE_MARKER_003:;
			//marker Oz
		}
	}

	std::cout << "apostoriory refinement: inx=" << inx <<" iny=" << iny << " inz="<< inz << "\n";
}// patch_mesh_refinement_21_11_2019



// проверяет принадлежит ли контрольный объём
// гидродинамической модели || hollow block.
// Возвращает параметр ib равный номеру блока
// которому принадлежит контрольный объём.
bool in_model_fluid_gap_1(TOCHKA p, integer &ib, BLOCK* b, integer lb) {
	int i = 0, k = 0;
	bool ret = true;// по умолчанию принадлежит модели
					// цикл по всем блокам
	for (i = 0; i<lb; i++) {

		if (b[i].g.itypegeom == PRISM) {
			// Prism
			if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {
				k = i;
			}
		}
		else if (b[i].g.itypegeom == CAD_STL) {
			// CAD_STL
			// 15.11.2020

			int k_loc = -1;
			if (b[i].g.in_CAD_STL_check(p,k_loc,i)) {
				k = i;
			}
		}
		else if (b[i].g.itypegeom == POLYGON) {
			// Мы проверяем принадлежность полигону только в случае если точка p находидся
			// строго внутри прямоугольной призмы окаймляющей полигон.
			if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {
				// определяет принадлежность точки полигону.
				in_polygon(p, b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2, k, i);
			}

		}
		else {
		/*if (b[i].g.itypegeom == 1) {*/

			// Cylinder
			switch (b[i].g.iPlane) {
			case XY_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.z > b[i].g.zC) && (p.z < b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							k = i;
						}
					}
				}
				else {
					if ((p.z > b[i].g.zC) && (p.z < b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) > b[i].g.R_in_cyl) {
								k = i;
							}
						}
					}
				}
				break;
			case XZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.y > b[i].g.yC) && (p.y < b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
						}
					}
				}
				else {
					if ((p.y > b[i].g.yC) && (p.y < b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
							}
						}
					}
				}
				break;
			case YZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.x > b[i].g.xC) && (p.x < b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
						}
					}
				}
				else {
					if ((p.x > b[i].g.xC) && (p.x < b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
							}
						}
					}
				}
				break;
			}
		}
	}
	if ((b[k].itype == PHYSICS_TYPE_IN_BODY::SOLID)) ret = false;
	ib = k;

	return ret;

} // in_model_fluid_gap_1

  // проверяет принадлежит ли контрольный объём
  // гидродинамической модели.
  // Возвращает параметр ib равный номеру блока
  // которому принадлежит контрольный объём.
bool in_model_fluid_gap(TOCHKA p, int &ib, BLOCK* b, int lb) {
	int i = 0, k = 0;
	bool ret = true;// по умолчанию принадлежит модели

	// 27_12_2017.
	// Блоки перечисляются согласно приоритетам.
	// Сначала идут блоки с низким приоритетом начиная с нуля.
	// Блок  большим номером перезаписывает свойства блока с меньшим номером.
	// Поэтому если сканировать блоки с конца (начиная с самомого большого номера lb и далее
	// в сторону уменьшения номеров, то первый найденный блок как раз и будет нужным и можно
	// досрочно прервать цикл for с помощью break и сэкономить существенную долю времени.
	
	
	// цикл по всем блокам
	for (i = lb - 1; i >= 0; i--) {

		if (b[i].g.itypegeom == PRISM) {
			// Prism
			if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {
				k = i;
				// Только нашли и сразу закончили проверку.
				goto OUTOF_IN_MODEL_FLOW_1;
			}
		}
		if (b[i].g.itypegeom == CAD_STL) {
			// CAD_STL
			// 15.11.2020

			int k_loc = -1;
			if (b[i].g.in_CAD_STL_check(p, k_loc, i)) {
				k = i;
				// Только нашли и сразу закончили проверку.
				goto OUTOF_IN_MODEL_FLOW_1;
			}
		}
		if (b[i].g.itypegeom == POLYGON) {

			if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {

				bool found = false;
				// определяет принадлежность точки полигону.
				found = in_polygon(p, b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2, k, i);
				if (found) {
					// Только нашли и сразу закончили проверку.
					goto OUTOF_IN_MODEL_FLOW_1;
				}

			}

		}

		if (b[i].g.itypegeom == CYLINDER) {
			// Cylinder
			switch (b[i].g.iPlane) {
			case XY_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.z > b[i].g.zC) && (p.z < b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							k = i;
							// Только нашли и сразу закончили проверку.
							goto OUTOF_IN_MODEL_FLOW_1;
						}
					}
				}
				else {
					if ((p.z > b[i].g.zC) && (p.z < b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) > b[i].g.R_in_cyl) {
								k = i;
								// Только нашли и сразу закончили проверку.
								goto OUTOF_IN_MODEL_FLOW_1;
							}
						}
					}
				}
				break;
			case XZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.y > b[i].g.yC) && (p.y < b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Только нашли и сразу закончили проверку.
							goto OUTOF_IN_MODEL_FLOW_1;
						}
					}
				}
				else {
					if ((p.y > b[i].g.yC) && (p.y < b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Только нашли и сразу закончили проверку.
								goto OUTOF_IN_MODEL_FLOW_1;
							}
						}
					}
				}
				break;
			case YZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.x > b[i].g.xC) && (p.x < b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Только нашли и сразу закончили проверку.
							goto OUTOF_IN_MODEL_FLOW_1;
						}
					}
				}
				else {
					if ((p.x > b[i].g.xC) && (p.x < b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Только нашли и сразу закончили проверку.
								goto OUTOF_IN_MODEL_FLOW_1;
							}
						}
					}
				}
				break;
			}
		}
	}

OUTOF_IN_MODEL_FLOW_1:

	if ((b[k].itype == PHYSICS_TYPE_IN_BODY::SOLID)) ret = false;
	ib = k;

	return ret;

} // in_model_fluid_gap

// улучшает качество аппроксимации сеткой цилиндров вводя дополнительно для
// каждого цилиндра четыре сеточных линии на расстоянии шага сетки от вшешнего
// контура основания цилиндра. Это требуется для сильно неравномерных сеток.
// 25.04.2018 Асемблесы.
void calc_minimum_fluid_gap3(int &inumboundaryx, doublereal* &rxboundary,
	int &inumboundaryy, doublereal* &ryboundary,
	int &inumboundaryz, doublereal* &rzboundary,
	int lb, int ls, int lw, BLOCK* b, SOURCE* &s, WALL* w,
	int lu, UNION* &my_union, integer &iunion_id_p1)
{
	doublereal dm_start = 1.0 / dcount_diametr_cylinder;
	int i;

	bool *bactive_cyl = new bool[lb];
	for (i = 0; i < lb; i++) {
		bactive_cyl[i] = true; // Поначалу все активны.
	}

	// доработка АМПЛИТРОН.
	//const doublereal dopusk = 1.0e-20; // Для определения совпадения центров.

	for (i = 1; i < lb; i++) {
		if (bactive_cyl[i]) {
			if ((b[i].g.itypegeom == CYLINDER) && (b[i].g.iPlane == XY_PLANE)) {
				integer b_big = i; // Совпадающий цилиндр с самым большим внешним радиусом.
				bool bincomming = false;
				for (integer j = 1; j < lb; j++) {
					if (i != j) {
						if ((b[j].g.itypegeom == CYLINDER) && (b[j].g.iPlane == XY_PLANE)) {
							if ((fabs(b[i].g.xC - b[j].g.xC) < /*dopusk*/ shorter_length_for_simplificationX(b[j].g.xC, b, lb, w, lw, s, ls)) &&
								((fabs(b[i].g.yC - b[j].g.yC) < /*dopusk*/	shorter_length_for_simplificationY(b[j].g.yC, b, lb, w, lw, s, ls)))) {
								// Центры и плоскости совпали.
								if (b[j].g.R_out_cyl > b[b_big].g.R_out_cyl) {
									b_big = j;
								}
								bincomming = true;
							}
						}
					}
				}
				if (bincomming) {
					if (b_big != i) {
						bactive_cyl[i] = false;// не добавляем линий для цилиндров совпадающими центрами, сетка и так очень подробна.
					}
					for (integer j = 1; j < lb; j++) {
						if (i != j) {
							if ((b[j].g.itypegeom == CYLINDER) && (b[j].g.iPlane == XY_PLANE)) {
								if ((fabs(b[i].g.xC - b[j].g.xC) < /*dopusk*/ shorter_length_for_simplificationX(b[j].g.xC, b, lb, w, lw, s, ls)) &&
									((fabs(b[i].g.yC - b[j].g.yC) < /*dopusk*/ shorter_length_for_simplificationY(b[j].g.yC, b, lb, w, lw, s, ls)))) {
									// Центры и плоскости совпали.
									if (b_big != j) {
										bactive_cyl[j] = false;// не добавляем линий для цилиндров совпадающими центрами, сетка и так очень подробна.
									}									
								}
							}
						}
					}
				}
			}
			if ((b[i].g.itypegeom == CYLINDER) && (b[i].g.iPlane == XZ_PLANE)) {
				integer b_big = i; // Совпадающий цилиндр с самым большим внешним радиусом.
				bool bincomming = false;
				for (integer j = 1; j < lb; j++) {
					if (i != j) {
						if ((b[j].g.itypegeom == CYLINDER) && (b[j].g.iPlane == XZ_PLANE)) {
							if ((fabs(b[i].g.xC - b[j].g.xC) < /*dopusk*/ shorter_length_for_simplificationX(b[j].g.xC, b, lb, w, lw, s, ls)) &&
								((fabs(b[i].g.zC - b[j].g.zC) < /*dopusk*/ shorter_length_for_simplificationZ(b[j].g.zC, b, lb, w, lw, s, ls)))) {
								// Центры и плоскости совпали.
								if (b[j].g.R_out_cyl > b[b_big].g.R_out_cyl) {
									b_big = j;
								}
								bincomming = true;
							}
						}
					}
				}
				if (bincomming) {
					if (b_big != i) {
						bactive_cyl[i] = false;// не добавляем линий для цилиндров совпадающими центрами, сетка и так очень подробна.
					}
					for (integer j = 1; j < lb; j++) {
						if (i != j) {
							if ((b[j].g.itypegeom == CYLINDER) && (b[j].g.iPlane == XZ_PLANE)) {
								if ((fabs(b[i].g.xC - b[j].g.xC) < /*dopusk*/shorter_length_for_simplificationX(b[j].g.xC, b, lb, w, lw, s, ls)) &&
									((fabs(b[i].g.zC - b[j].g.zC) < /*dopusk*/ shorter_length_for_simplificationZ(b[j].g.zC, b, lb, w, lw, s, ls)))) {
									// Центры и плоскости совпали.
									if (b_big != j) {
										bactive_cyl[j] = false;// не добавляем линий для цилиндров совпадающими центрами, сетка и так очень подробна.
									}
								}
							}
						}
					}
				}
			}
			if ((b[i].g.itypegeom == CYLINDER) && (b[i].g.iPlane == YZ_PLANE)) {
				integer b_big = i; // Совпадающий цилиндр с самым большим внешним радиусом.
				bool bincomming = false;
				for (integer j = 1; j < lb; j++) {
					if (i != j) {
						if ((b[j].g.itypegeom == CYLINDER) && (b[j].g.iPlane == YZ_PLANE)) {
							if ((fabs(b[i].g.yC - b[j].g.yC) < /*dopusk*/ shorter_length_for_simplificationY(b[j].g.yC, b, lb, w, lw, s, ls)) &&
								((fabs(b[i].g.zC - b[j].g.zC) < /*dopusk*/ shorter_length_for_simplificationZ(b[j].g.zC, b, lb, w, lw, s, ls)))) {
								// Центры и плоскости совпали.
								if (b[j].g.R_out_cyl > b[b_big].g.R_out_cyl) {
									b_big = j;
								}
								bincomming = true;
							}
						}
					}
				}
				if (bincomming) {
					if (b_big != i) {
						bactive_cyl[i] = false;// не добавляем линий для цилиндров совпадающими центрами, сетка и так очень подробна.
					}
					for (integer j = 1; j < lb; j++) {
						if (i != j) {
							if ((b[j].g.itypegeom == CYLINDER) && (b[j].g.iPlane == YZ_PLANE)) {
								if ((fabs(b[i].g.yC - b[j].g.yC) < /*dopusk*/ shorter_length_for_simplificationY(b[j].g.yC, b, lb, w, lw, s, ls)) &&
									((fabs(b[i].g.zC - b[j].g.zC) < /*dopusk*/ shorter_length_for_simplificationZ(b[j].g.zC, b, lb, w, lw, s, ls)))) {
									// Центры и плоскости совпали.
									if (b_big != j) {
										bactive_cyl[j] = false;// не добавляем линий для цилиндров совпадающими центрами, сетка и так очень подробна.
									}
								}
							}
						}
					}
				}
			}
		}
	}


	for (i = 1; i < lb; i++) {
		if (b[i].iunion_id == iunion_id_p1) {
			if ((b[i].g.itypegeom != POLYGON)&&(b[i].g.itypegeom != CAD_STL)) {

				if (bcylinder_meshing && (b[i].g.itypegeom == CYLINDER) && (bactive_cyl[i])) {
					// Cylinder

					switch (b[i].g.iPlane) {
					case XY_PLANE: case XZ_PLANE:

						// 24.01.2018
						// Небольшой запас в бока для того чтобы форма окружности цилиндра не искажалась.
						doublereal x_1 = b[i].g.xC - b[i].g.R_out_cyl - dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((x_1>=b[0].g.xS)&&(x_1<=b[0].g.xE)) {
							addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
						}
						x_1 = b[i].g.xC + b[i].g.R_out_cyl + dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
							addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
						}
						break;
					}

				}
			}
		}
	}


	for (i = 1; i < lb; i++) {
		if (b[i].iunion_id == iunion_id_p1) {
			if ((b[i].g.itypegeom != POLYGON) && (b[i].g.itypegeom != CAD_STL)) {

				if (bcylinder_meshing && (b[i].g.itypegeom == CYLINDER) && (bactive_cyl[i])) {
					// Cylinder
					switch (b[i].g.iPlane) {
					case XY_PLANE: case YZ_PLANE:

						// 24.01.2018
						// Небольшой запас в бока для того чтобы форма окружности цилиндра не искажалась.
						doublereal y_1 = b[i].g.yC - b[i].g.R_out_cyl - dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
							addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
						}
						y_1 = b[i].g.yC + b[i].g.R_out_cyl + dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
							addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
						}
						break;
					}

				}
			}
		}
	}

	for (i = 1; i < lb; i++) {
		if (b[i].iunion_id == iunion_id_p1) {
			if ((b[i].g.itypegeom != POLYGON) && (b[i].g.itypegeom != CAD_STL)) {


				if (bcylinder_meshing && (b[i].g.itypegeom == CYLINDER && (bactive_cyl[i]))) {
					// Cylinder
					switch (b[i].g.iPlane) {
					case XZ_PLANE: case YZ_PLANE:

						// 24.01.2018
						// Небольшой запас в бока для того чтобы форма окружности цилиндра не искажалась.
						doublereal z_1 = b[i].g.zC - b[i].g.R_out_cyl - dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
							addboundary(rzboundary, inumboundaryz, z_1, XY_PLANE, b, lb, w, lw, s, ls);
						}
						z_1 = b[i].g.zC + b[i].g.R_out_cyl + dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
							addboundary(rzboundary, inumboundaryz, z_1, XY_PLANE, b, lb, w, lw, s, ls);
						}
						break;
					}

				}
			}
		}
	}

	delete[] bactive_cyl;

	if (iunion_id_p1 == 0) {
		// Добавляем границы юнионов.
		for (i = 0; i < lu; i++) {
			if ((fabs(my_union[i].xE - my_union[i].xS) >= fabs(my_union[i].zE - my_union[i].zS)) && (fabs(my_union[i].xE - my_union[i].xS) >= fabs(my_union[i].yE - my_union[i].yS))) {
				// больше всех
				for (integer i_65 = 0; i_65 < 17; i_65++) {
					doublereal x_1 = my_union[i].xS + i_65 * 0.0625 * fabs(my_union[i].xE - my_union[i].xS);
					if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
						addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}
			else if ((fabs(my_union[i].xE - my_union[i].xS) <= fabs(my_union[i].zE - my_union[i].zS)) && (fabs(my_union[i].xE - my_union[i].xS) <= fabs(my_union[i].yE - my_union[i].yS))) {
				// меньше всех
				for (integer i_65 = 0; i_65 < 6; i_65++) {
					doublereal x_1 = my_union[i].xS + i_65 * 0.2 * fabs(my_union[i].xE - my_union[i].xS);
					if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
						addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}
			else {
				// Серединка.
				for (integer i_65 = 0; i_65 < 11; i_65++) {
					doublereal x_1 = my_union[i].xS + i_65 * 0.1 * fabs(my_union[i].xE - my_union[i].xS);
					if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
						addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}
			addboundary(rxboundary, inumboundaryx, my_union[i].xS,YZ_PLANE, b, lb, w, lw, s, ls);
			addboundary(rxboundary, inumboundaryx, my_union[i].xE,YZ_PLANE, b, lb, w, lw, s, ls);
		}
	

	
		// Добавляем границы юнионов.
		for (i = 0; i < lu; i++) {
			if ((fabs(my_union[i].yE - my_union[i].yS) >= fabs(my_union[i].zE - my_union[i].zS)) && (fabs(my_union[i].yE - my_union[i].yS) >= fabs(my_union[i].xE - my_union[i].xS))) {
				// больше всех
				for (integer i_65 = 0; i_65 < 17; i_65++) {
					doublereal y_1 = my_union[i].yS + i_65 * 0.0625 * fabs(my_union[i].yE - my_union[i].yS);
					if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
						addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}
			else if ((fabs(my_union[i].yE - my_union[i].yS) <= fabs(my_union[i].zE - my_union[i].zS)) && (fabs(my_union[i].yE - my_union[i].yS) <= fabs(my_union[i].xE - my_union[i].xS))) {
				// меньше всех
				for (integer i_65 = 0; i_65 < 6; i_65++) {
					doublereal y_1 = my_union[i].yS + i_65 * 0.2 * fabs(my_union[i].yE - my_union[i].yS);
					if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
						addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}
			else {
				// Серединка.
				for (integer i_65 = 0; i_65 < 11; i_65++) {
					doublereal y_1 = my_union[i].yS + i_65 * 0.1 * fabs(my_union[i].yE - my_union[i].yS);
					if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
						addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}
			doublereal y_1 = my_union[i].yS;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
			y_1 = my_union[i].yE;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	

	
		// Добавляем границы юнионов.
		for (i = 0; i < lu; i++) {
			if ((fabs(my_union[i].zE - my_union[i].zS) >= fabs(my_union[i].yE - my_union[i].yS)) && (fabs(my_union[i].zE - my_union[i].zS) >= fabs(my_union[i].xE - my_union[i].xS))) {
			// больше всех
				for (integer i_65 = 0; i_65 < 17; i_65++) {
					doublereal z_1 = my_union[i].zS + i_65 * 0.0625 * fabs(my_union[i].zE - my_union[i].zS);
					if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
						addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}
			else if ((fabs(my_union[i].zE - my_union[i].zS) <= fabs(my_union[i].yE - my_union[i].yS)) && (fabs(my_union[i].zE - my_union[i].zS) <= fabs(my_union[i].xE - my_union[i].xS))) {
				// меньше всех
				for (integer i_65 = 0; i_65 < 6; i_65++) {
					doublereal z_1 = my_union[i].zS + i_65 * 0.2 * fabs(my_union[i].zE - my_union[i].zS);
					if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
						addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}
			else {
				// Серединка.
				for (integer i_65 = 0; i_65 < 11; i_65++) {
					doublereal z_1 = my_union[i].zS + i_65 * 0.1 * fabs(my_union[i].zE - my_union[i].zS);
					if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
						addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}
			doublereal z_1 = my_union[i].zS;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
			z_1 = my_union[i].zE;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
		}
	} // (iunion_id_p1 == 0)

}// calc_minimum_fluid_gap3

// Вычисляет minimum fluid gap. 
// ЧАСТЬ 1.
// 12.03.2017
// 16.08.2017 Polygon
// 25.04.2018 Асемблесы.
// 15.11.2020 CAD_STL
void calc_minimum_fluid_gap1(int &inumboundaryx, doublereal* &rxboundary,
	int &inumboundaryy, doublereal* &ryboundary,
	int &inumboundaryz, doublereal* &rzboundary,
	int lb, int ls, int lw, BLOCK* b, SOURCE* &s, WALL* w, 
	int lu, UNION* &my_union, integer &iunion_id_p1)
{

	// Влючить для обнаружения щелей в модели.
	bool bdiagnostic_analysys = false;

	doublereal dm_start = 1.0/ dcount_diametr_cylinder;
	integer i;

	//***************** x **************************

	inumboundaryx = 1;
	rxboundary = new doublereal[inumboundaryx + 1]; // число границ по оси x

	if (iunion_id_p1 == 0) {
		// добавить границы кабинета
		rxboundary[0] = b[0].g.xS; // начало области
		rxboundary[inumboundaryx] = b[0].g.xE; // конец области
	}
	else {
		// добавить границы юниона.
		rxboundary[0] = my_union[iunion_id_p1-1].xS; // начало области
		rxboundary[inumboundaryx] = my_union[iunion_id_p1 - 1].xE; // конец области
	}

	// блоки
	for (i = 1; i<lb; i++) {
		if (b[i].iunion_id == iunion_id_p1) {
			if ((b[i].g.itypegeom != POLYGON) && (b[i].g.itypegeom != CAD_STL)) {
				doublereal x_1 = b[i].g.xS;
				if ((x_1>=b[0].g.xS) && (x_1<=b[0].g.xE)) {
					addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
				}
				x_1 = b[i].g.xE;
				if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
					addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
				}
				if (bcylinder_meshing && (b[i].g.itypegeom == 1)) {
					// Cylinder

					switch (b[i].g.iPlane) {
					case XY_PLANE: case XZ_PLANE:
						for (doublereal dm = dm_start; dm < 0.98; dm = dm + dm_start) {
							x_1 = b[i].g.xC - b[i].g.R_out_cyl + dm * 2.0 * b[i].g.R_out_cyl;
							if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
								addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
							}
						}
						// 24.01.2018
						// Небольшой запас в бока для того чтобы форма окружности цилиндра не искажалась.
						x_1 = b[i].g.xC - b[i].g.R_out_cyl - 0.5 * dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
							addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
						}
						x_1 = b[i].g.xC + b[i].g.R_out_cyl + 0.5 * dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
							addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
						}
						break;
					}

				}
				else if (0&&(b[i].arr_Sc[0] > 0.0) && (b[i].g.itypegeom == 0)) {
					// 27.10.2018
					// Измельчаем сетку на блоках в которых выделяется тепловая мощность. 
					for (doublereal dm = dm_start; dm < 0.98; dm = dm + dm_start) {
						doublereal x_1 = b[i].g.xS - fabs(b[i].g.xE - b[i].g.xS) + dm * 3.0 * fabs(b[i].g.xE - b[i].g.xS);
						if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
							addboundary(rxboundary, inumboundaryx, x_1, YZ_PLANE, b, lb, w, lw, s, ls);
						}
					}
				}
			}
			else if (b[i].g.itypegeom == POLYGON) {
				//Polygon
				doublereal x_1;
				switch (b[i].g.iPlane_obj2) {
				case XY_PLANE:
					for (integer i_65 = 0; i_65 < b[i].g.nsizei; i_65++) {
						x_1 = b[i].g.xi[i_65];
						if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
							addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
						}
					}
					break;
				case XZ_PLANE:
					for (integer i_65 = 0; i_65 < b[i].g.nsizei; i_65++) {
						x_1 = b[i].g.xi[i_65];
						if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
							addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
						}
					}
					break;
				case YZ_PLANE:
					x_1 = b[i].g.xi[0];
					if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
						addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
					}
					x_1 = b[i].g.xi[0] + b[i].g.hi[0];
					if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
						addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
					}
					break;
				}
			}
			else {
				// CAD STL

				gCAD_STL* tmp = b[i].g.root_CAD_STL;

				while (tmp != nullptr) {

					doublereal x_1=tmp->pa.x;

					addboundary(rxboundary, inumboundaryx, x_1, YZ_PLANE, b, lb, w, lw, s, ls);

					x_1 = tmp->pb.x;

					addboundary(rxboundary, inumboundaryx, x_1, YZ_PLANE, b, lb, w, lw, s, ls);

					x_1 = tmp->pc.x;

					addboundary(rxboundary, inumboundaryx, x_1, YZ_PLANE, b, lb, w, lw, s, ls);

					tmp = tmp->next;

				}

			}
		}
	}

	if (iunion_id_p1 == 0) {
		// Добавляем границы юнионов.
		for (i = 0; i < lu; i++) {
			doublereal x_1 = my_union[i].xS;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
			x_1 = my_union[i].xE;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// источники
	for (i = 0; i<ls; i++) {
		if (s[i].iunion_id == iunion_id_p1) {
			doublereal x_1 = s[i].g.xS;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
			x_1 = s[i].g.xE;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// стенки
	for (i = 0; i<lw; i++) {
		if (w[i].iunion_id == iunion_id_p1) {
			doublereal x_1 = w[i].g.xS;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
			x_1 = w[i].g.xE;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// упорядочивание по возрастанию
	//BubbleEnhSort(rxboundary, inumboundaryx);
	Sort_method(rxboundary, inumboundaryx);

	if (bdiagnostic_analysys) {
		printf("diagnostic X:");
		for (integer i73 = 0; i73 < inumboundaryx; i73++) {
			//printf("%lld %e %e\n", i73, rxboundary[i73], rxboundary[i73 + 1] - rxboundary[i73]);
			std::cout << i73 << " " << rxboundary[i73] << " " << rxboundary[i73 + 1] - rxboundary[i73] << std::endl;
		}
		//printf("%lld %e\n", inumboundaryx, rxboundary[inumboundaryx]);
		std::cout << inumboundaryx << " " << rxboundary[inumboundaryx] << std::endl;
		system("PAUSE");

	}

	//*********** y  **************************

	ryboundary = new doublereal[inumboundaryy + 1]; // число границ по оси y
	

	if (iunion_id_p1 == 0) {
		// добавить границы кабинета
		ryboundary[0] = b[0].g.yS; // начало области 
		ryboundary[inumboundaryy] = b[0].g.yE; // конец области 
	}
	else {
		// добавить границы юниона.
		ryboundary[0] = my_union[iunion_id_p1 - 1].yS; // начало области
		ryboundary[inumboundaryy] = my_union[iunion_id_p1 - 1].yE; // конец области
	}

	// блоки
	for (i = 1; i<lb; i++) {
		if (b[i].iunion_id == iunion_id_p1) {
			if ((b[i].g.itypegeom != POLYGON)&&(b[i].g.itypegeom != CAD_STL)) {
				doublereal y_1 = b[i].g.yS;
				if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
					addboundary(ryboundary, inumboundaryy, y_1, XZ_PLANE, b, lb, w, lw, s, ls);
				}
				y_1 = b[i].g.yE;
				if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
					addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
				}
				if (bcylinder_meshing && (b[i].g.itypegeom == 1)) {
					// Cylinder
					switch (b[i].g.iPlane) {
					case XY_PLANE: case YZ_PLANE:
						for (doublereal dm = dm_start; dm < 0.98; dm = dm + dm_start) {
							y_1 = b[i].g.yC - b[i].g.R_out_cyl + dm * 2.0 * b[i].g.R_out_cyl;
							if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
								addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
							}
						}
						// 24.01.2018
						// Небольшой запас в бока для того чтобы форма окружности цилиндра не искажалась.
						y_1 = b[i].g.yC - b[i].g.R_out_cyl - 0.5 * dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
							addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
						}
						y_1 = b[i].g.yC + b[i].g.R_out_cyl + 0.5 * dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
							addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
						}
						break;
					}

				}
				else if (0 && (b[i].arr_Sc[0] > 0.0) && (b[i].g.itypegeom == 0)) {
					// 27.10.2018
					// Измельчаем сетку на блоках в которых выделяется тепловая мощность. 
					for (doublereal dm = dm_start; dm < 0.98; dm = dm + dm_start) {
						y_1 = b[i].g.yS - fabs(b[i].g.yE - b[i].g.yS) + dm * 3.0 * fabs(b[i].g.yE - b[i].g.yS);
						if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
							addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
						}
					}
				}
			}
			else if (b[i].g.itypegeom == POLYGON) {
				//Polygon
				doublereal y_1;
				switch (b[i].g.iPlane_obj2) {
				case XY_PLANE:
					for (integer i_65 = 0; i_65 < b[i].g.nsizei; i_65++) {
						y_1 = b[i].g.yi[i_65];
						if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
							addboundary(ryboundary, inumboundaryy, y_1, XZ_PLANE, b, lb, w, lw, s, ls);
						}
					}
					break;
				case XZ_PLANE:
					y_1 = b[i].g.yi[0];
					if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
						addboundary(ryboundary, inumboundaryy, y_1, XZ_PLANE, b, lb, w, lw, s, ls);
					}
					y_1 = b[i].g.yi[0] + b[i].g.hi[0];
					if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
						addboundary(ryboundary, inumboundaryy, y_1, XZ_PLANE, b, lb, w, lw, s, ls);
					}
					break;
				case YZ_PLANE:
					for (integer i_65 = 0; i_65 < b[i].g.nsizei; i_65++) {
						y_1 = b[i].g.yi[i_65];
						if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
							addboundary(ryboundary, inumboundaryy, y_1, XZ_PLANE, b, lb, w, lw, s, ls);
						}
					}
					break;
				}
			}
			else {
				// CAD STL

				gCAD_STL* tmp = b[i].g.root_CAD_STL;

				while (tmp != nullptr) {

					doublereal y_1 = tmp->pa.y;

					addboundary(ryboundary, inumboundaryy, y_1, XZ_PLANE, b, lb, w, lw, s, ls);

					y_1 = tmp->pb.y;

					addboundary(ryboundary, inumboundaryy, y_1, XZ_PLANE, b, lb, w, lw, s, ls);

					y_1 = tmp->pc.y;

					addboundary(ryboundary, inumboundaryy, y_1, XZ_PLANE, b, lb, w, lw, s, ls);

					tmp = tmp->next;

				}
			}
		}
	}

	if (iunion_id_p1 == 0) {
		// Добавляем границы юнионов.
		for (i = 0; i < lu; i++) {
			doublereal y_1 = my_union[i].yS;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
			y_1 = my_union[i].yE;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// источники
	for (i = 0; i<ls; i++) {
		if (s[i].iunion_id == iunion_id_p1) {
			doublereal y_1 = s[i].g.yS;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
			y_1 = s[i].g.yE;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// стенки
	for (i = 0; i<lw; i++) {
		if (w[i].iunion_id == iunion_id_p1) {
			doublereal y_1 = w[i].g.yS;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
			y_1 = w[i].g.yE;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// упорядочивание по возрастанию
	//BubbleEnhSort(ryboundary, inumboundaryy);
	Sort_method(ryboundary, inumboundaryy);

	if (bdiagnostic_analysys) {
		printf("diagnostic Y:");
		for (integer i73 = 0; i73 < inumboundaryy; i73++) {
			//printf("%lld %e %e\n", i73, ryboundary[i73], ryboundary[i73 + 1] - ryboundary[i73]);
			std::cout << i73 << " " << ryboundary[i73] << " " << ryboundary[i73 + 1] - ryboundary[i73] << std::endl;
		}
		//printf("%lld %e\n", inumboundaryy, ryboundary[inumboundaryy]);
		std::cout << inumboundaryy << " " << ryboundary[inumboundaryy] << std::endl;
		system("PAUSE");
	}

	//*******************  z **************************************

	rzboundary = new doublereal[1 + inumboundaryz]; // число границ по оси y

	
	if (iunion_id_p1 == 0) {
		// добавить границы кабинета
		rzboundary[0] = b[0].g.zS; // начало области 
		rzboundary[inumboundaryz] = b[0].g.zE; // конец области 
	}
	else {
		// добавить границы юниона.
		rzboundary[0] = my_union[iunion_id_p1 - 1].zS; // начало области
		rzboundary[inumboundaryz] = my_union[iunion_id_p1 - 1].zE; // конец области
	}


	// блоки
	for (i = 1; i<lb; i++) {
		if (b[i].iunion_id == iunion_id_p1) {
			if ((b[i].g.itypegeom != POLYGON) && (b[i].g.itypegeom != CAD_STL)) {
				doublereal z_1 = b[i].g.zS;
				if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
					addboundary(rzboundary, inumboundaryz,z_1,XY_PLANE, b, lb, w, lw, s, ls);
				}
				z_1 = b[i].g.zE;
				if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
					addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
				}

				if (bcylinder_meshing && (b[i].g.itypegeom == 1)) {
					// Cylinder
					switch (b[i].g.iPlane) {
					case XZ_PLANE: case YZ_PLANE:
						for (doublereal dm = dm_start; dm < 0.98; dm = dm + dm_start) {
							doublereal z_1 = b[i].g.zC - b[i].g.R_out_cyl + dm * 2.0 * b[i].g.R_out_cyl;
							if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
								addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
							}
						}
						// 24.01.2018
						// Небольшой запас в бока для того чтобы форма окружности цилиндра не искажалась.
						doublereal z_1 = b[i].g.zC - b[i].g.R_out_cyl - 0.5 * dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
							addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
						}
						z_1 = b[i].g.zC + b[i].g.R_out_cyl + 0.5 * dm_start * 2.0 * b[i].g.R_out_cyl;
						if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
							addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
						}
						break;
					}

				}
				else if (0 && (b[i].arr_Sc[0] > 0.0) && (b[i].g.itypegeom == 0)) {
					// 27.10.2018
					// Измельчаем сетку на блоках в которых выделяется тепловая мощность. 

					for (doublereal dm = dm_start; dm < 0.98; dm = dm + dm_start) {
						doublereal z_1 = b[i].g.zS - fabs(b[i].g.zE - b[i].g.zS) + dm * 3.0 * fabs(b[i].g.zE - b[i].g.zS);
						if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
							addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
						}
					}
				}
			}
			else if (b[i].g.itypegeom == POLYGON) {
				//Polygon
				doublereal z_1;
				switch (b[i].g.iPlane_obj2) {
				case XY_PLANE:
					z_1 = b[i].g.zi[0];
					if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
						addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
					}
					z_1 = b[i].g.zi[0] + b[i].g.hi[0];
					if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
						addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
					}
					break;
				case XZ_PLANE:
					for (integer i_65 = 0; i_65 < b[i].g.nsizei; i_65++) {
						 z_1 = b[i].g.zi[i_65];
						if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
							addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
						}
					}
					break;
				case YZ_PLANE:
					for (integer i_65 = 0; i_65 < b[i].g.nsizei; i_65++) {
						 z_1 = b[i].g.zi[i_65];
						if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
							addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
						}
					}
					break;
				}
			}
			else {
				// CAD STL

				gCAD_STL* tmp = b[i].g.root_CAD_STL;

				while (tmp != nullptr) {

					doublereal z_1 = tmp->pa.z;

					addboundary(rzboundary, inumboundaryz, z_1, XY_PLANE, b, lb, w, lw, s, ls);

					z_1 = tmp->pb.z;

					addboundary(rzboundary, inumboundaryz, z_1, XY_PLANE, b, lb, w, lw, s, ls);

					z_1 = tmp->pc.z;

					addboundary(rzboundary, inumboundaryz, z_1, XY_PLANE, b, lb, w, lw, s, ls);

					tmp = tmp->next;

				}
			}
		}
	}

	if (iunion_id_p1 == 0) {
		// Добавляем границы юнионов.
		for (i = 0; i < lu; i++) {
			doublereal z_1 = my_union[i].zS;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz,z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
			z_1 = my_union[i].zE;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// источники
	for (i = 0; i<ls; i++) {
		if (s[i].iunion_id == iunion_id_p1) {
			doublereal z_1 = s[i].g.zS;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
			z_1 = s[i].g.zE;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// стенки
	for (i = 0; i<lw; i++) {
		if (w[i].iunion_id == iunion_id_p1) {
			doublereal z_1 = w[i].g.zS;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz,z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
			z_1 = w[i].g.zE;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// упорядочивание по возрастанию
	//BubbleEnhSort(rzboundary, inumboundaryz);
	Sort_method(rzboundary, inumboundaryz);

	if (bdiagnostic_analysys) {
		printf("diagnostic Z:");
		for (int i73 = 0; i73 < inumboundaryz; i73++) {
			//printf("%d %e %e\n", i73, rzboundary[i73], rzboundary[i73 + 1] - rzboundary[i73]);
			std::cout << i73 << " " << rzboundary[i73] << " " << rzboundary[i73 + 1] - rzboundary[i73] << std::endl;
		}
		//printf("%lld %e\n", inumboundaryz, rzboundary[inumboundaryz]);
		std::cout << inumboundaryz << " " << rzboundary[inumboundaryz] << std::endl;
		system("PAUSE");
	}

}


// Вычисляет minimum fluid gap. 
// ЧАСТЬ 2.
// 12.03.2017
void calc_minimum_fluid_gap2_1(integer &inumboundaryx, doublereal* &rxboundary,
	integer &inumboundaryy, doublereal* &ryboundary,
	integer &inumboundaryz, doublereal* &rzboundary,
	doublereal &minimum_fluid_gap_x,
	doublereal &minimum_fluid_gap_y,
	doublereal &minimum_fluid_gap_z,
	int lb, int ls, int lw, BLOCK* b, SOURCE* &s, WALL* w)
{
	// minimum fluid gap in Ox direction

	bool bincomming = false; // true начался жидкостный блок.
	doublereal startpos = -1.0e36;
	for (integer i7 = 0; i7 < inumboundaryy; i7++) {
		for (integer i8 = 0; i8 < inumboundaryz; i8++) {
			doublereal yc = 0.5*(ryboundary[i7] + ryboundary[i7 + 1]);
			doublereal zc = 0.5*(rzboundary[i8] + rzboundary[i8 + 1]);
			for (integer i9 = 0; i9 < inumboundaryx; i9++) {
				doublereal xc = 0.5*(rxboundary[i9] + rxboundary[i9 + 1]);
				int ib;
				TOCHKA p;
				p.x = xc;
				p.y = yc;
				p.z = zc;
				if ((bincomming) && (!in_model_fluid_gap(p, ib, b, lb))) {
					bincomming = false;
					if (fabs(rxboundary[i9] - startpos) < minimum_fluid_gap_x) {
						minimum_fluid_gap_x = fabs(rxboundary[i9] - startpos);
					}
				}
				if ((!bincomming) && (in_model_fluid_gap(p, ib, b, lb))) {
					startpos = rxboundary[i9];
					bincomming = true;
				}
			}
		}
	}


	// minimum fluid gap in Oy direction

	bincomming = false; // true начался жидкостный блок.
	startpos = -1.0e36;
	for (integer i7 = 0; i7 < inumboundaryx; i7++) {
		for (integer i8 = 0; i8 < inumboundaryz; i8++) {
			doublereal xc = 0.5*(rxboundary[i7] + rxboundary[i7 + 1]);
			doublereal zc = 0.5*(rzboundary[i8] + rzboundary[i8 + 1]);
			for (integer i9 = 0; i9 < inumboundaryy; i9++) {
				doublereal yc = 0.5*(ryboundary[i9] + ryboundary[i9 + 1]);
				int ib;
				TOCHKA p;
				p.x = xc;
				p.y = yc;
				p.z = zc;
				if ((bincomming) && (!in_model_fluid_gap(p, ib, b, lb))) {
					bincomming = false;
					if (fabs(ryboundary[i9] - startpos) < minimum_fluid_gap_y) {
						minimum_fluid_gap_y = fabs(ryboundary[i9] - startpos);
					}
				}
				if ((!bincomming) && (in_model_fluid_gap(p, ib, b, lb))) {
					startpos = ryboundary[i9];
					bincomming = true;
				}
			}
		}
	}

	// minimum fluid gap in Oz direction

	bincomming = false; // true начался жидкостный блок.
	startpos = -1.0e36;
	for (integer i7 = 0; i7 < inumboundaryx; i7++) {
		for (integer i8 = 0; i8 < inumboundaryy; i8++) {
			doublereal xc = 0.5*(rxboundary[i7] + rxboundary[i7 + 1]);
			doublereal yc = 0.5*(ryboundary[i8] + ryboundary[i8 + 1]);
			for (integer i9 = 0; i9 < inumboundaryz; i9++) {
				doublereal zc = 0.5*(rzboundary[i9] + rzboundary[i9 + 1]);
				int ib;
				TOCHKA p;
				p.x = xc;
				p.y = yc;
				p.z = zc;
				if ((bincomming) && (!in_model_fluid_gap(p, ib, b, lb))) {
					bincomming = false;
					if (fabs(rzboundary[i9] - startpos) < minimum_fluid_gap_z) {
						minimum_fluid_gap_z = fabs(rzboundary[i9] - startpos);
					}
				}
				if ((!bincomming) && (in_model_fluid_gap(p, ib, b, lb))) {
					startpos = rzboundary[i9];
					bincomming = true;
				}
			}
		}
	}
} //calc_minimum_fluid_gap2_1





// Вычисляет minimum fluid gap. 
// ЧАСТЬ 2.
// 12.03.2017
// 25.04.2018 готов для работы с АСЕМБЛЕСАМИ
void calc_minimum_fluid_gap2(int &inumboundaryx, doublereal* &rxboundary,
	int &inumboundaryy, doublereal* &ryboundary,
	int &inumboundaryz, doublereal* &rzboundary,
	doublereal &minimum_fluid_gap_x,
	doublereal &minimum_fluid_gap_y,
	doublereal &minimum_fluid_gap_z,
	integer lb, integer ls, integer lw, BLOCK* b, SOURCE* &s, WALL* w, 
	integer lu, UNION* &my_union, integer &iunion_id_p1)
{
	
	
	// preprocessing
	integer* ib_marker = new integer[inumboundaryx*inumboundaryy*inumboundaryz];
	integer* ib_marker_yxz = new integer[inumboundaryx*inumboundaryy*inumboundaryz];
	integer* ib_marker_zxy = new integer[inumboundaryx*inumboundaryy*inumboundaryz];
	/*
	for (integer i9 = 0; i9 < inumboundaryx; i9++) {
		for (integer i7 = 0; i7 < inumboundaryy; i7++) {
			for (integer i8 = 0; i8 < inumboundaryz; i8++) {
				doublereal yc = 0.5*(ryboundary[i7] + ryboundary[i7 + 1]);
				doublereal zc = 0.5*(rzboundary[i8] + rzboundary[i8 + 1]);
				doublereal xc = 0.5*(rxboundary[i9] + rxboundary[i9 + 1]);
				int ib;
				TOCHKA p;
				p.x = xc;
				p.y = yc;
				p.z = zc;
				in_model_fluid_gap(p, ib, b, lb);
				// x + y*dimx+z*dimx*dimy
				ib_marker[i9+ inumboundaryx*i7+ inumboundaryx*inumboundaryy*i8] = ib;
			}
		}
	}
	*/

	
	Block_indexes* block_indexes = nullptr;
	block_indexes = new Block_indexes[lb];
	//if (block_indexes == nullptr) {
		//printf("error in allocation memory for block_indexes in enumerate_volume_improved.\n");
		//system("pause");
		//exit(1);
	//}
	

	// 08.04.2018
	for (integer i = 0; i < lb; i++) {
		// инициализация, на случай если блоки не будут распознаны.
		block_indexes[i].iL = -1;
		block_indexes[i].iR = -2;
		block_indexes[i].jL = -1;
		block_indexes[i].jR = -2;
		block_indexes[i].kL = -1;
		block_indexes[i].kR = -2;
	}

	/*
	for (integer i = 0; i < lb; i++) {
		doublereal x4 = b[i].g.xS;
		doublereal distmax;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryx; j++) {
			if (fabs(rxboundary[j] - x4) < distmax) {
				block_indexes[i].iL = j;
				distmax = fabs(rxboundary[j] - x4);
			}
		}
		x4 = b[i].g.xE;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryx; j++) {
			if (fabs(rxboundary[j] - x4) < distmax) {
				block_indexes[i].iR = j;
				distmax = fabs(rxboundary[j] - x4);
			}
		}
		x4 = b[i].g.yS;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryy; j++) {
			if (fabs(ryboundary[j] - x4) < distmax) {
				block_indexes[i].jL = j;
				distmax = fabs(ryboundary[j] - x4);
			}
		}
		x4 = b[i].g.yE;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryy; j++) {
			if (fabs(ryboundary[j] - x4) < distmax) {
				block_indexes[i].jR = j;
				distmax = fabs(ryboundary[j] - x4);
			}
		}
		x4 = b[i].g.zS;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryz; j++) {
			if (fabs(rzboundary[j] - x4) < distmax) {
				block_indexes[i].kL = j;
				distmax = fabs(rzboundary[j] - x4);
			}
		}
		x4 = b[i].g.zE;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryz; j++) {
			if (fabs(rzboundary[j] - x4) < distmax) {
				block_indexes[i].kR = j;
				distmax = fabs(rzboundary[j] - x4);
			}
		}
	}
	*/
	
	for (integer i = 0; i < lb; i++) {
		//if (b[i].iunion_id == iunion_id_p1) {
		{	
		    doublereal x4 = b[i].g.xS;
			for (integer j = 0; j <= inumboundaryx; j++) {
				if (fabs(rxboundary[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].iL = j;
					break;
				}
			}
			x4 = b[i].g.xE;
			for (integer j = 0; j <= inumboundaryx; j++) {
				if (fabs(rxboundary[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].iR = j;
					break;
				}
			}
			x4 = b[i].g.yS;
			for (integer j = 0; j <= inumboundaryy; j++) {
				if (fabs(ryboundary[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].jL = j;
					break;
				}
			}
			x4 = b[i].g.yE;
			for (integer j = 0; j <= inumboundaryy; j++) {
				if (fabs(ryboundary[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].jR = j;
					break;
				}
			}
			x4 = b[i].g.zS;
			for (integer j = 0; j <= inumboundaryz; j++) {
				if (fabs(rzboundary[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].kL = j;
					break;
				}
			}
			x4 = b[i].g.zE;
			for (integer  j = 0; j <= inumboundaryz; j++) {
				if (fabs(rzboundary[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].kR = j;
					break;
				}
			}
		}
	}

	
	integer ismarker = 0;
	if (iunion_id_p1 > 0) {
		for (integer m7 = 0; m7 < lb; m7++) {

			if (block_indexes[m7].iL < 0) {
				block_indexes[m7].iL = 0;
				ismarker++;
			}
			if (block_indexes[m7].jL < 0) {
				block_indexes[m7].jL = 0;
				ismarker++;
			}
			if (block_indexes[m7].kL < 0) {
				block_indexes[m7].kL = 0;
				ismarker++;
			}
			if (block_indexes[m7].iR < 0) {
				block_indexes[m7].iR = inumboundaryx;
				ismarker++;
			}
			if (block_indexes[m7].jR < 0) {
				block_indexes[m7].jR = inumboundaryy;
				ismarker++;
			}
			if (block_indexes[m7].kR < 0) {
				block_indexes[m7].kR = inumboundaryz;
				ismarker++;
			}
		}
	}
	//printf("ismarker =%d\n",ismarker);
	//getchar();
	

	// Количество проходов существенно сократилось и в итоге это приводит к существенному
	// увеличению быстродействия.
	integer m7;
	integer ib_stub = -1;
	// Мы найдем самый большой по размеру Hollow block, иначе будет просто кабинет.
	ib_stub = 0;
	doublereal vol_stub = -1.0;
	
	for (integer i = 0; i < lb; i++) {
		//if (b[i].iunion_id == iunion_id_p1) {
		{
			if (b[i].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
				if (fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS) > vol_stub) {
					ib_stub = i;
					vol_stub = fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS);
				}
			}
		}
	}

	
#pragma omp parallel for
	for (integer k1 = 0; k1 < inumboundaryz; k1++) 
	{
		integer iPk1 = inumboundaryx * inumboundaryy * k1;
		for (integer j1 = 0; j1 < inumboundaryy; j1++) 
		{
			integer iPj1 = iPk1 + inumboundaryx * j1;
			for (integer i1 = 0; i1 < inumboundaryx; i1++)
			{
				integer iPi1 = i1 + iPj1;
				ib_marker[iPi1 ] = ib_stub; //-1
			}
		}
	}
	for (m7 = 0; m7 < lb; m7++) {
		if (b[m7].iunion_id == iunion_id_p1) {

#pragma omp parallel for
			for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++)
			{
				integer iPk1 = inumboundaryx * inumboundaryy * k1;
				for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++)
				{
					integer iPj1 = iPk1 + inumboundaryx * j1;
					for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++) {
						integer iPi1 = i1 + iPj1;
						ib_marker[iPi1] = m7;
					}
				}
			}
		}
	}

	//yxz
#pragma omp parallel for
	for (integer j1 = 0; j1 < inumboundaryy; j1++)
	{
		integer iPj1 = inumboundaryx * inumboundaryz * j1;
		for (integer i1 = 0; i1 < inumboundaryx; i1++)
		{
			integer iPi1 = iPj1 + inumboundaryz * i1;
			for (integer k1 = 0; k1 < inumboundaryz; k1++)
			{
				integer iPk1 = k1 + iPi1;
				ib_marker_yxz[iPk1] = ib_stub; //-1
			}
		}
	}
	for (m7 = 0; m7 < lb; m7++) {
		if (b[m7].iunion_id == iunion_id_p1) {

#pragma omp parallel for
			for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++)
			{
				integer iPj1 = inumboundaryx * inumboundaryz * j1;
				for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++)
				{
					integer iPi1 = iPj1 + inumboundaryz * i1;
					for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++) {
						integer iPk1 = k1 + iPi1;
						ib_marker_yxz[iPk1] = m7;
					}
				}
			}
		}
	}

	// zxy
#pragma omp parallel for
	for (integer k1 = 0; k1 < inumboundaryz; k1++)
	{
		integer iPk1 = inumboundaryx * inumboundaryy * k1;
		for (integer i1 = 0; i1 < inumboundaryx; i1++)
		{
			integer iPi1 = iPk1 + inumboundaryy * i1;
			for (integer j1 = 0; j1 < inumboundaryy; j1++)
			{
				integer iPj1 = j1 + iPi1;
				ib_marker_zxy[iPj1] = ib_stub; //-1
			}
		}
	}
	for (m7 = 0; m7 < lb; m7++) {
		if (b[m7].iunion_id == iunion_id_p1) {

#pragma omp parallel for
			for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++)
			{
				integer iPk1 = inumboundaryx * inumboundaryy * k1;
				for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++)
				{
					integer iPi1 = iPk1 + inumboundaryy * i1;
					for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++) {
						integer iPj1 = j1 + iPi1;
						ib_marker_zxy[iPj1] = m7;
					}
				}
			}
		}
	}

	delete[] block_indexes;
	//printf("identifire blocks number 80 procent.\n");


	// minimum fluid gap in Ox direction

	bool bincomming = false; // true начался жидкостный блок.
	doublereal startpos = -1.0e36;
	for (integer i8 = 0; i8 < inumboundaryz; i8++) {
		integer i8_ = inumboundaryx * inumboundaryy * i8;
	     for (integer i7 = 0; i7 < inumboundaryy; i7++) {
			 integer i7_ = inumboundaryx * i7 + i8_;
			//doublereal yc = 0.5*(ryboundary[i7] + ryboundary[i7 + 1]);
			//doublereal zc = 0.5*(rzboundary[i8] + rzboundary[i8 + 1]);
			for (integer i9 = 0; i9 < inumboundaryx; i9++) {
				//doublereal xc = 0.5*(rxboundary[i9] + rxboundary[i9 + 1]);
				//int ib;
				//TOCHKA p;
				//p.x = xc;
				//p.y = yc;
				//p.z = zc;
				bool b_this_is_SOLID_block = false;
				if (b[ib_marker[i9 +  i7_]].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
					b_this_is_SOLID_block = true;
				}
				if ((bincomming) && ((b_this_is_SOLID_block))) {
					bincomming = false;
					if (fabs(rxboundary[i9] - startpos) < minimum_fluid_gap_x) {
						minimum_fluid_gap_x = fabs(rxboundary[i9] - startpos);
					}
				}
				if ((!bincomming) && (!(b_this_is_SOLID_block))) {
					startpos = rxboundary[i9];
					bincomming = true;
				}
			}
		}
	}

	delete[] ib_marker;
	// minimum fluid gap in Oy direction

	bincomming = false; // true начался жидкостный блок.
	startpos = -1.0e36;
	for (integer i8 = 0; i8 < inumboundaryz; i8++) {
		integer i8_ = inumboundaryx * inumboundaryy * i8;
		for (integer i7 = 0; i7 < inumboundaryx; i7++) {		
			//doublereal xc = 0.5*(rxboundary[i7] + rxboundary[i7 + 1]);
			//doublereal zc = 0.5*(rzboundary[i8] + rzboundary[i8 + 1]);
			integer i7_ = i7*inumboundaryy + i8_;
			for (integer i9 = 0; i9 < inumboundaryy; i9++) {
				//doublereal yc = 0.5*(ryboundary[i9] + ryboundary[i9 + 1]);
				//int ib;
				//TOCHKA p;
				//p.x = xc;
				//p.y = yc;
				//p.z = zc;
				bool b_this_is_SOLID_block = false;
				if ((b[ib_marker_zxy[i9+ i7_]].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
					b_this_is_SOLID_block = true;
				}
				if ((bincomming) && (b_this_is_SOLID_block)) {
					bincomming = false;
					if (fabs(ryboundary[i9] - startpos) < minimum_fluid_gap_y) {
						minimum_fluid_gap_y = fabs(ryboundary[i9] - startpos);
					}
				}
				if ((!bincomming) && (!(b_this_is_SOLID_block))) {
					startpos = ryboundary[i9];
					bincomming = true;
				}
			}
		}
	}

	delete[] ib_marker_zxy;
	// minimum fluid gap in Oz direction

	bincomming = false; // true начался жидкостный блок.
	startpos = -1.0e36;
	integer imultxz = inumboundaryx * inumboundaryz;
	for (integer i8 = 0; i8 < inumboundaryy; i8++) {
		integer i8_ = imultxz*i8;
		for (integer i7 = 0; i7 < inumboundaryx; i7++) {
			integer i7_ = i7*inumboundaryz + i8_;
			//doublereal xc = 0.5*(rxboundary[i7] + rxboundary[i7 + 1]);
			//doublereal yc = 0.5*(ryboundary[i8] + ryboundary[i8 + 1]);
			for (integer i9 = 0; i9 < inumboundaryz; i9++) {
				//doublereal zc = 0.5*(rzboundary[i9] + rzboundary[i9 + 1]);
				//int ib;
				//TOCHKA p;
				//p.x = xc;
				//p.y = yc;
				//p.z = zc;
				bool b_this_is_SOLID_block = false;
				if (b[ib_marker_yxz[i7_  +  i9]].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
					b_this_is_SOLID_block = true;
				}
				if ((bincomming) && (b_this_is_SOLID_block)) {
					bincomming = false;
					if (fabs(rzboundary[i9] - startpos) < minimum_fluid_gap_z) {
						minimum_fluid_gap_z = fabs(rzboundary[i9] - startpos);
					}
				}
				if ((!bincomming) && (!(b_this_is_SOLID_block))) {
					startpos = rzboundary[i9];
					bincomming = true;
				}
			}
		}
	}

	delete[] ib_marker_yxz;
	
}



// 12.03.2017
// реализация snap to
// уменьшающает размерность сеточной модели на 33%.
// 25.04.2018 Готова к АССЕМБЛЕСАМ
void snap_to_moving(bool* &source_indexpopadaniqnagranYZ,
	bool* &source_indexpopadaniqnagranXY,
	bool* &source_indexpopadaniqnagranXZ,
	doublereal* &rxboundary, doublereal* &ryboundary, doublereal* &rzboundary,
	int &inumboundaryx, int &inumboundaryy, int &inumboundaryz,
	doublereal &minimum_fluid_gap_x, doublereal &minimum_fluid_gap_y, doublereal &minimum_fluid_gap_z,
	int lb, int ls, int lw, BLOCK* &b, SOURCE* &s, WALL* &w,
	int lu, UNION* &my_union, integer &iunion_id_p1)
{

	doublereal dm_start = 1.0/ dcount_diametr_cylinder;

	// 0.	none
    // 1.	Snap to grid
    // 2.	Snap to grid ALICE
    // 3.	Snap to grid ++
	bool bsnap_TO = false;  // snap to grid (избавляет от щелей в сложных моделях).
	if ((bsnap_TO_global == 1) || (bsnap_TO_global == 3)) {
		bsnap_TO = true;
	}
	doublereal snap_to_multiplyer = 0.3;// константа лежит в диапазоне (0..1) и регулирует какой зазор игнорировать.

	bool brepeat = true;

	integer i;

RESTARTX_SNAPTO:

	if (source_indexpopadaniqnagranYZ != nullptr) {
		delete[] source_indexpopadaniqnagranYZ;
		source_indexpopadaniqnagranYZ = nullptr;
	}
	if (rxboundary != nullptr) {
		delete[] rxboundary;
		rxboundary = nullptr;
	}

	inumboundaryx = 1;
	rxboundary = new doublereal[inumboundaryx + 1]; // число границ по оси x

	

	if (iunion_id_p1 == 0) {
		// добавить границы кабинета
		rxboundary[0] = b[0].g.xS; // начало области
		rxboundary[inumboundaryx] = b[0].g.xE; // конец области
	}
	else {
		// добавить границы юниона.
		rxboundary[0] = my_union[iunion_id_p1 - 1].xS; // начало области
		rxboundary[inumboundaryx] = my_union[iunion_id_p1 - 1].xE; // конец области
	}

	// блоки
	for (i = 1; i<lb; i++) {
		if (b[i].iunion_id == iunion_id_p1) {
			doublereal x_1 = b[i].g.xS;
			if ((x_1>=b[0].g.xS) && (x_1<=b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
			x_1 = b[i].g.xE;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}

			if (bcylinder_meshing && (b[i].g.itypegeom == 1)) {
				// Cylinder

				switch (b[i].g.iPlane) {
				case XY_PLANE: case XZ_PLANE:
					for (doublereal dm = dm_start; dm < 0.98; dm = dm + dm_start) {
						//addboundary(rxboundary, inumboundaryx, dm*(b[i].g.xS+ b[i].g.xE));
						x_1 = (b[i].g.xC - b[i].g.R_out_cyl + dm * 2.0 * b[i].g.R_out_cyl);
						if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
							addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
						}
					}
					break;
				}
			}
		}
	}

	// Надо выяснить попал ли источник на границу двух блоков,
	// если попал то это плохо и его придётся сдвигать.
	source_indexpopadaniqnagranYZ = new bool[ls];
	for (i = 0; i < ls; i++) {
		source_indexpopadaniqnagranYZ[i] = false;
	}
	for (i = 0; i < ls; i++) {
		if (s[i].iunion_id == iunion_id_p1) {
			if (s[i].iPlane == YZ_PLANE) {
				for (integer i1 = 0; i1 <= inumboundaryx; i1++) {
					s[i].g.xS = s[i].g.xE;
					if (fabs(s[i].g.xS - rxboundary[i1]) < 1.0e-36) {
						source_indexpopadaniqnagranYZ[i] = true;
					}
				}
			}
		}
	}

	if (iunion_id_p1 == 0) {
		// Добавляем границы юнионов.
		for (i = 0; i < lu; i++) {
			doublereal x_1 = my_union[i].xS;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
			x_1 = my_union[i].xE;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// источники
	for (i = 0; i<ls; i++) {
		if (s[i].iunion_id == iunion_id_p1) {
			doublereal x_1 = s[i].g.xS;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
			x_1 = s[i].g.xE;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// стенки
	for (i = 0; i<lw; i++) {
		if (w[i].iunion_id == iunion_id_p1) {
			doublereal x_1 = w[i].g.xS;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
			x_1 = w[i].g.xE;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(rxboundary, inumboundaryx, x_1,YZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// упорядочивание по возрастанию
	//BubbleEnhSort(rxboundary, inumboundaryx);
	Sort_method(rxboundary, inumboundaryx);

	// snap to grid
	if (bsnap_TO) {
		doublereal rmindisteps = 1.0e30;

		for (i = 0; i <inumboundaryx; i++) {
			if (fabs(rxboundary[i + 1] - rxboundary[i]) < rmindisteps) {
				rmindisteps = fabs(rxboundary[i + 1] - rxboundary[i]);
			}
		}


		rmindisteps *= 0.25;

		doublereal rmindist = 1.0e30;

		for (i = 0; i < lb; i++) {
			if (b[i].iunion_id == iunion_id_p1) {
				if (bcylinder_meshing && (b[i].g.itypegeom == CYLINDER)) {
					// Cylinder

					switch (b[i].g.iPlane) {
					case XY_PLANE: case XZ_PLANE:
						if (fabs(2.0*b[i].g.R_out_cyl) / dcount_diametr_cylinder < rmindist) {
							rmindist = fabs(2.0*b[i].g.R_out_cyl) / dcount_diametr_cylinder;
						}
						break;
					}
				}
				else if (b[i].g.itypegeom == POLYGON) {
					doublereal dist = 1.0e30;
					integer i_7;
					switch (b[i].g.iPlane_obj2) {
					case XY_PLANE:
						for (i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
							dist = sqrt((b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) * (b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) + (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) * (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						{
							dist = sqrt((b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) * (b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) + (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) * (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						break;
					case XZ_PLANE:
						for (i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
							dist = sqrt((b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) * (b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) + (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]) * (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						{
							dist = sqrt((b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) * (b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) + (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]) * (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						break;
					case YZ_PLANE:
						if (fabs(b[i].g.hi[0]) < rmindist) {
							rmindist = fabs(b[i].g.hi[0]);
						}
						break;
					}
				}
				else if (b[i].g.itypegeom == CAD_STL) {

					doublereal dist47 = b[i].g.min_size_edge();
					if (dist47 < rmindist) {
						rmindist = dist47;
					}
				}
				else if ((b[i].g.itypegeom == PRISM) && (fabs(b[i].g.xE - b[i].g.xS) < rmindist)) {
					rmindist = fabs(b[i].g.xE - b[i].g.xS);
				}
			}
		}
		for (i = 0; i < ls; i++) {
			if (s[i].iunion_id == iunion_id_p1) {
				if ((s[i].iPlane == XY_PLANE) || (s[i].iPlane == XZ_PLANE)) {
					if (fabs(s[i].g.xE - s[i].g.xS) < rmindist) {
						rmindist = fabs(s[i].g.xE - s[i].g.xS);
					}
				}
			}
		}
		for (i = 0; i < lw; i++) {
			if (w[i].iunion_id == iunion_id_p1) {
				if ((w[i].iPlane == XY_PLANE) || (w[i].iPlane == XZ_PLANE)) {
					if (fabs(w[i].g.xE - w[i].g.xS) < rmindist) {
						rmindist = fabs(w[i].g.xE - w[i].g.xS);
					}
				}
			}
		}
		if (iunion_id_p1 == 0) {
			// Добавляем границы юнионов.
			for (i = 0; i < lu; i++) {
				if (fabs(my_union[i].xE - my_union[i].xS) < rmindist) {
					rmindist = fabs(my_union[i].xE - my_union[i].xS);
				}
			}
		}

		// 28.02.2017.
		// Здесь нужно учитывать minimum fluid gap
		// Смысл в том что нельзя сливать блоки если между ними больше либо равно minimum fluid gap*snap_to_multiplyer.
		// Причем minimum fluid gap определяется не только по FLUID доменам. Если мы имеем дело с задачей чистой теплопередачи
		// то minimum fluid gap должен быть определён также по hollow блокам.
		// Вышесказанное приводит к тому что Snap to не работает корректно во всех случаях.

		if (minimum_fluid_gap_x < rmindist) rmindist = minimum_fluid_gap_x;


		rmindist *= snap_to_multiplyer;
		doublereal movetopos;
		doublereal changepos;
		bool bmove = false;
		bool bmove2 = false;
		for (i = 0; i < inumboundaryx; i++) {
			if (bmove) bmove2 = true;
			bmove = false;
			if (fabs(rxboundary[i + 1] - rxboundary[i]) < rmindist) {
				if (i > 0) {
					if (i + 2 <= inumboundaryx) {

						if (fabs(rxboundary[i - 1] - rxboundary[i]) <= fabs(rxboundary[i + 2] - rxboundary[i + 1])) {
							if (fabs(rxboundary[i - 1] - rxboundary[i]) < fabs(rxboundary[i + 1] - rxboundary[i])) {
								// такого быть скорее всего не может.
								movetopos = rxboundary[i - 1];
								changepos = rxboundary[i];
							}
							else {
								// наиболее вероятно.
								movetopos = rxboundary[i + 1];
								changepos = rxboundary[i];
							}
							bmove = true;
						}
						else {
							if (fabs(rxboundary[i + 2] - rxboundary[i + 1]) < fabs(rxboundary[i + 1] - rxboundary[i])) {
								// такого быть скорее всего не может.
								movetopos = rxboundary[i + 2];
								changepos = rxboundary[i + 1];
							}
							else {
								// наиболее вероятно.
								movetopos = rxboundary[i];
								changepos = rxboundary[i + 1];
							}

							bmove = true;
						}

					}
					else {
						movetopos = rxboundary[inumboundaryx];
						changepos = rxboundary[inumboundaryx - 1];
						bmove = true;
					}
				}
				else {
					movetopos = rxboundary[0];
					changepos = rxboundary[1];
					bmove = true;
				}
			}
			if (bmove) {
				//printf(" X changepos=%e to movetopos=%e \n", changepos, movetopos);
				std::cout << " X changepos=" << changepos << " to movetopos=" << movetopos << std::endl;
				for (integer i_2 = 0; i_2 < lb; i_2++) {
					if (b[i_2].iunion_id == iunion_id_p1) {
						if (fabs(b[i_2].g.xE - changepos) < rmindisteps) {
							if (b[i_2].g.itypegeom == PRISM) {
								b[i_2].g.xE = movetopos;
							}
						}
						if (fabs(b[i_2].g.xS - changepos) < rmindisteps) {
							if (b[i_2].g.itypegeom == PRISM) {
								b[i_2].g.xS = movetopos;
							}
						}
					}
				}
				for (integer i_2 = 0; i_2 < lw; i_2++) {
					if (w[i_2].iunion_id == iunion_id_p1) {
						if (fabs(w[i_2].g.xE - changepos) < rmindisteps) {
							w[i_2].g.xE = movetopos;
						}
						if (fabs(w[i_2].g.xS - changepos) < rmindisteps) {
							w[i_2].g.xS = movetopos;
						}
					}
				}
				for (integer i_2 = 0; i_2 < ls; i_2++) {
					if (s[i_2].iunion_id == iunion_id_p1) {
						if (fabs(s[i_2].g.xE - changepos) < rmindisteps) {
							s[i_2].g.xE = movetopos;
						}
						if (fabs(s[i_2].g.xS - changepos) < rmindisteps) {
							s[i_2].g.xS = movetopos;
						}
					}
				}
			}
		}

		if (bmove2) {
			if (brepeat) {
				brepeat = false;
				goto RESTARTX_SNAPTO;
			}

		}

	}


	brepeat = true;

	// По оси Oy
	inumboundaryy = 1;


RESTARTY_SNAPTO:

	if (source_indexpopadaniqnagranXZ != nullptr) {
		delete[] source_indexpopadaniqnagranXZ;
		source_indexpopadaniqnagranXZ = nullptr;
	}
	if (ryboundary != nullptr) {
		delete[] ryboundary;
		ryboundary = nullptr;
	}


	inumboundaryy = 1;
	ryboundary = new doublereal[inumboundaryy + 1]; // число границ по оси y
	

	if (iunion_id_p1 == 0) {
		// добавить границы кабинета
		ryboundary[0] = b[0].g.yS; // начало области 
		ryboundary[inumboundaryy] = b[0].g.yE; // конец области 
	}
	else {
		// добавить границы юниона.
		ryboundary[0] = my_union[iunion_id_p1 - 1].yS; // начало области
		ryboundary[inumboundaryy] = my_union[iunion_id_p1 - 1].yE; // конец области
	}

										   // блоки
	for (i = 1; i<lb; i++) {
		if (b[i].iunion_id == iunion_id_p1) {
			doublereal y_1 = b[i].g.yS;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
			y_1 = b[i].g.yE;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
			if (bcylinder_meshing && (b[i].g.itypegeom == CYLINDER)) {
				// Cylinder
				switch (b[i].g.iPlane) {
				case XY_PLANE: case YZ_PLANE:
					for (doublereal dm = dm_start; dm < 0.98; dm = dm + dm_start) {
						y_1 = (b[i].g.yC - b[i].g.R_out_cyl + dm * 2.0 * b[i].g.R_out_cyl);
						if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
							addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
						}
					}
					break;
				}
			}
		}
	}

	// Надо выяснить попал ли источник на границу двух блоков,
	// если попал то это плохо и его придётся сдвигать.
	source_indexpopadaniqnagranXZ = new bool[ls];
	for (i = 0; i < ls; i++) {
		source_indexpopadaniqnagranXZ[i] = false;
	}
	for (i = 0; i < ls; i++) {
		if (s[i].iunion_id == iunion_id_p1) {
			if (s[i].iPlane == XZ_PLANE) {
				for (integer i1 = 0; i1 <= inumboundaryy; i1++) {
					s[i].g.yS = s[i].g.yE;
					if (fabs(s[i].g.yS - ryboundary[i1]) < 1.0e-36) {
						source_indexpopadaniqnagranXZ[i] = true;
					}
				}
			}
		}
	}

	if (iunion_id_p1 == 0) {
		// Добавляем границы юнионов.
		for (i = 0; i < lu; i++) {
			doublereal y_1 = my_union[i].yS;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
			y_1 = my_union[i].yE;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// источники
	for (i = 0; i<ls; i++) {
		if (s[i].iunion_id == iunion_id_p1) {
			doublereal y_1 = s[i].g.yS;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
			y_1 = s[i].g.yE;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// стенки
	for (i = 0; i<lw; i++) {
		if (w[i].iunion_id == iunion_id_p1) {
			doublereal y_1 = w[i].g.yS;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
			y_1 = w[i].g.yE;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(ryboundary, inumboundaryy, y_1,XZ_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// упорядочивание по возрастанию
	//BubbleEnhSort(ryboundary, inumboundaryy);
	Sort_method(ryboundary, inumboundaryy);

	// snap to grid
	if (bsnap_TO) {
		doublereal rmindist = 1.0e30;

		for (i = 0; i < lb; i++) {
			if (b[i].iunion_id == iunion_id_p1) {
				if (bcylinder_meshing && (b[i].g.itypegeom == CYLINDER)) {
					// Cylinder
					switch (b[i].g.iPlane) {
					case XY_PLANE: case YZ_PLANE:
						if (fabs(2.0*b[i].g.R_out_cyl) / dcount_diametr_cylinder < rmindist) {
							rmindist = fabs(2.0*b[i].g.R_out_cyl) / dcount_diametr_cylinder;
						}
						break;
					}
				}
				else if (b[i].g.itypegeom == POLYGON) {
					doublereal dist = 1.0e30;
					integer i_7;
					switch (b[i].g.iPlane_obj2) {
					case XY_PLANE:
						for (i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
							dist = sqrt((b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) * (b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) + (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) * (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						{
							dist = sqrt((b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) * (b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) + (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) * (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						break;
					case YZ_PLANE:
						for (i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
							dist = sqrt((b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) * (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) + (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]) * (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						{
							dist = sqrt((b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) * (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) + (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]) * (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						break;
					case XZ_PLANE:
						if (fabs(b[i].g.hi[0]) < rmindist) {
							rmindist = fabs(b[i].g.hi[0]);
						}
						break;
					}
				}
				else if (b[i].g.itypegeom == CAD_STL) {

					doublereal dist47 = b[i].g.min_size_edge();
					if (dist47 < rmindist) {
						rmindist = dist47;
					}
				}
				else if ((b[i].g.itypegeom == PRISM) && (fabs(b[i].g.yE - b[i].g.yS) < rmindist)) {
					rmindist = fabs(b[i].g.yE - b[i].g.yS);
				}
			}
		}

		for (i = 0; i < ls; i++) {
			if (s[i].iunion_id == iunion_id_p1) {
				if ((s[i].iPlane == XY_PLANE) || (s[i].iPlane == YZ_PLANE)) {
					if (fabs(s[i].g.yE - s[i].g.yS) < rmindist) {
						rmindist = fabs(s[i].g.yE - s[i].g.yS);
					}
				}
			}
		}
		for (i = 0; i < lw; i++) {
			if (w[i].iunion_id == iunion_id_p1) {
				if ((w[i].iPlane == XY_PLANE) || (w[i].iPlane == YZ_PLANE)) {
					if (fabs(w[i].g.yE - w[i].g.yS) < rmindist) {
						rmindist = fabs(w[i].g.yE - w[i].g.yS);
					}
				}
			}
		}
		if (iunion_id_p1 == 0) {
			// Добавляем границы юнионов.
			for (i = 0; i < lu; i++) {
				if (fabs(my_union[i].yE - my_union[i].yS) < rmindist) {
					rmindist = fabs(my_union[i].yE - my_union[i].yS);
				}
			}
		}

		if (minimum_fluid_gap_y < rmindist) rmindist = minimum_fluid_gap_y;

		rmindist *= snap_to_multiplyer;
		doublereal movetopos;
		doublereal changepos;
		bool bmove = false;
		bool bmove2 = false;
		for (i = 0; i < inumboundaryy; i++) {
			if (bmove) bmove2 = true;
			bmove = false;
			if (fabs(ryboundary[i + 1] - ryboundary[i]) < rmindist) {
				if (i > 0) {
					if (i + 2 <= inumboundaryy) {
						if (fabs(ryboundary[i - 1] - ryboundary[i]) <= fabs(ryboundary[i + 2] - ryboundary[i + 1])) {
							if (fabs(ryboundary[i - 1] - ryboundary[i]) < fabs(ryboundary[i + 1] - ryboundary[i])) {
								// такого быть скорее всего не может.
								movetopos = ryboundary[i - 1];
								changepos = ryboundary[i];
							}
							else {
								// наиболее вероятно.
								movetopos = ryboundary[i + 1];
								changepos = ryboundary[i];
							}
							bmove = true;
						}
						else {
							if (fabs(ryboundary[i + 2] - ryboundary[i + 1]) < fabs(ryboundary[i + 1] - ryboundary[i])) {
								// такого быть скорее всего не может.
								movetopos = ryboundary[i + 2];
								changepos = ryboundary[i + 1];
							}
							else {
								// наиболее вероятно.
								movetopos = ryboundary[i];
								changepos = ryboundary[i + 1];
							}

							bmove = true;
						}



					}
					else {
						movetopos = ryboundary[inumboundaryy];
						changepos = ryboundary[inumboundaryy - 1];
						bmove = true;
					}
				}
				else {
					movetopos = ryboundary[0];
					changepos = ryboundary[1];
					bmove = true;
				}
			}
			if (bmove) {
				//printf("Y changepos=%e to movetopos=%e \n", changepos, movetopos);
				std::cout << " Y changepos=" << changepos << " to movetopos=" << movetopos << std::endl;
				for (integer i_2 = 0; i_2 < lb; i_2++) {
					if (b[i_2].iunion_id == iunion_id_p1) {
						if (fabs(b[i_2].g.yE - changepos) < 1.0e-33) {
							if (b[i_2].g.itypegeom == PRISM) {
								b[i_2].g.yE = movetopos;
							}
						}
						if (fabs(b[i_2].g.yS - changepos) < 1.0e-33) {
							if (b[i_2].g.itypegeom == PRISM) {
								b[i_2].g.yS = movetopos;
							}
						}
					}
				}
				for (integer i_2 = 0; i_2 < lw; i_2++) {
					if (w[i_2].iunion_id == iunion_id_p1) {
						if (fabs(w[i_2].g.yE - changepos) < 1.0e-33) {
							w[i_2].g.yE = movetopos;
						}
						if (fabs(w[i_2].g.yS - changepos) < 1.0e-33) {
							w[i_2].g.yS = movetopos;
						}
					}
				}
				for (integer i_2 = 0; i_2 < ls; i_2++) {
					if (s[i_2].iunion_id == iunion_id_p1) {
						if (fabs(s[i_2].g.yE - changepos) < 1.0e-33) {
							s[i_2].g.yE = movetopos;
						}
						if (fabs(s[i_2].g.yS - changepos) < 1.0e-33) {
							s[i_2].g.yS = movetopos;
						}
					}
				}
			}
		}

		if (bmove2) {
			if (brepeat) {
				brepeat = false;
				goto RESTARTY_SNAPTO;
			}
		}

	}


	// По оси Oz

	brepeat = true;

	inumboundaryz = 1;


RESTARTZ_SNAPTO:

	if (source_indexpopadaniqnagranXY != nullptr) {
		delete[] source_indexpopadaniqnagranXY;
		source_indexpopadaniqnagranXY = nullptr;
	}
	if (rzboundary != nullptr) {
		delete[] rzboundary;
		rzboundary = nullptr;
	}


	inumboundaryz = 1;
	rzboundary = new doublereal[1 + inumboundaryz]; // число границ по оси y

	
	if (iunion_id_p1 == 0) {
		// добавить границы кабинета
		rzboundary[0] = b[0].g.zS; // начало области 
		rzboundary[inumboundaryz] = b[0].g.zE; // конец области 
	}
	else {
		// добавить границы юниона.
		rzboundary[0] = my_union[iunion_id_p1 - 1].zS; // начало области
		rzboundary[inumboundaryz] = my_union[iunion_id_p1 - 1].zE; // конец области
	}


										   // блоки
	for (i = 1; i<lb; i++) {
		if (b[i].iunion_id == iunion_id_p1) {
			doublereal z_1 = b[i].g.zS;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz,z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
			z_1 = b[i].g.zE;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
			if (bcylinder_meshing && (b[i].g.itypegeom == CYLINDER)) {
				// Cylinder
				switch (b[i].g.iPlane) {
				case XZ_PLANE: case YZ_PLANE:
					for (doublereal dm = dm_start; dm < 0.98; dm = dm + dm_start) {
						z_1 = (b[i].g.zC - b[i].g.R_out_cyl + dm * 2.0 * b[i].g.R_out_cyl);
						if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
							addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
						}
					}
					break;
				}
			}
		}
	}

	// Надо выяснить попал ли источник на границу двух блоков,
	// если попал то это плохо и его придётся сдвигать.
	source_indexpopadaniqnagranXY = new bool[ls];
	for (i = 0; i < ls; i++) {
		source_indexpopadaniqnagranXY[i] = false;
	}
	for (i = 0; i < ls; i++) {
		if (s[i].iunion_id == iunion_id_p1) {
			if (s[i].iPlane == XY_PLANE) {
				for (integer i1 = 0; i1 <= inumboundaryz; i1++) {
					s[i].g.zS = s[i].g.zE;
					if (fabs(s[i].g.zS - rzboundary[i1]) < 1.0e-36) {
						source_indexpopadaniqnagranXY[i] = true;
					}
				}
			}
		}
	}

	if (iunion_id_p1 == 0) {
		// Добавляем границы юнионов.
		for (i = 0; i < lu; i++) {
			doublereal z_1 = my_union[i].zS;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
			z_1 = my_union[i].zE;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// источники
	for (i = 0; i<ls; i++) {
		if (s[i].iunion_id == iunion_id_p1) {
			doublereal z_1 = s[i].g.zS;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
			z_1 = s[i].g.zE;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// стенки
	for (i = 0; i<lw; i++) {
		if (w[i].iunion_id == iunion_id_p1) {
			doublereal z_1 = w[i].g.zS;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
			z_1 = w[i].g.zE;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(rzboundary, inumboundaryz, z_1,XY_PLANE, b, lb, w, lw, s, ls);
			}
		}
	}

	// упорядочивание по возрастанию
	//BubbleEnhSort(rzboundary, inumboundaryz);
	Sort_method(rzboundary, inumboundaryz);


	// snap to grid
	if (bsnap_TO) {
		doublereal rmindist = 1.0e30;

		for (i = 0; i < lb; i++) {
			if (b[i].iunion_id == iunion_id_p1) {
				if (bcylinder_meshing && (b[i].g.itypegeom == CYLINDER)) {
					// Cylinder
					switch (b[i].g.iPlane) {
					case XZ_PLANE: case YZ_PLANE:
						if (fabs(2.0*b[i].g.R_out_cyl) / dcount_diametr_cylinder < rmindist) {
							rmindist = fabs(2.0*b[i].g.R_out_cyl) / dcount_diametr_cylinder;
						}
						break;
					}
				}
				else if (b[i].g.itypegeom == POLYGON) {
					doublereal dist = 1.0e30;
					integer i_7;
					switch (b[i].g.iPlane_obj2) {
					case XZ_PLANE:
						for ( i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
							dist = sqrt((b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) * (b[i].g.xi[i_7 + 1] - b[i].g.xi[i_7]) + (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]) * (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						{
							dist = sqrt((b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) * (b[i].g.xi[b[i].g.nsizei - 1] - b[i].g.xi[0]) + (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]) * (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						break;
					case YZ_PLANE:
						for (i_7 = 0; i_7 < b[i].g.nsizei - 1; i_7++) {
							dist = sqrt((b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) * (b[i].g.yi[i_7 + 1] - b[i].g.yi[i_7]) + (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]) * (b[i].g.zi[i_7 + 1] - b[i].g.zi[i_7]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						{
							dist = sqrt((b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) * (b[i].g.yi[b[i].g.nsizei - 1] - b[i].g.yi[0]) + (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]) * (b[i].g.zi[b[i].g.nsizei - 1] - b[i].g.zi[0]));
							if (dist < rmindist) {
								rmindist = dist;
							}
						}
						break;
					case XY_PLANE:
				        if (fabs(b[i].g.hi[0]) < rmindist) {
					        rmindist = fabs(b[i].g.hi[0]);
			         	}
						break;
					}
				}
				else if (b[i].g.itypegeom == CAD_STL) {

					doublereal dist47 = b[i].g.min_size_edge();
					if (dist47 < rmindist) {
						rmindist = dist47;
					}
				}
				else if ((b[i].g.itypegeom == PRISM)&&(fabs(b[i].g.zE - b[i].g.zS) < rmindist)) {
					rmindist = fabs(b[i].g.zE - b[i].g.zS);
				}
			}
		}
		for (i = 0; i < ls; i++) {
			if (s[i].iunion_id == iunion_id_p1) {
				if ((s[i].iPlane == XZ_PLANE) || ((s[i].iPlane == YZ_PLANE))) {
					if (fabs(s[i].g.zE - s[i].g.zS) < rmindist) {
						rmindist = fabs(s[i].g.zE - s[i].g.zS);
					}
				}
			}
		}
		for (i = 0; i < lw; i++) {
			if (w[i].iunion_id == iunion_id_p1) {
				if ((w[i].iPlane == XZ_PLANE) || (w[i].iPlane == YZ_PLANE)) {
					if (fabs(w[i].g.zE - w[i].g.zS) < rmindist) {
						rmindist = fabs(w[i].g.zE - w[i].g.zS);
					}
				}
			}
		}

		if (iunion_id_p1 == 0) {
			// Добавляем границы юнионов.
			for (i = 0; i < lu; i++) {
				if (fabs(my_union[i].zE - my_union[i].zS) < rmindist) {
					rmindist = fabs(my_union[i].zE - my_union[i].zS);
				}
			}
		}

		if (minimum_fluid_gap_z < rmindist) rmindist = minimum_fluid_gap_z;
		
		rmindist *= snap_to_multiplyer;
		doublereal movetopos;
		doublereal changepos;
		bool bmove = false;
		bool bmove2 = false;
		for (i = 0; i < inumboundaryz; i++) {
			if (bmove) bmove2 = true;
			bmove = false;
			if (fabs(rzboundary[i + 1] - rzboundary[i]) < rmindist) {
				if (i > 0) {
					if (i + 2 <= inumboundaryz) {
						if (fabs(rzboundary[i - 1] - rzboundary[i]) <= fabs(rzboundary[i + 2] - rzboundary[i + 1])) {
							if (fabs(rzboundary[i - 1] - rzboundary[i]) < fabs(rzboundary[i + 1] - rzboundary[i])) {
								// такого быть скорее всего не может.
								movetopos = rzboundary[i - 1];
								changepos = rzboundary[i];
							}
							else {
								// наиболее вероятно.
								movetopos = rzboundary[i + 1];
								changepos = rzboundary[i];
							}
							bmove = true;
						}
						else {
							if (fabs(rzboundary[i + 2] - rzboundary[i + 1]) < fabs(rzboundary[i + 1] - rzboundary[i])) {
								// такого быть скорее всего не может.
								movetopos = rzboundary[i + 2];
								changepos = rzboundary[i + 1];
							}
							else {
								// наиболее вероятно.
								movetopos = rzboundary[i];
								changepos = rzboundary[i + 1];
							}

							bmove = true;
						}
					}
					else {
						movetopos = rzboundary[inumboundaryz];
						changepos = rzboundary[inumboundaryz - 1];
						bmove = true;
					}
				}
				else {
					movetopos = rzboundary[0];
					changepos = rzboundary[1];
					bmove = true;
				}
			}
			if (bmove) {
				//printf(" Z changepos=%e to movetopos=%e \n", changepos, movetopos);
				std::cout << " Z changepos=" << changepos << " to movetopos=" << movetopos << std::endl;
				for (integer i_2 = 0; i_2 < lb; i_2++) {
					if (b[i_2].iunion_id == iunion_id_p1) {
						if (fabs(b[i_2].g.zE - changepos) < 1.0e-33) {
							if (b[i_2].g.itypegeom == PRISM) {
								//printf("block zE=%e zS=%e movetopos%e\n", b[i_2].g.zE, b[i_2].g.zS, movetopos);
								std::cout << "block zE=" << b[i_2].g.zE << " zS=" << b[i_2].g.zS << " movetopos " << movetopos << std::endl;
								b[i_2].g.zE = movetopos;
							}
						}
						if (fabs(b[i_2].g.zS - changepos) < 1.0e-33) {
							if (b[i_2].g.itypegeom == PRISM) {
								//printf("block zS=%e zE=%e movetopos=%e\n", b[i_2].g.zS, b[i_2].g.zE, movetopos);
								std::cout << "block zS=" << b[i_2].g.zS << " zE=" << b[i_2].g.zE << " movetopos " << movetopos << std::endl;
								b[i_2].g.zS = movetopos;
							}
						}
					}
				}
				for (integer i_2 = 0; i_2 < lw; i_2++) {
					if (w[i_2].iunion_id == iunion_id_p1) {
						if (fabs(w[i_2].g.zE - changepos) < 1.0e-33) {
							printf("wall zE\n");
							w[i_2].g.zE = movetopos;
						}
						if (fabs(w[i_2].g.zS - changepos) < 1.0e-33) {
							printf("wall zS\n");
							w[i_2].g.zS = movetopos;
						}
					}
				}
				for (integer i_2 = 0; i_2 < ls; i_2++) {
					if (s[i_2].iunion_id == iunion_id_p1) {
						if (fabs(s[i_2].g.zE - changepos) < 1.0e-33) {
							printf("source zE\n");
							s[i_2].g.zE = movetopos;
						}
						if (fabs(s[i_2].g.zS - changepos) < 1.0e-33) {
							printf("source zS\n");
							s[i_2].g.zS = movetopos;
						}
					}
				}
			}
		}

		if (bmove2) {
			if (brepeat) {
				brepeat = false;
				goto RESTARTZ_SNAPTO;
			}
		}
		
	}

	brepeat = false;

}


/* генерирует простейшую расчётную сетку hex cartesian - 
 * структурированная прямоугольная сетка в трёхмерной расчётной 
 * области в виде кирпича. 
*/
void simplemeshgen(doublereal* &xpos, doublereal* &ypos, doublereal* &zpos, int &inx, int &iny, int &inz,
				   int lb, int ls, int lw,   BLOCK* b, SOURCE* &s, WALL* w, int lu, UNION* &my_union, TPROP* matlist,
	doublereal* &xposadd, doublereal* &yposadd, doublereal* &zposadd,
	int &inxadd, int &inyadd, int &inzadd, integer &iunion_id_p1)
{

	

	// Значение 0.1 подобрано и его лучше не трогать.
	doublereal deltavolkov = 0.1; // Неравномерная сетка в кубе по К.Н. Волкову.

	//bool bsnap_TO = bsnap_TO_global; // snap to grid (избавляет от щелей в сложных моделях).
	//doublereal snap_to_multiplyer = 0.3;

	bool bgeom=true; // если true, то используется неравномерная сетка по закону геометрической прогрессии.
	if (lb <= 2) {
		bgeom = false;
	}
	doublereal q=1.02; // 1.05 - очень слабо неравномерная, 1.1, 1.2 - умеренно неравномерная, 1.25, 1.5, 2 - сильно неравномерная. 
	// для низкорейнольдсовых моделей турбулентности знаменатель геометрической прогрессии должен быть не больше 1.3.

	//bool brepeat = true;

	// По Оси Ox
	doublereal *rxboundary = nullptr; // массив обязательных границ
	int inumboundaryx = 1;
	
	doublereal *ryboundary = nullptr; // массив обязательных границ
	int inumboundaryy = 1;

	doublereal *rzboundary = nullptr; // массив обязательных границ
	int inumboundaryz = 1;

	// подготовка к вычислению зазоров minimum fluid gap.
	calc_minimum_fluid_gap1(inumboundaryx, rxboundary, inumboundaryy, ryboundary, inumboundaryz, rzboundary,
		lb, ls, lw, b, s, w, lu, my_union, iunion_id_p1);

	doublereal minimum_fluid_gap_x = 1.0e36;
	doublereal minimum_fluid_gap_y = 1.0e36;
	doublereal minimum_fluid_gap_z = 1.0e36;

	

	// Непосредтвенное вычисление зазоров minimum fluid gap.
	calc_minimum_fluid_gap2(inumboundaryx, rxboundary, inumboundaryy, ryboundary,
		inumboundaryz, rzboundary, minimum_fluid_gap_x, minimum_fluid_gap_y, minimum_fluid_gap_z, 
		lb,  ls,  lw,  b, s, w, lu, my_union, iunion_id_p1);


	

	bool *source_indexpopadaniqnagranYZ = nullptr;
	bool *source_indexpopadaniqnagranXY = nullptr;
	bool *source_indexpopadaniqnagranXZ = nullptr;

	// 12.03.2017
	// реализация snap to
	// уменьшающает размерность сеточной модели на 33%.
	if (b_adhesion_Mesh) {
		snap_to_moving(source_indexpopadaniqnagranYZ,
			source_indexpopadaniqnagranXY,
			source_indexpopadaniqnagranXZ,
			rxboundary, ryboundary, rzboundary,
			inumboundaryx, inumboundaryy, inumboundaryz,
			minimum_fluid_gap_x, minimum_fluid_gap_y, minimum_fluid_gap_z,
			lb, ls, lw, b, s, w, lu, my_union, iunion_id_p1);
	}

	integer i;

	

    int *ixintervalcount; // число интервалов
	ixintervalcount = new int [inumboundaryx]; // на один меньше чем число границ.
	doublereal alphascal=1.0;
	int inowintervalcount;
	for (i=0; i<(inumboundaryx); i++) {
         alphascal=(rxboundary[i+1]-rxboundary[i])/(rxboundary[inumboundaryx]-rxboundary[0]);
         inowintervalcount=(int)(alphascal*inx);
		 if (inowintervalcount < min_elem_in_x_element) inowintervalcount=min_elem_in_x_element;


		 // FLUID
		 //bool b2div = false;
		 for (integer i_3 = 0; i_3 < (inumboundaryy); i_3++) {
			 for (integer i_4 = 0; i_4 < (inumboundaryz); i_4++) {
				 doublereal yp_3 = 0.5*(ryboundary[i_3 + 1] + ryboundary[i_3]);
				 doublereal zp_3 = 0.5*(rzboundary[i_4 + 1] + rzboundary[i_4]);
				 doublereal xp_1 = 0.5*(rxboundary[i + 1] + rxboundary[i]);
				 doublereal xp_2 = xp_1;
				 doublereal xp_3 = xp_1;
				 if (i < inumboundaryx - 1) {
					 xp_2 = 0.5*(rxboundary[i + 2] + rxboundary[i + 1]);
				 }
				 if (i > 0) {
					 xp_3 = 0.5*(rxboundary[i - 1] + rxboundary[i]);
				 }
				 // определяет номер блока по координате точки.
				 //myisblock_id(integer lb, BLOCK* &b, doublereal x11, doublereal y11, doublereal z11)
				 if (!((((b[myisblock_id(lb, b, xp_1, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||
					 (b[myisblock_id(lb, b, xp_1, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) && 
					 ((b[myisblock_id(lb, b, xp_2, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||
					 (b[myisblock_id(lb, b, xp_2, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) &&
						 ((b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||
					 (b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID))) || 
							 ((b[myisblock_id(lb, b, xp_1, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID) &&
					 (b[myisblock_id(lb, b, xp_2, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID) && 
								 (b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) || 
								 (((b[myisblock_id(lb, b, xp_1, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||
					 (b[myisblock_id(lb, b, xp_1, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) && 
									 ((b[myisblock_id(lb, b, xp_2, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||
									 (b[myisblock_id(lb, b, xp_2, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) &&
										 ((b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||
									 (b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)))))
				 {
					 // нету подряд трех кубиков вдоль линии Ох принадлежащих одновременно solid или fluid.
					 //b2div = true;
				 }
			 }
		 }
		 /*
		 if (b2div) {
			// ixintervalcount[i] = 2;
			ixintervalcount[i] = inowintervalcount;
		 }
		 else {
			 printf("inc X");
			 //getchar();
			 ixintervalcount[i] = 1;
		 }
		 */
		 // 20 02 2017
		 ixintervalcount[i] = inowintervalcount;

         
	}
    /* // debug
	for (i=0; i<(inumboundaryx); i++) {
	#if doubleintprecision == 1
		printf("%lld  ",ixintervalcount[i]);
	#else
		printf("%d  ",ixintervalcount[i]);
	#endif
         
	} //*/
    //getchar();

	

	// заполнение масива положения узлов.
	int iposmark = 1;
	doublereal dx;
	integer k;
	integer ixoldsize=0;
	SetLength(xpos,ixoldsize,1);
    ixoldsize=1;
	doublereal max_size_ratio_x = 1.0;
	for (i=0; i<(inumboundaryx); i++) 
	{
		if ((ixintervalcount[i]-2)%2==0) ixintervalcount[i]++; // чтобы число узлов было всегда чётное.
		integer n2=(integer)((ixintervalcount[i]-1)/2); // округление в меньшую сторону
		doublereal qn2=q;
		for (integer i1=1; i1<n2; i1++) qn2*=q; // возведение в степень.
		doublereal b1=(rxboundary[i+1]-rxboundary[i])*(q-1.0)/(2.0*(qn2-1.0));
		//printf("length=%e\n",(rxboundary[i+1]-rxboundary[i]));
		//getchar(); // OK

        dx=(rxboundary[i+1]-rxboundary[i])/(ixintervalcount[i]-1);
		SetLength(xpos,ixoldsize,(iposmark+ixintervalcount[i]));
        ixoldsize=(iposmark+ixintervalcount[i]);
#if doubleintprecision == 1
		//printf("%lld  ",ixoldsize);// debug
#else
		//printf("%d  ",ixoldsize);// debug
#endif
       
        for (k=iposmark; k<=iposmark+ixintervalcount[i]-2; k++)
		{
			if (!bgeom) {
				if (lb <= 2 ) {
					
					doublereal betavolkov = 1.0 / (sqrt(1.0 - deltavolkov));
					doublereal nvolkow = ((doublereal)(k - iposmark)) / ((doublereal)(ixintervalcount[i] - 1.0));
					xpos[k] = (rxboundary[i+1]- rxboundary[i]) * fmax(0.0,((betavolkov + 1.0)*pow((betavolkov + 1.0) / (betavolkov - 1.0), 2.0*nvolkow - 1.0) - (betavolkov - 1.0)));
					xpos[k] /= 2.0*(1.0+pow((betavolkov + 1.0) / (betavolkov - 1.0), 2.0*nvolkow - 1.0));
					xpos[k] += rxboundary[i];					
				}
				else {
					// равномерная сетка
					xpos[k] = rxboundary[i] + (k - iposmark)*dx;
				}
			}
			else
			{
				// неравномерная сетка
				integer ic1=k-iposmark;
				doublereal madd=b1;
				if (ic1<=n2) {
					// сначала левая сторона.
					for (integer i1=1; i1<ic1; i1++) madd*=q; 
				} else 
				{
					// потом правая сторона.
					for (integer i1=2*n2-ic1; i1>0; i1--) madd*=q;
				}
				if (k==iposmark) xpos[k]=rxboundary[i];
				else xpos[k]=xpos[k-1]+madd;
				
			}
		}
		iposmark=iposmark+ixintervalcount[i]-1;
	}
	SetLength(xpos,ixoldsize,iposmark+1);
	xpos[iposmark]=rxboundary[inumboundaryx];
	inx=iposmark;
    for (i=0; i<inx; i++) xpos[i]=xpos[i+1]; // сдвиг влево на 1
    SetLength(xpos,inx+1,inx); 
	inx--; // индексация начинается с нуля и заканчивается значением inx
	
   // 3 сентября 2017.
   // Добавление новых сеточных линий.
   // необходимость их добавления обусловлена требованнем АЛИС сеток.
	for (int i_28 = 0; i_28 <= inxadd; i_28++) {
		//SetLength(xpos, inx + 1, inx + 2);
		//xpos[inx + 1] = xposadd[i_28];
		//inx++;
		// Это добавление умное, препятствует возникновению одинаковых значений в массиве.
		addboundary(xpos, inx, xposadd[i_28],YZ_PLANE, b, lb, w, lw, s, ls);
	}
	Sort_method<doublereal>(xpos,inx);

	

	for (i=0; i<adapt_x; i++) simplecorrect_meshgen_x(xpos, inx, lb, ls, lw, b, s, w);

	//const doublereal etalon_max_size_ratio=2.0;
	integer inum_iter_ratio_good=0;
	while (1) {

	// Вычисление max size ratio x axis.
	max_size_ratio_x=1.0;
	for (i=0; i<inx-1; i++) {
		
		doublereal dmax=0.0;
		doublereal dmin=0.0;
		if (fabs(xpos[i+1]-xpos[i])>fabs(xpos[i+2]-xpos[i+1])) {
			dmax=fabs(xpos[i+1]-xpos[i]);
			dmin=fabs(xpos[i+2]-xpos[i+1]);
		}
		else {
			dmax=fabs(xpos[i+2]-xpos[i+1]);
			dmin=fabs(xpos[i+1]-xpos[i]);
		}
		if (dmax/dmin>max_size_ratio_x) {
			max_size_ratio_x=dmax/dmin;
		}
	}
	//printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	//printf("x axis max size ratio is equal = %1.4f",max_size_ratio_x);
	if (max_size_ratio_x != max_size_ratio_x) {
		//printf("x axis max size ratio is equal = %1.4f\n", max_size_ratio_x);
		std::cout << "x axis max size ratio is equal =" << max_size_ratio_x << std::endl;
		system("PAUSE");
	}
	//getchar();

	if (max_size_ratio_x < (etalon_max_size_ratio*1.1)) {
		break;
	}

	// Корректировка max size ratio x axis.
	// сделано 9 апреля 2013 года.
	
	bool bplus=false;
	max_size_ratio_x=1.0;
	for (i=0; i<inx-1; i++) {
		
		bplus=false;
		doublereal dmax=0.0;
		doublereal dmin=0.0;
		if (fabs(xpos[i+1]-xpos[i])>fabs(xpos[i+2]-xpos[i+1])) {
			dmax=fabs(xpos[i+1]-xpos[i]);
			dmin=fabs(xpos[i+2]-xpos[i+1]);
		}
		else {
			bplus=true;
			dmax=fabs(xpos[i+2]-xpos[i+1]);
			dmin=fabs(xpos[i+1]-xpos[i]);
		}
		if (dmax/dmin>etalon_max_size_ratio) {
			doublereal pos_candidate;
			if (bplus) {
				pos_candidate=xpos[i+1]+0.5*dmax;
			}
			else {
				pos_candidate=xpos[i]+0.5*dmax;
			}
			SetLength(xpos,inx+1,inx+2);
			xpos[inx+1]=pos_candidate;
            inx=inx+1;
			//BubbleEnhSort<doublereal>(xpos, 0, inx);
			Sort_method<doublereal>(xpos,inx);
			break;
		}
	}

	inum_iter_ratio_good++;

	//getchar();

	}
	
	//printf("x axis max size ratio is equal = %1.4f\n", max_size_ratio_x);
	std::cout << "x axis max size ratio is equal = " << max_size_ratio_x << std::endl;

#if doubleintprecision == 1
	printf("inum_iter_ratio_good is %lld\n", inum_iter_ratio_good);
#else
	printf("inum_iter_ratio_good is %d\n", inum_iter_ratio_good);
#endif
	
	if (bwait) {
	   //getchar();
		system("pause");

	}

	/* // debug
	for (i=0; i<=inx; i++) {
        printf("%f  ",xpos[i]);   
	}
	getchar(); //*/ 


   


    int *iyintervalcount; // число интервалов
    iyintervalcount = new int [inumboundaryy]; // на один меньше чем число границ.
    for (i=0; i<(inumboundaryy); i++) {
         alphascal=(ryboundary[i+1]-ryboundary[i])/(ryboundary[inumboundaryy]-rxboundary[0]);
         inowintervalcount=(int)(alphascal*iny);
		 if (inowintervalcount < min_elem_in_y_element) inowintervalcount=min_elem_in_y_element; 

		 // FLUID
		// bool b2div = false;
		 for (integer i_3 = 0; i_3 < (inumboundaryx); i_3++) {
			 for (integer i_4 = 0; i_4 < (inumboundaryz); i_4++) {
				 doublereal xp_3 = 0.5*(rxboundary[i_3 + 1] + rxboundary[i_3]);
				 doublereal zp_3 = 0.5*(rzboundary[i_4 + 1] + rzboundary[i_4]);
				 doublereal yp_1 = 0.5*(ryboundary[i + 1] + ryboundary[i]);
				 doublereal yp_2 = yp_1;
				 doublereal yp_3 = yp_1;
				 if (i < inumboundaryy - 1) {
					 yp_2 = 0.5*(ryboundary[i + 2] + ryboundary[i + 1]);
				 }
				 if (i > 0) {
					 yp_3 = 0.5*(ryboundary[i - 1] + ryboundary[i]);
				 }
				 // определяет номер блока по координате точки.
				 //myisblock_id(integer lb, BLOCK* &b, doublereal x11, doublereal y11, doublereal z11)
				 if (!((((b[myisblock_id(lb, b, xp_3, yp_1, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[myisblock_id(lb, b, xp_3, yp_1, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) && ((b[myisblock_id(lb, b, xp_3, yp_2, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[myisblock_id(lb, b, xp_3, yp_2, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID))) || ((b[myisblock_id(lb, b, xp_3, yp_1, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[myisblock_id(lb, b, xp_3, yp_2, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) || (((b[myisblock_id(lb, b, xp_3, yp_1, zp_3)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[myisblock_id(lb, b, xp_3, yp_1, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) && ((b[myisblock_id(lb, b, xp_3, yp_2, zp_3)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[myisblock_id(lb, b, xp_3, yp_2, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)))))
				 {
					 // нету подряд трех кубиков вдоль линии Ох принадлежащих одновременно solid или fluid.
					// b2div = true;
				 }
			 }
		 }
		 /*
		 if (b2div) {
			// iyintervalcount[i] = 2;
			iyintervalcount[i] = inowintervalcount;
		 }
		 else {
			 printf("inc Y");
			 //getchar();
			 iyintervalcount[i] = 1;
		 }
		 */
		 // 20 02 2017
		 iyintervalcount[i] = inowintervalcount;
         
	}
    // заполнение массива узлов
    iposmark = 1;
	doublereal dy;
    integer iyoldsize=0;
    SetLength(ypos,iyoldsize,1);
    iyoldsize=1;
	doublereal max_size_ratio_y = 1.0;
	for (i=0; i<inumboundaryy; i++) 
	{

		// чтобы число узлов було всегда чётное
		if ((iyintervalcount[i]-2)%2==0) iyintervalcount[i]++;
		int n2=(int)((iyintervalcount[i]-1)/2); // округление в меньшую сторону
		doublereal qn2=q;
		for (integer i1=1; i1<n2; i1++) qn2*=q; // возведение в степень.
		doublereal b1=(ryboundary[i+1]-ryboundary[i])*(q-1.0)/(2.0*(qn2-1.0));

        dy=(ryboundary[i+1]-ryboundary[i])/(iyintervalcount[i]-1);
		SetLength(ypos,iyoldsize,iposmark+iyintervalcount[i]);
        iyoldsize=iposmark+iyintervalcount[i];
        for (k=iposmark; k<=iposmark+iyintervalcount[i]-2; k++) 
		{
			if (!bgeom) {
				if (lb <= 2) {

					doublereal betavolkov = 1.0 / (sqrt(1.0 - deltavolkov));
					doublereal nvolkow = ((doublereal)(k - iposmark)) / ((doublereal)(iyintervalcount[i] - 1.0));
					ypos[k] = (ryboundary[i+1] - ryboundary[i]) * fmax(0.0, (((betavolkov + 1.0)*(pow(((betavolkov + 1.0) / (betavolkov - 1.0)), 2.0*nvolkow - 1.0)) - (betavolkov - 1.0))));
					ypos[k] /= 2.0*(1.0+pow( (betavolkov + 1.0) / (betavolkov - 1.0), 2.0*nvolkow - 1.0));
					ypos[k] += ryboundary[i];

				}
				else {
					// равномерная сетка
					ypos[k] = ryboundary[i] + (k - iposmark)*dy;
				}
			}
			else 
			{
				// неравномерная сетка
				// на основе геометрической прогрессии с двумя сгущениями.
				integer ic1=k-iposmark;
				doublereal madd=b1;
				if (ic1<=n2) {
					for (integer i1=1; i1<ic1; i1++) madd*=q;
				} else 
				{
					for (integer i1=2*n2-ic1; i1>0; i1--) madd*=q;
				}
				if (k==iposmark) ypos[k]=ryboundary[i];
				else ypos[k]=ypos[k-1]+madd;
			}
		}
		iposmark=iposmark+iyintervalcount[i]-1;
	}
	SetLength(ypos,iyoldsize,iposmark+1);
	ypos[iposmark]=ryboundary[inumboundaryy];
	iny=iposmark;
    for (i=0; i<iny; i++) ypos[i]=ypos[i+1]; // сдвиг влево на 1
    SetLength(ypos,iny+1,iny); 
	iny--; // индексация начинается с нуля и заканчивается значением iny
    

	// 3 сентября 2017.
	// Добавление новых сеточных линий.
	// необходимость их добавления обусловлена требованнем АЛИС сеток.
	for (int i_28 = 0; i_28 <= inyadd; i_28++) {
		//SetLength(ypos, iny + 1, iny + 2);
		//ypos[iny + 1] = yposadd[i_28];
		//iny++;
		// Это добавление умное, препятствует возникновению одинаковых значений в массиве.
		addboundary(ypos, iny, yposadd[i_28],XZ_PLANE, b, lb, w, lw, s, ls);
	}
	Sort_method<doublereal>(ypos,iny);

    for (i=0; i<adapt_y; i++) simplecorrect_meshgen_y(ypos, iny, lb, ls, lw, b, s, w);

	inum_iter_ratio_good=0;
	while (1) {

	max_size_ratio_y=1.0;
	for (i=0; i<iny-1; i++) {
		
		doublereal dmax=0.0;
		doublereal dmin=0.0;
		if (fabs(ypos[i+1]-ypos[i])>fabs(ypos[i+2]-ypos[i+1])) {
			dmax=fabs(ypos[i+1]-ypos[i]);
			dmin=fabs(ypos[i+2]-ypos[i+1]);
		}
		else {
			dmax=fabs(ypos[i+2]-ypos[i+1]);
			dmin=fabs(ypos[i+1]-ypos[i]);
		}
		if (dmax/dmin>max_size_ratio_y) {
			max_size_ratio_y=dmax/dmin;
		}
	}
	
	if (max_size_ratio_y != max_size_ratio_y) {
		//printf("y axis max size ratio is equal = %1.4f\n", max_size_ratio_y);
		std::cout << "y axis max size ratio is equal = " << max_size_ratio_y << std::endl;
		system("PAUSE");
	}
	//getchar();
	
	if (max_size_ratio_y < (etalon_max_size_ratio*1.1)) {
		break;
	}

	// Корректировка max size ratio y axis.
	// сделано 9 апреля 2013 года.
	
	bool bplus=false;
	max_size_ratio_y=1.0;
	for (i=0; i<iny-1; i++) {
		
		bplus=false;
		doublereal dmax=0.0;
		doublereal dmin=0.0;
		if (fabs(ypos[i+1]-ypos[i])>fabs(ypos[i+2]-ypos[i+1])) {
			dmax=fabs(ypos[i+1]-ypos[i]);
			dmin=fabs(ypos[i+2]-ypos[i+1]);
		}
		else {
			bplus=true;
			dmax=fabs(ypos[i+2]-ypos[i+1]);
			dmin=fabs(ypos[i+1]-ypos[i]);
		}
		if (dmax/dmin>etalon_max_size_ratio) {
			doublereal pos_candidate;
			if (bplus) {
				pos_candidate=ypos[i+1]+0.5*dmax;
			}
			else {
				pos_candidate=ypos[i]+0.5*dmax;
			}
			SetLength(ypos,iny+1,iny+2);
			ypos[iny+1]=pos_candidate;
            iny=iny+1;
			//BubbleEnhSort<doublereal>(ypos, 0, iny);
			Sort_method<doublereal>(ypos,iny);
			break;
		}
	}

	inum_iter_ratio_good++;

	//getchar();

	}
	
	//printf("y axis max size ratio is equal = %1.4f\n", max_size_ratio_y);
	std::cout << "y axis max size ratio is equal = " << max_size_ratio_y << std::endl;

#if doubleintprecision == 1
	printf("inum_iter_ratio_good is %lld\n", inum_iter_ratio_good);
#else
	printf("inum_iter_ratio_good is %d\n", inum_iter_ratio_good);
#endif
	
	if (bwait) {
	   //getchar();
		system("pause");
	}

	/* // debug
	for (i=0; i<=iny; i++) {
        printf("%f  ",ypos[i]);   
	}//*/

    // По оси Oz
     
   
	//brepeat = false;


    int *izintervalcount; // число интервалов
    izintervalcount = new int [inumboundaryz]; // на один меньше чем число границ.
    for (i=0; i<(inumboundaryz); i++) {
         alphascal=(rzboundary[i+1]-rzboundary[i])/(rzboundary[inumboundaryz]-rzboundary[0]);
         inowintervalcount=(int)(alphascal*inz);
		 if (inowintervalcount < min_elem_in_z_element) inowintervalcount=min_elem_in_z_element;

		 // FLUID
		 bool b2div = false;
		 for (integer i_3 = 0; i_3 < (inumboundaryx); i_3++) {
			 for (integer i_4 = 0; i_4 < (inumboundaryy); i_4++) {
				 doublereal xp_3 = 0.5*(rxboundary[i_3 + 1] + rxboundary[i_3]);
				 doublereal yp_3 = 0.5*(ryboundary[i_4 + 1] + ryboundary[i_4]);
				 doublereal zp_1 = 0.5*(rzboundary[i + 1] + rzboundary[i]);
				 doublereal zp_2 = zp_1;
				 doublereal zp_3 = zp_1;
				 if (i < inumboundaryz - 1) {
					 zp_2 = 0.5*(rzboundary[i + 2] + rzboundary[i + 1]);
				 }
				 if (i > 0) {
					 zp_3 = 0.5*(rzboundary[i - 1] + rzboundary[i]);
				 }
				 // определяет номер блока по координате точки.
				 //myisblock_id(integer lb, BLOCK* &b, doublereal x11, doublereal y11, doublereal z11)
				 if (!((((b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID))) || ((b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) || (((b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)))))
				 {
					 // нету подряд трех кубиков вдоль линии Ох принадлежащих одновременно solid или fluid.
					 b2div = true;
				 }
			 }
		 }
		 if (b2div) {
			// izintervalcount[i] = 2;
			// izintervalcount[i] = inowintervalcount;
		 }
		 else {
			 //printf("inc Z");
			// getchar();
			 //izintervalcount[i] = 1;
			 if ((inowintervalcount==4)&&(inumboundaryz==1)) {
				 b_one_cell_z = true;
			 }
		 }
		 // 20 02 2017
		 izintervalcount[i] = inowintervalcount;
		 //printf("izintervalcount=%lld\n", izintervalcount[i]);
         
	}
	//getchar();
    // заполнение массива узлов
    iposmark = 1;
	doublereal dz;
    integer izoldsize=0;
	if (b_one_cell_z) {
		SetLength(zpos, izoldsize, 2);
		zpos[0] = rzboundary[0];
		zpos[1]= rzboundary[inumboundaryz];
		inz = 1;
	}
	else {
		SetLength(zpos, izoldsize, 1);
		izoldsize = 1;
		doublereal max_size_ratio_z = 1.0;
		for (i = 0; i < (inumboundaryz); i++)
		{
			// чтобы число узлов було всегда чётное
			if ((izintervalcount[i] - 2) % 2 == 0) izintervalcount[i]++;
			integer n2 = (integer)((izintervalcount[i] - 1) / 2); // округление в меньшую сторону
			doublereal qn2 = q;
			for (integer i1 = 1; i1 < n2; i1++) qn2 *= q; // возведение в степень.
			doublereal b1 = (rzboundary[i + 1] - rzboundary[i]) * (q - 1.0) / (2.0 * (qn2 - 1.0));

			dz = (rzboundary[i + 1] - rzboundary[i]) / (izintervalcount[i] - 1);
			SetLength(zpos, izoldsize, iposmark + izintervalcount[i]);
			izoldsize = iposmark + izintervalcount[i];
			for (k = iposmark; k <= iposmark + izintervalcount[i] - 2; k++)
			{
				if (!bgeom) {
					if (lb <= 2) {

						doublereal betavolkov = 1.0 / (sqrt(1.0 - deltavolkov));
						doublereal nvolkow = ((doublereal)(k - iposmark)) / ((doublereal)(izintervalcount[i] - 1.0));
						zpos[k] = (rzboundary[i + 1] - rzboundary[i]) * fmax(0.0, (((betavolkov + 1.0) * (pow(((betavolkov + 1.0) / (betavolkov - 1.0)), 2.0 * nvolkow - 1.0)) - (betavolkov - 1.0))));
						zpos[k] /= 2.0 * (1.0 + pow((betavolkov + 1.0) / (betavolkov - 1.0), 2.0 * nvolkow - 1.0));
						zpos[k] += rzboundary[i];
						//printf("%e\n", zpos[k]);
						//getchar();

					}
					else {
						// равномерная сетка
						zpos[k] = rzboundary[i] + (k - iposmark) * dz;
					}
				}
				else
				{
					// неравномерная сетка
					// на основе геометрической прогрессии с двумя сгущениями.
					integer ic1 = k - iposmark;
					doublereal madd = b1;
					if (ic1 <= n2) {
						for (integer i1 = 1; i1 < ic1; i1++) madd *= q;
					}
					else
					{
						for (integer i1 = 2 * n2 - ic1; i1 > 0; i1--) madd *= q;
					}
					if (k == iposmark) zpos[k] = rzboundary[i];
					else zpos[k] = zpos[k - 1] + madd;
				}
			}
			iposmark = iposmark + izintervalcount[i] - 1;
		}
		SetLength(zpos, izoldsize, iposmark + 1);
		zpos[iposmark] = rzboundary[inumboundaryz];
		inz = iposmark;
		for (i = 0; i < inz; i++) zpos[i] = zpos[i + 1]; // сдвиг влево на 1
		SetLength(zpos, inz + 1, inz);
		inz--; // индексация начинается с нуля и заканчивается значением inz

		// 3 сентября 2017.
		// Добавление новых сеточных линий.
		// необходимость их добавления обусловлена требованнем АЛИС сеток.
		for (int i_28 = 0; i_28 <= inzadd; i_28++) {
			//SetLength(zpos, inz + 1, inz + 2);
			//zpos[inz + 1] = zposadd[i_28];
			//inz++;
			// Это добавление умное, препятствует возникновению одинаковых значений в массиве.
			addboundary(zpos, inz, zposadd[i_28], XY_PLANE, b, lb, w, lw, s, ls);
		}
		Sort_method<doublereal>(zpos, inz);

		for (i = 0; i < adapt_z; i++) simplecorrect_meshgen_z(zpos, inz, lb, ls, lw, b, s, w);

		inum_iter_ratio_good = 0;
		while (1) {

			max_size_ratio_z = 1.0;
			for (i = 0; i < inz - 1; i++) {

				doublereal dmax = 0.0;
				doublereal dmin = 0.0;
				if (fabs(zpos[i + 1] - zpos[i]) > fabs(zpos[i + 2] - zpos[i + 1])) {
					dmax = fabs(zpos[i + 1] - zpos[i]);
					dmin = fabs(zpos[i + 2] - zpos[i + 1]);
				}
				else {
					dmax = fabs(zpos[i + 2] - zpos[i + 1]);
					dmin = fabs(zpos[i + 1] - zpos[i]);
				}
				if (dmax / dmin > max_size_ratio_z) {
					max_size_ratio_z = dmax / dmin;
				}
			}

			if (max_size_ratio_z != max_size_ratio_z) {
				//printf("z axis max size ratio is equal = %1.4f\n", max_size_ratio_z);
				std::cout << "z axis max size ratio is equal = " << max_size_ratio_z << std::endl;
				system("PAUSE");
			}
			//getchar();

			if (max_size_ratio_z < (etalon_max_size_ratio * 1.1)) {
				break;
			}

			// Корректировка max size ratio z axis.
			// сделано 9 апреля 2013 года.

			bool bplus = false;
			max_size_ratio_z = 1.0;
			for (i = 0; i < inz - 1; i++) {

				bplus = false;
				doublereal dmax = 0.0;
				doublereal dmin = 0.0;
				if (fabs(zpos[i + 1] - zpos[i]) > fabs(zpos[i + 2] - zpos[i + 1])) {
					dmax = fabs(zpos[i + 1] - zpos[i]);
					dmin = fabs(zpos[i + 2] - zpos[i + 1]);
				}
				else {
					bplus = true;
					dmax = fabs(zpos[i + 2] - zpos[i + 1]);
					dmin = fabs(zpos[i + 1] - zpos[i]);
				}
				if (dmax / dmin > etalon_max_size_ratio) {
					doublereal pos_candidate;
					if (bplus) {
						pos_candidate = zpos[i + 1] + 0.5 * dmax;
					}
					else {
						pos_candidate = zpos[i] + 0.5 * dmax;
					}
					SetLength(zpos, inz + 1, inz + 2);
					zpos[inz + 1] = pos_candidate;
					inz = inz + 1;
					//BubbleEnhSort<doublereal>(zpos, 0, inz);
					Sort_method<doublereal>(zpos, inz);
					break;
				}
			}

			inum_iter_ratio_good++;

			//getchar();

		}



		//printf("z axis max size ratio is equal = %1.4f\n", max_size_ratio_z);
		std::cout << "z axis max size ratio is equal = " << max_size_ratio_z << std::endl;

	}

#if doubleintprecision == 1
	printf("inum_iter_ratio_good is %lld\n", inum_iter_ratio_good);
#else
	printf("inum_iter_ratio_good is %d\n", inum_iter_ratio_good);
#endif
	
	if (bwait) {
	  // getchar();
	  system("pause");

	}

	//for (i = 0; i <= inz; i++) {
	//printf("%e  ", zpos[i]);
	//}
	//getchar();

	integer inz_fix = inz;
	if (b_adhesion_Mesh) {
		// Корректировка источника XY
		for (i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranXY[i]) {
				doublereal xc = 0.5 * (s[i].g.xS + s[i].g.xE);
				doublereal yc = 0.5 * (s[i].g.yS + s[i].g.yE);
				doublereal zg = s[i].g.zS;
				// найти позицию +Z на сетке zpos
				// найти центр кО.
				// Если он принадлежит Solid block то сместить истоник в цент этой клетки.
				// Адаптацию maxsize ratio 2 не делать.
				// Если неуспех ищемцентр КО по -Z и смещаем туда.
				integer i55_found = -2;
				for (integer i55 = 0; i55 <= inz_fix; i55++) {
					if (fabs(zpos[i55] - zg) < 1.0e-36) {
						i55_found = i55;
						break;
					}
				}
				if (i55_found >= 0) {
					if (i55_found < inz_fix) {
						doublereal zg1 = 0.5 * (zg + zpos[i55_found + 1]);
						//printf("zg1=%e\n", zg1);
						std::cout << "zg1=" << zg1 << std::endl;
						integer i56_found = -2;
						for (integer ib55 = 0; ib55 < lb; ib55++) {
							if ((xc > b[ib55].g.xS) && (xc < b[ib55].g.xE) && (yc > b[ib55].g.yS) && (yc < b[ib55].g.yE) && (zg1 > b[ib55].g.zS) && (zg1 < b[ib55].g.zE))
							{
								i56_found = ib55;
							}
						}

						bool bzero_pos = true;
						doublereal zg2 = zpos[0] - 0.5 * fabs(zpos[1] - zpos[0]);
						if (i55_found > 0) {
							zg2 = 0.5 * (zg + zpos[i55_found - 1]);
							bzero_pos = false;
						}
						//printf("zg2=%e\n", zg2);
						std::cout << "zg2=" << zg2 << std::endl;
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (zg2 > b[ib57].g.zS) && (zg2 < b[ib57].g.zE))
							{
								if (zg2 > zpos[0]) {
									i57_found = ib57;
								}
							}
						}

						if (i56_found >= 0) {
							if (b[i56_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								if ((i57_found >= 0) && (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
									// Мы помещаем источник тепла в блок с большей теплопроводностью.
									// comparison_lam выдаёт истину если теплопроводность блока i56 больше.
									if (comparison_lam(matlist, b, i56_found, i57_found, 25.0)) {
										// в блоке i56 теплопроводность выше.
										s[i].g.zS = zg1;
										s[i].g.zE = zg1;
										addboundary(zpos, inz, zg1, XY_PLANE, b, lb, w, lw, s, ls);
									}
									else {
										// в блоке i57 теплопроводность выше.
										s[i].g.zS = zg2;
										s[i].g.zE = zg2;
										addboundary(zpos, inz, zg2, XY_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {

									// Найден Солид Блок.
									//printf("zg1==%e\n", zg1);
									std::cout << "zg1==" << zg1 << std::endl;
									doublereal zgolg = zpos[1];
									if (inz == 1) {
										// Слишком малоразмерная расчётная сетка.
										zgolg = 0.5 * (zpos[0] + zg1);
										addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
										zgolg = 0.5 * (zpos[1] + zg1);
										addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
										zgolg = 0.5 * (zpos[0] + zg1);
										// Расстояние от края вглубь расчётной области две клетки.
										s[i].g.zS = zgolg;
										s[i].g.zE = zgolg;
									}
									else {
										if (bzero_pos) {
											// Источник в позиции 0.0.
											// Разбиваем поподробнее.
											zgolg = 0.5 * (zpos[0] + zg1);
											addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
											zgolg = 0.5 * (zpos[1] + zg1);
											addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
											zgolg = 0.5 * (zpos[0] + zg1);
											// Расстояние от края вглубь расчётной области две клетки.
											s[i].g.zS = zgolg;
											s[i].g.zE = zgolg;
										}
										else {
											s[i].g.zS = zg1;
											s[i].g.zE = zg1;
										}
									}
									addboundary(zpos, inz, zg1, XY_PLANE, b, lb, w, lw, s, ls);



								}
							}
							else {

								if (i57_found >= 0) {
									if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
										// Найден Солид Блок.
										s[i].g.zS = zg2;
										s[i].g.zE = zg2;
										addboundary(zpos, inz, zg2, XY_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
									system("PAUSE");
									exit(1);
								}
							}
						}
					}
					else {
						doublereal zg2 = 0.5 * (zg + zpos[i55_found - 1]);
						//printf("zg2=%e\n", zg2);
						std::cout << "zg2=" << zg2 << std::endl;
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (zg2 > b[ib57].g.zS) && (zg2 < b[ib57].g.zE))
							{
								i57_found = ib57;
							}
						}
						if (i57_found >= 0) {
							if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								// Найден Солид Блок.
								s[i].g.zS = zg2;
								s[i].g.zE = zg2;
								addboundary(zpos, inz, zg2, XY_PLANE, b, lb, w, lw, s, ls);
							}
						}
						else {
							printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
							system("PAUSE");
							exit(1);
						}
					}
				}
			}
		}
	}
	//BubbleEnhSort<doublereal>(zpos, 0, inz);
	Sort_method<doublereal>(zpos,inz);
	//SetLength(zpos, inz + 2, inz);
	//inz++;

	integer inx_fix = inx;
	if (b_adhesion_Mesh) {
		// Корректировка источника YZ
		for (i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranYZ[i]) {
				doublereal zc = 0.5 * (s[i].g.zS + s[i].g.zE);
				doublereal yc = 0.5 * (s[i].g.yS + s[i].g.yE);
				doublereal xg = s[i].g.xS;
				// найти позицию +Z на сетке zpos
				// найти центр кО.
				// Если он принадлежит Solid block то сместить истоник в цент этой клетки.
				// Адаптацию maxsize ratio 2 не делать.
				// Если неуспех ищемцентр КО по -Z и смещаем туда.
				integer i55_found = -2;
				for (integer i55 = 0; i55 <= inx_fix; i55++) {
					if (fabs(xpos[i55] - xg) < 1.0e-36) {
						i55_found = i55;
						break;
					}
				}
				if (i55_found >= 0) {
					if (i55_found < inx_fix) {
						doublereal xg1 = 0.5 * (xg + xpos[i55_found + 1]);
						//printf("xg1=%e\n", xg1);
						std::cout << "xg1=" << xg1 << std::endl;
						integer i56_found = -2;
						for (integer ib55 = 0; ib55 < lb; ib55++) {
							if ((zc > b[ib55].g.zS) && (zc < b[ib55].g.zE) && (yc > b[ib55].g.yS) && (yc < b[ib55].g.yE) && (xg1 > b[ib55].g.xS) && (xg1 < b[ib55].g.xE))
							{
								i56_found = ib55;
							}
						}

						bool bzero_pos = true;
						doublereal xg2 = xpos[0] - 0.5 * fabs(xpos[1] - xpos[0]);
						if (i55_found > 0) {
							xg2 = 0.5 * (xg + xpos[i55_found - 1]);
							bzero_pos = false;
						}

						//printf("xg2=%e\n", xg2);
						std::cout << "xg2=" << xg2 << std::endl;
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (xg2 > b[ib57].g.xS) && (xg2 < b[ib57].g.xE))
							{
								if (xg2 > xpos[0]) {
									i57_found = ib57;
								}
							}
						}


						if (i56_found >= 0) {
							if (b[i56_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								// TODO 11.07.2016
								if ((i57_found >= 0) && (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
									// Мы помещаем источник тепла в блок с большей теплопроводностью.
									// comparison_lam выдаёт истину если теплопроводность блока i56 больше.
									if (comparison_lam(matlist, b, i56_found, i57_found, 25.0)) {
										// в блоке i56 теплопроводность выше.
										s[i].g.xS = xg1;
										s[i].g.xE = xg1;
										addboundary(xpos, inz, xg1, YZ_PLANE, b, lb, w, lw, s, ls);
									}
									else {
										// в блоке i57 теплопроводность выше.
										s[i].g.xS = xg2;
										s[i].g.xE = xg2;
										addboundary(xpos, inz, xg2, YZ_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									// Найден Солид Блок.
									//printf("xg1==%e\n", xg1);
									std::cout << "xg1==" << xg1 << std::endl;
									doublereal xgolg = xpos[1];
									if (inx == 1) {
										// Слишком малоразмерная расчётная сетка.
										xgolg = 0.5 * (xpos[0] + xg1);
										addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
										xgolg = 0.5 * (xpos[1] + xg1);
										addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
										xgolg = 0.5 * (xpos[0] + xg1);
										// Расстояние от края вглубь расчётной области две клетки.
										s[i].g.xS = xgolg;
										s[i].g.xE = xgolg;
									}
									else {
										if (bzero_pos) {
											// Слишком малоразмерная расчётная сетка.
											xgolg = 0.5 * (xpos[0] + xg1);
											addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
											xgolg = 0.5 * (xpos[1] + xg1);
											addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
											xgolg = 0.5 * (xpos[0] + xg1);
											// Расстояние от края вглубь расчётной области две клетки.
											s[i].g.xS = xgolg;
											s[i].g.xE = xgolg;
										}
										else {
											s[i].g.xS = xg1;
											s[i].g.xE = xg1;
										}
									}
									addboundary(xpos, inx, xg1, YZ_PLANE, b, lb, w, lw, s, ls);

								}
							}
							else {

								if (i57_found >= 0) {
									if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
										// Найден Солид Блок.
										s[i].g.xS = xg2;
										s[i].g.xE = xg2;
										addboundary(xpos, inx, xg2, YZ_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
									system("PAUSE");
									exit(1);
								}
							}
						}
					}
					else {
						doublereal xg2 = 0.5 * (xg + xpos[i55_found - 1]);
						//printf("xg2=%e\n", xg2);
						std::cout << "xg2=" << xg2 << std::endl;
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (xg2 > b[ib57].g.xS) && (xg2 < b[ib57].g.xE))
							{
								i57_found = ib57;
							}
						}
						if (i57_found >= 0) {
							if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								// Найден Солид Блок.
								s[i].g.xS = xg2;
								s[i].g.xE = xg2;
								addboundary(xpos, inx, xg2, YZ_PLANE, b, lb, w, lw, s, ls);
							}
						}
						else {
							printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
							system("PAUSE");
							exit(1);
						}
					}
				}
			}
		}
	}
	//BubbleEnhSort<doublereal>(xpos, 0, inx);
	Sort_method<doublereal>(xpos,inx);

	integer iny_fix = iny;
	if (b_adhesion_Mesh) {
		// Корректировка источника XZ
		for (i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranXZ[i]) {
				doublereal xc = 0.5 * (s[i].g.xS + s[i].g.xE);
				doublereal zc = 0.5 * (s[i].g.zS + s[i].g.zE);
				doublereal yg = s[i].g.yS;
				// найти позицию +Z на сетке zpos
				// найти центр кО.
				// Если он принадлежит Solid block то сместить истоник в цент этой клетки.
				// Адаптацию maxsize ratio 2 не делать.
				// Если неуспех ищемцентр КО по -Z и смещаем туда.
				integer i55_found = -2;
				for (integer i55 = 0; i55 <= iny_fix; i55++) {
					if (fabs(ypos[i55] - yg) < 1.0e-36) {
						i55_found = i55;
						break;
					}
				}
				if (i55_found >= 0) {
					if (i55_found < iny_fix) {
						doublereal yg1 = 0.5 * (yg + ypos[i55_found + 1]);
						//printf("yg1=%e\n", yg1);
						std::cout << "yg1=" << yg1 << std::endl;
						integer i56_found = -2;
						for (integer ib55 = 0; ib55 < lb; ib55++) {
							if ((xc > b[ib55].g.xS) && (xc < b[ib55].g.xE) && (zc > b[ib55].g.zS) && (zc < b[ib55].g.zE) && (yg1 > b[ib55].g.yS) && (yg1 < b[ib55].g.yE))
							{
								i56_found = ib55;
							}
						}
						doublereal yg2 = 0.5 * (yg + ypos[i55_found - 1]);
						//printf("yg2=%e\n", yg2);
						std::cout << "yg2=" << yg2 << std::endl;
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yg2 > b[ib57].g.yS) && (yg2 < b[ib57].g.yE))
							{
								if (yg2 > ypos[0]) {
									i57_found = ib57;
								}
							}
						}


						if (i56_found >= 0) {
							if (b[i56_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {

								if ((i57_found >= 0) && (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
									// Мы помещаем источник тепла в блок с большей теплопроводностью.
									// comparison_lam выдаёт истину если теплопроводность блока i56 больше.
									if (comparison_lam(matlist, b, i56_found, i57_found, 25.0)) {
										// в блоке i56 теплопроводность выше.
										s[i].g.yS = yg1;
										s[i].g.yE = yg1;
										addboundary(ypos, iny, yg1, XZ_PLANE, b, lb, w, lw, s, ls);
									}
									else {
										// в блоке i57 теплопроводность выше.
										s[i].g.yS = yg2;
										s[i].g.yE = yg2;
										addboundary(ypos, iny, yg2, XZ_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									// Найден Солид Блок.
									s[i].g.yS = yg1;
									s[i].g.yE = yg1;
									addboundary(ypos, iny, yg1, XZ_PLANE, b, lb, w, lw, s, ls);
								}
							}
							else {

								if (i57_found >= 0) {
									if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
										// Найден Солид Блок.
										s[i].g.yS = yg2;
										s[i].g.yE = yg2;
										addboundary(ypos, iny, yg2, XZ_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
									system("PAUSE");
									exit(1);
								}
							}
						}
					}
					else {
						doublereal yg2 = 0.5 * (yg + ypos[i55_found - 1]);
						//printf("yg2=%e\n", yg2);
						std::cout << "yg2=" << yg2 << std::endl;
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yg2 > b[ib57].g.yS) && (yg2 < b[ib57].g.yE))
							{
								i57_found = ib57;
							}
						}
						if (i57_found >= 0) {
							if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								// Найден Солид Блок.
								s[i].g.yS = yg2;
								s[i].g.yE = yg2;
								addboundary(ypos, iny, yg2, XZ_PLANE, b, lb, w, lw, s, ls);
							}
						}
						else {
							printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
							system("PAUSE");
							exit(1);
						}
					}
				}
			}
		}
	}
	//BubbleEnhSort<doublereal>(ypos, 0, iny);
	Sort_method<doublereal>(ypos,iny);


	// улучшение качества сетки по второму критерию с значением 30.0
	// как во flowvision.
	quolite_refinement(inx, iny, inz, xpos, ypos, zpos);

	// debug

	//for (i=0; i<=inz; i++) {
	//printf("%e  ",zpos[i]);
	//}
	//getchar();//

  
	/*// debug
	for (i=0; i<=inz; i++) {
        printf("%f  ",zpos[i]);   
	}
	getchar();//*/

    // Освобождение оперативной памяти.
	delete[] rxboundary;
	delete[] ryboundary;
	delete[] rzboundary;
	delete[] ixintervalcount;
	delete[] iyintervalcount;
    delete[] izintervalcount;
	delete[] source_indexpopadaniqnagranXY;
	delete[] source_indexpopadaniqnagranYZ;
	delete[] source_indexpopadaniqnagranXZ;

} // simplemeshgen

/*
// возвращает наибольшее из двух чисел
doublereal fmax(doublereal fA, doublereal fB) {
	doublereal r=fA;
	if (fB>r) r=fB;
	return r;
} // fmax
*/

// Балансировка значения n2 так чтобы ячейки 
// прилегающие к середине интервала имели примерно одинаковые размеры.
integer ibalancen2(integer n2param, doublereal qL, doublereal qR, doublereal *rboundary,
	           integer intervalcount, integer i, integer iposmark) {

		// i - номер интервала который подлежит разбиению,

        bool bcontinue=true;
		integer n2=n2param;
		doublereal qgold=fmax(qL,qR);
		// Защита от зацикливания:
		integer ibound=0, imax=1000; // ограничение на максимальное число проходов.

		// итерационная балансировка:
		while ((ibound<imax) && (bcontinue)) 
		{	    

		     doublereal qn2L=qL;
		     for (integer i1=1; i1<n2; i1++) qn2L*=qL; // возведение в степень.
		     doublereal b1L=(rboundary[i+1]-rboundary[i])*(qL-1.0)/(2.0*(qn2L-1.0));
		     doublereal qn2R=qR;
		     for (integer i1=1; i1<(intervalcount-n2); i1++) qn2R*=qR; // возведение в степень.
             doublereal b1R=(rboundary[i+1]-rboundary[i])*(qR-1.0)/(2.0*(qn2R-1.0));


		     doublereal rpos1=0.0, rpos2=0.0, rpos3=0.0;
        
             for (integer k=iposmark; k<=iposmark+intervalcount-2; k++)
		     {
			      // неравномерная сетка
			      integer ic1=k-iposmark;
			      doublereal madd;
			      if (ic1<n2) {
				      // сначала левая сторона.
				      madd=b1L;
				      for (integer i1=0; i1<ic1; i1++) madd*=qL;
					  if (ic1==n2-1) madd=0.5*(rboundary[i+1]+rboundary[i])-rpos3; // чтобы избежать ошибок округления
			      } else 
			      {
				      // потом правая сторона.
				      madd=b1R;
				      for (integer i1=intervalcount-1-ic1; i1>0; i1--) madd*=qR;
			      }
			      if (k==iposmark) {
				     rpos1=rboundary[i];
				     rpos2=rpos1;
				     rpos3=rpos1;
			      }
			      else {
				     rpos1=rpos2;
				     rpos2=rpos3;
				     rpos3=rpos3+madd;
			      }
			      if (ic1==n2) {
				      if (((rpos3-rpos2)/(rpos2-rpos1)<=qgold) && ((rpos3-rpos2)/(rpos2-rpos1)>=1.0/qgold)) {
					      bcontinue=false;
				      } else if ((rpos3-rpos2)/(rpos2-rpos1)<1.0/qgold) {
					      // верхний существенно меньше нижнего
					      // надо увеличить n2;
					      n2+=2;
						  if (n2>intervalcount-4) bcontinue=false;
				      }
					  else if ((rpos3-rpos2)/(rpos2-rpos1)>qgold) {
						  // верхний существенно больше нижнего
						  // надо уменьшить n2;
						  n2-=2;
						  if (n2<5) bcontinue=false;
					  }
				      break; // выход из цикла for
			      }
		     }
			 ibound++; // ограничитель на максимальное число проходов.
		} // end while

    return n2;  
} // ibalancen2

/* генерирует простейшую расчётную сетку hex cartesian - 
 * структурированная прямоугольная сетка в трёхмерной расчётной 
 * области в виде кирпича. Особенность данного генератора в том
 * что он может генерировать как равномерную сетку так и неравномерную.
 * Данный автоматический сеточный генератор содержит больше управляющих параметров чем
 * автоматический сеточный генератор simplemeshgen. 
 * дата реализации 8 февраля 2012.
*/
void unevensimplemeshgen(doublereal* &xpos, doublereal* &ypos, doublereal* &zpos, int &inx, int &iny, int &inz,
	int lb, int ls, int lw, BLOCK* b, SOURCE* s, WALL* w, int lu, UNION* &my_union, 
	TPROP* matlist, doublereal dgx, doublereal dgy, doublereal dgz,
	doublereal* &xposadd, doublereal* &yposadd, doublereal* &zposadd,
	int &inxadd, int &inyadd, int &inzadd, integer &iunion_id_p1)
{

	

	//*****************************************************************************************************************
	// Изменяемые параметры автоматического сеточного генератора
	doublereal deltavolkov = 0.1; // Неравномерная сетка в кубе по Волкову.
	//printf("incomming\n");
	//system("pause");

	// если bgeom  , то используется неравномерная сетка по закону геометрической прогрессии.
	bool bgeomx=false; // по оси Ox
	bool bgeomy=false; // по оси Oy
	bool bgeomz=true; // по оси Oz
    // 1.05 - очень слабо неравномерная, 1.1, 1.2 - умеренно неравномерная, 1.25, 1.5, 2 - сильно неравномерная.
	// для низкорейнольдсовых моделей турбулентности знаменатель геометрической прогрессии должен быть не больше 1.3.
	doublereal qxW=1.0005; // параметр сгущения к левой West стороне. 
	doublereal qxE=1.0005; // параметр сгущения к правой East стороне
	doublereal qyS=1.0005; // South
	doublereal qyN=1.0005; // North
	doublereal qzB=1.35; // Bottom
	doublereal qzT=1.00005; // Top
	//doublereal qzB=1.00005; // Bottom
	//doublereal qzT=1.35; // Top
	//*****************************************************************************************************************
	
	// Распознаём тест Валь-Девиса или задачу Релея Бенара.
	if (lb == 1) {
		if (dgx*dgx + dgy*dgy + dgz*dgz > 1.0e-20) {

			// Не работает 18.02.2017.
			//printf("uneven Davis mesh\n");
			bgeomx = true; // по оси Ox
			bgeomy = true; // по оси Oy
			bgeomz = true; // по оси Oz
								// 1.05 - очень слабо неравномерная, 1.1, 1.2 - умеренно неравномерная, 1.25, 1.5, 2 - сильно неравномерная.
								// для низкорейнольдсовых моделей турбулентности знаменатель геометрической прогрессии должен быть не больше 1.3.
			// 1.35
			qxW = 1.2; // параметр сгущения к левой West стороне. 
			qxE = 1.2; // параметр сгущения к правой East стороне
			qyS = 1.2; // South
			qyN = 1.2; // North
			qzB = 1.2; // Bottom
			qzT = 1.2; // Top
		}
	}

	if (lb == 1) {
		// 20.05.2017 К.Н. Волков
		bgeomx = false; // по оси Ox
		bgeomy = false; // по оси Oy
		bgeomz = false; // по оси Oz
	}

	//bool brepeat = true;

	// По Оси Ox
	doublereal *rxboundary = nullptr; // массив обязательных границ
	int inumboundaryx = 1;

	doublereal *ryboundary = nullptr; // массив обязательных границ
	int inumboundaryy = 1;

	doublereal *rzboundary = nullptr; // массив обязательных границ
	int inumboundaryz = 1;

	// подготовка к вычислению зазоров minimum fluid gap.
	calc_minimum_fluid_gap1(inumboundaryx, rxboundary, inumboundaryy, ryboundary, inumboundaryz, rzboundary,
		lb, ls, lw, b, s, w, lu, my_union, iunion_id_p1);

	doublereal minimum_fluid_gap_x = 1.0e36;
	doublereal minimum_fluid_gap_y = 1.0e36;
	doublereal minimum_fluid_gap_z = 1.0e36;

	

	// Непосредтвенное вычисление зазоров minimum fluid gap.
	calc_minimum_fluid_gap2( inumboundaryx, rxboundary, inumboundaryy, ryboundary,
		inumboundaryz, rzboundary, minimum_fluid_gap_x, minimum_fluid_gap_y, minimum_fluid_gap_z,
		lb, ls, lw, b, s, w, lu, my_union, iunion_id_p1);

	


	bool *source_indexpopadaniqnagranYZ = nullptr;
	bool *source_indexpopadaniqnagranXY = nullptr;
	bool *source_indexpopadaniqnagranXZ = nullptr;

	// 12.03.2017
	// реализация snap to
	// уменьшающает размерность сеточной модели на 33%.
	if (b_adhesion_Mesh) {
		snap_to_moving(source_indexpopadaniqnagranYZ,
			source_indexpopadaniqnagranXY,
			source_indexpopadaniqnagranXZ,
			rxboundary, ryboundary, rzboundary,
			inumboundaryx, inumboundaryy, inumboundaryz,
			minimum_fluid_gap_x, minimum_fluid_gap_y, minimum_fluid_gap_z,
			lb, ls, lw, b, s, w, lu, my_union, iunion_id_p1);
	}

	integer i;

	

    int *ixintervalcount; // число интервалов
	ixintervalcount = new int [inumboundaryx]; // на один меньше чем число границ.
	doublereal alphascal=1.0;
	int inowintervalcount;
	for (i=0; i<(inumboundaryx); i++) {
         alphascal=(rxboundary[i+1]-rxboundary[i])/(rxboundary[inumboundaryx]-rxboundary[0]);
         inowintervalcount=(int)(alphascal*inx);
		 if (inowintervalcount < min_elem_in_x_element) inowintervalcount=min_elem_in_x_element;
         ixintervalcount[i]=inowintervalcount;
	}
    /* // debug
	for (i=0; i<(inumboundaryx); i++) {
	#if doubleintprecision == 1
		printf("%lld  ",ixintervalcount[i]);
	#else
		printf("%d  ",ixintervalcount[i]);
	#endif
        
	} //*/
    //getchar();

	// заполнение масива положения узлов.
	int iposmark = 1;
	doublereal dx;
	int k;
	integer ixoldsize=0;
	SetLength(xpos,ixoldsize,1);
    ixoldsize=1;
	for (i=0; i<(inumboundaryx); i++) 
	{
		if ((ixintervalcount[i]-2)%2==0) ixintervalcount[i]++; // чтобы число узлов было всегда чётное.
		integer n2=(integer)((ixintervalcount[i]-1)/2); // округление в меньшую сторону
		if (bgeomx) n2=ibalancen2(n2, qxW, qxE, rxboundary, ixintervalcount[i], i, iposmark); // балансировка 
		doublereal qn2W=qxW;
		for (integer i1=1; i1<n2; i1++) qn2W*=qxW; // возведение в степень.
		doublereal b1W=(rxboundary[i+1]-rxboundary[i])*(qxW-1.0)/(2.0*(qn2W-1.0));
		doublereal qn2E=qxE;
		for (integer i1=1; i1<(ixintervalcount[i]-n2); i1++) qn2E*=qxE; // возведение в степень.
		doublereal b1E=(rxboundary[i+1]-rxboundary[i])*(qxE-1.0)/(2.0*(qn2E-1.0));
		//printf("length=%e\n",(rxboundary[i+1]-rxboundary[i]));
		//getchar(); // OK

        dx=(rxboundary[i+1]-rxboundary[i])/(ixintervalcount[i]-1);
		SetLength(xpos,ixoldsize,(iposmark+ixintervalcount[i]));
        ixoldsize=(iposmark+ixintervalcount[i]);
#if doubleintprecision == 1
		//printf("%lld  ",ixoldsize);// debug
#else
		//printf("%d  ",ixoldsize);// debug
#endif
        
        for (k=iposmark; k<=iposmark+ixintervalcount[i]-2; k++)
		{
			if (!bgeomx) {
				// равномерная сетка
				if (lb == 1) {

					doublereal betavolkov = 1.0 / (sqrt(1.0 - deltavolkov));
					doublereal nvolkow = ((doublereal)(k - iposmark)) / ((doublereal)(ixintervalcount[i] - 1.0));
					xpos[k] = (rxboundary[1] - rxboundary[0]) * fmax(0.0, (((betavolkov + 1.0)*(pow(((betavolkov + 1.0) / (betavolkov - 1.0)), 2.0*nvolkow - 1.0)) - (betavolkov - 1.0))));
					xpos[k] /= 2.0*(1.0 + pow((betavolkov + 1.0) / (betavolkov - 1.0), 2.0*nvolkow - 1.0));
					xpos[k] += rxboundary[0];
				}
				else {
					xpos[k] = rxboundary[i] + (k - iposmark)*dx;
				}
			}
			else
			{
				// неравномерная сетка
				integer ic1=k-iposmark;
				doublereal madd;
				if (ic1<n2) {
					// сначала левая сторона.
					madd=b1W;
					for (integer i1=1; i1<ic1; i1++) madd*=qxW; 
					if (ic1==n2-1) madd=0.5*(rxboundary[i+1]+rxboundary[i])-xpos[k-1]; // чтобы избежать ошибок округления
				} else 
				{
					// потом правая сторона.
					madd=b1E;
					for (integer i1=ixintervalcount[i]-1-ic1; i1>0; i1--) madd*=qxE;
				}
				if (k==iposmark) xpos[k]=rxboundary[i];
				else xpos[k]=xpos[k-1]+madd;
				
			}
		}
		iposmark=iposmark+ixintervalcount[i]-1;
	}
	SetLength(xpos,ixoldsize,iposmark+1);
	xpos[iposmark]=rxboundary[inumboundaryx];
	inx=iposmark;
    for (i=0; i<inx; i++) xpos[i]=xpos[i+1]; // сдвиг влево на 1
    SetLength(xpos,inx+1,inx); 
	inx--; // индексация начинается с нуля и заканчивается значением inx
	
	// 3 сентября 2017.
	// Добавление новых сеточных линий.
	// необходимость их добавления обусловлена требованнем АЛИС сеток.
	for (int i_28 = 0; i_28 <= inxadd; i_28++) {
		//SetLength(xpos, inx + 1, inx + 2);
		//xpos[inx + 1] = xposadd[i_28];
		//inx++;
		// Это добавление умное, препятствует возникновению одинаковых значений в массиве.
		addboundary(xpos, inx, xposadd[i_28],YZ_PLANE, b, lb, w, lw, s, ls);
	}
	Sort_method<doublereal>(xpos,inx);

	for (i=0; i<adapt_x; i++) simplecorrect_meshgen_x(xpos, inx, lb, ls, lw, b, s, w);


	//*****************************************************************************
	// Выравниватель, такой чтобы maxsize_ratio не превышал 2.0.
	// Выравнивание необходимо для сходимости гидродинамических задач.
	// 11.07.2016
	//const doublereal etalon_max_size_ratio=2.0;
	integer inum_iter_ratio_good = 0;
	doublereal max_size_ratio_x = 1.0;
	while (1) {

		// Вычисление max size ratio x axis.
		max_size_ratio_x = 1.0;
		for (i = 0; i<inx - 1; i++) {

			doublereal dmax = 0.0;
			doublereal dmin = 0.0;
			if (fabs(xpos[i + 1] - xpos[i])>fabs(xpos[i + 2] - xpos[i + 1])) {
				dmax = fabs(xpos[i + 1] - xpos[i]);
				dmin = fabs(xpos[i + 2] - xpos[i + 1]);
			}
			else {
				dmax = fabs(xpos[i + 2] - xpos[i + 1]);
				dmin = fabs(xpos[i + 1] - xpos[i]);
			}
			if (dmax / dmin>max_size_ratio_x) {
				max_size_ratio_x = dmax / dmin;
			}
		}
		
		if (max_size_ratio_x != max_size_ratio_x) {
			//printf("x axis max size ratio is equal = %1.4f\n", max_size_ratio_x);
			std::cout << "x axis max size ratio is equal = " << max_size_ratio_x << std::endl;
			system("PAUSE");
		}
		//getchar();

		if (max_size_ratio_x < (etalon_max_size_ratio*1.1)) {
			break;
		}

		// Корректировка max size ratio x axis.
		// сделано 9 апреля 2013 года.

		bool bplus = false;
		max_size_ratio_x = 1.0;
		for (i = 0; i<inx - 1; i++) {

			bplus = false;
			doublereal dmax = 0.0;
			doublereal dmin = 0.0;
			if (fabs(xpos[i + 1] - xpos[i])>fabs(xpos[i + 2] - xpos[i + 1])) {
				dmax = fabs(xpos[i + 1] - xpos[i]);
				dmin = fabs(xpos[i + 2] - xpos[i + 1]);
			}
			else {
				bplus = true;
				dmax = fabs(xpos[i + 2] - xpos[i + 1]);
				dmin = fabs(xpos[i + 1] - xpos[i]);
			}
			if (dmax / dmin>etalon_max_size_ratio) {
				doublereal pos_candidate;
				if (bplus) {
					pos_candidate = xpos[i + 1] + 0.5*dmax;
				}
				else {
					pos_candidate = xpos[i] + 0.5*dmax;
				}
				SetLength(xpos, inx + 1, inx + 2);
				xpos[inx + 1] = pos_candidate;
				inx = inx + 1;
				//BubbleEnhSort<doublereal>(xpos, 0, inx);
				Sort_method<doublereal>(xpos,inx);
				break;
			}
		}

		inum_iter_ratio_good++;

		//getchar();

	}

	//printf("x axis max size ratio is equal = %1.4f\n", max_size_ratio_x);
	std::cout << "x axis max size ratio is equal = " << max_size_ratio_x << std::endl;

#if doubleintprecision == 1
	printf("inum_iter_ratio_good is %lld\n", inum_iter_ratio_good);
#else
	printf("inum_iter_ratio_good is %d\n", inum_iter_ratio_good);
#endif
	
	if (bwait) {
		//getchar();
		system("pause");

	}

	// Выравниватель, такой чтобы max size_ratio не превышал 2.0.
	//*****************************************************************************
	

	/* // debug
	for (i=0; i<=inx; i++) {
        printf("%f  ",xpos[i]);   
	}
	getchar(); //*/ 


    // По оси Oy
    


    int *iyintervalcount; // число интервалов
    iyintervalcount = new int [inumboundaryy]; // на один меньше чем число границ.
    for (i=0; i<(inumboundaryy); i++) {
         alphascal=(ryboundary[i+1]-ryboundary[i])/(ryboundary[inumboundaryy]-rxboundary[0]);
         inowintervalcount=(int)(alphascal*iny);
		 if (inowintervalcount < min_elem_in_y_element) inowintervalcount=min_elem_in_y_element; 
         iyintervalcount[i]=inowintervalcount;
	}
    // заполнение массива узлов
    iposmark = 1;
	doublereal dy;
    integer iyoldsize=0;
    SetLength(ypos,iyoldsize,1);
    iyoldsize=1;
	for (i=0; i<inumboundaryy; i++) 
	{

		// чтобы число узлов було всегда чётное
		if ((iyintervalcount[i]-2)%2==0) iyintervalcount[i]++;
		integer n2=(integer)((iyintervalcount[i]-1)/2); // округление в меньшую сторону
		if (bgeomy) n2=ibalancen2(n2, qyS, qyN, ryboundary, iyintervalcount[i], i, iposmark); // балансировка 
		doublereal qn2S=qyS;
		for (integer i1=1; i1<n2; i1++) qn2S*=qyS; // возведение в степень.
		doublereal b1S=(ryboundary[i+1]-ryboundary[i])*(qyS-1.0)/(2.0*(qn2S-1.0));
		doublereal qn2N=qyN;
		for (integer i1=1; i1<(iyintervalcount[i]-n2); i1++) qn2N*=qyN; // возведение в степень.
		doublereal b1N=(ryboundary[i+1]-ryboundary[i])*(qyN-1.0)/(2.0*(qn2N-1.0));

        dy=(ryboundary[i+1]-ryboundary[i])/(iyintervalcount[i]-1);
		SetLength(ypos,iyoldsize,iposmark+iyintervalcount[i]);
        iyoldsize=iposmark+iyintervalcount[i];
        for (k=iposmark; k<=iposmark+iyintervalcount[i]-2; k++) 
		{
			if (!bgeomy) {
				if (lb == 1) {

					doublereal betavolkov = 1.0 / (sqrt(1.0 - deltavolkov));
					doublereal nvolkow = ((doublereal)(k - iposmark)) / ((doublereal)(iyintervalcount[i] - 1.0));
					ypos[k] = (ryboundary[1] - ryboundary[0]) * fmax(0.0, (((betavolkov + 1.0)*(pow(((betavolkov + 1.0) / (betavolkov - 1.0)), 2.0*nvolkow - 1.0)) - (betavolkov - 1.0))));
					ypos[k] /= 2.0*(1.0 + pow((betavolkov + 1.0) / (betavolkov - 1.0), 2.0*nvolkow - 1.0));
					ypos[k] += ryboundary[0];

				}
				else {
					// равномерная сетка
					ypos[k] = ryboundary[i] + (k - iposmark)*dy;
				}
			}
			else 
			{
				// неравномерная сетка
				// на основе геометрической прогрессии с двумя сгущениями.
				integer ic1=k-iposmark;
				doublereal madd;
				if (ic1<n2) {
					madd=b1S;
					for (integer i1=1; i1<ic1; i1++) madd*=qyS;
					if (ic1==n2-1) madd=0.5*(ryboundary[i+1]+ryboundary[i])-ypos[k-1];
				} else 
				{
					madd=b1N;
					for (integer i1=iyintervalcount[i]-1-ic1; i1>0; i1--) madd*=qyN;
				}
				if (k==iposmark) ypos[k]=ryboundary[i];
				else ypos[k]=ypos[k-1]+madd;
			}
		}
		iposmark=iposmark+iyintervalcount[i]-1;
	}
	SetLength(ypos,iyoldsize,iposmark+1);
	ypos[iposmark]=ryboundary[inumboundaryy];
	iny=iposmark;
    for (i=0; i<iny; i++) ypos[i]=ypos[i+1]; // сдвиг влево на 1
    SetLength(ypos,iny+1,iny); 
	iny--; // индексация начинается с нуля и заканчивается значением iny
    
    // 3 сентября 2017.
	// Добавление новых сеточных линий.
	// необходимость их добавления обусловлена требованнем АЛИС сеток.
	for (int i_28 = 0; i_28 <= inyadd; i_28++) {
		//SetLength(ypos, iny + 1, iny + 2);
		//ypos[iny + 1] = yposadd[i_28];
		//iny++;
		// Это добавление умное, препятствует возникновению одинаковых значений в массиве.
		addboundary(ypos, iny, yposadd[i_28],XZ_PLANE, b, lb, w, lw, s, ls);
	}
	Sort_method<doublereal>(ypos,iny);

    for (i=0; i<adapt_y; i++) simplecorrect_meshgen_y(ypos, iny, lb, ls, lw, b, s, w);


	inum_iter_ratio_good = 0;
	doublereal max_size_ratio_y = 1.0;
	while (1) {

		max_size_ratio_y = 1.0;
		for (i = 0; i<iny - 1; i++) {

			doublereal dmax = 0.0;
			doublereal dmin = 0.0;
			if (fabs(ypos[i + 1] - ypos[i])>fabs(ypos[i + 2] - ypos[i + 1])) {
				dmax = fabs(ypos[i + 1] - ypos[i]);
				dmin = fabs(ypos[i + 2] - ypos[i + 1]);
			}
			else {
				dmax = fabs(ypos[i + 2] - ypos[i + 1]);
				dmin = fabs(ypos[i + 1] - ypos[i]);
			}
			if (dmax / dmin>max_size_ratio_y) {
				max_size_ratio_y = dmax / dmin;
			}
		}
		
		if (max_size_ratio_y != max_size_ratio_y) {
			//printf("y axis max size ratio is equal = %1.4f\n", max_size_ratio_y);
			std::cout << "y axis max size ratio is equal = " << max_size_ratio_y << std::endl;
			system("PAUSE");
		}
		//getchar();

		if (max_size_ratio_y < (etalon_max_size_ratio*1.1)) {
			break;
		}

		// Корректировка max size ratio y axis.
		// сделано 9 апреля 2013 года.

		bool bplus = false;
		max_size_ratio_y = 1.0;
		for (i = 0; i<iny - 1; i++) {

			bplus = false;
			doublereal dmax = 0.0;
			doublereal dmin = 0.0;
			if (fabs(ypos[i + 1] - ypos[i])>fabs(ypos[i + 2] - ypos[i + 1])) {
				dmax = fabs(ypos[i + 1] - ypos[i]);
				dmin = fabs(ypos[i + 2] - ypos[i + 1]);
			}
			else {
				bplus = true;
				dmax = fabs(ypos[i + 2] - ypos[i + 1]);
				dmin = fabs(ypos[i + 1] - ypos[i]);
			}
			if (dmax / dmin>etalon_max_size_ratio) {
				doublereal pos_candidate;
				if (bplus) {
					pos_candidate = ypos[i + 1] + 0.5*dmax;
				}
				else {
					pos_candidate = ypos[i] + 0.5*dmax;
				}
				SetLength(ypos, iny + 1, iny + 2);
				ypos[iny + 1] = pos_candidate;
				iny = iny + 1;
				//BubbleEnhSort<doublereal>(ypos, 0, iny);
				Sort_method<doublereal>(ypos,iny);
				break;
			}
		}

		inum_iter_ratio_good++;

		//getchar();

	}

	//printf("y axis max size ratio is equal = %1.4f\n", max_size_ratio_y);
	std::cout << "y axis max size ratio is equal = " << max_size_ratio_y << std::endl;

#if doubleintprecision == 1
	printf("inum_iter_ratio_good is %lld\n", inum_iter_ratio_good);
#else
	printf("inum_iter_ratio_good is %d\n", inum_iter_ratio_good);
#endif
	
	if (bwait) {
		//getchar();
		system("pause");
	}

	/* // debug
	for (i=0; i<=iny; i++) {
        printf("%f  ",ypos[i]);   
	}//*/

    // По оси Oz
     
   
      

    int *izintervalcount; // число интервалов
    izintervalcount = new int [inumboundaryz]; // на один меньше чем число границ.
    for (i=0; i<(inumboundaryz); i++) {
         alphascal=(rzboundary[i+1]-rzboundary[i])/(rzboundary[inumboundaryz]-rzboundary[0]);
         inowintervalcount=(int)(alphascal*inz);
		 if (inowintervalcount < min_elem_in_z_element) inowintervalcount=min_elem_in_z_element;
         izintervalcount[i]=inowintervalcount;
	}
	
    // заполнение массива узлов
    iposmark = 1;
	doublereal dz;
    integer izoldsize=0;
    SetLength(zpos,izoldsize,1);
    izoldsize=1;
	for (i=0; i<(inumboundaryz); i++) 
	{
		// чтобы число узлов було всегда чётное
		if ((izintervalcount[i]-2)%2==0) izintervalcount[i]++;
		integer n2=(integer)((izintervalcount[i]-1)/2); // округление в меньшую сторону 
		if (bgeomz) n2=ibalancen2(n2, qzB, qzT, rzboundary, izintervalcount[i], i, iposmark); // балансировка 
#if doubleintprecision == 1
		//printf("n2=%lld\n",n2); // debug
#else
		//printf("n2=%d\n",n2); // debug
#endif
		
		doublereal qn2B=qzB; // нижняя сторона
		// точного значения в 0.25 всё равно достигнуть не получится из-за возникающих ошибок округления.
		for (int i1=1; i1<n2; i1++) qn2B*=qzB; // возведение в степень.
		doublereal b1B=(rzboundary[i+1]-rzboundary[i])*(qzB-1.0)/(2.0*(qn2B-1.0));
		doublereal qn2T=qzT; // верхняя сторона
		for (int i1=1; i1<(izintervalcount[i]-n2); i1++) qn2T*=qzT; // возведение в степень.
		doublereal b1T=(rzboundary[i+1]-rzboundary[i])*(qzT-1.0)/(2.0*(qn2T-1.0)); 

        dz=(rzboundary[i+1]-rzboundary[i])/(izintervalcount[i]-1);
		SetLength(zpos,izoldsize,iposmark+izintervalcount[i]);
        izoldsize=iposmark+izintervalcount[i];
        for (k=iposmark; k<=iposmark+izintervalcount[i]-2; k++)
		{
			if (!bgeomz) {

				if (lb == 1) {

					doublereal betavolkov = 1.0 / (sqrt(1.0 - deltavolkov));
					doublereal nvolkow = ((doublereal)(k - iposmark)) / ((doublereal)(izintervalcount[i] - 1.0));
					zpos[k] = (rzboundary[1] - rzboundary[0]) * fmax(0.0, (((betavolkov + 1.0)*(pow(((betavolkov + 1.0) / (betavolkov - 1.0)), 2.0*nvolkow - 1.0)) - (betavolkov - 1.0))));
					zpos[k] /= 2.0*(1.0 + pow((betavolkov + 1.0) / (betavolkov - 1.0), 2.0*nvolkow - 1.0));
					zpos[k] += rzboundary[0];
					//printf("%e\n", zpos[k]);
					//getchar();

				}
				else {
					// равномерная сетка
					zpos[k] = rzboundary[i] + (k - iposmark)*dz;
				}
			}
			else 
			{
				// неравномерная сетка
				// на основе геометрической прогрессии с двумя сгущениями.
				int ic1=k-iposmark;
				doublereal madd;
				if (ic1<n2) {
					madd=b1B;
					for (integer i1=0; i1<ic1; i1++) madd*=qzB; 
					if (ic1==n2-1) madd=0.5*(rzboundary[i+1]+rzboundary[i])-zpos[k-1]; // чтобы избежать ошибок округления.
					//printf("first.:");
				} else 
				{
					madd=b1T;
					for (int i1=izintervalcount[i]-1-ic1; i1>0; i1--) madd*=qzT;
					// это тока предпоследний.
					//if (k==iposmark+izintervalcount[i]-2) madd*=rzboundary[i+1]-zpos[k-1]; // чтобы избежать ошибок округления.
					//printf("second.:");
				}
				if (k==iposmark) zpos[k]=rzboundary[i]; // чтобы избежать ошибок округления.
				else {
					// учитывается что это начиная со второго до предпоследнего
					zpos[k]=zpos[k-1]+madd;
				}
#if doubleintprecision == 1
		//printf("zpos[%lld]=%f, madd=%f\n",k,zpos[k],madd);
#else
		//printf("zpos[%d]=%f, madd=%f\n",k,zpos[k],madd);
#endif
				
				//getchar();
			}
		}
		iposmark=iposmark+izintervalcount[i]-1;
	}
	SetLength(zpos,izoldsize,iposmark+1);
	zpos[iposmark]=rzboundary[inumboundaryz];
#if doubleintprecision == 1
	//printf("finish zpos[%lld]=%f\n",iposmark,zpos[iposmark]);
#else
	//printf("finish zpos[%d]=%f\n",iposmark,zpos[iposmark]);
#endif
	
	//getchar();
	inz=iposmark;
	for (i=0; i<inz; i++) zpos[i]=zpos[i+1]; // сдвиг влево на 1
    SetLength(zpos,inz+1,inz); 
	inz--; // индексация начинается с нуля и заканчивается значением inz

	// 3 сентября 2017.
	// Добавление новых сеточных линий.
	// необходимость их добавления обусловлена требованнем АЛИС сеток.
	for (int i_28 = 0; i_28 <= inzadd; i_28++) {
		//SetLength(zpos, inz + 1, inz + 2);
		//zpos[inz + 1] = zposadd[i_28];
		//inz++;
		// Это добавление умное, препятствует возникновению одинаковых значений в массиве.
		addboundary(zpos, inz, zposadd[i_28],XY_PLANE, b, lb, w, lw, s, ls);
	}
	Sort_method<doublereal>(zpos,inz);

	for (i=0; i<adapt_z; i++) simplecorrect_meshgen_z(zpos, inz, lb, ls, lw, b, s, w);
	
	inum_iter_ratio_good = 0;
	doublereal max_size_ratio_z = 1.0;
	while (1) {

		max_size_ratio_z = 1.0;
		for (i = 0; i<inz - 1; i++) {

			doublereal dmax = 0.0;
			doublereal dmin = 0.0;
			if (fabs(zpos[i + 1] - zpos[i])>fabs(zpos[i + 2] - zpos[i + 1])) {
				dmax = fabs(zpos[i + 1] - zpos[i]);
				dmin = fabs(zpos[i + 2] - zpos[i + 1]);
			}
			else {
				dmax = fabs(zpos[i + 2] - zpos[i + 1]);
				dmin = fabs(zpos[i + 1] - zpos[i]);
			}
			if (dmax / dmin>max_size_ratio_z) {
				max_size_ratio_z = dmax / dmin;
			}
		}
		
		if (max_size_ratio_z != max_size_ratio_z) {
			//printf("z axis max size ratio is equal = %1.4f\n", max_size_ratio_z);
			std::cout << "z axis max size ratio is equal = " << max_size_ratio_z << std::endl;
			system("PAUSE");
		}
		//getchar();

		if (max_size_ratio_z < (etalon_max_size_ratio*1.1)) {
			break;
		}

		// Корректировка max size ratio z axis.
		// сделано 9 апреля 2013 года.

		bool bplus = false;
		max_size_ratio_z = 1.0;
		for (i = 0; i<inz - 1; i++) {

			bplus = false;
			doublereal dmax = 0.0;
			doublereal dmin = 0.0;
			if (fabs(zpos[i + 1] - zpos[i])>fabs(zpos[i + 2] - zpos[i + 1])) {
				dmax = fabs(zpos[i + 1] - zpos[i]);
				dmin = fabs(zpos[i + 2] - zpos[i + 1]);
			}
			else {
				bplus = true;
				dmax = fabs(zpos[i + 2] - zpos[i + 1]);
				dmin = fabs(zpos[i + 1] - zpos[i]);
			}
			if (dmax / dmin>etalon_max_size_ratio) {
				doublereal pos_candidate;
				if (bplus) {
					pos_candidate = zpos[i + 1] + 0.5*dmax;
				}
				else {
					pos_candidate = zpos[i] + 0.5*dmax;
				}
				SetLength(zpos, inz + 1, inz + 2);
				zpos[inz + 1] = pos_candidate;
				inz = inz + 1;
				//BubbleEnhSort<doublereal>(zpos, 0, inz);
				Sort_method<doublereal>(zpos,inz);
				break;
			}
		}

		inum_iter_ratio_good++;

		//getchar();

	}

	//printf("z axis max size ratio is equal = %1.4f\n", max_size_ratio_z);
	std::cout << "z axis max size ratio is equal = " << max_size_ratio_z << std::endl;

#if doubleintprecision == 1
	printf("inum_iter_ratio_good is %lld\n", inum_iter_ratio_good);
#else
	printf("inum_iter_ratio_good is %d\n", inum_iter_ratio_good);
#endif
	
	if (bwait) {
		// getchar();
		system("pause");

	}

	//for (i = 0; i <= inz; i++) {
	//printf("%e  ", zpos[i]);
	//}
	//getchar();

	integer inz_fix = inz;
	if (b_adhesion_Mesh) {
		// Корректировка источника XY
		for (i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranXY[i]) {
				doublereal xc = 0.5 * (s[i].g.xS + s[i].g.xE);
				doublereal yc = 0.5 * (s[i].g.yS + s[i].g.yE);
				doublereal zg = s[i].g.zS;
				// найти позицию +Z на сетке zpos
				// найти центр кО.
				// Если он принадлежит Solid block то сместить истоник в цент этой клетки.
				// Адаптацию maxsize ratio 2 не делать.
				// Если неуспех ищемцентр КО по -Z и смещаем туда.
				integer i55_found = -2;
				for (integer i55 = 0; i55 <= inz_fix; i55++) {
					if (fabs(zpos[i55] - zg) < 1.0e-36) {
						i55_found = i55;
						break;
					}
				}
				if (i55_found >= 0) {
					if (i55_found < inz_fix) {
						doublereal zg1 = 0.5 * (zg + zpos[i55_found + 1]);
						printf("zg1=%e\n", zg1);
						integer i56_found = -2;
						for (integer ib55 = 0; ib55 < lb; ib55++) {
							if ((xc > b[ib55].g.xS) && (xc < b[ib55].g.xE) && (yc > b[ib55].g.yS) && (yc < b[ib55].g.yE) && (zg1 > b[ib55].g.zS) && (zg1 < b[ib55].g.zE))
							{
								i56_found = ib55;
							}
						}

						bool bzero_pos = true;
						doublereal zg2 = zpos[0] - 0.5 * fabs(zpos[1] - zpos[0]);
						if (i55_found > 0) {
							zg2 = 0.5 * (zg + zpos[i55_found - 1]);
							bzero_pos = false;
						}

						printf("zg2=%e\n", zg2);
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (zg2 > b[ib57].g.zS) && (zg2 < b[ib57].g.zE))
							{
								i57_found = ib57;
							}
						}

						if (i56_found >= 0) {
							if (b[i56_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								if ((i57_found >= 0) && (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
									// Мы помещаем источник тепла в блок с большей теплопроводностью.
									// comparison_lam выдаёт истину если теплопроводность блока i56 больше.
									if (comparison_lam(matlist, b, i56_found, i57_found, 25.0)) {
										// в блоке i56 теплопроводность выше.
										s[i].g.zS = zg1;
										s[i].g.zE = zg1;
										addboundary(zpos, inz, zg1, XY_PLANE, b, lb, w, lw, s, ls);
									}
									else {
										// в блоке i57 теплопроводность выше.
										s[i].g.zS = zg2;
										s[i].g.zE = zg2;
										addboundary(zpos, inz, zg2, XY_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									// Найден Солид Блок.
									printf("zg1==%e\n", zg1);
									doublereal zgolg = zpos[1];
									if (inz == 1) {
										// Слишком малоразмерная расчётная сетка.
										zgolg = 0.5 * (zpos[0] + zg1);
										addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
										zgolg = 0.5 * (zpos[1] + zg1);
										addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
										zgolg = 0.5 * (zpos[0] + zg1);
										// Расстояние от края вглубь расчётной области две клетки.
										s[i].g.zS = zgolg;
										s[i].g.zE = zgolg;
									}
									else {
										if (bzero_pos) {
											// Слишком малоразмерная расчётная сетка.
											zgolg = 0.5 * (zpos[0] + zg1);
											addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
											zgolg = 0.5 * (zpos[1] + zg1);
											addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
											zgolg = 0.5 * (zpos[0] + zg1);
											// Расстояние от края вглубь расчётной области две клетки.
											s[i].g.zS = zgolg;
											s[i].g.zE = zgolg;
										}
										else {
											s[i].g.zS = zg1;
											s[i].g.zE = zg1;
										}
									}
									addboundary(zpos, inz, zg1, XY_PLANE, b, lb, w, lw, s, ls);

								}
							}
							else {

								if (i57_found >= 0) {
									if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
										// Найден Солид Блок.
										s[i].g.zS = zg2;
										s[i].g.zE = zg2;
										addboundary(zpos, inz, zg2, XY_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
									system("PAUSE");
									exit(1);
								}
							}
						}
					}
					else {
						doublereal zg2 = 0.5 * (zg + zpos[i55_found - 1]);
						printf("zg2=%e\n", zg2);
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (zg2 > b[ib57].g.zS) && (zg2 < b[ib57].g.zE))
							{
								i57_found = ib57;
							}
						}
						if (i57_found >= 0) {
							if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								// Найден Солид Блок.
								s[i].g.zS = zg2;
								s[i].g.zE = zg2;
								addboundary(zpos, inz, zg2, XY_PLANE, b, lb, w, lw, s, ls);
							}
						}
						else {
							printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
							system("PAUSE");
							exit(1);
						}
					}
				}
			}
		}
	}
	//BubbleEnhSort<doublereal>(zpos, 0, inz);
	Sort_method<doublereal>(zpos,inz);
	//SetLength(zpos, inz + 2, inz);
	//inz++;

	integer inx_fix = inx;
	if (b_adhesion_Mesh) {
		// Корректировка источника YZ
		for (i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranYZ[i]) {
				doublereal zc = 0.5 * (s[i].g.zS + s[i].g.zE);
				doublereal yc = 0.5 * (s[i].g.yS + s[i].g.yE);
				doublereal xg = s[i].g.xS;
				// найти позицию +Z на сетке zpos
				// найти центр кО.
				// Если он принадлежит Solid block то сместить истоник в цент этой клетки.
				// Адаптацию maxsize ratio 2 не делать.
				// Если неуспех ищемцентр КО по -Z и смещаем туда.
				integer i55_found = -2;
				for (integer i55 = 0; i55 <= inx_fix; i55++) {
					if (fabs(xpos[i55] - xg) < 1.0e-36) {
						i55_found = i55;
						break;
					}
				}
				if (i55_found >= 0) {
					if (i55_found < inx_fix) {
						doublereal xg1 = 0.5 * (xg + xpos[i55_found + 1]);
						printf("xg1=%e\n", xg1);
						integer i56_found = -2;
						for (integer ib55 = 0; ib55 < lb; ib55++) {
							if ((zc > b[ib55].g.zS) && (zc < b[ib55].g.zE) && (yc > b[ib55].g.yS) && (yc < b[ib55].g.yE) && (xg1 > b[ib55].g.xS) && (xg1 < b[ib55].g.xE))
							{
								i56_found = ib55;
							}
						}

						bool bzero_pos = true;
						doublereal xg2 = xpos[0] - 0.5 * fabs(xpos[1] - xpos[0]);
						if (i55_found > 0) {
							xg2 = 0.5 * (xg + xpos[i55_found - 1]);
							bzero_pos = false;
						}


						printf("xg2=%e\n", xg2);
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (xg2 > b[ib57].g.xS) && (xg2 < b[ib57].g.xE))
							{
								i57_found = ib57;
							}
						}


						if (i56_found >= 0) {
							if (b[i56_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								// TODO 11.07.2016
								if ((i57_found >= 0) && (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
									// Мы помещаем источник тепла в блок с большей теплопроводностью.
									// comparison_lam выдаёт истину если теплопроводность блока i56 больше.
									if (comparison_lam(matlist, b, i56_found, i57_found, 25.0)) {
										// в блоке i56 теплопроводность выше.
										s[i].g.xS = xg1;
										s[i].g.xE = xg1;
										addboundary(xpos, inz, xg1, YZ_PLANE, b, lb, w, lw, s, ls);
									}
									else {
										// в блоке i57 теплопроводность выше.
										s[i].g.xS = xg2;
										s[i].g.xE = xg2;
										addboundary(xpos, inz, xg2, YZ_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									// Найден Солид Блок.
									printf("xg1==%e\n", xg1);
									doublereal xgolg = xpos[1];
									if (inx == 1) {
										// Слишком малоразмерная расчётная сетка.
										xgolg = 0.5 * (xpos[0] + xg1);
										addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
										xgolg = 0.5 * (xpos[1] + xg1);
										addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
										xgolg = 0.5 * (xpos[0] + xg1);
										// Расстояние от края вглубь расчётной области две клетки.
										s[i].g.xS = xgolg;
										s[i].g.xE = xgolg;
									}
									else {
										if (bzero_pos) {
											// Слишком малоразмерная расчётная сетка.
											xgolg = 0.5 * (xpos[0] + xg1);
											addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
											xgolg = 0.5 * (xpos[1] + xg1);
											addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
											xgolg = 0.5 * (xpos[0] + xg1);
											// Расстояние от края вглубь расчётной области две клетки.
											s[i].g.xS = xgolg;
											s[i].g.xE = xgolg;
										}
										else {
											s[i].g.xS = xg1;
											s[i].g.xE = xg1;
										}
									}
									addboundary(xpos, inx, xg1, YZ_PLANE, b, lb, w, lw, s, ls);


								}
							}
							else {

								if (i57_found >= 0) {
									if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
										// Найден Солид Блок.
										s[i].g.xS = xg2;
										s[i].g.xE = xg2;
										addboundary(xpos, inx, xg2, YZ_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
									system("PAUSE");
									exit(1);
								}
							}
						}
					}
					else {
						doublereal xg2 = 0.5 * (xg + xpos[i55_found - 1]);
						printf("xg2=%e\n", xg2);
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (xg2 > b[ib57].g.xS) && (xg2 < b[ib57].g.xE))
							{
								i57_found = ib57;
							}
						}
						if (i57_found >= 0) {
							if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								// Найден Солид Блок.
								s[i].g.xS = xg2;
								s[i].g.xE = xg2;
								addboundary(xpos, inx, xg2, YZ_PLANE, b, lb, w, lw, s, ls);
							}
						}
						else {
							printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
							system("PAUSE");
							exit(1);
						}
					}
				}
			}
		}
	}
	//BubbleEnhSort<doublereal>(xpos, 0, inx);
	Sort_method<doublereal>(xpos,inx);

	integer iny_fix = iny;
	if (b_adhesion_Mesh) {
		// Корректировка источника XZ
		for (i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranXZ[i]) {
				doublereal xc = 0.5 * (s[i].g.xS + s[i].g.xE);
				doublereal zc = 0.5 * (s[i].g.zS + s[i].g.zE);
				doublereal yg = s[i].g.yS;
				// найти позицию +Z на сетке zpos
				// найти центр кО.
				// Если он принадлежит Solid block то сместить истоник в цент этой клетки.
				// Адаптацию maxsize ratio 2 не делать.
				// Если неуспех ищемцентр КО по -Z и смещаем туда.
				integer i55_found = -2;
				for (integer i55 = 0; i55 <= iny_fix; i55++) {
					if (fabs(ypos[i55] - yg) < 1.0e-36) {
						i55_found = i55;
						break;
					}
				}
				if (i55_found >= 0) {
					if (i55_found < iny_fix) {
						doublereal yg1 = 0.5 * (yg + ypos[i55_found + 1]);
						printf("yg1=%e\n", yg1);
						integer i56_found = -2;
						for (integer ib55 = 0; ib55 < lb; ib55++) {
							if ((xc > b[ib55].g.xS) && (xc < b[ib55].g.xE) && (zc > b[ib55].g.zS) && (zc < b[ib55].g.zE) && (yg1 > b[ib55].g.yS) && (yg1 < b[ib55].g.yE))
							{
								i56_found = ib55;
							}
						}
						doublereal yg2 = 0.5 * (yg + ypos[i55_found - 1]);
						printf("yg2=%e\n", yg2);
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yg2 > b[ib57].g.yS) && (yg2 < b[ib57].g.yE))
							{
								i57_found = ib57;
							}
						}


						if (i56_found >= 0) {
							if (b[i56_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {

								if ((i57_found >= 0) && (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
									// Мы помещаем источник тепла в блок с большей теплопроводностью.
									// comparison_lam выдаёт истину если теплопроводность блока i56 больше.
									if (comparison_lam(matlist, b, i56_found, i57_found, 25.0)) {
										// в блоке i56 теплопроводность выше.
										s[i].g.yS = yg1;
										s[i].g.yE = yg1;
										addboundary(ypos, iny, yg1, XZ_PLANE, b, lb, w, lw, s, ls);
									}
									else {
										// в блоке i57 теплопроводность выше.
										s[i].g.yS = yg2;
										s[i].g.yE = yg2;
										addboundary(ypos, iny, yg2, XZ_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									// Найден Солид Блок.
									s[i].g.yS = yg1;
									s[i].g.yE = yg1;
									addboundary(ypos, iny, yg1, XZ_PLANE, b, lb, w, lw, s, ls);
								}
							}
							else {

								if (i57_found >= 0) {
									if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
										// Найден Солид Блок.
										s[i].g.yS = yg2;
										s[i].g.yE = yg2;
										addboundary(ypos, iny, yg2, XZ_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {
									printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
									system("PAUSE");
									exit(1);
								}
							}
						}
					}
					else {
						doublereal yg2 = 0.5 * (yg + ypos[i55_found - 1]);
						printf("yg2=%e\n", yg2);
						integer i57_found = -2;
						for (integer ib57 = 0; ib57 < lb; ib57++) {
							if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yg2 > b[ib57].g.yS) && (yg2 < b[ib57].g.yE))
							{
								i57_found = ib57;
							}
						}
						if (i57_found >= 0) {
							if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
								// Найден Солид Блок.
								s[i].g.yS = yg2;
								s[i].g.yE = yg2;
								addboundary(ypos, iny, yg2, XZ_PLANE, b, lb, w, lw, s, ls);
							}
						}
						else {
							printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
							system("PAUSE");
							exit(1);
						}
					}
				}
			}
		}
	}
	//BubbleEnhSort<doublereal>(ypos, 0, iny);
	Sort_method<doublereal>(ypos,iny);

	
	// улучшение качества сетки по второму критерию с значением 30.0
	// как во flowvision.
	//quolite_refinement(inx, iny, inz, xpos, ypos, zpos);

	/*// debug
	for (i=0; i<=inz; i++) {
	printf("%f  ",zpos[i]);
	}
	getchar();//*/

	// debug
    /*
	for (i=0; i<=inz; i++) {
        printf("%e  ",zpos[i]);   
	}
	getchar();//
    */

    // Освобождение оперативной памяти.
	delete[] rxboundary;
	delete[] ryboundary;
	delete[] rzboundary;
	delete[] ixintervalcount;
	delete[] iyintervalcount;
    delete[] izintervalcount;
	delete[] source_indexpopadaniqnagranXY;
	delete[] source_indexpopadaniqnagranYZ;
	delete[] source_indexpopadaniqnagranXZ;

} // unevensimplemeshgen



/* генерирует простейшую расчётную сетку hex cartesian -
* структурированная прямоугольная сетка в трёхмерной расчётной
* области в виде кирпича. Особенность данного генератора в том
* что он может генерировать как равномерную сетку так и неравномерную.
* Данный автоматический сеточный генератор содержит больше управляющих параметров чем
* автоматический сеточный генератор simplemeshgen.
* дата реализации 8 февраля 2012.
* Изучение сеточного генератора ANSYS Icepak а также требование чрезвычайно большеразмерных моделей 
* выдвинуло требование производить расчёт на самой грубой из всех возможных сеток по причине быстродействия.
* Часто нельзя ждать пока посчитается большеразмерная модель на подробной сетке (счёт длится несколько суток)
* а требуется быстро оценить (максимум один час) качественное поведение большеразмерной модели.
* 3 мая 2016.
*/
// 25.04.2018 готов к АССЕМБЛЕСАМ.
void coarsemeshgen2(doublereal* &xpos, doublereal* &ypos, doublereal* &zpos, int &inx, int &iny, int &inz,
	int lb, int ls, int lw, BLOCK* &b, SOURCE* &s, WALL* &w, int lu, UNION* &my_union, TPROP* matlist,
	doublereal* &xposadd, doublereal* &yposadd, doublereal* &zposadd,
	int &inxadd, int &inyadd, int &inzadd, integer &iunion_id_p1)
{

	//23.12.2018
	bool brefinement = false;// (AMG1R6_LABEL == 1); // Сгущение расчетной сетки на источниках тепла.

	//*****************************************************************************************************************
	// Изменяемые параметры автоматического сеточного генератора

	// если bgeom  , то используется неравномерная сетка по закону геометрической прогрессии.
	bool bgeomx = false; // по оси Ox
	bool bgeomy = false; // по оси Oy
	bool bgeomz = false; // по оси Oz
	// 1.05 - очень слабо неравномерная, 1.1, 1.2 - умеренно неравномерная, 1.25, 1.5, 2 - сильно неравномерная.
	// для низкорейнольдсовых моделей турбулентности знаменатель геометрической прогрессии должен быть не больше 1.3.
	//doublereal qxW = 10.0; // параметр сгущения к левой West стороне. 
	//doublereal qxE = 10.0; // параметр сгущения к правой East стороне
	//doublereal qyS = 10.0; // South
	//doublereal qyN = 10.0; // North
	//doublereal qzB = 10.0; // Bottom
	//doublereal qzT = 10.0; // Top
	//doublereal qzB=1.00005; // Bottom
	//doublereal qzT=1.35; // Top

	//bool bsnap_TO = bsnap_TO_global; // snap to grid (избавляет от щелей в сложных моделях).
	//doublereal snap_to_multiplyer = 0.3;// константа лежит в диапазоне (0..1) и регулирует какой зазор игнорировать.
	//*****************************************************************************************************************

	//bool brepeat = true;

	
	// По Оси Ox
	doublereal *rxboundary=nullptr; // массив обязательных границ
	int inumboundaryx = 1;
	
	doublereal *ryboundary = nullptr; // массив обязательных границ
	int inumboundaryy = 1;

	doublereal *rzboundary = nullptr; // массив обязательных границ
	int inumboundaryz = 1;

	// подготовка к вычислению зазоров minimum fluid gap.
	calc_minimum_fluid_gap1(inumboundaryx, rxboundary, inumboundaryy, ryboundary, inumboundaryz, rzboundary,
		lb, ls, lw, b, s, w, lu, my_union, iunion_id_p1);
	//printf("%d %d %d", inumboundaryx, inumboundaryy, inumboundaryz);

	doublereal minimum_fluid_gap_x = 1.0e36;
	doublereal minimum_fluid_gap_y = 1.0e36;
	doublereal minimum_fluid_gap_z = 1.0e36;

	

	// Непосредтвенное вычисление зазоров minimum fluid gap.
	calc_minimum_fluid_gap2(inumboundaryx, rxboundary, inumboundaryy, ryboundary,
		inumboundaryz, rzboundary, minimum_fluid_gap_x, minimum_fluid_gap_y, minimum_fluid_gap_z,
		lb, ls, lw, b, s, w, lu, my_union, iunion_id_p1);
	

	bool *source_indexpopadaniqnagranYZ = nullptr;
	bool *source_indexpopadaniqnagranXY = nullptr;
	bool *source_indexpopadaniqnagranXZ = nullptr;

	// 0.	none
	// 1.	Snap to grid
	// 2.	Snap to grid ALICE
	// 3.	Snap to grid ++
	

	// 12.03.2017
	// реализация snap to
	// уменьшает размерность сеточной модели на 33%.
	if ((bsnap_TO_global==1)||(bsnap_TO_global==3)) {
		if (b_adhesion_Mesh) {
			snap_to_moving(source_indexpopadaniqnagranYZ,
				source_indexpopadaniqnagranXY,
				source_indexpopadaniqnagranXZ,
				rxboundary, ryboundary, rzboundary,
				inumboundaryx, inumboundaryy, inumboundaryz,
				minimum_fluid_gap_x, minimum_fluid_gap_y, minimum_fluid_gap_z,
				lb, ls, lw, b, s, w, lu, my_union, iunion_id_p1);
		}
	}
	//printf("%d %d %d", inumboundaryx, inumboundaryy, inumboundaryz);
	//getchar();

	// улучшает качество аппроксимации сеткой цилиндров вводя дополнительно для
	// каждого цилиндра четыре сеточных линии на расстоянии шага сетки от вшешнего
	// контура основания цилиндра. Это требуется для сильно неравномерных сеток.
	calc_minimum_fluid_gap3(inumboundaryx, rxboundary, inumboundaryy, ryboundary, inumboundaryz, rzboundary,
		lb, ls, lw, b, s, w, lu, my_union, iunion_id_p1);

	//integer i;


	// preprocessing
	integer* ib_marker = new integer[inumboundaryx*inumboundaryy*inumboundaryz];
	bool* ib_marker_flag_power = new bool[inumboundaryx*inumboundaryy*inumboundaryz];
	bool* ib_marker_flag_fluid = new bool[inumboundaryx*inumboundaryy*inumboundaryz];
	/*
	for (integer i9 = 0; i9 < inumboundaryx; i9++) {
	for (integer i7 = 0; i7 < inumboundaryy; i7++) {
	for (integer i8 = 0; i8 < inumboundaryz; i8++) {
	doublereal yc = 0.5*(ryboundary[i7] + ryboundary[i7 + 1]);
	doublereal zc = 0.5*(rzboundary[i8] + rzboundary[i8 + 1]);
	doublereal xc = 0.5*(rxboundary[i9] + rxboundary[i9 + 1]);
	int ib;
	TOCHKA p;
	p.x = xc;
	p.y = yc;
	p.z = zc;
	in_model_fluid_gap(p, ib, b, lb);
	// x + y*dimx+z*dimx*dimy
	ib_marker[i9+ inumboundaryx*i7+ inumboundaryx*inumboundaryy*i8] = ib;
	}
	}
	}
	*/

	// Быстрый препроцессинг за линейное время.
	Block_indexes* block_indexes = new Block_indexes[lb];
	if (block_indexes == nullptr) {
		printf("error in allocation memory for block_indexes in enumerate_volume_improved.\n");
		system("pause");
		exit(1);
	}
	//integer  j, k;


	// 08.04.2018
	
#pragma omp parallel for
	for (integer i = 0; i < lb; i++) {
		// инициализация, на случай если блоки не будут распознаны.
		block_indexes[i].iL = -1;
		block_indexes[i].iR = -2;
		block_indexes[i].jL = -1;
		block_indexes[i].jR = -2;
		block_indexes[i].kL = -1;
		block_indexes[i].kR = -2;
	}
	
/*
	for (integer i = 0; i < lb; i++) {
		doublereal x4 = b[i].g.xS;
		doublereal distmax;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryx; j++) {
			if (fabs(rxboundary[j] - x4) < distmax) {
				block_indexes[i].iL = j;
				distmax = fabs(rxboundary[j] - x4);
			}
		}
		x4 = b[i].g.xE;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryx; j++) {
			if (fabs(rxboundary[j] - x4) < distmax) {
				block_indexes[i].iR = j;
				distmax = fabs(rxboundary[j] - x4);
			}
		}
		x4 = b[i].g.yS;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryy; j++) {
			if (fabs(ryboundary[j] - x4) < distmax) {
				block_indexes[i].jL = j;
				distmax = fabs(ryboundary[j] - x4);
			}
		}
		x4 = b[i].g.yE;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryy; j++) {
			if (fabs(ryboundary[j] - x4) < distmax) {
				block_indexes[i].jR = j;
				distmax = fabs(ryboundary[j] - x4);
			}
		}
		x4 = b[i].g.zS;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryz; j++) {
			if (fabs(rzboundary[j] - x4) < distmax) {
				block_indexes[i].kL = j;
				distmax = fabs(rzboundary[j] - x4);
			}
		}
		x4 = b[i].g.zE;
		distmax = 1.0e36;
		for (integer j = 0; j <= inumboundaryz; j++) {
			if (fabs(rzboundary[j] - x4) < distmax) {
				block_indexes[i].kR = j;
				distmax = fabs(rzboundary[j] - x4);
			}
		}
	}
	*/
	
#pragma omp parallel for
	for (integer i = 0; i < lb; i++) {
		//if (b[i].iunion_id == iunion_id_p1) {
		{
			

			doublereal x4 = b[i].g.xS;
			for (integer j = 0; j <= inumboundaryx; j++) {
				if (fabs(rxboundary[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].iL = j;
					break;
				}
			}
			x4 = b[i].g.xE;
			for (integer j = 0; j <= inumboundaryx; j++) {
				if (fabs(rxboundary[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].iR = j;
					break;
				}
			}
			x4 = b[i].g.yS;
			for (integer j = 0; j <= inumboundaryy; j++) {
				if (fabs(ryboundary[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].jL = j;
					break;
				}
			}
			x4 = b[i].g.yE;
			for (integer j = 0; j <= inumboundaryy; j++) {
				if (fabs(ryboundary[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].jR = j;
					break;
				}
			}
			x4 = b[i].g.zS;
			for (integer j = 0; j <= inumboundaryz; j++) {
				if (fabs(rzboundary[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].kL = j;
					break;
				}
			}
			x4 = b[i].g.zE;
			for (integer j = 0; j <= inumboundaryz; j++) {
				if (fabs(rzboundary[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i].kR = j;
					break;
				}
			}
		}
	}

	

	integer ismarker = 0;
	if (iunion_id_p1 > 0) {
#pragma omp parallel for reduction(+ : ismarker)
		for (integer m7 = 0; m7 < lb; m7++) {
			if (block_indexes[m7].iL < 0) {
				block_indexes[m7].iL = 0;
				ismarker++;
			}
			if (block_indexes[m7].jL < 0) {
				block_indexes[m7].jL = 0;
				ismarker++;
			}
			if (block_indexes[m7].kL < 0) {
				block_indexes[m7].kL = 0;
				ismarker++;
			}
			if (block_indexes[m7].iR < 0) {
				block_indexes[m7].iR = inumboundaryx;
				ismarker++;
			}
			if (block_indexes[m7].jR < 0) {
				block_indexes[m7].jR = inumboundaryy;
				ismarker++;
			}
			if (block_indexes[m7].kR < 0) {
				block_indexes[m7].kR = inumboundaryz;
				ismarker++;
			}
		}
	}
	
	//printf("ismarker =%d\n", ismarker);
	//getchar();
	

	// Количество проходов существенно сократилось и в итоге это приводит к существенному
	// увеличению быстродействия.
	integer m7;
	integer ib_stub = -1;
	// Мы найдем самый большой по размеру Hollow block, иначе будет просто кабинет.
	ib_stub = 0;
	doublereal vol_stub = -1.0;
	for (integer i = 0; i < lb; i++) {
		//if (b[i].iunion_id == iunion_id_p1) {
		{
			if (b[i].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
				if (fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS) > vol_stub) {
					ib_stub = i;
					vol_stub = fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS);
				}
			}
		}
	}

	integer mult_xy = inumboundaryx*inumboundaryy;
#pragma omp parallel for
	for (integer k1 = 0; k1 < inumboundaryz; k1++) {
		integer ik_1 = mult_xy*k1;
		for (integer j1 = 0; j1 < inumboundaryy; j1++) {
			integer ij_1 = inumboundaryx*j1 + ik_1;
			for (integer i1 = 0; i1 < inumboundaryx; i1++) {
				ib_marker[i1 + ij_1] = ib_stub;//-1
				ib_marker_flag_power[i1 + ij_1] = false;
				ib_marker_flag_fluid[i1 + ij_1] = false;
			}
		}
	}
	for (m7 = 0; m7 < lb; m7++) {

		bool bpowerON_loc = false;
		if ((b[m7].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (fabs(b[m7].arr_Sc[0]) > 0.0)) {
			// Обнаружена неудовлетворительная точность на источниках тепла,
			// Около них сетку необходимо сгущать.
			bpowerON_loc = true;
		}
		bool bFluidON = false;
		if (b[m7].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
			bFluidON = true;
		}

#pragma omp parallel for
		for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++) {
			integer ik_1 = mult_xy*k1;
			for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++) {
				integer ij_1 = inumboundaryx*j1 + ik_1;
				for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++) {
					ib_marker[i1 + ij_1] = m7;
					ib_marker_flag_power[i1 + ij_1] = bpowerON_loc;
					ib_marker_flag_fluid[i1 + ij_1] = bFluidON;
				}
			}
		}
	}

	delete[] block_indexes;
	//printf("identifire blocks number 80 procent.\n");


	int *ixintervalcount; // число интервалов
	ixintervalcount = new int[inumboundaryx]; // на один меньше чем число границ.
	//doublereal alphascal;
	//integer inowintervalcount;
//	for (i = 0; i<(inumboundaryx); i++) {
	//	alphascal = (rxboundary[i + 1] - rxboundary[i]) / (rxboundary[inumboundaryx] - rxboundary[0]);
		//inowintervalcount = (integer)(alphascal*inx);
	//	if (inowintervalcount < min_elem_in_x_element) inowintervalcount = min_elem_in_x_element;
		//ixintervalcount[i] = inowintervalcount;
	//}

#pragma omp parallel for
	for (integer i = 0; i < (inumboundaryx); i++) {

		bool bpowerON = false;

		// Если мы в solide то одна клетка.
		// Если во Fluide то две клетки.
		doublereal cpos = 0.5 * (rxboundary[i + 1] + rxboundary[i]);

		// В плоскости найден одиночный fluid контрольный объем имющий 
		// по E и W сторонам только Solid и hollow соседей.
		bool bfound_onex_fluid_cv = false;
		// применим нисходящее проектирование сверху вниз.
		//integer ibcur = 0; // номер текущего блока (кабинет по умолчанию).
		// сканируем плоскость.
		for (integer iy = 0; (bfound_onex_fluid_cv==false)&&(iy < (inumboundaryy)); iy++) {
			for (integer iz = 0; (bfound_onex_fluid_cv == false)&&(iz < (inumboundaryz)); iz++) {
				

				integer iP = i + inumboundaryx*iy + inumboundaryx*inumboundaryy*iz;

				//ibcur = ib_marker[iP];
				
				if (ib_marker_flag_power[iP]) {
					// Обнаружена неудовлетворительная точность на источниках тепла,
					// Около них сетку необходимо сгущать.
					bpowerON = true;
				}

				//if (bpowerON == false) {// так медленнее.
					//if ((b[ibcur].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (fabs(b[ibcur].arr_Sc[0]) > 0.0)) {
						// Обнаружена неудовлетворительная точность на источниках тепла,
						// Около них сетку необходимо сгущать.
						//bpowerON = true;
					//}
				//}

				//if (b[ibcur].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
				if (ib_marker_flag_fluid[iP]) {

					doublereal qgeom = 10.0;
					// Вычисляем координаты текущего КО.
					
					doublereal cposy = 0.5 * (ryboundary[iy + 1] + ryboundary[iy]);
					doublereal cposz = 0.5 * (rzboundary[iz + 1] + rzboundary[iz]);
					// Определяем номер блока к которому принадлежит данный КО.
					//for (integer i99 = 0; i99 < lb; i99++) {
						//if ((cpos>b[i99].g.xS) && (cpos < b[i99].g.xE)&&
							//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE)&&
							//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
							//ibcur = i99;
						//}
					//}

					// Делаем проверку: Есть хоть один сосед (E,W) тоже FLUID с учётом геометрической прогрессии 10.0 ?
					// Т.е. если окажутся два FLUID соседа но у них отношение сторон больше 10=qgeom то разбивать всё равно надо большего пополам.
					if (i < inumboundaryx - 1) {
						doublereal cpos_plus = 0.5*(rxboundary[i + 2] + rxboundary[i+1]);
						integer ibcur_plus = 0; 
						// Определяем номер блока к которому принадлежит данный КО.
						//for (integer i99 = 0; i99 < lb; i99++) {
							//if ((cpos_plus>b[i99].g.xS) && (cpos_plus < b[i99].g.xE) &&
								//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE) &&
								//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
								//ibcur_plus = i99;
							//}
						//}
						ibcur_plus = ib_marker[iP+1];
						if (b[ibcur_plus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
							if (i > 0) {
								doublereal cpos_minus = 0.5*(rxboundary[i] + rxboundary[i - 1]);
								integer ibcur_minus = 0;
								// Определяем номер блока к которому принадлежит данный КО.
								//for (integer i99 = 0; i99 < lb; i99++) {
									//if ((cpos_minus>b[i99].g.xS) && (cpos_minus < b[i99].g.xE) &&
										//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE) &&
										//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
										//ibcur_minus = i99;
									//}
								//}
								ibcur_minus = ib_marker[iP - 1];
								if (b[ibcur_minus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
									// ОДИНОЧНЫЙ FLUID блок
									bfound_onex_fluid_cv = true;
									break;
								}
								else {
									// Нижний сосед есть и он тоже FLUID
									// Проверка на qgeom
									doublereal q1=fabs((rxboundary[i] - rxboundary[i - 1]) / (rxboundary[i + 1] - rxboundary[i]));
									doublereal qgeom_inv = 1.0 / qgeom;
									if (q1 < qgeom_inv) {
										bfound_onex_fluid_cv = true;
										break;
									}
									if (q1 > qgeom) {
										// W сторону надо бить пополам, но этот случай обрабатывается при обработке W чейки. 
									}
								}
							}
							else {
								// А нижнего соседа просто нету !!!
								// ОДИНОЧНЫЙ FLUID блок
								bfound_onex_fluid_cv = true;
								break;
							}
						}
						else {
							// Верхняя ячейка точно FLUID
							// Дальнейшее рассмотрение имеет смысл если толщина верхней ячейки не удовлетворяет qgeom критерию
							// т.е. она слишком узкая
							doublereal q2 = fabs((rxboundary[i+2] - rxboundary[i+1]) / (rxboundary[i + 1] - rxboundary[i]));
							doublereal qgeom_inv = 1.0 / qgeom;
							if (q2 < qgeom_inv) {
							     // Верхняя ячейка неудовлетворительна,
								// переходим к рассмотрению нижней.
								if (i > 0) {
									doublereal cpos_minus = 0.5*(rxboundary[i] + rxboundary[i - 1]);
									integer ibcur_minus = 0;
									// Определяем номер блока к которому принадлежит данный КО.
									//for (integer i99 = 0; i99 < lb; i99++) {
										//if ((cpos_minus>b[i99].g.xS) && (cpos_minus < b[i99].g.xE) &&
											//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE) &&
											//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
											//ibcur_minus = i99;
									//	}
									//}
									//integer iP = i + inumboundaryx*iy + inumboundaryx*inumboundaryy*iz;
									ibcur_minus = ib_marker[iP-1];
									if (b[ibcur_minus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
										// ОДИНОЧНЫЙ FLUID блок
										bfound_onex_fluid_cv = true;
										break;
									}
									else {
										// Нижний сосед есть и он тоже FLUID
										// Проверка на qgeom
										doublereal q1 = fabs((rxboundary[i] - rxboundary[i - 1]) / (rxboundary[i + 1] - rxboundary[i]));
										doublereal qgeom_inv = 1.0 / qgeom;
										if (q1 < qgeom_inv) {
											bfound_onex_fluid_cv = true;
											break;
										}
										if (q1 > qgeom) {
											// W сторону надо бить пополам, но этот случай обрабатывается при обработке W чейки. 
										}
									}
								}
								else {
									// А нижнего соседа просто нету !!!
									// ОДИНОЧНЫЙ FLUID блок
									bfound_onex_fluid_cv = true;
									break;
								}
							}
						}
					}
					else {
						// проверяем нижнего соседа
						// Верхняя чейка просто отсутствует
						if (i > 0) {
							doublereal cpos_minus = 0.5*(rxboundary[i] + rxboundary[i - 1]);
							integer ibcur_minus = 0;
							// Определяем номер блока к которому принадлежит данный КО.
							//for (integer i99 = 0; i99 < lb; i99++) {
								//if ((cpos_minus>b[i99].g.xS) && (cpos_minus < b[i99].g.xE) &&
									//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE) &&
									//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
									//ibcur_minus = i99;
								//}
							//}
							//integer iP = i + inumboundaryx*iy + inumboundaryx*inumboundaryy*iz;
							ibcur_minus = ib_marker[iP - 1];
							if (b[ibcur_minus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
								// ОДИНОЧНЫЙ FLUID блок
								bfound_onex_fluid_cv = true;
								break;
							}
							else {
								// Нижний сосед есть и он тоже FLUID
								// Проверка на qgeom
								doublereal q1 = fabs((rxboundary[i] - rxboundary[i - 1]) / (rxboundary[i + 1] - rxboundary[i]));
								doublereal qgeom_inv = 1.0 / qgeom;
								if (q1 < qgeom_inv) {
									bfound_onex_fluid_cv = true;
									break;
								}
								if (q1 > qgeom) {
									// W сторону надо бить пополам, но этот случай обрабатывается при обработке W чейки. 
								}
							}
						}
						else {
							// Это означает что у нас всего вообще одна клетка.
							// Значит надо её разбить.
							bfound_onex_fluid_cv = true;
							break;
						}
					}
				}
			}
		}
		
		if (bfound_onex_fluid_cv==false) {
			// SOLID
			if (brefinement&&bpowerON) {
				ixintervalcount[i] = 3; // 23.12.2018
			}
			else {
				ixintervalcount[i] = 1;
			}
		}
		else {
			// FLUID
			bool b2div = false;
			for (integer i_3 = 0; (b2div==false)&&(i_3 < (inumboundaryy)); i_3++) {
				for (integer i_4 = 0; (b2div == false) && (i_4 < (inumboundaryz)); i_4++) {

					doublereal yp_3= 0.5*(ryboundary[i_3 + 1] + ryboundary[i_3]);
					doublereal zp_3 = 0.5*(rzboundary[i_4 + 1] + rzboundary[i_4]);
					doublereal xp_1= 0.5*(rxboundary[i + 1] + rxboundary[i]);
					doublereal xp_2 = xp_1;
					doublereal xp_3 = xp_1;
					if (i < inumboundaryx - 1) {
					    xp_2 = 0.5*(rxboundary[i + 2] + rxboundary[i+1]);
					}
					if (i > 0) {
						xp_3 = 0.5*(rxboundary[i - 1] + rxboundary[i]);
					}
					// определяет номер блока по координате точки.
					//myisblock_id(integer lb, BLOCK* &b, doublereal x11, doublereal y11, doublereal z11)
					// Функция myisblock_id ресурсоемкая. Она содержит вложенный цикл.
					integer ib1 = myisblock_id(lb, b, xp_1, yp_3, zp_3);
					integer ib2 = myisblock_id(lb, b, xp_2, yp_3, zp_3);
					integer ib3 = myisblock_id(lb, b, xp_3, yp_3, zp_3);
					if (!((((b[ib1].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[ib1].itype == PHYSICS_TYPE_IN_BODY::SOLID))
						&& ((b[ib2].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[ib2].itype == PHYSICS_TYPE_IN_BODY::SOLID)) 
						&& ((b[ib3].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[ib3].itype == PHYSICS_TYPE_IN_BODY::SOLID))) 
						|| ((b[ib1].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[ib2].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[ib3].itype == PHYSICS_TYPE_IN_BODY::SOLID)) 
						|| (((b[ib1].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[ib1].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) &&
						   ((b[ib2].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[ib2].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) && 
						   ((b[ib3].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[ib3].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)))))
					{
						// нету подряд трех кубиков вдоль линии Ох принадлежащих одновременно solid или fluid.
						b2div = true;
					}
				}
			}
			if (b2div) {
				ixintervalcount[i] = 2;
			}
			else {
				printf("inc X");
				//getchar();
				ixintervalcount[i] = 1;
			}
		}

		if (SPARSE_MESHER) ixintervalcount[i] = 1;// Одна клетка. И так очень много геометрических объектов.
	}

	

	 // debug
	/*
	for (i=0; i<(inumboundaryx); i++) {
	#if doubleintprecision == 1
		printf("%lld  ",ixintervalcount[i]);
	#else
		printf("%d  ",ixintervalcount[i]);
	#endif
		getchar();
	} //*/
	//getchar();

	

	// заполнение масива положения узлов.
	int iposmark = 1;
	doublereal dx;
	//integer k;
	integer ixoldsize = 0;
	SetLength(xpos, ixoldsize, 1);
	ixoldsize = 1;
	for (integer i = 0; i<(inumboundaryx); i++)
	{
		//if ((ixintervalcount[i] - 2) % 2 == 0) ixintervalcount[i]++; // чтобы число узлов было всегда чётное.
		//integer n2 = (integer)((ixintervalcount[i] - 1) / 2); // округление в меньшую сторону
		//if (bgeomx) n2 = ibalancen2(n2, qxW, qxE, rxboundary, ixintervalcount[i], i, iposmark); // балансировка 
		//doublereal qn2W = qxW;
		//for (integer i1 = 1; i1<n2; i1++) qn2W *= qxW; // возведение в степень.
		//doublereal b1W = (rxboundary[i + 1] - rxboundary[i])*(qxW - 1.0) / (2.0*(qn2W - 1.0));
		//doublereal qn2E = qxE;
		//for (integer i1 = 1; i1<(ixintervalcount[i] - n2); i1++) qn2E *= qxE; // возведение в степень.
		//doublereal b1E = (rxboundary[i + 1] - rxboundary[i])*(qxE - 1.0) / (2.0*(qn2E - 1.0));
		//printf("length=%e\n",(rxboundary[i+1]-rxboundary[i]));
		//getchar(); // OK

		dx = (rxboundary[i + 1] - rxboundary[i]) / (ixintervalcount[i]);
		SetLength(xpos, ixoldsize, (iposmark + ixintervalcount[i]+1));//+1
		ixoldsize = (iposmark + ixintervalcount[i]+1);//+1
#if doubleintprecision == 1
		//printf("%lld  ",ixoldsize);// debug
#else
		//printf("%d  ",ixoldsize);// debug
#endif
		
		for (integer k = iposmark; k <= iposmark + ixintervalcount[i] - 2+1; k++)
		{
			if (!bgeomx) {
				// равномерная сетка
				xpos[k] = rxboundary[i] + (k - iposmark)*dx;
			}
			/*else
			{
				// неравномерная сетка
				integer ic1 = k - iposmark;
				doublereal madd;
				if (ic1<n2) {
					// сначала левая сторона.
					madd = b1W;
					for (integer i1 = 1; i1<ic1; i1++) madd *= qxW;
					if (ic1 == n2 - 1) madd = 0.5*(rxboundary[i + 1] + rxboundary[i]) - xpos[k - 1]; // чтобы избежать ошибок округления
				}
				else
				{
					// потом правая сторона.
					madd = b1E;
					for (integer i1 = ixintervalcount[i] - 1 - ic1; i1>0; i1--) madd *= qxE;
				}
				if (k == iposmark) xpos[k] = rxboundary[i];
				else xpos[k] = xpos[k - 1] + madd;

			}
			*/
		}
		iposmark = iposmark + ixintervalcount[i]-1+1; //+1
	}
	SetLength(xpos, ixoldsize, iposmark + 1);
	xpos[iposmark] = rxboundary[inumboundaryx];
	inx = iposmark;
	for (integer i = 0; i<inx; i++) xpos[i] = xpos[i + 1]; // сдвиг влево на 1
	SetLength(xpos, inx + 1, inx);
	inx--; // индексация начинается с нуля и заканчивается значением inx

	// Нам нужна самая  грубая среди существующих сеток. 
	//for (i = 0; i<adapt_x; i++) simplecorrect_meshgen_x(xpos, inx, lb, ls, lw, b, s, w);

	/*
	// debug
	for (i=0; i<=inx; i++) {
	printf("%f  ",xpos[i]);
	}
	getchar(); //
	*/

	// 3 сентября 2017.
	// Добавление новых сеточных линий.
	// необходимость их добавления обусловлена требованнем АЛИС сеток.
	for (int i_28 = 0; i_28 <= inxadd; i_28++) {
		//SetLength(xpos, inx+1, inx+2);
		//xpos[inx+1] = xposadd[i_28];
		//inx++;
		// Это добавление умное, препятствует возникновению одинаковых значений в массиве.
		addboundary(xpos, inx, xposadd[i_28],YZ_PLANE, b, lb, w, lw, s, ls);
	}
	Sort_method<doublereal>(xpos,inx);

	if (1) {
		//const doublereal etalon_max_size_ratio=2.0;
		integer inum_iter_ratio_good = 0;
		doublereal max_size_ratio_x = 1.0;
		while (1) {

			// Вычисление max size ratio x axis.
			max_size_ratio_x = 1.0;
			for (integer i = 0; i<inx - 1; i++) {

				doublereal dmax = 0.0;
				doublereal dmin = 0.0;
				if (fabs(xpos[i + 1] - xpos[i])>fabs(xpos[i + 2] - xpos[i + 1])) {
					dmax = fabs(xpos[i + 1] - xpos[i]);
					dmin = fabs(xpos[i + 2] - xpos[i + 1]);
				}
				else {
					dmax = fabs(xpos[i + 2] - xpos[i + 1]);
					dmin = fabs(xpos[i + 1] - xpos[i]);
				}
				if (dmax / dmin>max_size_ratio_x) {
					max_size_ratio_x = dmax / dmin;
				}
			}
			
			if (max_size_ratio_x != max_size_ratio_x) {
				//printf("x axis max size ratio is equal = %1.4f\n", max_size_ratio_x);
				std::cout << "x axis max size ratio is equal = " << max_size_ratio_x << std::endl;
				system("PAUSE");
			}
			//getchar();

			if (max_size_ratio_x < (etalon_max_size_ratio*1.1)) {
				break;
			}

			// Корректировка max size ratio x axis.
			// сделано 9 апреля 2013 года.

			bool bplus = false;
			max_size_ratio_x = 1.0;
			for (integer i = 0; i<inx - 1; i++) {

				bplus = false;
				doublereal dmax = 0.0;
				doublereal dmin = 0.0;
				if (fabs(xpos[i + 1] - xpos[i])>fabs(xpos[i + 2] - xpos[i + 1])) {
					dmax = fabs(xpos[i + 1] - xpos[i]);
					dmin = fabs(xpos[i + 2] - xpos[i + 1]);
				}
				else {
					bplus = true;
					dmax = fabs(xpos[i + 2] - xpos[i + 1]);
					dmin = fabs(xpos[i + 1] - xpos[i]);
				}
				if (dmax / dmin>etalon_max_size_ratio) {
					doublereal pos_candidate;
					if (bplus) {
						pos_candidate = xpos[i + 1] + 0.5*dmax;
					}
					else {
						pos_candidate = xpos[i] + 0.5*dmax;
					}
					SetLength(xpos, inx + 1, inx + 2);
					xpos[inx + 1] = pos_candidate;
					inx = inx + 1;
					//BubbleEnhSort<doublereal>(xpos, 0, inx);
					Sort_method<doublereal>(xpos,inx);
					break;
				}
			}

			inum_iter_ratio_good++;

			//getchar();

		}

		//printf("x axis max size ratio is equal = %1.4f\n", max_size_ratio_x);
		std::cout << "x axis max size ratio is equal = " << max_size_ratio_x;
		std::cout << std::endl;

#if doubleintprecision == 1
		printf("inum_iter_ratio_good is %lld\n", inum_iter_ratio_good);
#else
		printf("inum_iter_ratio_good is %d\n", inum_iter_ratio_good);
#endif
		
	}


	int *iyintervalcount; // число интервалов
	iyintervalcount = new int[inumboundaryy]; // на один меньше чем число границ.
	//for (integer i = 0; i<(inumboundaryy); i++) {
		//alphascal = (ryboundary[i + 1] - ryboundary[i]) / (ryboundary[inumboundaryy] - rxboundary[0]);
		//inowintervalcount = (integer)(alphascal*iny);
		//if (inowintervalcount < min_elem_in_y_element) inowintervalcount = min_elem_in_y_element;
		//iyintervalcount[i] = inowintervalcount;
	//}

	
#pragma omp parallel for
	for (integer i = 0; i < (inumboundaryy); i++) {

		bool bpowerON = false;


		// Если мы в solide то одна клетка.
		// Если во Fluide то две клетки.
		doublereal cpos = 0.5 * (ryboundary[i + 1] + ryboundary[i]);
		// В плоскости найден одиночный fluid контрольный объем имющий 
		// по E и W сторонам только Solid и hollow соседей.
		bool bfound_onex_fluid_cv = false;
		// применим нисходящее проектирование сверху вниз.
		//integer ibcur = 0; // номер текущего блока (кабинет по умолчанию).
		// сканируем плоскость.
		for (integer ix = 0; (bfound_onex_fluid_cv == false)&&(ix < (inumboundaryx)); ix++) {
			for (integer iz = 0;  (bfound_onex_fluid_cv == false) && (iz < (inumboundaryz)); iz++) {
				doublereal qgeom = 10.0;
				// Вычисляем координаты текущего КО.
				
				doublereal cposx = 0.5*(rxboundary[ix + 1] + rxboundary[ix]);
				doublereal cposz = 0.5*(rzboundary[iz + 1] + rzboundary[iz]);
				// Определяем номер блока к которому принадлежит данный КО.
				//for (integer i99 = 0; i99 < lb; i99++) {
					//if ((cpos>b[i99].g.yS) && (cpos < b[i99].g.yE) &&
						//(cposx>b[i99].g.xS) && (cposx < b[i99].g.xE) &&
						//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
						//ibcur = i99;
					//}
				//}
				integer iP = ix + inumboundaryx*i + inumboundaryx*inumboundaryy*iz;
				//ibcur = ib_marker[iP];

				if (ib_marker_flag_power[iP]) {
					// Обнаружена неудовлетворительная точность на источниках тепла,
					// Около них сетку необходимо сгущать.
					bpowerON = true;
				}

				//if ((b[ibcur].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (fabs(b[ibcur].arr_Sc[0]) > 0.0)) {
					// Обнаружена неудовлетворительная точность на источниках тепла,
					// Около них сетку необходимо сгущать.
					//bpowerON = true;
				//}

				//if (b[ibcur].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
				if (ib_marker_flag_fluid[iP]) {
					// Делаем проверочку: Есть хоть один сосед (N,S) тоже FLUID с учётом геометрической прогрессии 10.0 ?
					// Т.е. если окажутся два FLUID соседа но у них отношение сторон больше 10=qgeom то разбивать всё равно надо большего пополам.
					if (i < inumboundaryy - 1) {
						doublereal cpos_plus = 0.5*(ryboundary[i + 2] + ryboundary[i + 1]);
						integer ibcur_plus = 0;
						// Определяем номер блока к которому принадлежит данный КО.
						//for (integer i99 = 0; i99 < lb; i99++) {
							//if ((cpos_plus>b[i99].g.yS) && (cpos_plus < b[i99].g.yE) &&
								//(cposx>b[i99].g.xS) && (cposx < b[i99].g.xE) &&
								//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
								//ibcur_plus = i99;
							//}
						//}
						ibcur_plus = ib_marker[ix + inumboundaryx*(i+1) + inumboundaryx*inumboundaryy*iz];
						if (b[ibcur_plus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
							if (i > 0) {
								doublereal cpos_minus = 0.5*(ryboundary[i] + ryboundary[i - 1]);
								integer ibcur_minus = 0;
								// Определяем номер блока к которому принадлежит данный КО.
								//for (integer i99 = 0; i99 < lb; i99++) {
									//if ((cpos_minus>b[i99].g.yS) && (cpos_minus < b[i99].g.yE) &&
										//(cposx>b[i99].g.xS) && (cposx < b[i99].g.xE) &&
										//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
										//ibcur_minus = i99;
									//}
								//}
								ibcur_minus = ib_marker[ix + inumboundaryx*(i - 1) + inumboundaryx*inumboundaryy*iz];
								if (b[ibcur_minus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
									// ОДИНОЧНЫЙ FLUID блок
									bfound_onex_fluid_cv = true;
									break;
								}
								else {
									// Нижний сосед есть и он тоже FLUID
									// Проверка на qgeom
									doublereal q1 = fabs((ryboundary[i] - ryboundary[i - 1]) / (ryboundary[i + 1] - ryboundary[i]));
									doublereal qgeom_inv = 1.0 / qgeom;
									if (q1 < qgeom_inv) {
										bfound_onex_fluid_cv = true;
										break;
									}
									if (q1 > qgeom) {
										// W сторону надо бить пополам, но этот случай обрабатывается при обработке W чейки. 
									}
								}
							}
							else {
								// А нижнего соседа просто нету !!!
								// ОДИНОЧНЫЙ FLUID блок
								bfound_onex_fluid_cv = true;
								break;
							}
						}
						else {
							// Верхняя ячейка точно FLUID
							// Дальнейшее рассмотрение имеет смысл если толщина верхней ячейки не удовлетворяет qgeom критерию
							// т.е. она слишком узкая
							doublereal q2 = fabs((ryboundary[i + 2] - ryboundary[i + 1]) / (ryboundary[i + 1] - ryboundary[i]));
							doublereal qgeom_inv = 1.0 / qgeom;
							if (q2 < qgeom_inv) {
								// Верхняя ячейка неудовлетворительна,
								// переходим к рассмотрению нижней.
								if (i > 0) {
									doublereal cpos_minus = 0.5*(ryboundary[i] + ryboundary[i - 1]);
									integer ibcur_minus = 0;
									// Определяем номер блока к которому принадлежит данный КО.
									//for (integer i99 = 0; i99 < lb; i99++) {
										//if ((cpos_minus>b[i99].g.yS) && (cpos_minus < b[i99].g.yE) &&
											//(cposx>b[i99].g.xS) && (cposx < b[i99].g.xE) &&
											//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
											//ibcur_minus = i99;
										//}
									//}
									ibcur_minus = ib_marker[ix + inumboundaryx*(i - 1) + inumboundaryx*inumboundaryy*iz];
									if (b[ibcur_minus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
										// ОДИНОЧНЫЙ FLUID блок
										bfound_onex_fluid_cv = true;
										break;
									}
									else {
										// Нижний сосед есть и он тоже FLUID
										// Проверка на qgeom
										doublereal q1 = fabs((ryboundary[i] - ryboundary[i - 1]) / (ryboundary[i + 1] - ryboundary[i]));
										doublereal qgeom_inv = 1.0 / qgeom;
										if (q1 < qgeom_inv) {
											bfound_onex_fluid_cv = true;
											break;
										}
										if (q1 > qgeom) {
											// W сторону надо бить пополам, но этот случай обрабатывается при обработке W чейки. 
										}
									}
								}
								else {
									// А нижнего соседа просто нету !!!
									// ОДИНОЧНЫЙ FLUID блок
									bfound_onex_fluid_cv = true;
									break;
								}
							}
						}
					}
					else {
						// проверяем нижнего соседа
						// Верхняя чейка просто отсутствует
						if (i > 0) {
							doublereal cpos_minus = 0.5*(ryboundary[i] + ryboundary[i - 1]);
							integer ibcur_minus = 0;
							// Определяем номер блока к которому принадлежит данный КО.
							//for (integer i99 = 0; i99 < lb; i99++) {
								//if ((cpos_minus>b[i99].g.yS) && (cpos_minus < b[i99].g.yE) &&
									//(cposx>b[i99].g.xS) && (cposx < b[i99].g.xE) &&
									//(cposz>b[i99].g.zS) && (cposz < b[i99].g.zE)) {
									//ibcur_minus = i99;
								//}
							//}
							ibcur_minus = ib_marker[ix + inumboundaryx*(i - 1) + inumboundaryx*inumboundaryy*iz];
							if (b[ibcur_minus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
								// ОДИНОЧНЫЙ FLUID блок
								bfound_onex_fluid_cv = true;
								break;
							}
							else {
								// Нижний сосед есть и он тоже FLUID
								// Проверка на qgeom
								doublereal q1 = fabs((ryboundary[i] - ryboundary[i - 1]) / (ryboundary[i + 1] - ryboundary[i]));
								doublereal qgeom_inv = 1.0 / qgeom;
								if (q1 < qgeom_inv) {
									bfound_onex_fluid_cv = true;
									break;
								}
								if (q1 > qgeom) {
									// W сторону надо бить пополам, но этот случай обрабатывается при обработке W чейки. 
								}
							}
						}
						else {
							// Это означает что у нас всего вообще одна клетка.
							// Значит надо её разбить.
							bfound_onex_fluid_cv = true;
							break;
						}
					}
				}
			}
		}

		if (bfound_onex_fluid_cv == false) {

			if (brefinement&&bpowerON) {
				iyintervalcount[i] = 3; // 23.12.2018
			}
			else {
				iyintervalcount[i] = 1;
			}

		}
		else {
			//iyintervalcount[i] = 2;


			// FLUID
			bool b2div = false;
			for (integer i_3 = 0; (b2div == false)&&(i_3 < (inumboundaryx)); i_3++) {
				for (integer i_4 = 0; (b2div == false) && (i_4 < (inumboundaryz)); i_4++) {
					doublereal xp_3 = 0.5*(rxboundary[i_3 + 1] + rxboundary[i_3]);
					doublereal zp_3 = 0.5*(rzboundary[i_4 + 1] + rzboundary[i_4]);
					doublereal yp_1 = 0.5*(ryboundary[i + 1] + ryboundary[i]);
					doublereal yp_2 = yp_1;
					doublereal yp_3 = yp_1;
					if (i < inumboundaryy - 1) {
						yp_2 = 0.5*(ryboundary[i + 2] + ryboundary[i + 1]);
					}
					if (i > 0) {
						yp_3 = 0.5*(ryboundary[i - 1] + ryboundary[i]);
					}
					// определяет номер блока по координате точки.
					//myisblock_id(integer lb, BLOCK* &b, doublereal x11, doublereal y11, doublereal z11)
					// Функция myisblock_id ресурсоемкая. Она содержит вложенный цикл.
					integer ib1 = myisblock_id(lb, b, xp_3, yp_1, zp_3);
					integer ib2 = myisblock_id(lb, b, xp_3, yp_2, zp_3);
					integer ib3 = myisblock_id(lb, b, xp_3, yp_3, zp_3);
					if (!((((b[ib1].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[ib1].itype == PHYSICS_TYPE_IN_BODY::SOLID))
						&& ((b[ib2].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[ib2].itype == PHYSICS_TYPE_IN_BODY::SOLID)) 
						&& ((b[ib3].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[ib3].itype == PHYSICS_TYPE_IN_BODY::SOLID))) 
						|| ((b[ib1].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[ib2].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[ib3].itype == PHYSICS_TYPE_IN_BODY::SOLID)) 
						|| (((b[ib1].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[ib1].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) &&
						((b[ib2].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[ib2].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) &&
							((b[ib3].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[ib3].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)))))
					{
						// нету подряд трех кубиков вдоль линии Ох принадлежащих одновременно solid или fluid.
						b2div = true;
					}
				}
			}
			if (b2div) {
				iyintervalcount[i] = 2;
			}
			else {
				printf("inc Y");
				//getchar();
				iyintervalcount[i] = 1;
			}

		}
		
		if (SPARSE_MESHER) iyintervalcount[i] = 1;// Одна клетка. И так очень много геометрических объектов.
	}


	// заполнение массива узлов
	iposmark = 1;
	doublereal dy;
	integer iyoldsize = 0;
	SetLength(ypos, iyoldsize, 1);
	iyoldsize = 1;
	for (integer i = 0; i<inumboundaryy; i++)
	{

		// чтобы число узлов було всегда чётное
	//	if ((iyintervalcount[i] - 2) % 2 == 0) iyintervalcount[i]++;
		//integer n2 = (integer)((iyintervalcount[i] - 1) / 2); // округление в меньшую сторону
		//if (bgeomy) n2 = ibalancen2(n2, qyS, qyN, ryboundary, iyintervalcount[i], i, iposmark); // балансировка 
		//doublereal qn2S = qyS;
		//for (integer i1 = 1; i1<n2; i1++) qn2S *= qyS; // возведение в степень.
		//doublereal b1S = (ryboundary[i + 1] - ryboundary[i])*(qyS - 1.0) / (2.0*(qn2S - 1.0));
		//doublereal qn2N = qyN;
		//for (integer i1 = 1; i1<(iyintervalcount[i] - n2); i1++) qn2N *= qyN; // возведение в степень.
		//doublereal b1N = (ryboundary[i + 1] - ryboundary[i])*(qyN - 1.0) / (2.0*(qn2N - 1.0));

		dy = (ryboundary[i + 1] - ryboundary[i]) / (iyintervalcount[i]);
		SetLength(ypos, iyoldsize, iposmark + iyintervalcount[i]+1);
		iyoldsize = iposmark + iyintervalcount[i]+1;
		for (integer k = iposmark; k <= iposmark + iyintervalcount[i] - 2+1; k++)
		{
			if (!bgeomy) {
				// равномерная сетка
				ypos[k] = ryboundary[i] + (k - iposmark)*dy;
			}
			/*else
			{
				// неравномерная сетка
				// на основе геометрической прогрессии с двумя сгущениями.
				integer ic1 = k - iposmark;
				doublereal madd;
				if (ic1<n2) {
					madd = b1S;
					for (integer i1 = 1; i1<ic1; i1++) madd *= qyS;
					if (ic1 == n2 - 1) madd = 0.5*(ryboundary[i + 1] + ryboundary[i]) - ypos[k - 1];
				}
				else
				{
					madd = b1N;
					for (integer i1 = iyintervalcount[i] - 1 - ic1; i1>0; i1--) madd *= qyN;
				}
				if (k == iposmark) ypos[k] = ryboundary[i];
				else ypos[k] = ypos[k - 1] + madd;
			}*/
		}
		iposmark = iposmark + iyintervalcount[i] - 1+1;
	}
	SetLength(ypos, iyoldsize, iposmark + 1);
	ypos[iposmark] = ryboundary[inumboundaryy];
	iny = iposmark;
	for (integer i = 0; i<iny; i++) ypos[i] = ypos[i + 1]; // сдвиг влево на 1
	SetLength(ypos, iny + 1, iny);
	iny--; // индексация начинается с нуля и заканчивается значением iny

	//for (i = 0; i<adapt_y; i++) simplecorrect_meshgen_y(ypos, iny, lb, ls, lw, b, s, w);


	/*
	// debug
	for (i=0; i<=iny; i++) {
	printf("%f  ",ypos[i]);
	}//
	*/

	

	// 3 сентября 2017.
	// Добавление новых сеточных линий.
	// необходимость их добавления обусловлена требованнем АЛИС сеток.
	for (int i_28 = 0; i_28 <= inyadd; i_28++) {
		//SetLength(ypos, iny+1, iny + 2);
		//ypos[iny+1] = yposadd[i_28];
		//iny++;
		// Это добавление умное, препятствует возникновению одинаковых значений в массиве.
		addboundary(ypos, iny, yposadd[i_28],XZ_PLANE, b, lb, w, lw, s, ls);
	}
	Sort_method<doublereal>(ypos,iny);
	
	if (1) {
		integer inum_iter_ratio_good = 0;
		doublereal max_size_ratio_y = 1.0;
		while (1) {

			max_size_ratio_y = 1.0;
			for (integer i = 0; i<iny - 1; i++) {

				doublereal dmax = 0.0;
				doublereal dmin = 0.0;
				if (fabs(ypos[i + 1] - ypos[i])>fabs(ypos[i + 2] - ypos[i + 1])) {
					dmax = fabs(ypos[i + 1] - ypos[i]);
					dmin = fabs(ypos[i + 2] - ypos[i + 1]);
				}
				else {
					dmax = fabs(ypos[i + 2] - ypos[i + 1]);
					dmin = fabs(ypos[i + 1] - ypos[i]);
				}
				if (dmax / dmin>max_size_ratio_y) {
					max_size_ratio_y = dmax / dmin;
				}
			}
			
			if (max_size_ratio_y != max_size_ratio_y) {
				//printf("y axis max size ratio is equal = %1.4f\n", max_size_ratio_y);
				std::cout << "y axis max size ratio is equal = " << max_size_ratio_y << std::endl;
				system("PAUSE");
			}
			//getchar();

			if (max_size_ratio_y < (etalon_max_size_ratio*1.1)) {
				break;
			}

			// Корректировка max size ratio y axis.
			// сделано 9 апреля 2013 года.

			bool bplus = false;
			max_size_ratio_y = 1.0;
			for (integer i = 0; i<iny - 1; i++) {

				bplus = false;
				doublereal dmax = 0.0;
				doublereal dmin = 0.0;
				if (fabs(ypos[i + 1] - ypos[i])>fabs(ypos[i + 2] - ypos[i + 1])) {
					dmax = fabs(ypos[i + 1] - ypos[i]);
					dmin = fabs(ypos[i + 2] - ypos[i + 1]);
				}
				else {
					bplus = true;
					dmax = fabs(ypos[i + 2] - ypos[i + 1]);
					dmin = fabs(ypos[i + 1] - ypos[i]);
				}
				if (dmax / dmin>etalon_max_size_ratio) {
					doublereal pos_candidate;
					if (bplus) {
						pos_candidate = ypos[i + 1] + 0.5*dmax;
					}
					else {
						pos_candidate = ypos[i] + 0.5*dmax;
					}
					SetLength(ypos, iny + 1, iny + 2);
					ypos[iny + 1] = pos_candidate;
					iny = iny + 1;
					//BubbleEnhSort<doublereal>(ypos, 0, iny);
					Sort_method<doublereal>(ypos,iny);
					break;
				}
			}

			inum_iter_ratio_good++;

			//getchar();

		}

		//printf("y axis max size ratio is equal = %1.4f\n", max_size_ratio_y);
		std::cout << "y axis max size ratio is equal = " << max_size_ratio_y << std::endl;

#if doubleintprecision == 1
		printf("inum_iter_ratio_good is %lld\n", inum_iter_ratio_good);
#else
		printf("inum_iter_ratio_good is %d\n", inum_iter_ratio_good);
#endif
		
	}

	int *izintervalcount; // число интервалов
	izintervalcount = new int[inumboundaryz]; // на один меньше чем число границ.
	//for (i = 0; i<(inumboundaryz); i++) {
		//alphascal = (rzboundary[i + 1] - rzboundary[i]) / (rzboundary[inumboundaryz] - rzboundary[0]);
		//inowintervalcount = (integer)(alphascal*inz);
		//if (inowintervalcount < min_elem_in_z_element) inowintervalcount = min_elem_in_z_element;
	//	izintervalcount[i] = inowintervalcount;
	//}

	/*Обязательно оставить для возможного отката назад.
	for (i = 0; i < (inumboundaryz); i++) {
		// Если мы в solide то одна клетка.
		// Если во Fluide то две клетки.
		doublereal cpos = 0.5*(rzboundary[i + 1] + rzboundary[i]);
		integer ibcur = 0; // номер текущего блока (кабинет по умолчанию).
		for (integer i99 = 0; i99 < lb; i99++) {
			if ((cpos>b[i99].g.zS) && (cpos < b[i99].g.zE)) {
				ibcur = i99;
			}
		}
		if (b[ibcur].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
			izintervalcount[i] = 1;
		}
		if (b[ibcur].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
			izintervalcount[i] = 2;
		}
		if (b[ibcur].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
			izintervalcount[i] = 1;
		}
	}
	*/

#pragma omp parallel for
	for (integer i = 0; i < (inumboundaryz); i++) {

		bool bpowerON = false;


		// Если мы в solide то одна клетка.
		// Если во Fluide то две клетки.
		doublereal cpos = 0.5*(rzboundary[i + 1] + rzboundary[i]);
		// В плоскости найден одиночный fluid контрольный объем имющий 
		// по E и W сторонам только Solid и hollow соседей.
		bool bfound_onex_fluid_cv = false;
		// применим нисходящее проектирование сверху вниз.
		//integer ibcur = 0; // номер текущего блока (кабинет по умолчанию).
		// сканируем плоскость.
		for (integer ix = 0; (bfound_onex_fluid_cv == false)&&( ix < (inumboundaryx)); ix++) {
			for (integer iy = 0; (bfound_onex_fluid_cv == false) && (iy < (inumboundaryy)); iy++) {
				doublereal qgeom = 10.0;
				// Вычисляем координаты текущего КО.
				doublereal cposx = 0.5*(rxboundary[ix + 1] + rxboundary[ix]);
				doublereal cposy = 0.5*(ryboundary[iy + 1] + ryboundary[iy]);
				// Определяем номер блока к которому принадлежит данный КО.
				//for (integer i99 = 0; i99 < lb; i99++) {
					//if ((cpos>b[i99].g.zS) && (cpos < b[i99].g.zE) &&
						//(cposx>b[i99].g.xS) && (cposx < b[i99].g.xE) &&
						//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE)) {
						//ibcur = i99;
					//}
				//}

				
				integer iP = ix + inumboundaryx*iy + inumboundaryx*inumboundaryy*i;
				//ibcur = ib_marker[iP];

				if (ib_marker_flag_power[iP]) {
					// Обнаружена неудовлетворительная точность на источниках тепла,
					// Около них сетку необходимо сгущать.
					bpowerON = true;
				}
				//if ((b[ibcur].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (fabs(b[ibcur].arr_Sc[0]) > 0.0)) {
					// Обнаружена неудовлетворительная точность на источниках тепла,
					// Около них сетку необходимо сгущать.
					//bpowerON = true;
				//}


				//if (b[ibcur].itype == PHYSICS_TYPE_IN_BODY::FLUID) {
				if (ib_marker_flag_fluid[iP]) {
					// Делаем проверочку: Есть хоть один сосед (E,W) тоже FLUID с учётом геометрической прогрессии 10.0 ?
					// Т.е. если окажутся два FLUID соседа но у них отношение сторон больше 10=qgeom то разбивать всё равно надо большего пополам.
					if (i < inumboundaryz - 1) {
						doublereal cpos_plus = 0.5*(rzboundary[i + 2] + rzboundary[i + 1]);
						integer ibcur_plus = 0;
						// Определяем номер блока к которому принадлежит данный КО.
						//for (integer i99 = 0; i99 < lb; i99++) {
							//if ((cpos_plus>b[i99].g.zS) && (cpos_plus < b[i99].g.zE) &&
								//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE) &&
								//(cposx>b[i99].g.xS) && (cposx < b[i99].g.xE)) {
								//ibcur_plus = i99;
							//}
						//}
						ibcur_plus = ib_marker[ix + inumboundaryx*iy + inumboundaryx*inumboundaryy*(i+1)];
						if (b[ibcur_plus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
							if (i > 0) {
								doublereal cpos_minus = 0.5*(rzboundary[i] + rzboundary[i - 1]);
								integer ibcur_minus = 0;
								// Определяем номер блока к которому принадлежит данный КО.
								//for (integer i99 = 0; i99 < lb; i99++) {
									//if ((cpos_minus>b[i99].g.zS) && (cpos_minus < b[i99].g.zE) &&
										//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE) &&
										//(cposx>b[i99].g.xS) && (cposx < b[i99].g.xE)) {
										//ibcur_minus = i99;
									//}
								//}
								ibcur_minus = ib_marker[ix + inumboundaryx*iy + inumboundaryx*inumboundaryy*(i-1)];
								if (b[ibcur_minus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
									// ОДИНОЧНЫЙ FLUID блок
									bfound_onex_fluid_cv = true;
									break;
								}
								else {
									// Нижний сосед есть и он тоже FLUID
									// Проверка на qgeom
									doublereal q1 = fabs((rzboundary[i] - rzboundary[i - 1]) / (rzboundary[i + 1] - rzboundary[i]));
									doublereal qgeom_inv = 1.0 / qgeom;
									if (q1 < qgeom_inv) {
										bfound_onex_fluid_cv = true;
										break;
									}
									if (q1 > qgeom) {
										// W сторону надо бить пополам, но этот случай обрабатывается при обработке W чейки. 
									}
								}
							}
							else {
								// А нижнего соседа просто нету !!!
								// ОДИНОЧНЫЙ FLUID блок
								bfound_onex_fluid_cv = true;
								break;
							}
						}
						else {
							// Верхняя ячейка точно FLUID
							// Дальнейшее рассмотрение имеет смысл если толщина верхней ячейки не удовлетворяет qgeom критерию
							// т.е. она слишком узкая
							doublereal q2 = fabs((rzboundary[i + 2] - rzboundary[i + 1]) / (rzboundary[i + 1] - rzboundary[i]));
							doublereal qgeom_inv = 1.0 / qgeom;
							if (q2 < qgeom_inv) {
								// Верхняя ячейка неудовлетворительна,
								// переходим к рассмотрению нижней.
								if (i > 0) {
									doublereal cpos_minus = 0.5*(rzboundary[i] + rzboundary[i - 1]);
									integer ibcur_minus = 0;
									// Определяем номер блока к которому принадлежит данный КО.
									//for (integer i99 = 0; i99 < lb; i99++) {
										//if ((cpos_minus>b[i99].g.zS) && (cpos_minus < b[i99].g.zE) &&
											//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE) &&
											//(cposx>b[i99].g.xS) && (cposx < b[i99].g.xE)) {
											//ibcur_minus = i99;
										//}
									//}
									ibcur_minus = ib_marker[ix + inumboundaryx*iy + inumboundaryx*inumboundaryy*(i - 1)];
									if (b[ibcur_minus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
										// ОДИНОЧНЫЙ FLUID блок
										bfound_onex_fluid_cv = true;
										break;
									}
									else {
										// Нижний сосед есть и он тоже FLUID
										// Проверка на qgeom
										doublereal q1 = fabs((rzboundary[i] - rzboundary[i - 1]) / (rzboundary[i + 1] - rzboundary[i]));
										doublereal qgeom_inv = 1.0 / qgeom;
										if (q1 < qgeom_inv) {
											bfound_onex_fluid_cv = true;
											break;
										}
										if (q1 > qgeom) {
											// W сторону надо бить пополам, но этот случай обрабатывается при обработке W чейки. 
										}
									}
								}
								else {
									// А нижнего соседа просто нету !!!
									// ОДИНОЧНЫЙ FLUID блок
									bfound_onex_fluid_cv = true;
									break;
								}
							}
						}
					}
					else {
						// проверяем нижнего соседа
						// Верхняя чейка просто отсутствует
						if (i > 0) {
							doublereal cpos_minus = 0.5*(rzboundary[i] + rzboundary[i - 1]);
							integer ibcur_minus = 0;
							// Определяем номер блока к которому принадлежит данный КО.
							//for (integer i99 = 0; i99 < lb; i99++) {
								//if ((cpos_minus>b[i99].g.zS) && (cpos_minus < b[i99].g.zE) &&
									//(cposy>b[i99].g.yS) && (cposy < b[i99].g.yE) &&
									//(cposx > b[i99].g.xS) && (cposx < b[i99].g.xE)) {
									//ibcur_minus = i99;
								//}
							//}
							ibcur_minus = ib_marker[ix + inumboundaryx*iy + inumboundaryx*inumboundaryy*(i - 1)];
							if (b[ibcur_minus].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
								// ОДИНОЧНЫЙ FLUID блок
								bfound_onex_fluid_cv = true;
								break;
							}
							else {
								// Нижний сосед есть и он тоже FLUID
								// Проверка на qgeom
								doublereal q1 = fabs((rzboundary[i] - rzboundary[i - 1]) / (rzboundary[i + 1] - rzboundary[i]));
								doublereal qgeom_inv = 1.0 / qgeom;
								if (q1 < qgeom_inv) {
									bfound_onex_fluid_cv = true;
									break;
								}
								if (q1 > qgeom) {
									// W сторону надо бить пополам, но этот случай обрабатывается при обработке W чейки. 
								}
							}
						}
						else {
							// Это означает что у нас всего вообще одна клетка.
							// Значит надо её разбить.
							bfound_onex_fluid_cv = true;
							break;
						}
					}
				}
			}
		}

		if (bfound_onex_fluid_cv == false) {
			

			if (brefinement&&bpowerON) {
				izintervalcount[i] = 3; // 23.12.2018
			}
			else {
				izintervalcount[i] = 1;
			}

		}
		else {
			//izintervalcount[i] = 2;

			// FLUID
			bool b2div = false;
			for (integer i_3 = 0; (b2div == false)&&(i_3 < (inumboundaryx)); i_3++) {
				for (integer i_4 = 0; (b2div == false) && (i_4 < (inumboundaryy)); i_4++) {
					doublereal xp_3 = 0.5*(rxboundary[i_3 + 1] + rxboundary[i_3]);
					doublereal yp_3 = 0.5*(ryboundary[i_4 + 1] + ryboundary[i_4]);
					doublereal zp_1 = 0.5*(rzboundary[i + 1] + rzboundary[i]);
					doublereal zp_2 = zp_1;
					doublereal zp_3 = zp_1;
					if (i < inumboundaryz - 1) {
						zp_2 = 0.5*(rzboundary[i + 2] + rzboundary[i + 1]);
					}
					if (i > 0) {
						zp_3 = 0.5*(rzboundary[i - 1] + rzboundary[i]);
					}
					// определяет номер блока по координате точки.
					//myisblock_id(integer lb, BLOCK* &b, doublereal x11, doublereal y11, doublereal z11)
					if (!((((b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID))) || ((b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::SOLID) && (b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::SOLID)) || (((b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_1)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_2)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) && ((b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::FLUID)||(b[myisblock_id(lb, b, xp_3, yp_3, zp_3)].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)))))
					{
						// нету подряд трех кубиков вдоль линии Ох принадлежащих одновременно solid или fluid.
						b2div = true;
					}
				}
			}
			if (b2div) {
				izintervalcount[i] = 2;
			}
			else {
				printf("inc Z");
				//getchar();
				izintervalcount[i] = 1;
			}

		}
		if (SPARSE_MESHER) izintervalcount[i] = 1;// Одна клетка. И так очень много геометрических объектов.
	}

	// заполнение массива узлов
	iposmark = 1;
	doublereal dz;
	integer izoldsize = 0;
	SetLength(zpos, izoldsize, 1);
	izoldsize = 1;
	for (integer i = 0; i<(inumboundaryz); i++)
	{
		// чтобы число узлов було всегда чётное
		//if ((izintervalcount[i] - 2) % 2 == 0) izintervalcount[i]++;
		//integer n2 = (integer)((izintervalcount[i] - 1) / 2); // округление в меньшую сторону 
		//if (bgeomz) n2 = ibalancen2(n2, qzB, qzT, rzboundary, izintervalcount[i], i, iposmark); // балансировка 
#if doubleintprecision == 1
		//printf("n2=%lld\n",n2); // debug
#else
		//printf("n2=%d\n",n2); // debug
#endif
		
		//doublereal qn2B = qzB; // нижняя сторона
		// точного значения в 0.25 всё равно достигнуть не получится из-за возникающих ошибок округления.
		//for (integer i1 = 1; i1<n2; i1++) qn2B *= qzB; // возведение в степень.
		//doublereal b1B = (rzboundary[i + 1] - rzboundary[i])*(qzB - 1.0) / (2.0*(qn2B - 1.0));
		//doublereal qn2T = qzT; // верхняя сторона
		//for (integer i1 = 1; i1<(izintervalcount[i] - n2); i1++) qn2T *= qzT; // возведение в степень.
		//doublereal b1T = (rzboundary[i + 1] - rzboundary[i])*(qzT - 1.0) / (2.0*(qn2T - 1.0));

		dz = (rzboundary[i + 1] - rzboundary[i]) / (izintervalcount[i]);
		SetLength(zpos, izoldsize, iposmark + izintervalcount[i]+1);
		izoldsize = iposmark + izintervalcount[i]+1;
		for (integer k = iposmark; k <= iposmark + izintervalcount[i] - 2+1; k++)
		{
			if (!bgeomz) {
				// равномерная сетка
				zpos[k] = rzboundary[i] + (k - iposmark)*dz;
			}
			/*else
			{
				// неравномерная сетка
				// на основе геометрической прогрессии с двумя сгущениями.
				integer ic1 = k - iposmark;
				doublereal madd;
				if (ic1<n2) {
					madd = b1B;
					for (integer i1 = 0; i1<ic1; i1++) madd *= qzB;
					if (ic1 == n2 - 1) madd = 0.5*(rzboundary[i + 1] + rzboundary[i]) - zpos[k - 1]; // чтобы избежать ошибок округления.
					//printf("first.:");
				}
				else
				{
					madd = b1T;
					for (integer i1 = izintervalcount[i] - 1 - ic1; i1>0; i1--) madd *= qzT;
					// это тока предпоследний.
					//if (k==iposmark+izintervalcount[i]-2) madd*=rzboundary[i+1]-zpos[k-1]; // чтобы избежать ошибок округления.
					//printf("second.:");
				}
				if (k == iposmark) zpos[k] = rzboundary[i]; // чтобы избежать ошибок округления.
				else {
					// учитывается что это начиная со второго до предпоследнего
					zpos[k] = zpos[k - 1] + madd;
				}
				#if doubleintprecision == 1
					//printf("zpos[%lld]=%f, madd=%f\n",k,zpos[k],madd);
				#else
					//printf("zpos[%d]=%f, madd=%f\n",k,zpos[k],madd);
				#endif
				
				//getchar();
			}*/
		}
		iposmark = iposmark + izintervalcount[i] - 1+1;
	}
	SetLength(zpos, izoldsize, iposmark + 1);
	zpos[iposmark] = rzboundary[inumboundaryz];
#if doubleintprecision == 1
	//printf("finish zpos[%lld]=%f\n",iposmark,zpos[iposmark]);
#else
	//printf("finish zpos[%d]=%f\n",iposmark,zpos[iposmark]);
#endif

	
	//getchar();
	inz = iposmark;
	for (integer i = 0; i<inz; i++) zpos[i] = zpos[i + 1]; // сдвиг влево на 1
	SetLength(zpos, inz + 1, inz);
	inz--; // индексация начинается с нуля и заканчивается значением inz

//	for (i = 0; i<adapt_z; i++) simplecorrect_meshgen_z(zpos, inz, lb, ls, lw, b, s, w);

	//for (i = 0; i <= inz; i++) {
		//printf("%e  ", zpos[i]);
	//}
	//getchar();

	integer inz_fix = inz;
	if (b_adhesion_Mesh) {
		// Корректировка источника XY
		for (integer i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranXY != nullptr) {
				if (source_indexpopadaniqnagranXY[i]) {
					doublereal xc = 0.5 * (s[i].g.xS + s[i].g.xE);
					doublereal yc = 0.5 * (s[i].g.yS + s[i].g.yE);
					doublereal zg = s[i].g.zS;
					// найти позицию +Z на сетке zpos
					// найти центр кО.
					// Если он принадлежит Solid block то сместить истоник в цент этой клетки.
					// Адаптацию maxsize ratio 2 не делать.
					// Если неуспех ищемцентр КО по -Z и смещаем туда.
					integer i55_found = -2;
					for (integer i55 = 0; i55 <= inz_fix; i55++) {
						if (fabs(zpos[i55] - zg) < 1.0e-36) {
							i55_found = i55;
							break;
						}
					}
					if (i55_found >= 0) {
						if (i55_found < inz_fix) {
							doublereal zg1 = 0.5 * (zg + zpos[i55_found + 1]);
							printf("zg1=%e\n", zg1);
							integer i56_found = -2;
							for (integer ib55 = 0; ib55 < lb; ib55++) {

								if ((xc > b[ib55].g.xS) && (xc < b[ib55].g.xE) && (yc > b[ib55].g.yS) && (yc < b[ib55].g.yE) && (zg1 > b[ib55].g.zS) && (zg1 < b[ib55].g.zE))
								{
									i56_found = ib55;
								}
							}

							bool bzero_pos = true;
							doublereal zg2 = zpos[0] - 0.5 * fabs(zpos[1] - zpos[0]);
							if (i55_found > 0) {
								zg2 = 0.5 * (zg + zpos[i55_found - 1]);
								bzero_pos = false;
							}
							printf("zg2=%e\n", zg2);
							integer i57_found = -2;
							for (integer ib57 = 0; ib57 < lb; ib57++) {

								if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (zg2 > b[ib57].g.zS) && (zg2 < b[ib57].g.zE))
								{
									if (zg2 > zpos[0]) {
										i57_found = ib57;
									}
								}
							}

							if (i56_found >= 0) {
								if (b[i56_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
									if ((i57_found >= 0) && (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
										// Мы помещаем источник тепла в блок с большей теплопроводностью.
										// comparison_lam выдаёт истину если теплопроводность блока i56 больше.
										if (comparison_lam(matlist, b, i56_found, i57_found, 25.0)) {
											// в блоке i56 теплопроводность выше.
											s[i].g.zS = zg1;
											s[i].g.zE = zg1;
											addboundary(zpos, inz, zg1, XY_PLANE, b, lb, w, lw, s, ls);
										}
										else {
											// в блоке i57 теплопроводность выше.
											s[i].g.zS = zg2;
											s[i].g.zE = zg2;
											addboundary(zpos, inz, zg2, XY_PLANE, b, lb, w, lw, s, ls);
										}
									}
									else {
										// Найден Солид Блок.
										s[i].g.zS = zg1;
										s[i].g.zE = zg1;
										addboundary(zpos, inz, zg1, XY_PLANE, b, lb, w, lw, s, ls);
									}
								}
								else {

									if (i57_found >= 0) {
										if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {

											// Найден Солид Блок.
											printf("zg1==%e\n", zg1);
											doublereal zgolg = zpos[1];
											if (inz == 1) {
												// Слишком малоразмерная расчётная сетка.
												zgolg = 0.5 * (zpos[0] + zg1);
												addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
												zgolg = 0.5 * (zpos[1] + zg1);
												addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
												zgolg = 0.5 * (zpos[0] + zg1);
												// Расстояние от края вглубь расчётной области две клетки.
												s[i].g.zS = zgolg;
												s[i].g.zE = zgolg;
											}
											else {
												if (bzero_pos) {
													// Источник в позиции 0.0 разбиваем поподробнее.
													zgolg = 0.5 * (zpos[0] + zg1);
													addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
													zgolg = 0.5 * (zpos[1] + zg1);
													addboundary(zpos, inz, zgolg, XY_PLANE, b, lb, w, lw, s, ls);
													zgolg = 0.5 * (zpos[0] + zg1);
													// Расстояние от края вглубь расчётной области две клетки.
													s[i].g.zS = zgolg;
													s[i].g.zE = zgolg;
												}
												else {
													s[i].g.zS = zg1;
													s[i].g.zE = zg1;
												}
											}
											addboundary(zpos, inz, zg1, XY_PLANE, b, lb, w, lw, s, ls);


										}
									}
									else {
										printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
										system("PAUSE");
										exit(1);
									}
								}
							}
						}
						else {
							doublereal zg2 = 0.5 * (zg + zpos[i55_found - 1]);
							printf("zg2=%e\n", zg2);
							integer i57_found = -2;
							for (integer ib57 = 0; ib57 < lb; ib57++) {
								if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (zg2 > b[ib57].g.zS) && (zg2 < b[ib57].g.zE))
								{
									i57_found = ib57;
								}
							}
							if (i57_found >= 0) {
								if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
									// Найден Солид Блок.
									s[i].g.zS = zg2;
									s[i].g.zE = zg2;
									addboundary(zpos, inz, zg2, XY_PLANE, b, lb, w, lw, s, ls);
								}
							}
							else {
								printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
								system("PAUSE");
								exit(1);
							}
						}
					}
				}
			}
		}
	}
	//BubbleEnhSort<doublereal>(zpos, 0, inz);
	Sort_method<doublereal>(zpos,inz);
	//SetLength(zpos, inz + 2, inz);
	//inz++;

	/*
	if (source_indexpopadaniqnagranYZ!=nullptr) {
		for (i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranYZ[i] == false) {
				printf("error\n");
				getchar();
			}
		}
	}
	*/

	integer inx_fix = inx;
	if (b_adhesion_Mesh) {
		// Корректировка источника YZ
		for (integer i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranYZ != nullptr) {
				if (source_indexpopadaniqnagranYZ[i]) {
					doublereal zc = 0.5 * (s[i].g.zS + s[i].g.zE);
					doublereal yc = 0.5 * (s[i].g.yS + s[i].g.yE);
					doublereal xg = s[i].g.xS;
					// найти позицию +Z на сетке zpos
					// найти центр кО.
					// Если он принадлежит Solid block то сместить истоник в цент этой клетки.
					// Адаптацию maxsize ratio 2 не делать.
					// Если неуспех ищемцентр КО по -Z и смещаем туда.
					integer i55_found = -2;
					for (integer i55 = 0; i55 <= inx_fix; i55++) {
						if (fabs(xpos[i55] - xg) < 1.0e-36) {
							i55_found = i55;
							break;
						}
					}
					if (i55_found >= 0) {
						if (i55_found < inx_fix) {
							doublereal xg1 = 0.5 * (xg + xpos[i55_found + 1]);
							printf("xg1=%e\n", xg1);
							integer i56_found = -2;
							for (integer ib55 = 0; ib55 < lb; ib55++) {

								if ((zc > b[ib55].g.zS) && (zc < b[ib55].g.zE) && (yc > b[ib55].g.yS) && (yc < b[ib55].g.yE) && (xg1 > b[ib55].g.xS) && (xg1 < b[ib55].g.xE))
								{
									i56_found = ib55;
								}
							}
							bool bzero_pos = true;
							doublereal xg2 = xpos[0] - 0.5 * fabs(xpos[1] - xpos[0]);
							if (i55_found > 0) {
								xg2 = 0.5 * (xg + xpos[i55_found - 1]);
								bzero_pos = false;
							}
							printf("xg2=%e\n", xg2);
							integer i57_found = -2;
							for (integer ib57 = 0; ib57 < lb; ib57++) {
								if ((zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (xg2 > b[ib57].g.xS) && (xg2 < b[ib57].g.xE))
								{
									if (xg2 > xpos[0]) {
										i57_found = ib57;
									}
								}

							}



							if (i56_found >= 0) {
								if (b[i56_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
									// TODO 11.07.2016
									if ((i57_found >= 0) && (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
										// Мы помещаем источник тепла в блок с большей теплопроводностью.
										// comparison_lam выдаёт истину если теплопроводность блока i56 больше.
										if (comparison_lam(matlist, b, i56_found, i57_found, 25.0)) {
											// в блоке i56 теплопроводность выше.
											s[i].g.xS = xg1;
											s[i].g.xE = xg1;
											addboundary(xpos, inz, xg1, YZ_PLANE, b, lb, w, lw, s, ls);
										}
										else {
											// в блоке i57 теплопроводность выше.
											s[i].g.xS = xg2;
											s[i].g.xE = xg2;
											addboundary(xpos, inz, xg2, YZ_PLANE, b, lb, w, lw, s, ls);
										}
									}
									else {
										// Найден Солид Блок.
										printf("xg1==%e\n", xg1);
										doublereal xgolg = xpos[1];
										if (inx == 1) {
											// Слишком малоразмерная расчётная сетка.
											xgolg = 0.5 * (xpos[0] + xg1);
											addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
											xgolg = 0.5 * (xpos[1] + xg1);
											addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
											xgolg = 0.5 * (xpos[0] + xg1);
											// Расстояние от края вглубь расчётной области две клетки.
											s[i].g.xS = xgolg;
											s[i].g.xE = xgolg;
										}
										else {
											if (bzero_pos) {
												// Источник в позиции 0.0.
												// Разбиваем поподробнее.
												xgolg = 0.5 * (xpos[0] + xg1);
												addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
												xgolg = 0.5 * (xpos[1] + xg1);
												addboundary(xpos, inx, xgolg, YZ_PLANE, b, lb, w, lw, s, ls);
												xgolg = 0.5 * (xpos[0] + xg1);
												// Расстояние от края вглубь расчётной области две клетки.
												s[i].g.xS = xgolg;
												s[i].g.xE = xgolg;
											}
											else {
												s[i].g.xS = xg1;
												s[i].g.xE = xg1;
											}
										}
										addboundary(xpos, inx, xg1, YZ_PLANE, b, lb, w, lw, s, ls);


										//printf("ok");
									}
								}
								else {

									if (i57_found >= 0) {
										if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
											// Найден Солид Блок.
											s[i].g.xS = xg2;
											s[i].g.xE = xg2;
											addboundary(xpos, inx, xg2, YZ_PLANE, b, lb, w, lw, s, ls);
										}
									}
									else {
										printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
										system("PAUSE");
										exit(1);
									}
								}
							}
						}
						else {
							doublereal xg2 = 0.5 * (xg + xpos[i55_found - 1]);
							printf("xg2=%e\n", xg2);
							integer i57_found = -2;
							for (integer ib57 = 0; ib57 < lb; ib57++) {
								if ((zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yc > b[ib57].g.yS) && (yc < b[ib57].g.yE) && (xg2 > b[ib57].g.xS) && (xg2 < b[ib57].g.xE))
								{
									i57_found = ib57;
								}
							}
							if (i57_found >= 0) {
								if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
									// Найден Солид Блок.
									s[i].g.xS = xg2;
									s[i].g.xE = xg2;
									addboundary(xpos, inx, xg2, YZ_PLANE, b, lb, w, lw, s, ls);
								}
							}
							else {
								printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
								system("PAUSE");
								exit(1);
							}
						}
					}
				}
			}
		}
	}
	//getchar();
	//BubbleEnhSort<doublereal>(xpos, 0, inx);
	Sort_method<doublereal>(xpos,inx);

	/*
	for (i = 0; i <= inx; i++) {
		printf("%e\n", xpos[i]);
	}
	printf("%e %e %e\n", s[0].g.xS, s[1].g.xS, s[2].g.xS);
	getchar();
	*/

	integer iny_fix = iny;
	if (b_adhesion_Mesh) {
		// Корректировка источника XZ
		for (integer i = 0; i < ls; i++) {
			if (source_indexpopadaniqnagranXZ != nullptr) {
				if (source_indexpopadaniqnagranXZ[i]) {
					doublereal xc = 0.5 * (s[i].g.xS + s[i].g.xE);
					doublereal zc = 0.5 * (s[i].g.zS + s[i].g.zE);
					doublereal yg = s[i].g.yS;
					// найти позицию +Z на сетке zpos
					// найти центр кО.
					// Если он принадлежит Solid block то сместить истоник в цент этой клетки.
					// Адаптацию maxsize ratio 2 не делать.
					// Если неуспех ищемцентр КО по -Z и смещаем туда.
					integer i55_found = -2;
					for (integer i55 = 0; i55 <= iny_fix; i55++) {
						if (fabs(ypos[i55] - yg) < shorter_length_for_simplificationY(yg, b, lb, w, lw, s, ls)) {
							i55_found = i55;
							break;
						}
					}
					if (i55_found >= 0) {
						if (i55_found < iny_fix) {
							doublereal yg1 = 0.5 * (yg + ypos[i55_found + 1]);
							//printf("yg1=%e\n", yg1);
							std::cout << "yg1=" << yg1 << std::endl;
							integer i56_found = -2;
							for (integer ib55 = 0; ib55 < lb; ib55++) {
								if ((xc > b[ib55].g.xS) && (xc < b[ib55].g.xE) && (zc > b[ib55].g.zS) && (zc < b[ib55].g.zE) && (yg1 > b[ib55].g.yS) && (yg1 < b[ib55].g.yE))
								{
									i56_found = ib55;
								}
							}

							bool bzero_pos = true;
							doublereal yg2 = ypos[0] - 0.5 * fabs(ypos[1] - ypos[0]);
							if (i55_found > 0) {
								yg2 = 0.5 * (yg + ypos[i55_found - 1]);
								bzero_pos = false;
							}

							//printf("yg2=%e\n", yg2);
							std::cout << "yg2=" << yg2 << std::endl;
							integer i57_found = -2;
							for (integer ib57 = 0; ib57 < lb; ib57++) {
								if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yg2 > b[ib57].g.yS) && (yg2 < b[ib57].g.yE))
								{
									if (yg2 > ypos[0]) {
										i57_found = ib57;
									}
								}
							}


							if (i56_found >= 0) {
								if (b[i56_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {

									if ((i57_found >= 0) && (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID)) {
										// Мы помещаем источник тепла в блок с большей теплопроводностью.
										// comparison_lam выдаёт истину если теплопроводность блока i56 больше.
										if (comparison_lam(matlist, b, i56_found, i57_found, 25.0)) {
											// в блоке i56 теплопроводность выше.
											s[i].g.yS = yg1;
											s[i].g.yE = yg1;
											addboundary(ypos, iny, yg1, XZ_PLANE, b, lb, w, lw, s, ls);
										}
										else {
											// в блоке i57 теплопроводность выше.
											s[i].g.yS = yg2;
											s[i].g.yE = yg2;
											addboundary(ypos, iny, yg2, XZ_PLANE, b, lb, w, lw, s, ls);
										}
									}
									else {

										// Найден Солид Блок.
										//printf("yg1==%e\n", yg1);
										std::cout << "yg1=" << yg1 << std::endl;
										doublereal ygolg = ypos[1];
										if (iny == 1) {
											// Слишком малоразмерная расчётная сетка.
											ygolg = 0.5 * (ypos[0] + yg1);
											addboundary(ypos, iny, ygolg, XZ_PLANE, b, lb, w, lw, s, ls);
											ygolg = 0.5 * (ypos[1] + yg1);
											addboundary(ypos, iny, ygolg, XZ_PLANE, b, lb, w, lw, s, ls);
											ygolg = 0.5 * (ypos[0] + yg1);
											// Расстояние от края вглубь расчётной области две клетки.
											s[i].g.yS = ygolg;
											s[i].g.yE = ygolg;
										}
										else {
											if (bzero_pos) {
												// Источник в позиции 0.0.
												// Разбиваем поподробнее.
												ygolg = 0.5 * (ypos[0] + yg1);
												addboundary(ypos, iny, ygolg, XZ_PLANE, b, lb, w, lw, s, ls);
												ygolg = 0.5 * (ypos[1] + yg1);
												addboundary(ypos, iny, ygolg, XZ_PLANE, b, lb, w, lw, s, ls);
												ygolg = 0.5 * (ypos[0] + yg1);
												s[i].g.yS = ygolg;
												s[i].g.yE = ygolg;
											}
											else {
												s[i].g.yS = yg1;
												s[i].g.yE = yg1;
											}
										}
										addboundary(ypos, iny, yg1, XZ_PLANE, b, lb, w, lw, s, ls);


									}
								}
								else {

									if (i57_found >= 0) {
										if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
											// Найден Солид Блок.
											s[i].g.yS = yg2;
											s[i].g.yE = yg2;
											addboundary(ypos, iny, yg2, XZ_PLANE, b, lb, w, lw, s, ls);
										}
									}
									else {
										printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
										system("PAUSE");
										exit(1);
									}
								}
							}
						}
						else {
							doublereal yg2 = 0.5 * (yg + ypos[i55_found - 1]);
							//printf("yg2=%e\n", yg2);
							std::cout << "yg2=" << yg2 << std::endl;
							integer i57_found = -2;
							for (integer ib57 = 0; ib57 < lb; ib57++) {
								if ((xc > b[ib57].g.xS) && (xc < b[ib57].g.xE) && (zc > b[ib57].g.zS) && (zc < b[ib57].g.zE) && (yg2 > b[ib57].g.yS) && (yg2 < b[ib57].g.yE))
								{
									i57_found = ib57;
								}
							}
							if (i57_found >= 0) {
								if (b[i57_found].itype == PHYSICS_TYPE_IN_BODY::SOLID) {
									// Найден Солид Блок.
									s[i].g.yS = yg2;
									s[i].g.yE = yg2;
									addboundary(ypos, iny, yg2, XZ_PLANE, b, lb, w, lw, s, ls);
								}
							}
							else {
								printf("ERROR: sourse na granice dvus hollow or fluid blockov.");
								system("PAUSE");
								exit(1);
							}
						}
					}
				}
			}
		}
	}
	//BubbleEnhSort<doublereal>(ypos, 0, iny);
	Sort_method<doublereal>(ypos,iny);
	// debug
	
	//for (i=0; i<=inz; i++) {
	  // printf("%e  ",zpos[i]);
	//}
	//getchar();//
	
	// 3 сентября 2017.
	// Добавление новых сеточных линий.
	// необходимость их добавления обусловлена требованнем АЛИС сеток.
	for (int i_28 = 0; i_28 <= inzadd; i_28++) {
		//SetLength(zpos, inz+1, inz + 2);
		//zpos[inz+1] = zposadd[i_28];
		//inz++;
		// Это добавление умное, препятствует возникновению одинаковых значений в массиве.
		addboundary(zpos, inz, zposadd[i_28],XY_PLANE, b, lb, w, lw, s, ls);
	}
	Sort_method<doublereal>(zpos,inz);


	if (1) {
		integer inum_iter_ratio_good = 0;
		doublereal max_size_ratio_z = 1.0;
		while (1) {

			max_size_ratio_z = 1.0;
			for (integer i = 0; i<inz - 1; i++) {

				doublereal dmax = 0.0;
				doublereal dmin = 0.0;
				if (fabs(zpos[i + 1] - zpos[i])>fabs(zpos[i + 2] - zpos[i + 1])) {
					dmax = fabs(zpos[i + 1] - zpos[i]);
					dmin = fabs(zpos[i + 2] - zpos[i + 1]);
				}
				else {
					dmax = fabs(zpos[i + 2] - zpos[i + 1]);
					dmin = fabs(zpos[i + 1] - zpos[i]);
				}
				if (dmax / dmin>max_size_ratio_z) {
					max_size_ratio_z = dmax / dmin;
				}
			}
			
			if (max_size_ratio_z != max_size_ratio_z) {
				//printf("z axis max size ratio is equal = %1.4f\n", max_size_ratio_z);
				std::cout << "z axis max size ratio is equal = " << max_size_ratio_z << std::endl;
				system("PAUSE");
			}
			//getchar();

			if (max_size_ratio_z < (etalon_max_size_ratio*1.1)) {
				break;
			}

			// Корректировка max size ratio z axis.
			// сделано 9 апреля 2013 года.

			bool bplus = false;
			max_size_ratio_z = 1.0;
			for (integer i = 0; i<inz - 1; i++) {

				bplus = false;
				doublereal dmax = 0.0;
				doublereal dmin = 0.0;
				if (fabs(zpos[i + 1] - zpos[i])>fabs(zpos[i + 2] - zpos[i + 1])) {
					dmax = fabs(zpos[i + 1] - zpos[i]);
					dmin = fabs(zpos[i + 2] - zpos[i + 1]);
				}
				else {
					bplus = true;
					dmax = fabs(zpos[i + 2] - zpos[i + 1]);
					dmin = fabs(zpos[i + 1] - zpos[i]);
				}
				if (dmax / dmin>etalon_max_size_ratio) {
					doublereal pos_candidate;
					if (bplus) {
						pos_candidate = zpos[i + 1] + 0.5*dmax;
					}
					else {
						pos_candidate = zpos[i] + 0.5*dmax;
					}
					SetLength(zpos, inz + 1, inz + 2);
					zpos[inz + 1] = pos_candidate;
					inz = inz + 1;
					//BubbleEnhSort<doublereal>(zpos, 0, inz);
					Sort_method<doublereal>(zpos,inz);
					break;
				}
			}

			inum_iter_ratio_good++;

			//getchar();

		}

		//printf("z axis max size ratio is equal = %1.4f\n", max_size_ratio_z);
		std::cout << "z axis max size ratio is equal = " << max_size_ratio_z << std::endl;

#if doubleintprecision == 1
		printf("inum_iter_ratio_good is %lld\n", inum_iter_ratio_good);
#else
		printf("inum_iter_ratio_good is %d\n", inum_iter_ratio_good);
#endif

		
	}

	// улучшение качества сетки по второму критерию с значением 30.0
	// как во flowvision.
	quolite_refinement(inx, iny, inz, xpos, ypos, zpos);

	patch_mesh_refinement_21_11_2019(inx, iny, inz,
		xpos, ypos, zpos, lb, b);

	


	// Освобождение оперативной памяти.
	delete[] ib_marker;
	delete[] ib_marker_flag_power;
	delete[] ib_marker_flag_fluid;
	delete[] rxboundary;
	delete[] ryboundary;
	delete[] rzboundary;
	delete[] ixintervalcount;
	delete[] iyintervalcount;
	delete[] izintervalcount;
	delete[] source_indexpopadaniqnagranXY;
	delete[] source_indexpopadaniqnagranYZ;
	delete[] source_indexpopadaniqnagranXZ;

	
} // coarsemeshgen2


// 28.06.2020
void cad_geometry_octree_meshgen(doublereal*& xpos, doublereal*& ypos, doublereal*& zpos,
	int& inx, int& iny, int& inz,
	int lb, int ls, int lw, 
	BLOCK*& b, SOURCE*& s, WALL*& w, int lu, UNION*& my_union, TPROP* matlist,
	doublereal*& xposadd, doublereal*& yposadd, doublereal*& zposadd,
	int& inxadd, int& inyadd, int& inzadd, integer& iunion_id_p1)
{

	doublereal base_min_side = fmin(b[0].g.xE - b[0].g.xS,fmin(b[0].g.yE - b[0].g.yS, b[0].g.zE - b[0].g.zS));
	LINE_DIRECTIONAL side_id;
	if (fabs(base_min_side - fabs(b[0].g.xE - b[0].g.xS)) < 1.0e-30) side_id = LINE_DIRECTIONAL::X_LINE_DIRECTIONAL;
	if (fabs(base_min_side - fabs(b[0].g.yE - b[0].g.yS)) < 1.0e-30) side_id = LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL;
	if (fabs(base_min_side - fabs(b[0].g.zE - b[0].g.zS)) < 1.0e-30) side_id = LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL;

	int ing = (int)(fmax(inx, fmax(iny, inz)));
	doublereal cell_size = base_min_side / ing;

	switch (side_id) {
	case LINE_DIRECTIONAL::X_LINE_DIRECTIONAL:
		std::cout << "X" << std::endl;
		inx = ing;
		xpos = new doublereal[ing + 1];
		xpos[0] = b[0].g.xS;
		for (integer i = 1; i <= ing; i++) xpos[i] = b[0].g.xS + i * cell_size;
		xpos[ing] = b[0].g.xE;
		ing = (int)(fabs(b[0].g.yE - b[0].g.yS)/cell_size);
		cell_size = fabs(b[0].g.yE - b[0].g.yS) / ing;
		iny = ing;
		ypos = new doublereal[ing + 1];
		ypos[0] = b[0].g.yS;
		for (integer i = 1; i <= ing; i++) ypos[i] = b[0].g.yS + i * cell_size;
		ypos[ing] = b[0].g.yE;
		ing = (int)(fabs(b[0].g.zE - b[0].g.zS) / cell_size);
		cell_size = fabs(b[0].g.zE - b[0].g.zS) / ing;
		iny = ing;
		zpos = new doublereal[ing + 1];
		zpos[0] = b[0].g.zS;
		for (integer i = 1; i <= ing; i++) zpos[i] = b[0].g.zS + i * cell_size;
		zpos[ing] = b[0].g.zE;
		break;
	case LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL:
		std::cout << "Y" << std::endl;
		iny = ing;
		ypos = new doublereal[ing + 1];
		ypos[0] = b[0].g.yS;
		for (integer i = 1; i <= ing; i++) ypos[i] = b[0].g.yS + i * cell_size;
		ypos[ing] = b[0].g.yE;
		ing = (int)(fabs(b[0].g.xE - b[0].g.xS) / cell_size);
		cell_size = fabs(b[0].g.xE - b[0].g.xS) / ing;
		inx = ing;
		xpos = new doublereal[ing + 1];
		xpos[0] = b[0].g.xS;
		for (integer i = 1; i <= ing; i++) xpos[i] = b[0].g.xS + i * cell_size;
		xpos[ing] = b[0].g.xE;
		ing = (int)(fabs(b[0].g.zE - b[0].g.zS) / cell_size);
		cell_size = fabs(b[0].g.zE - b[0].g.zS) / ing;
		iny = ing;
		zpos = new doublereal[ing + 1];
		zpos[0] = b[0].g.zS;
		for (integer i = 1; i <= ing; i++) zpos[i] = b[0].g.zS + i * cell_size;
		zpos[ing] = b[0].g.zE;
		break;
	case LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL:
		std::cout << "Z" << std::endl;
		inz = ing;
		zpos = new doublereal[ing + 1];
		zpos[0] = b[0].g.zS;
		for (integer i = 1; i <= ing; i++) zpos[i] = b[0].g.zS + i * cell_size;
		if (zpos[ing - 1] > b[0].g.zE) {
			std::cout << zpos[ing - 1] << " " << b[0].g.zE << " " << cell_size << std::endl;
		}
		zpos[ing] = b[0].g.zE;
		ing = (int)(fabs(b[0].g.xE - b[0].g.xS) / cell_size);
		cell_size = fabs(b[0].g.xE - b[0].g.xS) / ing;
		inx = ing;
		xpos = new doublereal[ing + 1];
		xpos[0] = b[0].g.xS;
		for (integer i = 1; i <= ing; i++) xpos[i] = b[0].g.xS + i * cell_size;
		if (xpos[ing - 1] > b[0].g.xE) {
			std::cout << xpos[ing - 1] << " " << b[0].g.xE << " " << cell_size << std::endl;
		}
		xpos[ing] = b[0].g.xE;
		ing = (int)(fabs(b[0].g.yE - b[0].g.yS) / cell_size);
		cell_size = fabs(b[0].g.yE - b[0].g.yS) / ing;
		iny = ing;
		ypos = new doublereal[ing + 1];
		ypos[0] = b[0].g.yS;
		for (integer i = 1; i <= ing; i++) ypos[i] = b[0].g.yS + i * cell_size;
		if (ypos[ing - 1] > b[0].g.yE) {
			std::cout << ypos[ing - 1] << " " << b[0].g.yE << " " << cell_size << std::endl;
		}
		ypos[ing] = b[0].g.yE;
		break;
	}

	// Добавление новых сеточных линий не производится.

}// cad_geometry_octree_meshgen


void coarsemeshgen(doublereal*& xpos, doublereal*& ypos, doublereal*& zpos,
	int& inx, int& iny, int& inz,
	int lb, int ls, int lw,
	BLOCK*& b, SOURCE*& s, WALL*& w, int lu, UNION*& my_union, TPROP* matlist,
	doublereal*& xposadd, doublereal*& yposadd, doublereal*& zposadd,
	int& inxadd, int& inyadd, int& inzadd, integer& iunion_id_p1)
{
	if (CAD_GEOMETRY_OCTREE_MESHGEN) {
		cad_geometry_octree_meshgen(xpos, ypos, zpos,
			inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
			xposadd, yposadd, zposadd,
			inxadd, inyadd, inzadd, iunion_id_p1);
	}
	else {
		coarsemeshgen2(xpos, ypos, zpos,
			inx, iny, inz, lb, ls, lw, b, s, w, lu, my_union, matlist,
			xposadd, yposadd, zposadd,
			inxadd, inyadd, inzadd, iunion_id_p1);
	}
}// coarsemeshgen

#endif
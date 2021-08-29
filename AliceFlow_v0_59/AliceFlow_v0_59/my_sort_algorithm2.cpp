// Алгоритмы сортировки
// Используются в алгебраическом многосеточном методе.
#pragma once
#ifndef _MY_SORT_ALGORITHM2_CPP_
#define _MY_SORT_ALGORITHM2_CPP_ 1

void swap(Ak2& Amat, integer i, integer j)
{
	Ak1 A_temp;

	// change Amat[i] <-> Amat[j]
	A_temp.i = Amat.i[i];
	Amat.i[i] = Amat.i[j];
	Amat.i[j] = A_temp.i;

	A_temp.j = Amat.j[i];
	Amat.j[i] = Amat.j[j];
	Amat.j[j] = A_temp.j;

	A_temp.aij = Amat.aij[i];
	Amat.aij[i] = Amat.aij[j];
	Amat.aij[j] = A_temp.aij;

} // swap


bool comparei(const Ak2& Amat, integer index, const Ak1& Amat2) {
	return ((Amat.i[index] < Amat2.i) || ((Amat.i[index] == Amat2.i) && (Amat.j[index] < Amat2.j)));
}

bool comparej(const Ak2& Amat, integer index, const Ak1& Amat2) {
	return ((Amat.j[index] < Amat2.j) || ((Amat.j[index] == Amat2.j) && (Amat.i[index] < Amat2.i)));
}

// PivotList содержит ошибку.Обнаружено 12 декабря 2015.
// Вот алгоритм PivotList
integer PivotList(Ak2& Amat, integer first, integer last,
	bool (*compare)(Ak2& Amat, integer index, Ak1& Amat2)) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	Ak1 Afirst;
	Afirst.i= Amat.i[first];
	Afirst.j = Amat.j[first]; 
	integer PivotPoint = first;

	for (integer index_1 = (first + 1); index_1 <= last; index_1++) {
		if (compare(Amat, index_1, Afirst)) {
			PivotPoint++;
			swap(Amat, PivotPoint, index_1);
		}
	}

	swap(Amat, first, PivotPoint);

	return PivotPoint;
} // PivotList


// Сортировка Тима (Тим Петерсон).
// см. файл uniformsimplemeshgen.cpp
//const integer RUN = 32;

// данная функция сортирует массив начиная с левого индекса left
// до правого индекса right который имеет размер не более RUN
void insertionSortTim_amg(Ak2& Amat, integer left, integer right)
{
	for (integer i = left + 1; i <= right; i++)
	{
		Ak1 temp;
		temp.i = Amat.i[i];
		temp.j = Amat.j[i];
		temp.aij = Amat.aij[i];

		integer j = i - 1;
		while ((j >= left) && (Amat.i[j] > temp.i))
		{
			Amat.i[j + 1] = Amat.i[j];
			Amat.j[j + 1] = Amat.j[j];
			Amat.aij[j + 1] = Amat.aij[j];
			j--;
		}
		Amat.i[j + 1] = temp.i;
		Amat.j[j + 1] = temp.j;
		Amat.aij[j + 1] = temp.aij;
	}
}// insertion sort по i.

// данная функция сортирует массив начиная с левого индекса left
// до правого индекса right который имеет размер не более RUN
void insertionSortTim_amg_j(Ak2& Amat, integer left, integer right)
{
	for (integer i = left + 1; i <= right; i++)
	{
		Ak1 temp;
		temp.i = Amat.i[i];
		temp.j = Amat.j[i];
		temp.aij = Amat.aij[i];

		integer j = i - 1;
		while ((j >= left) && (Amat.j[j] > temp.j))
		{
			Amat.i[j + 1] = Amat.i[j];
			Amat.j[j + 1] = Amat.j[j];
			Amat.aij[j + 1] = Amat.aij[j];
			j--;
		}
		Amat.i[j + 1] = temp.i;
		Amat.j[j + 1] = temp.j;
		Amat.aij[j + 1] = temp.aij;
	}
} // insertion sort по j.

// данная функция сортирует массив начиная с левого индекса left
// до правого индекса right который имеет размер не более RUN
void insertionSortTim_amg(Ak1*& Amat, integer left, integer right)
{
	for (integer i = left + 1; i <= right; i++)
	{
		Ak1 temp;
		temp = Amat[i];


		integer j = i - 1;
		while ((j >= left) && (Amat[j].i > temp.i))
		{
			Amat[j + 1] = Amat[j];
			j--;
		}
		Amat[j + 1] = temp;
	}
}// insertion sort по i.

// данная функция сортирует массив начиная с левого индекса left
// до правого индекса right который имеет размер не более RUN
void insertionSortTim_amg_j(Ak1*& Amat, integer left, integer right)
{
	for (integer i = left + 1; i <= right; i++)
	{
		Ak1 temp;
		temp = Amat[i];


		integer j = i - 1;
		while ((j >= left) && (Amat[j].j > temp.j))
		{
			Amat[j + 1] = Amat[j];
			j--;
		}
		Amat[j + 1] = temp;
	}
} // insertion sort по j.

// merge function merges the sorted runs
void mergeTim_amg(Ak2& Amat, integer l, integer m, integer r)
{
	// original array is broken in two parts
	// left and right array
	integer len1 = m - l + 1, len2 = r - m;

	//myARRT left[len1], right[len2];
	Ak2 left;
	left.i = nullptr;
	left.j = nullptr;
	left.aij = nullptr;
	if (len1 >= 0) {
		left.i = new integer_mix_precision[len1 + 1];
		left.j = new integer_mix_precision[len1 + 1];
		left.aij = new real_mix_precision[len1 + 1];
	}
	Ak2 right;
	right.i = nullptr;
	right.j = nullptr;
	right.aij = nullptr;
	if (len2 >= 0) {
		right.i = new integer_mix_precision[len2 + 1];
		right.j = new integer_mix_precision[len2 + 1];
		right.aij = new real_mix_precision[len2 + 1];
	}

	for (integer i = 0; i < len1; i++) {
		left.i[i] = Amat.i[l + i];
		left.j[i] = Amat.j[l + i];
		left.aij[i] = Amat.aij[l + i];

	}
	for (integer i = 0; i < len2; i++) {
		right.i[i] = Amat.i[m + 1 + i];
		right.j[i] = Amat.j[m + 1 + i];
		right.aij[i] = Amat.aij[m + 1 + i];
	}

	integer i = 0;
	integer j = 0;
	integer k = l;

	// after comparing, we merge those two array
	// in larger sub array
	while (i < len1 && j < len2)
	{
		if (left.i[i] <= right.i[j])
		{
			Amat.i[k] = left.i[i];
			Amat.j[k] = left.j[i];
			Amat.aij[k] = left.aij[i];
			i++;
		}
		else
		{
			Amat.i[k] = right.i[j];
			Amat.j[k] = right.j[j];
			Amat.aij[k] = right.aij[j];
			j++;
		}
		k++;
	}

	// copy remaining elements of left, if any
	while (i < len1)
	{
		Amat.i[k] = left.i[i];
		Amat.j[k] = left.j[i];
		Amat.aij[k] = left.aij[i];
		k++;
		i++;
	}

	// copy remaining element of right, if any
	while (j < len2)
	{
		Amat.i[k] = right.i[j];
		Amat.j[k] = right.j[j];
		Amat.aij[k] = right.aij[j];
		k++;
		j++;
	}

	if (left.i != nullptr) {
		delete[] left.i;
	}
	if (left.j != nullptr) {
		delete[] left.j;
	}
	if (left.aij != nullptr) {
		delete[] left.aij;
	}
	if (right.i != nullptr) {
		delete[] right.i;
	}
	if (right.j != nullptr) {
		delete[] right.j;
	}
	if (right.aij != nullptr) {
		delete[] right.aij;
	}
}


// merge function merges the sorted runs
void mergeTim_amg_j(Ak2& Amat, integer l, integer m, integer r)
{
	// original array is broken in two parts
	// left and right array
	integer len1 = m - l + 1, len2 = r - m;

	//myARRT left[len1], right[len2];
	Ak2 left;
	left.i = nullptr;
	left.j = nullptr;
	left.aij = nullptr;
	if (len1 >= 0) {
		left.i = new integer_mix_precision[len1 + 1];
		left.j = new integer_mix_precision[len1 + 1];
		left.aij = new real_mix_precision[len1 + 1];
	}
	Ak2 right;
	right.i = nullptr;
	right.j = nullptr;
	right.aij = nullptr;
	if (len2 >= 0) {
		right.i = new integer_mix_precision[len2 + 1];
		right.j = new integer_mix_precision[len2 + 1];
		right.aij = new real_mix_precision[len2 + 1];
	}

	for (integer i = 0; i < len1; i++) {
		left.i[i] = Amat.i[l + i];
		left.j[i] = Amat.j[l + i];
		left.aij[i] = Amat.aij[l + i];

	}
	for (integer i = 0; i < len2; i++) {
		right.i[i] = Amat.i[m + 1 + i];
		right.j[i] = Amat.j[m + 1 + i];
		right.aij[i] = Amat.aij[m + 1 + i];
	}

	integer i = 0;
	integer j = 0;
	integer k = l;

	// after comparing, we merge those two array
	// in larger sub array
	while (i < len1 && j < len2)
	{
		if (left.j[i] <= right.j[j])
		{
			Amat.i[k] = left.i[i];
			Amat.j[k] = left.j[i];
			Amat.aij[k] = left.aij[i];
			i++;
		}
		else
		{
			Amat.i[k] = right.i[j];
			Amat.j[k] = right.j[j];
			Amat.aij[k] = right.aij[j];
			j++;
		}
		k++;
	}

	// copy remaining elements of left, if any
	while (i < len1)
	{
		Amat.i[k] = left.i[i];
		Amat.j[k] = left.j[i];
		Amat.aij[k] = left.aij[i];
		k++;
		i++;
	}

	// copy remaining element of right, if any
	while (j < len2)
	{
		Amat.i[k] = right.i[j];
		Amat.j[k] = right.j[j];
		Amat.aij[k] = right.aij[j];
		k++;
		j++;
	}

	if (left.i != nullptr) {
		delete[] left.i;
	}
	if (left.j != nullptr) {
		delete[] left.j;
	}
	if (left.aij != nullptr) {
		delete[] left.aij;
	}
	if (right.i != nullptr) {
		delete[] right.i;
	}
	if (right.j != nullptr) {
		delete[] right.j;
	}
	if (right.aij != nullptr) {
		delete[] right.aij;
	}
}

// merge function merges the sorted runs
void mergeTim_amg(Ak1*& Amat, integer l, integer m, integer r,
	integer (*indx_compare)(Ak1& Amat))
{
	// original array is broken in two parts
	// left and right array
	integer len1 = m - l + 1, len2 = r - m;

	//myARRT left[len1], right[len2];
	Ak1* left = nullptr;
	if (len1 >= 0) {
		left = new Ak1[len1 + 1];
	}
	Ak1* right = nullptr;
	if (len2 >= 0) {
		right = new Ak1[len2 + 1];
	}

	for (integer i = 0; i < len1; i++) {
		left[i] = Amat[l + i];
	}
	for (integer i = 0; i < len2; i++) {
		right[i] = Amat[m + 1 + i];
	}

	integer i = 0;
	integer j = 0;
	integer k = l;

	// after comparing, we merge those two array
	// in larger sub array
	while (i < len1 && j < len2)
	{
		if (indx_compare(left[i]) <= indx_compare(right[j]))
		{
			Amat[k] = left[i];
			i++;
		}
		else
		{
			Amat[k] = right[j];
			j++;
		}
		k++;
	}

	// copy remaining elements of left, if any
	while (i < len1)
	{
		Amat[k] = left[i];
		k++;
		i++;
	}

	// copy remaining element of right, if any
	while (j < len2)
	{
		Amat[k] = right[j];
		k++;
		j++;
	}

	if (left != nullptr) {
		delete[] left;
	}
	if (right != nullptr) {
		delete[] right;
	}

}

//С помощью неполной НЕУБЫВАЮЩЕЙ кучи 
//крупные элементы закидываем поближе к концу массива 
void reheap(Ak1*& a, integer length, integer i,
	integer (*indx_compare)(Ak1& Amat)) {

	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	Ak1 Temp = a[i];
	integer parent = i;
	integer child = 2 * (i + 1) - 1;
	//Просматриваем потомков, а также потомков потомков 
	//и сравниваем их с родителем (если что - передвигаем потомков влево) //Цикл продолжается пока не выпадем за пределы массива
	 //или пока не обменяем какого-нибудь потомка на родителя. 
	while ((child < length) && (!done)) {
		//Если правый потомок в пределах массива 
		if (child < length - 1) {
			//То из левого и правого потомка выбираем наименьшего 
			if (indx_compare(a[child]) >= indx_compare(a[child + 1])) { child += 1; }
		}
		//Родитель меньше потомков? 
		if (indx_compare(Temp) < indx_compare(a[child])) {
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


//С помощью неполной НЕУБЫВАЮЩЕЙ кучи 
//крупные элементы закидываем поближе к концу массива 
void reheap_i(Ak2& a, integer length, integer i) {
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	Ak1 Temp;
	Temp.i = a.i[i];
	Temp.j = a.j[i];
	Temp.aij = a.aij[i];

	integer parent = i;
	integer child = 2 * (i + 1) - 1;
	//Просматриваем потомков, а также потомков потомков 
	//и сравниваем их с родителем (если что - передвигаем потомков влево) //Цикл продолжается пока не выпадем за пределы массива
	 //или пока не обменяем какого-нибудь потомка на родителя. 
	while ((child < length) && (!done)) {
		//Если правый потомок в пределах массива 
		if (child < length - 1) {
			//То из левого и правого потомка выбираем наименьшего 
			if (a.i[child] >= a.i[child + 1]) { child += 1; }
		}
		//Родитель меньше потомков? 
		if (Temp.i < a.i[child]) {
			//Тогда с этим родителем и его потомками разобрались
			done = true;
			//Родитель НЕ меньше чем наименьший из его потомков. 
			//Перемещаем потомка на место родителя 
			//и с родителем в цикле сравниваем уже потомков этого потомка 
		}
		else
		{
			a.i[parent] = a.i[child];
			a.j[parent] = a.j[child];
			a.aij[parent] = a.aij[child];
			parent = child;
			child = 2 * (parent + 1) - 1;
		}
	}
	//Родитель, с которого всё начиналось 
	//передвигается ближе к концу массива 
	//(или остаётся на месте если не повезло) 
	a.i[parent] = Temp.i;
	a.j[parent] = Temp.j;
	a.aij[parent] = Temp.aij;
}

//С помощью неполной НЕУБЫВАЮЩЕЙ кучи 
//крупные элементы закидываем поближе к концу массива 
void reheap_j(Ak2& a, integer length, integer i) {
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	Ak1 Temp;
	Temp.i = a.i[i];
	Temp.j = a.j[i];
	Temp.aij = a.aij[i];

	integer parent = i;
	integer child = 2 * (i + 1) - 1;
	//Просматриваем потомков, а также потомков потомков 
	//и сравниваем их с родителем (если что - передвигаем потомков влево) //Цикл продолжается пока не выпадем за пределы массива
	 //или пока не обменяем какого-нибудь потомка на родителя. 
	while ((child < length) && (!done)) {
		//Если правый потомок в пределах массива 
		if (child < length - 1) {
			//То из левого и правого потомка выбираем наименьшего 
			if (a.j[child] >= a.j[child + 1]) { child += 1; }
		}
		//Родитель меньше потомков? 
		if (Temp.j < a.j[child]) {
			//Тогда с этим родителем и его потомками разобрались
			done = true;
			//Родитель НЕ меньше чем наименьший из его потомков. 
			//Перемещаем потомка на место родителя 
			//и с родителем в цикле сравниваем уже потомков этого потомка 
		}
		else
		{
			a.i[parent] = a.i[child];
			a.j[parent] = a.j[child];
			a.aij[parent] = a.aij[child];
			parent = child;
			child = 2 * (parent + 1) - 1;
		}
	}
	//Родитель, с которого всё начиналось 
	//передвигается ближе к концу массива 
	//(или остаётся на месте если не повезло) 
	a.i[parent] = Temp.i;
	a.j[parent] = Temp.j;
	a.aij[parent] = Temp.aij;
}

//С помощью неполной НЕВОЗРАСТАЮЩЕЙ кучи 
//мелкие элементы закидываем поближе к началу массива 
void invreheap(Ak1*& a, integer length, integer i,
	integer (*indx_compare)(Ak1& Amat))
{
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	Ak1 Temp = a[length - 1 - i];
	integer parent = i;
	integer child = 2 * (i + 1) - 1;
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
			if (indx_compare(a[length - 1 - child]) <= indx_compare(a[length - 1 - (child + 1)]))
			{
				child += 1;
			}
		}
		//Родитель больше потомков? 
		if (indx_compare(Temp) > indx_compare(a[length - 1 - child]))
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

//С помощью неполной НЕВОЗРАСТАЮЩЕЙ кучи 
//мелкие элементы закидываем поближе к началу массива 
void invreheap_i(Ak2& a, integer length, integer i)
{
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	Ak1 Temp;
	Temp.i = a.i[length - 1 - i];
	Temp.j = a.j[length - 1 - i];
	Temp.aij = a.aij[length - 1 - i];

	integer parent = i;
	integer child = 2 * (i + 1) - 1;
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
			if (a.i[length - 1 - child] <= a.i[length - 1 - (child + 1)])
			{
				child += 1;
			}
		}
		//Родитель больше потомков? 
		if (Temp.i > a.i[length - 1 - child])
		{
			//Тогда с этим родителем и его потомками разобрались 
			done = true;
		}
		else
		{
			//Родитель НЕ больше чем наибольший из его потомков. 
			//Перемещаем потомка на место родителя 
			//и с родителем в цикле сравниваем уже потомков этого потомка 
			integer itemp = length - 1 - parent;
			integer jtemp = length - 1 - child;
			a.i[itemp] = a.i[jtemp];
			a.j[itemp] = a.j[jtemp];
			a.aij[itemp] = a.aij[jtemp];
			parent = child;
			child = 2 * (parent + 1) - 1;
		}
	}
	//Родитель, с которого всё начиналось 
	//передвигается ближе к началу массива 
	//(или остаётся на месте если не повезло)
	integer itemp = length - 1 - parent;
	a.i[itemp] = Temp.i;
	a.j[itemp] = Temp.j;
	a.aij[itemp] = Temp.aij;
}

// С помощью неполной НЕВОЗРАСТАЮЩЕЙ кучи 
// мелкие элементы закидываем поближе к началу массива 
void invreheap_j(Ak2& a, integer length, integer i)
{
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	Ak1 Temp;
	Temp.i = a.i[length - 1 - i];
	Temp.j = a.j[length - 1 - i];
	Temp.aij = a.aij[length - 1 - i];

	integer parent = i;
	integer child = 2 * (i + 1) - 1;
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
			if (a.j[length - 1 - child] <= a.j[length - 1 - (child + 1)])
			{
				child += 1;
			}
		}
		//Родитель больше потомков? 
		if (Temp.j > a.j[length - 1 - child])
		{
			//Тогда с этим родителем и его потомками разобрались 
			done = true;
		}
		else
		{
			//Родитель НЕ больше чем наибольший из его потомков. 
			//Перемещаем потомка на место родителя 
			//и с родителем в цикле сравниваем уже потомков этого потомка 
			integer itemp = length - 1 - parent;
			integer jtemp = length - 1 - child;
			a.i[itemp] = a.i[jtemp];
			a.j[itemp] = a.j[jtemp];
			a.aij[itemp] = a.aij[jtemp];
			parent = child;
			child = 2 * (parent + 1) - 1;
		}
	}
	//Родитель, с которого всё начиналось 
	//передвигается ближе к началу массива 
	//(или остаётся на месте если не повезло)
	integer itemp = length - 1 - parent;
	a.i[itemp] = Temp.i;
	a.j[itemp] = Temp.j;
	a.aij[itemp] = Temp.aij;
}


/** * Демонстраицонный алгоритм для J-сортировки (JSort).
* Автор алгоритма - Джейсон Моррисон (Jason Morrison)
* <http://www.scs.carleton.ca/~morrison>
* * JSortAlgorithm.java
* * Автор реализации - Патрик Морин
* @author Patrick Morin */

//Основная процедура сортировки 
void sortJ_amg(Ak1*& a, integer first, integer last,
	integer (*indx_compare)(Ak1& Amat))
{
	// Строим неубывающую кучу 
	// Большие элементы из начала массива 
	// закидываем поближе к концу 
	for (integer i = last; i >= first; i--) reheap(a, last + 1, i, indx_compare);
	// Строим невозрастающую кучу
	// Меньшие элементы из конца массива 
	// закидываем поближе к началу 
	for (integer i = last; i >= first; i--) invreheap(a, last + 1, i, indx_compare);
	// Массив ПОЧТИ упорядочен 
	// Досортировываем вставками 
	for (integer j = first + 1; j <= last; j++)
	{
		Ak1 Temp = a[j];
		integer i = j - 1;
		while (i >= first && indx_compare(a[i]) > indx_compare(Temp))
		{
			a[i + 1] = a[i];
			i -= 1;
		}
		a[i + 1] = Temp;
	}
} // sortJ_amg

//Основная процедура сортировки 
void sortJ_amg(Ak2& a, integer first, integer last)
{
	//Строим неубывающую кучу 
	//Большие элементы из начала массива 
	//закидываем поближе к концу 
	for (integer i = last; i >= first; i--) reheap_i(a, last + 1, i);
	//Строим невозрастающую кучу
	//Меньшие элементы из конца массива 
	//закидываем поближе к началу 
	for (integer i = last; i >= first; i--) invreheap_i(a, last + 1, i);
	//Массив ПОЧТИ упорядочен 
	//Досортировываем вставками 
	for (integer j = first + 1; j <= last; j++)
	{
		Ak1 Temp;
		Temp.i = a.i[j];
		Temp.j = a.j[j];
		Temp.aij = a.aij[j];

		integer i = j - 1;
		while (i >= first && a.i[i] > Temp.i)
		{
			a.i[i + 1] = a.i[i];
			a.j[i + 1] = a.j[i];
			a.aij[i + 1] = a.aij[i];
			i -= 1;
		}
		a.i[i + 1] = Temp.i;
		a.j[i + 1] = Temp.j;
		a.aij[i + 1] = Temp.aij;
	}
} // sortJ_amg


void sortJ_amg_j(Ak2& a, integer first, integer last)
{
	//Строим неубывающую кучу 
	//Большие элементы из начала массива 
	//закидываем поближе к концу 
	for (integer i = last; i >= first; i--) reheap_j(a, last + 1, i);
	//Строим невозрастающую кучу
	//Меньшие элементы из конца массива 
	//закидываем поближе к началу 
	for (integer i = last; i >= first; i--) invreheap_j(a, last + 1, i);
	//Массив ПОЧТИ упорядочен 
	//Досортировываем вставками 
	for (integer j = first + 1; j <= last; j++)
	{
		Ak1 Temp;
		Temp.i = a.i[j];
		Temp.j = a.j[j];
		Temp.aij = a.aij[j];

		integer i = j - 1;
		while (i >= first && a.j[i] > Temp.j)
		{
			a.i[i + 1] = a.i[i];
			a.j[i + 1] = a.j[i];
			a.aij[i + 1] = a.aij[i];
			i -= 1;
		}
		a.i[i + 1] = Temp.i;
		a.j[i + 1] = Temp.j;
		a.aij[i + 1] = Temp.aij;
	}
} // sortJ_amg_j

//Очень коротко суть алгоритма TimSort можно объяснить так:
// 1. По специальному алгоритму разделяем входной массив на подмассивы.
// 2. Сортируем каждый подмассив обычной сортировкой вставками.
// 3. Собираем отсортированные подмассивы в единый массив с помощью модифицированной сортировки слиянием.
// Дьявол, как всегда, скрывается в деталях, а именно в алгоритме 
// из пункта 1 и модификации сортировки слиянием из пункта 3.

// Требует дополнительной памяти O(n).
// iterative Timsort function to sort the
// Amat[first...last] (similar to merge sort)
void timSort_amg(Ak2& Amat, integer first, integer last)
{
	// Sort individual subarrays of size RUN
	for (integer i = first; i <= last; i += RUN) {
		//insertionSortTim_amg(Amat, i, min(i + RUN - 1, last));
		sortJ_amg(Amat, i, (i + RUN - 1 < last ? i + RUN - 1 : last));
	}

	// start merging from size RUN (or 32). It will merge
	// to form size 64, then 128, 256 and so on ....
	for (integer size = RUN; size < last - first + 1; size = 2 * size)
	{
		// pick starting point of left sub array. We
		// are going to merge arr[left..left+size-1]
		// and arr[left+size, left+2*size-1]
		// After every merge, we increase left by 2*size
		for (integer left = first; left <= last; left += 2 * size)
		{
			// find ending point of left sub array
			// mid+1 is starting point of right sub array
			integer mid = left + size - 1;
			integer right = (left + 2 * size - 1 < last ? left + 2 * size - 1 : last);

			// слить подмассивы arr[left.....mid] &
			// arr[mid+1....right]
			if ((mid < right) && (mid >= left)) {
				mergeTim_amg(Amat, left, mid, right);
			}
		}
	}
} //  timSort_amg


// Требует дополнительной памяти O(n).
// iterative Timsort function to sort the
// Amat[first...last] (similar to merge sort)
void timSort_amg_j(Ak2& Amat, integer first, integer last)
{
	// Sort individual subarrays of size RUN
	for (integer i = first; i <= last; i += RUN) {
		//insertionSortTim_amg_j(Amat, i, min(i + RUN - 1, last));
		sortJ_amg_j(Amat, i, (i + RUN - 1 < last ? i + RUN - 1 : last));
	}

	// start merging from size RUN (or 32). It will merge
	// to form size 64, then 128, 256 and so on ....
	for (integer size = RUN; size < last - first + 1; size = 2 * size)
	{
		// pick starting point of left sub array. We
		// are going to merge arr[left..left+size-1]
		// and arr[left+size, left+2*size-1]
		// After every merge, we increase left by 2*size
		for (integer left = first; left <= last; left += 2 * size)
		{
			// find ending point of left sub array
			// mid+1 is starting point of right sub array
			integer mid = left + size - 1;
			integer right = (left + 2 * size - 1 < last ? left + 2 * size - 1 : last);

			// merge sub array arr[left.....mid] &
			// arr[mid+1....right]
			if ((mid < right) && (mid >= left)) {
				mergeTim_amg_j(Amat, left, mid, right);
			}
		}
	}
} //  timSort_amg_j


// Требует дополнительной памяти O(n).
// iterative Timsort function to sort the
// Amat[first...last] (similar to merge sort)
void timSort_amg(Ak1*& Amat, integer first, integer last, 
	integer (*indx_compare)(Ak1& Amat))
{
	// Sort individual subarrays of size RUN
	for (integer i = first; i <= last; i += RUN) {
		//insertionSortTim_amg(Amat, i, min((i + RUN-1), (last)));
		sortJ_amg(Amat, i, (i + RUN - 1 < last ? i + RUN - 1 : last), indx_compare);
	}

	// start merging from size RUN (or 32). It will merge
	// to form size 64, then 128, 256 and so on ....
	for (integer size = RUN; size < last - first + 1; size = 2 * size)
	{
		// pick starting point of left sub array. We
		// are going to merge arr[left..left+size-1]
		// and arr[left+size, left+2*size-1]
		// After every merge, we increase left by 2*size
		for (integer left = first; left <= last; left += 2 * size)
		{
			// find ending point of left sub array
			// mid+1 is starting point of right sub array
			integer mid = left + size - 1;
			integer right = (left + 2 * size - 1 < last ? left + 2 * size - 1 : last);

			// merge sub array arr[left.....mid] &
			// arr[mid+1....right]
			if ((mid < right) && (mid >= left)) {
				mergeTim_amg(Amat, left, mid, right, indx_compare);
			}
		}
	}
} //  timSort_amg


// Переформировать пирамиду
void FixHeap_j(Ak2& Amat,
	integer root,
	Ak1 m,
	integer bound,
	integer iadd)
{
	integer vacant;
	

	// list сортируемый список пирамида
	// root номер корня пирамиды
	// m ключевое значение вставляемое в пирамиду
	// bound правая граница (номер) в пирамиде
	vacant = root;
	while (2 * vacant <= bound)
	{
		integer largerChild;
		largerChild = 2 * vacant;
		integer lCadd = largerChild + iadd;
		integer lCadd1 = lCadd + 1;

		// поиск наибольшего из двух непосредственных потомков
		bool compare_result = false;
		if (Amat.j[lCadd1] > Amat.j[lCadd]) {
			compare_result = true;
		}
		else if (Amat.j[lCadd1] == Amat.j[lCadd]) {
			if (Amat.i[lCadd1] > Amat.i[lCadd]) {
				compare_result = true;
			}
		}
		if ((largerChild < bound) && (compare_result /*Amat.j[largerChild + 1+iadd]*n+ Amat.i[largerChild + 1 + iadd]> Amat.j[largerChild + iadd]*n+ Amat.i[largerChild + iadd]*/))
		{
			largerChild = largerChild + 1;
		}

		lCadd = largerChild + iadd;
		// находится ли ключ выше текущего потомка ?
		compare_result = false;
		if (m.j > Amat.j[lCadd]) {
			compare_result = true;
		}
		else if (m.j == Amat.j[lCadd]) {
			if (m.i > Amat.i[lCadd]) {
				compare_result = true;
			}
		}
		if (compare_result /*m.j*n+m.i >  Amat.j[largerChild + iadd]*n+ Amat.i[largerChild + iadd]*/)
		{
			// да, цикл завершается
			break;
		}
		else
		{
			// нет, большего непосредственного потомка
			// следует поднять
			Amat.i[vacant + iadd] = Amat.i[lCadd];
			Amat.j[vacant + iadd] = Amat.j[lCadd];
			Amat.aij[vacant + iadd] = Amat.aij[lCadd];
			vacant = largerChild;
		}
	}
	Amat.i[vacant + iadd] = m.i;
	Amat.j[vacant + iadd] = m.j;
	Amat.aij[vacant + iadd] = m.aij;
} // FixHeap_j

// Переформировать пирамиду
void FixHeap(Ak2& Amat,
	integer root,
	Ak1 m,
	integer bound,
	integer iadd)
{
	integer vacant;
	

	// list сортируемый список пирамида
	// root номер корня пирамиды
	// m ключевое значение вставляемое в пирамиду
	// bound правая граница (номер) в пирамиде
	vacant = root;
	while (2 * vacant <= bound)
	{
		integer largerChild;
		largerChild = 2 * vacant;
		integer lCadd = largerChild + iadd;
		integer lCadd1 = lCadd + 1;

		// поиск наибольшего из двух непосредственных потомков
		//integer key1 = Amat.i[largerChild + 1 + iadd]*n + Amat.j[largerChild + 1 + iadd];
		//integer key2 = Amat.i[largerChild + iadd]*n + Amat.j[largerChild + iadd];
		bool compare_result = false;
		if (Amat.i[lCadd1] > Amat.i[lCadd]) {
			compare_result = true;
		}
		else if (Amat.i[lCadd1] == Amat.i[lCadd]) {
			if (Amat.j[lCadd1] > Amat.j[lCadd]) {
				compare_result = true;
			}
		}
		if ((largerChild < bound) && compare_result/*(key1>key2)*/)
		{
			largerChild = largerChild + 1;
		}

		lCadd = largerChild + iadd;
		// находится ли ключ выше текущего потомка ?
		//integer key5 = m.i*n + m.j;
		//integer key6 = Amat.i[largerChild + iadd]*n + Amat.j[largerChild + iadd];
		compare_result = false;
		if (m.i > Amat.i[lCadd]) {
			compare_result = true;
		}
		else if (m.i == Amat.i[lCadd]) {
			if (m.j > Amat.j[lCadd]) {
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
			Amat.i[vacant + iadd] = Amat.i[lCadd];
			Amat.j[vacant + iadd] = Amat.j[lCadd];
			Amat.aij[vacant + iadd] = Amat.aij[lCadd];
			vacant = largerChild;
		}
	}
	Amat.i[vacant + iadd] = m.i;
	Amat.j[vacant + iadd] = m.j;
	Amat.aij[vacant + iadd] = m.aij;
} // FixHeap

// Пирамидальная сортировка оптимальна как
// по памяти, так и по быстродействию, к тому же её алгоритм
// очень интересен.
// Ограничение состоит в том, что нумерация массива должна начинаться с 1.
void HeapSort(Ak2& Amat, integer first, integer last)
{

	Ak1 maxelm; // элемент с наибольшим значением ключа

	// конструирование пирамиды
	for (integer i = ((last - first + 1) / 2); i >= 1; i--)
	{
		Ak1 temp;
		temp.i = Amat.i[i + first - 1];
		temp.j = Amat.j[i + first - 1];
		temp.aij = Amat.aij[i + first - 1];
		FixHeap(Amat, i, temp, last - first + 1, first - 1);
	}
	for (integer i = last - first + 1; i >= 2; i--)
	{
		// скопировать корень пирамиды в список
		// переформировать пирамиду
		maxelm.i = Amat.i[first];
		maxelm.j = Amat.j[first];
		maxelm.aij = Amat.aij[first];
		Ak1 temp;
		temp.i = Amat.i[i + first - 1];
		temp.j = Amat.j[i + first - 1];
		temp.aij = Amat.aij[i + first - 1];
		FixHeap(Amat, 1, temp, i - 1, first - 1);
		Amat.i[i + first - 1] = maxelm.i;
		Amat.j[i + first - 1] = maxelm.j;
		Amat.aij[i + first - 1] = maxelm.aij;
	}
} // HeapSort

// Пирамидальная сортировка оптимальна как
// по памяти, так и по быстродействию, к тому же её алгоритм
// очень интересен.
// Ограничение состоит в том, что нумерация массива должна начинаться с 1.
void HeapSort_j(Ak2& Amat, integer first, integer last)
{

	Ak1 maxelm; // элемент с наибольшим значением ключа

	integer iadd = first - 1;
	// конструирование пирамиды
	for (integer i = ((last - first + 1) / 2); i >= 1; i--)
	{
		Ak1 temp;
		temp.i = Amat.i[i + iadd];
		temp.j = Amat.j[i + iadd];
		temp.aij = Amat.aij[i + iadd];
		FixHeap_j(Amat, i, temp, last - first + 1, iadd);
	}
	for (integer i = last - first + 1; i >= 2; i--)
	{
		// скопировать корень пирамиды в список
		// переформировать пирамиду
		maxelm.i = Amat.i[1 + iadd];
		maxelm.j = Amat.j[1 + iadd];
		maxelm.aij = Amat.aij[1 + iadd];
		Ak1 temp;
		temp.i = Amat.i[i + iadd];
		temp.j = Amat.j[i + iadd];
		temp.aij = Amat.aij[i + iadd];
		FixHeap_j(Amat, 1, temp, i - 1, iadd);
		Amat.i[i + iadd] = maxelm.i;
		Amat.j[i + iadd] = maxelm.j;
		Amat.aij[i + iadd] = maxelm.aij;
	}
} // HeapSort_j

// Сортировка подсчётом.
// За время O(n)
// Томас Кормен стр. 224.
// Внимание: Алгоритм потребляет очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sort(Ak2& Amat, integer first, integer last, bool bmemo)
{
	// смена на malloc и calloc 7 января 2016.
	//если bmemo==true то запоминаем первоначальный порядок значений.
	integer* the_original_order_of_values_buf = nullptr;

	integer k = -1;
	for (integer j = first; j <= last; j++) {
		if (Amat.i[j] > k) k = Amat.i[j];
	}
	//integer* C = new integer[k + 1];
	integer* C = (integer*)malloc((k + 1) * sizeof(integer));
	char c1[2] = "C";
	char c2[14] = "Counting_Sort";
	handle_error(C, c1, c2, (k + 1));

	the_original_order_of_values_buf = (integer*)malloc((last + 1) * sizeof(integer));
	char c7[34] = "the_original_order_of_values_buf";
	char c6[14] = "Counting_Sort";
	handle_error(the_original_order_of_values_buf, c7, c6, (last + 1));

	if (bmemo) {
		the_original_order_of_values = (integer*)malloc((last + 1) * sizeof(integer));
		char c5[29] = "the_original_order_of_values";

		handle_error(the_original_order_of_values, c5, c6, (last + 1));

		the_original_order_of_values_reverse = (integer*)malloc((last + 1) * sizeof(integer));
		char c8[38] = "the_original_order_of_values_reverse";
		handle_error(the_original_order_of_values_reverse, c8, c6, (last + 1));

	}

#pragma omp parallel for
	for (integer i = 0; i <= k; i++) {
		C[i] = 0; // инициализация.
	}
	for (integer j = first; j <= last; j++) {
		C[Amat.i[j]]++;
	}
	// В C[i] хранится количество элементов равных i.
//НИ в коем случае !!! #pragma omp parallel for
	for (integer i = 1; i <= k; i++) {
		C[i] += C[i - 1];
	}
	// В C[i] количество элементов не превышающих i
	//Ak1* Bm = new Ak1[last - first + 2];
	Ak1* Bm = (Ak1*)malloc((last - first + 2) * sizeof(Ak1));
	char c3[3] = "Bm";
	char c4[14] = "Counting_Sort";
	handle_error(Bm, c3, c4, (last - first + 2));

	
	for (integer j = last; j >= first; j--) {
		integer ind;
		ind = Amat.i[j];
		Bm[C[ind]].i = Amat.i[j];
		Bm[C[ind]].j = Amat.j[j];
		Bm[C[ind]].aij = Amat.aij[j];
		if (bmemo) {
			// j стал C[ind]
			if (the_original_order_of_values_buf != nullptr) {
				the_original_order_of_values_buf[C[ind]] = j;
			}
		}
		C[ind]--;
	}
	// Обратное копирование.
	for (integer jnew = first, i = 1; jnew <= last; jnew++, i++) {
		//Amat.i[jnew] = B[jnew - first + 1].i;
		//Amat.j[jnew] = B[jnew - first + 1].j;
		//Amat.aij[jnew] = B[jnew - first + 1].aij;
		// i стал jnew. i ассоциируется с C[ind].
		Amat.i[jnew] = Bm[i].i;
		Amat.j[jnew] = Bm[i].j;
		Amat.aij[jnew] = Bm[i].aij;
		if (bmemo) {
			if (the_original_order_of_values != nullptr) {
				the_original_order_of_values[the_original_order_of_values_buf[i]] = jnew;
				the_original_order_of_values_reverse[jnew] = the_original_order_of_values_buf[i];
			}
		}
	}
	//delete[] Bm;
	if (Bm != nullptr) {
		free(Bm);
		Bm = nullptr;
	}
	//delete[] C;
	if (C != nullptr) {
		free(C);
		C = nullptr;
	}

	if (the_original_order_of_values_buf != nullptr) {
		free(the_original_order_of_values_buf);
		the_original_order_of_values_buf = nullptr;
	}


} // Counting_Sort

// Сортировка подсчётом.
// За время O(n)
// Томас Кормен стр. 224.
// Внимание: Алгоритм потребляет очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sort(Ak2& Amat, integer first, integer last, bool bmemo, integer bucket_len)
{
	// смена на malloc и calloc 7 января 2016.
	//если bmemo==true то запоминаем первоначальный порядок значений.
	integer* the_original_order_of_values_buf = nullptr;

	//integer bucket_len = -1;
	//for (integer j = first; j <= last; j++) {
		//if (Amat.i[j] > bucket_len) bucket_len = Amat.i[j];
	//}
	//integer* C = new integer[k + 1];
	integer* C = (integer*)malloc((bucket_len + 1) * sizeof(integer));
	char c1[2] = "C";
	char c2[14] = "Counting_Sort";
	handle_error<integer>(C, c1, c2, (bucket_len + 1));

	the_original_order_of_values_buf = (integer*)malloc((last + 1) * sizeof(integer));
	char c7[34] = "the_original_order_of_values_buf";
	char c6[14] = "Counting_Sort";
	handle_error<integer>(the_original_order_of_values_buf, c7, c6, (last + 1));

	if (bmemo) {
		the_original_order_of_values = (integer*)malloc((last + 1) * sizeof(integer));
		char c5[29] = "the_original_order_of_values";

		handle_error<integer>(the_original_order_of_values, c5, c6, (last + 1));

		the_original_order_of_values_reverse = (integer*)malloc((last + 1) * sizeof(integer));
		char c8[38] = "the_original_order_of_values_reverse";
		handle_error<integer>(the_original_order_of_values_reverse, c8, c6, (last + 1));

	}

#pragma omp parallel for
	for (integer i = 0; i <= bucket_len; i++) {
		C[i] = 0; // инициализация.
	}
	//memset(C,0,sizeof(integer)*(bucket_len + 1));


	for (integer j = first; j <= last; j++) {
		C[Amat.i[j]]++;
	}
	// В C[i] хранится количество элементов равных i.
//НИ в коем случае !!! #pragma omp parallel for
	for (integer i = 1; i <= bucket_len; i++) {
		C[i] += C[i - 1];
	}
	// В C[i] количество элементов не превышающих i
	//Ak1* Bm = new Ak1[last - first + 2];
	Ak1* Bm = (Ak1*)malloc((last - first + 2) * sizeof(Ak1));
	char c3[3] = "Bm";
	char c4[14] = "Counting_Sort";
	handle_error<Ak1>(Bm, c3, c4, (last - first + 2));

	integer ind;
	for (integer j = last; j >= first; j--) {
		ind = Amat.i[j];
		Bm[C[ind]].i = Amat.i[j];
		Bm[C[ind]].j = Amat.j[j];
		Bm[C[ind]].aij = Amat.aij[j];
		if (bmemo) {
			// j стал C[ind]
			if (the_original_order_of_values_buf != nullptr) {
				the_original_order_of_values_buf[C[ind]] = j;
			}
		}
		C[ind]--;
	}
	// Обратное копирование.
	for (integer jnew = first, i = 1; jnew <= last; jnew++, i++) {
		//Amat.i[jnew] = B[jnew - first + 1].i;
		//Amat.j[jnew] = B[jnew - first + 1].j;
		//Amat.aij[jnew] = B[jnew - first + 1].aij;
		// i стал jnew. i ассоциируется с C[ind].
		Amat.i[jnew] = Bm[i].i;
		Amat.j[jnew] = Bm[i].j;
		Amat.aij[jnew] = Bm[i].aij;
		if (bmemo) {
			if (the_original_order_of_values != nullptr) {
				the_original_order_of_values[the_original_order_of_values_buf[i]] = jnew;
				the_original_order_of_values_reverse[jnew] = the_original_order_of_values_buf[i];
			}
		}
	}
	//delete[] Bm;
	if (Bm != nullptr) {
		free(Bm);
		Bm = nullptr;
	}
	//delete[] C;
	if (C != nullptr) {
		free(C);
		C = nullptr;
	}

	if (the_original_order_of_values_buf != nullptr) {
		free(the_original_order_of_values_buf);
		the_original_order_of_values_buf = nullptr;
	}


} // Counting_Sort

// Сортировка подсчётом.
// За время O(n)
// Томас Кормен стр. 224.
// Внимание: Алгоритм потребляет очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sortj(Ak2& Amat, integer first, integer last)
{

	integer k = -1;
	for (integer j = first; j <= last; j++) {
		if (Amat.j[j] > k) k = Amat.j[j];
	}
	//integer* C = new integer[k + 1];
	integer* C = (integer*)malloc((k + 1) * sizeof(integer));
	char c1[2] = "C";
	char c2[14] = "Counting_Sort";
	handle_error(C, c1, c2, (k + 1));

#pragma omp parallel for
	for (integer i = 0; i <= k; i++) {
		C[i] = 0; // инициализация.
	}
	for (integer j = first; j <= last; j++) {
		C[Amat.j[j]]++;
	}
	// В C[i] хранится количество элементов равных i.
	for (integer i = 1; i <= k; i++) {
		C[i] += C[i - 1];
	}
	// В C[i] количество элементов не превышающих i
	//Ak1* Bm = new Ak1[last - first + 2];
	Ak1* Bm = (Ak1*)malloc((last - first + 2) * sizeof(Ak1));
	char c3[3] = "Bm";
	char c4[14] = "Counting_Sort";
	handle_error(Bm, c3, c4, (last - first + 2));

	
	for (integer j = last; j >= first; j--) {
		integer ind;
		ind = Amat.j[j];
		Bm[C[ind]].i = Amat.i[j];
		Bm[C[ind]].j = Amat.j[j];
		Bm[C[ind]].aij = Amat.aij[j];
		C[ind]--;
	}
	// Обратное копирование.
	for (integer j = first, i = 1; j <= last; j++, i++) {
		//Amat.i[j] = B[j - first + 1].i;
		//Amat.j[j] = B[j - first + 1].j;
		//Amat.aij[j] = B[j - first + 1].aij;
		Amat.i[j] = Bm[i].i;
		Amat.j[j] = Bm[i].j;
		Amat.aij[j] = Bm[i].aij;
	}
	//delete[] Bm;
	free(Bm);
	//delete[] C;
	free(C);

}// Counting_Sortj

// Сортировка подсчётом.
// За время O(n)
// Томас Кормен стр. 224.
// Внимание: Алгоритм потребляет очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sortj(Ak2& Amat, integer first, integer last, integer bucket_len)
{

	//integer bucket_len = -1;
	//for (integer j = first; j <= last; j++) {
		//if (Amat.j[j] > bucket_len) bucket_len = Amat.j[j];
	//}
	//integer* C = new integer[k + 1];
	integer* C = (integer*)malloc((bucket_len + 1) * sizeof(integer));
	char c1[2] = "C";
	char c2[14] = "Counting_Sort";
	handle_error<integer>(C, c1, c2, (bucket_len + 1));

#pragma omp parallel for
	for (integer i = 0; i <= bucket_len; i++) {
		C[i] = 0; // инициализация.
	}
	//memset(C, 0, sizeof(integer)*(bucket_len + 1));
	for (integer j = first; j <= last; j++) {
		C[Amat.j[j]]++;
	}
	// В C[i] хранится количество элементов равных i.
	for (integer i = 1; i <= bucket_len; i++) {
		C[i] += C[i - 1];
	}
	// В C[i] количество элементов не превышающих i
	//Ak1* Bm = new Ak1[last - first + 2];
	Ak1* Bm = (Ak1*)malloc((last - first + 2) * sizeof(Ak1));
	char c3[3] = "Bm";
	char c4[14] = "Counting_Sort";
	handle_error<Ak1>(Bm, c3, c4, (last - first + 2));

	
	for (integer j = last; j >= first; j--) {
		integer ind;
		ind = Amat.j[j];
		Bm[C[ind]].i = Amat.i[j];
		Bm[C[ind]].j = Amat.j[j];
		Bm[C[ind]].aij = Amat.aij[j];
		C[ind]--;
	}
	// Обратное копирование.
	for (integer j = first, i = 1; j <= last; j++, i++) {
		//Amat.i[j] = B[j - first + 1].i;
		//Amat.j[j] = B[j - first + 1].j;
		//Amat.aij[j] = B[j - first + 1].aij;
		Amat.i[j] = Bm[i].i;
		Amat.j[j] = Bm[i].j;
		Amat.aij[j] = Bm[i].aij;
	}
	//delete[] Bm;
	free(Bm);
	//delete[] C;
	free(C);

}// Counting_Sortj



// Сортировка подсчётом.
// За время O(n)
// Томас Кормен стр. 224.
// Внимание: Алгоритм потребляет очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
// Однопоточный вариант, т.к. распараллеливание будет сделано на уровень выше
// путём разбиения на два подмассива сортировки их по отдельности разными потоками,
// а на заключительной части используется алгоритм слияния.
// 08.06.2021 Не содержит выделений и уничтожений памяти. Оперативная память выделяется централизованно 
// один раз для всех многосеточных уровней внутри setup_phase.
void Counting_Sort_bmemo_false(Ak1*& Amat, integer first, integer last, integer bucket_len,
	integer(*indx_compare)(Ak1& Amat), integer*& C, Ak1* &Bm)
{
	// смена на malloc и calloc 7 января 2016.


	//integer bucket_len = -1;
	//for (integer j = first; j <= last; j++) {
		//if (indx_compare(Amat[j]) > bucket_len) bucket_len = indx_compare(Amat[j]);
	//}
	//integer* C = new integer[bucket_len + 1];
	/*
	// Вынесен во внешний код, чтобы избежать частых alloc() и free().
	integer* C = (integer*)malloc((bucket_len + 1) * sizeof(integer));
	char c1[2] = "C";
	char c2[14] = "Counting_Sort";
	handle_error<integer>(C, c1, c2, (bucket_len + 1));

	// Инициализация также выполняется параллельно во внешнем коде.
//#pragma omp parallel for
	for (integer i = 0; i <= bucket_len; i++) {
		C[i] = 0; // инициализация.
	}*/
	//memset(C, 0, sizeof(integer)*(bucket_len + 1));
	for (integer j = first; j <= last; j++) {
		C[indx_compare(Amat[j])]++;
	}
	// В C[i] хранится количество элементов равных i.
//НИ в коем случае !!! #pragma omp parallel for
	for (integer i = 1; i <= bucket_len; i++) {
		C[i] += C[i - 1];
	}
	// В C[i] количество элементов не превышающих i
	//Ak1* Bm = new Ak1[last - first + 2];
	/*Ak1* Bm = (Ak1*)malloc((last - first + 2) * sizeof(Ak1));
	char c3[3] = "Bm";
	char c4[14] = "Counting_Sort";
	handle_error<Ak1>(Bm, c3, c4, (last - first + 2));
	*/

	for (integer j = last; j >= first; j--) {
		integer ind;
		ind = indx_compare(Amat[j]);
		Bm[C[ind]] = Amat[j];
		C[ind]--;
	}

	//delete[] C;
	/*if (C != nullptr) {
		free(C);
		C = nullptr;
	}*/

	// Обратное копирование.
	for (integer jnew = first, i = 1; jnew <= last; jnew++, i++) {
		//Amat[jnew] = B[jnew - first + 1];
		// i стал jnew. i ассоциируется с C[ind].
		Amat[jnew] = Bm[i];
	}
	//delete[] Bm;
	/*if (Bm != nullptr) {
		free(Bm);
		Bm = nullptr;
	}*/

}


// Сортировка подсчётом.
// За время O(n)
// Томас Кормен стр. 224.
// Внимание: Алгоритм потребляет очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
// Однопоточный вариант, т.к. распараллеливание будет сделано на уровень выше
// путём разбиения на два подмассива сортировки их по отдельности разными потоками,
// а на заключительной части используется алгоритм слияния.
// 08.06.2021
void Counting_Sort_bmemo_false(Ak1*& Amat, integer first, integer last,  integer bucket_len,
	integer(*indx_compare)(Ak1& Amat), integer* &C)
{
	// смена на malloc и calloc 7 января 2016.
	

	//integer bucket_len = -1;
	//for (integer j = first; j <= last; j++) {
		//if (indx_compare(Amat[j]) > bucket_len) bucket_len = indx_compare(Amat[j]);
	//}
	//integer* C = new integer[bucket_len + 1];
	/*
	// Вынесен во внешний код, чтобы избежать частых alloc() и free().
	integer* C = (integer*)malloc((bucket_len + 1) * sizeof(integer));
	char c1[2] = "C";
	char c2[14] = "Counting_Sort";
	handle_error<integer>(C, c1, c2, (bucket_len + 1));	

	// Инициализация также выполняется параллельно во внешнем коде.
//#pragma omp parallel for
	for (integer i = 0; i <= bucket_len; i++) {
		C[i] = 0; // инициализация.
	}*/
	//memset(C, 0, sizeof(integer)*(bucket_len + 1));
	for (integer j = first; j <= last; j++) {
		C[indx_compare(Amat[j])]++;
	}
	// В C[i] хранится количество элементов равных i.
//НИ в коем случае !!! #pragma omp parallel for
	for (integer i = 1; i <= bucket_len; i++) {
		C[i] += C[i - 1];
	}
	// В C[i] количество элементов не превышающих i
	//Ak1* Bm = new Ak1[last - first + 2];
	Ak1* Bm = (Ak1*)malloc((last - first + 2) * sizeof(Ak1));
	char c3[3] = "Bm";
	char c4[14] = "Counting_Sort";
	handle_error<Ak1>(Bm, c3, c4, (last - first + 2));


	for (integer j = last; j >= first; j--) {
		integer ind;
		ind = indx_compare(Amat[j]);
		Bm[C[ind]] = Amat[j];
		C[ind]--;
	}

	//delete[] C;
	/*if (C != nullptr) {
		free(C);
		C = nullptr;
	}*/

	// Обратное копирование.
	for (integer jnew = first, i = 1; jnew <= last; jnew++, i++) {
		//Amat[jnew] = B[jnew - first + 1];
		// i стал jnew. i ассоциируется с C[ind].
		Amat[jnew] = Bm[i];
	}
	//delete[] Bm;
	if (Bm != nullptr) {
		free(Bm);
		Bm = nullptr;
	}
	
}

// Сортировка подсчётом.
// За время O(n)
// Томас Кормен стр. 224.
// Внимание: Алгоритм потребляет очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sort(Ak1*& Amat, integer first, integer last, bool bmemo, integer bucket_len,
	integer (*indx_compare)(Ak1& Amat))
{
	// смена на malloc и calloc 7 января 2016.
	//если bmemo==true то запоминаем первоначальный порядок значений.
	integer* the_original_order_of_values_buf = nullptr;

	//integer bucket_len = -1;
	//for (integer j = first; j <= last; j++) {
		//if (indx_compare(Amat[j]) > bucket_len) bucket_len = indx_compare(Amat[j]);
	//}
	//integer* C = new integer[bucket_len + 1];
	integer* C = (integer*)malloc((bucket_len + 1) * sizeof(integer));
	char c1[2] = "C";
	char c2[14] = "Counting_Sort";
	handle_error<integer>(C, c1, c2, (bucket_len + 1));

	the_original_order_of_values_buf = (integer*)malloc((last + 1) * sizeof(integer));
	char c7[34] = "the_original_order_of_values_buf";
	char c6[14] = "Counting_Sort";
	handle_error<integer>(the_original_order_of_values_buf, c7, c6, (last + 1));

	if (bmemo) {
		the_original_order_of_values = (integer*)malloc((last + 1) * sizeof(integer));
		char c5[29] = "the_original_order_of_values";

		handle_error<integer>(the_original_order_of_values, c5, c6, (last + 1));

		the_original_order_of_values_reverse = (integer*)malloc((last + 1) * sizeof(integer));
		char c8[38] = "the_original_order_of_values_reverse";
		handle_error<integer>(the_original_order_of_values_reverse, c8, c6, (last + 1));

	}

#pragma omp parallel for
	for (integer i = 0; i <= bucket_len; i++) {
		C[i] = 0; // инициализация.
	}
	//memset(C, 0, sizeof(integer)*(bucket_len + 1));
	for (integer j = first; j <= last; j++) {
		C[indx_compare(Amat[j])]++;
	}
	// В C[i] хранится количество элементов равных i.
//НИ в коем случае !!! #pragma omp parallel for
	for (integer i = 1; i <= bucket_len; i++) {
		C[i] += C[i - 1];
	}
	// В C[i] количество элементов не превышающих i
	//Ak1* Bm = new Ak1[last - first + 2];
	Ak1* Bm = (Ak1*)malloc((last - first + 2) * sizeof(Ak1));
	char c3[3] = "Bm";
	char c4[14] = "Counting_Sort";
	handle_error<Ak1>(Bm, c3, c4, (last - first + 2));

	
	for (integer j = last; j >= first; j--) {
		integer ind;
		ind = indx_compare(Amat[j]);
		Bm[C[ind]] = Amat[j];
		if (bmemo) {
			// j стал C[ind]
			if (the_original_order_of_values_buf != nullptr) {
				the_original_order_of_values_buf[C[ind]] = j;
			}
		}
		C[ind]--;
	}
	// Обратное копирование.
	for (integer jnew = first, i = 1; jnew <= last; jnew++, i++) {
		//Amat[jnew] = B[jnew - first + 1];
		// i стал jnew. i ассоциируется с C[ind].
		Amat[jnew] = Bm[i];
		if (bmemo) {
			if (the_original_order_of_values != nullptr) {
				the_original_order_of_values[the_original_order_of_values_buf[i]] = jnew;
				the_original_order_of_values_reverse[jnew] = the_original_order_of_values_buf[i];
			}
		}
	}
	//delete[] Bm;
	if (Bm != nullptr) {
		free(Bm);
		Bm = nullptr;
	}
	//delete[] C;
	if (C != nullptr) {
		free(C);
		C = nullptr;
	}

	if (the_original_order_of_values_buf != nullptr) {
		free(the_original_order_of_values_buf);
		the_original_order_of_values_buf = nullptr;
	}
}

#include "leftist_heap.cpp" // Левосторонняя куча(пирамида). Алгоритм сортировки с её помощью.

integer compAi(integer a, integer b) {
	if (a > b) return (1);
	if (a < b) return (-1);
	return (0);
}// compAi

integer compAj(integer a, integer b) {
	if (a > b) return (1);
	if (a < b) return (-1);
	return (0);
}// compAj

// Версия быстрой сортировки Чарльза Хоара упорядочивание по строкам i.
// Скорость возрастает: heapsort -> quicksort -> timsort.
void qs(Ak2& Amat, integer first, integer last) {
	integer i = first, j = last;
	Ak1 tmp;

	const integer INSERTION_SIZE = 32;// 64;//my_amg_manager.memory_size_Stress;
	if (last - first < INSERTION_SIZE) {
		//insertionSortTim_amg(Amat, first, last);
		sortJ_amg(Amat, first, last);
	}
	else {
		/*
		В случае явной рекурсии, как в программе выше, в стеке сохраняются не только границы подмассивов, но и ряд совершенно ненужных параметров, таких как локальные переменные. Если эмулировать стек программно, его размер можно уменьшить в несколько раз.
		Чем на более равные части будет делиться массив - тем лучше. Потому в качестве опорного целесообразно брать средний из трех, а если массив достаточно велик - то из девяти произвольных элементов.
		Пусть входные последовательности очень плохи для алгоритма. Например, их специально подбирают, чтобы средний элемент оказывался каждый раз минимумом. Как сделать QuickSort устойчивой к такому "саботажу" ? Очень просто - выбирать в качестве опорного случайный элемент входного массива. Тогда любые неприятные закономерности во входном потоке будут нейтрализованы. Другой вариант - переставить перед сортировкой элементы массива случайным образом.
		Быструю сортировку можно использовать и для двусвязных списков. Единственная проблема при этом - отсутствие непосредственного доступа к случайному элементу. Так что в качестве опорного приходится выбирать первый элемент, и либо надеяться на хорошие исходные данные, либо случайным образом переставить элементы перед сортировкой.
		*/
		integer pivot;
		/*
		if (1) {
			if (last - first < 3000) {
				pivot = Amat.i[(integer)((first + last) / 2)];
			}
			else if (last - first < 800000) {
				pivot = (integer)((Amat.i[first + 100] + Amat.i[(first + last) / 2] + Amat.i[last - 100]) / 3.0);
			}
			else {
				pivot = 0;
				pivot = (integer)((Amat.i[first] + Amat.i[first + 100000] + Amat.i[first + 200000] + Amat.i[first + 300000] + Amat.i[first + 400000] + Amat.i[first + 500000] + Amat.i[first + 600000] + Amat.i[first + 700000] + Amat.i[last]) / 9.0);
			}

		}
		else {
		*/
		pivot = Amat.i[(first + last) / 2];
		//}

		// partition
		while (i <= j) {
			//while ((Amat.i[i] < pivot) || ((Amat.i[i] == pivot) && (Amat.j[i] < pivot.j)))
			//i++;
			//while ((Amat.i[j] > pivot) || ((Amat.i[j] == pivot) && (Amat.j[j] > pivot.j)))
			//j--;
			while (Amat.i[i] < pivot)
				i++;
			while (Amat.i[j] > pivot)
				j--;
			if (i <= j) {
				tmp.i = Amat.i[i];
				tmp.j = Amat.j[i];
				tmp.aij = Amat.aij[i];
				Amat.i[i] = Amat.i[j];
				Amat.j[i] = Amat.j[j];
				Amat.aij[i] = Amat.aij[j];
				Amat.i[j] = tmp.i;
				Amat.j[j] = tmp.j;
				Amat.aij[j] = tmp.aij;
				i++;
				j--;
			}
		}

		// recursion
		/*
		if (1) {
			if (first < j) {
				if (j - first <= INSERTION_SIZE)
				{
					insertionSortTim_amg(Amat, first, j);
				}
				else {
					qs(Amat, first, j);
				}
			}
			if (i < last) {
				if (last - i <= INSERTION_SIZE)
				{
					insertionSortTim_amg(Amat, i, last);
				}
				else {
					qs(Amat, i, last);
				}
			}
		}
		else {
			*/
		if (first < j)
			qs(Amat, first, j);
		if (i < last)
			qs(Amat, i, last);
		//}
	}

} //qs

// Версия быстрой сортировки Чарльза Хоара упорядочивание по строкам j.
// Скорость возрастает: heapsort -> quicksort -> timsort.
void qsj(Ak2& Amat, integer first, integer last) {

	const integer INSERTION_SIZE = 32;// 64;// my_amg_manager.memory_size_Stress;

	if (last - first < INSERTION_SIZE) {
		//insertionSortTim_amg_j(Amat, first, last);
		sortJ_amg_j(Amat, first, last);
	}
	else {

		integer i = first, j = last;
		Ak1 tmp;
		//Ak1 pivot = Amat[(first + last) / 2];
		/*
		В случае явной рекурсии, как в программе выше, в стеке сохраняются не только границы подмассивов, но и ряд совершенно ненужных параметров, таких как локальные переменные. Если эмулировать стек программно, его размер можно уменьшить в несколько раз.
		Чем на более равные части будет делиться массив - тем лучше. Потому в качестве опорного целесообразно брать средний из трех, а если массив достаточно велик - то из девяти произвольных элементов.
		Пусть входные последовательности очень плохи для алгоритма. Например, их специально подбирают, чтобы средний элемент оказывался каждый раз минимумом. Как сделать QuickSort устойчивой к такому "саботажу" ? Очень просто - выбирать в качестве опорного случайный элемент входного массива. Тогда любые неприятные закономерности во входном потоке будут нейтрализованы. Другой вариант - переставить перед сортировкой элементы массива случайным образом.
		Быструю сортировку можно использовать и для двусвязных списков. Единственная проблема при этом - отсутствие непосредственного доступа к случайному элементу. Так что в качестве опорного приходится выбирать первый элемент, и либо надеяться на хорошие исходные данные, либо случайным образом переставить элементы перед сортировкой.
		*/
		integer pivot;
		/*
		if (1) {
			if (last - first < 3000) {
				pivot = Amat.j[(integer)((first + last) / 2)];
			}
			else if (last - first < 800000) {
				pivot = (integer)((Amat.j[first + 100] + Amat.j[(first + last) / 2] + Amat.j[last - 100]) / 3.0);
			}
			else {
				pivot = 0;
				pivot = (integer)((Amat.j[first] + Amat.j[first + 100000] + Amat.j[first + 200000] + Amat.j[first + 300000] + Amat.j[first + 400000] + Amat.j[first + 500000] + Amat.j[first + 600000] + Amat.j[first + 700000] + Amat.j[last]) / 9.0);
			}

		}
		else {
		*/
		pivot = Amat.j[(first + last) / 2];
		//}



		// partition
		while (i <= j) {
			//while ((Amat.j[i] < pivot.j) || ((Amat.j[i] == pivot.j) && (Amat.i[i] < pivot.i)))
			//i++;
			//while ((Amat.j[j] > pivot.j) || ((Amat.j[j] == pivot.j) && (Amat.i[j] > pivot.i)))
			//j--;

			while (Amat.j[i] < pivot)
				i++;
			while (Amat.j[j] > pivot)
				j--;
			if (i <= j) {
				tmp.i = Amat.i[i];
				tmp.j = Amat.j[i];
				tmp.aij = Amat.aij[i];
				Amat.i[i] = Amat.i[j];
				Amat.j[i] = Amat.j[j];
				Amat.aij[i] = Amat.aij[j];
				Amat.i[j] = tmp.i;
				Amat.j[j] = tmp.j;
				Amat.aij[j] = tmp.aij;
				i++;
				j--;
			}
		}

		// recursion
		/*
		if (0) {
			if (first < j) {
				if (j - first <= INSERTION_SIZE)
				{
					insertionSortTim_amg_j(Amat, first, j);
				}
				else {
					qsj(Amat, first, j);
				}
			}
			if (i < last) {
				if (last - i <= INSERTION_SIZE)
				{
					insertionSortTim_amg_j(Amat, i, last);
				}
				else {
					qsj(Amat, i, last);
				}
			}
		}
		else {
		*/
		if (first < j)
			qsj(Amat, first, j);
		if (i < last)
			qsj(Amat, i, last);
		//}
	}

}

#endif
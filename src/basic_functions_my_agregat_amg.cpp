
#pragma once
#ifndef BASIC_FUNCTIONS_MY_AGREGAT_AMG_CPP
#define BASIC_FUNCTIONS_MY_AGREGAT_AMG_CPP 1


void swap(Ak2 &Amat, integer i, integer j)
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


// PivotList содержит ошибку.Обнаружено 12 декабря 2015.
// Вот алгоритм PivotList
integer PivotList(Ak2 &Amat, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue_j = Amat.j[first];
	integer PivotValue_i = Amat.i[first];
	integer PivotPoint = first;

	for (integer index = (first + 1); index <= last; index++) {
		if (Amat.i[index] < PivotValue_i) {
			PivotPoint++;
			swap(Amat, PivotPoint, index);
		}
		else if ((Amat.i[index] == PivotValue_i) && (Amat.j[index] < PivotValue_j)) {
			PivotPoint++;
			swap(Amat, PivotPoint, index);
		}
	}

	swap(Amat, first, PivotPoint);

	return PivotPoint;
} // PivotList

// Вот алгоритм PivotList
integer PivotList_j(Ak2 &Amat, integer first, integer last) {
	// list==jptr and altr обрабатываемый список
	// first номер первого элемента
	// last номер последнего элемента

	integer PivotValue_j = Amat.j[first];
	integer PivotValue_i = Amat.i[first];
	integer PivotPoint = first;

	for (integer index = (first + 1); index <= last; index++) {
		if (Amat.j[index] < PivotValue_j) {
			PivotPoint++;
			swap(Amat, PivotPoint, index);
		}
		else if ((Amat.j[index] == PivotValue_j) && (Amat.i[index] < PivotValue_i)) {
			PivotPoint++;
			swap(Amat, PivotPoint, index);
		}
	}

	swap(Amat, first, PivotPoint);

	return PivotPoint;
} // PivotListamg_j


// Сортировка Тима (Тим Петерсон).
// см. файл uniformsimplemeshgen.cpp
//const integer RUN = 32;

// this function sorts array from left index to
// to right index which is of size atmost RUN
void insertionSortTim_amg(Ak2 &Amat, integer left, integer right)
{
	for (integer i = left + 1; i <= right; i++)
	{
		Ak1 temp;
		temp.i= Amat.i[i];
		temp.j = Amat.j[i];
		temp.aij = Amat.aij[i];

		integer j = i - 1;
		while ((j >= left)&&(Amat.i[j] > temp.i))
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

// this function sorts array from left index to
// to right index which is of size atmost RUN
void insertionSortTim_amg_j(Ak2 &Amat, integer left, integer right)
{
	for (integer i = left + 1; i <= right; i++)
	{
		Ak1 temp;
		temp.i = Amat.i[i];
		temp.j = Amat.j[i];
		temp.aij = Amat.aij[i];

		integer j = i - 1;
		while ((j >= left)&&(Amat.j[j] > temp.j))
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

// this function sorts array from left index to
// to right index which is of size atmost RUN
void insertionSortTim_amg(Ak1* &Amat, integer left, integer right)
{
	for (integer i = left + 1; i <= right; i++)
	{
		Ak1 temp;
		temp = Amat[i];
		

		integer j = i - 1;
		while ((j >= left)&&(Amat[j].i> temp.i))
		{
			Amat[j + 1] = Amat[j];
			j--;
		}
		Amat[j + 1] = temp;
	}
}// insertion sort по i.

// this function sorts array from left index to
// to right index which is of size atmost RUN
void insertionSortTim_amg_j(Ak1* &Amat, integer left, integer right)
{
	for (integer i = left + 1; i <= right; i++)
	{
		Ak1 temp;
		temp = Amat[i];
		

		integer j = i - 1;
		while ((j >= left)&&(Amat[j].j > temp.j))
		{
			Amat[j + 1] = Amat[j];
			j--;
		}
		Amat[j + 1] = temp;
	}
} // insertion sort по j.

// merge function merges the sorted runs
void mergeTim_amg(Ak2 &Amat, integer l, integer m, integer r)
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
		left.i = new integer[len1 + 1];
		left.j = new integer[len1 + 1];
		left.aij = new doublereal[len1 + 1];
	}
	Ak2 right;
	right.i = nullptr;
	right.j = nullptr;
	right.aij = nullptr;
	if (len2 >= 0) {
		right.i = new integer[len2 + 1];
		right.j = new integer[len2 + 1];
		right.aij = new doublereal[len2 + 1];
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
void mergeTim_amg_j(Ak2 &Amat, integer l, integer m, integer r)
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
		left.i = new integer[len1 + 1];
		left.j = new integer[len1 + 1];
		left.aij = new doublereal[len1 + 1];
	}
	Ak2 right;
	right.i = nullptr;
	right.j = nullptr;
	right.aij = nullptr;
	if (len2 >= 0) {
		right.i = new integer[len2 + 1];
		right.j = new integer[len2 + 1];
		right.aij = new doublereal[len2 + 1];
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
void mergeTim_amg(Ak1* &Amat, integer l, integer m, integer r)
{
	// original array is broken in two parts
	// left and right array
	integer len1 = m - l + 1, len2 = r - m;

	//myARRT left[len1], right[len2];
	Ak1* left=nullptr;
	if (len1 >= 0) {
		left = new Ak1[len1 + 1];
	}
	Ak1* right=nullptr;
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
		if (left[i].i <= right[j].i)
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

// merge function merges the sorted runs
void mergeTim_amg_j(Ak1* &Amat, integer l, integer m, integer r)
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
		if (left[i].j <= right[j].j)
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

//С помошью неполной НЕУБЫВАЮЩЕЙ кучи 
//крупные элементы закидываем поближе к концу массива 
void reheap_i(Ak1* &a, integer length, integer i) {
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
			if (a[child].i >= a[child + 1].i) { child += 1; }
		}
		//Родитель меньше потомков? 
		if (Temp.i < a[child].i) {
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

//С помошью неполной НЕУБЫВАЮЩЕЙ кучи 
//крупные элементы закидываем поближе к концу массива 
void reheap_i(Ak2 &a, integer length, integer i) {
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	Ak1 Temp;
	Temp.i= a.i[i];
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

//С помошью неполной НЕУБЫВАЮЩЕЙ кучи 
//крупные элементы закидываем поближе к концу массива 
void reheap_j(Ak1* &a, integer length, integer i) {
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
			if (a[child].j >= a[child + 1].j) { child += 1; }
		}
		//Родитель меньше потомков? 
		if (Temp.j < a[child].j) {
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

//С помошью неполной НЕУБЫВАЮЩЕЙ кучи 
//крупные элементы закидываем поближе к концу массива 
void reheap_j(Ak2 &a, integer length, integer i) {
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

//С помошью неполной НЕВОЗРАСТАЮЩЕЙ кучи 
//мелкие элементы закидываем поближе к началу массива 
void invreheap_i(Ak1* &a, integer length, integer i)
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
			if (a[length - 1 - child].i <= a[length - 1 - (child + 1)].i)
			{
				child += 1;
			}
		}
		//Родитель больше потомков? 
		if (Temp.i > a[length - 1 - child].i)
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

//С помошью неполной НЕВОЗРАСТАЮЩЕЙ кучи 
//мелкие элементы закидываем поближе к началу массива 
void invreheap_i(Ak2 &a, integer length, integer i)
{
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	Ak1 Temp;
	Temp.i= a.i[length - 1 - i];
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


//С помошью неполной НЕВОЗРАСТАЮЩЕЙ кучи 
//мелкие элементы закидываем поближе к началу массива 
void invreheap_j(Ak1* &a, integer length, integer i)
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
			if (a[length - 1 - child].j <= a[length - 1 - (child + 1)].j)
			{
				child += 1;
			}
		}
		//Родитель больше потомков? 
		if (Temp.j > a[length - 1 - child].j)
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

// С помошью неполной НЕВОЗРАСТАЮЩЕЙ кучи 
// мелкие элементы закидываем поближе к началу массива 
void invreheap_j(Ak2 &a, integer length, integer i)
{
	//С этим родителем ещё не разобрались 
	bool done = false;
	//Запоминаем отдельно родителя 
	//и смотрим на его потомка слева 
	Ak1 Temp;
	Temp.i= a.i[length - 1 - i];
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
void sortJ_amg(Ak1* &a, integer first, integer last)
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
	for (integer j = first+1; j <= last; j++)
	{
		Ak1 Temp = a[j];
		integer i = j - 1;
		while (i >= first && a[i].i > Temp.i)
		{
			a[i + 1] = a[i];
			i -= 1;
		}
		a[i + 1] = Temp;
	}
} // sortJ_amg


void sortJ_amg_j(Ak1* &a, integer first, integer last)
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
		Ak1 Temp = a[j];
		integer i = j - 1;
		while (i >= first && a[i].j > Temp.j)
		{
			a[i + 1] = a[i];
			i -= 1;
		}
		a[i + 1] = Temp;
	}
} // sortJ_amg_j

//Основная процедура сортировки 
void sortJ_amg(Ak2 &a, integer first, integer last)
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
		Temp.i= a.i[j];
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


void sortJ_amg_j(Ak2 &a, integer first, integer last)
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

//Очень коротко суть алгоритма TimSort можно объяснить так :
// 1. По специальному алгоритму разделяем входной массив на подмассивы.
// 2. Сортируем каждый подмассив обычной сортировкой вставками.
// 3. Собираем отсортированные подмассивы в единый массив с помощью модифицированной сортировки слиянием.
// Дьявол, как всегда, скрывается в деталях, а именно в алгоритме 
// из пункта 1 и модификации сортировки слиянием из пункта 3.

// Требует дополнительной памяти O(n).
// iterative Timsort function to sort the
// Amat[first...last] (similar to merge sort)
void timSort_amg(Ak2 &Amat, integer first, integer last)
{
	// Sort individual subarrays of size RUN
	for (integer i = first; i <= last; i += RUN)
		insertionSortTim_amg(Amat, i, min(i + RUN - 1, last));

	// start merging from size RUN (or 32). It will merge
	// to form size 64, then 128, 256 and so on ....
	for (integer size = RUN; size < last-first+1; size = 2 * size)
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
			integer right = min(left + 2 * size - 1, last);

			// merge sub array arr[left.....mid] &
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
void timSort_amg_j(Ak2 &Amat, integer first, integer last)
{
	// Sort individual subarrays of size RUN
	for (integer i = first; i <= last; i += RUN)
		insertionSortTim_amg_j(Amat, i, min(i + RUN - 1, last));

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
			integer right = min(left + 2 * size - 1, last);

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
void timSort_amg(Ak1* &Amat, integer first, integer last)
{
	// Sort individual subarrays of size RUN
	for (integer i = first; i <= last; i += RUN) {
	//	insertionSortTim_amg(Amat, i, min((i + RUN-1), (last)));
		sortJ_amg(Amat, i, min(i + RUN - 1, last));
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
			integer right = min(left + 2 * size - 1, last);

			// merge sub array arr[left.....mid] &
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
void timSort_amg_j(Ak1* &Amat, integer first, integer last)
{
	// Sort individual subarrays of size RUN
	for (integer i = first; i <= last; i += RUN) {
		//insertionSortTim_amg_j(Amat, i, min((i + RUN-1), (last)));
		sortJ_amg_j(Amat, i, min(i + RUN - 1, last));
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
			integer right = min(left + 2 * size - 1, last);

			// merge sub array arr[left.....mid] &
			// arr[mid+1....right]
			if ((mid < right) && (mid >= left)) {
				mergeTim_amg_j(Amat, left, mid, right);
			}
		}
	}
} //  timSort_amg_j

// Переформировать пирамиду
void FixHeap_j(Ak2 &Amat,
	integer root,
	Ak1 m,
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
void FixHeap(Ak2 &Amat,
	integer root,
	Ak1 m,
	integer bound,
	integer iadd)
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
void HeapSort(Ak2 &Amat, integer first, integer last)
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
		maxelm.j= Amat.j[first];
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
void HeapSort_j(Ak2 &Amat, integer first, integer last)
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
// Внимание : Алгоритм жрёт очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sort(Ak2 &Amat, integer first, integer last, bool bmemo)
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
// Внимание : Алгоритм жрёт очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sort(Ak2 &Amat, integer first, integer last, bool bmemo, integer bucket_len)
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
// Внимание : Алгоритм жрёт очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sortj(Ak2 &Amat, integer first, integer last)
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

	integer ind;
	for (integer j = last; j >= first; j--) {
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
// Внимание : Алгоритм жрёт очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sortj(Ak2 &Amat, integer first, integer last, integer bucket_len)
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

	integer ind;
	for (integer j = last; j >= first; j--) {
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
// Внимание : Алгоритм жрёт очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sort(Ak1* &Amat, integer first, integer last, bool bmemo, integer bucket_len)
{
	// смена на malloc и calloc 7 января 2016.
	//если bmemo==true то запоминаем первоначальный порядок значений.
	integer* the_original_order_of_values_buf = nullptr;

	//integer bucket_len = -1;
	//for (integer j = first; j <= last; j++) {
		//if (Amat[j].i > bucket_len) bucket_len = Amat[j].i;
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
		C[Amat[j].i]++;
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
		ind = Amat[j].i;
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

// Сортировка подсчётом.
// За время O(n)
// Томас Кормен стр. 224.
// Внимание : Алгоритм жрёт очень много оперативной памяти.
// Впервые предложена Севардом (H.H.Seward) в 1954 году.
void Counting_Sortj(Ak1* &Amat, integer first, integer last, integer bucket_len)
{

	//integer bucket_len = -1;
	//for (integer j = first; j <= last; j++) {
		//if (Amat[j].j > bucket_len) bucket_len = Amat[j].j;
	//}
	//integer* C = new integer[bucket_len + 1];
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
		C[Amat[j].j]++;
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

	integer ind;
	for (integer j = last; j >= first; j--) {
		ind = Amat[j].j;
		Bm[C[ind]] = Amat[j];
		C[ind]--;
	}
	// Обратное копирование.
	for (integer j = first, i = 1; j <= last; j++, i++) {
		//Amat[j] = B[j - first + 1];
		Amat[j] = Bm[i];
	}
	//delete[] Bm;
	free(Bm);
	//delete[] C;
	free(C);

}


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

// Правильная версия сортировки Чарльза Хоара которая раз в 5 быстрее чем,
// пирамидальная сортировка. Но ещё быстрее обещает быть TimSort (Futures).
void qs(Ak2 &Amat, integer first, integer last) {
	integer i = first, j = last;
	Ak1 tmp;

	const integer INSERTION_SIZE = 64;//my_amg_manager.memory_size_Stress;
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

// Правильная версия сортировки Чарльза Хоара которая раз в 5 быстрее чем,
// пирамидальная сортировка. Но ещё быстрее обещает быть TimSort (Futures).
void qsj(Ak2 &Amat, integer first, integer last) {

	const integer INSERTION_SIZE = 64;// my_amg_manager.memory_size_Stress;

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




// 15 11 2016 метод полностью на malloc с контролем через handle_error.
// Мы усилим РУМБА0.14 алгоритм ILU2 предобуславливанием на каждом уровне вложенности.
// Преобразует equation3D  формат хранения в CRS формат.
// Цель написания этого преобразователя: экономия оперативной памяти компьютера.
// Т.к. формат SIMPLESPARSE требует слишком много памяти.
void equation3DtoCRSRUMBA1(LEVEL_ADDITIONAL_DATA &milu2,
	bool ballocmemory, Ak2 &Amat, integer istartq, integer iendq,
	integer* &row_ptr_start, integer* &row_ptr_end, integer iadd, integer ilevel) {

	// Если ballocmemory равен true то происходит выделение памяти.
	const bool bonly_negative_connections = false;
	bool flag = true;
	integer n = 0; // число ненулевых элементов
	integer maxelm_plus_maxbound = 0; // Число уравнений в СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;


	// подсчёт количества ненулевых элементов
	// во всех точках расчётной области как внутренних так и граничных.
	for (integer ii = startpos; ii <= endpos; ii++) {
		maxelm_plus_maxbound++;
		integer istr = ii - iadd;
		if (fabs(Amat.aij[row_ptr_start[ii]]) < 1.0e-20)
		{
#if doubleintprecision == 1
			printf("zero diagonal coefficient in level = %lld, istr=%lld\n", ilevel, istr);
#else
			printf("zero diagonal coefficient in level = %d, istr=%d\n", ilevel, istr);
#endif

			//getchar();
			system("PAUSE");
			exit(1);
		}
		doublereal ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

		//if (fabs(ap_now)> 1e10*nonzeroEPS) n++; // диагональный элемент
		if (ap_now > nonzeroEPS) n++; // Диагональный элемент.
		else {
			// 5 августа 2016. 
			flag = false;
			printf("internal zero diagonal element.\n");
			//printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e sum_nb=%e", sl[k].ae, sl[k].aw, sl[k].an, sl[k].as, sl[k].at, sl[k].ab, sl[k].ae + sl[k].aw + sl[k].an + sl[k].as + sl[k].at + sl[k].ab);
			if (ap_now < 0.0) {
				printf("found negativ diagonal coefficient=%e...\n", ap_now);
			}
			printf("fatal error equation3DtoCRS in RUMBA...\n");
			//getchar();
			system("PAUSE");
			exit(1);
			//n++;
			//sl[k].ap = fabs(sl[k].ae) + fabs(sl[k].aw) + fabs(sl[k].an) + fabs(sl[k].as) + fabs(sl[k].at) + fabs(sl[k].ab);
		}

		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];

		for (integer ii1 = is1; ii1 <= is2; ii1++)
		{
			if (bonly_negative_connections) {
				if (Amat.aij[ii1] < -nonzeroEPS) {
					n++;
				}
			}
			else {
				if (fabs(Amat.aij[ii1]) > nonzeroEPS) n++;
			}
		}
	}


	if (flag) {
		// memory +15N
		// Теперь выделение памяти будет происходить централизованно, вне данного кода.
		// Это сделано для кода BICGSTAB_internal3. дата изменения 12 апреля 2013.
		// Другой код, использующий equation3dtoCRS может оказаться неработоспособным после этого изменения.
		if (ballocmemory) {
			// Важно выделить память с запасом, т.к. одна и таже память используется и для компонент скорости и для попрапвки давления.
			//val = new doublereal[7 * (maxelm + maxbound) + 2 * maxbound + 2];
			//col_ind = new integer[7 * (maxelm + maxbound) + 2 * maxbound + 2];

			//milu2.val = new doublereal[n + 2];
			milu2.val = (doublereal*)malloc((n + 2) * sizeof(doublereal));
			char c1[11] = "milu2.val";
			char c2[23] = "equation3DtoCRSRUMBA1";
			handle_error(milu2.val, c1, c2, (n + 2));
			//milu2.col_ind = new integer[n + 2];
			milu2.col_ind = (integer*)malloc((n + 2) * sizeof(integer));
			char c3[14] = "milu2.col_ind";
			handle_error(milu2.col_ind, c3, c2, (n + 2));
			//row_ptr = new integer[(maxelm + maxbound) + 1];
			//milu2.row_ptr = new integer[maxelm_plus_maxbound + 1];
			milu2.row_ptr = (integer*)malloc((maxelm_plus_maxbound + 1) * sizeof(integer));
			char c4[14] = "milu2.row_ptr";
			handle_error(milu2.row_ptr, c4, c2, (maxelm_plus_maxbound + 1));

			if ((milu2.val == nullptr) || (milu2.col_ind == nullptr) || (milu2.row_ptr == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}


		// инициализация
#pragma omp parallel for
		for (integer k = 0; k<(n); k++) {
			milu2.val[k] = 0.0;
			milu2.col_ind[k] = -1;
		}
#pragma omp parallel for
		for (integer k = 0; k <= (maxelm_plus_maxbound); k++) {
			milu2.row_ptr[k] = n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}

		// Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (integer k=0; k<M.n; k++) {
		val[k]=M.a[k].aij;
		col_ind[k]=M.a[k].j;
		row_ptr[M.a[k].i]=min(k,row_ptr[M.a[k].i]);
		}
		*/
		integer ik = 0; // счётчик ненулевых элементов СЛАУ

						// для всех узлов расчётной области как внутренних так и внешних:
		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			integer k = istr - 1; // нумерация с нуля.
			doublereal ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];


			// Диагональный коэффициент.
			if (fabs(ap_now) > nonzeroEPS) {
				milu2.val[ik] = ap_now;
				milu2.col_ind[ik] = istr - 1; // Нуменрация начинается с нуля.
				milu2.row_ptr[k] = min(ik, milu2.row_ptr[k]);
				ik++;
			}

			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];

			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				if (bonly_negative_connections) {
					if (Amat.aij[ii1] < -nonzeroEPS) {
						milu2.val[ik] = Amat.aij[ii1]; // уже с учётом знака.
						milu2.col_ind[ik] = Amat.j[ii1] - 1; // Нумерация начинается с нуля.
						milu2.row_ptr[k] = min(ik, milu2.row_ptr[k]);
						ik++;
					}
				}
				else {
					if (fabs(Amat.aij[ii1]) > nonzeroEPS) {
						milu2.val[ik] = Amat.aij[ii1]; // уже с учётом знака.
						milu2.col_ind[ik] = Amat.j[ii1] - 1; // Нумерация начинается с нуля.
						milu2.row_ptr[k] = min(ik, milu2.row_ptr[k]);
						ik++;
					}
				}
			}
		}

		// в каждой строке элементы отсортированы по номерам столбцов:
#pragma omp parallel for
		for (integer k = 0; k < (maxelm_plus_maxbound); k++) {
			QuickSortCSIR(milu2.col_ind, milu2.val, milu2.row_ptr[k] + 1, milu2.row_ptr[k + 1] - 1);
		}

#if doubleintprecision == 1
		//printf("n==%lld row_ptr=%lld\n", n, milu2.row_ptr[maxelm_plus_maxbound]);
#else
		//printf("n==%d row_ptr=%d\n", n, milu2.row_ptr[maxelm_plus_maxbound]);
#endif

		//getchar();

		milu2.maxelm_plus_maxbound = maxelm_plus_maxbound;
#pragma omp parallel for
		for (integer i = 0; i<milu2.row_ptr[maxelm_plus_maxbound]; i++) {
			milu2.col_ind[i] = milu2.col_ind[i] + 1;
		}
#pragma omp parallel for
		for (integer i = 0; i<maxelm_plus_maxbound + 1; i++) {
			milu2.row_ptr[i] = milu2.row_ptr[i] + 1;
		}

		/*
		FILE *fp;
		errno_t err;
		// создание файла для записи.
		if ((err = fopen_s( &fp, "matr.txt", "w")) != 0) {
		printf("Create File Error\n");
		}
		else {

		// debug
		for (k=0; k<=maxelm+maxbound; k++) {
		#if doubleintprecision == 1
		fprintf(fp,"%lld ",milu2.row_ptr[k]);
		#else
		fprintf(fp,"%d ",milu2.row_ptr[k]);
		#endif

		}
		fprintf(fp,"\n");
		for (k=0; k<milu2.row_ptr[maxelm+maxbound]; k++) {
		#if doubleintprecision == 1
		fprintf(fp, "%e %lld\n",milu2.val[k],milu2.col_ind[k]);
		#else
		fprintf(fp, "%e %d\n",milu2.val[k],milu2.col_ind[k]);
		#endif

		}

		fclose(fp);
		}
		printf("ready");
		getchar();
		*/
		doublereal radd = 0.0;
		/*
		if (ilevel == 0) {
		radd = 1.8;
		}
		else if (ilevel == 1) {
		radd = 4.2;
		}
		else if ((ilevel >=2)&&(ilevel<=5)) {
		radd = 7.0;
		}
		else if (ilevel==6) {
		radd = 7.0;
		}
		else {
		radd = 8.0;
		}
		*/
		if (1) {

			// Параметры подбирал на АЛИС сетке.
			integer ilevel_loc = ilevel;
			if (ilevel == 0) ilevel_loc = 1;
			//radd = 3.3*(n / (5.0*maxelm_plus_maxbound));
			//milu2.iwk = (integer)((milu2.lfil + 1 + radd) * n + 4 * maxelm_plus_maxbound);
			//milu2.iwk = n + 4; // ilu0 не требует лишней памяти.
			integer lfil = my_amg_manager.lfil;
			if (lfil == 0) {
				milu2.iwk = (integer)((lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil==1) {
				// 2.4 +; 2.2 -; 2.3 + opt; 
				doublereal m_1 = 2.3;
				if (ilevel_loc == 10) m_1 = 4.0;
				if (ilevel_loc >= 11) m_1 = 5.0;
				if (ilevel_loc >= 15) m_1 = 8.0;
				milu2.iwk = (integer)((m_1 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 2) {
				// milu2.iwk = (40 * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound;
				// ilevel
				// lfil==1 3 is work Ok
				// milu2.iwk = (3* ilevel * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound;
				// 2.5 -; 2.6 + opt;
				doublereal m_1 = 2.6;
				if (ilevel_loc == 10) m_1 = 4.3;
				if (ilevel_loc >= 11) m_1 = 5.3;
				if (ilevel_loc >= 15) m_1 = 8.0;
				milu2.iwk = (integer)((m_1 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 3) {
				// 3 -; 3.2 + opt;
				milu2.iwk = (integer)((3.2 * ilevel_loc*lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 4) {
				// 3.8 -; 3.9 + opt;
				milu2.iwk = (integer)((3.9 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 5) {
				// 4.4 -; 
				milu2.iwk = (integer)((4.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 6) {
				milu2.iwk = (integer)((5.2 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil >= 7) {
				milu2.iwk = (integer)((5.8 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
		}

		// размерность памяти под матрицу предобуславливания.
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.alu != nullptr) {
			//delete[] milu2.alu;
			free(milu2.alu);
			milu2.alu = nullptr;
		}
		//milu2.alu = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
		milu2.alu = (doublereal*)malloc((milu2.iwk + 2) * sizeof(doublereal));
		char c1[23] = "equation3DtoCRSRUMBA1";
		char c2[11] = "milu2.alu";
		handle_error(milu2.alu, c2, c1, (milu2.iwk + 2));
		for (integer i_63 = 0; i_63 < (milu2.iwk + 2); i_63++) {
			milu2.alu[i_63] = 0.0; // Инициализация матрицы.
		}
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.jlu != nullptr) {
			//delete[] milu2.jlu;
			free(milu2.jlu);
			milu2.jlu = nullptr;
		}
		//milu2.jlu = new integer[milu2.iwk + 2];
		milu2.jlu = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
		char c3[11] = "milu2.jlu";
		handle_error(milu2.jlu, c3, c1, (milu2.iwk + 2));
		for (integer i_63 = 0; i_63 < (milu2.iwk + 2); i_63++) {
			milu2.jlu[i_63] = -1; // Инициализация матрицы. Не существующее значение столбца.
		}
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.ju != nullptr) {
			//delete[] milu2.ju;
			free(milu2.ju);
			milu2.ju = nullptr;
		}
		//milu2.ju = new integer[maxelm_plus_maxbound + 2];
		milu2.ju = (integer*)malloc((maxelm_plus_maxbound + 2) * sizeof(integer));
		char c4[10] = "milu2.ju";
		handle_error(milu2.ju, c4, c1, (maxelm_plus_maxbound + 2));
		for (integer i_63 = 0; i_63 < (maxelm_plus_maxbound + 2); i_63++) {
			milu2.ju[i_63] = -1; // Инициализация матрицы. Не существующее значение.
		}
		//milu2.levs = new integer[milu2.iwk + 2]; // уровень.
		milu2.levs = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
		char c5[12] = "milu2.levs";
		handle_error(milu2.levs, c5, c1, (milu2.iwk + 2));
		//milu2.w = new doublereal[maxelm_plus_maxbound + 2]; // +2 запас по памяти.
		milu2.w = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c6[9] = "milu2.w";
		handle_error(milu2.w, c6, c1, (maxelm_plus_maxbound + 2));
		char c7[10] = "milu2.jw";
		if (1) {
			//milu2.jw = new integer[3 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
			milu2.jw = (integer*)malloc((4 * maxelm_plus_maxbound + 2) * sizeof(integer));
			handle_error(milu2.jw, c7, c1, (4 * maxelm_plus_maxbound + 2));
		}
		else {
			if (ilevel < 6) {
				///milu2.jw = new integer[5 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
				milu2.jw = (integer*)malloc((5 * maxelm_plus_maxbound + 2) * sizeof(integer));
				handle_error(milu2.jw, c7, c1, (5 * maxelm_plus_maxbound + 2));
			}
			else {
				//milu2.jw = new integer[12 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
				milu2.jw = (integer*)malloc((12 * maxelm_plus_maxbound + 2) * sizeof(integer));
				handle_error(milu2.jw, c7, c1, (12 * maxelm_plus_maxbound + 2));
			}
		}
		if ((milu2.alu == nullptr) || (milu2.jlu == nullptr) || (milu2.levs == nullptr) || (milu2.ju == nullptr) || (milu2.w == nullptr) || (milu2.jw == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		// копии объектов
		// 9 ноября 2016 Копии объектов теперь хранятся централизованно :
		// см. LEVEL_ADDITIONAL_DATA_BUFER milu_gl_buffer; in ilut.c module.
		//milu2.alu_copy = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
		//milu2.jlu_copy = new integer[milu2.iwk + 2];
		//milu2.ju_copy = new integer[maxelm_plus_maxbound + 2];
		//milu2.b_copy = new doublereal[maxelm_plus_maxbound + 2];
		char c8[13] = "milu2.b_copy";
		milu2.b_copy = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		handle_error(milu2.b_copy, c8, c1, (maxelm_plus_maxbound + 2));
		//milu2.x_copy = new doublereal[maxelm_plus_maxbound + 2];
		milu2.x_copy = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c9[13] = "milu2.x_copy";
		handle_error(milu2.x_copy, c9, c1, (maxelm_plus_maxbound + 2));
		//milu2.zbuf = new doublereal[maxelm_plus_maxbound + 2];
		milu2.zbuf = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c10[11] = "milu2.zbuf";
		handle_error(milu2.zbuf, c10, c1, (maxelm_plus_maxbound + 2));
		//milu2.zbuf2 = new doublereal[maxelm_plus_maxbound + 2];
		milu2.zbuf2 = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c11[12] = "milu2.zbuf2";
		handle_error(milu2.zbuf2, c11, c1, (maxelm_plus_maxbound + 2));
		for (integer i_63 = 0; i_63 < (maxelm_plus_maxbound + 2); i_63++) {
			milu2.zbuf2[i_63] = 0.0;
		}

		//if ((milu2.alu_copy == nullptr) || (milu2.jlu_copy == nullptr) || (milu2.ju_copy == nullptr) || (milu2.b_copy == nullptr) || (milu2.x_copy == nullptr) || (milu2.zbuf==nullptr)||(milu2.zbuf2==nullptr)) {
		if ((milu2.b_copy == nullptr) || (milu2.x_copy == nullptr) || (milu2.zbuf == nullptr) || (milu2.zbuf2 == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		//milu2.lfil = 2;
		integer ierr = 0;
		iluk_(maxelm_plus_maxbound, milu2.val, milu2.col_ind, milu2.row_ptr, my_amg_manager.lfil/*milu2.lfil*/, milu2.alu, milu2.jlu, milu2.ju, milu2.levs, milu2.iwk, milu2.w, milu2.jw, ierr);
		if (ierr != 0) {
#if doubleintprecision == 1
			printf("ierr=%lld\n", ierr);
#else
			printf("ierr=%d\n", ierr);
#endif

			//getchar();
			system("PAUSE");
			exit(1);
		}

		if ((ierr == -2) || (ierr == -3)) {

			integer ipassage = 1;
			do {
				printf("\nPlease WAIT... ... ...\n");

				// задаче не хватило памяти, значит нужно перевыделить !
				//if (milu2.alu != nullptr) delete[] milu2.alu;
				//if (milu2.jlu != nullptr) delete[] milu2.jlu;
				//if (milu2.levs != nullptr) delete[] milu2.levs;
				if (milu2.alu != nullptr) free(milu2.alu);
				if (milu2.jlu != nullptr) free(milu2.jlu);
				if (milu2.levs != nullptr) free(milu2.levs);
				//if (milu2.alu_copy != nullptr) delete[] milu2.alu_copy;
				//if (milu2.jlu_copy != nullptr) delete[] milu2.jlu_copy;

				// инициализация !
				milu2.alu = nullptr;
				milu2.jlu = nullptr;
				milu2.levs = nullptr;
				//milu2.alu_copy = nullptr;
				//milu2.jlu_copy = nullptr;
				integer ilevel_loc = ilevel;
				if (ilevel == 0) ilevel_loc = 1;
				integer lfil = my_amg_manager.lfil;
				if (lfil == 0) {
					milu2.iwk = (integer)((lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 1) {
					milu2.iwk = (integer)((2.3* ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 2) {
					milu2.iwk = (integer)((2.6* ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 3) {
					milu2.iwk = (integer)((3.2* ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 4) {
					milu2.iwk = (integer)((3.9 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 5) {
					milu2.iwk = (integer)((4.4 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 6) {
					milu2.iwk = (integer)((5.0 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil >= 7) {
					milu2.iwk = (integer)((5.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				
				// Выделение оперативной памяти может быть 
		        // повторным, поэтому необходимо освободить оперативную память,
		        // перед повторным выделением.
				if (milu2.alu != nullptr) {
					//delete[] milu2.alu;
					free(milu2.alu);
					milu2.alu = nullptr;
				}
				//milu2.alu = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
				milu2.alu = (doublereal*)malloc((milu2.iwk + 2) * sizeof(doublereal));
				char c12[11] = "milu2.alu";
				handle_error(milu2.alu, c12, c1, (milu2.iwk + 2));
				// Выделение оперативной памяти может быть 
		        // повторным, поэтому необходимо освободить оперативную память,
	        	// перед повторным выделением.
				if (milu2.jlu != nullptr) {
					//delete[] milu2.jlu;
					free(milu2.jlu);
					milu2.jlu = nullptr;
				}
				//milu2.jlu = new integer[milu2.iwk + 2];
				milu2.jlu = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
				char c13[11] = "milu2.jlu";
				handle_error(milu2.jlu, c13, c1, (milu2.iwk + 2));
				//milu2.levs = new integer[milu2.iwk + 2]; // уровень.
				milu2.levs = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
				char c14[11] = "milu2.levs";
				handle_error(milu2.levs, c14, c1, (milu2.iwk + 2));
				//milu2.alu_copy = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
				//milu2.jlu_copy = new integer[milu2.iwk + 2];

				//if ((milu2.alu_copy != nullptr) && (milu2.jlu_copy != nullptr) && (milu2.alu != nullptr) && (milu2.jlu != nullptr) && (milu2.levs != nullptr)) {
				if ((milu2.alu != nullptr) && (milu2.jlu != nullptr) && (milu2.levs != nullptr)) {
					iluk_(maxelm_plus_maxbound, milu2.val, milu2.col_ind, milu2.row_ptr, lfil, milu2.alu, milu2.jlu, milu2.ju, milu2.levs, milu2.iwk, milu2.w, milu2.jw, ierr);
				}
				else {
					// недостаточно памяти на данном оборудовании.
					ipassage = 4;
					printf("Problem : not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

				ipassage++;
				// Пока просто контролируем, в добавок хорошо бы выводить информацию на каком уровне это произошло.
#if doubleintprecision == 1
				printf("control memory allocation for ilu2 in RUMBA0.14 in level=%lld\n", ilevel);
#else
				printf("control memory allocation for ilu2 in RUMBA0.14 in level=%d\n", ilevel);
#endif

				system("PAUSE");
			} while ((ierr != 0) && (ipassage<4));

			if (ipassage == 4) {
				printf("Error memory alloc !!!\n");
				printf("failed to obtain an expansion for the 4 approaches...\n");
				printf("Please any key to exit...\n");
				//getchar();
				system("PAUSE");
				exit(1);
			}
		}

		for (integer k = 0; k<n; k++) if (milu2.col_ind[k] == (-1)) {
			printf("Error equation3D to CRS.\n");
			//getchar();
			system("pause");
		}

		// Нумерация начинается с единицы (это важно).
		for (integer k = 0; k<maxelm_plus_maxbound; k++) {
			if (milu2.val[milu2.row_ptr[k] - 1]<nonzeroEPS) {
#if doubleintprecision == 1
				printf("negativ diagonal element equation3DtoCRS %lld\n", k);
#else
				printf("negativ diagonal element equation3DtoCRS %d\n", k);
#endif

				//getchar();
				system("pause");
			}
		}

		if (1 || milu2.lfil == 0) {
			// экономим оперативную память.
			//if (milu2.val != nullptr) delete[] milu2.val;
			//if (milu2.col_ind != nullptr) delete[] milu2.col_ind;
			//if (milu2.row_ptr != nullptr) delete[] milu2.row_ptr;
			if (milu2.val != nullptr) free(milu2.val);
			if (milu2.col_ind != nullptr) free(milu2.col_ind);
			if (milu2.row_ptr != nullptr) free(milu2.row_ptr);
			milu2.val = nullptr;
			milu2.col_ind = nullptr;
			milu2.row_ptr = nullptr;
		}

	}

	if (!flag) {
		printf("Error equation 3D to CRS: zero diagonal element...\n");
		//getchar();
		system("pause");
	}


} // equation3DtoCRSRUMBA1


// 31.12.2019 iluk сглаживатель для amg1r5 алгоритма. Окончание 1.01.2020.
// 15.11.2016 метод полностью на malloc с контролем через handle_error.
// Мы усилим РУМБА0.14 алгоритм ILU2 предобуславливанием на каждом уровне вложенности.
// Преобразует equation3D  формат хранения в CRS формат.
// Цель написания этого преобразователя: экономия оперативной памяти компьютера.
// Т.к. формат SIMPLESPARSE требует слишком много памяти.
void equation3DtoCRSRUMBA1_amg1r5(LEVEL_ADDITIONAL_DATA& milu2,
	doublereal* a, integer* ia, integer* ja,
	integer* imin, integer* imax, integer ilevel) {

	// Если ballocmemory равен true то происходит выделение памяти.
	const bool bonly_negative_connections = false;
	
	integer n = 0; // число ненулевых элементов
	integer maxelm_plus_maxbound = 0; // Число уравнений в СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	integer istartq = imin[ilevel + 1];
	integer iendq = imax[ilevel + 1];
	maxelm_plus_maxbound = 1 + imax[ilevel + 1] - imin[ilevel + 1];
	n = ia[imax[ilevel + 1] + 1] - 1 - ia[imin[ilevel + 1]] + 1;
	

	integer startpos = istartq;
	integer endpos = iendq;

	bool ballocmemory = false;
	//if (milu2.jw == nullptr) ballocmemory = true;
	if (bflag_visit_amg1r5[ilevel] == false) {
		ballocmemory = true;
		bflag_visit_amg1r5[ilevel] = true;
	}
	
	// memory +15N
	// Теперь выделение памяти будет происходить централизованно, вне данного кода.
	// Это сделано для кода BICGSTAB_internal3. дата изменения 12 апреля 2013.
	// Другой код, использующий equation3dtoCRS может оказаться неработоспособным после этого изменения.
	if (ballocmemory) {

		// Замер времени.
		unsigned int calculation_main_start_time = 0; // начало счёта мс.
		unsigned int calculation_main_end_time = 0; // окончание счёта мс.


		calculation_main_start_time = clock(); // момент начала счёта.

		printf("apply ilu%lld smoother for number ", my_amg_manager.lfil);
		printf("%2lld level %3d ", ilevel, (int)(n/(maxelm_plus_maxbound)));

		// Важно выделить память с запасом, т.к. одна и таже память используется и для компонент скорости и для попрапвки давления.
		//val = new doublereal[7 * (maxelm + maxbound) + 2 * maxbound + 2];
		//col_ind = new integer[7 * (maxelm + maxbound) + 2 * maxbound + 2];
		
		//milu2.val = new doublereal[n + 2];
		milu2.val = (doublereal*)malloc((n + 2) * sizeof(doublereal));
		char c1[11] = "milu2.val";
		char c2[23] = "equation3DtoCRSRUMBA1";
		handle_error(milu2.val, c1, c2, (n + 2));
		//milu2.col_ind = new integer[n + 2];
		milu2.col_ind = (integer*)malloc((n + 2) * sizeof(integer));
		char c3[14] = "milu2.col_ind";
		handle_error(milu2.col_ind, c3, c2, (n + 2));
		//row_ptr = new integer[(maxelm + maxbound) + 1];
		//milu2.row_ptr = new integer[maxelm_plus_maxbound + 1];
		milu2.row_ptr = (integer*)malloc((maxelm_plus_maxbound + 1) * sizeof(integer));
		char c4[14] = "milu2.row_ptr";
		handle_error(milu2.row_ptr, c4, c2, (maxelm_plus_maxbound + 1));

		if ((milu2.val == nullptr) || (milu2.col_ind == nullptr) || (milu2.row_ptr == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		integer iminus = 0;
		if (ilevel > 0) iminus = ja[ia[imin[ilevel + 1]]]-1;

		integer k = 0, k_1 = 0;
		integer i__1 = imax[ilevel + 1];
		for (integer i__ = imin[ilevel + 1]; i__ <= i__1; ++i__) {
			milu2.row_ptr[k_1] = k;
			k_1++;

			integer i__2 = ia[i__ + 1] - 1;
			for (integer j = ia[i__]; j <= i__2; ++j) {
				milu2.val[k] = a[j];
				milu2.col_ind[k] = ja[j] -1 -iminus;
				k++;
			}
		}
		milu2.row_ptr[k_1] = k; // n==k гарантированно.
	

		// в каждой строке элементы отсортированы по номерам столбцов:
#pragma omp parallel for
		for (integer k = 0; k < (maxelm_plus_maxbound); k++) {
			QuickSortCSIR(milu2.col_ind, milu2.val, milu2.row_ptr[k] + 1, milu2.row_ptr[k + 1] - 1);
		}

		milu2.maxelm_plus_maxbound = maxelm_plus_maxbound;
		
#pragma omp parallel for
		for (integer i = 0; i < milu2.row_ptr[maxelm_plus_maxbound]; i++) {
			milu2.col_ind[i] = milu2.col_ind[i] + 1;
		}
#pragma omp parallel for
		for (integer i = 0; i < maxelm_plus_maxbound + 1; i++) {
			milu2.row_ptr[i] = milu2.row_ptr[i] + 1;
		}
		
			
		doublereal radd = 0.0;
			
		if (1) {

			// Параметры подбирал на АЛИС сетке.
			integer ilevel_loc = ilevel;
			if (ilevel == 0) ilevel_loc = 1;
		   	//radd = 3.3*(n / (5.0*maxelm_plus_maxbound));
		    //milu2.iwk = (integer)((milu2.lfil + 1 + radd) * n + 4 * maxelm_plus_maxbound);
			//milu2.iwk = n + 4; // ilu0 не требует лишней памяти.
			integer lfil = my_amg_manager.lfil;
			if (lfil == 0) {
				milu2.iwk = (integer)((lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 1) {
				// 2.4 +; 2.2 -; 2.3 + opt; 
				doublereal m_1 = 2.3;
				if (ilevel_loc == 10) m_1 = 4.0;
				if (ilevel_loc >= 11) m_1 = 5.0;
				if (ilevel_loc >= 15) m_1 = 8.0;
				milu2.iwk = (integer)((m_1 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 2) {
				// milu2.iwk = (40 * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound;
				// ilevel
				// lfil==1 3 is work Ok
				// milu2.iwk = (3* ilevel * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound;
				// 2.5 -; 2.6 + opt;
				doublereal m_1 = 2.6;
				if (ilevel_loc == 10) m_1 = 4.3;
				if (ilevel_loc >= 11) m_1 = 5.3;
				if (ilevel_loc >= 15) m_1 = 8.0;
				milu2.iwk = (integer)((m_1 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 3) {
				// 3 -; 3.2 + opt;
				milu2.iwk = (integer)((3.2 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 4) {
				// 3.8 -; 3.9 + opt;
				milu2.iwk = (integer)((3.9 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 5) {
				// 4.4 -; 
				milu2.iwk = (integer)((4.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil == 6) {
				milu2.iwk = (integer)((5.2 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
			else if (lfil >= 7) {
				milu2.iwk = (integer)((5.8 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + 4 * maxelm_plus_maxbound);
			}
		}

		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.alu != nullptr) {
			//delete[] milu2.alu;
			free(milu2.alu);
			milu2.alu = nullptr;
		}
		// размерность памяти под матрицу предобуславливания.
		//milu2.alu = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
		milu2.alu = (doublereal*)malloc((milu2.iwk + 2) * sizeof(doublereal));
		char c11[23] = "equation3DtoCRSRUMBA1";
		char c21[11] = "milu2.alu";
		handle_error(milu2.alu, c21, c11, (milu2.iwk + 2));
		// Инициализация матрицы нулевыми значениями.
		for (integer i_63 = 0; i_63 < (milu2.iwk + 2); i_63++) {
			milu2.alu[i_63] = 0.0;
		}
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.jlu != nullptr) {
			//delete[] milu2.jlu;
			free(milu2.jlu);
			milu2.jlu = nullptr;
		}
		//milu2.jlu = new integer[milu2.iwk + 2];
		milu2.jlu = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
		char c31[11] = "milu2.jlu";
		handle_error(milu2.jlu, c31, c11, (milu2.iwk + 2));
		// Инициализация матрицы несуществующими номерами столбцов.
		for (integer i_63 = 0; i_63 < (milu2.iwk + 2); i_63++) {
			milu2.jlu[i_63] = -1; // Столбец с номером -1 не существует в матрице.
		}
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.ju != nullptr) {
			//delete[] milu2.ju;
			free(milu2.ju);
			milu2.ju = nullptr;
		}
		//milu2.ju = new integer[maxelm_plus_maxbound + 2];
		milu2.ju = (integer*)malloc((maxelm_plus_maxbound + 2) * sizeof(integer));
		char c41[10] = "milu2.ju";
		handle_error(milu2.ju, c41, c11, (maxelm_plus_maxbound + 2));
		// Инициализация матрицы несуществующими номерами столбцов.
		for (integer i_63 = 0; i_63 < (maxelm_plus_maxbound + 2); i_63++) {
			milu2.ju[i_63] = -1; // Значение -1 не существует в матрице.
		}
		//milu2.levs = new integer[milu2.iwk + 2]; // уровень.
		milu2.levs = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
		char c5[12] = "milu2.levs";
		handle_error(milu2.levs, c5, c1, (milu2.iwk + 2));
		//milu2.w = new doublereal[maxelm_plus_maxbound + 2]; // +2 запас по памяти.
		milu2.w = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c6[9] = "milu2.w";
		handle_error(milu2.w, c6, c1, (maxelm_plus_maxbound + 2));
		char c7[10] = "milu2.jw";
		if (1) {
			//milu2.jw = new integer[3 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
			milu2.jw = (integer*)malloc((4 * maxelm_plus_maxbound + 2) * sizeof(integer));
			handle_error(milu2.jw, c7, c1, (4 * maxelm_plus_maxbound + 2));
		}
		else {
			if (ilevel < 6) {
				///milu2.jw = new integer[5 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
				milu2.jw = (integer*)malloc((5 * maxelm_plus_maxbound + 2) * sizeof(integer));
				handle_error(milu2.jw, c7, c1, (5 * maxelm_plus_maxbound + 2));
			}
			else {
				//milu2.jw = new integer[12 * maxelm_plus_maxbound + 2]; // +2 запас по памяти.
				milu2.jw = (integer*)malloc((12 * maxelm_plus_maxbound + 2) * sizeof(integer));
				handle_error(milu2.jw, c7, c1, (12 * maxelm_plus_maxbound + 2));
			}
		}
		if ((milu2.alu == nullptr) || (milu2.jlu == nullptr) || (milu2.levs == nullptr) || (milu2.ju == nullptr) || (milu2.w == nullptr) || (milu2.jw == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		// копии объектов
		// 9 ноября 2016 Копии объектов теперь хранятся централизованно :
		// см. LEVEL_ADDITIONAL_DATA_BUFER milu_gl_buffer; in ilut.c module.
		//milu2.alu_copy = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
		//milu2.jlu_copy = new integer[milu2.iwk + 2];
		//milu2.ju_copy = new integer[maxelm_plus_maxbound + 2];
		//milu2.b_copy = new doublereal[maxelm_plus_maxbound + 2];
		char c8[13] = "milu2.b_copy";
		milu2.b_copy = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		handle_error(milu2.b_copy, c8, c1, (maxelm_plus_maxbound + 2));
		//milu2.x_copy = new doublereal[maxelm_plus_maxbound + 2];
		milu2.x_copy = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c9[13] = "milu2.x_copy";
		handle_error(milu2.x_copy, c9, c1, (maxelm_plus_maxbound + 2));
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.zbuf != nullptr) {
			//delete[] milu2.zbuf;
			free(milu2.zbuf);
			milu2.zbuf = nullptr;
		}
		//milu2.zbuf = new doublereal[maxelm_plus_maxbound + 2];
		milu2.zbuf = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c10[11] = "milu2.zbuf";
		handle_error(milu2.zbuf, c10, c1, (maxelm_plus_maxbound + 2));
		// Выделение оперативной памяти может быть 
		// повторным, поэтому необходимо освободить оперативную память,
		// перед повторным выделением.
		if (milu2.zbuf2 != nullptr) {
			//delete[] milu2.zbuf2;
			free(milu2.zbuf2);
			milu2.zbuf2 = nullptr;
		}
		//milu2.zbuf2 = new doublereal[maxelm_plus_maxbound + 2];
		milu2.zbuf2 = (doublereal*)malloc((maxelm_plus_maxbound + 2) * sizeof(doublereal));
		char c111[12] = "milu2.zbuf2";
		handle_error(milu2.zbuf2, c111, c1, (maxelm_plus_maxbound + 2));
		for (integer i_63 = 0; i_63 < (maxelm_plus_maxbound + 2); i_63++) {
			milu2.zbuf2[i_63] = 0.0;
		}

		//if ((milu2.alu_copy == nullptr) || (milu2.jlu_copy == nullptr) || (milu2.ju_copy == nullptr) || (milu2.b_copy == nullptr) || (milu2.x_copy == nullptr) || (milu2.zbuf==nullptr)||(milu2.zbuf2==nullptr)) {
		if ((milu2.b_copy == nullptr) || (milu2.x_copy == nullptr) || (milu2.zbuf == nullptr) || (milu2.zbuf2 == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		//milu2.lfil = 2;
		integer ierr = 0;
		iluk_(maxelm_plus_maxbound, milu2.val, milu2.col_ind, milu2.row_ptr, my_amg_manager.lfil/*milu2.lfil*/, milu2.alu, milu2.jlu, milu2.ju, milu2.levs, milu2.iwk, milu2.w, milu2.jw, ierr);
		if (ierr != 0) {
#if doubleintprecision == 1
		    printf("ierr=%lld\n", ierr);
#else
		    printf("ierr=%d\n", ierr);
#endif
			//getchar();
			system("PAUSE");
			exit(1);
		}

		if ((ierr == -2) || (ierr == -3)) {

			integer ipassage = 1;
			do {
				printf("\nPlease WAIT... ... ...\n");

				// задаче не хватило памяти, значит нужно перевыделить !
				//if (milu2.alu != nullptr) delete[] milu2.alu;
				//if (milu2.jlu != nullptr) delete[] milu2.jlu;
		     	//if (milu2.levs != nullptr) delete[] milu2.levs;
				if (milu2.alu != nullptr) free(milu2.alu);
				if (milu2.jlu != nullptr) free(milu2.jlu);
				if (milu2.levs != nullptr) free(milu2.levs);
				//if (milu2.alu_copy != nullptr) delete[] milu2.alu_copy;
				//if (milu2.jlu_copy != nullptr) delete[] milu2.jlu_copy;

				// инициализация !
				milu2.alu = nullptr;
				milu2.jlu = nullptr;
				milu2.levs = nullptr;
				//milu2.alu_copy = nullptr;
				//milu2.jlu_copy = nullptr;
				integer ilevel_loc = ilevel;
				if (ilevel == 0) ilevel_loc = 1;
				integer lfil = my_amg_manager.lfil;
				if (lfil == 0) {
					milu2.iwk = (integer)((lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 1) {
					milu2.iwk = (integer)((2.3 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 2) {
					milu2.iwk = (integer)((2.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 3) {
					milu2.iwk = (integer)((3.2 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 4) {
					milu2.iwk = (integer)((3.9 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 5) {
					milu2.iwk = (integer)((4.4 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil == 6) {
					milu2.iwk = (integer)((5.0 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}
				else if (lfil >= 7) {
					milu2.iwk = (integer)((5.6 * ilevel_loc * lfil + 1) * (n + 2 * maxelm_plus_maxbound + 2) + (4 + 3 * ipassage) * maxelm_plus_maxbound);
				}

				//milu2.alu = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
				milu2.alu = (doublereal*)malloc((milu2.iwk + 2) * sizeof(doublereal));
				char c12[11] = "milu2.alu";
				handle_error(milu2.alu, c12, c1, (milu2.iwk + 2));
				//milu2.jlu = new integer[milu2.iwk + 2];
				milu2.jlu = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
				char c13[11] = "milu2.jlu";
				handle_error(milu2.jlu, c13, c1, (milu2.iwk + 2));
				//milu2.levs = new integer[milu2.iwk + 2]; // уровень.
				milu2.levs = (integer*)malloc((milu2.iwk + 2) * sizeof(integer));
				char c14[11] = "milu2.levs";
				handle_error(milu2.levs, c14, c1, (milu2.iwk + 2));
				//milu2.alu_copy = new doublereal[milu2.iwk + 2]; // +2 запас по памяти.
				//milu2.jlu_copy = new integer[milu2.iwk + 2];

				//if ((milu2.alu_copy != nullptr) && (milu2.jlu_copy != nullptr) && (milu2.alu != nullptr) && (milu2.jlu != nullptr) && (milu2.levs != nullptr)) {
				if ((milu2.alu != nullptr) && (milu2.jlu != nullptr) && (milu2.levs != nullptr)) {
					iluk_(maxelm_plus_maxbound, milu2.val, milu2.col_ind, milu2.row_ptr, lfil, milu2.alu, milu2.jlu, milu2.ju, milu2.levs, milu2.iwk, milu2.w, milu2.jw, ierr);
				}
				else {
					// недостаточно памяти на данном оборудовании.
					ipassage = 4;
					printf("Problem : not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

				ipassage++;
				// Пока просто контролируем, в добавок хорошо бы выводить информацию на каком уровне это произошло.
#if doubleintprecision == 1
				printf("control memory allocation for ilu2 in RUMBA0.14 in level=%lld\n", ilevel);
#else
				printf("control memory allocation for ilu2 in RUMBA0.14 in level=%d\n", ilevel);
#endif

				system("PAUSE");
				
			} while ((ierr != 0) && (ipassage < 4));

			if (ipassage == 4) {
				printf("Error memory alloc !!!\n");
				printf("failed to obtain an expansion for the 4 approaches...\n");
			    printf("Please any key to exit...\n");
				//getchar();
				system("PAUSE");
				exit(1);
			}
		}

		for (integer k = 0; k < n; k++) if (milu2.col_ind[k] == (-1)) {
			printf("Error equation3D to CRS.\n");
			//getchar();
			system("pause");
		}

		// Нумерация начинается с единицы (это важно).
		for (integer k = 0; k < maxelm_plus_maxbound; k++) {
			if (milu2.val[milu2.row_ptr[k] - 1] < nonzeroEPS) {
#if doubleintprecision == 1
		    	printf("negativ diagonal element equation3DtoCRS %lld\n", k);
#else
				printf("negativ diagonal element equation3DtoCRS %d\n", k);
#endif

				//getchar();
				system("pause");
		    }
	    }

		if (1 || milu2.lfil == 0) {
			// экономим оперативную память.
			//if (milu2.val != nullptr) delete[] milu2.val;
			//if (milu2.col_ind != nullptr) delete[] milu2.col_ind;
			//if (milu2.row_ptr != nullptr) delete[] milu2.row_ptr;
			if (milu2.val != nullptr) free(milu2.val);
			if (milu2.col_ind != nullptr) free(milu2.col_ind);
			if (milu2.row_ptr != nullptr) free(milu2.row_ptr);
			milu2.val = nullptr;
			milu2.col_ind = nullptr;
			milu2.row_ptr = nullptr;
		}

		calculation_main_end_time = clock(); // момент окончания счёта.
		printf("%6d \n", calculation_main_end_time - calculation_main_start_time);
    }

} // equation3DtoCRSRUMBA1_amg1r5

  // Мы усилим РУМБА0.14 алгоритм ILU0 предобуславливанием (сглаживателем) на каждом уровне вложенности.
  // Преобразует equation3D  формат хранения в CRS формат.
  // Цель написания этого преобразователя: экономия оперативной памяти компьютера.
  // Т.к. формат SIMPLESPARSE требует слишком много памяти.
  // Для многих задач сходимость мульгрида на базе Зейдель сглаживателя недостаточна. 
  // Судя по документации ILU0 smoother more effective по сравнению с Зейдель сглаживателем.
  // ILU0 - универсален и не накладывает ограничений на структуру матрицы СЛАУ. ILU0 потребляет мало пользователей.
void equation3DtoCRSRUMBA0(LEVEL_ADDITIONAL_DATA0 &milu0,
	bool ballocmemory, Ak2 &Amat, integer istartq, integer iendq,
	integer* &row_ptr_start, integer* &row_ptr_end, integer iadd, integer ilevel) {


	// iadd_now=n_a[0]+...+n_a[ilevel_detector-1];
	//equation3DtoCRSRUMBA0(milu0[ilevel_detector], true,
	//Amat, 1, n_a[ilevel_detector], row_ptr_start, row_ptr_end, iadd_now, ilevel_detector);

	// Если ballocmemory равен true то происходит выделение памяти.
	const bool bonly_negative_connections = false;
	bool flag = true;
	integer n = 0; // число ненулевых элементов
	integer maxelm_plus_maxbound = 0; // Число уравнений в СЛАУ.

	const doublereal nonzeroEPS = 1e-37; // для отделения вещественного нуля

	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;


	// подсчёт количества ненулевых элементов
	// во всех точках расчётной области как внутренних так и граничных.
	for (integer ii = startpos; ii <= endpos; ii++) {
		maxelm_plus_maxbound++;
		integer istr = ii - iadd;
		if (fabs(Amat.aij[row_ptr_start[ii]]) < 1.0e-20)
		{
#if doubleintprecision == 1
			printf("zero diagonal coefficient in level = %lld, istr=%lld\n", ilevel, istr);
#else
			printf("zero diagonal coefficient in level = %d, istr=%d\n", ilevel, istr);
#endif

			//getchar();
			system("PAUSE");
			exit(1);
		}
		doublereal ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

		//if (fabs(ap_now)> 1e10*nonzeroEPS) n++; // диагональный элемент
		if (ap_now > nonzeroEPS) n++; // Диагональный элемент.
		else {
			// 5 августа 2016. 
			flag = false;
			printf("internal zero diagonal element.\n");
			//printf("ae=%e aw=%e an=%e as=%e at=%e ab=%e sum_nb=%e", sl[k].ae, sl[k].aw, sl[k].an, sl[k].as, sl[k].at, sl[k].ab, sl[k].ae + sl[k].aw + sl[k].an + sl[k].as + sl[k].at + sl[k].ab);
			if (ap_now < 0.0) {
				printf("found negativ diagonal coefficient=%e...\n", ap_now);
			}
			printf("fatal error equation3DtoCRS in RUMBA...\n");
			//getchar();
			system("PAUSE");
			exit(1);
			//n++;
			//sl[k].ap = fabs(sl[k].ae) + fabs(sl[k].aw) + fabs(sl[k].an) + fabs(sl[k].as) + fabs(sl[k].at) + fabs(sl[k].ab);
		}

		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];

		for (integer ii1 = is1; ii1 <= is2; ii1++)
		{
			if (bonly_negative_connections) {
				if (Amat.aij[ii1] < -nonzeroEPS) {
					n++;
				}
			}
			else {
				if (fabs(Amat.aij[ii1]) > nonzeroEPS) {
					n++;
				}
				else {
#if doubleintprecision == 1
					printf("igmoring a[%lld][%lld]=%e istr=%lld\n", Amat.i[ii1], Amat.j[ii1], fabs(Amat.aij[ii1]), istr);
#else
					printf("igmoring a[%d][%d]=%e istr=%d\n", Amat.i[ii1], Amat.j[ii1], fabs(Amat.aij[ii1]), istr);
#endif

				}
			}
		}
	}

#if doubleintprecision == 1
	printf("nnz=%lld, %lld startpos=%lld endpos=%lld\n", n, row_ptr_end[endpos] - row_ptr_start[startpos] + 1, startpos, endpos);
#else
	printf("nnz=%d, %d startpos=%d endpos=%d\n", n, row_ptr_end[endpos] - row_ptr_start[startpos] + 1, startpos, endpos);
#endif

	system("PAUSE");

	if (flag) {
		// memory +15N
		// Теперь выделение памяти будет происходить централизованно, вне данного кода.
		// Это сделано для кода BICGSTAB_internal3. дата изменения 12 апреля 2013.
		// Другой код, использующий equation3dtoCRS может оказаться неработоспособным после этого изменения.
		if (ballocmemory) {
			// Важно выделить память с запасом, т.к. одна и таже память используется и для компонент скорости и для попрапвки давления.
			//val = new doublereal[7 * (maxelm + maxbound) + 2 * maxbound + 2];
			//col_ind = new integer[7 * (maxelm + maxbound) + 2 * maxbound + 2];
			milu0.val = new doublereal[n + 2];
			milu0.col_ind = new integer[n + 2];
			//row_ptr = new integer[(maxelm + maxbound) + 1];
			milu0.row_ptr = new integer[maxelm_plus_maxbound + 1];
			if ((milu0.val == nullptr) || (milu0.col_ind == nullptr) || (milu0.row_ptr == nullptr)) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for ILU0 decomposition...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
		}


		// инициализация
#pragma omp parallel for
		for (integer k = 0; k<(n); k++) {
			milu0.val[k] = 0.0;
			milu0.col_ind[k] = -1;
		}
#pragma omp parallel for
		for (integer k = 0; k <= (maxelm_plus_maxbound); k++) {
			milu0.row_ptr[k] = n; // присваиваем количество ненулевых элементов плюс 1 с учётом того что нумерация массива начинается с 0
		}

		// Быстрая Сортировка Хоара.
		// упорядочивание по строкам
		//QuickSort(...); не требуется,
		// т.к. сама структура хранения 
		// подразумевает упорядочивание по строкам.

		/*
		// заполнение разреженной матрицы
		for (integer k=0; k<M.n; k++) {
		val[k]=M.a[k].aij;
		col_ind[k]=M.a[k].j;
		row_ptr[M.a[k].i]=min(k,row_ptr[M.a[k].i]);
		}
		*/
		integer ik = 0; // счётчик ненулевых элементов СЛАУ

						// для всех узлов расчётной области как внутренних так и внешних:
		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			integer k = istr - 1; // нумерация с нуля.
			doublereal ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];


			// Диагональный коэффициент.
			if (fabs(ap_now) > nonzeroEPS) {
				milu0.val[ik] = ap_now;
				milu0.col_ind[ik] = istr - 1; // Нуменрация начинается с нуля.
				milu0.row_ptr[k] = min(ik, milu0.row_ptr[k]);
				ik++;
			}

			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];

			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				if (bonly_negative_connections) {
					if (Amat.aij[ii1] < -nonzeroEPS) {
						milu0.val[ik] = Amat.aij[ii1]; // уже с учётом знака.
						milu0.col_ind[ik] = Amat.j[ii1] - 1; // Нумерация начинается с нуля.
						milu0.row_ptr[k] = min(ik, milu0.row_ptr[k]);
						ik++;
					}
				}
				else {
					if (fabs(Amat.aij[ii1]) > nonzeroEPS) {
						milu0.val[ik] = Amat.aij[ii1]; // уже с учётом знака.
						milu0.col_ind[ik] = Amat.j[ii1] - 1; // Нумерация начинается с нуля.
						milu0.row_ptr[k] = min(ik, milu0.row_ptr[k]);
						ik++;
					}
				}
			}
		}

		// в каждой строке элементы отсортированы по номерам столбцов:
#pragma omp parallel for
		for (integer k = 0; k<(maxelm_plus_maxbound); k++) QuickSortCSIR(milu0.col_ind, milu0.val, milu0.row_ptr[k] + 1, milu0.row_ptr[k + 1] - 1);

#if doubleintprecision == 1
		//printf("n==%lld row_ptr=%lld\n", n, milu0.row_ptr[maxelm_plus_maxbound]);
#else
		//printf("n==%d row_ptr=%d\n", n, milu0.row_ptr[maxelm_plus_maxbound]);
#endif

		//getchar();
#pragma omp parallel for
		for (integer i = 0; i<milu0.row_ptr[maxelm_plus_maxbound]; i++) {
			milu0.col_ind[i] = milu0.col_ind[i] + 1;
		}
#pragma omp parallel for
		for (integer i = 0; i<maxelm_plus_maxbound + 1; i++) {
			milu0.row_ptr[i] = milu0.row_ptr[i] + 1;
		}

		/*
		FILE *fp;
		errno_t err;
		// создание файла для записи.
		if ((err = fopen_s( &fp, "matr.txt", "w")) != 0) {
		printf("Create File Error\n");
		}
		else {
		#if doubleintprecision == 1
		// debug
		for (k=0; k<=maxelm+maxbound; k++) {
		fprintf(fp,"%lld ",milu0.row_ptr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<milu0.row_ptr[maxelm+maxbound]; k++) {
		fprintf(fp, "%e %lld\n",milu0.val[k],milu0.col_ind[k]);
		}
		#else
		// debug
		for (k=0; k<=maxelm+maxbound; k++) {
		fprintf(fp,"%d ",milu0.row_ptr[k]);
		}
		fprintf(fp,"\n");
		for (k=0; k<milu0.row_ptr[maxelm+maxbound]; k++) {
		fprintf(fp, "%e %d\n",milu0.val[k],milu0.col_ind[k]);
		}
		#endif



		fclose(fp);
		}
		printf("ready");
		getchar();
		*/
		doublereal radd = 0.0;
		/*
		if (ilevel == 0) {
		radd = 1.8;
		}
		else if (ilevel == 1) {
		radd = 4.2;
		}
		else if ((ilevel >= 2) && (ilevel <= 5)) {
		radd = 6.5;
		}
		else if (ilevel == 6) {
		radd = 7.0;
		}
		else {
		radd = 8.0;
		}
		*/

		// ILU0 не требует памяти больше чем размер исходных матриц. 
		// Новых ненулевых элементов для него не появляется.
		milu0.iwk = milu0.row_ptr[maxelm_plus_maxbound];
		milu0.alu = new doublereal[milu0.row_ptr[maxelm_plus_maxbound] + 2];
		milu0.jlu = new integer[milu0.row_ptr[maxelm_plus_maxbound] + 2];

		milu0.ju = new integer[maxelm_plus_maxbound + 2];


		milu0.iw = new integer[maxelm_plus_maxbound + 2]; // +2 запас по памяти.

		if ((milu0.alu == nullptr) || (milu0.jlu == nullptr) || (milu0.ju == nullptr) || (milu0.iw == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for ILU0 decomposition: alu, jlu, ju, iw...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		// копии объектов

		// ILU0 не требует памяти больше чем размер исходных матриц. 
		// Новых ненулевых элементов для него не появляется.
		milu0.alu_copy = new doublereal[milu0.row_ptr[maxelm_plus_maxbound] + 2];
		milu0.jlu_copy = new integer[milu0.row_ptr[maxelm_plus_maxbound] + 2];

		milu0.ju_copy = new integer[maxelm_plus_maxbound + 2];

		milu0.b_copy = new doublereal[maxelm_plus_maxbound + 2];
		milu0.x_copy = new doublereal[maxelm_plus_maxbound + 2];
		milu0.zbuf = new doublereal[maxelm_plus_maxbound + 2];
		milu0.zbuf2 = new doublereal[maxelm_plus_maxbound + 2];

		if ((milu0.alu_copy == nullptr) || (milu0.jlu_copy == nullptr) || (milu0.ju_copy == nullptr) || (milu0.b_copy == nullptr) || (milu0.x_copy == nullptr) || (milu0.zbuf == nullptr) || (milu0.zbuf2 == nullptr)) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for ILU0 decomposition: alu_copy, jlu_copy, ju_copy etc...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		// Построение ILU0 декомпозиции.
		integer ierr = 0;
		ilu0_(maxelm_plus_maxbound, milu0.val, milu0.col_ind, milu0.row_ptr, milu0.alu, milu0.jlu, milu0.ju, milu0.iw, ierr);
		if (ierr != 0) {
#if doubleintprecision == 1
			printf("ierr=%lld\n", ierr);
#else
			printf("ierr=%d\n", ierr);
#endif

			//getchar();
			system("PAUSE");
			exit(1);
		}

		if ((ierr == -2) || (ierr == -3)) {

			integer ipassage = 1;
			do {
				printf("\nPlease WAIT... ... ...\n");

				// задаче не хватило памяти, значит нужно перевыделить !
				if (milu0.alu != nullptr) delete milu0.alu;
				if (milu0.jlu != nullptr) delete milu0.jlu;
				if (milu0.alu_copy != nullptr) delete milu0.alu_copy;
				if (milu0.jlu_copy != nullptr) delete milu0.jlu_copy;

				// инициализация !
				milu0.alu = nullptr;
				milu0.jlu = nullptr;
				milu0.alu_copy = nullptr;
				milu0.jlu_copy = nullptr;

				// При использовании ILU0 новой памяти не требуется должно хватать размера исходной матрицы.
				integer iwk = milu0.row_ptr[maxelm_plus_maxbound] + 2; // ИЗМЕНИТЬ ЭТУ СТРОКУ для увеличения памяти.
				milu0.iwk = iwk - 2;
				milu0.alu = new doublereal[iwk + 2]; // +2 запас по памяти.
				milu0.jlu = new integer[iwk + 2];
				milu0.alu_copy = new doublereal[iwk + 2]; // +2 запас по памяти.
				milu0.jlu_copy = new integer[iwk + 2];

				if ((milu0.alu_copy != nullptr) && (milu0.jlu_copy != nullptr) && (milu0.alu != nullptr) && (milu0.jlu != nullptr)) {
					ilu0_(maxelm_plus_maxbound, milu0.val, milu0.col_ind, milu0.row_ptr, milu0.alu, milu0.jlu, milu0.ju, milu0.iw, ierr);
				}
				else {
					// недостаточно памяти на данном оборудовании.
					ipassage = 4;
					printf("Problem : not enough memory on your equipment...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}

				ipassage++;
				// Пока просто контролируем, в добавок хорошо бы выводить информацию на каком уровне это произошло.
#if doubleintprecision == 1
				printf("control memory allocation for ilu0 in RUMBA0.14 in level=%lld\n", ilevel);
#else
				printf("control memory allocation for ilu0 in RUMBA0.14 in level=%d\n", ilevel);
#endif

				system("PAUSE");
			} while ((ierr != 0) && (ipassage<4));

			if (ipassage == 4) {
				printf("Error memory alloc !!!\n");
				printf("failed to obtain an expansion for the 4 approaches...\n");
				printf("Please any key to exit...\n");
				//getchar();
				system("PAUSE");
				exit(1);
			}
		}

		// По идее здесь уже можно освободить ОЗУ от milu0.val, milu0.col_ind, milu0.row_ptr.
		// 7 ноября 2016.

	}

	if (!flag) {
		printf("Error equation 3D to CRS: zero diagonal element...\n");
		//getchar();
		system("pause");
	}

	for (integer k = 0; k<n; k++) if (milu0.col_ind[k] == (-1)) {
		printf("Error equation3D to CRS.\n");
		//getchar();
		system("pause");
	}

	// Нумерация начинается с единицы (это важно).
	for (integer k = 0; k<maxelm_plus_maxbound; k++) {
		if (milu0.val[milu0.row_ptr[k] - 1]<nonzeroEPS) {
#if doubleintprecision == 1
			printf("negativ diagonal element equation3DtoCRS %lld\n", k);
#else
			printf("negativ diagonal element equation3DtoCRS %d\n", k);
#endif

			//getchar();
			system("pause");
		}
	}


} // equation3DtoCRSRUMBA0

void nested_desection_patch(Ak2 &Amat, integer isize_na, bool* &nested_desection,
	integer* &row_ptr_start, integer* &row_ptr_end, integer iadd_now) {

	// Данный метод самописный и очень хреновый. Его не рекомендуется использовать в частности из-за граничных узлов.
	// Лучше производить разбиения по геометрическому признаку. Т.е. если есть координаты узлов, то слева выделить всё с координатами меньше
	// чем центральная координата, а справа всё что больше. Помоему так делают в Open Foam X.

	// Описание способа работы.
	// Метод требует подачи на вход матрицы в Ak1 формате.
	// Данная функция выдаёт булевый вектор nested desection [1..isize_na]:
	// Ячейки равные false относятся к границе сепаратора (разделителя двух подобластей).
	// Введем число integer middle = (integer)(0.5*(1 + isize_na));
	// Ячейки помеченные true индекс которых [1..middle] Могут выполнятся первым потоком.
	// Ячейки помеченные true индекс которых [middle+1..isize_na] Могут выполнятся вторым потоком независимо от первого.

	// Готовим nested desection
	// для двух потоков.
	// Самая подробная матрица 0.
	// nested_desection0
	bool *nested_desection1 = nullptr;
	nested_desection1 = (bool*)malloc((isize_na + 1) * sizeof(bool));
	char c1[19] = "nested_desection1";
	char c2[24] = "nested_desection_patch";
	char c3[19] = "nested_desection2";
	handle_error(nested_desection1, c1, c2, (isize_na + 1));
	bool *nested_desection2 = nullptr;
	nested_desection2 = (bool*)malloc((isize_na + 1) * sizeof(bool));
	handle_error(nested_desection2, c3, c2, (isize_na + 1));

	// инициализация.
	for (integer i_73 = 1; i_73 <= isize_na; i_73++) {
		nested_desection1[i_73] = false;
		nested_desection2[i_73] = false;
		nested_desection[i_73] = true;
	}
	integer middle = (1 + isize_na)/2;
	for (integer i_72 = 1; i_72 <= middle; i_72++) {
		for (integer i_82 = row_ptr_start[i_72 + iadd_now]; i_82 <= row_ptr_end[i_72 + iadd_now]; i_82++) {
			nested_desection1[Amat.j[i_82]] = true;
		}
	}
	for (integer i_72 = middle + 1; i_72 <= isize_na; i_72++) {
		for (integer i_82 = row_ptr_start[i_72 + iadd_now]; i_82 <= row_ptr_end[i_72 + iadd_now]; i_82++) {
			nested_desection2[Amat.j[i_82]] = true;
		}
	}
	// На пересечении двух множеств находятся узлы сепаратора.
	for (integer i_73 = 1; i_73 <= isize_na; i_73++) {
		if ((nested_desection1[i_73]) && (nested_desection2[i_73])) {
			nested_desection[i_73] = false;
		}
	}
	// запоминаем узлов из сепаратора.
	for (integer i_73 = 1; i_73 <= isize_na; i_73++) {
		nested_desection1[i_73] = nested_desection[i_73];
	}
	// Если какой-либо сосед принадлежит сепаратору то и сам узел принадлежит сепаратору.
	for (integer i_72 = 1; i_72 <= isize_na; i_72++) {
		for (integer i_82 = row_ptr_start[i_72 + iadd_now]; i_82 <= row_ptr_end[i_72 + iadd_now]; i_82++) {
			if (!nested_desection1[Amat.j[i_82]]) nested_desection[Amat.i[i_82]] = false;
		}
	}

	free(nested_desection1);
	free(nested_desection2);
} // nested_desection_patch

  // residual.
  // 3 jan 2016 ещё быстрее.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT1, typename doublerealT2>
void residualq2(Ak2 &Amat, integer istartq, integer iendq, doublerealT1* &x, doublerealT1* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, doublerealT2* &residual, doublerealT2* &my_diag)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	doublerealT2 dsum = 0.0;

#pragma omp parallel for 
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		residual[istr] = b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++) {
			if (1 || (Amat.aij[ii1] < 0.0)) {
				residual[istr] += -Amat.aij[ii1] *x[Amat.j[ii1]];
			}
			else {
				// 16 августа 2016.
				// не работает.
				//dsum += Amat.aij[ii1];
			}
		}
		residual[istr] -= (1.0 / Amat.aij[row_ptr_start[ii]])*x[istr];  // 1 april 2017
		//residual[istr] -= (my_diag[Amat.i[row_ptr_start[ii]]] + dsum)* x[istr]; // 3 jan 2016
	}


} // residualq2

 // residual.
  // 3 jan 2016 ещё быстрее.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT1, typename doublerealT2>
void residualq2(Ak2& Amat, integer istartq, integer iendq, doublerealT1*& x, doublerealT1*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublerealT2*& residual, doublerealT2*& my_diag, doublerealT2*& diag_minus_one)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	doublerealT2 dsum = 0.0;

#pragma omp parallel for 
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		// домножаем на минус 1.0 если требуется.
		residual[istr] = diag_minus_one[istr]*b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++) {
			if (1 || (Amat.aij[ii1] < 0.0)) {
				residual[istr] += -Amat.aij[ii1] * x[Amat.j[ii1]];
			}
			else {
				// 16 августа 2016.
				// не работает.
				//dsum += Amat.aij[ii1];
			}
		}
		residual[istr] -= (1.0 / Amat.aij[row_ptr_start[ii]]) * x[istr];  // 1 april 2017
		//residual[istr] -= (my_diag[Amat.i[row_ptr_start[ii]]] + dsum)* x[istr]; // 3 jan 2016
	}


} // residualq2

  // residual.
  // Анализ структуры невязки в случае проблемы с солвером.
  // 3 jan 2016 ещё быстрее.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT1, typename doublerealT2>
void residualq2_analysys(Ak2 &Amat, integer istartq, integer iendq, doublerealT1* &x, doublerealT1* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, doublerealT2* &residual, doublerealT2* &my_diag)
{
	// Анализ структуры невязки.
	// Мы делим вектор невязки на 20 равных частей и смотрим где невязка наиболее велика.



	printf("residual2 analysys start:\n");
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	doublerealT2 dsum = 0.0;
	doublerealT2 statistics_sum = 0.0;
	integer icount_log = 1;
	integer i_srez = ((endpos - startpos) / 20);
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		residual[istr] = b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++) {
			if (1 || (Amat.aij[ii1] < 0.0)) {
				residual[istr] += -Amat.aij[ii1] *x[Amat.j[ii1]];
			}
			else {
				// 16 августа 2016.
				// не работает.
				//dsum += Amat.aij[ii1];
			}
		}
		//residual[istr] += (-1.0 / Amat.aij[row_ptr_start[ii]])*x[istr];
		//residual[istr] -= (my_diag[Amat.i[row_ptr_start[ii]]] + dsum)* x[istr]; // 3 jan 2016
		residual[istr] -= (my_diag[istr] + dsum)* x[istr]; // 8 feb 2019
		statistics_sum += fabs(residual[istr]);
		if ((ii - startpos + 1) % i_srez == 0) {
#if doubleintprecision == 1
			printf("%lld %e ", icount_log++, statistics_sum);
#else
			printf("%d %e ", icount_log++, statistics_sum);
#endif

			statistics_sum = 0.0;
		}
	}
	printf("\n");

} // residualq2_analysys


// residual.
  // Анализ структуры невязки в случае проблемы с солвером.
  // 3 jan 2016 ещё быстрее.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT1, typename doublerealT2>
void residualq2_analysys(Ak2& Amat, integer istartq, integer iendq, doublerealT1*& x, doublerealT1*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublerealT2*& residual, doublerealT2*& my_diag, doublerealT2*& diag_minus_one)
{
	// Анализ структуры невязки.
	// Мы делим вектор невязки на 20 равных частей и смотрим где невязка наиболее велика.



	printf("residual2 analysys start:\n");
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	doublerealT2 dsum = 0.0;
	doublerealT2 statistics_sum = 0.0;
	integer icount_log = 1;
	integer i_srez = ((endpos - startpos) / 20);
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		residual[istr] = diag_minus_one[istr]*b[istr];
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++) {
			if (1 || (Amat.aij[ii1] < 0.0)) {
				residual[istr] += -Amat.aij[ii1] * x[Amat.j[ii1]];
			}
			else {
				// 16 августа 2016.
				// не работает.
				//dsum += Amat.aij[ii1];
			}
		}
		//residual[istr] += (-1.0 / Amat.aij[row_ptr_start[ii]])*x[istr];
		//residual[istr] -= (my_diag[Amat.i[row_ptr_start[ii]]] + dsum)* x[istr]; // 3 jan 2016
		residual[istr] -= (my_diag[istr] + dsum) * x[istr]; // 8 feb 2019
		statistics_sum += fabs(residual[istr]);
		if ((ii - startpos + 1) % i_srez == 0) {
#if doubleintprecision == 1
			printf("%lld %e ", icount_log++, statistics_sum);
#else
			printf("%d %e ", icount_log++, statistics_sum);
#endif

			statistics_sum = 0.0;
		}
	}
	printf("\n");

} // residualq2_analysys

  // smoother.
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

							 // За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega); //28.07.2016
							  //omega = 0.7;

							  //if (isorintmemo == iadd) {
							  // Это точно не первый раз
							  //bfirst = false;
							  //}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	if (omega < 1.0) {
		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy

		if (x_jacoby_buffer == nullptr) {
			printf("ERROR : x_jacoby_buffer == nullptr.\n");
			system("PAUSE");
			exit(1);
		}

		//#pragma loop(hint_parallel(8))
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			x_jacoby_buffer[istr] = x[istr];
		}


		//#pragma loop(hint_parallel(8))
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			doublerealT rold = x_jacoby_buffer[istr];

			// 13.07.2016
			doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// Распараллеливание почемуто тормозит очень сильно.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				// 13.07.2016
				// игнорирование positive connections.
				//if ((Amat.aij[ii1] < 0.0)) {
				rsum += -Amat.aij[ii1] *x_jacoby_buffer[ipos];
				//}
				//else {
				// не рабтает.
				//	ap_now += Amat.aij[ii1];
				//}
			}
			x[istr] += rsum;
			//x[istr] *= Amat.aij[row_ptr_start[ii]];
			// 13.07.2016
			x[istr] /= ap_now;

			// Запускается только после первого раза сразу запускать нельзя
			// т.к. начальное значение не является адекватным.
			x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
		}
	}
	else {

		// 3 ноября 2016. Симметричный метод Гаусса-Зейделя.
		if (isimmetricGS_switch == 0) {
			// 3 ноября 2016 было экспериментально доказано на BSKDmitrii что симметричный метод Гаусса - Зейделя РАСХОДИТСЯ.
			// НЕЛЬЗЯ ИСПОЛЬЗОВАТЬ СИММЕТРИЧНЫЙ МЕТОД ГАУССА -ЗЕЙДЕЛЯ. ИСПОЛЬЗУЙТЕ ОБЫЧНЫЙ МЕТОД ЗЕЙДЕЛЯ.


#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
				// 28.01.2017
				//ap_now = 0.0;

				x[istr] = b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ii1++)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					if (1) {
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						// Здесь всё в каше и позитив и негатив коннектионшс.
						rsum += -Amat.aij[ii1] *x[ipos];
						//}
						//else {
						// не работает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					else {
						// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

						// игнорирование positive connections.
						if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] *x[ipos];
							//ap_now += fabs(Amat.aij[ii1]);
						}
						else {
							// не работает.
							// Вероятно из-за того что так нарушен баланс.
							//ap_now += fabs(Amat.aij[ii1]);
							if (fabs(x[ipos]) > fabs(x[istr])) {
								ap_now += fabs(Amat.aij[ii1]);
							}
							else rsum += -Amat.aij[ii1] *x[ipos];
						}
					}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;

				// Запускается только после первого раза сразу запускать нельзя
				// т.к. начальное значение не является адекватным.
				x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
			}
			// Ни в коем случае не переключать. 3 ноября 2016.
			//isimmetricGS_switch = 1;

		}
		else {

#pragma omp parallel for
			for (integer ii = endpos; ii >= startpos; ii--) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

				x[istr] = b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ii1++)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					// игнорирование positive connections.
					//if ((Amat.aij[ii1] < 0.0)) {
					rsum += -Amat.aij[ii1] *x[ipos];
					//}
					//else {
					// не рабтает.
					//	ap_now += Amat.aij[ii1];
					//}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;

				// Запускается только после первого раза сразу запускать нельзя
				// т.к. начальное значение не является адекватным.
				x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
			}
			//isimmetricGS_switch = 0;
		}
	}


} // seidelqsor2

// smoother.
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2(Ak2& Amat, integer istartq, integer iendq,
	doublerealT*& x, doublerealT*& b,
	integer*& row_ptr_start, integer*& row_ptr_end,
	integer iadd, doublereal* &diag_minus_one)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

							 // За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega); //28.07.2016
							  //omega = 0.7;

							  //if (isorintmemo == iadd) {
							  // Это точно не первый раз
							  //bfirst = false;
							  //}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	if (omega < 1.0) {
		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy

		if (x_jacoby_buffer == nullptr) {
			printf("ERROR : x_jacoby_buffer == nullptr.\n");
			system("PAUSE");
			exit(1);
		}

		//#pragma loop(hint_parallel(8))
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			x_jacoby_buffer[istr] = x[istr];
		}


		//#pragma loop(hint_parallel(8))
#pragma omp parallel for
		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			doublerealT rold = x_jacoby_buffer[istr];

			// 13.07.2016
			doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

			x[istr] = diag_minus_one[istr]*b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// Распараллеливание почемуто тормозит очень сильно.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				// 13.07.2016
				// игнорирование positive connections.
				//if ((Amat.aij[ii1] < 0.0)) {
				rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
				//}
				//else {
				// не рабтает.
				//	ap_now += Amat.aij[ii1];
				//}
			}
			x[istr] += rsum;
			//x[istr] *= Amat.aij[row_ptr_start[ii]];
			// 13.07.2016
			x[istr] /= ap_now;

			// Запускается только после первого раза сразу запускать нельзя
			// т.к. начальное значение не является адекватным.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}
	}
	else {

		// 3 ноября 2016. Симметричный метод Гаусса-Зейделя.
		if (isimmetricGS_switch == 0) {
			// 3 ноября 2016 было экспериментально доказано на BSKDmitrii что симметричный метод Гаусса - Зейделя РАСХОДИТСЯ.
			// НЕЛЬЗЯ ИСПОЛЬЗОВАТЬ СИММЕТРИЧНЫЙ МЕТОД ГАУССА -ЗЕЙДЕЛЯ. ИСПОЛЬЗУЙТЕ ОБЫЧНЫЙ МЕТОД ЗЕЙДЕЛЯ.


#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
				// 28.01.2017
				//ap_now = 0.0;

				x[istr] = diag_minus_one[istr]*b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ii1++)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					if (1) {
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						// Здесь всё в каше и позитив и негатив коннектионшс.
						rsum += -Amat.aij[ii1] * x[ipos];
						//}
						//else {
						// не работает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					else {
						// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

						// игнорирование positive connections.
						if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//ap_now += fabs(Amat.aij[ii1]);
						}
						else {
							// не работает.
							// Вероятно из-за того что так нарушен баланс.
							//ap_now += fabs(Amat.aij[ii1]);
							if (fabs(x[ipos]) > fabs(x[istr])) {
								ap_now += fabs(Amat.aij[ii1]);
							}
							else rsum += -Amat.aij[ii1] * x[ipos];
						}
					}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;

				// Запускается только после первого раза сразу запускать нельзя
				// т.к. начальное значение не является адекватным.
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
			// Ни в коем случае не переключать. 3 ноября 2016.
			//isimmetricGS_switch = 1;

		}
		else {

#pragma omp parallel for
			for (integer ii = endpos; ii >= startpos; ii--) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

				x[istr] = diag_minus_one[istr]*b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ii1++)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					// игнорирование positive connections.
					//if ((Amat.aij[ii1] < 0.0)) {
					rsum += -Amat.aij[ii1] * x[ipos];
					//}
					//else {
					// не рабтает.
					//	ap_now += Amat.aij[ii1];
					//}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;

				// Запускается только после первого раза сразу запускать нельзя
				// т.к. начальное значение не является адекватным.
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
			//isimmetricGS_switch = 0;
		}
	}


} // seidelqsor2


  // smoother.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{

	//seidelqstable<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

	// SOR!!!
	// 3 jan 2016
	//seidelqsor(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// лучший выбор
	seidelqsor2(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// 14 января 2016 смена направлений сканирования.
	// оказалось значительно хуже.
	//seidelqsor3(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	//seidelqsor2Pcpu(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// early_naive_relaxation_method2.
    // Этот код (данная функция) больше не используется. Это одна из ранних
    // наивных попыток написать метод релаксации.
	//early_naive_relaxation_method2(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// Экспериментальный вариант слишком медленный:
	//bool* flag = new bool[iendq - istartq + 2];
	//early_naive_relaxation_method(Amat, istartq, iendq, x, b, flag, iendq - istartq + 1);
	//delete[] flag;

} // seidelq

// smoother.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublereal* &diag_minus_one)
{

	//seidelqstable<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

	// SOR!!!
	// 3 jan 2016
	//seidelqsor(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// лучший выбор
	seidelqsor2(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
	// 14 января 2016 смена направлений сканирования.
	// оказалось значительно хуже.
	//seidelqsor3(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	//seidelqsor2Pcpu(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// early_naive_relaxation_method2.
	//early_naive_relaxation_method2(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	// Экспериментальный вариант слишком медленный:
	//bool* flag = new bool[iendq - istartq + 2];
	//early_naive_relaxation_method(Amat, istartq, iendq, x, b, flag, iendq - istartq + 1);
	//delete[] flag;

} // seidelq

  // smoother.
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
  // 9 июня 2016. 
// early_naive_relaxation_method2.
// Этот код (данная функция) больше не используется. Это одна из ранних
// наивных попыток написать метод релаксации.
template <typename doublerealT>
void early_naive_relaxation_method2(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

							 // За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);
	omega = 1.0;

	// это нужне сделать единожды на фазе препроцессинга.
	doublerealT* diagonal = nullptr;
	diagonal = new doublerealT[iendq - istartq + 2];
	if (diagonal == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for diagonal in classical_ilu2...\n");
		printf("Please any key to exit...\n");
		//getchar();
		system("pause");
		exit(1);
	}
	// initialization:
	for (integer i_6 = 0; i_6 < iendq - istartq + 2; i_6++) diagonal[i_6] = 1.0;
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		{
			integer istr = ii - iadd;
			integer ic = ii;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			doublerealT ap = 0.0;
			doublerealT mult = 0.0;
			for (ic = is1; ic <= is2; ic++) {
				if (Amat.j[ic] != istr) {
					if (Amat.j[ic] < istr) {
						// aij*aji/djj;
						// найти aji;
						doublerealT aji = 0.0;
						// линейный поиск очень медленный особенно на каждой итерации.
						//for (integer ii1 = startpos; ii1 <= endpos; ii1++) {
						integer ii1 = iadd + Amat.j[ic];
						if (Amat.i[ii1] == Amat.j[ic]) {
							integer ic1 = ii1;
							integer is11 = row_ptr_start[ii1] + 1;
							integer is21 = row_ptr_end[ii1];
							for (ic1 = is11; ic1 <= is21; ic1++) {
								if (Amat.j[ic1] != Amat.j[ic]) {
									if (Amat.j[ic1] == istr) {
										aji = Amat.aij[ic1];
									}
								}
							}
							break;
						}
						//}
						mult += (Amat.aij[ic] *aji) / diagonal[Amat.j[ic]];
					}
				}
			}
			ap = 1.0 / Amat.aij[is1 - 1];
			diagonal[istr] = ap - mult;
		}
	}

	//if (isorintmemo == iadd) {
	// Это точно не первый раз
	//bfirst = false;
	//}
	startpos = istartq + iadd;
	endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		doublerealT rold = x[istr];

		x[istr] = b[istr];

		doublerealT rsum = 0.0;
		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];
		// Распараллеливание почемуто тормозит очень сильно.
		//#pragma omp parallel for reduction(+:rsum)
		for (integer ii1 = is1; ii1 <= is2; ii1++)
		{
			//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
			integer ipos = Amat.j[ii1];
			rsum += -Amat.aij[ii1] *x[ipos];
		}
		x[istr] += rsum;
		//x[istr] *= Amat.aij[row_ptr_start[ii]];
		// incomplete lower upwer decomposition.
		x[istr] /= diagonal[istr];

		// Запускается только после первого раза сразу запускать нельзя
		// т.к. начальное значение не является адекватным.
		//x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
	}

	delete[] diagonal;

} // early_naive_relaxation_method2

  // smoother.
  // 14 января 2015 каждый раз осуществляется смена направления сканирования.
  // Работает только с ноой логикой : nFinestSweeps=2, nPreSweeps=0, nPostSweeps=2.
  // Смена направлений оказалась значительно хуже по быстродействию.
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor3(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

	// За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// Это точно не первый раз
	//bfirst = false;
	//}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	if (bswitch_direct_seidelqsor3) {
		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// Распараллеливание почемуто тормозит очень сильно.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] *x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// Запускается только после первого раза сразу запускать нельзя
			// т.к. начальное значение не является адекватным.
			x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
		}
	}
	else {
		// другое направление сканирования.

		for (integer ii = endpos; ii >= startpos; ii--) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// Распараллеливание почемуто тормозит очень сильно.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] *x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// Запускается только после первого раза сразу запускать нельзя
			// т.к. начальное значение не является адекватным.
			x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
		}
	}

	// смена направления сканирования на обратное.
	bswitch_direct_seidelqsor3 = !bswitch_direct_seidelqsor3;

} // seidelqsor3


  // smoother.
  // 16 января 2016 распараллеливание на центральном процессоре.
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2Pcpu(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

							 // За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// Это точно не первый раз
	//bfirst = false;
	//}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	integer inumcore_loc = 2;

	if (inumcore_loc == 1) {

		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// Распараллеливание почемуто тормозит очень сильно.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] *x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// Запускается только после первого раза сразу запускать нельзя
			// т.к. начальное значение не является адекватным.
			x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
		}

	}

	if (inumcore_loc == 2) {

		// Здесь отсутствует информация о nested_desection
		// поэтому куча времени тратится на барьерную синхронизацию.

		integer middle = (startpos + endpos)/2;
		integer iadd1 = iadd;
		integer iadd2 = iadd;
		integer middle1 = middle;

#pragma omp parallel sections
		{
#pragma omp section
			{
				// работа первого потока.
				for (integer ii = startpos; ii <= middle; ii++) {
					integer istr = ii - iadd1;
					doublerealT rold = x[istr];

					doublerealT x1buf = 0.0;
					x1buf = b[istr];

					doublerealT rsum1 = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x1buf += -Amat.aij[ii1]*x[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						rsum1 = rsum1 - Amat.aij[ii1] *x[ipos];
					}
					x1buf = x1buf + rsum1;
					x1buf = x1buf * Amat.aij[row_ptr_start[ii]];

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					x1buf = omega*x1buf + (1.0 - omega)*rold; // this is SOR
					x[istr] = x1buf;
				}
			}

#pragma omp section
			{
				// работа второго потока. 
				for (integer ii_1 = middle1 + 1; ii_1 <= endpos; ii_1++) {
					integer istr1 = ii_1 - iadd2;
					doublerealT rold1 = x[istr1];

					doublerealT x2buf = 0.0;
					x2buf = b[istr1];

					doublerealT rsum2 = 0.0;
					integer is3 = row_ptr_start[ii_1] + 1;
					integer is4 = row_ptr_end[ii_1];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii2 = is3; ii2 <= is4; ii2++)
					{
						//x[istr1] += -Amat.aij[ii2]*x[Amat.j[ii2]];
						integer ipos = Amat.j[ii2];
						rsum2 = rsum2 - Amat.aij[ii2] *x[ipos];
					}
					x2buf = x2buf + rsum2;
					x2buf = x2buf*Amat.aij[row_ptr_start[ii_1]];

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					x2buf = omega*x2buf + (1.0 - omega)*rold1; // this is SOR
					x[istr1] = x2buf;
				}
			}
		}

	}


} // seidelqsor2Pcpu



  // smoother.
  // 16 января 2016 распараллеливание на центральном процессоре.
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2Pcpu(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, bool* &bnested_desection, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

	// За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// Это точно не первый раз
	//bfirst = false;
	//}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	integer inumcore_loc = 1;

	if (inumcore_loc == 1) {

		// Однопоточный вариант программы.
		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// Распараллеливание почемуто тормозит очень сильно.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] *x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// Запускается только после первого раза сразу запускать нельзя
			// т.к. начальное значение не является адекватным.
			x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
		}

	}

	if (inumcore_loc == 2) {

		integer middle = (integer)(0.5*(startpos + endpos));
		doublerealT omega1 = omega;
		doublerealT omega2 = omega;
		integer iadd1 = iadd;
		integer iadd2 = iadd;
		integer middle1 = middle;
		doublerealT ic_all = 0.0;
		doublerealT ic_separator = 0.0;
		for (integer i_51 = startpos - iadd1; i_51 <= endpos - iadd1; i_51++) {
			bnested_desection_global_amg[i_51] = bnested_desection[i_51];
			ic_all += 1.0;
			if (!bnested_desection_global_amg[i_51]) ic_separator += 1.0;
		}
		//printf("all=%e separator=%e\n",ic_all,ic_separator);
		//system("pause");

		if (2.0*ic_separator < ic_all) {
			// параллельно только при хорошем делении.



			//#pragma omp parallel sections num_threads(4)
			//{

			//printf_s("Hello from thread %d\n", omp_get_thread_num());
			//#pragma omp section
			//	printf_s("Hello from thread %d\n", omp_get_thread_num());
			///}

			// Версия с nesteddesection избавляет от барьерных синхронизаций и мы
			// получаем линейное ускорение производительности при увеличении числа ядер.

			//default(shared)

#pragma omp parallel  shared(bnested_desection_global_amg, bnested_desection,x,row_ptr_start,row_ptr_end,b,Amat)
			{
#pragma omp	sections
				{
#pragma omp section
					{
						//#pragma omp parallel sections num_threads(2)
						//	{

						// работа первого потока.
						for (integer ii = startpos; ii <= middle; ii++) {
							integer istr = ii - iadd1;
							if (bnested_desection_global_amg[istr]) {
								doublerealT rold = x[istr];

								doublerealT x1buf = 0.0;
								x1buf = b[istr];

								doublerealT rsum1 = 0.0;
								integer is1 = row_ptr_start[ii] + 1;
								integer is2 = row_ptr_end[ii];
								// Распараллеливание почемуто тормозит очень сильно.
								//#pragma omp parallel for reduction(+:rsum)
								for (integer ii1 = is1; ii1 <= is2; ii1++)
								{
									//x1buf += -Amat.aij[ii1]*x[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									rsum1 = rsum1 - Amat.aij[ii1] *x[ipos];
								}
								x1buf = x1buf + rsum1;
								x1buf = x1buf * Amat.aij[row_ptr_start[ii]];

								// Запускается только после первого раза сразу запускать нельзя
								// т.к. начальное значение не является адекватным.
								x1buf = omega1*x1buf + (1.0 - omega1)*rold; // this is SOR
								x[istr] = x1buf;
							}
						}
					}

#pragma omp section
					{
						// работа второго потока. 
						for (integer ii_1 = middle1 + 1; ii_1 <= endpos; ii_1++) {
							integer istr1 = ii_1 - iadd2;
							if (bnested_desection[istr1]) {
								doublerealT rold1 = x[istr1];

								doublerealT x2buf = 0.0;
								x2buf = b[istr1];

								doublerealT rsum2 = 0.0;
								integer is3 = row_ptr_start[ii_1] + 1;
								integer is4 = row_ptr_end[ii_1];
								// Распараллеливание почемуто тормозит очень сильно.
								//#pragma omp parallel for reduction(+:rsum)
								for (integer ii2 = is3; ii2 <= is4; ii2++)
								{
									//x[istr1] += -Amat.aij[ii2]*x[Amat.j[ii2]];
									integer ipos = Amat.j[ii2];
									rsum2 = rsum2 - Amat.aij[ii2] *x[ipos];
								}
								x2buf = x2buf + rsum2;
								x2buf = x2buf*Amat.aij[row_ptr_start[ii_1]];

								// Запускается только после первого раза сразу запускать нельзя
								// т.к. начальное значение не является адекватным.
								x2buf = omega2*x2buf + (1.0 - omega2)*rold1; // this is SOR
								x[istr1] = x2buf;
							}
						}
					}
				}
			}

			// Однопоточный смыкающий кусок.
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				if (!bnested_desection[istr]) {
					doublerealT rold = x[istr];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						rsum += -Amat.aij[ii1] *x[ipos];
					}
					x[istr] += rsum;
					x[istr] *= Amat.aij[row_ptr_start[ii]];

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
				}
			}
		}
		else {
			// серийно
			// Однопоточный вариант программы.
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				x[istr] = b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ii1++)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					rsum += -Amat.aij[ii1] *x[ipos];
				}
				x[istr] += rsum;
				x[istr] *= Amat.aij[row_ptr_start[ii]];

				// Запускается только после первого раза сразу запускать нельзя
				// т.к. начальное значение не является адекватным.
				x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
			}
		}

	}


} // seidelqsor2Pcpu+nested desection


// smoother.
  // 16 января 2016 распараллеливание на центральном процессоре.
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2Pcpu(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublereal* &diag_minus_one)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

							 // За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// Это точно не первый раз
	//bfirst = false;
	//}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	integer inumcore_loc = 2;

	if (inumcore_loc == 1) {

		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = diag_minus_one[istr]*b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// Распараллеливание почемуто тормозит очень сильно.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] * x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// Запускается только после первого раза сразу запускать нельзя
			// т.к. начальное значение не является адекватным.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}

	}

	if (inumcore_loc == 2) {

		// Здесь отсутствует информация о nested_desection
		// поэтому куча времени тратится на барьерную синхронизацию.

		integer middle = (startpos + endpos) / 2;
		integer iadd1 = iadd;
		integer iadd2 = iadd;
		integer middle1 = middle;

#pragma omp parallel sections
		{
#pragma omp section
			{
				// работа первого потока.
				for (integer ii = startpos; ii <= middle; ii++) {
					integer istr = ii - iadd1;
					doublerealT rold = x[istr];

					doublerealT x1buf = 0.0;
					x1buf = diag_minus_one[istr]*b[istr];

					doublerealT rsum1 = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x1buf += -Amat.aij[ii1]*x[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						rsum1 = rsum1 - Amat.aij[ii1] * x[ipos];
					}
					x1buf = x1buf + rsum1;
					x1buf = x1buf * Amat.aij[row_ptr_start[ii]];

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					x1buf = omega * x1buf + (1.0 - omega) * rold; // this is SOR
					x[istr] = x1buf;
				}
			}

#pragma omp section
			{
				// работа второго потока. 
				for (integer ii_1 = middle1 + 1; ii_1 <= endpos; ii_1++) {
					integer istr1 = ii_1 - iadd2;
					doublerealT rold1 = x[istr1];

					doublerealT x2buf = 0.0;
					x2buf = diag_minus_one[istr1]*b[istr1];

					doublerealT rsum2 = 0.0;
					integer is3 = row_ptr_start[ii_1] + 1;
					integer is4 = row_ptr_end[ii_1];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii2 = is3; ii2 <= is4; ii2++)
					{
						//x[istr1] += -Amat.aij[ii2]*x[Amat.j[ii2]];
						integer ipos = Amat.j[ii2];
						rsum2 = rsum2 - Amat.aij[ii2] * x[ipos];
					}
					x2buf = x2buf + rsum2;
					x2buf = x2buf * Amat.aij[row_ptr_start[ii_1]];

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					x2buf = omega * x2buf + (1.0 - omega) * rold1; // this is SOR
					x[istr1] = x2buf;
				}
			}
		}

	}


} // seidelqsor2Pcpu



  // smoother.
  // 16 января 2016 распараллеливание на центральном процессоре.
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelqsor2Pcpu(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, bool*& bnested_desection, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublereal*& diag_minus_one)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

	// За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega);

	//if (isorintmemo == iadd) {
	// Это точно не первый раз
	//bfirst = false;
	//}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	integer inumcore_loc = 1;

	if (inumcore_loc == 1) {

		// Однопоточный вариант программы.
		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			doublerealT rold = x[istr];

			x[istr] = diag_minus_one[istr] * b[istr];

			doublerealT rsum = 0.0;
			integer is1 = row_ptr_start[ii] + 1;
			integer is2 = row_ptr_end[ii];
			// Распараллеливание почемуто тормозит очень сильно.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
				integer ipos = Amat.j[ii1];
				rsum += -Amat.aij[ii1] * x[ipos];
			}
			x[istr] += rsum;
			x[istr] *= Amat.aij[row_ptr_start[ii]];

			// Запускается только после первого раза сразу запускать нельзя
			// т.к. начальное значение не является адекватным.
			x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
		}

	}

	if (inumcore_loc == 2) {

		integer middle = (integer)(0.5 * (startpos + endpos));
		doublerealT omega1 = omega;
		doublerealT omega2 = omega;
		integer iadd1 = iadd;
		integer iadd2 = iadd;
		integer middle1 = middle;
		doublerealT ic_all = 0.0;
		doublerealT ic_separator = 0.0;
		for (integer i_51 = startpos - iadd1; i_51 <= endpos - iadd1; i_51++) {
			bnested_desection_global_amg[i_51] = bnested_desection[i_51];
			ic_all += 1.0;
			if (!bnested_desection_global_amg[i_51]) ic_separator += 1.0;
		}
		//printf("all=%e separator=%e\n",ic_all,ic_separator);
		//system("pause");

		if (2.0 * ic_separator < ic_all) {
			// параллельно только при хорошем делении.



			//#pragma omp parallel sections num_threads(4)
			//{

			//printf_s("Hello from thread %d\n", omp_get_thread_num());
			//#pragma omp section
			//	printf_s("Hello from thread %d\n", omp_get_thread_num());
			///}

			// Версия с nesteddesection избавляет от барьерных синхронизаций и мы
			// получаем линейное ускорение производительности при увеличении числа ядер.

			//default(shared)

#pragma omp parallel  shared(bnested_desection_global_amg, bnested_desection,x,row_ptr_start,row_ptr_end,b,Amat)
			{
#pragma omp	sections
				{
#pragma omp section
					{
						//#pragma omp parallel sections num_threads(2)
						//	{

						// работа первого потока.
						for (integer ii = startpos; ii <= middle; ii++) {
							integer istr = ii - iadd1;
							if (bnested_desection_global_amg[istr]) {
								doublerealT rold = x[istr];

								doublerealT x1buf = 0.0;
								x1buf = diag_minus_one[istr] * b[istr];

								doublerealT rsum1 = 0.0;
								integer is1 = row_ptr_start[ii] + 1;
								integer is2 = row_ptr_end[ii];
								// Распараллеливание почемуто тормозит очень сильно.
								//#pragma omp parallel for reduction(+:rsum)
								for (integer ii1 = is1; ii1 <= is2; ii1++)
								{
									//x1buf += -Amat.aij[ii1]*x[Amat.j[ii1]];
									integer ipos = Amat.j[ii1];
									rsum1 = rsum1 - Amat.aij[ii1] * x[ipos];
								}
								x1buf = x1buf + rsum1;
								x1buf = x1buf * Amat.aij[row_ptr_start[ii]];

								// Запускается только после первого раза сразу запускать нельзя
								// т.к. начальное значение не является адекватным.
								x1buf = omega1 * x1buf + (1.0 - omega1) * rold; // this is SOR
								x[istr] = x1buf;
							}
						}
					}

#pragma omp section
					{
						// работа второго потока. 
						for (integer ii_1 = middle1 + 1; ii_1 <= endpos; ii_1++) {
							integer istr1 = ii_1 - iadd2;
							if (bnested_desection[istr1]) {
								doublerealT rold1 = x[istr1];

								doublerealT x2buf = 0.0;
								x2buf = diag_minus_one[istr1] * b[istr1];

								doublerealT rsum2 = 0.0;
								integer is3 = row_ptr_start[ii_1] + 1;
								integer is4 = row_ptr_end[ii_1];
								// Распараллеливание почемуто тормозит очень сильно.
								//#pragma omp parallel for reduction(+:rsum)
								for (integer ii2 = is3; ii2 <= is4; ii2++)
								{
									//x[istr1] += -Amat.aij[ii2]*x[Amat.j[ii2]];
									integer ipos = Amat.j[ii2];
									rsum2 = rsum2 - Amat.aij[ii2] * x[ipos];
								}
								x2buf = x2buf + rsum2;
								x2buf = x2buf * Amat.aij[row_ptr_start[ii_1]];

								// Запускается только после первого раза сразу запускать нельзя
								// т.к. начальное значение не является адекватным.
								x2buf = omega2 * x2buf + (1.0 - omega2) * rold1; // this is SOR
								x[istr1] = x2buf;
							}
						}
					}
				}
			}

			// Однопоточный смыкающий кусок.
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				if (!bnested_desection[istr]) {
					doublerealT rold = x[istr];

					x[istr] = diag_minus_one[istr] * b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						rsum += -Amat.aij[ii1] * x[ipos];
					}
					x[istr] += rsum;
					x[istr] *= Amat.aij[row_ptr_start[ii]];

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				}
			}
		}
		else {
			// серийно
			// Однопоточный вариант программы.
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				doublerealT rold = x[istr];

				x[istr] = diag_minus_one[istr] * b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ii1++)
				{
					//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					rsum += -Amat.aij[ii1] * x[ipos];
				}
				x[istr] += rsum;
				x[istr] *= Amat.aij[row_ptr_start[ii]];

				// Запускается только после первого раза сразу запускать нельзя
				// т.к. начальное значение не является адекватным.
				x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
			}
		}

	}


} // seidelqsor2Pcpu+nested desection

// smoother.
// 1 september 2015.
template <typename doublerealT>
void seidel1(Ak2 &Amat, integer istart, integer iend, doublerealT* &x, doublerealT* &b, bool* &flag, integer n)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[Amat.i[ii]] == false) {
			integer istr = Amat.i[ii];
			integer ic = ii;
			doublerealT ap = 0.0;
			x[istr] = b[istr];
			while ((ic <= iend) && (Amat.i[ic] == istr)) {
				if (Amat.j[ic] != istr) {
					x[istr] += -Amat.aij[ic] *x[Amat.j[ic]];
				}
				else ap = Amat.aij[ic];
				ic++;
			}
			/*
			if (fabs(ap) < 1.0e-30) {
			#if doubleintprecision == 1
				printf("zero diagonal elements in string %lld",istr);
			#else
				printf("zero diagonal elements in string %d",istr);
			#endif

			getchar();
			exit(1);
			}
			else */ {
				x[istr] /= ap;
			}
			flag[Amat.i[ii]] = true;
		}
	}


} // seidel1

// smoother.
// 5 jan 2016.
// 1 september 2015.
template <typename doublerealT>
void seidelsor(Ak2 &Amat, integer istart, integer iend, doublerealT* &x, doublerealT* &b, bool* &flag, integer n)
{

	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	doublerealT rn = 1.0;
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[Amat.i[ii]] == false) {
			flag[Amat.i[ii]] = true;
			rn += 1.0;
		}
	}

	doublerealT omega = 1.0; // initialize.

	// За подробностями смотри книгу Патрика Роуча стр. 183.
	optimal_omega(rn, omega);//28.07.2016
	//omega = 0.7;

	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[Amat.i[ii]] == false) {
			integer istr = Amat.i[ii];
			integer ic = ii;
			doublerealT ap = 0.0;
			doublerealT rold = x[istr];
			x[istr] = b[istr];
			while ((ic <= iend) && (Amat.i[ic] == istr)) {
				if (Amat.j[ic] != istr) {
					if (Amat.aij[ic] < 0.0) {
						x[istr] += -Amat.aij[ic] *x[Amat.j[ic]];
					}
					else {
						ap += Amat.aij[ic];
					}
				}
				else ap += Amat.aij[ic];
				ic++;
			}
			/*
			if (fabs(ap) < 1.0e-30) {
			#if doubleintprecision == 1
					printf("zero diagonal elements in string %lld",istr);
			#else
					printf("zero diagonal elements in string %d",istr);
			#endif

			getchar();
			exit(1);
			}
			else */ {
				x[istr] /= ap;
				x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
			}
			flag[Amat.i[ii]] = true;
		}
	}


} // seidelsor


// smoother.
// 5 jan 2016.
// 1 september 2015.
// 9 июня 2016.
// early_naive_relaxation_method.
// Этот код (данная функция) больше не используется. Это одна из ранних
// наивных попыток написать метод релаксации.
template <typename doublerealT>
void early_naive_relaxation_method(Ak2 &Amat, integer istart, integer iend,
	doublerealT* &x, doublerealT* &b, bool* &flag, integer n)
{

	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	doublerealT rn = 1.0;
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[Amat.i[ii]] == false) {
			flag[Amat.i[ii]] = true;
			rn += 1.0;
		}
	}

	doublerealT omega = 1.0; // initialize.
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}

	doublerealT* diagonal = new doublerealT[(integer)(rn)];
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[Amat.i[ii]] == false) {
			integer istr = Amat.i[ii];
			integer ic = ii;
			doublerealT ap = 0.0;
			doublerealT mult = 0.0;
			x[istr] = b[istr];
			while ((ic <= iend) && (Amat.i[ic] == istr)) {
				if (Amat.j[ic] != istr) {
					if (Amat.j[ic] < istr) {
						// aij*aji/djj;
						// найти aji;
						doublerealT aji = 0.0;
						// Здесь нужно сделать хеш таблицу которая позволяет избавиться от линейного поиска.
						// линейный поиск очень медленный особенно на каждой итерации.
						for (integer ii1 = istart; ii1 <= iend; ii1++) {
							if (Amat.i[ii1] == Amat.j[ic]) {
								integer ic1 = ii1;
								while ((ic1 <= iend) && (Amat.i[ic1] == Amat.j[ic])) {
									if (Amat.j[ic1] != Amat.j[ic]) {
										if (Amat.j[ic1] == istr) {
											aji = Amat.aij[ic1];
										}
									}
								}
								break;
							}
						}
						mult += (Amat.aij[ic] *aji) / diagonal[Amat.j[ic]];
					}
				}
				else ap = Amat.aij[ic];
				ic++;
			}
			diagonal[istr] = ap - mult;
			flag[Amat.i[ii]] = true;
		}
	}


	// За подробностями смотри книгу Патрика Роуча стр. 183.
	optimal_omega(rn, omega);
	omega = 1.0;

	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	for (integer i = 1; i <= n; i++) {
		flag[i] = false;
	}
	for (integer ii = istart; ii <= iend; ii++) {
		if (flag[Amat.i[ii]] == false) {
			integer istr = Amat.i[ii];
			integer ic = ii;
			doublerealT ap = 0.0;
			doublerealT rold = x[istr];
			x[istr] = b[istr];
			while ((ic <= iend) && (Amat.i[ic] == istr)) {
				if (Amat.j[ic] != istr) {
					x[istr] += -Amat.aij[ic] *x[Amat.j[ic]];
				}
				else ap = Amat.aij[ic];
				ic++;
			}
			/*
			if (fabs(ap) < 1.0e-30) {
			#if doubleintprecision == 1
				printf("zero diagonal elements in string %lld",istr);
			#else
				printf("zero diagonal elements in string %d",istr);
			#endif

			getchar();
			exit(1);
			}
			else */ {
			//x[istr] /= ap;
				x[istr] /= diagonal[istr];
				//x[istr] = omega*x[istr] + (1.0 - omega)*rold; // this is SOR
			}
			flag[Amat.i[ii]] = true;
		}
	}

	delete[] diagonal;

} // early_naive_relaxation_method


// smoother.
// 5 jan 2016.
// 1 september 2015.
template <typename doublerealT>
void seidel(Ak2 &Amat, integer istart, integer iend, doublerealT* &x, doublerealT* &b, bool* &flag, integer n)
{
	//seidel1<doublerealT>(Amat,  istart,  iend, x, b, flag, n);
	seidelsor<doublerealT>(Amat, istart, iend, x, b, flag, n);
	//early_naive_relaxation_method<doublerealT>(Amat, istart, iend, x, b, flag, n);
}

// smoother.
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelqstable(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		x[istr] = b[istr];

		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++)
		{
			x[istr] += -Amat.aij[ii1] *x[Amat.j[ii1]];
		}
		x[istr] *= Amat.aij[row_ptr_start[ii]];
	}


} // seidelq


template <typename doublerealT>
void seidelqsor(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.15; // лучший выбор.
	bool bfirst = false;
	//if (isorintmemo == iadd) {
	// Это точно не первый раз
	//bfirst = false;
	//}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		doublerealT rold = x[istr];

		x[istr] = b[istr];

		doublerealT rsum = 0.0;
		// Параллельность почему-то тормозит очень сильно.
		//#pragma omp parallel for reduction(+:rsum)
		for (integer ii1 = row_ptr_start[ii] + 1; ii1 <= row_ptr_end[ii]; ii1++)
		{
			rsum += -Amat.aij[ii1] *x[Amat.j[ii1]];
			//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
		}
		x[istr] += rsum;

		x[istr] *= Amat.aij[row_ptr_start[ii]];
		if (bfirst) {
			bfirst = false;
		}
		else {
			// Запускается только после первого раза сразу запускать нельзя
			// т.к. начальное значение не является адекватным.
			x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
		}
	}


} // seidelqsor


 // smoother.
  // 6 июня 2017 Трёхшаговый метод Рунге-Кутты, параметры взяты из литературы.
// Достоинство в том что оптимальные параметры известны из литературы и не надо
// ничего вручную подбирать.
// Методы Рунге-Кутты в качестве сглаживателей рекомендованы для очень плохообусловленных задач,
// что вызвано отчасти плохими сетками (например, АЛИС сетками).
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void Runge_Kutt_3or5(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, integer iorder)
{
	// iorder == 3 or 5. Трёхшаговый или пятишаговый методы Рунге - Кутты.

	if ((iorder == 3) || (iorder == 5)) {

		// Трёхшаговый метод Рунге - Кутты.

		// Методы ускорения газодинамических расчётов на неструктурированных сетках. К.Н.Волков, под редакцией проф. В.Н.Емельянова
		// Москва ФИЗМАТЛИТ 2014.
		doublerealT m[5];
		if (iorder == 3) {
			if (1) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.2075;
				m[1] = 0.5915;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
			else {
				// Направленное огрубление.
				m[0] = 0.2239;
				m[1] = 0.5653;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
		}

		if (iorder == 5) {
			if (0) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.0962;
				m[1] = 0.2073;
				m[2] = 0.3549;
				m[3] = 0.6223;
				m[4] = 1.0;
			}
			else {
				// Направленное огрубление.
				m[0] = 0.0870;
				m[1] = 0.1892;
				m[2] = 0.3263;
				m[3] = 0.5558;
				m[4] = 1.0;
			}
		}

		// istart - начальная позиция ненулевых элементов в матрице А.
		// iend - конечная позиция ненулевых элементов в матрице А.


		doublerealT rn = (doublerealT)(iendq - istartq + 1);

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;


		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy


		if (x_jacoby_buffer == nullptr) {
			printf("ERROR : x_jacoby_buffer == nullptr.\n");
			system("PAUSE");
			exit(1);
		}



		for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; inumber_step_Runge_Kutt++) {

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {

				integer istr = ii - iadd;
				doublerealT rold = x_jacoby_buffer[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

				x[istr] = b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ii1++)
				{
					//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					// игнорирование positive connections.
					//if ((Amat.aij[ii1] < 0.0)) {
					rsum += -Amat.aij[ii1] *x_jacoby_buffer[ipos];

					//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
					//}
					//else {
					// не рабтает.
					//	ap_now += Amat.aij[ii1];
					//}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;


				x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
				//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

				//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
				//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
			}

		}

	}
	else {
		// Неправильный порядок метода Рунге - Кутты. 
		// Предусмотрены только третий и пятый порядки.
		printf("order Runge Kutt method is bad...\n ");
		system("pause");
		exit(1);
	}


} // Runge_Kutt_3or5

  // smoother.
  // 6 июня 2017 Трёхшаговый метод Рунге-Кутты, параметры взяты из литературы.
  // Достоинство в том что оптимальные параметры известны из литературы и не надо
  // ничего вручную подбирать.
  // Методы Рунге-Кутты в качестве сглаживателей рекомендованы для очень плохообусловленных задач,
  // что вызвано отчасти плохими сетками (например, АЛИС сетками).
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void Runge_Kutt_3or5(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start,
	integer * &row_ptr_end, integer iadd, integer iorder, bool* &F_false_C_true, integer idirect)
{
	// iorder == 3 or 5. Трёхшаговый или пятишаговый методы Рунге - Кутты.

	if ((iorder == 3) || (iorder == 5)) {

		// Трёхшаговый метод Рунге - Кутты.

		// Методы ускорения газодинамических расчётов на неструктурированных сетках. К.Н.Волков, под редакцией проф. В.Н.Емельянова
		// Москва ФИЗМАТЛИТ 2014.
		doublerealT m[5];
		if (iorder == 3) {
			if (1) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.2075;
				m[1] = 0.5915;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
			else {
				// Направленное огрубление.
				m[0] = 0.2239;
				m[1] = 0.5653;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
		}

		if (iorder == 5) {
			if (0) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.0962;
				m[1] = 0.2073;
				m[2] = 0.3549;
				m[3] = 0.6223;
				m[4] = 1.0;
			}
			else {
				// Направленное огрубление.
				m[0] = 0.0870;
				m[1] = 0.1892;
				m[2] = 0.3263;
				m[3] = 0.5558;
				m[4] = 1.0;
			}
		}

		// istart - начальная позиция ненулевых элементов в матрице А.
		// iend - конечная позиция ненулевых элементов в матрице А.


		doublerealT rn = (doublerealT)(iendq - istartq + 1);

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;


		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy





		for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; inumber_step_Runge_Kutt++) {

			if (idirect == 1) {

				// Восходящая ветвь : сначала F потом C.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == false) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] *x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
					}
				}


				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == true) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] *x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

					}
				}

			}
			else {
				// Нисходящая ветвь : сначала C потом F.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == true) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] *x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
					}
				}


				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == false) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] *x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

					}
				}


			}

		}

	}
	else {
		// Неправильный порядок метода Рунге - Кутты. 
		// Предусмотрены только третий и пятый порядки.
		printf("order Runge Kutt method is bad...\n ");
		system("pause");
		exit(1);
	}


} // Runge_Kutt_3or5

// Sparse approximate inverse relaxation scheme.
// n – число неизвестных,
// x – вектор решения (текущее приближение к решению)
// val, col_ind, row_ptr матрица в CRS формате.
// rthdsd – правая часть СЛАУ.
template <typename doublerealT>
void spai0_smoother(integer n, doublerealT*& val,
	integer*& col_ind, integer*& row_ptr,
	doublerealT*& x, doublerealT*& rthdsd)
{

#pragma omp parallel for 
	for (integer i = 0; i < n; i++) {
		doublereal num = 0.0;
		doublereal den = 0.0;
		doublereal sum = 0.0;

		for (integer j = row_ptr[i]; j <= row_ptr[i + 1] - 1; j++) {

			den += val[j] * val[j];

			if (col_ind[j] == i) {
				num += val[j]; // диагональ
			}
			else {
				sum += val[j] * x[col_ind[j]];
			}
		}
		doublereal M_precond =  num / den;
		x[i] = M_precond * (rthdsd[i] - sum) / num;
	}
} // spai0_smoother

// 21.12.2019
// Sparse approximate inverse relaxation scheme.
// n – число неизвестных,
// x – вектор решения (текущее приближение к решению)
// val, col_ind, row_ptr матрица в CRS формате.
// rthdsd – правая часть СЛАУ.
template <typename doublerealT>
void spai0_smoother(Ak2 &Amat, integer istartq, integer iendq,
	doublerealT* &x, doublerealT* &b, integer* &row_ptr_start, integer* &row_ptr_end, integer iadd)
{

	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

//#pragma omp parallel for 
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
	
		doublereal num = 0.0;
		doublereal den = 0.0;
		doublereal sum = 0.0;
		doublereal E1 = 1.0 / Amat.aij[row_ptr_start[ii]];
		E1 *= E1;// Для нормировки иначе происходит неправильная обработка.

		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];

		num =  1.0 / Amat.aij[row_ptr_start[ii]];
		den += 1.0;// num* num;

		// Распараллеливание почемуто тормозит очень сильно.
		//#pragma omp parallel for reduction(+:rsum)
		for (integer ii1 = is1; ii1 <= is2; ii1++)
		{

			den += (Amat.aij[ii1] * Amat.aij[ii1])/ E1;

			sum += Amat.aij[ii1] * x[Amat.j[ii1]];
			
		}
		doublereal M_precond = 1.0 / den; //num / den;
		x[istr] = M_precond * (b[istr] - sum) / num;
		/*// debug
		if (x[istr] != x[istr]) {
			printf("%e 1\n", x[istr]);
			system("pause");
		}
		*/
	}
} // spai0_smoother

// 21.12.2019
// Sparse approximate inverse relaxation scheme.
// n – число неизвестных,
// x – вектор решения (текущее приближение к решению)
// val, col_ind, row_ptr матрица в CRS формате.
// rthdsd – правая часть СЛАУ.
template <typename doublerealT>
void spai0_smoother(Ak2 &Amat, integer istartq, integer iendq,
	doublerealT* &x, doublerealT* &b, integer* &row_ptr_start,
	integer* &row_ptr_end, integer iadd, doublereal* &diag_minus_one)
{

	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	//#pragma omp parallel for 
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;

		doublereal num = 0.0;
		doublereal den = 0.0;
		doublereal sum = 0.0;
		doublereal E1 = 1.0 / Amat.aij[row_ptr_start[ii]];
		E1 *= E1;// Для нормировки иначе происходит неправильная обработка.

		integer is1 = row_ptr_start[ii] + 1;
		integer is2 = row_ptr_end[ii];

		num=1.0 / Amat.aij[row_ptr_start[ii]];
		//num = 1.0;
		den += 1.0;// num* num;

		// Распараллеливание почемуто тормозит очень сильно.
		//#pragma omp parallel for reduction(+:rsum)
		for (integer ii1 = is1; ii1 <= is2; ii1++)
		{

			integer ipos = Amat.j[ii1];

			den += (Amat.aij[ii1] * Amat.aij[ii1])/ E1;

			sum += Amat.aij[ii1] * x[ipos];
			
		}
		doublereal M_precond = 1.0 / den;//num/den;
		x[istr] = M_precond * (diag_minus_one[istr] * b[istr] - sum) / num;
		/* // debug 
		if (x[istr] != x[istr]) {
			printf("%e M=%e num=%e den=%e b=%e sum=%e 2\n", x[istr], M_precond,num,den, b[istr],sum);
			system("pause");
		}
		*/
	}
} // spai0_smoother

// smoother.
  // 6 июня 2017 Трёхшаговый метод Рунге-Кутты, параметры взяты из литературы.
// Достоинство в том что оптимальные параметры известны из литературы и не надо
// ничего вручную подбирать.
// Методы Рунге-Кутты в качестве сглаживателей рекомендованы для очень плохообусловленных задач,
// что вызвано отчасти плохими сетками (например, АЛИС сетками).
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void Runge_Kutt_3or5(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, 
	integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, integer iorder, doublereal* &diag_minus_one)
{
	// iorder == 3 or 5. Трёхшаговый или пятишаговый методы Рунге - Кутты.

	if ((iorder == 3) || (iorder == 5)) {

		// Трёхшаговый метод Рунге - Кутты.

		// Методы ускорения газодинамических расчётов на неструктурированных сетках. К.Н.Волков, под редакцией проф. В.Н.Емельянова
		// Москва ФИЗМАТЛИТ 2014.
		doublerealT m[5];
		if (iorder == 3) {
			if (1) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.2075;
				m[1] = 0.5915;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
			else {
				// Направленное огрубление.
				m[0] = 0.2239;
				m[1] = 0.5653;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
		}

		if (iorder == 5) {
			if (0) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.0962;
				m[1] = 0.2073;
				m[2] = 0.3549;
				m[3] = 0.6223;
				m[4] = 1.0;
			}
			else {
				// Направленное огрубление.
				m[0] = 0.0870;
				m[1] = 0.1892;
				m[2] = 0.3263;
				m[3] = 0.5558;
				m[4] = 1.0;
			}
		}

		// istart - начальная позиция ненулевых элементов в матрице А.
		// iend - конечная позиция ненулевых элементов в матрице А.


		doublerealT rn = (doublerealT)(iendq - istartq + 1);

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;


		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy


		if (x_jacoby_buffer == nullptr) {
			printf("ERROR : x_jacoby_buffer == nullptr.\n");
			system("PAUSE");
			exit(1);
		}



		for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; inumber_step_Runge_Kutt++) {

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {

				integer istr = ii - iadd;
				doublerealT rold = x_jacoby_buffer[istr];

				// 13.07.2016
				doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

				x[istr] = diag_minus_one[istr]*b[istr];

				doublerealT rsum = 0.0;
				integer is1 = row_ptr_start[ii] + 1;
				integer is2 = row_ptr_end[ii];
				// Распараллеливание почемуто тормозит очень сильно.
				//#pragma omp parallel for reduction(+:rsum)
				for (integer ii1 = is1; ii1 <= is2; ii1++)
				{
					//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
					integer ipos = Amat.j[ii1];
					// 13.07.2016
					// игнорирование positive connections.
					//if ((Amat.aij[ii1] < 0.0)) {
					rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

					//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
					//}
					//else {
					// не рабтает.
					//	ap_now += Amat.aij[ii1];
					//}
				}
				x[istr] += rsum;
				//x[istr] *= Amat.aij[row_ptr_start[ii]];
				// 13.07.2016
				x[istr] /= ap_now;


				x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
				//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

				//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
				//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
			}

		}

	}
	else {
		// Неправильный порядок метода Рунге - Кутты. 
		// Предусмотрены только третий и пятый порядки.
		printf("order Runge Kutt method is bad...\n ");
		system("pause");
		exit(1);
	}


} // Runge_Kutt_3or5

  // smoother.
  // 6 июня 2017 Трёхшаговый метод Рунге-Кутты, параметры взяты из литературы.
  // Достоинство в том что оптимальные параметры известны из литературы и не надо
  // ничего вручную подбирать.
  // Методы Рунге-Кутты в качестве сглаживателей рекомендованы для очень плохообусловленных задач,
  // что вызвано отчасти плохими сетками (например, АЛИС сетками).
  // 5 января 2016 с использованием формулы из книги Патрика Роуча.
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void Runge_Kutt_3or5(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start,
	integer*& row_ptr_end, integer iadd, integer iorder, bool*& F_false_C_true, integer idirect, doublereal* &diag_minus_one)
{
	// iorder == 3 or 5. Трёхшаговый или пятишаговый методы Рунге - Кутты.

	if ((iorder == 3) || (iorder == 5)) {

		// Трёхшаговый метод Рунге - Кутты.

		// Методы ускорения газодинамических расчётов на неструктурированных сетках. К.Н.Волков, под редакцией проф. В.Н.Емельянова
		// Москва ФИЗМАТЛИТ 2014.
		doublerealT m[5];
		if (iorder == 3) {
			if (1) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.2075;
				m[1] = 0.5915;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
			else {
				// Направленное огрубление.
				m[0] = 0.2239;
				m[1] = 0.5653;
				m[2] = 1.0;
				m[3] = 0.0; // не используется в трёхшаговом методе.
				m[4] = 0.0; // не используется в трёхшаговом методе.
			}
		}

		if (iorder == 5) {
			if (0) {
				// ПИОНЕР - 1 лучший выбор.
				// Полное огрубление.
				m[0] = 0.0962;
				m[1] = 0.2073;
				m[2] = 0.3549;
				m[3] = 0.6223;
				m[4] = 1.0;
			}
			else {
				// Направленное огрубление.
				m[0] = 0.0870;
				m[1] = 0.1892;
				m[2] = 0.3263;
				m[3] = 0.5558;
				m[4] = 1.0;
			}
		}

		// istart - начальная позиция ненулевых элементов в матрице А.
		// iend - конечная позиция ненулевых элементов в матрице А.


		doublerealT rn = (doublerealT)(iendq - istartq + 1);

		integer startpos = istartq + iadd;
		integer endpos = iendq + iadd;


		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy





		for (integer inumber_step_Runge_Kutt = 0; inumber_step_Runge_Kutt < iorder - 1; inumber_step_Runge_Kutt++) {

			if (idirect == 1) {

				// Восходящая ветвь : сначала F потом C.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == false) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr]*b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
					}
				}


				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == true) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr]*b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

					}
				}

			}
			else {
				// Нисходящая ветвь : сначала C потом F.

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == true) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr]*b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
					}
				}


				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					integer istr = ii - iadd;
					x_jacoby_buffer[istr] = x[istr];
				}

				//#pragma loop(hint_parallel(8))
#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == false) {

						integer istr = ii - iadd;
						doublerealT rold = x_jacoby_buffer[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr]*b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];

							//rsum += -Amat.aij[ii1]*x[ipos]; // experiment
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;


						x[istr] = (1.0 - m[inumber_step_Runge_Kutt]) * x[istr] + m[inumber_step_Runge_Kutt] * rold; // 21
						//x[istr] = m[inumber_step_Runge_Kutt]*x[istr] + (1.0-m[inumber_step_Runge_Kutt]) * rold; // 23

						//x[istr] = (1.0 - m[inumber_step_Runge_Kutt])*x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];
						//x[istr] = x_jacoby_buffer[istr] + m[inumber_step_Runge_Kutt] * x[istr];

					}
				}


			}

		}

	}
	else {
		// Неправильный порядок метода Рунге - Кутты. 
		// Предусмотрены только третий и пятый порядки.
		printf("order Runge Kutt method is bad...\n ");
		system("pause");
		exit(1);
	}


} // Runge_Kutt_3or5


typedef struct TDATA_BASE_GMRES {
	doublereal* val = nullptr;
	integer* col_ind = nullptr;
	integer* row_ptr = nullptr;
} DATA_BASE_GMRES;

DATA_BASE_GMRES* dbgmres_smoother = nullptr;

// 13.09.2017.
template <typename doublerealT>
void gmres_smoother(Ak2 &Amat, integer istartq, integer iendq, 
	doublerealT* &x, doublerealT* &b, 
	integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, integer ilevel)
{
	integer n = iendq - istartq + 1;

	// Разреженная матрица СЛАУ
	// в CRS формате.

	doublerealT* val=nullptr;
	integer* col_ind = nullptr, *row_ptr = nullptr;
	doublerealT* dX0 = new doublerealT[n];
	doublerealT* dV = new doublerealT[n];

	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		dX0[ii - startpos] = x[istr];
		dV[ii - startpos] = b[istr];
	}

	integer nnz = row_ptr_end[endpos] - row_ptr_start[startpos] + 1;
	if (dbgmres_smoother[ilevel].val == nullptr) {
		dbgmres_smoother[ilevel].val = new doublerealT[nnz];
		dbgmres_smoother[ilevel].col_ind = new integer[nnz];
		dbgmres_smoother[ilevel].row_ptr = new integer[n + 1];
		for (integer i_1 = 0; i_1 <= n; i_1++) dbgmres_smoother[ilevel].row_ptr[i_1] = nnz;

		for (integer ii = startpos; ii <= endpos; ii++) {
			integer istr = ii - iadd;
			integer is1 = row_ptr_start[ii];// +1;
			integer is2 = row_ptr_end[ii];


			// Распараллеливание почемуто тормозит очень сильно.
			//#pragma omp parallel for reduction(+:rsum)
			for (integer ii1 = is1; ii1 <= is2; ii1++)
			{
				///printf("is1=%d is2=%d\n",is1,is2);
				if (iadd > 20) {
					//printf("i==%d j==%d\n", Amat.i[ii1], Amat.j[ii1]); getchar();

					//printf("%d \n", ii1 - row_ptr_start[startpos]);
					//getchar();
				}
				dbgmres_smoother[ilevel].col_ind[ii1 - row_ptr_start[startpos]] = Amat.j[ii1] - 1;
				dbgmres_smoother[ilevel].val[ii1 - row_ptr_start[startpos]] = Amat.aij[ii1];
				if (Amat.j[ii1] != istr/*Amat.i[ii1]*/) {
					dbgmres_smoother[ilevel].val[ii1 - row_ptr_start[startpos]] *= 1.0;
				}
				else {
					dbgmres_smoother[ilevel].val[ii1 - row_ptr_start[startpos]] = 1.0 / dbgmres_smoother[ilevel].val[ii1 - row_ptr_start[startpos]];
				}
				dbgmres_smoother[ilevel].row_ptr[istr/*Amat.i[ii1]*/ - 1] = min(dbgmres_smoother[ilevel].row_ptr[istr/*Amat.i[ii1]*/ - 1],
					ii1 - row_ptr_start[startpos]);
			}
		}
	}

	/*
	if (0) {
		// Матрица col_ind, row_ptr, val собрана верно!!! 14 сентября 2017.

		doublerealT omega = 1.0; // initialize.

		// За подробностями смотри книгу Патрика Роуча стр. 183.
		doublerealT rn = (doublerealT)(iendq - istartq + 1);
		optimal_omega(rn, omega);

		for (integer i_1 = 0; i_1 < n; i_1++) {
			doublerealT r = dV[i_1];
			doublerealT ap = 0.0;
			for (integer j_1 = row_ptr[i_1]; j_1 <= row_ptr[i_1 + 1] - 1; j_1++) {
				if (i_1 != col_ind[j_1]) {
					r += -val[j_1] * dX0[col_ind[j_1]];
					//printf("%e ",);
				}
				else {
					ap = val[j_1];
				}
			}
			dX0[i_1] = (1.0 - omega)*dX0[i_1] + ((omega)*(r)) / ap;
		}
	}
	*/

	// Методы GMRES и BiCGStab расходятся в качестве сглаживателей и 
	// не обеспечивают не то что сглаживания а даже дают расходимость.
	// GMRES Саад и Шульц. [1986]
	integer maxit = 20;
	integer m_restart = my_amg_manager.m_restart;// 20; // Так рекомендовано в статье.
	gmres(n, dbgmres_smoother[ilevel].val, dbgmres_smoother[ilevel].col_ind, dbgmres_smoother[ilevel].row_ptr, dV, dX0, maxit, m_restart);
	//Bi_CGStabCRS_smoother(n, val, col_ind, row_ptr, dV, dX0, maxit);


	// Возвращение результата расчёта.
	for (integer ii = startpos; ii <= endpos; ii++) {
		integer istr = ii - iadd;
		x[istr] = dX0[ii - startpos];
	}

	delete[] dX0;
	delete[] dV;
	delete[] val;
	delete[] col_ind;
	delete[] row_ptr;
}

// smoother.
// 5 июня 2017 добавлено CF-Jacobi smoothing (F - smoothing).
// 5 января 2016 с использованием формулы из книги Патрика Роуча.
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelqsor2(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, bool* &F_false_C_true, integer idirect)
{
	// F_false_C_true - нумерация начинается с 1.
	// idirect==0 douwn
	// idirect==1 up

	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

							// За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega); //28.07.2016
							  //omega = 0.7;

							  //if (isorintmemo == iadd) {
							  // Это точно не первый раз
							  //bfirst = false;
							  //}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	if (omega < 1.0) {
		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy



		if (idirect == 1) {
			// Восходящая ветвь.

			// Сначала F потом C.

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				if (F_false_C_true[ii] == false) { // F nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] *x_jacoby_buffer[ipos];
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					//if (is1 <= is2) {
						// Только если это не условие Дирихле применяем релаксацию.
					x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
				//}
				}
			}

			// update.
			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				if (F_false_C_true[ii] == true) { // C nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] *x_jacoby_buffer[ipos];
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					//if (is1 <= is2) {
						// Только если это не условие Дирихле применяем релаксацию.
					x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
				//}
				}
			}



		}
		else {
			// idirect==0
			// Сначала С потом F.

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				if (F_false_C_true[ii] == true) { // C nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] *x_jacoby_buffer[ipos];
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					//if (is1 <= is2) {
						// Только если это не условие Дирихле применяем релаксацию.
					x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
				//}
				}
			}


			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				if (F_false_C_true[ii] == false) { // F nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1]*x_jacoby_buffer[ipos];
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					//if (is1 <= is2) {
						// Только если это не условие Дирихле применяем релаксацию.
					x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
				//}
				}
			}

		}


	}
	else {

		// 3 ноября 2016. Симметричный метод Гаусса-Зейделя.
		if (isimmetricGS_switch == 0) {
			// 3 ноября 2016 было экспериментально доказано на BSKDmitrii что симметричный метод Гаусса - Зейделя РАСХОДИТСЯ.
			// НЕЛЬЗЯ ИСПОЛЬЗОВАТЬ СИММЕТРИЧНЫЙ МЕТОД ГАУССА -ЗЕЙДЕЛЯ. ИСПОЛЬЗУЙТЕ ОБЫЧНЫЙ МЕТОД ЗЕЙДЕЛЯ.

			if (idirect == 1) {
				// Восходящая ветвь.



				// Сначала F потом C.

//----->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == false) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// Здесь всё в каше и позитив и негатив коннектионшс.
								rsum += -Amat.aij[ii1]*x[ipos];
								//}
								//else {
								// не работает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

								// игнорирование positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1]*x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// не работает.
									// Вероятно из-за того что так нарушен баланс.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1]*x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
					//}
					}
				}

				//---->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == true) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// Здесь всё в каше и позитив и негатив коннектионшс.
								rsum += -Amat.aij[ii1]*x[ipos];
								//}
								//else {
								// не работает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

								// игнорирование positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1]*x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// не работает.
									// Вероятно из-за того что так нарушен баланс.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1]*x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
					//}
					}
				}

			}
			else {

				// idirect==0
				// Сначала С потом F.

//---->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == true) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// Здесь всё в каше и позитив и негатив коннектионшс.
								rsum += -Amat.aij[ii1]*x[ipos];
								//}
								//else {
								// не работает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

								// игнорирование positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1]*x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// не работает.
									// Вероятно из-за того что так нарушен баланс.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1]*x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
					//}
					}
				}

				//---->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == false) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// Здесь всё в каше и позитив и негатив коннектионшс.
								rsum += -Amat.aij[ii1]*x[ipos];
								//}
								//else {
								// не работает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

								// игнорирование positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1]*x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// не работает.
									// Вероятно из-за того что так нарушен баланс.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1]*x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
					//}
					}
				}

			}

			// Ни в коем случае не переключать. 3 ноября 2016.
			//isimmetricGS_switch = 1;

		}
		else {

			if (idirect == 1) {
				// Восходящая ветвь.

				// Сначала F потом C.

//---->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if (F_false_C_true[ii] == false) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1]*x[ipos];
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
					//}
					}
				}


				//----->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if (F_false_C_true[ii] == true) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1]*x[ipos];
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
					//}
					}
				}



			}
			else {
				// idirect==0
				// Сначала С потом F.

//--->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if (F_false_C_true[ii] == true) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1]*x[ipos];
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
					//}
					}
				}


				//--->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if (F_false_C_true[ii] == false) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1]*x[ipos];
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega)*rold; // this is SOR
					//}
					}
				}

			}

			//isimmetricGS_switch = 0;  // смена направления.
		}
	}


} // seidelqsor2

// smoother.
// 5 июня 2017 добавлено CF-Jacobi smoothing (F - smoothing).
// 5 января 2016 с использованием формулы из книги Патрика Роуча.
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelqsor2(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, bool*& F_false_C_true, integer idirect, doublereal* &diag_minus_one)
{
	// F_false_C_true - нумерация начинается с 1.
	// idirect==0 douwn
	// idirect==1 up

	// istart - начальная позиция ненулевых элементов в матрице А.
	// iend - конечная позиция ненулевых элементов в матрице А.
	// sor 1.855 использовать нельзя мы имеем расходимость вычислительного процесса.
	// Попробуем нижнюю релаксацию.
	// Расходимость имеем и при нижней релаксации. 0.8

	// BSKDmitrii
	// omega   iter  time,s
	// 1.0 106 43
	// 1.1 98 42
	// 1.15 94 40 best
	// 1.2 90 40
	// 1.225 413 1min 37s
	// 1.25 divergence detected
	// 1.3 divergence detected

	doublerealT omega = 1.0; // initialize.

							// За подробностями смотри книгу Патрика Роуча стр. 183.
	doublerealT rn = (doublerealT)(iendq - istartq + 1);
	optimal_omega(rn, omega); //28.07.2016
							  //omega = 0.7;

							  //if (isorintmemo == iadd) {
							  // Это точно не первый раз
							  //bfirst = false;
							  //}
	integer startpos = istartq + iadd;
	integer endpos = iendq + iadd;

	if (omega < 1.0) {
		if (bfirst_jacoby_start) {
			x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
			i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
			bfirst_jacoby_start = false;
		}
		else {
			// Перевыделение оперативной памяти в случае nu1==0.
			if (i_x_jacoby_buffer_pool_size < 3 * (endpos - startpos + 1)) {
				if (x_jacoby_buffer != nullptr) {
					delete[] x_jacoby_buffer;
					x_jacoby_buffer = nullptr;
					x_jacoby_buffer = new doublereal[3 * (endpos - startpos + 1)];
					i_x_jacoby_buffer_pool_size = 3 * (endpos - startpos + 1);
					bfirst_jacoby_start = false;
				}
			}
		}
		// copy



		if (idirect == 1) {
			// Восходящая ветвь.

			// Сначала F потом C.

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				if (F_false_C_true[ii] == false) { // F nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = diag_minus_one[istr]*b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					//if (is1 <= is2) {
						// Только если это не условие Дирихле применяем релаксацию.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				//}
				}
			}

			// update.
			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				if (F_false_C_true[ii] == true) { // C nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = diag_minus_one[istr]*b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					//if (is1 <= is2) {
						// Только если это не условие Дирихле применяем релаксацию.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				//}
				}
			}



		}
		else {
			// idirect==0
			// Сначала С потом F.

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				if (F_false_C_true[ii] == true) { // C nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = diag_minus_one[istr]*b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					//if (is1 <= is2) {
						// Только если это не условие Дирихле применяем релаксацию.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				//}
				}
			}


			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				integer istr = ii - iadd;
				x_jacoby_buffer[istr] = x[istr];
			}

			//#pragma loop(hint_parallel(8))
#pragma omp parallel for
			for (integer ii = startpos; ii <= endpos; ii++) {
				if (F_false_C_true[ii] == false) { // F nodes

					integer istr = ii - iadd;
					doublerealT rold = x_jacoby_buffer[istr];

					// 13.07.2016
					doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

					x[istr] = diag_minus_one[istr]*b[istr];

					doublerealT rsum = 0.0;
					integer is1 = row_ptr_start[ii] + 1;
					integer is2 = row_ptr_end[ii];
					// Распараллеливание почемуто тормозит очень сильно.
					//#pragma omp parallel for reduction(+:rsum)
					for (integer ii1 = is1; ii1 <= is2; ii1++)
					{
						//x[istr] += -Amat.aij[ii1]*x_jacoby_buffer[Amat.j[ii1]];
						integer ipos = Amat.j[ii1];
						// 13.07.2016
						// игнорирование positive connections.
						//if ((Amat.aij[ii1] < 0.0)) {
						rsum += -Amat.aij[ii1] * x_jacoby_buffer[ipos];
						//}
						//else {
						// не рабтает.
						//	ap_now += Amat.aij[ii1];
						//}
					}
					x[istr] += rsum;
					//x[istr] *= Amat.aij[row_ptr_start[ii]];
					// 13.07.2016
					x[istr] /= ap_now;

					// Запускается только после первого раза сразу запускать нельзя
					// т.к. начальное значение не является адекватным.
					//if (is1 <= is2) {
						// Только если это не условие Дирихле применяем релаксацию.
					x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
				//}
				}
			}

		}


	}
	else {

		// 3 ноября 2016. Симметричный метод Гаусса-Зейделя.
		if (isimmetricGS_switch == 0) {
			// 3 ноября 2016 было экспериментально доказано на BSKDmitrii что симметричный метод Гаусса - Зейделя РАСХОДИТСЯ.
			// НЕЛЬЗЯ ИСПОЛЬЗОВАТЬ СИММЕТРИЧНЫЙ МЕТОД ГАУССА -ЗЕЙДЕЛЯ. ИСПОЛЬЗУЙТЕ ОБЫЧНЫЙ МЕТОД ЗЕЙДЕЛЯ.

			if (idirect == 1) {
				// Восходящая ветвь.



				// Сначала F потом C.

//----->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == false) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = diag_minus_one[istr]*b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// Здесь всё в каше и позитив и негатив коннектионшс.
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// не работает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

								// игнорирование positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// не работает.
									// Вероятно из-за того что так нарушен баланс.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1] * x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

				//---->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == true) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// Здесь всё в каше и позитив и негатив коннектионшс.
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// не работает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

								// игнорирование positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// не работает.
									// Вероятно из-за того что так нарушен баланс.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1] * x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

			}
			else {

				// idirect==0
				// Сначала С потом F.

//---->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == true) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// Здесь всё в каше и позитив и негатив коннектионшс.
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// не работает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

								// игнорирование positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// не работает.
									// Вероятно из-за того что так нарушен баланс.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1] * x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

				//---->#pragma omp parallel for
				for (integer ii = startpos; ii <= endpos; ii++) {
					if (F_false_C_true[ii] == false) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];
						// 28.01.2017
						//ap_now = 0.0;

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							if (1) {
								// игнорирование positive connections.
								//if ((Amat.aij[ii1] < 0.0)) {
								// Здесь всё в каше и позитив и негатив коннектионшс.
								rsum += -Amat.aij[ii1] * x[ipos];
								//}
								//else {
								// не работает.
								//	ap_now += Amat.aij[ii1];
								//}
							}
							else {
								// Здесь я перемудрил. По видимому нужен Крукиер : TKM, TKM1, TKM2.

								// игнорирование positive connections.
								if ((Amat.aij[ii1] < 0.0)) {
									rsum += -Amat.aij[ii1] * x[ipos];
									//ap_now += fabs(Amat.aij[ii1]);
								}
								else {
									// не работает.
									// Вероятно из-за того что так нарушен баланс.
									//ap_now += fabs(Amat.aij[ii1]);
									if (fabs(x[ipos]) > fabs(x[istr])) {
										ap_now += fabs(Amat.aij[ii1]);
									}
									else rsum += -Amat.aij[ii1] * x[ipos];
								}
							}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

			}

			// Ни в коем случае не переключать. 3 ноября 2016.
			//isimmetricGS_switch = 1;

		}
		else {

			if (idirect == 1) {
				// Восходящая ветвь.

				// Сначала F потом C.

//---->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if (F_false_C_true[ii] == false) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}


				//----->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if (F_false_C_true[ii] == true) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}



			}
			else {
				// idirect==0
				// Сначала С потом F.

//--->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if (F_false_C_true[ii] == true) { // C nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}


				//--->#pragma omp parallel for
				for (integer ii = endpos; ii >= startpos; ii--) {
					if (F_false_C_true[ii] == false) { // F nodes

						integer istr = ii - iadd;
						doublerealT rold = x[istr];

						// 13.07.2016
						doublerealT ap_now = 1.0 / Amat.aij[row_ptr_start[ii]];

						x[istr] = diag_minus_one[istr] * b[istr];

						doublerealT rsum = 0.0;
						integer is1 = row_ptr_start[ii] + 1;
						integer is2 = row_ptr_end[ii];
						// Распараллеливание почемуто тормозит очень сильно.
						//#pragma omp parallel for reduction(+:rsum)
						for (integer ii1 = is1; ii1 <= is2; ii1++)
						{
							//x[istr] += -Amat.aij[ii1]*x[Amat.j[ii1]];
							integer ipos = Amat.j[ii1];
							// 13.07.2016
							// игнорирование positive connections.
							//if ((Amat.aij[ii1] < 0.0)) {
							rsum += -Amat.aij[ii1] * x[ipos];
							//}
							//else {
							// не рабтает.
							//	ap_now += Amat.aij[ii1];
							//}
						}
						x[istr] += rsum;
						//x[istr] *= Amat.aij[row_ptr_start[ii]];
						// 13.07.2016
						x[istr] /= ap_now;

						// Запускается только после первого раза сразу запускать нельзя
						// т.к. начальное значение не является адекватным.
						//if (is1 <= is2) {
							// Только если это не условие Дирихле применяем релаксацию.
						x[istr] = omega * x[istr] + (1.0 - omega) * rold; // this is SOR
					//}
					}
				}

			}

			//isimmetricGS_switch = 0;  // смена направления.
		}
	}


} // seidelqsor2

// smoother.
  // 9 september 2015 and 4 june 2017.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2 &Amat, integer istartq, integer iendq,
	doublerealT* &x, doublerealT* &b, 
	integer * &row_ptr_start, integer * &row_ptr_end,
	integer iadd, bool* &F_false_C_true, integer idirect, integer ilevel)
{
	
	if (my_amg_manager.b_spai0) {
		spai0_smoother<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
	}
	else if (my_amg_manager.b_gmres) {
		if (my_amg_manager.bCFJacoby == true) {
			seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, F_false_C_true, idirect);
		}
		else {
			seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
		}
		// gmres smoother
		gmres_smoother<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd,ilevel);
	}
	else {
		// лучший выбор
		if (my_amg_manager.bCFJacoby == true) {
			if ((my_amg_manager.iRunge_Kutta_smoother == 3) || (my_amg_manager.iRunge_Kutta_smoother == 5)) {
				// Трёхшаговый метод Рунге-Кутты.
				integer iorder = my_amg_manager.iRunge_Kutta_smoother;
				Runge_Kutt_3or5<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iorder, F_false_C_true, idirect);
			}
			else {

				seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, F_false_C_true, idirect);
			}
		}
		else {
			if ((my_amg_manager.iRunge_Kutta_smoother == 3) || (my_amg_manager.iRunge_Kutta_smoother == 5)) {
				// Трёхшаговый метод Рунге-Кутты.
				integer iorder = my_amg_manager.iRunge_Kutta_smoother;
				Runge_Kutt_3or5<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iorder);
			}
			else {

				seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);
			}

		}
	}

} // seidelq

// smoother.
  // 9 september 2015 and 4 june 2017.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq,
	doublerealT*& x, doublerealT*& b,
	integer*& row_ptr_start, integer*& row_ptr_end, 
	integer iadd, bool*& F_false_C_true, integer idirect,
	doublereal* &diag_minus_one, integer ilevel)
{
	if (my_amg_manager.b_spai0) {
		spai0_smoother<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
	}
	else if (my_amg_manager.b_gmres) {
		if (my_amg_manager.bCFJacoby == true) {
			seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, F_false_C_true, idirect, diag_minus_one);
		}
		else {
			seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
		}
		// gmres smoother
		gmres_smoother<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, ilevel);
	}
	else {
		// лучший выбор
		if (my_amg_manager.bCFJacoby == true) {
			if ((my_amg_manager.iRunge_Kutta_smoother == 3) || (my_amg_manager.iRunge_Kutta_smoother == 5)) {
				// Трёхшаговый метод Рунге-Кутты.
				integer iorder = my_amg_manager.iRunge_Kutta_smoother;
				Runge_Kutt_3or5<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iorder, F_false_C_true, idirect, diag_minus_one);
			}
			else {

				seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, F_false_C_true, idirect, diag_minus_one);
			}
		}
		else {
			if ((my_amg_manager.iRunge_Kutta_smoother == 3) || (my_amg_manager.iRunge_Kutta_smoother == 5)) {
				// Трёхшаговый метод Рунге-Кутты.
				integer iorder = my_amg_manager.iRunge_Kutta_smoother;
				Runge_Kutt_3or5<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, iorder, diag_minus_one);
			}
			else {

				seidelqsor2<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
			}

		}
	}

} // seidelq

// smoother.
// 16 jan 2016.  Seidel q -quick SOR + parallel
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelq(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, bool* &bnested_desection, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd)
{
	seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, bnested_desection, row_ptr_start, row_ptr_end, iadd);
	//seidelqsor2Pcpu(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

} // seidelq

// smoother.
// 16 jan 2016.  Seidel q -quick SOR + parallel
// 9 september 2015.
// q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, bool*& bnested_desection, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, doublereal*& diag_minus_one)
{
	seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, bnested_desection, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
	//seidelqsor2Pcpu(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

} // seidelq

  // smoother.
  // 16 jan 2016.  Seidel q -quick SOR + parallel
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2 &Amat, integer istartq, integer iendq, doublerealT* &x, doublerealT* &b, bool* &bnested_desection, integer * &row_ptr_start, integer * &row_ptr_end, integer iadd, bool* &F_false_C_true, integer idirect)
{
	// , bool* &F_false_C_true, integer idirect Заглушка, параметры не используются.
	// Внимание обратная совместимость.

	seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, bnested_desection, row_ptr_start, row_ptr_end, iadd);
	//seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

} // seidelq

 // smoother.
  // 16 jan 2016.  Seidel q -quick SOR + parallel
  // 9 september 2015.
  // q - quick.
template <typename doublerealT>
void seidelq(Ak2& Amat, integer istartq, integer iendq, doublerealT*& x, doublerealT*& b, bool*& bnested_desection, integer*& row_ptr_start, integer*& row_ptr_end, integer iadd, bool*& F_false_C_true, integer idirect, doublereal*& diag_minus_one)
{
	// , bool* &F_false_C_true, integer idirect Заглушка, параметры не используются.
	// Внимание обратная совместимость.

	seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, bnested_desection, row_ptr_start, row_ptr_end, iadd, diag_minus_one);
	//seidelqsor2Pcpu<doublerealT>(Amat, istartq, iendq, x, b, row_ptr_start, row_ptr_end, iadd);

} // seidelq



// 08.01.2018 Перенесено в отдельную функцию, т.к. используется неоднократно. Передаётся большое число параметров.
//  A*z76=s76; 891 строка.
template <typename doublerealT>
void V_cycle_solve(Ak2 &Amat, doublereal* &z76, doublereal* &s76, bool process_flow_logic, integer* &row_ptr_start,
	integer* &row_ptr_end, doublerealT** &residual_fine, doublerealT** &diag, integer* n_a, bool bonly_serial,
	doublerealT process_flow_beta, bool* &F_false_C_true, integer &nu1, integer &nu2, integer bILU2smoother,
	integer ilevel, integer inumberVcyclelocbicgstab, integer imyinit, const integer idim_diag,
	LEVEL_ADDITIONAL_DATA* &milu2, LEVEL_ADDITIONAL_DATA0* milu0, bool** &nested_desection, 
	Ak1* &P, // prolongation он же restriction (метод не чувствителен к сортировке).
	integer* nnz_aRP, doublerealT** &residual_coarse, integer igam, integer* nnz_a,
	doublerealT** &error_approx_coarse, doublerealT dapply_ilu_max_pattern_size,
	doublerealT process_flow_alpha, doublerealT** &error_approx_fine,
	integer nFinestSweeps) {

	// Один V - цикл алгебраического многосеточного метода.
	// A*z76=s76;

	const integer ZERO_INIT = 0;
	const integer RANDOM_INIT = 1;// надо увеличивать nu1, nu2 с 1,2 до 5 наверно.

	for (integer i_13 = 0; i_13 < inumberVcyclelocbicgstab; i_13++)
	{


		doublerealT R0_0 = 0.0;
		doublerealT Rprev_0 = 0.0, Rnext_0 = 0.0;
		if (process_flow_logic) {
			// calculate initial residual.
			//residualq(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
			R0_0 = norma(residual_fine[0], n_a[0]);
			Rprev_0 = R0_0;

			// smother
			integer iter = 0;
			for (iter = 0; iter < nu1; iter++) {
				//quick seidel
				if (bonly_serial) {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, 0);
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 0);
				}
				//residualq(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
				Rnext_0 = norma(residual_fine[0], n_a[0]);
				// this is process flow logic
				if (Rnext_0 > process_flow_beta*Rprev_0) {
					// Смысл модификации в том что мы экономим итерации на пресмутере.
					break; // досрочно опускаемся на следующий уровень если он есть конечно.
				}
				else {
					Rprev_0 = Rnext_0;
				}
			}
			if (iter == nu1) {
				printf("level 0 limit presmother iteration is reached\n");
			}

		}
		else {
			// smoother
			for (integer iter = 0; iter < nu1; iter++) {
				//seidel(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0]);
				//quick seidel
				if (bonly_serial) {
					if (bILU2smoother == 1) {
						// ILU0
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0,0);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							milu0[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
						}
						lusol_1patchforRUMBA(n_a[0], milu0[0].zbuf, milu0[0].zbuf2, milu0[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							z76[i43 + 1] += milu0[0].zbuf2[i43 + 1];
						}
					}
					else if ((bILU2smoother == 2) || (my_amg_manager.iFinnest_ilu == 1)) {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0,0);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							milu2[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
						}
						lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							z76[i43 + 1] += milu2[0].zbuf2[i43 + 1];
						}
					}
					else {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0,0);
					}
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 0);
				}
			}
		}

		//exporttecplot(x, n);

		move_down(nu1, nu2);

		if (!process_flow_logic) {
			// residual_r
			//doublerealT *residual_fine[0] = new doublerealT[n_a[0] + 1];
			//residual<doublereal>(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0], residual_fine[0]);
			//residualq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
		}



		//if (iprohod%5==0) getchar();
		if (ilevel > 1) {

			//doublerealT *residual_coarse = new doublerealT[n_a[1] + 1];
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[1]; ii++) {
				residual_coarse[0][ii] = 0.0;
			}

			// restriction
			restriction(P, 1, nnz_aRP[0],  residual_fine[0], residual_coarse[0],  n_a[1]);


			// Amat*e=r;
			//doublerealT* error_approx_coarse = new doublerealT[n_a[1] + 1];
			if ((imyinit == ZERO_INIT)) {
#pragma omp parallel for
				for (integer ii = 1; ii <= n_a[1]; ii++) {
					error_approx_coarse[0][ii] = 0.0;
				}
			}

			if ((imyinit == RANDOM_INIT)) {
				// (1,110); (0.8, 37); (0.7, 29); (0.6, 25); (0.5, 20); (0.4, 17); (0.3, 18); (0.0, 19);
#pragma omp parallel for
				for (integer ii = 1; ii <= n_a[1]; ii++) {
					error_approx_coarse[0][ii] = 0.4*fabs(residual_coarse[0][ii])*(rand() / ((doublerealT)RAND_MAX));
				}
			}


			for (integer i_37 = 1; i_37 <= igam; i_37++) {
				doublerealT R0_1 = 0.0;
				doublerealT Rprev_1 = 0.0, Rnext_1 = 0.0;
				if (process_flow_logic) {
					// calculate initial residual.
					//residualq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1]);
					residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);

					R0_1 = norma(residual_fine[1], n_a[1]);
					Rprev_1 = R0_1;

					// smother
					integer iter = 0;
					for (iter = 0; iter < nu1; iter++) {
						//quick seidel
						if (bonly_serial) {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0,1);
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0);
						}

						//residualq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1]);
						residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);

						Rnext_1 = norma(residual_fine[1], n_a[1]);
						// this is process flow logic
						if (Rnext_1 > process_flow_beta*Rprev_1) {
							// Смысл модификации в том что мы экономим итерации на пресмутере.
							break; // досрочно опускаемся на следующий уровень если он есть конечно.
						}
						else {
							Rprev_1 = Rnext_1;
						}
					}

					if (iter == nu1) {
						printf("level 1 limit presmother iteration is reached\n");
					}

				}
				else {

					// pre smothing
					for (integer iter = 0; iter < nu1; iter++) {
						//seidel(Amat, 1 + 2 * nnz_a[0], 2 * nnz_a[0] + nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1]);
						if (bonly_serial) {
							if (bILU2smoother == 1) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0,1);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									milu0[1].zbuf[i43 + 1] = residual_fine[1][i43 + 1];
								}
								lusol_1patchforRUMBA(n_a[1], milu0[1].zbuf, milu0[1].zbuf2, milu0[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									error_approx_coarse[0][i43 + 1] += milu0[1].zbuf2[i43 + 1];
								}


							}
							else if (1 && bILU2smoother == 2) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0,1);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									milu2[1].zbuf[i43 + 1] = residual_fine[1][i43 + 1];
								}
								lusol_1patchforRUMBA(n_a[1], milu2[1].zbuf, milu2[1].zbuf2, milu2[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									error_approx_coarse[0][i43 + 1] += milu2[1].zbuf2[i43 + 1];
								}


							}
							else {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0,1);
							}
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0);
						}
					}
				}





				doublerealT* R0_21 = new doublerealT[idim_diag];
				doublerealT* Rprev_21 = new doublerealT[idim_diag];
				doublerealT* Rnext_21 = new doublerealT[idim_diag];

				//3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
				for (integer i_id_level_local = 2; i_id_level_local < idim_diag; i_id_level_local++) {

					move_down(nu1, nu2);


					if (ilevel > i_id_level_local) {
						// Впервые создан универсальный програмный код на произвольное количество уровней вложенности 
						// в solution phase.
						// Нужно только написать универсальный обработчик для R0_21.
						// 4 декабря 2016. 

						// residual
						if (!process_flow_logic) {
							integer in_a_loc = 0;
							for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
								// sum(n_a[0]+n_a[19];
								in_a_loc += n_a[i_72];
							}
							//doublerealT *residual_fine[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local-1] + 1];
							//residual(Amat, 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6]+ nnz_a[7]+ nnz_a[8]+ nnz_a[9]+ nnz_a[10]+ nnz_a[11]+ nnz_a[12]+ nnz_a[13]+ nnz_a[14]+ nnz_a[15]+ nnz_a[16]+nnz_a[17]+nnz_a[18]+nnz_a[19], nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7]+ nnz_a[8]+ nnz_a[9]+ nnz_a[10]+ nnz_a[11]+ nnz_a[12]+ nnz_a[13]+ nnz_a[14]+ nnz_a[15]+ nnz_a[16]+ nnz_a[17]+nnz_a[18]+nnz_a[19]+nnz_a[20], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], flag, n_a[i_id_level_local - 1], residual_fine[i_id_level_local - 1]);
							//residualq(Amat, 1, n_a[i_id_level_local-1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end,  in_a_loc , residual_fine[i_id_level_local - 1]);
							//residualq2(Amat, 1, n_a[i_id_level_local-1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end, in_a_loc, residual_fine[i_id_level_local - 1], diag[i_id_level_local - 1]);
							residualq2(Amat, 1, n_a[i_id_level_local - 1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end, in_a_loc, residual_fine[i_id_level_local - 1], diag[i_id_level_local - 1]);

						}

						//doublerealT *residual_coarse[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local] + 1];

						// restriction

						integer innz_aRP_loc = 1;
						for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
							innz_aRP_loc += nnz_aRP[i_72];
						}
						//restriction(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16] + nnz_aRP[17] + nnz_aRP[18] + nnz_aRP[19], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16] + nnz_aRP[17] + nnz_aRP[18] + nnz_aRP[19] + nnz_aRP[20], flag, residual_fine[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], n_a[i_id_level_local - 1], n_a[i_id_level_local]);
						restriction(P, innz_aRP_loc, innz_aRP_loc - 1 + nnz_aRP[i_id_level_local - 1],  residual_fine[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], n_a[i_id_level_local]);



						// Amat*e=r;
						//doublerealT* error_approx_coarse[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local] + 1];

						if (0) {
							//if (icount_V_cycle == 1) {
							//	for (integer ii = 1; ii <= n_a[i_id_level_local]; ii++) {
							// обнуление только на первом проходе.
							//		error_approx_coarse[i_id_level_local - 1][ii] = 0.0;
							//error_approx_coarse[i_id_level_local - 1][ii] = (rand() / ((doublerealT)RAND_MAX));
							//	}
							//}
							// На задаче BSK_Dmitrii random оптимизация оказалась гораздо хуже, но сходимость тоже присутствовала.
#pragma omp parallel for
							for (integer ii = 1; ii <= n_a[i_id_level_local]; ii++) {
								// 0.4*fabs(residual_coarse[i_id_level_local - 1][ii]) - амплитуда.
								// 0.4 - демпфирующий множитель.
								// diag[i_id_level_local]
								// (1,110); (0.8, 37); (0.7, 29); (0.6, 25); (0.5, 20); (0.4, 17); (0.3, 18); (0.0, 19);
								error_approx_coarse[i_id_level_local - 1][ii] = 0.4*fabs(residual_coarse[i_id_level_local - 1][ii])*(rand() / ((doublerealT)RAND_MAX));
							}
						}
						else {
							if (imyinit == ZERO_INIT) {
#pragma omp parallel for
								for (integer ii = 1; ii <= n_a[i_id_level_local]; ii++) {
									error_approx_coarse[i_id_level_local - 1][ii] = 0.0;
								}
							}
							if (imyinit == RANDOM_INIT) {
#pragma omp parallel for
								for (integer ii = 1; ii <= n_a[i_id_level_local]; ii++) {
									// (1,110); (0.8, 37); (0.7, 29); (0.6, 25); (0.5, 20); (0.4, 17); (0.3, 18); (0.0, 19);
									error_approx_coarse[i_id_level_local - 1][ii] = 0.4*fabs(residual_coarse[i_id_level_local - 1][ii])*(rand() / ((doublerealT)RAND_MAX));
								}
							}
						}


						//for (integer i_37 = 1; i_37 <= igam; i_37++)
						{
							// pre smothing
							//doublerealT R0_21 = 0.0;
							//doublerealT Rprev_21 = 0.0, Rnext_21 = 0.0;
							R0_21[i_id_level_local] = 0.0;
							Rprev_21[i_id_level_local] = 0.0;
							Rnext_21[i_id_level_local] = 0.0;
							if (process_flow_logic) {

								integer in_a_loc = 0;
								for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
									// sum(n_a[0]+n_a[19];
									in_a_loc += n_a[i_72];
								}

								// calculate initial residual.
								//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end,in_a_loc + n_a[i_id_level_local-1] , residual_fine[i_id_level_local]);
								residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);



								R0_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
								Rprev_21[i_id_level_local] = R0_21[i_id_level_local];



								// smother
								integer iter = 0;
								integer nu1_count = nu1;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu1_count = 100;

									printf("Direct method is start.../n");

									IMatrix sparseS; // разреженная матрица в формате IMatrix
									initIMatrix(&sparseS, n_a[i_id_level_local]);

									integer istartpos = 1 + in_a_loc + n_a[i_id_level_local - 1];
									integer iendpos = n_a[i_id_level_local] + in_a_loc + n_a[i_id_level_local - 1];

									for (integer i56 = istartpos; i56 <= iendpos; i56++) {

										//setValueIMatrix(&sparseS, Amat.i[row_ptr_start[i56]] - 1, Amat.i[row_ptr_start[i56]] - 1, 1.0 / Amat.aij[row_ptr_start[i56]]);
										setValueIMatrix(&sparseS, i56 - istartpos, i56 - istartpos, 1.0 / Amat.aij[row_ptr_start[i56]]);//3.02.2019
										//if (Amat.i[row_ptr_start[i56]] != i56 - istartpos + 1) {
											//printf("Amat.i[row_ptr_start[i56]]=%lld i56 - istartpos + 1=%lld\n", Amat.i[row_ptr_start[i56]],i56 - istartpos + 1);
											//getchar();
										//}
										const doublerealT nonzeroEPS = 1e-37; // для отделения вещественного нуля

										integer is15 = row_ptr_start[i56] + 1;
										integer is25 = row_ptr_end[i56];

										for (integer ii17 = is15; ii17 <= is25; ii17++)
										{

											if ((fabs(Amat.aij[ii17]) > nonzeroEPS)) {
												//setValueIMatrix(&sparseS, Amat.i[ii17] - 1, Amat.j[ii17] - 1, Amat.aij[ii17]);
												setValueIMatrix(&sparseS, ii17 - is15, Amat.j[ii17] - 1, Amat.aij[ii17]);
											}

										}
									}

									doublereal* dX025 = new doublereal[n_a[i_id_level_local]];
									doublereal* dV25 = new doublereal[n_a[i_id_level_local]];
#pragma omp parallel for
									for (integer i57 = 0; i57 < n_a[i_id_level_local]; i57++) {
										dX025[i57] = error_approx_coarse[i_id_level_local - 1][i57 + 1];
										dV25[i57] = residual_coarse[i_id_level_local - 1][i57 + 1];
									}

									// главный метод, возвращающий решение x,
									// принимает вектор свободных членов b и 
									// квадратную матрицу xO в специальном разреженном формате.
									// реализация без барьера и итерационного уточнения.
									calculateSPARSEgaussArray(&sparseS, dX025, dV25);
#pragma omp parallel for
									for (integer i57 = 0; i57 < n_a[i_id_level_local]; i57++) {
										error_approx_coarse[i_id_level_local - 1][i57 + 1] = dX025[i57];
									}

									delete[] dX025;
									delete[] dV25;


									freeIMatrix(&sparseS);

									//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);
									residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);


									Rnext_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
									// this is process flow logic
									if (Rnext_21[i_id_level_local] > process_flow_beta*Rprev_21[i_id_level_local]) {
										// Смысл модификации в том что мы экономим итерации на пресмутере.
										break; // досрочно опускаемся на следующий уровень если он есть конечно.
									}
									else {
										Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
									}
								}
								else {
									for (iter = 0; iter < nu1_count; iter++) {
										//quick seidel
										if (bonly_serial) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, i_id_level_local);
										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0);
										}

										//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);
										residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);


										Rnext_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
										// this is process flow logic
										if (Rnext_21[i_id_level_local] > process_flow_beta*Rprev_21[i_id_level_local]) {
											// Смысл модификации в том что мы экономим итерации на пресмутере.
											break; // досрочно опускаемся на следующий уровень если он есть конечно.
										}
										else {
											Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
										}
									}
								}

								if (iter == nu1) {
#if doubleintprecision == 1
									printf("level %lld limit presmother iteration is reached\n", i_id_level_local);
#else
									printf("level %d limit presmother iteration is reached\n", i_id_level_local);
#endif

								}

							}
							else {
								integer nu1_count = nu1;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu1_count = 100;
								}
								for (integer iter = 0; iter < nu1_count; iter++) {

									integer in_a_loc = 0;
									for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
										// sum(n_a[0]+n_a[19];
										in_a_loc += n_a[i_72];
									}

									if (bonly_serial) {

										bool bflag56 = false;
										if (my_amg_manager.iFinnest_ilu == 1) {
											if (my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6) {
												doublerealT dn = 1.0*n_a[i_id_level_local];
												doublerealT dnnz = 1.0*nnz_a[i_id_level_local];
												if (dnnz / dn <= dapply_ilu_max_pattern_size) {
													bflag56 = true;
												}
											}
										}

										if (bILU2smoother == 1) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, i_id_level_local);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												milu0[i_id_level_local].zbuf[i43 + 1] = residual_fine[i_id_level_local][i43 + 1];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu0[i_id_level_local].zbuf, milu0[i_id_level_local].zbuf2, milu0[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												error_approx_coarse[i_id_level_local - 1][i43 + 1] += milu0[i_id_level_local].zbuf2[i43 + 1];
											}

										}
										else if (1 && ((bILU2smoother == 2) || bflag56)) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, i_id_level_local);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												milu2[i_id_level_local].zbuf[i43 + 1] = residual_fine[i_id_level_local][i43 + 1];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu2[i_id_level_local].zbuf, milu2[i_id_level_local].zbuf2, milu2[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												error_approx_coarse[i_id_level_local - 1][i43 + 1] += milu2[i_id_level_local].zbuf2[i43 + 1];
											}

										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, i_id_level_local);
										}
									}
									else {
										seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0);
									}
								}
							}
						}

					}

				}







				//center
				// ЭТО сердцевина SOLUTION PHASE.

				// TODO нижний 3.12.2016 (осталось один нижний и один верхний).


				// 21
				for (integer i_id_level_local = idim_diag - 1; i_id_level_local >= 2; i_id_level_local--) {

					if (ilevel > i_id_level_local) {

						{

							integer in_a_loc = 0;
							for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
								// sum(n_a[0]+n_a[19];
								in_a_loc += n_a[i_72];
							}

							// post smoothing
							// doublerealT R0_20 = 0.0;
							///doublerealT Rprev_20 = 0.0, Rnext_20 = 0.0;
							if (process_flow_logic) {
								// calculate initial residual.
								//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc  + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);



								//Rprev_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);


								// smother
								integer iter = 0;
								integer nu2_count = nu2;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu2_count = 100;																												
								}
								else {
									for (iter = 0; iter < nu2_count; iter++) {
										//quick seidel
										if (bonly_serial) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1);
										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1);
										}

										//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local-1], residual_fine[i_id_level_local]);
										residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);


										Rnext_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
										// this is process flow logic
										if (Rnext_21[i_id_level_local] < process_flow_alpha*R0_21[i_id_level_local]) {
											// Смысл модификации в том что мы экономим итерации на пресмутере.
											break; // досрочно опускаемся на следующий уровень если он есть конечно.
										}
										else {
											Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
										}
									}
								}

								if (iter == nu2) {
#if doubleintprecision == 1
									printf("level %lld limit postsmother iteration is reached\n", i_id_level_local);
#else
									printf("level %d limit postsmother iteration is reached\n", i_id_level_local);
#endif

								}

							}
							else {
								integer nu2_count = nu2;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu2_count = 100;																																					
								}
								for (integer iter = 0; iter < nu2_count; iter++) {
									if (bonly_serial) {

										bool bflag56 = false;
										if (my_amg_manager.iFinnest_ilu == 1) {
											if (my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6) {
												doublerealT dn = 1.0*n_a[i_id_level_local];
												doublerealT dnnz = 1.0*nnz_a[i_id_level_local];
												if (dnnz / dn <= dapply_ilu_max_pattern_size) {
													bflag56 = true;
												}
											}
										}

										if (bILU2smoother == 1) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												milu0[i_id_level_local].zbuf[i43 + 1] = residual_fine[i_id_level_local][i43 + 1];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu0[i_id_level_local].zbuf, milu0[i_id_level_local].zbuf2, milu0[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												error_approx_coarse[i_id_level_local - 1][i43 + 1] += milu0[i_id_level_local].zbuf2[i43 + 1];
											}
										}
										else if (1 && ((bILU2smoother == 2) || bflag56)) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												milu2[i_id_level_local].zbuf[i43 + 1] = residual_fine[i_id_level_local][i43 + 1];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu2[i_id_level_local].zbuf, milu2[i_id_level_local].zbuf2, milu2[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												error_approx_coarse[i_id_level_local - 1][i43 + 1] += milu2[i_id_level_local].zbuf2[i43 + 1];
											}
										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1);
										}
									}
									else {
										seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1);
									}
								}
							}


						}

						// prolongation
						// residual_r
						//doublerealT *error_approx_fine[i_id_level_local - 1] = new doublerealT[n_a[i_id_level_local - 1] + 1];
#pragma omp parallel for
						for (integer ii = 1; ii <= n_a[i_id_level_local - 1]; ii++) {
							error_approx_fine[i_id_level_local - 1][ii] = 0.0;
						}

						integer innz_aRP_loc = 1;
						for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
							innz_aRP_loc += nnz_aRP[i_72];
						}
						prolongation(P, innz_aRP_loc, innz_aRP_loc - 1 + nnz_aRP[i_id_level_local - 1],  error_approx_fine[i_id_level_local - 1], error_approx_coarse[i_id_level_local - 1], n_a[i_id_level_local - 1]);



						// correction
#pragma omp parallel for
						for (integer ii = 1; ii <= n_a[i_id_level_local - 1]; ii++) {
							error_approx_coarse[i_id_level_local - 2][ii] += error_approx_fine[i_id_level_local - 1][ii];
						}

						// free
						//delete[] error_approx_fine[i_id_level_local - 1];
						//delete[] error_approx_coarse[i_id_level_local - 1];
						//delete[] residual_coarse[i_id_level_local - 1];
						//delete[] residual_fine[i_id_level_local - 1];

					} // 21
				}






				// post smothing
				if (process_flow_logic) {


					// smother
					integer iter = 0;
					for (iter = 0; iter < nu2; iter++) {
						//quick seidel
						if (bonly_serial) {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1);
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1);
						}

						//residualq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0] , residual_fine[1]);
						residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);


						Rnext_1 = norma(residual_fine[1], n_a[1]);
						// this is process flow logic
						if (Rnext_1 < process_flow_alpha*R0_1) {
							// Смысл модификации в том что мы экономим итерации на пресмутере.
							break; // досрочно опускаемся на следующий уровень если он есть конечно.
						}
						else {
							Rprev_1 = Rnext_1;
						}
					}

					if (iter == nu2) {
						printf("level 1 limit postsmother iteration is reached\n");
					}

				}
				else {
					for (integer iter = 0; iter < nu2; iter++) {
						if (bonly_serial) {
							if (bILU2smoother == 1) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									milu0[1].zbuf[i43 + 1] = residual_fine[1][i43 + 1];
								}
								lusol_1patchforRUMBA(n_a[1], milu0[1].zbuf, milu0[1].zbuf2, milu0[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									error_approx_coarse[0][i43 + 1] += milu0[1].zbuf2[i43 + 1];
								}
							}
							else if (1 && bILU2smoother == 2) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									milu2[1].zbuf[i43 + 1] = residual_fine[1][i43 + 1];
								}
								lusol_1patchforRUMBA(n_a[1], milu2[1].zbuf, milu2[1].zbuf2, milu2[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									error_approx_coarse[0][i43 + 1] += milu2[1].zbuf2[i43 + 1];
								}
							}
							else {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1);
							}
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1);
						}
					}
				}

				delete[] R0_21;
				delete[] Rprev_21;
				delete[] Rnext_21;

				R0_21 = nullptr;
				Rprev_21 = nullptr;
				Rnext_21 = nullptr;

			}

			move_up(nu1, nu2);

			// prolongation
			// residual_r
			//doublerealT *error_approx_fine[0] = new doublerealT[n_a[0] + 1];
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[0]; ii++) {
				error_approx_fine[0][ii] = 0.0;
			}

			prolongation(P, 1, nnz_aRP[0],  error_approx_fine[0], error_approx_coarse[0], n_a[0]);

			// correction
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[0]; ii++) {
				//if (row_ptr_start[ii] != row_ptr_end[ii]) {
					// Не условие Дирихле.
				z76[ii] += error_approx_fine[0][ii];
				//	}
			}

			// free
			//delete[] error_approx_fine[0];
			//delete[] error_approx_coarse[0];
			//delete[] residual_coarse[0];
			//delete[] residual_fine[0];
		}



		//doublerealT R0_0 = 0.0;
		//doublerealT Rprev_0 = 0.0, Rnext_0 = 0.0;
		// post smothing
		if (process_flow_logic) {
			// calculate initial residual.
			//residualq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
			Rprev_0 = norma(residual_fine[0], n_a[0]);

			// smother
			integer iter = 0;
			for (iter = 0; iter < nFinestSweeps; iter++) {
				//quick seidel
				if (bonly_serial) {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1);
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 1);
				}
				//residualq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
				Rnext_0 = norma(residual_fine[0], n_a[0]);
				// this is process flow logic
				if (Rnext_0 < process_flow_alpha*Rprev_0) {
					// Смысл модификации в том что мы экономим итерации на пресмутере.
					break; // досрочно опускаемся на следующий уровень если он есть конечно.
				}
				else {
					Rprev_0 = Rnext_0;
				}
			}
			if (iter == nFinestSweeps) {
				printf("level 0 limit postsmother iteration is reached\n");
			}

		}
		else {
			// nFinnestSweeps new logic 14 jan 2016.
			// smother
			for (integer iter = 0; iter < nFinestSweeps; iter++) {
				//seidel<doublereal>(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0]);
				//quick seidel
				if (bonly_serial) {
					if (bILU2smoother == 1) {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							milu0[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
						}
						lusol_1patchforRUMBA(n_a[0], milu0[0].zbuf, milu0[0].zbuf2, milu0[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							z76[i43 + 1] += milu0[0].zbuf2[i43 + 1];
						}
					}
					else if ((bILU2smoother == 2) || (my_amg_manager.iFinnest_ilu == 1)) {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							milu2[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
						}
						lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							z76[i43 + 1] += milu2[0].zbuf2[i43 + 1];
						}
					}
					else {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1);
					}
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 1);
				}
			}
		}

	}
} // V_cycle_solve

// 08.01.2018 Перенесено в отдельную функцию, т.к. используется неоднократно. Передаётся большое число параметров.
//  A*z76=s76; 891 строка.
// Для плохих строк при их наличии мы домножаем правую часть на минус 1.0
template <typename doublerealT>
void V_cycle_solve(Ak2& Amat, doublereal*& z76, doublereal*& s76, bool process_flow_logic, integer*& row_ptr_start,
	integer*& row_ptr_end, doublerealT**& residual_fine, doublerealT**& diag, 
	doublerealT**& diag_minus_one, integer* n_a, bool bonly_serial,
	doublerealT process_flow_beta, bool*& F_false_C_true, integer& nu1, integer& nu2, integer bILU2smoother,
	integer ilevel, integer inumberVcyclelocbicgstab, integer imyinit, const integer idim_diag,
	LEVEL_ADDITIONAL_DATA*& milu2, LEVEL_ADDITIONAL_DATA0* milu0, bool**& nested_desection,
	Ak1*& P, // prolongation он же restriction (метод не чувствителен к сортировке).
	integer* nnz_aRP, doublerealT**& residual_coarse, integer igam, integer* nnz_a,
	doublerealT**& error_approx_coarse, doublerealT dapply_ilu_max_pattern_size,
	doublerealT process_flow_alpha, doublerealT**& error_approx_fine,
	integer nFinestSweeps) {

	// Один V - цикл алгебраического многосеточного метода.
	// A*z76=s76;

	const integer ZERO_INIT = 0;
	const integer RANDOM_INIT = 1;// надо увеличивать nu1, nu2 с 1,2 до 5 наверно.

	for (integer i_13 = 0; i_13 < inumberVcyclelocbicgstab; i_13++)
	{


		doublerealT R0_0 = 0.0;
		doublerealT Rprev_0 = 0.0, Rnext_0 = 0.0;
		if (process_flow_logic) {
			// calculate initial residual.
			//residualq(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
			R0_0 = norma(residual_fine[0], n_a[0]);
			Rprev_0 = R0_0;

			// smother
			integer iter = 0;
			for (iter = 0; iter < nu1; iter++) {
				//quick seidel
				if (bonly_serial) {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0],0);
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0]);
				}
				//residualq(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
				Rnext_0 = norma(residual_fine[0], n_a[0]);
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
				printf("level 0 limit presmother iteration is reached\n");
			}

		}
		else {
			// smoother
			for (integer iter = 0; iter < nu1; iter++) {
				//seidel(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0]);
				//quick seidel
				if (bonly_serial) {
					if (bILU2smoother == 1) {
						// ILU0
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0],0);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							milu0[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
						}
						lusol_1patchforRUMBA(n_a[0], milu0[0].zbuf, milu0[0].zbuf2, milu0[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							z76[i43 + 1] += milu0[0].zbuf2[i43 + 1];
						}
					}
					else if ((bILU2smoother == 2) || (my_amg_manager.iFinnest_ilu == 1)) {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0],0);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							milu2[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
						}
						lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							z76[i43 + 1] += milu2[0].zbuf2[i43 + 1];
						}
					}
					else {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0],0);
					}
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 0, diag_minus_one[0]);
				}
			}
		}

		//exporttecplot(x, n);

		move_down(nu1, nu2);

		if (!process_flow_logic) {
			// residual_r
			//doublerealT *residual_fine[0] = new doublerealT[n_a[0] + 1];
			//residual<doublereal>(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0], residual_fine[0]);
			//residualq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
		}



		//if (iprohod%5==0) getchar();
		if (ilevel > 1) {

			//doublerealT *residual_coarse = new doublerealT[n_a[1] + 1];
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[1]; ii++) {
				residual_coarse[0][ii] = 0.0;
			}

			// restriction
			restriction(P, 1, nnz_aRP[0], residual_fine[0], residual_coarse[0], n_a[1]);


			// Amat*e=r;
			//doublerealT* error_approx_coarse = new doublerealT[n_a[1] + 1];
			if ((imyinit == ZERO_INIT)) {
#pragma omp parallel for
				for (integer ii = 1; ii <= n_a[1]; ii++) {
					error_approx_coarse[0][ii] = 0.0;
				}
			}

			if ((imyinit == RANDOM_INIT)) {
				// (1,110); (0.8, 37); (0.7, 29); (0.6, 25); (0.5, 20); (0.4, 17); (0.3, 18); (0.0, 19);
#pragma omp parallel for
				for (integer ii = 1; ii <= n_a[1]; ii++) {
					error_approx_coarse[0][ii] = 0.4 * fabs(residual_coarse[0][ii]) * (rand() / ((doublerealT)RAND_MAX));
				}
			}


			for (integer i_37 = 1; i_37 <= igam; i_37++) {
				doublerealT R0_1 = 0.0;
				doublerealT Rprev_1 = 0.0, Rnext_1 = 0.0;
				if (process_flow_logic) {
					// calculate initial residual.
					//residualq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1]);
					residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);

					R0_1 = norma(residual_fine[1], n_a[1]);
					Rprev_1 = R0_1;

					// smother
					integer iter = 0;
					for (iter = 0; iter < nu1; iter++) {
						//quick seidel
						if (bonly_serial) {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1],1);
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1]);
						}

						//residualq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1]);
						residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);

						Rnext_1 = norma(residual_fine[1], n_a[1]);
						// this is process flow logic
						if (Rnext_1 > process_flow_beta * Rprev_1) {
							// Смысл модификации в том что мы экономим итерации на пресмутере.
							break; // досрочно опускаемся на следующий уровень если он есть конечно.
						}
						else {
							Rprev_1 = Rnext_1;
						}
					}

					if (iter == nu1) {
						printf("level 1 limit presmother iteration is reached\n");
					}

				}
				else {

					// pre smothing
					for (integer iter = 0; iter < nu1; iter++) {
						//seidel(Amat, 1 + 2 * nnz_a[0], 2 * nnz_a[0] + nnz_a[1], error_approx_coarse, residual_coarse, flag, n_a[1]);
						if (bonly_serial) {
							if (bILU2smoother == 1) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1],1);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									milu0[1].zbuf[i43 + 1] = residual_fine[1][i43 + 1];
								}
								lusol_1patchforRUMBA(n_a[1], milu0[1].zbuf, milu0[1].zbuf2, milu0[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									error_approx_coarse[0][i43 + 1] += milu0[1].zbuf2[i43 + 1];
								}


							}
							else if (1 && bILU2smoother == 2) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1],1);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									milu2[1].zbuf[i43 + 1] = residual_fine[1][i43 + 1];
								}
								lusol_1patchforRUMBA(n_a[1], milu2[1].zbuf, milu2[1].zbuf2, milu2[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									error_approx_coarse[0][i43 + 1] += milu2[1].zbuf2[i43 + 1];
								}


							}
							else {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1],1);
							}
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 0, diag_minus_one[1]);
						}
					}
				}





				doublerealT* R0_21 = new doublerealT[idim_diag];
				doublerealT* Rprev_21 = new doublerealT[idim_diag];
				doublerealT* Rnext_21 = new doublerealT[idim_diag];

				//3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
				for (integer i_id_level_local = 2; i_id_level_local < idim_diag; i_id_level_local++) {

					move_down(nu1, nu2);


					if (ilevel > i_id_level_local) {
						// Впервые создан универсальный програмный код на произвольное количество уровней вложенности 
						// в solution phase.
						// Нужно только написать универсальный обработчик для R0_21.
						// 4 декабря 2016. 

						// residual
						if (!process_flow_logic) {
							integer in_a_loc = 0;
							for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
								// sum(n_a[0]+n_a[19];
								in_a_loc += n_a[i_72];
							}
							//doublerealT *residual_fine[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local-1] + 1];
							//residual(Amat, 1 + nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6]+ nnz_a[7]+ nnz_a[8]+ nnz_a[9]+ nnz_a[10]+ nnz_a[11]+ nnz_a[12]+ nnz_a[13]+ nnz_a[14]+ nnz_a[15]+ nnz_a[16]+nnz_a[17]+nnz_a[18]+nnz_a[19], nnz_a[0] + nnz_a[1] + nnz_a[2] + nnz_a[3] + nnz_a[4] + nnz_a[5] + nnz_a[6] + nnz_a[7]+ nnz_a[8]+ nnz_a[9]+ nnz_a[10]+ nnz_a[11]+ nnz_a[12]+ nnz_a[13]+ nnz_a[14]+ nnz_a[15]+ nnz_a[16]+ nnz_a[17]+nnz_a[18]+nnz_a[19]+nnz_a[20], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], flag, n_a[i_id_level_local - 1], residual_fine[i_id_level_local - 1]);
							//residualq(Amat, 1, n_a[i_id_level_local-1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end,  in_a_loc , residual_fine[i_id_level_local - 1]);
							//residualq2(Amat, 1, n_a[i_id_level_local-1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end, in_a_loc, residual_fine[i_id_level_local - 1], diag[i_id_level_local - 1]);
							residualq2(Amat, 1, n_a[i_id_level_local - 1], error_approx_coarse[i_id_level_local - 2], residual_coarse[i_id_level_local - 2], row_ptr_start, row_ptr_end, in_a_loc, residual_fine[i_id_level_local - 1], diag[i_id_level_local - 1], diag_minus_one[i_id_level_local - 1]);

						}

						//doublerealT *residual_coarse[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local] + 1];

						// restriction

						integer innz_aRP_loc = 1;
						for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
							innz_aRP_loc += nnz_aRP[i_72];
						}
						//restriction(P, 1 + nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16] + nnz_aRP[17] + nnz_aRP[18] + nnz_aRP[19], nnz_aRP[0] + nnz_aRP[1] + nnz_aRP[2] + nnz_aRP[3] + nnz_aRP[4] + nnz_aRP[5] + nnz_aRP[6] + nnz_aRP[7] + nnz_aRP[8] + nnz_aRP[9] + nnz_aRP[10] + nnz_aRP[11] + nnz_aRP[12] + nnz_aRP[13] + nnz_aRP[14] + nnz_aRP[15] + nnz_aRP[16] + nnz_aRP[17] + nnz_aRP[18] + nnz_aRP[19] + nnz_aRP[20], flag, residual_fine[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], n_a[i_id_level_local - 1], n_a[i_id_level_local]);
						restriction(P, innz_aRP_loc, innz_aRP_loc - 1 + nnz_aRP[i_id_level_local - 1], residual_fine[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], n_a[i_id_level_local]);



						// Amat*e=r;
						//doublerealT* error_approx_coarse[i_id_level_local-1] = new doublerealT[n_a[i_id_level_local] + 1];

						if (0) {
							//if (icount_V_cycle == 1) {
							//	for (integer ii = 1; ii <= n_a[i_id_level_local]; ii++) {
							// обнуление только на первом проходе.
							//		error_approx_coarse[i_id_level_local - 1][ii] = 0.0;
							//error_approx_coarse[i_id_level_local - 1][ii] = (rand() / ((doublerealT)RAND_MAX));
							//	}
							//}
							// На задаче BSK_Dmitrii random оптимизация оказалась гораздо хуже, но сходимость тоже присутствовала.
#pragma omp parallel for
							for (integer ii = 1; ii <= n_a[i_id_level_local]; ii++) {
								// 0.4*fabs(residual_coarse[i_id_level_local - 1][ii]) - амплитуда.
								// 0.4 - демпфирующий множитель.
								// diag[i_id_level_local]
								// (1,110); (0.8, 37); (0.7, 29); (0.6, 25); (0.5, 20); (0.4, 17); (0.3, 18); (0.0, 19);
								error_approx_coarse[i_id_level_local - 1][ii] = 0.4 * fabs(residual_coarse[i_id_level_local - 1][ii]) * (rand() / ((doublerealT)RAND_MAX));
							}
						}
						else {
							if (imyinit == ZERO_INIT) {
#pragma omp parallel for
								for (integer ii = 1; ii <= n_a[i_id_level_local]; ii++) {
									error_approx_coarse[i_id_level_local - 1][ii] = 0.0;
								}
							}
							if (imyinit == RANDOM_INIT) {
#pragma omp parallel for
								for (integer ii = 1; ii <= n_a[i_id_level_local]; ii++) {
									// (1,110); (0.8, 37); (0.7, 29); (0.6, 25); (0.5, 20); (0.4, 17); (0.3, 18); (0.0, 19);
									error_approx_coarse[i_id_level_local - 1][ii] = 0.4 * fabs(residual_coarse[i_id_level_local - 1][ii]) * (rand() / ((doublerealT)RAND_MAX));
								}
							}
						}


						//for (integer i_37 = 1; i_37 <= igam; i_37++)
						{
							// pre smothing
							//doublerealT R0_21 = 0.0;
							//doublerealT Rprev_21 = 0.0, Rnext_21 = 0.0;
							R0_21[i_id_level_local] = 0.0;
							Rprev_21[i_id_level_local] = 0.0;
							Rnext_21[i_id_level_local] = 0.0;
							if (process_flow_logic) {

								integer in_a_loc = 0;
								for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
									// sum(n_a[0]+n_a[19];
									in_a_loc += n_a[i_72];
								}

								// calculate initial residual.
								//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end,in_a_loc + n_a[i_id_level_local-1] , residual_fine[i_id_level_local]);
								residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local - 1]);



								R0_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
								Rprev_21[i_id_level_local] = R0_21[i_id_level_local];



								// smother
								integer iter = 0;
								integer nu1_count = nu1;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu1_count = 100;

									printf("Direct method is start.../n");

									IMatrix sparseS; // разреженная матрица в формате IMatrix
									initIMatrix(&sparseS, n_a[i_id_level_local]);

									integer istartpos = 1 + in_a_loc + n_a[i_id_level_local - 1];
									integer iendpos = n_a[i_id_level_local] + in_a_loc + n_a[i_id_level_local - 1];

									for (integer i56 = istartpos; i56 <= iendpos; i56++) {

										//setValueIMatrix(&sparseS, Amat.i[row_ptr_start[i56]] - 1, Amat.i[row_ptr_start[i56]] - 1, 1.0 / Amat.aij[row_ptr_start[i56]]);
										setValueIMatrix(&sparseS, i56 - istartpos, i56 - istartpos, 1.0 / Amat.aij[row_ptr_start[i56]]);//3.02.2019
										//if (Amat.i[row_ptr_start[i56]] != i56 - istartpos + 1) {
											//printf("Amat.i[row_ptr_start[i56]]=%lld i56 - istartpos + 1=%lld\n", Amat.i[row_ptr_start[i56]],i56 - istartpos + 1);
											//getchar();
										//}
										const doublerealT nonzeroEPS = 1e-37; // для отделения вещественного нуля

										integer is15 = row_ptr_start[i56] + 1;
										integer is25 = row_ptr_end[i56];

										for (integer ii17 = is15; ii17 <= is25; ii17++)
										{

											if ((fabs(Amat.aij[ii17]) > nonzeroEPS)) {
												//setValueIMatrix(&sparseS, Amat.i[ii17] - 1, Amat.j[ii17] - 1, Amat.aij[ii17]);
												setValueIMatrix(&sparseS, ii17 - is15, Amat.j[ii17] - 1, Amat.aij[ii17]);
											}

										}
									}

									doublereal* dX025 = new doublereal[n_a[i_id_level_local]];
									doublereal* dV25 = new doublereal[n_a[i_id_level_local]];
#pragma omp parallel for
									for (integer i57 = 0; i57 < n_a[i_id_level_local]; i57++) {
										dX025[i57] = error_approx_coarse[i_id_level_local - 1][i57 + 1];
										dV25[i57] = diag_minus_one[i_id_level_local - 1][i57 + 1] *residual_coarse[i_id_level_local - 1][i57 + 1];
									}

									// главный метод, возвращающий решение x,
									// принимает вектор свободных членов b и 
									// квадратную матрицу xO в специальном разреженном формате.
									// реализация без барьера и итерационного уточнения.
									calculateSPARSEgaussArray(&sparseS, dX025, dV25);
#pragma omp parallel for
									for (integer i57 = 0; i57 < n_a[i_id_level_local]; i57++) {
										error_approx_coarse[i_id_level_local - 1][i57 + 1] = dX025[i57];
									}

									delete[] dX025;
									delete[] dV25;


									freeIMatrix(&sparseS);

									//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);
									residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);


									Rnext_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
									// this is process flow logic
									if (Rnext_21[i_id_level_local] > process_flow_beta * Rprev_21[i_id_level_local]) {
										// Смысл модификации в том что мы экономим итерации на пресмутере.
										break; // досрочно опускаемся на следующий уровень если он есть конечно.
									}
									else {
										Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
									}
								}
								else {
									for (iter = 0; iter < nu1_count; iter++) {
										//quick seidel
										if (bonly_serial) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local], i_id_level_local);
										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local]);
										}

										//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local-1], residual_coarse[i_id_level_local-1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);
										residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);


										Rnext_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
										// this is process flow logic
										if (Rnext_21[i_id_level_local] > process_flow_beta * Rprev_21[i_id_level_local]) {
											// Смысл модификации в том что мы экономим итерации на пресмутере.
											break; // досрочно опускаемся на следующий уровень если он есть конечно.
										}
										else {
											Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
										}
									}
								}

								if (iter == nu1) {
#if doubleintprecision == 1
									printf("level %lld limit presmother iteration is reached\n", i_id_level_local);
#else
									printf("level %d limit presmother iteration is reached\n", i_id_level_local);
#endif

								}

							}
							else {
								integer nu1_count = nu1;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu1_count = 100;
								}
								for (integer iter = 0; iter < nu1_count; iter++) {

									integer in_a_loc = 0;
									for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
										// sum(n_a[0]+n_a[19];
										in_a_loc += n_a[i_72];
									}

									if (bonly_serial) {

										bool bflag56 = false;
										if (my_amg_manager.iFinnest_ilu == 1) {
											if (my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6) {
												doublerealT dn = 1.0 * n_a[i_id_level_local];
												doublerealT dnnz = 1.0 * nnz_a[i_id_level_local];
												if (dnnz / dn <= dapply_ilu_max_pattern_size) {
													bflag56 = true;
												}
											}
										}

										if (bILU2smoother == 1) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local], i_id_level_local);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												milu0[i_id_level_local].zbuf[i43 + 1] = residual_fine[i_id_level_local][i43 + 1];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu0[i_id_level_local].zbuf, milu0[i_id_level_local].zbuf2, milu0[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												error_approx_coarse[i_id_level_local - 1][i43 + 1] += milu0[i_id_level_local].zbuf2[i43 + 1];
											}

										}
										else if (1 && ((bILU2smoother == 2) || bflag56)) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local], i_id_level_local);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												milu2[i_id_level_local].zbuf[i43 + 1] = residual_fine[i_id_level_local][i43 + 1];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu2[i_id_level_local].zbuf, milu2[i_id_level_local].zbuf2, milu2[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												error_approx_coarse[i_id_level_local - 1][i43 + 1] += milu2[i_id_level_local].zbuf2[i43 + 1];
											}

										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local], i_id_level_local);
										}
									}
									else {
										seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 0, diag_minus_one[i_id_level_local]);
									}
								}
							}
						}

					}

				}







				//center
				// ЭТО сердцевина SOLUTION PHASE.

				// TODO нижний 3.12.2016 (осталось один нижний и один верхний).


				// 21
				for (integer i_id_level_local = idim_diag - 1; i_id_level_local >= 2; i_id_level_local--) {

					if (ilevel > i_id_level_local) {

						{

							integer in_a_loc = 0;
							for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
								// sum(n_a[0]+n_a[19];
								in_a_loc += n_a[i_72];
							}

							// post smoothing
							// doublerealT R0_20 = 0.0;
							///doublerealT Rprev_20 = 0.0, Rnext_20 = 0.0;
							if (process_flow_logic) {
								// calculate initial residual.
								//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc  + n_a[i_id_level_local - 1], residual_fine[i_id_level_local]);



								//Rprev_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);


								// smother
								integer iter = 0;
								integer nu2_count = nu2;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu2_count = 100;																												
								}
								else {
									for (iter = 0; iter < nu2_count; iter++) {
										//quick seidel
										if (bonly_serial) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local], i_id_level_local);
										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local]);
										}

										//residualq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local-1], residual_fine[i_id_level_local]);
										residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);


										Rnext_21[i_id_level_local] = norma(residual_fine[i_id_level_local], n_a[i_id_level_local]);
										// this is process flow logic
										if (Rnext_21[i_id_level_local] < process_flow_alpha * R0_21[i_id_level_local]) {
											// Смысл модификации в том что мы экономим итерации на пресмутере.
											break; // досрочно опускаемся на следующий уровень если он есть конечно.
										}
										else {
											Rprev_21[i_id_level_local] = Rnext_21[i_id_level_local];
										}
									}
								}

								if (iter == nu2) {
#if doubleintprecision == 1
									printf("level %lld limit postsmother iteration is reached\n", i_id_level_local);
#else
									printf("level %d limit postsmother iteration is reached\n", i_id_level_local);
#endif

								}

							}
							else {
								integer nu2_count = nu2;
								if (i_id_level_local == ilevel - 1) {
									// на самом грубом уровне нам требуется точное решение.
									//nu2_count = 100;																																					
								}
								for (integer iter = 0; iter < nu2_count; iter++) {
									if (bonly_serial) {

										bool bflag56 = false;
										if (my_amg_manager.iFinnest_ilu == 1) {
											if (my_amg_manager.b_ilu_smoothers_in_nnz_n_LE_6) {
												doublerealT dn = 1.0 * n_a[i_id_level_local];
												doublerealT dnnz = 1.0 * nnz_a[i_id_level_local];
												if (dnnz / dn <= dapply_ilu_max_pattern_size) {
													bflag56 = true;
												}
											}
										}

										if (bILU2smoother == 1) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local], i_id_level_local);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												milu0[i_id_level_local].zbuf[i43 + 1] = residual_fine[i_id_level_local][i43 + 1];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu0[i_id_level_local].zbuf, milu0[i_id_level_local].zbuf2, milu0[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												error_approx_coarse[i_id_level_local - 1][i43 + 1] += milu0[i_id_level_local].zbuf2[i43 + 1];
											}
										}
										else if (1 && ((bILU2smoother == 2) || bflag56)) {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local], i_id_level_local);
											residualq2(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], residual_fine[i_id_level_local], diag[i_id_level_local], diag_minus_one[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												milu2[i_id_level_local].zbuf[i43 + 1] = residual_fine[i_id_level_local][i43 + 1];
											}
											lusol_1patchforRUMBA(n_a[i_id_level_local], milu2[i_id_level_local].zbuf, milu2[i_id_level_local].zbuf2, milu2[i_id_level_local]);
#pragma omp parallel for
											for (integer i43 = 0; i43 < n_a[i_id_level_local]; i43++) {
												error_approx_coarse[i_id_level_local - 1][i43 + 1] += milu2[i_id_level_local].zbuf2[i43 + 1];
											}
										}
										else {
											seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local], i_id_level_local);
										}
									}
									else {
										seidelq(Amat, 1, n_a[i_id_level_local], error_approx_coarse[i_id_level_local - 1], residual_coarse[i_id_level_local - 1], nested_desection[i_id_level_local], row_ptr_start, row_ptr_end, in_a_loc + n_a[i_id_level_local - 1], F_false_C_true, 1, diag_minus_one[i_id_level_local]);
									}
								}
							}


						}

						// prolongation
						// residual_r
						//doublerealT *error_approx_fine[i_id_level_local - 1] = new doublerealT[n_a[i_id_level_local - 1] + 1];
#pragma omp parallel for
						for (integer ii = 1; ii <= n_a[i_id_level_local - 1]; ii++) {
							error_approx_fine[i_id_level_local - 1][ii] = 0.0;
						}

						integer innz_aRP_loc = 1;
						for (integer i_72 = 0; i_72 < i_id_level_local - 1; i_72++) {
							innz_aRP_loc += nnz_aRP[i_72];
						}
						prolongation(P, innz_aRP_loc, innz_aRP_loc - 1 + nnz_aRP[i_id_level_local - 1], error_approx_fine[i_id_level_local - 1], error_approx_coarse[i_id_level_local - 1], n_a[i_id_level_local - 1]);



						// correction
#pragma omp parallel for
						for (integer ii = 1; ii <= n_a[i_id_level_local - 1]; ii++) {
							error_approx_coarse[i_id_level_local - 2][ii] += error_approx_fine[i_id_level_local - 1][ii];
						}

						// free
						//delete[] error_approx_fine[i_id_level_local - 1];
						//delete[] error_approx_coarse[i_id_level_local - 1];
						//delete[] residual_coarse[i_id_level_local - 1];
						//delete[] residual_fine[i_id_level_local - 1];

					} // 21
				}






				// post smothing
				if (process_flow_logic) {


					// smother
					integer iter = 0;
					for (iter = 0; iter < nu2; iter++) {
						//quick seidel
						if (bonly_serial) {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1],1);
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1]);
						}

						//residualq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0] , residual_fine[1]);
						residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);


						Rnext_1 = norma(residual_fine[1], n_a[1]);
						// this is process flow logic
						if (Rnext_1 < process_flow_alpha * R0_1) {
							// Смысл модификации в том что мы экономим итерации на пресмутере.
							break; // досрочно опускаемся на следующий уровень если он есть конечно.
						}
						else {
							Rprev_1 = Rnext_1;
						}
					}

					if (iter == nu2) {
						printf("level 1 limit postsmother iteration is reached\n");
					}

				}
				else {
					for (integer iter = 0; iter < nu2; iter++) {
						if (bonly_serial) {
							if (bILU2smoother == 1) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1],1);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									milu0[1].zbuf[i43 + 1] = residual_fine[1][i43 + 1];
								}
								lusol_1patchforRUMBA(n_a[1], milu0[1].zbuf, milu0[1].zbuf2, milu0[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									error_approx_coarse[0][i43 + 1] += milu0[1].zbuf2[i43 + 1];
								}
							}
							else if (1 && bILU2smoother == 2) {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1],1);
								residualq2(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], residual_fine[1], diag[1], diag_minus_one[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									milu2[1].zbuf[i43 + 1] = residual_fine[1][i43 + 1];
								}
								lusol_1patchforRUMBA(n_a[1], milu2[1].zbuf, milu2[1].zbuf2, milu2[1]);
#pragma omp parallel for
								for (integer i43 = 0; i43 < n_a[1]; i43++) {
									error_approx_coarse[0][i43 + 1] += milu2[1].zbuf2[i43 + 1];
								}
							}
							else {
								seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1],1);
							}
						}
						else {
							seidelq(Amat, 1, n_a[1], error_approx_coarse[0], residual_coarse[0], nested_desection[1], row_ptr_start, row_ptr_end, n_a[0], F_false_C_true, 1, diag_minus_one[1]);
						}
					}
				}

				delete[] R0_21;
				delete[] Rprev_21;
				delete[] Rnext_21;

				R0_21 = nullptr;
				Rprev_21 = nullptr;
				Rnext_21 = nullptr;

			}

			move_up(nu1, nu2);

			// prolongation
			// residual_r
			//doublerealT *error_approx_fine[0] = new doublerealT[n_a[0] + 1];
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[0]; ii++) {
				error_approx_fine[0][ii] = 0.0;
			}

			prolongation(P, 1, nnz_aRP[0], error_approx_fine[0], error_approx_coarse[0], n_a[0]);

			// correction
#pragma omp parallel for
			for (integer ii = 1; ii <= n_a[0]; ii++) {
				//if (row_ptr_start[ii] != row_ptr_end[ii]) {
					// Не условие Дирихле.
				z76[ii] += error_approx_fine[0][ii];
				//	}
			}

			// free
			//delete[] error_approx_fine[0];
			//delete[] error_approx_coarse[0];
			//delete[] residual_coarse[0];
			//delete[] residual_fine[0];
		}



		//doublerealT R0_0 = 0.0;
		//doublerealT Rprev_0 = 0.0, Rnext_0 = 0.0;
		// post smothing
		if (process_flow_logic) {
			// calculate initial residual.
			//residualq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
			residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
			Rprev_0 = norma(residual_fine[0], n_a[0]);

			// smother
			integer iter = 0;
			for (iter = 0; iter < nFinestSweeps; iter++) {
				//quick seidel
				if (bonly_serial) {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0],0);
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0]);
				}
				//residualq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0]);
				residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
				Rnext_0 = norma(residual_fine[0], n_a[0]);
				// this is process flow logic
				if (Rnext_0 < process_flow_alpha * Rprev_0) {
					// Смысл модификации в том что мы экономим итерации на пресмутере.
					break; // досрочно опускаемся на следующий уровень если он есть конечно.
				}
				else {
					Rprev_0 = Rnext_0;
				}
			}
			if (iter == nFinestSweeps) {
				printf("level 0 limit postsmother iteration is reached\n");
			}

		}
		else {
			// nFinnestSweeps new logic 14 jan 2016.
			// smother
			for (integer iter = 0; iter < nFinestSweeps; iter++) {
				//seidel<doublereal>(Amat, 1, nnz_a[0], z76, s76, flag, n_a[0]);
				//quick seidel
				if (bonly_serial) {
					if (bILU2smoother == 1) {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0],0);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							milu0[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
						}
						lusol_1patchforRUMBA(n_a[0], milu0[0].zbuf, milu0[0].zbuf2, milu0[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							z76[i43 + 1] += milu0[0].zbuf2[i43 + 1];
						}
					}
					else if ((bILU2smoother == 2) || (my_amg_manager.iFinnest_ilu == 1)) {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0],0);
						residualq2(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, residual_fine[0], diag[0], diag_minus_one[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							milu2[0].zbuf[i43 + 1] = residual_fine[0][i43 + 1];
						}
						lusol_1patchforRUMBA(n_a[0], milu2[0].zbuf, milu2[0].zbuf2, milu2[0]);
#pragma omp parallel for
						for (integer i43 = 0; i43 < n_a[0]; i43++) {
							z76[i43 + 1] += milu2[0].zbuf2[i43 + 1];
						}
					}
					else {
						seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0],0);
					}
				}
				else {
					seidelq<doublereal>(Amat, 1, n_a[0], z76, s76, nested_desection[0], row_ptr_start, row_ptr_end, 0, F_false_C_true, 1, diag_minus_one[0]);
				}
			}
		}

	}
} // V_cycle_solve

#endif /*BASIC_FUNCTIONS_MY_AGREGAT_AMG_CPP*/





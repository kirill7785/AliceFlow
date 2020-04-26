// leftist_heap.cpp
// Левосторонняя куча. 10.04.2020 По материалам лекции Андрея Станкевича.
// Куча относится к сливаемым кучам.
// Robert E. Tarjan (1983).

// Операции над leftist heap:
// singlton,
// merge
// insert,
// delete_min, get_min.
// Всё выполняется за двоичный логарифм. Куча сливаемая.

#pragma once
#ifndef MY_LEFTIST_HEAP_CPP
#define MY_LEFTIST_HEAP_CPP 1

typedef struct Tleftist_heap {
	integer key=0; // ключ сортировки.
	integer link=-1; // ссылка на номер массива.
	// Расстояние до ближайшего вакантного места в правом дереве.
	int rank=-1;
	Tleftist_heap* left = nullptr, * right = nullptr;
	//Tleftist_heap* parent=nullptr; // для операции decrease_key
} leftist_heap;

// меняет местами кучи x и y.
void swap(leftist_heap* &x, leftist_heap* &y) {
	//leftist_heap* p = x->parent; x->parent = y->parent; y->parent = p; p = nullptr;
	leftist_heap* z = x; x = y; y = z; z = nullptr;
} // swap

// Создаёт кучу с ключом key.
leftist_heap* singlton(integer key, integer link) {
	leftist_heap* x = new leftist_heap;
	x->rank = 1; // поскольку правое поддерево равно нулю, кратчайший путь к дочернему листу из узла x равен 1
	x->key = key;
	x->link = link;
	x->left = x->right = nullptr;
	//x->parent = nullptr;
	return x;
}

// Слияние двух левосторонних куч.
leftist_heap* merge(leftist_heap* x, leftist_heap* y) {
	if (x == nullptr) return y;
	if (y == nullptr) return x;
	if (x->key > y->key) {
		swap(x, y);
	}
	x->right = merge(x->right,y);// результат не может быть нулевым, так как y ненулевой
	if (x->left == nullptr) {
		swap(x->left, x->right);
		x->rank = 1; // поскольку правое поддерево равно нулю, кратчайший путь к дочернему листу из узла x равен 1
		return x;
	}
	if (x->right->rank > x->left->rank) {
		swap(x->left,x->right);
	}
	x->rank = x->right->rank + 1;
	return x;
} // merge

// вставка узла в левостороннюю кучу.
leftist_heap* insert(leftist_heap* heap, integer key, integer link) {
	leftist_heap* z = singlton(key,link);
	heap = merge(heap, z);
	return heap;
} // insert

// Минимальный элемент в левосторонней куче является корнем. 
// Таким образом, чтобы удалить минимальный элемент, корень удаляется,
// и его поддеревья объединяются для формирования
// новой левосторонней кучи с новым минимумом.
leftist_heap* delete_min(leftist_heap* heap) {
	//heap->left->parent = heap->right->parent = nullptr;
	leftist_heap* z = merge(heap->left,heap->right);
	heap->left = heap->right = nullptr;
	delete heap;
	heap = nullptr;
	return z;
} // delete_min

// посмотреть минимум
integer get_min(leftist_heap* heap) {
	return heap->key;
} // get_min

// посмотреть ссылку на минимум
integer get_min_link(leftist_heap* heap) {
	return heap->link;
} // get_min

/*
Действительно, правильное решение  по реализации функции DecreaseKey 
для левосторонних куч(деревьев), кажется,
скрыто в Интернете. (Моя книга Вайса имеет решение, но в главе об
амортизированном анализе, и я не быстро понимаю аргументацию.)
Во-первых, мы не можем рассматривать DecreaseKey как в двоичных кучах.
Там каждый раз меняются местами узел с его родителем, пока не будет
найден родитель, который на самом деле меньше, чем новое значение. Эта
операция работает в левосторонних деревьях, но может занимать линейное время.
В левых деревьях нет ограничений на длину путей слева, и действительно,
такое дерево может быть только одним длинным линейным крайним левым деревом-списком.
Решение состоит в том, чтобы вырезать уменьшенный узел q (вместе с его 
поддеревом), чтобы изменить значение q на новое значение, и объединить его
с остатком T исходного дерева. Прежде чем мы сможем это сделать, мы должны 
вернуть T в левую форму и восстановить значения длины нулевого пути, которые 
были изменены путем обрезания поддерева. Поскольку мы вырезали q из дерева,
должны быть проверены только узлы от родительского p вверх до корня. 
(К сожалению, это может быть линейное число узлов, как мы уже видели.)
Решение состоит в том, что новые значения длины нулевого пути, которые мы
выравниваем вверх по дереву, всегда на единицу больше, чем значение на
предыдущем уровне. Этот процесс останавливается, когда мы находимся в левом
дочернем элементе, который имеет новое значение ℓ, а его правый брат имеет 
значение r с value≥r. Тогда у родителя уже было значение 
r + 1 = min {ℓ, r} +1, и все готово.
Мы знаем (из базового свойства, которое делает левые деревья эффективными),
что значение нулевой длины пути в ℓ означает, что в этом поддереве есть
2ℓ − 1 узлов, что делает операцию логарифмической.
Тот факт, что значения вверху в новом дереве всегда на единицу больше, 
имеет здесь важное значение, что не было (обязательно) верно по пути в
исходном дереве.
*/
/* Не реализовано 11.04.2020
leftist_heap* decrease_key(leftist_heap* heap, leftist_heap* q, integer new_key) {
	if (new_key >= q->key) return heap; // Увеличивать ключ мы не умеем.
	leftist_heap* z = q->parent;
	if (z->left == q) {
		z->left = nullptr;
	}
	else {
		// z->right==x
		z->right = nullptr;
	}
	// Мы должны вернуть дерево heap в левую форму
	// на пути от z к корню heap.
	// TODO
	q->parent = nullptr;
	q->key = new_key; // уменьшенный ключ.
	heap = merge(heap, q);
	return heap;
} // decrease_key
*/

// Медленная сортировка по времени:
/* CountingSort 21m 6s 160ms
*  TimSort 21min 47s 740ms
*  QuickSort 20min 49s 700ms
*  HeapSort 21min 30s 20ms
*  LeftistHeapSort 49min 27s 730ms
*/
// Сортировка по возрастанию по индексу i
// с помощью leftistheap пирамиды.
void LeftistHeapSort(Ak2& Amat, integer first, integer last) {
	leftist_heap* heap = singlton(Amat.i[first],0);
	for (integer i = first + 1; i <= last; i++) heap = insert(heap, Amat.i[i], i-first);
	doublereal* a = new doublereal[last-first+1];
	integer* j_a = new integer[last - first + 1];
	for (integer i = first; i <= last; i++) {
		a[i - first] = Amat.aij[i];
		j_a[i - first] = Amat.j[i];
	}
	integer i = first;
	while (heap != nullptr) {
		integer i_ = get_min_link(heap);
		Amat.i[i] = get_min(heap);
		Amat.j[i] = j_a[i_];
		Amat.aij[i] = a[i_];
		i++;
		heap = delete_min(heap);
	}
	delete[] a;
	delete[] j_a;
}

// Сортировка по возрастанию по индексу j
// с помощью leftistheap пирамиды.
void LeftistHeapSort_j(Ak2& Amat, integer first, integer last) {
	leftist_heap* heap = singlton(Amat.j[first], 0);
	for (integer i = first + 1; i <= last; i++) heap = insert(heap, Amat.j[i], i - first);
	doublereal* a = new doublereal[last - first + 1];
	integer* i_a = new integer[last - first + 1];
	for (integer i = first; i <= last; i++) {
		a[i - first] = Amat.aij[i];
		i_a[i - first] = Amat.i[i];
	}
	integer i = first;
	while (heap != nullptr) {
		integer i_ = get_min_link(heap);
		Amat.j[i] = get_min(heap);
		Amat.i[i] = i_a[i_];
		Amat.aij[i] = a[i_];
		i++;
		heap = delete_min(heap);
	}
	delete[] a;
	delete[] i_a;
}

// Сортировка по возрастанию по индексу i
// с помощью leftistheap пирамиды.
void LeftistHeapSort(Ak1* &Amat, integer first, integer last) {
	leftist_heap* heap = singlton(Amat[first].i, 0);
	for (integer i = first + 1; i <= last; i++) heap = insert(heap, Amat[i].i, i - first);
	doublereal* a = new doublereal[last - first + 1];
	integer* j_a = new integer[last - first + 1];
	for (integer i = first; i <= last; i++) {
		a[i - first] = Amat[i].aij;
		j_a[i - first] = Amat[i].j;
	}
	integer i = first;
	while (heap != nullptr) {
		integer i_ = get_min_link(heap);
		Amat[i].i = get_min(heap);
		Amat[i].j = j_a[i_];
		Amat[i].aij = a[i_];
		i++;
		heap = delete_min(heap);
	}
	delete[] a;
	delete[] j_a;
}

// Сортировка по возрастанию по индексу j
// с помощью leftistheap пирамиды.
void LeftistHeapSort_j(Ak1* &Amat, integer first, integer last) {
	leftist_heap* heap = singlton(Amat[first].j, 0);
	for (integer i = first + 1; i <= last; i++) heap = insert(heap, Amat[i].j, i - first);
	doublereal* a = new doublereal[last - first + 1];
	integer* i_a = new integer[last - first + 1];
	for (integer i = first; i <= last; i++) {
		a[i - first] = Amat[i].aij;
		i_a[i - first] = Amat[i].i;
	}
	integer i = first;
	while (heap != nullptr) {
		integer i_ = get_min_link(heap);
		Amat[i].j = get_min(heap);
		Amat[i].i = i_a[i_];
		Amat[i].aij = a[i_];
		i++;
		heap = delete_min(heap);
	}
	delete[] a;
	delete[] i_a;
}

#endif // !MY_LEFTIST_HEAP_CPP
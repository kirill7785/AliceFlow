// Очередь по приоритетам на основе сортирующего дерева
// 12.06.2017 (июнь).
// priority_queue.cpp

// Быстро добавляем элементы, быстро читаем максимальный элемент среди
// имеющихся, быстро удаляем максимальный элемент.
// Всё на массиве и очень быстро потому что без рекурсии.
// Время выполнения логарифмическое.

// Закомментированная реализация это неполная 
// (нет операции удаления элемента и поиска по ключу) сильно упрощённая, но 
// правильно отражающая основную идею реализация, а вот рабочий код
// приведённый в конце это полноценный 100% рабочий гибрид этой идеи 
// с быстродействующей хеш-таблицей. 
// Это должен быть очень быстрый код т.к. 1. нет постоянных new и delete -
// память выделяется лишь один раз в момент инициализации и освобождается автоматическим
// деструктором. 2. Рекурсия заменена сверхбыстродействующей итерацией.
// Критика: А где балансировка обеспечивающая логарифмическую длину поиска ?
// Ответ: Балансировка не требуется т.к. по своему определению двоичная куча это полностью заполненное двоичное
// дерево в котором лишь последний уровень содержит позиции для добавления новых элементов. Тем самым
// высота дерева двоичной кучи всегда логарифм по основанию 2 от числа элементов в куче.

/*

// Восходящая установка структуры сортирующего дерева.
// Роберт Седжвик с. 366 в книге 2002 года.
template <class Item>
void fixUp(Item a[], integer k)
{
	while (k > 1 && a[k / 2] < a[k])
	{
		integer kdiv2 = k / 2;
		Item buf = a[k];
		a[k] = a[kdiv2];
		a[kdiv2] = buf;

		k = kdiv2;
	}
}

// Нисходящая установка структуры сортирующего дерева.
template <class Item>
void fixDown(Item a[], integer k, integer N)
{
	while (2 * k <= N)
	{
		integer j = 2 * k;
		if (j < N&&a[j] < a[j + 1]) j++;
		if (!(a[k] < a[j])) break;

		Item buf = a[k];
		a[k] = a[j];
		a[j] = buf;

		k = j;
	}
}

template <class Item>
class PQ
{
private:
	Item *pq;
	integer N;
	integer isize;
public:
	PQ(integer maxN)
	{
		isize = maxN;
		pq = new Item[maxN + 1]; 
		N = 0;
	}
	~PQ()
	{
		if (pq != nullptr) delete[] pq;
	    N = 0;
	}
	integer empty() const
	{
		return N == 0;
	}
	// Вставить элемент в очередь по 
	// приоритетам.
	void insert(Item item)
	{
		if (N + 1 > isize) {
			printf("ERROR!!! priority_queue stack overflow...\n");
			system("pause");
			exit(1);
		}
		else {
			pq[++N] = item;
			fixUp(pq, N);
		}
	}
	// Просто прочитать максимальный элемент.
	Item readmax()
	{
		return pq[1];
	}
	// Возвратить максимальный элемент
	// и удалить его.
	Item getmax()
	{
		Item buf = pq[1];
		pq[1] = pq[N];
		pq[N] = buf;
		fixDown(pq, 1, N - 1);
		return pq[N--];
	}
};

*/

// Используется в алгебраическом многосеточном методе.
// Используются только следующие функции:
// clear, remove, insert, readkeymaxelm.

// Соединяем с быстродействующей хеш таблицей.
template <class Item>
void exch(integer i, integer j, Item* &pq, integer* &qp, integer* &hash) {
	// exchange
	/*
	// begin
	integer t;

	Item buf1 = pq[qp[i]];
	Item buf2 = pq[qp[j]];
	t = qp[i];

	qp[i] = qp[j];
	qp[j] = t;

	pq[qp[i]] = buf1;
	pq[qp[j]] = buf2;

	// end
	*/
	//printf("exchange\n");

	Item t;

	t = pq[j];
	pq[j] = pq[i];
	pq[i] = t;

	integer p;

	p = hash[qp[i]];
	hash[qp[i]] = hash[qp[j]];
	hash[qp[j]] = p;

	p = qp[j];
	qp[j] = qp[i];
	qp[i] = p;

	

}

// Восходящая установка структуры сортирующего дерева.
// Роберт Седжвик с. 366 в книге 2002 года.
template <class Item>
void fixUp(Item* &a, integer* &inda, integer* &hash, integer k)
{
	while (k > 1 && a[k / 2] < a[k])
	{
		integer kdiv2 = k / 2;

		exch(k, kdiv2, a, inda, hash);

		k = kdiv2;
	}
}

// Нисходящая установка структуры сортирующего дерева.
template <class Item>
void fixDown(Item* &a, integer* &inda, integer* &hash, integer k, integer N)
{
	while (2 * k <= N)
	{
		

		integer j = 2 * k;
		if (j < N&&a[j] < a[j + 1]) j++;
		if (!(a[k] < a[j])) break;


		exch(k, j, a, inda, hash);

		k = j;
	}
}

//PQ(integer maxN, integer max_key_size);
//~PQ();
// Есть ли элемент с данным ключом в таблице ?
//bool isfound(integer key);
//bool empty() const;
// Очищаем содержимое и она снова готова к использованию.
//void clear();
// Вернуть элемент с заданным ключом:
// Обязательно предполагается что ключ существует внутри таблицы.
//Item getPQ(integer key);
//Item readmax();
//integer readkeymaxelm();
// Вставить элемент item в очередь по 
// приоритетам если элемент item имеет ключ key.
//template <class Item>
//void insert(Item item, integer key);
// Возвратить максимальный элемент
// и удалить его.
//Item getmax();
// Заменяет элемент с ключом key на элемент val с тем же ключом key.
// При этом ключ key должен быть уникальным.
//void modify(integer key, Item val);
// Удаление элемента с заданным значением ключа.
//void remove(integer key);
// У элемента изменить значение старого ключа на новый ключ
// при этом меняется и само содержимое элемента.
//void change(integer key_serch, integer key_new, integer item_new);


// Ключи должны быть уникальны, целочисленны и различны.
// Двух одинаковых ключей быть не должно, иначе коллизия в хеш-таблице.

template <class Item>
class PQ
{
private:
	// Хранение binary heap.
	Item *pq;
	// Обратный доступ по номеру в qp на ячейку в hash.
	integer *qp; // Ссылка на хеш-таблицу.
	// Доступ по ключу к полю в pq.
	integer *hash; // Хеш-таблица !!!
	integer N;
	integer isize;
	integer ihash_size;

public:
	PQ(integer maxN, integer max_key_size)
	{
		pq=nullptr;
		qp=nullptr;
		hash=nullptr;
		this->isize = maxN;
		if (this->pq != nullptr) delete[] this->pq;
		this->pq = new Item[maxN + 1];
		if (this->qp != nullptr) delete[] this->qp;
		this->qp = new integer[maxN + 1];
		for (integer i_1 = 0; i_1 < maxN + 1; i_1++) {
			// Инициализация: таблица пуста т.к. поле 0
			// в массиве pq никогда не используется.
			this->qp[i_1] = 0;
		}
		this->N = 0;
		// Хеш таблица !!!
		this->ihash_size = max_key_size;
		if (this->hash != nullptr) delete[] this->hash;
		this->hash = new integer[max_key_size+2];
		for (integer i_1 = 0; i_1 < max_key_size + 2; i_1++) {
			// Инициализация: таблица пуста т.к. поле 0
			// в массиве pq никогда не используется.
			this->hash[i_1] = 0;
		}
	}
	~PQ()
	{
		if (this->pq != nullptr) delete[] this->pq;
		this->N = 0;
		if (this->qp != nullptr) delete[] this->qp;
		if (this->hash != nullptr) delete[] this->hash;
	}
	void print_log(char ch) {
		printf("%c\n",ch);
		for (integer i_1 = 1; i_1 <= this->N; i_1++) {
#if doubleintprecision == 1
			printf("[%lld %lld] ", this->pq[i_1], this->qp[i_1]);
#else
			printf("[%d %d] ", this->pq[i_1], this->qp[i_1]);
#endif
		}
		printf("\n");
		system("PAUSE");
	}

	// Меняет местами значения элементов с ключами i и j. 
	// Сохраняет порядок кучи.
	void exchange_speshial_for_Saad(integer i, integer j) {

		//Item t;



		//t = pq[hash[j]];
		//pq[hash[j]] = pq[hash[i]];
		//pq[hash[i]] = t;

		// Этот обмен местами сохраняет порядок кучи.
		Item t1 = this->getPQ(i);
		//Item t2 = getPQ(j);

		//this->remove(i);
		//this->remove(j);
		//this->insert(t1, j);
		//this->insert(t2, i);

		// быстродействущая модификация
		this->remove(j);
		this->insert(t1, j);

	}

	// Очищаем содержимое и она снова готова к использованию.
	void clear()
	{
		for (integer i_1 = 0; i_1 < this->N + 1; i_1++) {
			// Ускоренная очистка хеш таблицы.
			this->hash[this->qp[i_1]] = 0;
		}
		/*
		// Стабильно верный вариант.
		for (integer i_1 = 0; i_1 < this->isize + 1; i_1++) {
			// Инициализация: таблица пуста т.к. поле 0
			// в массиве pq никогда не используется.
			this->qp[i_1] = 0;
		}
		*/
		// 25.02.2018
		// Быстродействие выше т.к. мы просматриваем только 
		// ключи лежащие в очереди а не весь диапазон с запасом.
		for (integer i_1 = 0; i_1 < this->N + 1; i_1++) {
			this->qp[i_1] = 0;
		}
		this->N = 0;
		/*
		for (integer i_1 = 0; i_1 < ihash_size + 2; i_1++) {
			// Инициализация: таблица пуста т.к. поле 0
			// в массиве pq никогда не используется.
			this->hash[i_1] = 0;
		}
		*/
	}
	bool empty() const
	{
		return this->N == 0;
	}
	// Есть ли элемент с данным ключом в таблице ?
	bool isfound(integer key) {
		if (this->hash[key] == 0) {
			// Элемент отсутствует в хеш-таблице.
			return false;
		}
		return true;
	}
	// Вернуть элемент с заданным ключом:
	// Обязательно предполагается что ключ существует внутри таблицы.
	Item getPQ(integer key) {
		if (this->hash[key] == 0) {
			// Элемент отсутствует в хеш таблице.
			printf("priority queue get ERROR: get element not found.\n");
			system("pause");
			exit(1);
		}
		return this->pq[this->hash[key]];
	}
	// Просто прочитать максимальный элемент.
	Item readmax()
	{
		return this->pq[1];
	}
	integer readkeymaxelm() {
		return this->qp[1];
	}

	

	// Вставить элемент item в очередь по 
	// приоритетам если элемент item имеет ключ key.
	void insert(Item item, integer key)
	{
		if (N + 1 > isize) {
			printf("ERROR!!! priority_queue stack overflow...\n");
#if doubleintprecision == 1
			printf("N=%lld\n",N);
#else
			printf("N=%d\n", N);
#endif
			system("pause");
			exit(1);
		}
		else {
			pq[++N] = item;
			hash[key] = N;
			qp[N] = key;
			fixUp(pq, qp, hash, N);
		}
		//print_log('i');
	}
	
	


	// Возвратить максимальный элемент
	// и удалить его.
	Item getmax()
	{
		exch(1, N, pq, qp, hash);

		fixDown(pq, qp, hash, 1, N - 1);
		return pq[N--];
	}

	

	// Заменяет элемент с ключом key на элемент val с тем же ключом key.
	// При этом ключ key должен быть уникальным.
	void modify(integer key, Item val)
	{
		if (hash[key] == 0) {
			// Элемент отсутствует в хеш таблице.
			printf("priority queue modify ERROR: get element not found.\n");
			system("pause");
			exit(1);
		}

		pq[hash[key]] = val;
		// Теперь необходимо восстановить порядок кучи.
		integer i = hash[key];
		fixUp(pq, qp, hash, i);
		fixDown(pq, qp, hash, i, N);
	}


	

	// Удаление элемента с заданным значением ключа.
	void remove(integer key)
	{
		if (N > 0) {
			if (hash[key] == 0) {
				// Элемент отсутствует в хеш-таблице.
				// Ничего не делаем т.к. элемента уже нет.
			}
			else {

				// Удаление.
				if (hash[key] == N) {
					N--;
					hash[key] = 0;
					qp[N + 1] = 0;
					// Ключ исключён из таблицы.
				}
				else {

					integer i = hash[key];

					exch(hash[key], N, pq, qp, hash);

					hash[qp[N]] = 0;
					qp[N] = 0;
					N--;

					// Теперь необходимо восстановить порядок кучи.
					fixUp(pq, qp, hash, i);
					fixDown(pq, qp, hash, i, N);

				}

			}
		}
		//print_log('r');

	}
	
	// У элемента изменить значение старого ключа на новый ключ
	// при этом меняется и само содержимое элемента.
	void change(integer key_search, integer key_new, integer item_new)
	{
		if (hash[key_search] == 0) {
			// Элемент отсутствует в хеш-таблице.
			if (hash[key_new] != 0) {
				// Элемент присутствует в хеш-таблице.
				pq[hash[key_new]] = item_new;
				// Теперь необходимо восстановить порядок кучи.
				fixUp(pq, qp, hash, hash[key_new]);
				fixDown(pq, qp, hash, hash[key_new], N);
			}
			else {
				// Вставка нового ключа с новыми данными.
				insert(item_new, key_new);
			}
		}
		else {
			if (hash[key_new] != 0) {
				// удаление старого ключа со всем его содержимым.
				remove(key_search);
				// Элемент присутствует в хеш-таблице.
				pq[hash[key_new]] = item_new;
				// Теперь необходимо восстановить порядок кучи.
				fixUp(pq, qp, hash[key_new]);
				fixDown(pq, qp, hash[key_new], N);
				
			}
			else {
				// key_new отсутствует.

				hash[key_new] = hash[key_search];
				hash[key_search] = 0; // исключение из дерева.
				pq[hash[key_new]] = item_new;
				qp[hash[key_new]] = key_new;
				// Теперь необходимо восстановить порядок кучи.
				fixUp(pq, qp, hash, hash[key_new]);
				fixDown(pq, qp, hash, hash[key_new], N);

			}
		}
	}

	
};

// Очередь по приоритетам на основе кучи Фибоначчи.
// 21.13.2018
// Все операции O(1) и только операция delete log2(N).

// Поле данных в АВЛ дереве.
class data_BalTree
{
public:
	// --> high priority --> for operation <,>
	//integer ii;
	integer  i, count_neighbour;
	// count_neighbour есть key.
	 // Перегруженные операции сравнения для составного ключа.
	data_BalTree() {
		this->i = -1;
		this->count_neighbour = -1;
	}
	data_BalTree(integer i_tmp, integer count_neighbour_tmp) {
		this->i = i_tmp;
		this->count_neighbour = count_neighbour_tmp;
	}
#if doubleintprecision == 1
	data_BalTree(int i_tmp, int count_neighbour_tmp) {
		this->i = (integer)(i_tmp);
		this->count_neighbour = (integer)(count_neighbour_tmp);
	}
#endif
	bool operator <(const data_BalTree&);
	bool operator >(const data_BalTree&);
	bool operator ==(const data_BalTree&);

};

bool  data_BalTree::operator <(const data_BalTree& b) {
	if (this->count_neighbour < b.count_neighbour) {
		return true;
	}
	else if (this->count_neighbour > b.count_neighbour) {
		return false;
	}
	else if (this->i < b.i) {
		return true;
	}
	else if (this->i > b.i) {
		return false;
	}
	else {
		return false;
	}

}

bool  data_BalTree::operator >(const data_BalTree& b) {
	if (this->count_neighbour < b.count_neighbour) {
		return false;
	}
	else if (this->count_neighbour > b.count_neighbour) {
		return true;
	}
	else if (this->i < b.i) {
		return false;
	}
	else if (this->i > b.i) {
		return true;
	}
	else {
		return false;
	}

}


bool  data_BalTree::operator ==(const data_BalTree& b) {
	if (this->count_neighbour < b.count_neighbour) {
		return false;
	}
	else if (this->count_neighbour > b.count_neighbour) {
		return false;
	}
	else if (this->i < b.i) {
		return false;
	}
	else if (this->i > b.i) {
		return false;
	}
	else {
		return true;
	}

}



/*Copyright (c) 2010, Robin Message <Robin.Message@cl.cam.ac.uk>
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
* Neither the name of the Univsersity of Cambridge nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF CAMBRIDGE OR ROBIN MESSAGE
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* Copyright(c) 2010, Robin Message <Robin.Message@cl.cam.ac.uk>
Все права защищены.
Распространение и использование в виде исходного и двоичного кода, с или без
модификация допускается при соблюдении следующих условий:
*При распространении исходного кода должны сохраняться вышеуказанные авторские права
обратите внимание, этот список условий и следующий отказ от ответственности.
* При повторном распространении двоичного кода должна сохраняться
обратите внимание, этот список условий и следующий отказ от ответственности в
документация и / или другие материалы, поставляемые вместе с дистрибутивом.
* Ни название University Кембриджа, ни
имена его участников могут использоваться для поддержки или продвижения продуктов
производные от данного программного обеспечения без предварительного письменного разрешения.
ЭТО ПРОГРАММНОЕ ОБЕСПЕЧЕНИЕ ПРЕДОСТАВЛЯЕТСЯ ВЛАДЕЛЬЦАМИ АВТОРСКИХ ПРАВ И УЧАСТНИКАМИ "КАК ЕСТЬ" И
ЛЮБЫЕ ЯВНЫЕ ИЛИ ПОДРАЗУМЕВАЕМЫЕ ГАРАНТИИ, ВКЛЮЧАЯ, НО НЕ ОГРАНИЧИВАЯСЬ, ПОДРАЗУМЕВАЕМЫЕ
ГАРАНТИИ ТОВАРНОЙ ПРИГОДНОСТИ И ПРИГОДНОСТИ ДЛЯ ОПРЕДЕЛЕННОЙ ЦЕЛИ
ОТКАЗАВШИЙСЯ.НИ В КОЕМ СЛУЧАЕ КЕМБРИДЖСКИЙ УНИВЕРСИТЕТ
*/
/*
template <class V> class FibonacciHeap;

template <class V> struct node {
private:
	// Указатель на левый сестринский узел.
	node<V>* prev;
	// указатель на правый сестринский узел.
	node<V>* next;
	// указатель на один из дочерних узлов.
	node<V>* child;
	// указатель на родительский узел.
	node<V>* parent;
	V value;

	// количество дочерних узлов.
	int degree;

	
	//логическое значение, которое указывает,
	//были ли потери узлом x дочерних узлов,
	//начиная с момента, когда  x стал дочерним
	//узлом какого-то другого узла.
	//FIBONNACCI_HEAP
	bool marked;
public:
	friend class FibonacciHeap<V>;
	node<V>* getPrev() { return prev; }
	node<V>* getNext() { return next; }
	node<V>* getChild() { return child; }
	node<V>* getParent() { return parent; }
	V getValue() { return value; }
	bool isMarked() { return marked; }

	bool hasChildren() { return child; }
	bool hasParent() { return parent; }
};

template <class V> class FibonacciHeap {
protected:
	node<V>* heap;
public:

	FibonacciHeap() {
		heap = _empty();
	}
	virtual ~FibonacciHeap() {
		if (heap) {
			_deleteAll(heap);
		}
	}
	node<V>* insert(V value) {
		node<V>* ret = _singleton(value);
		heap = _merge(heap, ret);
		return ret;
	}
	void merge(FibonacciHeap& other) {
		heap = _merge(heap, other.heap);
		other.heap = _empty();
	}

	bool isEmpty() {
		return heap == nullptr;
	}

	V getMinimum() {
		return heap->value;
	}

	V removeMinimum() {
		node<V>* old = heap;
		heap = _removeMinimum(heap);
		V ret = old->value;
		delete old;
		old = nullptr;
		return ret;
	}

	void decreaseKey(node<V>* n, V value) {
		heap = _decreaseKey(heap, n, value);
	}

	void deleteKey(V value) {
		node<V>* find_ = find(value);
		if (find_ != nullptr) {
#if doubleintprecision == 1
				decreaseKey(find_, -big_FIBO_integer_Value);
#else
				decreaseKey(find_, -2147483645);
#endif
			removeMinimum();
		}
	}

	node<V>* find(V value) {
		return _find(heap, value);
	}
private:
	node<V>* _empty() {
		return nullptr;
	}

	node<V>* _singleton(V value) {
		node<V>* n = new node<V>;
		n->value = value;
		n->prev = n->next = n;
		n->degree = 0;
		n->marked = false;
		n->child = nullptr;
		n->parent = nullptr;
		return n;
	}

	node<V>* _merge(node<V>* a, node<V>* b) {
		if (a == nullptr)return b;
		if (b == nullptr)return a;
		if (a->value>b->value) {
			node<V>* temp = a;
			a = b;
			b = temp;
		}
		node<V>* an = a->next;
		node<V>* bp = b->prev;
		a->next = b;
		b->prev = a;
		an->prev = bp;
		bp->next = an;
		return a;
	}

	void _deleteAll(node<V>* n) {
		if (n != nullptr) {
			node<V>* c = n;
			do {
				node<V>* d = c;
				c = c->next;
				_deleteAll(d->child);
				delete d;
				d = nullptr;
			} while (c != n);
		}
	}

	void _addChild(node<V>* parent, node<V>* child) {
		child->prev = child->next = child;
		child->parent = parent;
		parent->degree++;
		parent->child = _merge(parent->child, child);
	}

	void _unMarkAndUnParentAll(node<V>* n) {
		if (n == nullptr)return;
		node<V>* c = n;
		do {
			c->marked = false;
			c->parent = nullptr;
			c = c->next;
		} while (c != n);
	}

	node<V>* _removeMinimum(node<V>* n) {
		if (n == nullptr)return n;
		_unMarkAndUnParentAll(n->child);
		if (n->next == n) {
			n = n->child;
		}
		else {
			n->next->prev = n->prev;
			n->prev->next = n->next;
			n = _merge(n->next, n->child);
		}
		if (n == nullptr)return n;
		node<V>* trees[64] = { nullptr };

		while (true) {
			if (trees[n->degree] != nullptr) {
				node<V>* t = trees[n->degree];
				if (t == n)break;
				trees[n->degree] = nullptr;
				if (n->value<t->value) {
					t->prev->next = t->next;
					t->next->prev = t->prev;
					_addChild(n, t);
				}
				else {
					t->prev->next = t->next;
					t->next->prev = t->prev;
					if (n->next == n) {
						t->next = t->prev = t;
						_addChild(t, n);
						n = t;
					}
					else {
						n->prev->next = t;
						n->next->prev = t;
						t->next = n->next;
						t->prev = n->prev;
						_addChild(t, n);
						n = t;
					}
				}
				continue;
			}
			else {
				trees[n->degree] = n;
			}
			n = n->next;
		}
		node<V>* min = n;
		node<V>* start = n;
		do {
			if (n->value<min->value)min = n;
			n = n->next;
		} while (n != start);
		return min;
	}

	node<V>* _cut(node<V>* heap, node<V>* n) {
		if (n->next == n) {
			n->parent->child = nullptr;
		}
		else {
			n->next->prev = n->prev;
			n->prev->next = n->next;
			n->parent->child = n->next;
		}
		n->next = n->prev = n;
		n->marked = false;
		return _merge(heap, n);
	}

	node<V>* _decreaseKey(node<V>* heap, node<V>* n, V value) {
		if (n->value<value)return heap;
		n->value = value;
		if (n->parent) {
			if (n->value<n->parent->value) {
				heap = _cut(heap, n);
				node<V>* parent = n->parent;
				n->parent = nullptr;
				while (parent != nullptr && parent->marked) {
					heap = _cut(heap, parent);
					n = parent;
					parent = n->parent;
					n->parent = nullptr;
				}
				if (parent != nullptr && parent->parent != nullptr)parent->marked = true;
			}
		}
		else {
			if (n->value < heap->value) {
				heap = n;
			}
		}
		return heap;
	}

	// Неприемлемо медленный поиск _find.
	// В следующей версии мы заменим операцию
	//поиска на быстродействующую хеш-таблицу.
	node<V>* _find(node<V>* heap, V value) {
		node<V>* n = heap;
		if (n == nullptr)return nullptr;
		do {
			if (n->value == value)return n;
			node<V>* ret = _find(n->child, value);
			if (ret)return ret;
			n = n->next;
		} while (n != heap);
		return nullptr;
	}
};
*/

/* Copyright(c) 2010, Robin Message <Robin.Message@cl.cam.ac.uk>
Все права защищены.
Распространение и использование в виде исходного и двоичного кода, с или без
модификация допускается при соблюдении следующих условий:
*При распространении исходного кода должны сохраняться вышеуказанные авторские права
обратите внимание, этот список условий и следующий отказ от ответственности.
* При повторном распространении двоичного кода должна сохраняться
обратите внимание, этот список условий и следующий отказ от ответственности в
документация и / или другие материалы, поставляемые вместе с дистрибутивом.
* Ни название University Кембриджа, ни
имена его участников могут использоваться для поддержки или продвижения продуктов
производные от данного программного обеспечения без предварительного письменного разрешения.
ЭТО ПРОГРАММНОЕ ОБЕСПЕЧЕНИЕ ПРЕДОСТАВЛЯЕТСЯ ВЛАДЕЛЬЦАМИ АВТОРСКИХ ПРАВ И УЧАСТНИКАМИ "КАК ЕСТЬ" И
ЛЮБЫЕ ЯВНЫЕ ИЛИ ПОДРАЗУМЕВАЕМЫЕ ГАРАНТИИ, ВКЛЮЧАЯ, НО НЕ ОГРАНИЧИВАЯСЬ, ПОДРАЗУМЕВАЕМЫЕ
ГАРАНТИИ ТОВАРНОЙ ПРИГОДНОСТИ И ПРИГОДНОСТИ ДЛЯ ОПРЕДЕЛЕННОЙ ЦЕЛИ
ОТКАЗАВШИЙСЯ.НИ В КОЕМ СЛУЧАЕ КЕМБРИДЖСКИЙ УНИВЕРСИТЕТ
*/
/*
template <class V> class FibonacciHeap;
const integer size_fibonacci_cashe = 2047483647;

template <class V> struct node {
private:
// Указатель на левый сестринский узел.
node<V>* prev;
// указатель на правый сестринский узел.
node<V>* next;
// указатель на один из дочерних узлов.
node<V>* child;
// указатель на родительский узел.
node<V>* parent;
V value;

// количество дочерних узлов.
int degree;


//логическое значение, которое указывает,
//были ли потери узлом x дочерних узлов,
//начиная с момента, когда  x стал дочерним
//узлом какого-то другого узла.
//FIBONNACCI_HEAP
bool marked;
public:
	friend class FibonacciHeap<V>;
		node<V>* getPrev() { return prev; }
		node<V>* getNext() { return next; }
		node<V>* getChild() { return child; }
		node<V>* getParent() { return parent; }
		V getValue() { return value; }
		bool isMarked() { return marked; }

		bool hasChildren() { return child; }
		bool hasParent() { return parent; }
};

template <class V> class FibonacciHeap {
protected:
	node<V>* heap;
	node<V>** hash_index; // Хеш таблица !!!
public:

FibonacciHeap() {
	heap = _empty();

	hash_index = nullptr;
}

void WakeUp2() {
	hash_index = new node<V>*[size_fibonacci_cashe];
	for (integer i = 0; i < size_fibonacci_cashe; i++) hash_index[i] = nullptr;
}

void WakeUp() {
	heap = _empty();

	hash_index = new node<V>*[size_fibonacci_cashe];
	for (integer i = 0; i < size_fibonacci_cashe; i++) hash_index[i] = nullptr;
}

void Clear() {
	if (heap) {
		if (hash_index != nullptr) {
			for (integer i = 0; i < size_fibonacci_cashe; i++) hash_index[i] = nullptr;
		}
		_deleteAll(heap);
	}
}

virtual ~FibonacciHeap() {
	if (heap) {
		for (integer i = 0; i < size_fibonacci_cashe; i++) hash_index[i] = nullptr;
		delete[] hash_index;
		hash_index = nullptr;
		_deleteAll(heap);
	}
	else {
		if (hash_index != nullptr) {
			for (integer i = 0; i < size_fibonacci_cashe; i++) hash_index[i] = nullptr;
			delete[] hash_index;
			hash_index = nullptr;
		}
	}
}
node<V>* insert(V value) {
	node<V>* ret = _singleton(value);
	heap = _merge(heap, ret);
	hash_index[-value] = ret;
	return ret;
}
void merge(FibonacciHeap& other) {
	heap = _merge(heap, other.heap);
	other.heap = _empty();
}

bool isEmpty() {
	return heap == nullptr;
}

V getMinimum() {
	return heap->value;
}

V removeMinimum() {
	node<V>* old = heap;
	heap = _removeMinimum(heap);
	V ret = old->value;
	delete old;
	old = nullptr;
	return ret;
}

void decreaseKey(node<V>* n, V value) {
heap = _decreaseKey(heap, n, value);
}

void deleteKey(V value) {
//node<V>* find_ = find(value);
	if (hash_index != nullptr) {
		node<V>* find_ = hash_index[-value];
		if (find_ != nullptr) {
			hash_index[-value] = nullptr;
#if doubleintprecision == 1
				decreaseKey(find_, -big_FIBO_integer_Value);
#else
				decreaseKey(find_, -2147483645);
#endif
			removeMinimum();
		}
	}
}

//fibo_n = fibo_heap.find(-veb_dsearch_key);
//if (fibo_n == nullptr) {
//	fibo_heap.insert(-veb_dadd_key);
//}
//else {
//	fibo_heap.decreaseKey(fibo_n, -veb_dadd_key);
//}
//fibo_n = nullptr;
void insert_and_modify(V value_search, V value_add) {
	//node<V>* find_ = find(value_search);
	if (hash_index == nullptr) {
		insert(value_add);
	}
	else {
		node<V>* find_ = hash_index[-value_search];
		if (find_ == nullptr) {
			insert(value_add);
		}
		else {
			// меняем позицию указателя на find_.
			hash_index[-value_search] = nullptr;
			hash_index[-value_add] = find_;
			decreaseKey(find_, value_add);
		}
		find_ = nullptr;
	}
}

node<V>* find(V value) {
  //return _find(heap, value);
	return hash_index[-value];
}
private:
	node<V>* _empty() {
	return nullptr;
}

node<V>* _singleton(V value) {
	node<V>* n = new node<V>;
	n->value = value;
	n->prev = n->next = n;
	n->degree = 0;
	n->marked = false;
	n->child = nullptr;
	n->parent = nullptr;
	return n;
}

node<V>* _merge(node<V>* a, node<V>* b) {
	if (a == nullptr)return b;
	if (b == nullptr)return a;
	if (a->value>b->value) {
		node<V>* temp = a;
		a = b;
		b = temp;
	}
	node<V>* an = a->next;
	node<V>* bp = b->prev;
	a->next = b;
	b->prev = a;
	an->prev = bp;
	bp->next = an;
	return a;
}

void _deleteAll(node<V>* n) {
	if (n != nullptr) {
		node<V>* c = n;
		do {
			node<V>* d = c;
			c = c->next;
			_deleteAll(d->child);
			delete d;
			d = nullptr;
		} while (c != n);
	}
}

void _addChild(node<V>* parent, node<V>* child) {
	child->prev = child->next = child;
	child->parent = parent;
	parent->degree++;
	parent->child = _merge(parent->child, child);
}

void _unMarkAndUnParentAll(node<V>* n) {
	if (n == nullptr)return;
	node<V>* c = n;
	do {
		c->marked = false;
		c->parent = nullptr;
		c = c->next;
	} while (c != n);
}

node<V>* _removeMinimum(node<V>* n) {
	if (n == nullptr)return n;
	_unMarkAndUnParentAll(n->child);
	if (n->next == n) {
		n = n->child;
	}
	else {
		n->next->prev = n->prev;
		n->prev->next = n->next;
		n = _merge(n->next, n->child);
	}
	if (n == nullptr)return n;
	node<V>* trees[64] = { nullptr };

	while (true) {
		if (trees[n->degree] != nullptr) {
			node<V>* t = trees[n->degree];
			if (t == n)break;
			trees[n->degree] = nullptr;
			if (n->value<t->value) {
				t->prev->next = t->next;
				t->next->prev = t->prev;
				_addChild(n, t);
			}
			else {
				t->prev->next = t->next;
				t->next->prev = t->prev;
				if (n->next == n) {
					t->next = t->prev = t;
					_addChild(t, n);
					n = t;
				}
				else {
					n->prev->next = t;
					n->next->prev = t;
					t->next = n->next;
					t->prev = n->prev;
					_addChild(t, n);
					n = t;
				}
			}
			continue;
		}
		else {
			trees[n->degree] = n;
		}
		n = n->next;
	}
	node<V>* min = n;
	node<V>* start = n;
	do {
		if (n->value<min->value)min = n;
		n = n->next;
	} while (n != start);
	return min;
}

node<V>* _cut(node<V>* heap, node<V>* n) {
	if (n->next == n) {
		n->parent->child = nullptr;
	}
	else {
		n->next->prev = n->prev;
		n->prev->next = n->next;
		n->parent->child = n->next;
	}
	n->next = n->prev = n;
	n->marked = false;
	return _merge(heap, n);
}

node<V>* _decreaseKey(node<V>* heap, node<V>* n, V value) {
	if (n->value<value)return heap;
	n->value = value;
	if (n->parent) {
		if (n->value<n->parent->value) {
			heap = _cut(heap, n);
			node<V>* parent = n->parent;
			n->parent = nullptr;
			while (parent != nullptr && parent->marked) {
				heap = _cut(heap, parent);
				n = parent;
				parent = n->parent;
				n->parent = nullptr;
			}
			if (parent != nullptr && parent->parent != nullptr)parent->marked = true;
		}
	}
	else {
		if (n->value < heap->value) {
			heap = n;
		}
	}
	return heap;
}

// Неприемлемо медленный поиск _find.
// В следующей версии мы заменим операцию
//поиска на быстродействующую хеш-таблицу.
// Недостаток в том что хеш-таблица получается слишком большой,
// много времени на выделение памяти и на её освобождение.
node<V>* _find(node<V>* heap, V value) {
	node<V>* n = heap;
	if (n == nullptr)return nullptr;
	do {
		if (n->value == value)return n;
		node<V>* ret = _find(n->child, value);
		if (ret)return ret;
		n = n->next;
	} while (n != heap);
	return nullptr;
}
};
*/

// Наиболее адаптированная к программе AliceFlow версия фибоначчиевой кучи.
// 15.07.2018
/* Copyright(c) 2010, Robin Message <Robin.Message@cl.cam.ac.uk>
Все права защищены.
Распространение и использование в виде исходного и двоичного кода, с или без
модификация допускается при соблюдении следующих условий:
*При распространении исходного кода должны сохраняться вышеуказанные авторские права
обратите внимание, этот список условий и следующий отказ от ответственности.
* При повторном распространении двоичного кода должна сохраняться
обратите внимание, этот список условий и следующий отказ от ответственности в
документация и / или другие материалы, поставляемые вместе с дистрибутивом.
* Ни название University Кембриджа, ни
имена его участников могут использоваться для поддержки или продвижения продуктов
производные от данного программного обеспечения без предварительного письменного разрешения.
ЭТО ПРОГРАММНОЕ ОБЕСПЕЧЕНИЕ ПРЕДОСТАВЛЯЕТСЯ ВЛАДЕЛЬЦАМИ АВТОРСКИХ ПРАВ И УЧАСТНИКАМИ "КАК ЕСТЬ" И
ЛЮБЫЕ ЯВНЫЕ ИЛИ ПОДРАЗУМЕВАЕМЫЕ ГАРАНТИИ, ВКЛЮЧАЯ, НО НЕ ОГРАНИЧИВАЯСЬ, ПОДРАЗУМЕВАЕМЫЕ
ГАРАНТИИ ТОВАРНОЙ ПРИГОДНОСТИ И ПРИГОДНОСТИ ДЛЯ ОПРЕДЕЛЕННОЙ ЦЕЛИ
ОТКАЗАВШИЙСЯ.НИ В КОЕМ СЛУЧАЕ КЕМБРИДЖСКИЙ УНИВЕРСИТЕТ
*/

template <class V> class FibonacciHeap;


template <class V> struct node {
private:
	// Указатель на левый сестринский узел.
	node<V>* prev;
	// указатель на правый сестринский узел.
	node<V>* next;
	// указатель на один из дочерних узлов.
	node<V>* child;
	// указатель на родительский узел.
	node<V>* parent;
	V value;

	// количество дочерних узлов.
	int degree;


	//логическое значение, которое указывает,
	//были ли потери узлом x дочерних узлов,
	//начиная с момента, когда  x стал дочерним
	//узлом какого-то другого узла.
	//FIBONNACCI_HEAP
	bool marked;

	node<V>() {
		// Указатель на левый сестринский узел.
		prev=nullptr;
		// указатель на правый сестринский узел.
		next=nullptr;
		// указатель на один из дочерних узлов.
		child=nullptr;
		// указатель на родительский узел.
		parent=nullptr;
		value=0;

		// количество дочерних узлов.
		degree=0;


		//логическое значение, которое указывает,
		//были ли потери узлом x дочерних узлов,
		//начиная с момента, когда  x стал дочерним
		//узлом какого-то другого узла.
		//FIBONNACCI_HEAP
		marked = false;
	}

public:
	friend class FibonacciHeap<V>;
	node<V>* getPrev() { return prev; }
	node<V>* getNext() { return next; }
	node<V>* getChild() { return child; }
	node<V>* getParent() { return parent; }
	V getValue() { return value; }
	bool isMarked() { return marked; }

	bool hasChildren() { return child; }
	bool hasParent() { return parent; }
};

template <class V> struct FiboHashNode {
	node<V>* link;
	integer count_neighbour;

	FiboHashNode() {
		link=nullptr;
		count_neighbour=0;
	}
};

template <class V> class FibonacciHeap {
protected:
	node<V>* heap; // указатель на минимальный элемент в Фибоначчиевой пирамиде.
	FiboHashNode<V>* hash_index; // Хеш-таблица !!!
	integer isize; // рзмер хеш таблицы.
public:

	void put_out_a_link(integer i) {
		hash_index[i].link = nullptr;
		hash_index[i].count_neighbour = 0;
	}

	// Конструктор
	FibonacciHeap() {
		heap = _empty();

		isize = 0;//future n_a+1
		hash_index = nullptr;
	}

	// Инициализация хеш-таблицы. Выделение памяти.
	void WakeUp2(integer isize_loc) {
		isize = isize_loc;
		hash_index = new FiboHashNode<V>[isize];
		for (integer i = 0; i < isize; i++) {
			put_out_a_link(i);
		}
	}

	// Инициализация хеш-таблицы. Выделение памяти.
    // Инициализация указателя на кучу.
	void WakeUp(integer isize_loc) {
		heap = _empty();

		isize = isize_loc;
		hash_index = new FiboHashNode<V>[isize];
		for (integer i = 0; i < isize; i++) {
			put_out_a_link(i);
		}
	}

	// Очистка хеш-таблицы если куча не пуста.
	void Clear() {
		if (heap) {
			if (hash_index != nullptr) {
				// Память из под хеш таблицы не освобождается.
				for (integer i = 0; i < isize; i++) {
					put_out_a_link(i);
				}
			}
			_deleteAll(heap);
		}
	}

	// Присваиваем размеру хеш-таблицы новое значение.
    // Память не перевыделяем.
	void UpdateSize(integer isize_loc) {
		// Объём памяти выделенный под хеш таблицу остается неизменным.
		isize = isize_loc;
	}

	// Деструктор. Освобождение оперативной	 памяти из под пирамиды Fibonacci.
	virtual ~FibonacciHeap() {
		if (heap) {
			for (integer i = 0; i < isize; i++) {
				put_out_a_link(i);
			}
			delete[] hash_index;
			hash_index = nullptr;
			_deleteAll(heap);
		}
		else {
			if (hash_index != nullptr) {
				for (integer i = 0; i < isize; i++) {
					put_out_a_link(i);
				}
				delete[] hash_index;
				hash_index = nullptr;
			}
		}
	}

	// Вставка узла. Создаём пирамиду из одного узла и сливаем её 
    // с первоначальной пирамидой. Так как value длинное целое, то
    // она хранит как бы два ключа сразу: i – номер узла в сетке,
    // count_neighbour – число соседей. i и count_neighbour – нужны
    // для работы RS функции. Распаковка осуществляется с помощью 
    // операций деления на цело / и взятия остатка от деления %.
    // Ссылка на добавленный узел корректно заносится в хеш-таблицу. 
	node<V>* insert(V value) {
		node<V>* ret = _singleton(value);
		heap = _merge(heap, ret);

		integer i = ((abs(value)) % (isize));
		integer count_neighbour = ((abs(value)) / (isize));
		hash_index[i].link = ret;
		hash_index[i].count_neighbour = count_neighbour;

		return ret;
	}

	// Слияние двух куч. Освобождение памяти из под второй кучи.
	void merge(FibonacciHeap& other) {
		heap = _merge(heap, other.heap);
		other.heap = _empty();
	}

	// Пуста ли куча ?
	bool isEmpty() {
		return heap == nullptr;
	}

	// Указатель Фибоначчиевой кучи указывает на минимальный элемент.
    // Показать минимальный элемент.
	V getMinimum() {
		return heap->value;
	}

	// Удаление минимального элемента из кучи.
	V removeMinimum() {
		node<V>* old = heap;	
		integer i= -getMinimum();
		if (abs(i) != big_FIBO_integer_Value) {
			i = ((abs(i)) % (isize));
			put_out_a_link(i);
		}
		heap = _removeMinimum(heap);
		V ret = old->value;
		delete old;
		old = nullptr;
		return ret;
	}

	// На вход подаётся ссылка n на узел кучи значение ключа которого требуется 
    // уменьшить до нового меньшего значения value.
	void decreaseKey(node<V>* n, V value) {
		heap = _decreaseKey(heap, n, value);
	}

	// Удаление узла кучи значение ключа которого равно value.
	void deleteKey(V value) {
		//node<V>* find_ = find(value);
		// по значению ключа value получить 
        // параметры составного ключа (i, count_neighbour).
        // i – номер узла сетки.  
        // count_neighbour – число соседних узлов.
		integer i = ((abs(value)) % (isize));

		if (hash_index != nullptr) {
			node<V>* find_ = hash_index[i].link;
			if (find_ != nullptr) {
				// Удаляем узел из хеш-таблицы.
				put_out_a_link(i);

				// Уменьшаем текущее значение найденного
                // ключа до минус бесконечности.
                // В результате удаляемый ключ 
                // становится текущим минимумом в куче.
                // Удаляем минимальный элемент в куче,
                // тем самым удаляя запрашиваемый при вызове 
                // функции deleteKey узел.
				//if (isign == -1)
				{
#if doubleintprecision == 1
					decreaseKey(find_, -big_FIBO_integer_Value);
#else
					decreaseKey(find_, -2147483645);
#endif
				}
				
				removeMinimum();
			}
		}
	}

	// Удаление узла кучи значение ключа которого равно ddel.
	void deleteKey(data_BalTree ddel) {
		// ddel сразу ключ в терминах (i, count_neighbour).

		if (hash_index != nullptr) {
			node<V>* find_ = hash_index[ddel.i].link;
			if (find_ != nullptr) {
				// Удаляем узел из хеш-таблицы.
				put_out_a_link(ddel.i);

				// Уменьшаем текущее значение найденного
                // ключа до минус бесконечности.
                // В результате удаляемый ключ 
                // становится текущим минимумом в куче.
                // Удаляем минимальный элемент в куче,
                // тем самым удаляя запрашиваемый при вызове 
                // функции deleteKey узел.
#if doubleintprecision == 1
				decreaseKey(find_, -big_FIBO_integer_Value);
#else
				decreaseKey(find_, -2147483645);
#endif
				removeMinimum();
			}
		}
	}

	//fibo_n = fibo_heap.find(-veb_dsearch_key);
	//if (fibo_n == nullptr) {
	//	fibo_heap.insert(-veb_dadd_key);
	//}
	//else {
	//	fibo_heap.decreaseKey(fibo_n, -veb_dadd_key);
	//}
	//fibo_n = nullptr;



	// Ключ value_search ищется в пирамиде. Если ключа value_search в пирамиде
    // нет, то происходит добавление в пирамиду ключа с значением value_add.
    // Если ключ value_search в пирамиде есть, то узел с ним удаляется из пирамиды.
    // Далее происходит добавление в пирамиду строго меньшего ключа с значением value_add.
	void insert_and_modify(V value_search, V value_add) {
		
		//node<V>* find_ = find(value_search);
		// по значению ключа value_search получить 
        // параметры составного ключа (i, count_neighbour).
        // i – номер узла сетки.  
        // count_neighbour – число соседних узлов.
		
		integer i = ((abs(value_search)) % (isize));
		

		// Если пирамида пуста то просто добавление ключа.
		if (hash_index == nullptr) {
			insert(value_add);
		}
		else {
			node<V>* find_ = hash_index[i].link;
			// Ключа value_search в пирамиде нет.
            // Просто добавление ключа value_add.
			if (find_ == nullptr) {
				insert(value_add);
			}
			else {
				// Если добавляемый ключ не равен искомому ключу.
				if (value_search!= value_add) {
					// меняем позицию указателя на find_.
					put_out_a_link(i);

					// по значению ключа value_add получить 
                    // параметры составного ключа (i, count_neighbour).
                    // i – номер узла сетки.  
                    // count_neighbour – число соседних узлов.
					
					i = ((abs(value_add)) % (isize));
					integer count_neighbour = ((abs(value_add)) / (isize));

					hash_index[i].link = find_;
					hash_index[i].count_neighbour = count_neighbour;
					if (value_add > value_search) {
						// Экспериментально проверено что
						// такой ситуации быть не может 
						// при работе RS функции.
						printf("\nERROR!!! Fibonacci Heap insert_and_modify");
						printf("\nvalue_search=%lld", value_search);
						printf("value_add=%lld\n", value_add);
						system("PAUSE");
					}
					decreaseKey(find_, value_add);
				}
			}
			find_ = nullptr;
		}
	}

	// С помощью хеш-таблицы hash_index возвращает по значению ключа value
    // ссылку на узел пирамиды с этим ключом.
	node<V>* find(V value) {
		//return _find(heap, value);
		return hash_index[-value];
	}
private:

	// Возвращает нулевой указатель.
	node<V>* _empty() {
		return nullptr;
	}

	// Создаёт Фибоначчиеву пирамиду из одного узла.
    // Возвращает ссылку на корень пирамиды.
	node<V>* _singleton(V value) {
		node<V>* n = new node<V>;
		n->value = value;
		n->prev = n->next = n;
		n->degree = 0; // нет дочерних узлов.
		n->marked = false; // не было потери узлов.
		n->child = nullptr;
		n->parent = nullptr;
		return n;
	}

	// Слияние двух Фибоначчиевых пирамид a и b.
    // Возвращает ссылку на корень пирамиды результата слияния.
	node<V>* _merge(node<V>* a, node<V>* b) {
		if (a == nullptr)return b;
		if (b == nullptr)return a;
		if (a->value>b->value) {
			node<V>* temp = a;
			a = b;
			b = temp;
			temp = nullptr;
		}
		// Пирамида a хранит меньший ключ.
        // Слитие двух двухсвязных списков корней в один 
        // двухсвязный список.
		node<V>* an = a->next;
		node<V>* bp = b->prev;
		a->next = b;
		b->prev = a;
		an->prev = bp;
		bp->next = an;
		an = nullptr;
		bp = nullptr;
		return a;
	}

	// Рекурсивная процедура удаления – освобождения памяти
    // из под Фибоначчиевой пирамиды.
	void _deleteAll(node<V>* n) {
		if (n != nullptr) {
			node<V>* c = n;
			do {
				node<V>* d = c;
				c = c->next;
				_deleteAll(d->child);
				delete d;
				d = nullptr;
			} while (c != n);
		}
	}

	// Добавление дочернего child узла в Фибоначчиеву кучу parent.
	// Подвешивание к Фибоначчиеву дереву parent дочернего дерева 
    // child в результате чего parent указывает на дерево на единицу 
    // более высокого ранга. Данная операция требуется для «окучивания».
	void _addChild(node<V>* parent, node<V>* child) {
		child->prev = child->next = child;
		child->parent = parent;
		parent->degree++;
		parent->child = _merge(parent->child, child);
	}

	// Снятие пометки маркера и обнуление ссылки на родителя для
    // всех узлов в двухсвязном списке на который указывает передаваемый узел n.
	void _unMarkAndUnParentAll(node<V>* n) {
		if (n == nullptr)return;
		node<V>* c = n;
		do {
			c->marked = false;
			c->parent = nullptr;
			c = c->next;
		} while (c != n);
	}

	// Удаление минимального узла из Фибоначчиевой кучи n. После удаления минимума
    // производится окучивание Фибоначчиевой кучи, чтобы гарантировать логарифмическую 
    // стоимость операции удаления минимума из Фибоначчиевой кучи – не давать куче сильно 
    // вытягиваться в длину или высоту. Самая ресурсоёмкая операция.
	node<V>* _removeMinimum(node<V>* n) {
		if (n == nullptr)return n;
		// обнуляем ссылку на родительский узел.
        // Снимаем метку marked = false
		_unMarkAndUnParentAll(n->child);
		// Отсекаем детей узла n и добавляем детей узла n (n->child) в родительский 
        // корневой список деревьев. У нас в корневом списке много
        // Фибоначчиевых деревьев и ранги деревьев повторяются.
		if (n->next == n) {
			// Двухсвязный список корней состоит из 
            // одного узла n. Возвращаем двухсвязный список 
            // детей узла n в качестве корневого списка.
			n = n->child;
		}
		else {
			// Исключаем узел n из двухсвязного списка корней
			n->next->prev = n->prev;
			n->prev->next = n->next;
			// Сливаем двухсвязный список детей узла n (n->child) с
            // корневым списком.
			n = _merge(n->next, n->child);
		}
		if (n == nullptr)return n;
		// В массиве trees будем хранить ссылки на деревья из корневого
        // списка корней. Причем в ячейке degree будет хранится ссылка
        // на корень дерева из корневого списка с рангом равным degree.
        // В случае если в списке trees в одной ячейке будет два дерева
        // одинакового ранга degree то мы сольём эти два дерева и получим 
        // дерево ранга degree+1. Из двух деревьев ранга degree корнем 
        // дерева ранга (degree+1) станет то, у которого ключ наименьший.
        // Мы приследуем цель чтобы все деревья в корневом списке корней
        // были разных рангов. Данная операция называется «окучивание».
		node<V>* trees[64] = { nullptr };

		while (true) {
			if (n != nullptr) {
				if (trees[n->degree] != nullptr) {
					node<V>* t = trees[n->degree];
					if (t == n) break;
					// Делаем из двух деревьев ранга degree
                    // одно дерево рангом (degree+1).
					trees[n->degree] = nullptr;
					if (n->value < t->value) {
						// Исключаем узел t из двухсвязного 
                        // списка корней 
						t->prev->next = t->next;
						t->next->prev = t->prev;
						// Подвешиваем дерево с корнем t ребёнком к
                        // дереву с корнем n. n – корень дерева ранга
                        // (degree+1). 
                        // Ключ n меньше ключа t: n->value<t->value.
						_addChild(n, t);
					}
					else {
						// Исключаем узел t из двухсвязного 
                        // списка корней
						t->prev->next = t->next;
						t->next->prev = t->prev;
						// Корневой список деревьев на которые
                        // указывает указатель n состоит из
                        // одного дерева или из нескольких
                        // деревьев ?
						if (n->next == n) {
							// Дерево n единственное в корневом списке.

                            // Делаем из дерева t –
                            // единственное дерево в корневом списке.
							t->next = t->prev = t;
							// Ключ дерева t (t->value) меньше чем
                            // ключ дерева n (n->value). 
                            // t->value<n->value.
                            // К единственному дереву t в корневом списке
                            // подвешиваем единственное дерево n в 
                            // корневом списке.
							_addChild(t, n);
							// Возвращаем указатель на дерево t.
							n = t;
						}
						else {
							// Указатель n указывает на целый
							// двухсвязный список деревьев.
							// Исключение n узла из двухъсвязного
							// списка узлов.
							n->prev->next = t;
							n->next->prev = t;
							t->next = n->next;
							t->prev = n->prev;
							// Ключ дерева t (t->value) меньше чем
							// ключ дерева n (n->value). 
							// t->value<n->value.
							// К  дереву t в корневом списке,
							// подвешиваем единственное дерево n в 
							// корневом списке как сына.
							_addChild(t, n);
							// Возвращаем указатель на дерево t.
							n = t;
						}
					}
					continue;
				}
				else {
					// В массиве деревьев trees еще нет 
					// дерева n с рангом n->degree.
					// Помещаем дерево n в массив trees.
					trees[n->degree] = n;
				}
				// Переходим к следующему дереву в корневом списке.
				n = n->next;
			}
			else {
				break;
			}
		}
		// На данный момент все деревья в корневом списке
        // деревьев имеют разный ранг.
        // Деревьев в корневом списке не более логарифма 
        // по основанию 2 от количества элементов в Фибоначчиевой
        // куче.
        // Пробегаемся по двухсвязному списку деревьев в корневом
        // списке и ищем новый минимум. Т.к. деревьев не более логарифма
        // то операция поиска имеет сложность O(log_2(n)).
		node<V>* min_1 = n;
		node<V>* start = n;
		if (n != nullptr) {
			do {
				if (n->value < min_1->value) min_1 = n;
				n = n->next;
			} while (n != start);
		}
		// Возвращаем указатель на дерево в корневом списке
        // с минимальным значением ключа min->value – искомый 
		// минимальный ключ. min искомый указатель на результирующую 
		// Фибоначчиеву кучу, дерево в корневом списке с минимальным 
		// элементом.
		return min_1;
	}

	// Вырезаем из пирамиды heap вершину n.
	node<V>* _cut(node<V>* heap, node<V>* n) {
		if (n->next == n) {
			// Вершина n единственная в двухсвязном списке.
			// Из родителя убираем связь на вершину n.
			n->parent->child = nullptr;
		}
		else {
			// Исключаем узел n из двухсвязного 
			// списка корней 
			n->next->prev = n->prev;
			n->prev->next = n->next;
			// Из родителя убираем связь на вершину n,
			// делая ссылку на непустую вершину n->next.
			n->parent->child = n->next;
			// Теперь родитель указывает на новую 
			// вершину n->next, а вершина n вырезана
			// из кольцевого двухсвязного списка вершин.
		}
		// Исключаем дерево n из двухсязного
		// кольцевого списка. Теперь это одиночное
		// дерево.
		n->next = n->prev = n;
		// У которого не удаляли детей. 
		// Т.к. теперь n будет располагаться в 
		// корневом двухсвязном списке heap.
		// У члена корневого списка marked всегда равен false.
		n->marked = false;
		// Сливаем дерево n и исходную пирамиду heap,
		// путем связывания (объединения) двух
		// двухсвязных списков.
		return _merge(heap, n);
	}

	// В пирамиде heap ключ в узле с указателем n будет заменён
    // на меньшее значение ключа value. При этом данная операция не нарушит,
    // свойство кучи.
	node<V>* _decreaseKey(node<V>* heap, node<V>* n, V value) {
		// Априори value < n->value. 
		// Мы рассматриваем только уменьшение ключа
		// в узле n.
		if (n->value<value)return heap;
		n->value = value; // заменили ключ на меньший.
		// Проверяем не нарушен ли порядок кучи.
		// Восстанавливаем порядок кучи.
		if (n->parent) {
			// Вершина пирамиды еще не достигнута.
			if (n->value<n->parent->value) {
				// Порядок кучи нарушен.
				// Вырезаем из пирамиды heap вершину n.
				heap = _cut(heap, n);
				// запоминаем указатель на родителя 
				// у вырезанной вершины n.
				node<V>* parent = n->parent;
				// Гасим ссылку на родителя 
				// у вырезанной вершины n.
				n->parent = nullptr;
				// Пока корень пирамиды heap не достигнут
				// и текущий корень помечен что у него 
				// уже были потери детей
				while (parent != nullptr && parent->marked) {
					// Вырезаем из пирамиды heap вершину parent.
					heap = _cut(heap, parent);
					// запоминаем указатель на текущую
					// вырезанную вершину.
					n = parent;
					// Двигаемся вверх к корню пирамиды,
					// пока корень пирамиды не достигнут и
					// у текущего корневого узла уже
					// были потери детей.
					parent = n->parent;
					// Гасим ссылку на родителя 
					// у вырезанной вершины n.
					n->parent = nullptr;
				}
				// Цепочка вырезаний завершена на узле
				// parent у которого еще не отрезали детей. 
				// В случае существования у 
				// текущего узла parent родителя parent->parent!=NULL,
				// отмечаем у текущего узла parent
				// маркер parent->marked = true;
				// Это означает у узла parent отрезали
				// ребёнка.
				if (parent != nullptr && parent->parent != nullptr)parent->marked = true;
			}
		}
		else {
			if (n->value < heap->value) {
				// Выбрали меньший из двух ключей.
				// Теперь указатель пирамиды указывает на вершину
				// с меньшим значением ключа.
				heap = n;
			}
		}
		return heap;
	}

	// Неприемлемо медленный поиск _find.
	// В следующей версии мы заменим операцию
	//поиска на быстродействующую хеш-таблицу.
	// Недостаток в том что хеш-таблица получается слишком большой,
	// много времени на выделение памяти и на её освобождение.
	node<V>* _find(node<V>* heap, V value) {
		node<V>* n = heap;
		if (n == nullptr)return nullptr;
		do {
			if (n->value == value)return n;
			node<V>* ret = _find(n->child, value);
			if (ret)return ret;
			n = n->next;
		} while (n != heap);
		return nullptr;
	}

	


};

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
// привдённый в конце это полноценный 100% рабочий гибрид этой идеи 
// с быстродействующей хеш таблицей. 
// Это должен быть очень быстрый код т.к. 1. нет постоянных new и delete -
// память выделяется лишь один раз в момент инициализации и освобождается автоматическим
// деструктором. 2. Рекурсия заменена сверхбыстродействующей итерацией.
// Критика : А где балансировка обеспечивающая логарифмическую длину поиска ?
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
		if (pq != NULL) delete[] pq;
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
	// Просто прочитать масимальный элемент.
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
// Есть ли элемент с данным ключём в таблице ?
//bool isfound(integer key);
//bool empty() const;
// Очищаем содержимое и она снова готова к использованию.
//void clear();
// Вернуть элемент с заданным ключём:
// Обязательно предполагается что ключ существует внутри таблицы.
//Item get(integer key);
//Item readmax();
//integer readkeymaxelm();
// Вставить элемент item в очередь по 
// приоритетам если элемент item имеет ключ key.
//template <class Item>
//void insert(Item item, integer key);
// Возвратить максимальный элемент
// и удалить его.
//Item getmax();
// Заменяет элемент с ключём key на элемент val с тем же ключём key.
// При этом ключ key должен быть уникальным.
//void modify(integer key, Item val);
// Удаление элемента с заданным значением ключа.
//void remove(integer key);
// У элемента изменть значение старого ключа на новый ключ
// при этом меняется и само содержимое элемента.
//void change(integer key_serch, integer key_new, integer item_new);


// Ключи должны быть уникальны, целочисленны и различны.
// Двух одинаковых ключей быть недолжно, иначе коллизия в хеш таблице.

template <class Item>
class PQ
{
private:
	// Хранение binary heap.
	Item *pq;
	// Обратный доступ по номеру в qp на ячейку в hash.
	integer *qp; // Ссылка на хеш таблицу.
	// Доступ по ключу к полю в pq.
	integer *hash; // Хеш таблица !!!
	integer N;
	integer isize;
	integer ihash_size;

public:
	PQ(integer maxN, integer max_key_size)
	{
		isize = maxN;
		pq = new Item[maxN + 1];
		qp = new integer[maxN + 1];
		for (integer i_1 = 0; i_1 < maxN + 1; i_1++) {
			// Инициализация: таблица пуста т.к. поле 0
			// в массиве pq никогда не используется.
			qp[i_1] = 0;
		}
		N = 0;
		// Хеш таблица !!!
		ihash_size = max_key_size;
		hash = new integer[max_key_size+2];
		for (integer i_1 = 0; i_1 < max_key_size + 2; i_1++) {
			// Инициализация: таблица пуста т.к. поле 0
			// в массиве pq никогда не используется.
			hash[i_1] = 0;
		}
	}
	~PQ()
	{
		if (pq != NULL) delete[] pq;
		N = 0;
		if (qp != NULL) delete[] qp;
		if (hash != NULL) delete[] hash;
	}
	void print_log(char ch) {
		printf("%c\n",ch);
		for (integer i_1 = 1; i_1 <= N; i_1++) {
#if doubleintprecision == 1
			printf("[%lld %lld] ",pq[i_1],qp[i_1]);
#else
			printf("[%d %d] ", pq[i_1], qp[i_1]);
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
		Item t1 = get(i);
		//Item t2 = get(j);

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
		for (integer i_1 = 0; i_1 < N + 1; i_1++) {
			// Ускоренная очистка хеш таблицы.
			hash[qp[i_1]] = 0;
		}
		/*
		// Стабильно верный вариант.
		for (integer i_1 = 0; i_1 < isize + 1; i_1++) {
			// Инициализация: таблица пуста т.к. поле 0
			// в массиве pq никогда не используется.
			qp[i_1] = 0;
		}
		*/
		// 25.02.2018
		// Быстродействие выше т.к. мы просматриваем только 
		// ключи лежащие в очереди а не весь диапазон с запасом.
		for (integer i_1 = 0; i_1 < N + 1; i_1++) {
			qp[i_1] = 0;
		}
		N = 0;
		/*
		for (integer i_1 = 0; i_1 < ihash_size + 2; i_1++) {
			// Инициализация: таблица пуста т.к. поле 0
			// в массиве pq никогда не используется.
			hash[i_1] = 0;
		}
		*/
	}
	bool empty() const
	{
		return N == 0;
	}
	// Есть ли элемент с данным ключём в таблице ?
	bool isfound(integer key) {
		if (hash[key] == 0) {
			// Элемент отсутствует в хеш таблице.
			return false;
		}
		return true;
	}
	// Вернуть элемент с заданным ключём:
	// Обязательно предполагается что ключ существует внутри таблицы.
	Item get(integer key) {
		if (hash[key] == 0) {
			// Элемент отсутствует в хеш таблице.
			printf("priority queue get ERROR: get element not found.\n");
			system("pause");
			exit(1);
		}
		return pq[hash[key]];
	}
	// Просто прочитать масимальный элемент.
	Item readmax()
	{
		return pq[1];
	}
	integer readkeymaxelm() {
		return qp[1];
	}

	

	// Вставить элемент item в очередь по 
	// приоритетам если элемент item имеет ключ key.
	template <class Item>
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

	

	// Заменяет элемент с ключём key на элемент val с тем же ключём key.
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
				// Элемент отсутствует в хеш таблице.
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
	
	// У элемента изменть значение старого ключа на новый ключ
	// при этом меняется и само содержимое элемента.
	void change(integer key_serch, integer key_new, integer item_new)
	{
		if (hash[key_serch] == 0) {
			// Элемент отсутствует в хеш таблице.
			if (hash[key_new] != 0) {
				// Элемент присутствует в хеш таблице.
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
				// Элемент присутствует в хеш таблице.
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

// Все операции O(1) и только операция delete log2(N).

typedef struct T_Fibbonacci_node {
	integer key;

	// указатель на родительский узел.
	T_Fibbonacci_node* p;
	// указатель на один из дочерних узлов.
	T_Fibbonacci_node* child;
	// Указатель на левый сестринский узел.
	T_Fibbonacci_node* left;
	// указатель на правый сестринский узел.
	T_Fibbonacci_node* right;

	// количество дочерних узлов.
	integer degree;
	/*
	логическое значение, которое указывает,
	были ли потери узлом x дочерних узлов,
	начиная с момента, когда  x стал дочерним 
	узлом какого-то другого узла.
	*/
	bool mark;

} Fibbonacci_node;


// Сбалансированные деревья поиска.
// Используются в алгебраическом многосеточном методе.
#pragma once
#ifndef _AVL_BALANCED_TREE_CPP_
#define _AVL_BALANCED_TREE_CPP_ 1

// Для АВЛ дерева
// 12.12.2015
// По видимому в 3D очень большое количество соседей и простой 
// линейный двухсвязный список не справляется с таким большим 
// количеством элементов по параметру быстродействие. 
// Для сравнения в 2D при 1М неизвестных видно что линейный
// список справляется прекрасным образом со своим 2D числом соседей
// и даже по результатам 
// профайлинга не видно чтобы он испытывал сколько-нибудь ощутимую
// нагрузку. 

// АВЛ дерево это структура данных позволяющая производить все 
// операции за log_2 (number of neighbour count) операций. 

// Г.М. Адельсон-Вельский и Е.М. Ландис 1962 
// (взято из интернета Статья Николая Ершова на Хабре.).
// узел с максимальным значением ключа крайний правый.

// для поля данных надо определить перегруженные операции,
// меньше, больше и равно.

// Поле данных в АВЛ дереве.
// перенесена в файл priority_queue.cpp
//class data_BalTree
//{
	// public:	
	// --> high priority --> for operation <,>
	//integer ii;
	//integer  i, count_neighbour;
	// count_neighbour есть key.
	// Перегруженные операции сравнения для составного ключа.
	//bool operator <(const data_BalTree&);
	//bool operator >(const data_BalTree&);
	//bool operator ==(const data_BalTree&);
//};

//отношения порядка только в insert и remove

// Узел АВЛ дерева.
struct node_AVL
{
	data_BalTree key;
	// Высота поддерева с корнем в данном узле.
	unsigned char height;
	node_AVL* left;
	node_AVL* right;
	// Конструктор.
	node_AVL(data_BalTree k) { key = k; left = right = 0; height = 1; }
	node_AVL() {
		height=1;
		left = nullptr;
		right = nullptr;
	}
};

// Работает также и с пустыми деревьями.
// Обёртка поля height
unsigned char height(node_AVL* p)
{
	return p ? p->height : 0;
};

// Вычисляет balance factor заданного узла
// работает только с ненулевыми указателями.
integer bfactor(node_AVL* p)
{
	return height(p->right) - height(p->left);
};

// Восстанавливает корректное значение поля height
// заданного узла (при условии, что значения этого поля 
// в правом и левом дочерних узлах являются корректными).
void fixheight(node_AVL* p)
{
	unsigned char hl = height(p->left);
	unsigned char hr = height(p->right);
	p->height = (hl > hr ? hl : hr) + 1;
};

// Балансировка узлов.
node_AVL* rotateright(node_AVL* p)
{
	// правый поворот вокруг p
	node_AVL* q = p->left;
	p->left = q->right;
	q->right = p;
	fixheight(p);
	fixheight(q);
	return q;
};

// Левый поворот является симметричной копией правого:
node_AVL* rotateleft(node_AVL* q)
{
	// левый поворот вокруг q
	node_AVL* p = q->right;
	q->right = p->left;
	p->left = q;
	fixheight(q);
	fixheight(p);
	return p;
};

// Код выполняющий балансировку сводится к проверке условий и выполнению поворотов
node_AVL* balance(node_AVL* p) // балансировка узла p
{
	fixheight(p);
	if (bfactor(p) == 2)
	{
		if (bfactor(p->right) < 0)
			p->right = rotateright(p->right);
		return rotateleft(p);
	}
	if (bfactor(p) == -2)
	{
		if (bfactor(p->left) > 0)
			p->left = rotateleft(p->left);
		return rotateright(p);
	}
	return p; // балансировка не нужна
}

// Вставка ключей в дерево.
// Возвращает новое значение корня АВЛ дерева.
node_AVL* insert(node_AVL*& p, data_BalTree k)
{
	// Вставка ключа k в дерево с корнем p
	if (p == nullptr) {
		node_AVL* r1 = nullptr;
		r1 = new node_AVL(k);
		if (r1 == nullptr) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for r1 in insert my_agregat_amg...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
		return r1;
	}
	if (k < p->key)
		p->left = insert(p->left, k);
	else
		p->right = insert(p->right, k);

	return balance(p);
} // insert

  // Возвращает true если узел найден в дереве
bool isfound(node_AVL* p, data_BalTree k)
{
	if (p == 0) return false; // ненайден.
	if (k < p->key)
		return isfound(p->left, k);
	else if (k > p->key)
		return isfound(p->right, k);
	else return true; // найден.
}

// Полное удаление бинарного дерева.
void clear_AVL(node_AVL* p)
{
	if (p != 0) {
		clear_AVL(p->left);
		clear_AVL(p->right);
		// удаляем лист.
		delete p;
		p = 0;
	}
} // clear_AVL

  // Удаление узла с заданными свойства с сохранением сбалансированности.
node_AVL* findmin(node_AVL* p)
{
	// поиск узла с минимальным ключом в дереве p
	//if (!p) {
	return p->left ? findmin(p->left) : p;
	//}
	//else {
	// на поиск минимума подан нулевой указатель.
	//return 0;
	//}
} // findmin

node_AVL* findmax(node_AVL* p)
{
	// поиск узла с максимальным ключом в дереве p
	if (p != 0) {
		return p->right ? findmax(p->right) : p;
		/*
		#if doubleintprecision == 1
			if (p->right == 0) {
				//printf("%lld %lld %lld\n", p->key.count_neighbour, p->key.i, p->key.ii);
				//system("pause");
				return p;
			}
			else {
				//printf("%lld %lld %lld\n",p->key.count_neighbour,p->key.i,p->key.ii);
				findmax(p->right);
			}
		#else
			if (p->right == 0) {
				//printf("%d %d %d\n", p->key.count_neighbour, p->key.i, p->key.ii);
				//system("pause");
				return p;
			}
			else {
				//printf("%d %d %d\n",p->key.count_neighbour,p->key.i,p->key.ii);
				findmax(p->right);
			}
		#endif

		*/
	}
	else {
		// На поиск максимума подан нулевой указатель.
		return 0;
	}
} // findmax

data_BalTree get_max_AVL(node_AVL* p)
{
	// возвращение максимального узла в дереве.
	return p->right ? get_max_AVL(p->right) : p->key;
}

node_AVL* removemin(node_AVL* p)
{
	// удаление узла с минимальным ключом из дерева p
	if (p->left == 0)
		return p->right;
	p->left = removemin(p->left);
	return balance(p);
}

// Удаление заданного элемента из AVL дерева
// с полным сохранением балансировки.
// на возвращаемое значение можно не обращать внимания.
node_AVL* remove_AVL(node_AVL* p, data_BalTree k)
{
	// Отношение порядка определено 
	// для структуры из трёх целых чисел.
	// удаление ключа k из дерева p
	if (p == 0) return 0;
	// Двоичный поиск нужного элемента.
	if (k < p->key)
		p->left = remove_AVL(p->left, k);
	else if (k > p->key)
		p->right = remove_AVL(p->right, k);
	else // k==p->key
	{
		node_AVL* q = p->left;
		node_AVL* r = p->right;
		delete p;
		p = 0;
		if (r == 0) return q;
		node_AVL* min = findmin(r);
		min->right = removemin(r);
		min->left = q;
		return balance(min);
	}

	// При выходе из рекурсии делаем балансировку.
	return balance(p);
}

// Вставка ключа К в дерево если ключа k_search
// еще нет в дереве или модификация ключа К_search на К.
// Возвращает новое значение корня АВЛ дерева.
node_AVL* insert_and_modify(node_AVL* p, data_BalTree k, data_BalTree k_search)
{
	if (isfound(p, k_search) == false) {
		//   узла в дереве нет.
		p = insert(p, k);
		return p;
	}
	else {
		// удаление k_search
		p = remove_AVL(p, k_search); // необходимое действие
									 //remove_AVL(p, k_search); // приводит к ошибке.
		p = insert(p, k);							 // вставка к.
		return p;
	}
}

// print_AVL for debug.
void print_AVL(node_AVL* p)
{
	if (p != 0) {
		print_AVL(p->left);
		for (integer i = 0; i <= p->height; i++) {
			printf(" ");
		}
#if doubleintprecision == 1
		//printf("%lld %lld %lld\n", p->key.count_neighbour, p->key.i, p->key.ii);
		printf("%lld %lld\n", p->key.count_neighbour, p->key.i);
#else
		//printf("%d %d %d\n", p->key.count_neighbour, p->key.i, p->key.ii);
		printf("%d %d\n", p->key.count_neighbour, p->key.i);
#endif

		print_AVL(p->right);
	}
} // print_AVL

  // тестирование АВЛ дерева
  // 12 декабря 2015 удовлетворительно.
void test_AVL()
{
	node_AVL* root = 0;
	data_BalTree d3;
	d3.count_neighbour = rand();
	d3.i = rand();
	//d3.ii = rand();
	root = remove_AVL(root, d3);
	for (integer i = 0; i < 10; i++)
	{
		data_BalTree d;
		d.count_neighbour = rand();
		d.i = rand();
		//d.ii = rand();
		root = insert_and_modify(root, d, d);
		print_AVL(root);
		if (i == 5) {

			//41 18467 6334


			data_BalTree d1;
			d1.count_neighbour = 41;
			d1.i = 18467;
			//d1.ii = 6334;
			data_BalTree d2;
			d2.count_neighbour = rand();
			d2.i = rand();
			//d2.ii = rand();
			root = insert_and_modify(root, d2, d1);
			print_AVL(root);
			printf("remove 41 18467 6334\n");
			root = remove_AVL(root, d1);
			print_AVL(root);
			printf("found max\n");
			node_AVL* emax;
			emax = findmax(root);
#if doubleintprecision == 1
			printf("maximum id %lld\n", emax->key.count_neighbour);
#else
			printf("maximum id %d\n", emax->key.count_neighbour);
#endif

			emax = 0;
			//clear_AVL(root);
			//root = 0;
			print_AVL(root);
		}

		system("pause");
	}
	system("pause");
}



// Конец данных и методов используемых в АВЛ дереве.

#endif
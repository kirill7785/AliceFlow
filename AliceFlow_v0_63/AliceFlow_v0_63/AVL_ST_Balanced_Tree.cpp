// Сбалансированные деревья поиска.
// Используются в алгебраическом многосеточном методе.
#pragma once
#ifndef _AVL_ST_BALANCED_TREE_CPP_
#define _AVL_ST_BALANCED_TREE_CPP_ 1

#include "Bidirectional_Linear_List.cpp"

// 5.01.2017
// Для нахождения множества Strong Transpose нам 
// потребуется хеш-таблица АВЛ деревьев. (Лес АВЛ деревьев).
// В литературе пишут что при создании C/F разбиения используют
// именно Strong Transpose множества соседей.
// Поле данных в АВЛ дереве.
struct data_BalTreeST
{
	// --> high priority --> for operation <,>
	integer  i;
	// i есть key.

	data_BalTreeST() {
		i = -1;
	}
};

//отношения порядка только в insert и remove

// Узел АВЛ дерева.
struct node_AVLST
{
	data_BalTreeST key;
	// Высота поддерева с корнем в данном узле.
	unsigned char height;
	node_AVLST* left;
	node_AVLST* right;
	// Конструктор.
	node_AVLST(data_BalTreeST k) { key = k; left = right = nullptr; height = 1; }

	node_AVLST() {
		height = 1;
		left = nullptr;
		right = nullptr;
	}

};

// Работает также и с пустыми деревьями.
// Обёртка поля height
unsigned char height(node_AVLST* p)
{
	return p ? p->height : 0;
};

// Вычисляет balance factor заданного узла
// работает только с ненулевыми указателями.
integer bfactor(node_AVLST* p)
{
	return height(p->right) - height(p->left);
};

// Восстанавливает корректное значение поля height
// заданного узла (при условии, что значения этого поля 
// в правом и левом дочерних узлах являются корректными).
void fixheight(node_AVLST* p)
{
	unsigned char hl = height(p->left);
	unsigned char hr = height(p->right);
	p->height = (hl > hr ? hl : hr) + 1;
};

// Балансировка узлов.
node_AVLST* rotateright(node_AVLST* p)
{
	// правый поворот вокруг p
	node_AVLST* q = p->left;
	p->left = q->right;
	q->right = p;
	fixheight(p);
	fixheight(q);
	return q;
};

// Левый поворот является симметричной копией правого:
node_AVLST* rotateleft(node_AVLST* q)
{
	// левый поворот вокруг q
	node_AVLST* p = q->right;
	q->right = p->left;
	p->left = q;
	fixheight(q);
	fixheight(p);
	return p;
};

// Код выполняющий балансировку сводится к проверке условий и выполнению поворотов
node_AVLST* balance(node_AVLST* p) // балансировка узла p
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
node_AVLST* insert(node_AVLST*& p, data_BalTreeST k)
{
	// Вставка ключа k в дерево с корнем p
	if (p == nullptr) {
		node_AVLST* r1 = nullptr;
		r1 = new node_AVLST(k);
		if (r1 == nullptr) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for r1 in insert my_agregat_amg...\n");
			printf("Please any key to exit...\n");
			//system("pause");
			system("pause");
			exit(1);
		}
		return r1;
	}
	if (k.i >= p->key.i)
		p->right = insert(p->right, k);
	else
		p->left = insert(p->left, k);

	return balance(p);
} // insert

  // Возвращает true если узел найден в дереве
bool isfound(node_AVLST* p, data_BalTreeST k)
{
	if (p == nullptr) return false; // ненайден.
	if (k.i < p->key.i)
		return isfound(p->left, k);
	else if (k.i > p->key.i)
		return isfound(p->right, k);
	else return true; // найден.
}

// Полное удаление бинарного дерева.
void clear_AVL(node_AVLST* p)
{
	if (p != nullptr) {
		clear_AVL(p->left);
		clear_AVL(p->right);
		// удаляем лист.
		delete p;
		p = nullptr;
	}
} // clear_AVL

/*
// НЕ РАБОТАЕТ
// Ни в коем случае не раскомментировать. 6 июня 2017.
// Скорее всего нужна parent ссылка.
  // Полное удаление бинарного дерева.
void clear_AVL(node_AVLST* top)
{


	node_AVLST* stack[128];
	integer count = 0;

	while (!count == 0 || top != nullptr) {
		if (!count == 0) {
			top = stack[count];
			stack[count] = nullptr;
			count--;
			if (!count == 0 && top->right == stack[count]) {
				top = stack[count];
				stack[count] = nullptr;
				count--;
			}
			else {
				delete top;
				top = nullptr;
			}
		}

		while (top != nullptr) {
			stack[count++]=top;
			if (top->right != nullptr) {
				stack[count++]=(top->right);
				stack[count++] = top;
			}
			top = top->left;
		}
		*/
		/*
		if (p) {
			stack[count++] = p;
			p = p->left;
		}
		else {
			peekn = stack[count];
			if (peekn->right && lnp != peekn->right) {
				p = peekn->right;
			}
			else {

				stack[count] = nullptr;
				count--;
				p = stack[count];
				#if doubleintprecision == 1
					//printf("visited %lld\n", peekn->data);
				#else
					//printf("visited %d\n", peekn->data);
				#endif

				delete peekn;
				peekn = nullptr;
				lnp = peekn;
			}
		}
		*//*
	}



} // clear_AVL
*/

// Предварительный подсчёт размерности множества set при 
  // передачи списка Si Transpose из АВЛ дерева в множество set.
void formirate_F_SiTranspose_get_num_elements_set(node_AVLST*& p, integer& imarker75_scan)
{
	if (p != nullptr) {
		formirate_F_SiTranspose_get_num_elements_set(p->left, imarker75_scan);
		formirate_F_SiTranspose_get_num_elements_set(p->right, imarker75_scan);
		imarker75_scan++;

	}
} // formirate_F_SiTranspose_get_num_elements_set

// Передача списка Si Transpose из АВЛ дерева в множество set.
void formirate_F_SiTranspose(node_AVLST*& p, integer*& set, integer& imarker75_scan)
{
	if (p != nullptr) {
		formirate_F_SiTranspose(p->left, set, imarker75_scan);
		formirate_F_SiTranspose(p->right, set, imarker75_scan);
		set[imarker75_scan] = p->key.i;
		imarker75_scan++;

	}
}




// Удаление узла с заданными свойства с сохранением сбалансированности.
node_AVLST* findmin(node_AVLST* p)
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

node_AVLST* findmax(node_AVLST* p)
{
	// поиск узла с минимальным ключом в дереве p
	if (p != nullptr) {
		return p->right ? findmax(p->right) : p;
	}
	else {
		// На поиск максимума подан нулевой указатель.
		return nullptr;
	}
} // findmax

data_BalTreeST get_max_AVL(node_AVLST* p)
{
	// возвращение максимального узла в дереве.
	return p->right ? get_max_AVL(p->right) : p->key;
}

node_AVLST* removemin(node_AVLST* p)
{
	// удаление узла с минимальным ключом из дерева p
	if (p->left == nullptr)
		return p->right;
	p->left = removemin(p->left);
	return balance(p);
}

// Удаление заданного элемента из AVL дерева
// с полным сохранением балансировки.
// на возвращаемое значение можно не обращать внимания.
node_AVLST* remove_AVL(node_AVLST* p, data_BalTreeST k)
{
	// Отношение порядка определено 
	// для структуры из трёх целых чисел.
	// удаление ключа k из дерева p
	if (p == nullptr) return nullptr;
	// Двоичный поиск нужного элемента.
	if (k.i < p->key.i)
		p->left = remove_AVL(p->left, k);
	else if (k.i > p->key.i)
		p->right = remove_AVL(p->right, k);
	else // k==p->key
	{
		node_AVLST* q = p->left;
		node_AVLST* r = p->right;
		delete p;
		p = nullptr;
		if (r == nullptr) return q;
		node_AVLST* min = findmin(r);
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
node_AVLST* insert_and_modify(node_AVLST* p, data_BalTreeST k, data_BalTreeST k_search)
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
void print_AVL(node_AVLST* p)
{
	if (p != nullptr) {
		print_AVL(p->left);
		for (integer i = 0; i <= p->height; ++i) {
			printf(" ");
		}
#if doubleintprecision == 1
		//printf("%lld %lld %lld\n", p->key.count_neighbour, p->key.i, p->key.ii);
		printf("%lld\n", p->key.i);
#else
		//printf("%d %d %d\n", p->key.count_neighbour, p->key.i, p->key.ii);
		printf("%d\n", p->key.i);
#endif

		print_AVL(p->right);
	}
} // print_AVL

// static - локальная область видимости.
// Возвращает количество узлов в АВЛ дереве.
static integer inumber_AVL_node_gl = 0;

// getnumber_AVL_node_local.
static void getnumber_AVL_node_local(node_AVLST* p)
{
	if (p != nullptr) {
		inumber_AVL_node_gl++;
		getnumber_AVL_node_local(p->left);
		getnumber_AVL_node_local(p->right);
	}
} // getnumber_AVL_node_local.

static integer getnumber_AVL_node_global(node_AVLST* p)
{
	inumber_AVL_node_gl = 0;
	getnumber_AVL_node_local(p);
	integer ir = inumber_AVL_node_gl;
	inumber_AVL_node_gl = 0;
	return ir;
}


// getlist_AVL_node_local.
void getlist_AVL_node_local(node_AVLST* p, hashlist_i*& hl)
{
	if (p != nullptr) {
		inumber_AVL_node_gl++;
		insertion_list_i(hl, p->key.i);
		getlist_AVL_node_local(p->left, hl);
		getlist_AVL_node_local(p->right, hl);
	}
} // getlist_AVL_node_local.

integer getlist_AVL_node_global(node_AVLST* p, hashlist_i*& hl)
{
	inumber_AVL_node_gl = 0;
	getlist_AVL_node_local(p, hl);
	integer ir = inumber_AVL_node_gl;
	inumber_AVL_node_gl = 0;
	return ir;
} // getlist_AVL_node_global

  // тестирование АВЛ дерева
  // 12 декабря 2015 удовлетворительно.
void test_AVLST()
{
	node_AVLST* root = 0;
	data_BalTreeST d3;
	d3.i = rand();
	root = remove_AVL(root, d3);
	for (integer i = 0; i < 10; ++i)
	{
		data_BalTreeST d;
		d.i = rand();
		root = insert_and_modify(root, d, d);
		print_AVL(root);
		if (i == 5) {

			//41 18467 6334


			data_BalTreeST d1;
			d1.i = 18467;
			data_BalTreeST d2;
			d2.i = rand();
			root = insert_and_modify(root, d2, d1);
			print_AVL(root);
			printf("remove 41 18467 6334\n");
			root = remove_AVL(root, d1);
			print_AVL(root);
			printf("found max\n");
			node_AVLST* emax;
			emax = findmax(root);
#if doubleintprecision == 1
			printf("maximum id %lld\n", emax->key.i);
#else
			printf("maximum id %d\n", emax->key.i);
#endif

			emax = nullptr;
			//clear_AVL(root);
			//root = nullptr;
			print_AVL(root);
		}

		system("pause");
	}
	system("pause");
}
// 5.01.2017 Конец реализации леса АВЛ деревьев.

#endif
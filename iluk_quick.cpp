
#pragma once
#ifndef ILUK_QUICK_CPP
#define ILUK_QUICK_CPP 1

// Пока всё неверно. В следующие выходные 16-17 декабря нужно устранить линейный поиск с помощью
// очереди по приоритетам.
// АВЛ двоичное дерево не подходит. Требуется priority_queue 
// Пирамида соединённая с хеш таблицей.
// АВЛ дерево реализовано с ошибками. Нужно быть внимательней.
#include "priority_queue.cpp"


// Поле данных в АВЛ дереве.
struct data_BalTree_for_Saad
{
	// --> high priority --> for operation <,>
	integer  ind, val;
	// countsosed есть key.
};

// Узел АВЛ дерева.
struct node_AVL_for_Saad
{
	data_BalTree_for_Saad key;
	// Высота поддерева с корнем в данном узле.
	unsigned char height;
	node_AVL_for_Saad* left;
	node_AVL_for_Saad* right;
	// Конструктор.
	node_AVL_for_Saad(data_BalTree_for_Saad k) { key = k; left = right = 0; height = 1; }
};

// Работает также и с пустыми деревьями.
// Обёртка поля height
unsigned char height(node_AVL_for_Saad* p)
{
	return p ? p->height : 0;
};

// Вычисляет balance factor заданного узла
// работает только с ненулевыми указателями.
integer bfactor(node_AVL_for_Saad* p)
{
	return height(p->right) - height(p->left);
};

// Восстанавливает корректное значение поля height
// заданного узла (при условии, что значения этого поля 
// в правом и левом дочерних узлах являются корректными).
void fixheight(node_AVL_for_Saad* p)
{
	unsigned char hl = height(p->left);
	unsigned char hr = height(p->right);
	p->height = (hl > hr ? hl : hr) + 1;
};

// Балансировка узлов.
node_AVL_for_Saad* rotateright(node_AVL_for_Saad* p)
{
	// правый поворот вокруг p
	node_AVL_for_Saad* q = p->left;
	p->left = q->right;
	q->right = p;
	fixheight(p);
	fixheight(q);
	return q;
};

// Левый поворот является симметричной копией правого :
node_AVL_for_Saad* rotateleft(node_AVL_for_Saad* q)
{
	// левый поворот вокруг q
	node_AVL_for_Saad* p = q->right;
	q->right = p->left;
	p->left = q;
	fixheight(q);
	fixheight(p);
	return p;
};

// Код выполняющий балансировку сводится к проверке условий и выполнению поворотов
node_AVL_for_Saad* balance(node_AVL_for_Saad* p) // балансировка узла p
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

// Вставка ключей в дерево. Упорядлчение по полю val
// Возвращает новое значение корня АВЛ дерева.
node_AVL_for_Saad* insert(node_AVL_for_Saad* &p, data_BalTree_for_Saad k)
{
	// Вставка ключа k в дерево с корнем p
	if (p == NULL) {
		node_AVL_for_Saad* r1 = NULL;
		r1 = new node_AVL_for_Saad(k);
		if (r1 == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for r1 in insert my_agregat_amg...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
		return r1;
	}

	if (k.val < p->key.val)
		p->left = insert(p->left, k);
	else
		p->right = insert(p->right, k);


	return balance(p);
} // insert

  // Возвращает true если узел найден в дереве
  // поиск по val
  // Not Work
bool isfound(node_AVL_for_Saad* p, data_BalTree_for_Saad k)
{
	if (p == 0) return false; // ненайден.
							  /*
							  */
	if (k.val < p->key.val)
		return isfound(p->left, k);
	else if (k.val > p->key.val)
		return isfound(p->right, k);
	else if (k.ind < p->key.ind)
		return isfound(p->left, k);
	else if (k.ind > p->key.ind)
		return isfound(p->right, k);
	else

	return true; // найден.
}

// Полное удаление бинарного дерева.
void clear_AVL(node_AVL_for_Saad* p)
{
	if (p != 0) {
		clear_AVL(p->left);
		clear_AVL(p->right);
		// удаляем лист.
		delete p;
		p = 0;
	}
} // clear_AVL

  // Удаление узла с занными свойства с сохранением сбалансированности.
node_AVL_for_Saad* findmin(node_AVL_for_Saad* p)
{
	// поиск узла с минимальным ключём в дереве p
	if (p!=NULL) {
		return p->left ? findmin(p->left) : p;
	}
	else {
		// на поиск минимума подан нулевой указатель.
		return 0;
	}

	/*
	// Дерево упорядочено по полю val.
	if (p != NULL) {
	node_AVL_for_Saad* parent = p;
	while (p->left != 0) {
	parent = p;
	p = p->left;
	}

	return p;
	}
	*/
} // findmin

node_AVL_for_Saad* findmax(node_AVL_for_Saad* p)
{
	// поиск узла с минимальным ключём в дереве p
	if (p != NULL) {
		return p->right ? findmax(p->right) : p;
		/*
		#if doubleintprecision == 1
		if (p->right == 0) {
		//printf("%lld %lld %lld\n", p->key.countsosed, p->key.i, p->key.ii);
		//system("pause");
		return p;
		}
		else {
		//printf("%lld %lld %lld\n",p->key.countsosed,p->key.i,p->key.ii);
		findmax(p->right);
		}
		#else
		if (p->right == 0) {
		//printf("%d %d %d\n", p->key.countsosed, p->key.i, p->key.ii);
		//system("pause");
		return p;
		}
		else {
		//printf("%d %d %d\n",p->key.countsosed,p->key.i,p->key.ii);
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

data_BalTree_for_Saad get_max_AVL(node_AVL_for_Saad* p)
{
	// возвращение максимального узла в дереве.
	return p->right ? get_max_AVL(p->right) : p->key;
}

node_AVL_for_Saad* removemin(node_AVL_for_Saad* p)
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
node_AVL_for_Saad* remove_AVL(node_AVL_for_Saad* p, data_BalTree_for_Saad k)
{
	// Отношение порядка определено 
	// для структуры из трёх целых чисел.
	// удаление ключа k из дерева p
	if (p == 0) return 0;
	// Двоичный поиск нужного элемента.
	/*
	if (k.countsosed < p->key.countsosed)
	p->left = remove_AVL(p->left, k);
	else if (k.countsosed>p->key.countsosed)
	p->right = remove_AVL(p->right, k);
	else
	*/
	if (k.ind < p->key.ind)
		p->left = remove_AVL(p->left, k);
	else if (k.ind>p->key.ind)
		p->right = remove_AVL(p->right, k);
	//else if (k.ii < p->key.ii)
	//p->left = remove_AVL(p->left, k);
	//else if (k.ii>p->key.ii)
	//p->right = remove_AVL(p->right, k);
	else // k==p->key
	{
		node_AVL_for_Saad* q = p->left;
		node_AVL_for_Saad* r = p->right;
		delete p;
		p = 0;
		if (r == 0) return q;
		node_AVL_for_Saad* min = findmin(r);
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
node_AVL_for_Saad* insert_and_modify(node_AVL_for_Saad* p, data_BalTree_for_Saad k, data_BalTree_for_Saad k_search)
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
void print_AVL(node_AVL_for_Saad* p)
{
	if (p != 0) {
		print_AVL(p->left);
		for (integer i = 0; i <= p->height; i++) {
			printf(" ");
		}
#if doubleintprecision == 1
		//printf("%lld %lld %lld\n", p->key.countsosed, p->key.i, p->key.ii);
		printf("%lld %lld\n", p->key.val, p->key.ind);
#else
		//printf("%d %d %d\n", p->key.countsosed, p->key.i, p->key.ii);
		printf("%d %d\n", p->key.val, p->key.ind);
#endif

		print_AVL(p->right);
	}
} // print_AVL

  // тестирование АВЛ дерева
  // 12 декабря 2015 удовлетворительно.
void test_AVL_for_Saad()
{
	node_AVL_for_Saad* root = 0;
	data_BalTree_for_Saad d3;
	d3.val = rand();
	d3.ind = rand();
	//d3.ii = rand();
	root = remove_AVL(root, d3);
	for (integer i = 0; i < 10; i++)
	{
		data_BalTree_for_Saad d;
		d.val = rand();
		d.ind = rand();
		//d.ii = rand();
		root = insert_and_modify(root, d, d);
		print_AVL(root);
		if (i == 5) {

			//41 18467 6334


			data_BalTree_for_Saad d1;
			d1.val = 41;
			d1.ind = 18467;
			//d1.ii = 6334;
			data_BalTree_for_Saad d2;
			d2.val = rand();
			d2.ind = rand();
			//d2.ii = rand();
			root = insert_and_modify(root, d2, d1);
			print_AVL(root);
			printf("remove 41 18467 6334\n");
			root = remove_AVL(root, d1);
			print_AVL(root);
			printf("found max\n");
			node_AVL_for_Saad* emax;
			emax = findmax(root);
#if doubleintprecision == 1
			printf("maximum id %lld\n", emax->key.val);
#else
			printf("maximum id %d\n", emax->key.val);
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


// Очередь по приоритетам в соединении с хеш таблицей.
//PQ(integer maxN, integer max_key_size); // конструктор
//~PQ(); // деструктор
// Есть ли элемент с данным ключём в таблице ?
//bool isfound(integer key);
//bool empty() const; // проверка на пустоту.
// Очищаем содержимое и она снова готова к использованию.
//void clear();
// Вернуть элемент с заданным ключём:
// Обязательно предполагается что ключ существует внутри таблицы.
//Item get(integer key);
// прочитать максимальный элемент.
//Item readmax();
// Прочитать значение ключа максимального элемента.
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
//  Меняет местами значения элементов с ключами i и j. 
// Сохраняет порядок кучи.
//void exchange(integer i, integer j);



// Медленный линейный поиск ликвидирован.
/* ----------------------------------------------------------------------- */
/* Subroutine */ integer iluk_quick_stable(integer n, doublereal* &a, integer* &ja, integer* &ia,
	integer lfil, doublereal* &alu, integer* &jlu, integer* &ju,
	integer* &levs, integer iwk, doublereal* &w, integer* &jw, integer &ierr)
{
	/* System generated locals */
	integer i__1, i__2, i__3, i__4;

	/* Local variables */
	integer i__, j, k;
	doublereal s, t;
	integer j1, j2, n2, ii, jj, ju0;
	doublereal fact;
	integer lenl, jlev, lenu, jpos, jrow;

	// Переменные для проверки корректности исходного кода.
	integer jrow1=-1, jrow2=-1, k1=-1, k2=-1;

	/* ----------------------------------------------------------------------* */
	/*     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) * */
	/* ----------------------------------------------------------------------* */

	/* on entry: */
	/* ========== */
	/* n       = integer. The row dimension of the matrix A. The matrix */

	/* a,ja,ia = matrix stored in Compressed Sparse Row format. */

	/* lfil    = integer. The fill-in parameter. Each element whose */
	/*           leve-of-fill exceeds lfil during the ILU process is dropped. */
	/*           lfil must be .ge. 0 */

	/* tol     = real*8. Sets the threshold for dropping small terms in the */
	/*           factorization. See below for details on dropping strategy. */

	/* iwk     = integer. The minimum length of arrays alu, jlu, and levs. */

	/* On return: */
	/* =========== */

	/* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing */
	/*           the L and U factors together. The diagonal (stored in */
	/*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix */
	/*           contains the i-th row of L (excluding the diagonal entry=1) */
	/*           followed by the i-th row of U. */

	/* ju      = integer array of length n containing the pointers to */
	/*           the beginning of each row of U in the matrix alu,jlu. */

	/* levs    = integer (work) array of size iwk -- which contains the */
	/*           levels of each element in alu, jlu. */

	/* ierr    = integer. Error message with the following meaning. */
	/*           ierr  = 0    --> successful return. */
	/*           ierr .gt. 0  --> zero pivot encountered at step number ierr. */
	/*           ierr  = -1   --> Error. input matrix may be wrong. */
	/*                            (The elimination process has generated a */
	/*                            row in L or U whose length is .gt.  n.) */
	/*           ierr  = -2   --> The matrix L overflows the array al. */
	/*           ierr  = -3   --> The matrix U overflows the array alu. */
	/*           ierr  = -4   --> Illegal value for lfil. */
	/*           ierr  = -5   --> zero row encountered in A or U. */

	/* work arrays: */
	/* ============= */
	/* jw      = integer work array of length 3*n. */
	/* w       = real work array of length n */

	/* Notes/known bugs: This is not implemented efficiently storage-wise. */
	/*       For example: Only the part of the array levs(*) associated with */
	/*       the U-matrix is needed in the routine.. So some storage can */
	/*       be saved if needed. The levels of fills in the LU matrix are */
	/*       output for information only -- they are not needed by LU-solve. */

	/* ---------------------------------------------------------------------- */
	/* w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] */
	/* jw(n+1:2n)  stores the nonzero indicator. */

	/* Notes: */
	/* ------ */
	/* All the diagonal elements of the input matrix must be  nonzero. */

	/* ----------------------------------------------------------------------* */
	/*     locals */
	/* Parameter adjustments */
	--jw;
	--w;
	--ju;
	--ia;
	--a;
	--ja;
	--alu;
	--jlu;
	--levs;

	// Целочисленная очередь по приоритетам.
	PQ<integer>  pq(n + 1, n + 1);

	/* Function Body */
	if (lfil < 0) {
		goto L998;
	}
	/* ----------------------------------------------------------------------- */
	/*     initialize ju0 (points to next element to be added to alu,jlu) */
	/*     and pointer array. */
	/* ----------------------------------------------------------------------- */
	n2 = n + n;
	ju0 = n + 2;
	jlu[1] = ju0;

	/*     initialize nonzero indicator array + levs array -- */

	i__1 = n << 1;
	for (j = 1; j <= i__1; ++j) {
		jw[j] = 0;
		/* L1: */
	}

	

	/* ----------------------------------------------------------------------- */
	/*     beginning of main loop. */
	/* ----------------------------------------------------------------------- */
	i__1 = n;
	for (ii = 1; ii <= i__1; ++ii) {

		

		j1 = ia[ii];
		j2 = ia[ii + 1] - 1;

		/*     unpack L-part and U-part of row of A in arrays w */

		lenu = 1;
		lenl = 0;
		jw[ii] = ii;

		

		w[ii] = 0.f;
		jw[n + ii] = ii;

		i__2 = j2;
		for (j = j1; j <= i__2; ++j) {
			k = ja[j];
			t = a[j];
			if (t == 0.f) {
				goto L170;
			}
			if (k < ii) {
				++lenl;
				jw[lenl] = k;

				// push jw
				pq.insert(-k, lenl);
							

				w[lenl] = t;
				jw[n2 + lenl] = 0;
				jw[n + k] = lenl;
			}
			else if (k == ii) {
				w[ii] = t;
				jw[n2 + ii] = 0;
			}
			else {
				++lenu;
				jpos = ii + lenu - 1;
				jw[jpos] = k;
				w[jpos] = t;
				jw[n2 + jpos] = 0;
				jw[n + k] = jpos;
			}
		L170:
			;
		}

		jj = 0;

		/*     eliminate previous rows */

	L150:

		// remove in jw
		pq.remove(jj);

		++jj;
		if (jj > lenl) {
			goto L160;
		}
		/* ----------------------------------------------------------------------- */
		/*     in order to do the elimination in the correct order we must select */
		/*     the smallest column index among jw(k), k=jj+1, ..., lenl. */
		/* ----------------------------------------------------------------------- */
		jrow = jw[jj];
		k = jj;

		/*     determine smallest column index */

		i__2 = lenl;
		// Дьявольски медленный поиск минимума. Это линейный поиск.
		//printf("jj=%d\n",jj);// jj==1 далеко не всегда.
		// Это означает, что нужно поддерживать удаление элемента по ключу.
		if (0) {
			// Чрезвычайно медленный линейный поиск.
			for (j = jj + 1; j <= i__2; ++j) {
				if (jw[j] < jrow) {
					jrow = jw[j];
					k = j;
				}
			}
		}
		else {
			if (0) {

				// Проверочный участок кода.
				// Раскоментировать если нужно проверить код.

				jrow1 = jrow;
				k1 = k;

				for (j = jj + 1; j <= i__2; ++j) {
					if (jw[j] < jrow1) {
						jrow1 = jw[j];
						k1 = j;
					}
				}

				if (jj + 1 <= i__2) {
					// remove in jw
					pq.remove(jj);
					integer jrow3 = jrow, k3 = k;
					// на основе очереди по приоритетам в сочетании с хеш таблицей.
					if (jrow > (-pq.readmax())) {
						jrow = -pq.readmax();
						k = pq.readkeymaxelm();
					}
					pq.insert(-jrow3, k3);
				}


				jrow2 = jrow;
				k2 = k;

				if ((jrow1 != jrow2) || (k1 != k2)) {
					for (j = jj; j <= i__2; ++j) {
						printf("jw[%lld]=%lld ", j, jw[j]);
					}
					printf("\n");
					pq.print_log('s');
					printf("k1=%lld k2=%lld jrow1=%lld jrow2=%lld\n", k1, k2, jrow1, jrow2);
					getchar();
				}

			}
			else {
				jrow = -pq.readmax();
				k = pq.readkeymaxelm();
				/*
				if (jj + 1 <= i__2) {
					// remove in jw
					pq.remove(jj);
					integer jrow3 = jrow, k3 = k;
					// на основе очереди по приоритетам в сочетании с хеш таблицей.
					if (jrow > (-pq.readmax())) {
						jrow = -pq.readmax();
						k = pq.readkeymaxelm();
					}
					pq.insert(-jrow3, k3);
				}
				*/
			}

		}

		if (k != jj) {
			/*     exchange in jw */
			j = jw[jj];
			jw[jj] = jw[k];
			jw[k] = j;
			/*     exchange in jw(n+  (pointers/ nonzero indicator). */
			jw[n + jrow] = jj;
			jw[n + j] = k;
			/*     exchange in jw(n2+  (levels) */
			j = jw[n2 + jj];
			jw[n2 + jj] = jw[n2 + k];
			jw[n2 + k] = j;
			/*     exchange in w */
			s = w[jj];
			w[jj] = w[k];
			w[k] = s;

			
			pq.exchange_speshial_for_Saad(jj, k);
		}

		/*     zero out element in row by resetting jw(n+jrow) to zero. */

		jw[n + jrow] = 0;

		/*     get the multiplier for row to be eliminated (jrow) + its level */

		fact = w[jj] * alu[jrow];
		jlev = jw[n2 + jj];
		if (jlev > lfil) {
			goto L150;
		}

		/*     combine current row and row jrow */

		i__2 = jlu[jrow + 1] - 1;
		for (k = ju[jrow]; k <= i__2; ++k) {
			s = fact * alu[k];
			j = jlu[k];
			jpos = jw[n + j];
			if (j >= ii) {

				/*     dealing with upper part. */

				if (jpos == 0) {

					/*     this is a fill-in element */

					++lenu;
					if (lenu > n) {
						goto L995;
					}
					i__ = ii + lenu - 1;
					jw[i__] = j;
					jw[n + j] = i__;
					w[i__] = -s;
					jw[n2 + i__] = jlev + levs[k] + 1;
				}
				else {

					/*     this is not a fill-in element */

					w[jpos] -= s;
					/* Computing MIN */
					i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
					jw[n2 + jpos] = min(i__3, i__4);
				}
			}
			else {

				/*     dealing with lower part. */

				if (jpos == 0) {

					/*     this is a fill-in element */

					++lenl;
					if (lenl > n) {
						goto L995;
					}

					// push jw
					pq.insert(-j, lenl);

					jw[lenl] = j;
					jw[n + j] = lenl;
					w[lenl] = -s;
					jw[n2 + lenl] = jlev + levs[k] + 1;
				}
				else {

					/*     this is not a fill-in element */

					w[jpos] -= s;
					/* Computing MIN */
					i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
					jw[n2 + jpos] = min(i__3, i__4);
				}
			}
			/* L203: */
		}
		w[jj] = fact;

		//pq.remove(jj);
		//pq.insert(-jrow, jj);

		jw[jj] = jrow;
		goto L150;
	L160:

		/*     reset double-pointer to zero (U-part) */

		i__2 = lenu;
		for (k = 1; k <= i__2; ++k) {
			jw[n + jw[ii + k - 1]] = 0;
			/* L308: */
		}

		/*     update l-matrix */

		i__2 = lenl;
		for (k = 1; k <= i__2; ++k) {
			if (ju0 > iwk) {
				goto L996;
			}
			if (jw[n2 + k] <= lfil) {
				alu[ju0] = w[k];
				jlu[ju0] = jw[k];
				++ju0;
			}
			/* L204: */
		}

		/*     save pointer to beginning of row ii of U */

		ju[ii] = ju0;

		/*     update u-matrix */

		i__2 = ii + lenu - 1;
		for (k = ii + 1; k <= i__2; ++k) {
			if (jw[n2 + k] <= lfil) {
				jlu[ju0] = jw[k];
				alu[ju0] = w[k];
				levs[ju0] = jw[n2 + k];
				++ju0;
			}
			/* L302: */
		}
		if (fabs(w[ii]) < 1.0e-30) {
			printf("w[%lld]=%e\n", ii, w[ii]);
			w[ii] = 1.0;
			printf("k1=%lld k2=%lld jrow1=%lld jrow2=%lld\n",k1,k2,jrow1,jrow2);
			getchar();
			//goto L999;
		}

		alu[ii] = 1.0 / w[ii];

		/*     update pointer to beginning of next row of U. */

		jlu[ii + 1] = ju0;
		/* ----------------------------------------------------------------------- */
		/*     end main loop */
		/* ----------------------------------------------------------------------- */
		/* L500: */
	}

	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;
	

	ierr = 0;
	return 0;

	/*     incomprehensible error. Matrix must be wrong. */

L995:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;
	

	ierr = -1;
	return 0;

	/*     insufficient storage in L. */

L996:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;
	

	ierr = -2;
	return 0;

	/*     insufficient storage in U. */

	/* L997: */
	// ierr = -3;
	// return 0;

	/*     illegal lfil entered. */

L998:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;
	

	ierr = -4;
	return 0;

	/*     zero row encountered in A or U. */

	/*
L999:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;
	

	ierr = -5;
	return 0;
	*/
	/* ----------------end-of-iluk-------------------------------------------- */
	/* ----------------------------------------------------------------------- */
} /* iluk_quick */


  // Моя версия в которую добавлен поиск минимального элемента в АВЛ дереве. 10.12.2017
  // Медленный линейный поиск ликвидирован.
  /* ----------------------------------------------------------------------- */
/* Subroutine */ integer iluk_quickerr(integer n, doublereal* &a, integer* &ja, integer* &ia,
	integer lfil, doublereal* &alu, integer* &jlu, integer* &ju,
	integer* &levs, integer iwk, doublereal* &w, integer* &jw, integer &ierr)
{
	/* System generated locals */
	integer i__1, i__2, i__3, i__4;

	/* Local variables */
	integer i__, j, k;
	doublereal s, t;
	integer j1, j2, n2, ii, jj, ju0;
	doublereal fact;
	integer lenl, jlev, lenu, jpos, jrow;

	/* ----------------------------------------------------------------------* */
	/*     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) * */
	/* ----------------------------------------------------------------------* */

	/* on entry: */
	/* ========== */
	/* n       = integer. The row dimension of the matrix A. The matrix */

	/* a,ja,ia = matrix stored in Compressed Sparse Row format. */

	/* lfil    = integer. The fill-in parameter. Each element whose */
	/*           leve-of-fill exceeds lfil during the ILU process is dropped. */
	/*           lfil must be .ge. 0 */

	/* tol     = real*8. Sets the threshold for dropping small terms in the */
	/*           factorization. See below for details on dropping strategy. */

	/* iwk     = integer. The minimum length of arrays alu, jlu, and levs. */

	/* On return: */
	/* =========== */

	/* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing */
	/*           the L and U factors together. The diagonal (stored in */
	/*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix */
	/*           contains the i-th row of L (excluding the diagonal entry=1) */
	/*           followed by the i-th row of U. */

	/* ju      = integer array of length n containing the pointers to */
	/*           the beginning of each row of U in the matrix alu,jlu. */

	/* levs    = integer (work) array of size iwk -- which contains the */
	/*           levels of each element in alu, jlu. */

	/* ierr    = integer. Error message with the following meaning. */
	/*           ierr  = 0    --> successful return. */
	/*           ierr .gt. 0  --> zero pivot encountered at step number ierr. */
	/*           ierr  = -1   --> Error. input matrix may be wrong. */
	/*                            (The elimination process has generated a */
	/*                            row in L or U whose length is .gt.  n.) */
	/*           ierr  = -2   --> The matrix L overflows the array al. */
	/*           ierr  = -3   --> The matrix U overflows the array alu. */
	/*           ierr  = -4   --> Illegal value for lfil. */
	/*           ierr  = -5   --> zero row encountered in A or U. */

	/* work arrays: */
	/* ============= */
	/* jw      = integer work array of length 3*n. */
	/* w       = real work array of length n */

	/* Notes/known bugs: This is not implemented efficiently storage-wise. */
	/*       For example: Only the part of the array levs(*) associated with */
	/*       the U-matrix is needed in the routine.. So some storage can */
	/*       be saved if needed. The levels of fills in the LU matrix are */
	/*       output for information only -- they are not needed by LU-solve. */

	/* ---------------------------------------------------------------------- */
	/* w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] */
	/* jw(n+1:2n)  stores the nonzero indicator. */

	/* Notes: */
	/* ------ */
	/* All the diagonal elements of the input matrix must be  nonzero. */

	/* ----------------------------------------------------------------------* */
	/*     locals */
	/* Parameter adjustments */
	--jw;
	--w;
	--ju;
	--ia;
	--a;
	--ja;
	--alu;
	--jlu;
	--levs;

	node_AVL_for_Saad* root = 0; // АВЛ дерево для быстрого поиска минимума.

	/* Function Body */
	if (lfil < 0) {
		goto L998;
	}
	/* ----------------------------------------------------------------------- */
	/*     initialize ju0 (points to next element to be added to alu,jlu) */
	/*     and pointer array. */
	/* ----------------------------------------------------------------------- */
	n2 = n + n;
	ju0 = n + 2;
	jlu[1] = ju0;

	/*     initialize nonzero indicator array + levs array -- */

	i__1 = n << 1;
	for (j = 1; j <= i__1; ++j) {
		jw[j] = 0;
		/* L1: */
	}

	

								 /* ----------------------------------------------------------------------- */
								 /*     beginning of main loop. */
								 /* ----------------------------------------------------------------------- */
	i__1 = n;
	for (ii = 1; ii <= i__1; ++ii) {

		clear_AVL(root);

		j1 = ia[ii];
		j2 = ia[ii + 1] - 1;

		/*     unpack L-part and U-part of row of A in arrays w */

		lenu = 1;
		lenl = 0;
		jw[ii] = ii;

		// push jw
		data_BalTree_for_Saad ddel;
		ddel.ind = ii;
		ddel.val = ii;
		data_BalTree_for_Saad dsearch = ddel;
		root = insert_and_modify(root, ddel, dsearch);

		w[ii] = 0.f;
		jw[n + ii] = ii;

		i__2 = j2;
		for (j = j1; j <= i__2; ++j) {
			k = ja[j];
			t = a[j];
			if (t == 0.f) {
				goto L170;
			}
			if (k < ii) {
				++lenl;
				jw[lenl] = k;

				// push jw
				ddel.ind = lenl;
				ddel.val = k;
				dsearch = ddel;
				root = insert_and_modify(root, ddel, dsearch);

				w[lenl] = t;
				jw[n2 + lenl] = 0;
				jw[n + k] = lenl;
			}
			else if (k == ii) {
				w[ii] = t;
				jw[n2 + ii] = 0;
			}
			else {
				++lenu;
				jpos = ii + lenu - 1;
				jw[jpos] = k;
				w[jpos] = t;
				jw[n2 + jpos] = 0;
				jw[n + k] = jpos;
			}
		L170:
			;
		}

		jj = 0;

		/*     eliminate previous rows */

	L150:

		// remove in jw
		ddel.ind = jj;
		remove_AVL(root, ddel);

		node_AVL_for_Saad* emin = 0; // инициализация.

		++jj;
		if (jj > lenl) {
			goto L160;
		}
		/* ----------------------------------------------------------------------- */
		/*     in order to do the elimination in the correct order we must select */
		/*     the smallest column index among jw(k), k=jj+1, ..., lenl. */
		/* ----------------------------------------------------------------------- */
		// Присваивание для jrow и k будет ниже после линейного поиска с использованием очереди по приоритетам.
		// Раскоментировать при линейном поиске.
		//jrow = jw[jj];
		//k = jj;

		/*     determine smallest column index */

		i__2 = lenl;
		// Дьявольски медленный поиск минимума. Это линейный поиск.
		//printf("jj=%d\n",jj);// jj==1 далеко не всегда.
		// Это означает что нужно поддерживать удаление элемента по ключу.
		//getchar();
		/*
		for (j = jj + 1; j <= i__2; ++j) {
		if (jw[j] < jrow) {
		jrow = jw[j];
		k = j;
		}
		}
		*/
		
		emin = findmin(root);
		jrow = emin->key.val;
		k = emin->key.ind;
		if (k != jj) {
			ddel.ind = k;
			ddel.val = jrow;
			remove_AVL(root, ddel);
		}
		emin = 0;

		if (k != jj) {
			/*     exchange in jw */
			j = jw[jj];
			jw[jj] = jw[k];
			jw[k] = j;
			/*     exchange in jw(n+  (pointers/ nonzero indicator). */
			jw[n + jrow] = jj;
			jw[n + j] = k;
			/*     exchange in jw(n2+  (levels) */
			j = jw[n2 + jj];
			jw[n2 + jj] = jw[n2 + k];
			jw[n2 + k] = j;
			/*     exchange in w */
			s = w[jj];
			w[jj] = w[k];
			w[k] = s;
		}

		/*     zero out element in row by resetting jw(n+jrow) to zero. */

		jw[n + jrow] = 0;

		/*     get the multiplier for row to be eliminated (jrow) + its level */

		fact = w[jj] * alu[jrow];
		jlev = jw[n2 + jj];
		if (jlev > lfil) {
			goto L150;
		}

		/*     combine current row and row jrow */

		i__2 = jlu[jrow + 1] - 1;
		for (k = ju[jrow]; k <= i__2; ++k) {
			s = fact * alu[k];
			j = jlu[k];
			jpos = jw[n + j];
			if (j >= ii) {

				/*     dealing with upper part. */

				if (jpos == 0) {

					/*     this is a fill-in element */

					++lenu;
					if (lenu > n) {
						goto L995;
					}
					i__ = ii + lenu - 1;
					jw[i__] = j;
					jw[n + j] = i__;
					w[i__] = -s;
					jw[n2 + i__] = jlev + levs[k] + 1;
				}
				else {

					/*     this is not a fill-in element */

					w[jpos] -= s;
					/* Computing MIN */
					i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
					jw[n2 + jpos] = min(i__3, i__4);
				}
			}
			else {

				/*     dealing with lower part. */

				if (jpos == 0) {

					/*     this is a fill-in element */

					++lenl;
					if (lenl > n) {
						goto L995;
					}

					// push jw
					ddel.ind = lenl;
					ddel.val = j;
					dsearch = ddel;
					root = insert_and_modify(root, ddel, dsearch);

					jw[lenl] = j;
					jw[n + j] = lenl;
					w[lenl] = -s;
					jw[n2 + lenl] = jlev + levs[k] + 1;
				}
				else {

					/*     this is not a fill-in element */

					w[jpos] -= s;
					/* Computing MIN */
					i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
					jw[n2 + jpos] = min(i__3, i__4);
				}
			}
			/* L203: */
		}
		w[jj] = fact;

		// push jw
		ddel.ind = jj;
		ddel.val = jrow;
		dsearch = ddel;
		root = insert_and_modify(root, ddel, dsearch);


		jw[jj] = jrow;
		goto L150;
	L160:

		/*     reset double-pointer to zero (U-part) */

		i__2 = lenu;
		for (k = 1; k <= i__2; ++k) {
			jw[n + jw[ii + k - 1]] = 0;
			/* L308: */
		}

		/*     update l-matrix */

		i__2 = lenl;
		for (k = 1; k <= i__2; ++k) {
			if (ju0 > iwk) {
				goto L996;
			}
			if (jw[n2 + k] <= lfil) {
				alu[ju0] = w[k];
				jlu[ju0] = jw[k];
				++ju0;
			}
			/* L204: */
		}

		/*     save pointer to beginning of row ii of U */

		ju[ii] = ju0;

		/*     update u-matrix */

		i__2 = ii + lenu - 1;
		for (k = ii + 1; k <= i__2; ++k) {
			if (jw[n2 + k] <= lfil) {
				jlu[ju0] = jw[k];
				alu[ju0] = w[k];
				levs[ju0] = jw[n2 + k];
				++ju0;
			}
			/* L302: */
		}
		if (fabs(w[ii]) < 1.0e-30) {
			printf("w[%lld]=%e\n", ii, w[ii]);
			goto L999;
		}

		alu[ii] = 1.0 / w[ii];

		/*     update pointer to beginning of next row of U. */

		jlu[ii + 1] = ju0;
		/* ----------------------------------------------------------------------- */
		/*     end main loop */
		/* ----------------------------------------------------------------------- */
		/* L500: */
	}

	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;

	ierr = 0;
	return 0;

	/*     incomprehensible error. Matrix must be wrong. */

L995:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;

	ierr = -1;
	return 0;

	/*     insufficient storage in L. */

L996:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;

	ierr = -2;
	return 0;

	/*     insufficient storage in U. */

	/* L997: */
	// ierr = -3;
	// return 0;

	/*     illegal lfil entered. */

L998:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;

	ierr = -4;
	return 0;

	/*     zero row encountered in A or U. */

L999:
	++jw;
	++w;
	++ju;
	++ia;
	++a;
	++ja;
	++alu;
	++jlu;
	++levs;

	ierr = -5;
	return 0;
	/* ----------------end-of-iluk-------------------------------------------- */
	/* ----------------------------------------------------------------------- */
} /* iluk_quick */

#endif // !ILUK_QUICK_CPP
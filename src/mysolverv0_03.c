// Файл mysolverv0_03.c 
// решает одно уравнение.

#pragma once
#ifndef MY_SOLVER_v0_03_C
#define MY_SOLVER_v0_03_C 1

// Based on Renormalization Group Theory turbulence model.
#include "my_RNG_LES.cpp"
#include "test_filtr.cpp" // всё для применения двойной фильтрации.

// Проверено 08,08,2020
// index_of(i_1,'x');
integer index_of(integer ibase, char ch) {
	// 0 << ibase < t.maxnod
	// 0 <= return < 3 * t.maxnod
	switch (ch) {
	case 'x': return (3*ibase + 0); break; // X
	case 'y': return (3*ibase + 1); break; // Y
	case 'z': return (3*ibase + 2); break; // Z
	case 'X': return (3*ibase + 0); break; // X
	case 'Y': return (3*ibase + 1); break; // Y
	case 'Z': return (3*ibase + 2); break; // Z
	default: 
		printf("ERROR !!! index_of 08.08.2020\n");
		system("PAUSE");
		exit(1);
		break; 
	}
}// index_of

/*
// index_of(i_1,'x');
integer index_of(integer ibase, char ch) {
	// 0 <= return < 3 * t.maxnod
	switch (ch) {
	case 'x': return (3 * ibase + 0); break; // X
	case 'y': return (3 * ibase + 1); break; // Y
	case 'z': return (3 * ibase + 2); break; // Z
	case 'X': return (3 * ibase + 0); break; // X
	case 'Y': return (3 * ibase + 1); break; // Y
	case 'Z': return (3 * ibase + 2); break; // Z
	default: return (3 * ibase + 0); break; // X
	}
}// index_of
*/

/*
integer index_of(integer ibase, char ch) {
	switch (ch) {
	case 'x': return (3 * ibase); break; // X
	case 'y': return (3 * ibase+2); break; // Y
	case 'z': return (3 * ibase+1); break; // Z
	case 'X': return (3 * ibase); break; // X
	case 'Y': return (3 * ibase + 2); break; // Y
	case 'Z': return (3 * ibase + 1); break; // Z
	}
}// index_of
*/



// Данная функция используется для отладки.
// Она печатает матрицу СЛАУ для уравнения теплопередачи.
void print_temp_slau(TEMPER& t) {
	FILE* fptslau = nullptr; // файл в который будут записываться коэффициенты матрицы СЛАУ.

	// создание файла для записи значений собранной матрицы теплопроводности.
#ifdef MINGW_COMPILLER
	int err = 0;
	fptslau = fopen64("temperslau.txt", "w");
	if (fptslau == NULL) {
		err = 1;
    }
#else
	errno_t err = 0;
	err = fopen_s(&fptslau, "temperslau.txt", "w");
#endif
	
	if ((err == 0)&&(fptslau!=nullptr)) {
		if (fptslau != nullptr) {
#if doubleintprecision == 1
			fprintf(fptslau, "SLAU maxelm=%lld, maxbound=%lld, maxelm+maxbound=%lld\n", t.maxelm, t.maxbound, t.maxelm + t.maxbound);

			for (integer i = 0; i < t.maxelm; i++) {
				fprintf(fptslau, "id=%lld ap=%e ae=%e aw=%e an=%e as=%e at=%e ab=%e b=%e\n", i, t.slau[i].ap, t.slau[i].ae, t.slau[i].aw, t.slau[i].an, t.slau[i].as, t.slau[i].at, t.slau[i].ab, t.slau[i].b);
		    }

			fprintf(fptslau, "BOUNDARY CONDITION...\n");

			for (integer i = 0; i < t.maxbound; i++) {
				fprintf(fptslau, "id=%lld aw=%e ai=%e b=%e\n", i + t.maxelm, t.slau_bon[i].aw, t.slau_bon[i].ai, t.slau_bon[i].b);
			}
#else
			fprintf(fptslau, "SLAU maxelm=%d, maxbound=%d, maxelm+maxbound=%d\n", t.maxelm, t.maxbound, t.maxelm + t.maxbound);

			for (integer i = 0; i < t.maxelm; i++) {
				fprintf(fptslau, "id=%d ap=%e ae=%e aw=%e an=%e as=%e at=%e ab=%e b=%e\n", i, t.slau[i].ap, t.slau[i].ae, t.slau[i].aw, t.slau[i].an, t.slau[i].as, t.slau[i].at, t.slau[i].ab, t.slau[i].b);
			}

			fprintf(fptslau, "BOUNDARY CONDITION...\n");

			for (integer i = 0; i < t.maxbound; i++) {
				fprintf(fptslau, "id=%d aw=%e ai=%e b=%e\n", i + t.maxelm, t.slau_bon[i].aw, t.slau_bon[i].ai, t.slau_bon[i].b);
			}
#endif

			

			fclose(fptslau); // закрытие файла
		}
	}
	
} // print_temp_slau

typedef struct TSortElm {
	TOCHKA p;
	integer i;

	TSortElm() {
		//TOCHKA p;
	    i=-1;
	}
} SortElm;

// Сортировка Дональда Шелла из Роберта Седжвика 1959год.
void ShellSort(SortElm* &rb, integer in) {
	integer i, j;
	SortElm x;
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

			while (j >= h && ((x.p.x < rb[j - h].p.x)||((x.p.x == rb[j - h].p.x)&&(x.p.y < rb[j - h].p.y))||((x.p.x == rb[j - h].p.x) && (x.p.y == rb[j - h].p.y)&&(x.p.z < rb[j - h].p.z)))) {
				rb[j] = rb[j - h]; j -= h;
			}
			rb[j] = x;
		}
	}
} // ShellSort[1959]

// Заменяем линейный поиск на двоичный поиск и одну сортировку.
// Самой медленной операцией становится вставка в отсортированный массив.
bool BinarySearch(SortElm* &pa_global, TOCHKA &key, integer &i, integer n) {
	integer left = 0, right = n, mid;
	while (left <= right)
	{
		mid = left + (right - left) / 2;
		if ((key.x<pa_global[mid].p.x)||((key.x==pa_global[mid].p.x)&&(key.y<pa_global[mid].p.y))||((key.x == pa_global[mid].p.x)&&(key.y == pa_global[mid].p.y)&&(key.z < pa_global[mid].p.z))) right = mid - 1;
		else if ((key.x>pa_global[mid].p.x) || ((key.x == pa_global[mid].p.x) && (key.y>pa_global[mid].p.y)) || ((key.x == pa_global[mid].p.x) && (key.y == pa_global[mid].p.y) && (key.z > pa_global[mid].p.z)))  left = mid + 1;
		else {
			i = mid;
			return true;
		}
	}
	i = left;
	return false;
} // BinarySearch

// АВЛ дерево начало.
// Узел АВЛ дерева.
struct node_AVL1
{
	SortElm key;
	// Высота поддерева с корнем в данном узле.
	unsigned char height;
	node_AVL1* left;
	node_AVL1* right;
	// Конструктор.
	node_AVL1(SortElm k) { key = k; left = right = 0; height = 1; }
};

// Работает также и с пустыми деревьями.
// Обёртка поля height
unsigned char height(node_AVL1* p)
{
	return p ? p->height : 0;
};

// Вычисляет balance factor заданного узла
// работает только с ненулевыми указателями.
integer bfactor(node_AVL1* p)
{
	return height(p->right) - height(p->left);
};

// Восстанавливает корректное значение поля height
// заданного узла (при условии, что значения этого поля 
// в правом и левом дочерних узлах являются корректными).
void fixheight(node_AVL1* p)
{
	unsigned char hl = height(p->left);
	unsigned char hr = height(p->right);
	p->height = (hl > hr ? hl : hr) + 1;
};

// Балансировка узлов.
node_AVL1* rotateright(node_AVL1* p)
{
	// правый поворот вокруг p
	node_AVL1* q = p->left;
	p->left = q->right;
	q->right = p;
	fixheight(p);
	fixheight(q);
	return q;
};

// Левый поворот является симметричной копией правого:
node_AVL1* rotateleft(node_AVL1* q)
{
	// левый поворот вокруг q
	node_AVL1* p = q->right;
	q->right = p->left;
	p->left = q;
	fixheight(q);
	fixheight(p);
	return p;
};

// Код выполняющий балансировку сводится к проверке условий и выполнению поворотов
node_AVL1* balance(node_AVL1* p) // балансировка узла p
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
node_AVL1* insert(node_AVL1* &p, SortElm k)
{
	// Вставка ключа k в дерево с корнем p
	if (p == nullptr) {
		node_AVL1* r1 = nullptr;
		r1 = new node_AVL1(k);
		if (r1 == nullptr) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for r1 in insert AVL mysolverv0_03...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
		return r1;
	}
	if ((k.p.x<p->key.p.x) || ((k.p.x == p->key.p.x) && (k.p.y<p->key.p.y)) || ((k.p.x == p->key.p.x) && (k.p.y == p->key.p.y) && (k.p.z <  p->key.p.z)))
		p->left = insert(p->left, k);
	else 
		p->right = insert(p->right, k);
	
	return balance(p);
} // insert

  // Возвращает true если узел найден в дереве
bool isfound(node_AVL1* p, SortElm &k)
{
	if (p == 0) return false; // не найден.
	if ((k.p.x<p->key.p.x) || ((k.p.x == p->key.p.x) && (k.p.y<p->key.p.y)) || ((k.p.x == p->key.p.x) && (k.p.y == p->key.p.y) && (k.p.z <  p->key.p.z)))
		return isfound(p->left, k);
	else if ((k.p.x>p->key.p.x) || ((k.p.x == p->key.p.x) && (k.p.y>p->key.p.y)) || ((k.p.x == p->key.p.x) && (k.p.y == p->key.p.y) && (k.p.z >  p->key.p.z)))
		return isfound(p->right, k);
	else {
		k.i = p->key.i;
		return true; // найден.
	}
}

// Полное удаление бинарного дерева.
void clear_AVL(node_AVL1* p)
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
node_AVL1* findmin(node_AVL1* p)
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

node_AVL1* findmax(node_AVL1* p)
{
	// поиск узла с максимальным ключом в дереве p
	if (p != 0) {
		return p->right ? findmax(p->right) : p;
	}
	else {
		// На поиск максимума подан нулевой указатель.
		return 0;
	}
} // findmax

SortElm get_max_AVL(node_AVL1* p)
{
	// возвращение максимального узла в дереве.
	return p->right ? get_max_AVL(p->right) : p->key;
}

node_AVL1* removemin(node_AVL1* p)
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
node_AVL1* remove_AVL(node_AVL1* p, SortElm k)
{
	// Отношение порядка определено 
	// для структуры из трёх целых чисел.
	// удаление ключа k из дерева p
	if (p == 0) return 0;
	// Двоичный поиск нужного элемента.
	if ((k.p.x<p->key.p.x) || ((k.p.x == p->key.p.x) && (k.p.y<p->key.p.y)) || ((k.p.x == p->key.p.x) && (k.p.y == p->key.p.y) && (k.p.z <  p->key.p.z)))
		p->left = remove_AVL(p->left, k);
	else if ((k.p.x>p->key.p.x) || ((k.p.x == p->key.p.x) && (k.p.y>p->key.p.y)) || ((k.p.x == p->key.p.x) && (k.p.y == p->key.p.y) && (k.p.z >  p->key.p.z)))
		p->right = remove_AVL(p->right, k);
	else // k==p->key
	{
		node_AVL1* q = p->left;
		node_AVL1* r = p->right;
		p->left = 0;
		p->right = 0;
		delete p;
		p = 0;
		if (r == 0) return q;
		node_AVL1* min = findmin(r);
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
node_AVL1* insert_and_modify(node_AVL1* p, SortElm k, SortElm k_search)
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
void print_AVL(node_AVL1* p)
{
	if (p != 0) {
		print_AVL(p->left);
		for (integer i = 0; i <= p->height; i++) {
			printf(" ");
		}
#if doubleintprecision == 1		
		printf("%e %e %e %lld\n", p->key.p.x, p->key.p.y, p->key.p.z, p->key.i);
#else
		printf("%e %e %e %d\n", p->key.p.x, p->key.p.y, p->key.p.z, p->key.i);
#endif

		print_AVL(p->right);
	}
} // print_AVL 


// АВЛ дерево конец.

  // Решение уравнения теплопередачи в 3D методом конечных элементов.
  // 6 августа 2017. - 19.07.2019.
void solve_Thermal(TEMPER &t, FLOW* &fglobal, TPROP* matlist, 
	WALL* &w, integer lw, integer lu, BLOCK* b, integer lb, integer ls,
	QuickMemVorst& m,
	bool bThermalStress, doublereal operatingtemperature,
	// для нестационарного температурного моделирования 10.11.2018
	bool btimedep, doublereal timestep_seq_current_step,
	doublereal* &toldtimestep, doublereal* &tnewtimestep, integer &maxelm_global_ret,
	doublereal poweron_multiplier_sequence, integer irealesation_selector,
	doublereal* &t_for_Mechanical)
{
     // tnewtimestep - результат вычисления нестационарного моделирования.

	

	integer maxelm_global = t.maxnod;
	integer ncell_shadow_gl = t.maxelm;
	for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
		maxelm_global += my_union[iu_74].t.maxnod;
		ncell_shadow_gl += my_union[iu_74].t.maxelm;
	}

	// Вычисляем минимальный размер ячейки для поиска совпадения по вещественным координатам.
	doublereal epsx = 1.0e+30, epsy = 1.0e+30, epsz = 1.0e+30;
	for (integer ie = 0; ie < t.maxelm; ie++) {
		doublereal hx = 0.0, hy = 0.0, hz = 0.0;
		volume3D(ie, t.nvtx, t.pa, hx, hy, hz);
		if (0.3*hx < epsx) epsx = 0.3*hx;
		if (0.3*hy < epsy) epsy = 0.3*hy;
		if (0.3*hz < epsz) epsz = 0.3*hz;
	}

	for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
		for (integer ie = 0; ie < my_union[iu_74].t.maxelm; ie++) {
			doublereal hx = 0.0, hy = 0.0, hz = 0.0;
			volume3D(ie, my_union[iu_74].t.nvtx, my_union[iu_74].t.pa, hx, hy, hz);
			if (0.3*hx < epsx) epsx = 0.3*hx;
			if (0.3*hy < epsy) epsy = 0.3*hy;
			if (0.3*hz < epsz) epsz = 0.3*hz;
		}
	}

	
	// TODO Черногория 25.08.2018 На основе АВЛ дерева не реализовано.
	//integer irealesation_selector = bAVLrealesation;// bARRAYrealesation;

	node_AVL1* root=0;
	SortElm* pa_global1 = nullptr;
	if (irealesation_selector == bARRAYrealesation) {
		pa_global1 = new SortElm[maxelm_global];
	}
	TOCHKA* pa_global = new TOCHKA[maxelm_global];
	maxelm_global = t.maxnod;
	for (integer i = 0; i < maxelm_global; i++) {
		pa_global[i] = t.pa[i];
		if (irealesation_selector == bARRAYrealesation) {
			pa_global1[i].p = t.pa[i];
			pa_global1[i].i = i;
		}
		if (irealesation_selector == bAVLrealesation) {
			SortElm inskey;
			inskey.i = i;
			inskey.p = t.pa[i];
			root=insert(root, inskey);
		}
	}
	integer **hash_table_pa = new integer*[lu];
	for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
		hash_table_pa[iu_74] = new integer[my_union[iu_74].t.maxnod];
	}

	if (irealesation_selector == bARRAYrealesation) {
		ShellSort(pa_global1, maxelm_global - 1);
	}

	for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
		for (integer j = 0; j < my_union[iu_74].t.maxnod; j++) {
			bool bfound = false;
			/*
			for (integer i = 0; i < maxelm_global; i++) {
				if ((fabs(pa_global[i].x - my_union[iu_74].t.pa[j].x)<epsx)&& (fabs(pa_global[i].y - my_union[iu_74].t.pa[j].y)<epsy) && (fabs(pa_global[i].z - my_union[iu_74].t.pa[j].z)<epsz)) {
					bfound = true;
					hash_table_pa[iu_74][j] = i;
					break;
				}
			}
			*/
			integer i = -1;
			if (irealesation_selector == bARRAYrealesation) {
				bfound = BinarySearch(pa_global1, my_union[iu_74].t.pa[j], i, maxelm_global - 1);
			}
			if (irealesation_selector == bAVLrealesation) {
				SortElm searchkey;
				searchkey.i = i;// не участвует в поиске.
				searchkey.p = my_union[iu_74].t.pa[j];
				bfound = isfound(root, searchkey);
				if (bfound) {
					hash_table_pa[iu_74][j] = searchkey.i;
				}
			}
			if (bfound) {
				if (irealesation_selector == bARRAYrealesation) {
					hash_table_pa[iu_74][j] = pa_global1[i].i;
				}
			}
			if (!bfound) {
				// вставка сохраняющая порядок сортировки.
				// поиск места для вставки
				if (irealesation_selector == bARRAYrealesation) {
					integer mid;
					for (integer i75 = i - 1; i75 < maxelm_global; i75++) {
						mid = i75;
						if (((my_union[iu_74].t.pa[j].x < pa_global1[mid].p.x) || ((my_union[iu_74].t.pa[j].x == pa_global1[mid].p.x) && (my_union[iu_74].t.pa[j].y < pa_global1[mid].p.y)) || ((my_union[iu_74].t.pa[j].x == pa_global1[mid].p.x) && (my_union[iu_74].t.pa[j].y == pa_global1[mid].p.y) && (my_union[iu_74].t.pa[j].z < pa_global1[mid].p.z)))) {
							break;
						}
					}
					// вставка. 
					// сдвиг вправо (неэффективно по памяти.).
					for (integer i95 = maxelm_global; i95 > mid; i95--) {
						pa_global1[i95] = pa_global1[i95 - 1];
					}
					// вставка в освобожденную позицию. 
					pa_global1[mid].p = my_union[iu_74].t.pa[j];
					pa_global1[mid].i = maxelm_global;
				}
				if (irealesation_selector == bAVLrealesation) {
					SortElm inskey;
					inskey.i = maxelm_global;
					inskey.p = my_union[iu_74].t.pa[j];
					root = insert(root, inskey);
				}
				pa_global[maxelm_global] = my_union[iu_74].t.pa[j];

				hash_table_pa[iu_74][j] = maxelm_global;

				maxelm_global++;
			}
		}
	}

	maxelm_global_ret = maxelm_global;

	delete[] pa_global1;
	pa_global1 = nullptr;
	
	if (irealesation_selector == bAVLrealesation) {
		clear_AVL(root);
	}
	
	printf("pa ok\n");
	//getchar();

	// Нужен nvtx, pa, prop.
	// pa --- Ok.

	// Делаем nvtx
	bool* bcheck_visible = nullptr;
		bcheck_visible = new bool[ncell_shadow_gl];
		if ((bcheck_visible == nullptr)) {
			printf("problem allocate memory for bcheck_visible array in solve_Thermal function in module mysolverv0_03.c\n");
			system("PAUSE");
			exit(1);
		}
	    doublereal* Ux_arr = nullptr;
		Ux_arr = new doublereal[ncell_shadow_gl];
		doublereal* Uy_arr = nullptr;
		Uy_arr = new doublereal[ncell_shadow_gl];
		doublereal* Uz_arr = nullptr;
		Uz_arr = new doublereal[ncell_shadow_gl];
		doublereal* mut_arr = nullptr;
		mut_arr = new doublereal[ncell_shadow_gl];
		if ((Ux_arr == nullptr)) {
			printf("problem allocate memory for Ux_arr array in solve_Thermal function in module mysolverv0_03.c\n");
			system("PAUSE");
			exit(1);
		}
		if ((Uy_arr == nullptr)) {
			printf("problem allocate memory for Uy_arr array in solve_Thermal function in module mysolverv0_03.c\n");
			system("PAUSE");
			exit(1);
		}
		if ((Uz_arr == nullptr)) {
			printf("problem allocate memory for Uz_arr array in solve_Thermal function in module mysolverv0_03.c\n");
			system("PAUSE");
			exit(1);
		}
		if ((mut_arr == nullptr)) {
			printf("problem allocate memory for mut_arr array in solve_Thermal function in module mysolverv0_03.c\n");
			system("PAUSE");
			exit(1);
		}
	
	integer** nvtx_global = new integer*[8];
	for (integer i = 0; i < 8; i++) {
		nvtx_global[i] = new integer[ncell_shadow_gl];
	}
		
	// copy first
	for (integer i = 0; i < t.maxelm; i++) {
		for (integer j = 0; j < 8; j++) {
			nvtx_global[j][i] = t.nvtx[j][i];
		}
		integer ib = t.whot_is_block[i];
		bcheck_visible[i] = b[ib].bvisible;
		if ((t.ptr!=nullptr)&&(fglobal!=nullptr)&&(t.ptr[0][i] > -1) && (fglobal[t.ptr[1][i]].potent!=nullptr)) {
			Ux_arr[i] = fglobal[t.ptr[1][i]].potent[VELOCITY_X_COMPONENT][t.ptr[0][i]];
			Uy_arr[i] = fglobal[t.ptr[1][i]].potent[VELOCITY_Y_COMPONENT][t.ptr[0][i]];
			Uz_arr[i] = fglobal[t.ptr[1][i]].potent[VELOCITY_Z_COMPONENT][t.ptr[0][i]];
			mut_arr[i] = fglobal[t.ptr[1][i]].potent[MUT][t.ptr[0][i]];
			//printf("%e %e %e %e %lld %lld\n", Ux_arr[i], Uy_arr[i], Uz_arr[i], mut_arr[i], t.ptr[1][i], t.ptr[0][i]);
			//system("pause");
		}
		else {
			Ux_arr[i] = 0.0;
			Uy_arr[i] = 0.0;
			Uz_arr[i] = 0.0;
			mut_arr[i] = 0.0;
		}
	}
	integer ic_nvtx = t.maxelm;
	//assembles.
	// для каждого ассемблеса надо найти позиции значений nvtx в новом объединенном pa.
	for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
		for (integer j = 0; j < my_union[iu_74].t.maxelm; j++) {
			
			TOCHKA pa_loc;
			for (integer k = 0; k < 8; k++) {
				pa_loc = my_union[iu_74].t.pa[my_union[iu_74].t.nvtx[k][j]-1];
				bool bfound = false;
				// Только если lu==1
				if (lu == 1) {
					// Мы присоеденили второй кусок pa к первому но он может быть лишь короче из-за совпадений узлов двух сеток при сшитии,
					// поэтому искомая позиция близка к первоначальной плюс заход влево из-за совпадения точек.
					/*
					for (integer i= min(t.maxnod+ my_union[iu_74].t.nvtx[k][j] - 1, maxelm_global-1); i>=0; i--) {
					    if ((fabs(pa_global[i].x - pa_loc.x) < epsx) && (fabs(pa_global[i].y - pa_loc.y) < epsy) && (fabs(pa_global[i].z - pa_loc.z) < epsz)) {
						   bfound = true;
						   nvtx_global[k][ic_nvtx] = i + 1;
						   break; // Прерываем поиск как только находим то что надо.
					    }
				    }
				*/
					nvtx_global[k][ic_nvtx]=hash_table_pa[iu_74][my_union[iu_74].t.nvtx[k][j]-1]+1;
					integer ib = my_union[iu_74].t.whot_is_block[j];
					bcheck_visible[ic_nvtx] = b[ib].bvisible;
					if ((my_union[iu_74].t.ptr != nullptr) && (my_union[iu_74].f != nullptr) && (my_union[iu_74].t.ptr[0][j] > -1)&&(my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent!=nullptr)) {
						Ux_arr[ic_nvtx] = my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent[VELOCITY_X_COMPONENT][my_union[iu_74].t.ptr[0][j]];
						Uy_arr[ic_nvtx] = my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent[VELOCITY_Y_COMPONENT][my_union[iu_74].t.ptr[0][j]];
						Uz_arr[ic_nvtx] = my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent[VELOCITY_Z_COMPONENT][my_union[iu_74].t.ptr[0][j]];
						mut_arr[ic_nvtx] = my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent[MUT][my_union[iu_74].t.ptr[0][j]];
					}
					else {
						Ux_arr[ic_nvtx] = 0.0;
						Uy_arr[ic_nvtx] = 0.0;
						Uz_arr[ic_nvtx] = 0.0;
						mut_arr[ic_nvtx] = 0.0;
					}
				}
				else {
					for (integer i = 0; i < maxelm_global; i++) {
						if ((fabs(pa_global[i].x - pa_loc.x) < epsx) && (fabs(pa_global[i].y - pa_loc.y) < epsy) && (fabs(pa_global[i].z - pa_loc.z) < epsz)) {
							bfound = true;
							nvtx_global[k][ic_nvtx] = i + 1;
							integer ib = my_union[iu_74].t.whot_is_block[j];
							bcheck_visible[ic_nvtx] = b[ib].bvisible;
							if ((my_union[iu_74].t.ptr != nullptr) && (my_union[iu_74].f != nullptr) && (my_union[iu_74].t.ptr[0][j] > -1) && (my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent != nullptr)) {
								Ux_arr[ic_nvtx] = my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent[VELOCITY_X_COMPONENT][my_union[iu_74].t.ptr[0][j]];
								Uy_arr[ic_nvtx] = my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent[VELOCITY_Y_COMPONENT][my_union[iu_74].t.ptr[0][j]];
								Uz_arr[ic_nvtx] = my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent[VELOCITY_Z_COMPONENT][my_union[iu_74].t.ptr[0][j]];
								mut_arr[ic_nvtx] = my_union[iu_74].f[my_union[iu_74].t.ptr[1][j]].potent[MUT][my_union[iu_74].t.ptr[0][j]];
							}
							else {
								Ux_arr[ic_nvtx] = 0.0;
								Uy_arr[ic_nvtx] = 0.0;
								Uz_arr[ic_nvtx] = 0.0;
								mut_arr[ic_nvtx] = 0.0;
							}
							break;
						}
					}
				}
				

			}
			ic_nvtx++;
		}
	}

	integer* icount_number_visit = new integer[maxelm_global];
	for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
		icount_number_visit[i_1] = 0;
	}
	for (integer i_1 = 0; i_1 < ncell_shadow_gl; i_1++) {
		for (integer i = 0; i < 8; i++) {
			icount_number_visit[nvtx_global[i][i_1]-1]++;
		}
	}
	// Если icount_number_visit<=4 то это граничный узел.

	// Отмечаем в true узлы которые принадлежат пользовательским стенкам.
	// Здесь просто инициализация значением false.
	bool* bwall_active=new bool[maxelm_global]; 
	for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
		bwall_active[i_1] = false;
	}

	printf("nvtx ok\n");
	//getchar();
	// pa - Ok.
	// nvtx - Ok.
	// maxnod == maxelm_global
	doublereal* lam_export = new doublereal[maxelm_global];
	// для нестационарных моделирований.
	doublereal* rho_export = new doublereal[maxelm_global];
	doublereal* Cp_export = new doublereal[maxelm_global];
	doublereal* Vol_export= new doublereal[maxelm_global];
	
	for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
		Vol_export[i_1] = 0.0;		
	}

	

	// Делаем prop_global
	doublereal** prop_global = new doublereal*[8];
	for (integer i = 0; i < 8; i++) {
		prop_global[i] = new doublereal[ncell_shadow_gl];
	}
	// copy cabinet
	
	for (integer ie = 0; ie < t.maxelm; ie++) {
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
		volume3D(ie, nvtx_global, pa_global, dx, dy, dz);

		for (integer j = 0; j < 8; j++) {
			prop_global[j][ie] = t.prop[j][ie];

			ic_nvtx = nvtx_global[j][ie]-1; // maxelm->maxnod

			lam_export[ic_nvtx] = t.prop[LAM][ie];
			rho_export[ic_nvtx] = t.prop[RHO][ie];
			Cp_export[ic_nvtx] = t.prop[HEAT_CAPACITY][ie];

			Vol_export[ic_nvtx] += 0.125 * dx * dy * dz;
			
		}
		
		
	}
	ic_nvtx = t.maxelm;
	integer ic_nvtx1;
	// copy assembles
	for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
		for (integer j = 0; j < my_union[iu_74].t.maxelm; j++) {
			// вычисление размеров текущего контрольного объёма:
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
			volume3D(j, my_union[iu_74].t.nvtx, my_union[iu_74].t.pa, dx, dy, dz);

			for (integer k = 0; k < 8; k++) {
				prop_global[k][ic_nvtx] = my_union[iu_74].t.prop[k][j];

				integer ie = ic_nvtx;
				ic_nvtx1 = nvtx_global[k][ie] - 1; // maxelm->maxnod


				lam_export[ic_nvtx1] = my_union[iu_74].t.prop[LAM][j];
				rho_export[ic_nvtx1] = my_union[iu_74].t.prop[RHO][j];
				Cp_export[ic_nvtx1] = my_union[iu_74].t.prop[HEAT_CAPACITY][j];
				Vol_export[ic_nvtx1] += 0.125*dx*dy*dz;
				
			}
			ic_nvtx++;
		}
	}

	printf("prop ok\n");
	//getchar();
	// prop_global готов.

	printf("New temperature solver with all meshes n=%lld\n", maxelm_global);

	doublereal* rthdsd = new doublereal[maxelm_global + 2]; // Правая часть.
	doublereal* temp_potent = new doublereal[maxelm_global + 2]; // Температура.
	bool* constr = new bool[maxelm_global + 2]; // Фиксированная температура.

	doublereal told_iter = operatingtemperature;

	// инициализация.
	for (integer i_1 = 0; i_1 < maxelm_global + 2; i_1++) {
		rthdsd[i_1] = 0.0;
		temp_potent[i_1] = operatingtemperature;
		constr[i_1] = false; // По умолчанию все узлы свободны.
	}
	if (btimedep) {
		// Нестационарный температурный солвер.
		for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
			temp_potent[i_1] = toldtimestep[i_1];
		}
	}
	//for (integer i_1 = 0; i_1 < t.maxelm + t.maxnod; i_1++) {
		//t.potent[i_1]= operatingtemperature;
	//}
	//for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
		//for (integer i_1 = 0; i_1 < my_union[iu_74].t.maxelm + my_union[iu_74].t.maxnod; i_1++) {
			//my_union[iu_74].t.potent[i_1] = operatingtemperature;
		//}
	//}

	bool* boundary = new bool[maxelm_global];
	for (integer i32 = 0; i32 < maxelm_global; i32++) {
		boundary[i32] = false;
	}

	
	

	//integer ie = 0;
	//for (integer j = 0; j < 8; j++) {
	//	printf("%e %e %e\n", t.pa[t.nvtx[j][ie] - 1].x, t.pa[t.nvtx[j][ie] - 1].y, t.pa[t.nvtx[j][ie] - 1].z);
	//}
	//	getchar();
	// Учёт граничных условий.
	// На зафиксированных участках границы мы выставляем флаг true.
	for (integer i_1 = 0; i_1 < lw; i_1++) {

		//const doublereal eps1 = 1.0e-30;
		// pa нумеруется  нуля.
		// maxelm_global == maxnod
		for (integer j_1 = 0; j_1 < maxelm_global; j_1++) {
			bool bfound = false;
			switch (w[i_1].iPlane) {
			case XY_PLANE: 
				if ((fabs(pa_global[j_1].z - w[i_1].g.zS)<epsz) && 
					(pa_global[j_1].x<w[i_1].g.xE + epsx) &&
					(pa_global[j_1].x>w[i_1].g.xS - epsx) &&
					(pa_global[j_1].y>w[i_1].g.yS - epsy) &&
					(pa_global[j_1].y<w[i_1].g.yE + epsy)) {
				    bfound = true;
			    }
				break;
			case YZ_PLANE:
				if ((fabs(pa_global[j_1].x - w[i_1].g.xS)<epsx) &&
					(pa_global[j_1].z<w[i_1].g.zE + epsz) &&
					(pa_global[j_1].z>w[i_1].g.zS - epsz) && 
					(pa_global[j_1].y>w[i_1].g.yS - epsy) && 
					(pa_global[j_1].y<w[i_1].g.yE + epsy)) {
					bfound = true;
				}
				break;
			case XZ_PLANE:
				if ((fabs(pa_global[j_1].y - w[i_1].g.yS)<epsy) && 
					(pa_global[j_1].z<w[i_1].g.zE + epsz) &&
					(pa_global[j_1].z>w[i_1].g.zS - epsz) &&
					(pa_global[j_1].x>w[i_1].g.xS - epsx) && 
					(pa_global[j_1].x<w[i_1].g.xE + epsx)) {
					bfound = true;
				}
				break;
			}
			
			if (bfound) {

				boundary[j_1] = true;

				bwall_active[j_1] = true;// Пользовательская стенка.

				// Фиксированный потенциал.
				if (w[i_1].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) {
					constr[j_1] = true;
					temp_potent[j_1] = w[i_1].Tamb;
					rthdsd[j_1]= w[i_1].Tamb;
					//printf("%d Tamb=%e\n",i_1, w[i_1].Tamb);
					//getchar();
				}
				// граничное условие Стефана - Больцмана
				if (w[i_1].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) {
					// 16.12.2018

					breakRUMBAcalc_for_nonlinear_boundary_condition = true;
					doublereal alpha_relax1 = 0.25;

					if ((temp_potent[j_1] < -271.0) && (temp_potent[j_1] <= w[i_1].Tamb)) {
						// Защита от опускания температуры ниже абсолютного нуля.
						// 7.12.2019
						// qb=0.0;
					}
					else {

						if (temp_potent[j_1] < -272.15) {
							temp_potent[j_1] = -272.15;
						}

						// Температура с предыдущего шага по времени !!!. Работает только в нестационарном солвере.

						//rthdsd[j_1] = alpha_relax1 * (-w[i_1].emissivity*STEFAN_BOLCMAN_CONST*((273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1]) - (273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb))) +
							//(1.0 - alpha_relax1)*(-w[i_1].emissivity*STEFAN_BOLCMAN_CONST*((273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1]) - (273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)));

						rthdsd[j_1] = (-w[i_1].emissivity*w[i_1].ViewFactor*STEFAN_BOLCMAN_CONST*(
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1]) -
							(273.15 + w[i_1].Tamb)*
							(273.15 + w[i_1].Tamb)*
							(273.15 + w[i_1].Tamb)*
							(273.15 + w[i_1].Tamb)));


						doublereal dS = 0.0;
						integer ie1, j81 = 0;

						for (ie1 = 0; ie1 < t.maxelm; ie1++) {
							for (j81 = 0; j81 < 8; j81++) {
								if (t.nvtx[j81][ie1] - 1 == j_1) {
									doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
									volume3D(ie1, nvtx_global, pa_global, dx, dy, dz);
									switch (w[i_1].iPlane) {
									case XY_PLANE:
										dS += 0.25 * dx * dy;
										break;
									case XZ_PLANE:
										dS += 0.25 * dx * dz;
										break;
									case YZ_PLANE: dS += 0.25 * dy * dz;
										break;
									}

								}
							}
						}

						ic_nvtx = t.maxelm;
						
						// copy assembles
						for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
							for (integer j = 0; j < my_union[iu_74].t.maxelm; j++) {								
								
								for (integer k = 0; k < 8; k++) {
									prop_global[k][ic_nvtx] = my_union[iu_74].t.prop[k][j];

									integer ie = ic_nvtx;
									ic_nvtx1 = nvtx_global[k][ie] - 1; // maxelm->maxnod

									if (ic_nvtx1 == j_1) {
										// вычисление размеров текущего контрольного объёма:
										doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
										volume3D(j, my_union[iu_74].t.nvtx, my_union[iu_74].t.pa, dx, dy, dz);

										switch (w[i_1].iPlane) {
										case XY_PLANE:
											dS += 0.25 * dx * dy;
											break;
										case XZ_PLANE:
											dS += 0.25 * dx * dz;
											break;
										case YZ_PLANE: dS += 0.25 * dy * dz;
											break;
										}

									}

								}
								ic_nvtx++;
							}
						}


						rthdsd[j_1] *= dS;

						
					
					}
				}// Stefan-Bolcman

				// Граничное условие Ньютона-Рихмана.
				if (w[i_1].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) {
					// 16.12.2018

					breakRUMBAcalc_for_nonlinear_boundary_condition = true;

					rthdsd[j_1] = -w[i_1].film_coefficient*(temp_potent[j_1] - w[i_1].Tamb);


					doublereal dS = 0.0;
					integer ie1, j81 = 0;

					for (ie1 = 0; ie1 < t.maxelm; ie1++) {
						for (j81 = 0; j81 < 8; j81++) {
							if (t.nvtx[j81][ie1] - 1 == j_1) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
								volume3D(ie1, nvtx_global, pa_global, dx, dy, dz);
								switch (w[i_1].iPlane) {
								case XY_PLANE: 
									dS += 0.25 * dx * dy;									
									break;
								case XZ_PLANE:
									dS += 0.25 * dx * dz;									
									break;
								case YZ_PLANE: dS += 0.25 * dy * dz;									
									break;
								}
								
							}
						}
					}

					ic_nvtx = t.maxelm;

					// copy assembles
					for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
						for (integer j = 0; j < my_union[iu_74].t.maxelm; j++) {

							for (integer k = 0; k < 8; k++) {
								prop_global[k][ic_nvtx] = my_union[iu_74].t.prop[k][j];

								integer ie = ic_nvtx;
								ic_nvtx1 = nvtx_global[k][ie] - 1; // maxelm->maxnod

								if (ic_nvtx1 == j_1) {
									// вычисление размеров текущего контрольного объёма:
									doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
									volume3D(j, my_union[iu_74].t.nvtx, my_union[iu_74].t.pa, dx, dy, dz);

									switch (w[i_1].iPlane) {
									case XY_PLANE:
										dS += 0.25 * dx * dy;
										break;
									case XZ_PLANE:
										dS += 0.25 * dx * dz;
										break;
									case YZ_PLANE: dS += 0.25 * dy * dz;
										break;
									}

								}

							}
							ic_nvtx++;
						}
					}


					rthdsd[j_1] *= dS;
					

				} // Newton-Richman

			}
		}

	}
	
	doublereal* temp_potent_old = new doublereal[maxelm_global];



	for (integer j_1 = 0; j_1 < maxelm_global; j_1++) {
		temp_potent_old[j_1] = temp_potent[j_1];
	}

	doublereal* dsquare_default_non_linear = new doublereal[maxelm_global];
	for (integer i32 = 0; i32 < maxelm_global; i32++) {
		dsquare_default_non_linear[i32] = 0.0;
	}
	doublereal* emissivity_default_non_linear= new doublereal[maxelm_global];
	for (integer i32 = 0; i32 < maxelm_global; i32++) {
		emissivity_default_non_linear[i32] = 0.0;
	}

	if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
		(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
		(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {


		for (integer ie1 = 0; ie1 < t.maxelm; ie1++) {
			for (integer j81 = 0; j81 < 8; j81++) {
				integer j_1 = t.nvtx[j81][ie1] - 1;
				if (((icount_number_visit[j_1] >= 2)) && (icount_number_visit[j_1] <= 4) && (bwall_active[j_1] == false)) {

					integer iPlane = -1;

					if (t.neighbors_for_the_internal_node[E_SIDE][ie1].iNODE1 >= t.maxelm) {
						iPlane = YZ_PLANE;
						integer inumber = t.neighbors_for_the_internal_node[E_SIDE][ie1].iNODE1 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
					}
					else if (t.neighbors_for_the_internal_node[W_SIDE][ie1].iNODE1 >= t.maxelm) {
						iPlane = YZ_PLANE;
						integer inumber = t.neighbors_for_the_internal_node[W_SIDE][ie1].iNODE1 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
					}
					else if (t.neighbors_for_the_internal_node[N_SIDE][ie1].iNODE1 >= t.maxelm) {
						iPlane = XZ_PLANE;
						integer inumber = t.neighbors_for_the_internal_node[N_SIDE][ie1].iNODE1 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
					}
					else if (t.neighbors_for_the_internal_node[S_SIDE][ie1].iNODE1 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[S_SIDE][ie1].iNODE1 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[T_SIDE][ie1].iNODE1 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[T_SIDE][ie1].iNODE1 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XY_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[B_SIDE][ie1].iNODE1 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[B_SIDE][ie1].iNODE1 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XY_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[E_SIDE][ie1].iNODE2 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[E_SIDE][ie1].iNODE2 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = YZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[W_SIDE][ie1].iNODE2 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[W_SIDE][ie1].iNODE2 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = YZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[N_SIDE][ie1].iNODE2 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[N_SIDE][ie1].iNODE2 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[S_SIDE][ie1].iNODE2 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[S_SIDE][ie1].iNODE2 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[T_SIDE][ie1].iNODE2 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[T_SIDE][ie1].iNODE2 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XY_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[B_SIDE][ie1].iNODE2 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[B_SIDE][ie1].iNODE2 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XY_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[E_SIDE][ie1].iNODE3 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[E_SIDE][ie1].iNODE3 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = YZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[W_SIDE][ie1].iNODE3 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[W_SIDE][ie1].iNODE3 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = YZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[N_SIDE][ie1].iNODE3 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[N_SIDE][ie1].iNODE3 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[S_SIDE][ie1].iNODE3 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[S_SIDE][ie1].iNODE3 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[T_SIDE][ie1].iNODE3 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[T_SIDE][ie1].iNODE3 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XY_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[B_SIDE][ie1].iNODE3 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[B_SIDE][ie1].iNODE3 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XY_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[E_SIDE][ie1].iNODE4 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[E_SIDE][ie1].iNODE4 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = YZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[W_SIDE][ie1].iNODE4 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[W_SIDE][ie1].iNODE4 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = YZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[N_SIDE][ie1].iNODE4 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[N_SIDE][ie1].iNODE4 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[S_SIDE][ie1].iNODE4 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[S_SIDE][ie1].iNODE4 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XZ_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[T_SIDE][ie1].iNODE4 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[T_SIDE][ie1].iNODE4 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XY_PLANE;
					}
					else if (t.neighbors_for_the_internal_node[B_SIDE][ie1].iNODE4 >= t.maxelm) {
						integer inumber = t.neighbors_for_the_internal_node[B_SIDE][ie1].iNODE4 - t.maxelm;
						emissivity_default_non_linear[j_1] = t.border_neighbor[inumber].emissivity;
						iPlane = XY_PLANE;
					}
					else {
						//printf("Error!!! plane unknown!");
						//getchar();
					}

					// Нормаль определена. Вычисление площади.

					if (iPlane > -1) {
						doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
						volume3D(ie1, nvtx_global, pa_global, dx, dy, dz);
						switch (iPlane) {
						case XY_PLANE:
							dsquare_default_non_linear[j_1] += 0.25 * dx * dy;
							break;
						case XZ_PLANE:
							dsquare_default_non_linear[j_1] += 0.25 * dx * dz;
							break;
						case YZ_PLANE:
							dsquare_default_non_linear[j_1] += 0.25 * dy * dz;
							break;
						}
					}
				
				}
			}
		}

		{
			integer ie1 = t.maxelm;

			// Асемблесы не учитываются 21.07.2020.

			// copy assembles
			for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
				for (integer j = 0; j < my_union[iu_74].t.maxelm; j++) {

					for (integer j81 = 0; j81 < 8; j81++) {
						
						integer j_1 = nvtx_global[j81][ie1] - 1; // maxelm->maxnod

						
						if (((icount_number_visit[j_1] >= 2)) && (icount_number_visit[j_1] <= 4) && (bwall_active[j_1] == false)) {

							integer iPlane = -1;

							if (my_union[iu_74].t.neighbors_for_the_internal_node[E_SIDE][j].iNODE1 >= my_union[iu_74].t.maxelm) {
								iPlane = YZ_PLANE;
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[E_SIDE][j].iNODE1 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[W_SIDE][j].iNODE1 >= my_union[iu_74].t.maxelm) {
								iPlane = YZ_PLANE;
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[W_SIDE][j].iNODE1 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[N_SIDE][j].iNODE1 >= my_union[iu_74].t.maxelm) {
								iPlane = XZ_PLANE;
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[N_SIDE][j].iNODE1 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[S_SIDE][j].iNODE1 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[S_SIDE][j].iNODE1 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[T_SIDE][j].iNODE1 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[T_SIDE][j].iNODE1 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XY_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[B_SIDE][j].iNODE1 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[B_SIDE][j].iNODE1 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XY_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[E_SIDE][j].iNODE2 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[E_SIDE][j].iNODE2 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = YZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[W_SIDE][j].iNODE2 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[W_SIDE][j].iNODE2 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = YZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[N_SIDE][j].iNODE2 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[N_SIDE][j].iNODE2 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[S_SIDE][j].iNODE2 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[S_SIDE][j].iNODE2 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[T_SIDE][j].iNODE2 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[T_SIDE][j].iNODE2 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XY_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[B_SIDE][j].iNODE2 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[B_SIDE][j].iNODE2 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XY_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[E_SIDE][j].iNODE3 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[E_SIDE][j].iNODE3 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = YZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[W_SIDE][j].iNODE3 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[W_SIDE][j].iNODE3 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = YZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[N_SIDE][j].iNODE3 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[N_SIDE][j].iNODE3 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[S_SIDE][j].iNODE3 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[S_SIDE][j].iNODE3 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[T_SIDE][j].iNODE3 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[T_SIDE][j].iNODE3 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XY_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[B_SIDE][j].iNODE3 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[B_SIDE][j].iNODE3 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XY_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[E_SIDE][j].iNODE4 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[E_SIDE][j].iNODE4 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = YZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[W_SIDE][j].iNODE4 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[W_SIDE][j].iNODE4 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = YZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[N_SIDE][j].iNODE4 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[N_SIDE][j].iNODE4 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[S_SIDE][j].iNODE4 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[S_SIDE][j].iNODE4 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XZ_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[T_SIDE][j].iNODE4 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[T_SIDE][j].iNODE4 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XY_PLANE;
							}
							else if (my_union[iu_74].t.neighbors_for_the_internal_node[B_SIDE][j].iNODE4 >= my_union[iu_74].t.maxelm) {
								integer inumber = my_union[iu_74].t.neighbors_for_the_internal_node[B_SIDE][j].iNODE4 - my_union[iu_74].t.maxelm;
								emissivity_default_non_linear[j_1] = my_union[iu_74].t.border_neighbor[inumber].emissivity;
								iPlane = XY_PLANE;
							}
							else {
								//printf("Error!!! plane unknown!");
								//getchar();
							}

							// Нормаль определена. Вычисление площади.

							if (iPlane > -1) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
								volume3D(ie1, nvtx_global, pa_global, dx, dy, dz);
								switch (iPlane) {
								case XY_PLANE:
									dsquare_default_non_linear[j_1] += 0.25 * dx * dy;
									break;
								case XZ_PLANE:
									dsquare_default_non_linear[j_1] += 0.25 * dx * dz;
									break;
								case YZ_PLANE:
									dsquare_default_non_linear[j_1] += 0.25 * dy * dz;
									break;
								}
							}

						}						

					}
					ie1++;
				}
			}
		}

				
		for (integer j_1 = 0; j_1 < maxelm_global; j_1++) {
			if (((icount_number_visit[j_1] >= 2)) && (icount_number_visit[j_1] <= 4) && (bwall_active[j_1] == false)) {
					
				if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {
						

					breakRUMBAcalc_for_nonlinear_boundary_condition = true;

					boundary[j_1] = true;

					// film_coefficient - глобальный.
					// operating_temperature_for_film_coeff - глобальный.
					rthdsd[j_1] = -film_coefficient*(temp_potent[j_1] - operating_temperature_for_film_coeff);
					//printf("film_coefficient=%e operating_temperature_for_film_coeff=%e\n", film_coefficient, operating_temperature_for_film_coeff);
					//getchar();
					
					// Определение Нормали к стенке.
					// Работает на АЛИС сетке тоже.
					// Без асемблесов.						

					// Нормаль определена. Вычисление площади.
					rthdsd[j_1] *= dsquare_default_non_linear[j_1];

				} // Newton-Richman

				if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) {

					breakRUMBAcalc_for_nonlinear_boundary_condition = true;
					doublereal alpha_relax1 = 0.25;

					boundary[j_1] = true;

					if ((temp_potent[j_1] < -271.0) && (temp_potent[j_1] <= operating_temperature_for_film_coeff)) {
						// Защита от опускания температуры ниже абсолютного нуля.
						// 7.12.2019
						// qb=0.0;
					}
					else {

						if (temp_potent[j_1] < -272.15) {
							temp_potent[j_1] = -272.15;
						}

						// Температура с предыдущего шага по времени !!!. Работает только в нестационарном солвере.

						//rthdsd[j_1] = alpha_relax1 * (-emissivity_default_non_linear[j_1]*STEFAN_BOLCMAN_CONST*((273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1]) - (273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff))) +
						//(1.0 - alpha_relax1)*(-emissivity_default_non_linear[j_1]*STEFAN_BOLCMAN_CONST*((273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1]) - (273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)));
					
						

						rthdsd[j_1] = (-emissivity_default_non_linear[j_1] *STEFAN_BOLCMAN_CONST*(
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1]) -
							(273.15 + operating_temperature_for_film_coeff)*
							(273.15 + operating_temperature_for_film_coeff)*
							(273.15 + operating_temperature_for_film_coeff)*
							(273.15 + operating_temperature_for_film_coeff)));
							

						/*
						rthdsd[j_1] = alpha_relax1*(-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1]) -
							(273.15 + operating_temperature_for_film_coeff)*
							(273.15 + operating_temperature_for_film_coeff)*
							(273.15 + operating_temperature_for_film_coeff)*
							(273.15 + operating_temperature_for_film_coeff))) +
							(1.0 - alpha_relax1)*(-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
							(273.15 + temp_potent_old[j_1])*
								(273.15 + temp_potent_old[j_1])*
								(273.15 + temp_potent_old[j_1])*
								(273.15 + temp_potent_old[j_1]) -
								(273.15 + operating_temperature_for_film_coeff)*
								(273.15 + operating_temperature_for_film_coeff)*
								(273.15 + operating_temperature_for_film_coeff)*
								(273.15 + operating_temperature_for_film_coeff)));
								*/

						// Определение Нормали к стенке. 
						// Работает на АЛИС сетке тоже.
						// Без асемблесов.


						// Площадь была вычислена один раз и запомнена.
						rthdsd[j_1] *= dsquare_default_non_linear[j_1];

					}
				} // Stefan-Bolcman

				if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC) {

					boundary[j_1] = true;

					breakRUMBAcalc_for_nonlinear_boundary_condition = true;
					doublereal alpha_relax1 = 0.25;

					//if ((temp_potent[j_1] < -271.0) && (temp_potent[j_1] <= operating_temperature_for_film_coeff)) {
					// Защита от опускания температуры ниже абсолютного нуля.
					// 7.12.2019
					// qb=0.0;
					//}
					//else 
					{

						if (temp_potent[j_1] < -272.15) {
							temp_potent[j_1] = -272.15;
						}
						if (temp_potent_old[j_1] < -272.15) {
							temp_potent_old[j_1] = -272.15;
						}

						// Температура с предыдущего шага по времени !!!. Работает только в нестационарном солвере.

						//rthdsd[j_1] = alpha_relax1 * (-emissivity_default_non_linear[j_1]*STEFAN_BOLCMAN_CONST*((273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1]) - (273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff))) +
						//(1.0 - alpha_relax1)*(-emissivity_default_non_linear[j_1]*STEFAN_BOLCMAN_CONST*((273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1]) - (273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)));

						rthdsd[j_1] = -film_coefficient*(temp_potent[j_1] - operating_temperature_for_film_coeff);
						rthdsd[j_1] += (-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1])*
							(273.15 + temp_potent[j_1]) -
							(273.15 + operating_temperature_for_film_coeff)*
							(273.15 + operating_temperature_for_film_coeff)*
							(273.15 + operating_temperature_for_film_coeff)*
							(273.15 + operating_temperature_for_film_coeff)));

						/*
						rthdsd[j_1] = alpha_relax1*(-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
						(273.15 + temp_potent[j_1])*
						(273.15 + temp_potent[j_1])*
						(273.15 + temp_potent[j_1])*
						(273.15 + temp_potent[j_1]) -
						(273.15 + operating_temperature_for_film_coeff)*
						(273.15 + operating_temperature_for_film_coeff)*
						(273.15 + operating_temperature_for_film_coeff)*
						(273.15 + operating_temperature_for_film_coeff)))+
						(1.0-alpha_relax1)*(-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
						(273.15 + temp_potent_old[j_1])*
						(273.15 + temp_potent_old[j_1])*
						(273.15 + temp_potent_old[j_1])*
						(273.15 + temp_potent_old[j_1]) -
						(273.15 + operating_temperature_for_film_coeff)*
						(273.15 + operating_temperature_for_film_coeff)*
						(273.15 + operating_temperature_for_film_coeff)*
						(273.15 + operating_temperature_for_film_coeff)));
						*/

						// Определение Нормали к стенке. 
						// Работает на АЛИС сетке тоже.
						//Без асемблесов.


						// Площадь была вычислена один раз и запомнена.
						rthdsd[j_1] *= dsquare_default_non_linear[j_1];

					}
				} // Mix-Condition

			}
		} 
		
	}
	
	// Инициализация.
	for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
		bwall_active[i_1] = false;
	}


	doublereal* dsquare_default_non_linear1 = new doublereal[maxelm_global];
	for (integer i32 = 0; i32 < maxelm_global; i32++) {
		dsquare_default_non_linear1[i32] = 0.0;
	}

	for (integer i_1 = 0; i_1 < lw; i_1++) {

		//const doublereal eps1 = 1.0e-30;
		// pa нумеруется  нуля.
		for (integer j_1 = 0; j_1 < maxelm_global; j_1++) {
			bool bfound = false;
			switch (w[i_1].iPlane) {
			case XY_PLANE: if ((fabs(pa_global[j_1].z - w[i_1].g.zS) < epsz) &&
				(pa_global[j_1].x < w[i_1].g.xE + epsx) &&
				(pa_global[j_1].x > w[i_1].g.xS - epsx) &&
				(pa_global[j_1].y > w[i_1].g.yS - epsy) &&
				(pa_global[j_1].y < w[i_1].g.yE + epsy)) {
				bfound = true;
			}
						 break;
			case YZ_PLANE:
				if ((fabs(pa_global[j_1].x - w[i_1].g.xS) < epsx) &&
					(pa_global[j_1].z < w[i_1].g.zE + epsz) &&
					(pa_global[j_1].z > w[i_1].g.zS - epsz) &&
					(pa_global[j_1].y > w[i_1].g.yS - epsy) &&
					(pa_global[j_1].y < w[i_1].g.yE + epsy)) {
					bfound = true;
				}
				break;
			case XZ_PLANE:
				if ((fabs(pa_global[j_1].y - w[i_1].g.yS) < epsy) &&
					(pa_global[j_1].z < w[i_1].g.zE + epsz) &&
					(pa_global[j_1].z > w[i_1].g.zS - epsz) &&
					(pa_global[j_1].x > w[i_1].g.xS - epsx) &&
					(pa_global[j_1].x < w[i_1].g.xE + epsx)) {
					bfound = true;
				}
				break;
			}
			if (bfound) {

				
				

				// граничное условие Стефана - Больцмана
				if (w[i_1].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) {
					// 16.12.2018

					{

						

						doublereal dS = 0.0;
						integer ie1, j81 = 0;

						for (ie1 = 0; ie1 < t.maxelm; ie1++) {
							for (j81 = 0; j81 < 8; j81++) {
								if (t.nvtx[j81][ie1] - 1 == j_1) {
									doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
									volume3D(ie1, nvtx_global, pa_global, dx, dy, dz);
									switch (w[i_1].iPlane) {
									case XY_PLANE:
										dS += 0.25 * dx * dy;
										break;
									case XZ_PLANE:
										dS += 0.25 * dx * dz;
										break;
									case YZ_PLANE: dS += 0.25 * dy * dz;
										break;
									}

								}
							}
						}

						ic_nvtx = t.maxelm;

						// copy assembles
						for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
							for (integer j = 0; j < my_union[iu_74].t.maxelm; j++) {

								for (integer k = 0; k < 8; k++) {
									prop_global[k][ic_nvtx] = my_union[iu_74].t.prop[k][j];

									integer ie = ic_nvtx;
									ic_nvtx1 = nvtx_global[k][ie] - 1; // maxelm->maxnod

									if (ic_nvtx1 == j_1) {
										// вычисление размеров текущего контрольного объёма:
										doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
										volume3D(j, my_union[iu_74].t.nvtx, my_union[iu_74].t.pa, dx, dy, dz);

										switch (w[i_1].iPlane) {
										case XY_PLANE:
											dS += 0.25 * dx * dy;
											break;
										case XZ_PLANE:
											dS += 0.25 * dx * dz;
											break;
										case YZ_PLANE: dS += 0.25 * dy * dz;
											break;
										}

									}

								}
								ic_nvtx++;
							}
						}


						dsquare_default_non_linear1[j_1] = dS;

					}
				} // Stefan-Bolcman

				// Граничное условие Ньютона-Рихмана.
				if (w[i_1].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) {
					// 16.12.2018



					doublereal dS = 0.0;
					integer ie1, j81 = 0;

					for (ie1 = 0; ie1 < t.maxelm; ie1++) {
						for (j81 = 0; j81 < 8; j81++) {
							if (t.nvtx[j81][ie1] - 1 == j_1) {
								doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
								volume3D(ie1, nvtx_global, pa_global, dx, dy, dz);
								switch (w[i_1].iPlane) {
								case XY_PLANE:
									dS += 0.25 * dx * dy;
									break;
								case XZ_PLANE:
									dS += 0.25 * dx * dz;
									break;
								case YZ_PLANE: dS += 0.25 * dy * dz;
									break;
								}

							}
						}
					}

					ic_nvtx = t.maxelm;

					doublereal dS1 = dS;

					// copy assembles
					for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
						for (integer j = 0; j < my_union[iu_74].t.maxelm; j++) {

							for (integer k = 0; k < 8; k++) {
								prop_global[k][ic_nvtx] = my_union[iu_74].t.prop[k][j];

								integer ie = ic_nvtx;
								ic_nvtx1 = nvtx_global[k][ie] - 1; // maxelm->maxnod

								if (ic_nvtx1 == j_1) {
									// вычисление размеров текущего контрольного объёма:
									doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
									volume3D(j, my_union[iu_74].t.nvtx, my_union[iu_74].t.pa, dx, dy, dz);

									switch (w[i_1].iPlane) {
									case XY_PLANE:
										dS += 0.25 * dx * dy;
										break;
									case XZ_PLANE:
										dS += 0.25 * dx * dz;
										break;
									case YZ_PLANE: dS += 0.25 * dy * dz;
										break;
									}

								}

							}
							ic_nvtx++;
						}
					}

					//dS1 = dS - dS1;
					//if (dS1 > 1.0e-30) {
						//printf("%e \n", dS1); getchar();
				    //}

					dsquare_default_non_linear1[j_1] = dS;

				} // Newton-Richman
			}
		}
	}


	// для преобразования из центра в края.
	doublereal* lam_for_export = new doublereal[maxelm_global];
	doublereal* sum_vol = new doublereal[maxelm_global]; // суммарный объем.


	doublereal** Kmatrix_local = nullptr;
	Kmatrix_local = new doublereal * [8];
	for (integer i_1 = 0; i_1 < 8; i_1++) {
		Kmatrix_local[i_1] = new doublereal[8];
	}


	

	integer iprohod = 0;
	doublereal temp_max= -1.0e30;
	doublereal temp_min = +1.0e30;
	do {

		iprohod++;
		printf("nonlinear prohod number %lld\n", iprohod);

		// В случае мощности зависящей от температуры производит 
		// обновление рабочего значения мощности.
		//update_power_temperature_depend(s, ls, t, t.border_neighbor, gtdps, ltdp, toldtimestep1);
		//for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
			//update_power_temperature_depend(s, ls, my_union[iu_74].t, my_union[iu_74].t.border_neighbor, gtdps, ltdp, toldtimestep);
		//}
		
		update_temp_properties1(t, fglobal, b, lb, matlist, temp_potent,0,-1, lam_export, nvtx_global); // обновляем свойства твёрдых материалов
		for (integer iu_74 = 0; iu_74 < lu; iu_74++) {			
			integer iadd = t.maxelm;
			for (integer iu_75 = 0; iu_75 < iu_74; iu_75++) {
				iadd += my_union[iu_75].t.maxelm;
			}
			update_temp_properties1(my_union[iu_74].t, my_union[iu_74].f, b, lb, matlist, temp_potent,iadd,iu_74, lam_export, nvtx_global); // обновляем свойства твёрдых материалов
		}
		
		// copy cabinet
		for (integer i = 0; i < t.maxelm; i++) {
			for (integer j = 0; j < 8; j++) {
				prop_global[j][i] = t.prop[j][i];
			}
		}
		integer ic_nvtx1 = t.maxelm;
		// copy assembles
		for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
			for (integer j = 0; j < my_union[iu_74].t.maxelm; j++) {
				for (integer k = 0; k < 8; k++) {
					prop_global[k][ic_nvtx1] = my_union[iu_74].t.prop[k][j];
				}
				ic_nvtx1++;
			}
		}


		for (integer i_1 = 0; i_1 < maxelm_global + 2; i_1++) {
			rthdsd[i_1] = 0.0;
		}
		for (integer i_1 = 0; i_1 < lw; i_1++) {

			//const doublereal eps1 = 1.0e-30;
			// pa нумеруется  нуля.
			for (integer j_1 = 0; j_1 < maxelm_global; j_1++) {
				bool bfound = false;
				switch (w[i_1].iPlane) {
				case XY_PLANE: if ((fabs(pa_global[j_1].z - w[i_1].g.zS) < epsz) && 
					(pa_global[j_1].x < w[i_1].g.xE + epsx) && 
					(pa_global[j_1].x > w[i_1].g.xS - epsx) && 
					(pa_global[j_1].y > w[i_1].g.yS - epsy) && 
					(pa_global[j_1].y < w[i_1].g.yE + epsy)) {
					bfound = true;
				}
						 break;
				case YZ_PLANE:
					if ((fabs(pa_global[j_1].x - w[i_1].g.xS) < epsx) &&
						(pa_global[j_1].z < w[i_1].g.zE + epsz) &&
						(pa_global[j_1].z > w[i_1].g.zS - epsz) &&
						(pa_global[j_1].y > w[i_1].g.yS - epsy) && 
						(pa_global[j_1].y < w[i_1].g.yE + epsy)) {
						bfound = true;
					}
					break;
				case XZ_PLANE:
					if ((fabs(pa_global[j_1].y - w[i_1].g.yS) < epsy) &&
						(pa_global[j_1].z < w[i_1].g.zE + epsz) &&
						(pa_global[j_1].z > w[i_1].g.zS - epsz) && 
						(pa_global[j_1].x > w[i_1].g.xS - epsx) &&
						(pa_global[j_1].x < w[i_1].g.xE + epsx)) {
						bfound = true;
					}
					break;
				}
				if (bfound) {

					bwall_active[j_1] = true;

					boundary[j_1] = true;
					// Фиксированный потенциал.
					if (w[i_1].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) {
						//constr[j_1] = true;
						//temp_potent[j_1] = w[i_1].Tamb;
						rthdsd[j_1] = w[i_1].Tamb;
						//printf("%d Tamb=%e\n",i_1, w[i_1].Tamb);
						//getchar();
					}

					// граничное условие Стефана - Больцмана
					if (w[i_1].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) {
						// 16.12.2018

						breakRUMBAcalc_for_nonlinear_boundary_condition = true;
						doublereal alpha_relax1 = 0.25;

						if (temp_potent[j_1] > w[i_1].Tamb) {

							if (temp_potent[j_1] < -272.15) {
								temp_potent[j_1] = -272.15;
							}

							// Температура с предыдущего шага по времени !!!. Работает только в нестационарном солвере.

							//rthdsd[j_1] = alpha_relax1 * (-w[i_1].emissivity*STEFAN_BOLCMAN_CONST*((273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1]) - (273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb))) +
								//(1.0 - alpha_relax1)*(-w[i_1].emissivity*STEFAN_BOLCMAN_CONST*((273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1]) - (273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)*(273.15 + w[i_1].Tamb)));

							rthdsd[j_1] = (-w[i_1].emissivity*w[i_1].ViewFactor*STEFAN_BOLCMAN_CONST*
								   ((273.15 + temp_potent[j_1])*
								    (273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1]) -
									(273.15 + w[i_1].Tamb)*
									(273.15 + w[i_1].Tamb)*
									(273.15 + w[i_1].Tamb)*
									(273.15 + w[i_1].Tamb)));					


							rthdsd[j_1] *= dsquare_default_non_linear1[j_1];

						}
					} // Stefan-Bolcman

					// Граничное условие Ньютона-Рихмана.
					if (w[i_1].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY) {
						// 16.12.2018

						breakRUMBAcalc_for_nonlinear_boundary_condition = true;

						rthdsd[j_1] = -w[i_1].film_coefficient*(temp_potent[j_1] - w[i_1].Tamb);						

						rthdsd[j_1] *=  dsquare_default_non_linear1[j_1];

					} // Newton-Richman
				}
			}
		}

		
		
		if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {

				for (integer j_1 = 0; j_1 < maxelm_global; j_1++) {
					if (((icount_number_visit[j_1] >= 2)) && (icount_number_visit[j_1] <= 4) && (bwall_active[j_1] == false)) {
						
						if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {


							breakRUMBAcalc_for_nonlinear_boundary_condition = true;

							boundary[j_1] = true;

							// film_coefficient - глобальный.
							// operating_temperature_for_film_coeff - глобальный.
							rthdsd[j_1] = -film_coefficient*(temp_potent[j_1] - operating_temperature_for_film_coeff);


							// Определение Нормали к стенке. 
							// Работает на АЛИС сетке тоже.
							//Без асемблесов.

							
							// Площадь была вычислена один раз и запомнена.
							rthdsd[j_1] *= dsquare_default_non_linear[j_1];

						}  // Newton-Richman

						if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) {

							breakRUMBAcalc_for_nonlinear_boundary_condition = true;
							doublereal alpha_relax1 = 0.25;

							boundary[j_1] = true;

							//if ((temp_potent[j_1] < -271.0) && (temp_potent[j_1] <= operating_temperature_for_film_coeff)) {
								// Защита от опускания температуры ниже абсолютного нуля.
								// 7.12.2019
								// qb=0.0;
							//}
							//else 
							{

								if (temp_potent[j_1] < -272.15) {
									temp_potent[j_1] = -272.15;
								}
								if (temp_potent_old[j_1] < -272.15) {
									temp_potent_old[j_1] = -272.15;
								}

								// Температура с предыдущего шага по времени !!!. Работает только в нестационарном солвере.

								//rthdsd[j_1] = alpha_relax1 * (-emissivity_default_non_linear[j_1]*STEFAN_BOLCMAN_CONST*((273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1]) - (273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff))) +
								//(1.0 - alpha_relax1)*(-emissivity_default_non_linear[j_1]*STEFAN_BOLCMAN_CONST*((273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1]) - (273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)));
								
								rthdsd[j_1] = (-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1]) -
									(273.15 + operating_temperature_for_film_coeff)*
									(273.15 + operating_temperature_for_film_coeff)*
									(273.15 + operating_temperature_for_film_coeff)*
									(273.15 + operating_temperature_for_film_coeff)));
									
								/*
								rthdsd[j_1] = alpha_relax1*(-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1]) -
									(273.15 + operating_temperature_for_film_coeff)*
									(273.15 + operating_temperature_for_film_coeff)*
									(273.15 + operating_temperature_for_film_coeff)*
									(273.15 + operating_temperature_for_film_coeff)))+
									(1.0-alpha_relax1)*(-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
									    (273.15 + temp_potent_old[j_1])*
										(273.15 + temp_potent_old[j_1])*
										(273.15 + temp_potent_old[j_1])*
										(273.15 + temp_potent_old[j_1]) -
										(273.15 + operating_temperature_for_film_coeff)*
										(273.15 + operating_temperature_for_film_coeff)*
										(273.15 + operating_temperature_for_film_coeff)*
										(273.15 + operating_temperature_for_film_coeff)));
										*/

								// Определение Нормали к стенке. 
								// Работает на АЛИС сетке тоже.
								//Без асемблесов.


								// Площадь была вычислена один раз и запомнена.
								rthdsd[j_1] *= dsquare_default_non_linear[j_1];

							}
						} // Stefan-Bolcman

						if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC) {

							boundary[j_1] = true;

							breakRUMBAcalc_for_nonlinear_boundary_condition = true;
							doublereal alpha_relax1 = 0.25;

							//if ((temp_potent[j_1] < -271.0) && (temp_potent[j_1] <= operating_temperature_for_film_coeff)) {
							// Защита от опускания температуры ниже абсолютного нуля.
							// 7.12.2019
							// qb=0.0;
							//}
							//else 
							{

								if (temp_potent[j_1] < -272.15) {
									temp_potent[j_1] = -272.15;
								}
								if (temp_potent_old[j_1] < -272.15) {
									temp_potent_old[j_1] = -272.15;
								}

								// Температура с предыдущего шага по времени !!!. Работает только в нестационарном солвере.

								//rthdsd[j_1] = alpha_relax1 * (-emissivity_default_non_linear[j_1]*STEFAN_BOLCMAN_CONST*((273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1])*(273.15 + temp_potent[j_1]) - (273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff))) +
								//(1.0 - alpha_relax1)*(-emissivity_default_non_linear[j_1]*STEFAN_BOLCMAN_CONST*((273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1])*(273.15 + toldtimestep[j_1]) - (273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)*(273.15 + operating_temperature_for_film_coeff)));

								rthdsd[j_1] = -film_coefficient*(temp_potent[j_1] - operating_temperature_for_film_coeff);
								rthdsd[j_1] += (-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1])*
									(273.15 + temp_potent[j_1]) -
									(273.15 + operating_temperature_for_film_coeff)*
									(273.15 + operating_temperature_for_film_coeff)*
									(273.15 + operating_temperature_for_film_coeff)*
									(273.15 + operating_temperature_for_film_coeff)));

								/*
								rthdsd[j_1] = alpha_relax1*(-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
								(273.15 + temp_potent[j_1])*
								(273.15 + temp_potent[j_1])*
								(273.15 + temp_potent[j_1])*
								(273.15 + temp_potent[j_1]) -
								(273.15 + operating_temperature_for_film_coeff)*
								(273.15 + operating_temperature_for_film_coeff)*
								(273.15 + operating_temperature_for_film_coeff)*
								(273.15 + operating_temperature_for_film_coeff)))+
								(1.0-alpha_relax1)*(-emissivity_default_non_linear[j_1] * STEFAN_BOLCMAN_CONST*(
								(273.15 + temp_potent_old[j_1])*
								(273.15 + temp_potent_old[j_1])*
								(273.15 + temp_potent_old[j_1])*
								(273.15 + temp_potent_old[j_1]) -
								(273.15 + operating_temperature_for_film_coeff)*
								(273.15 + operating_temperature_for_film_coeff)*
								(273.15 + operating_temperature_for_film_coeff)*
								(273.15 + operating_temperature_for_film_coeff)));
								*/

								// Определение Нормали к стенке. 
								// Работает на АЛИС сетке тоже.
								//Без асемблесов.


								// Площадь была вычислена один раз и запомнена.
								rthdsd[j_1] *= dsquare_default_non_linear[j_1];

							}
						} // Mix-Condition

					}
				}
			
		}
		// Инициализация.
		for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
			bwall_active[i_1] = false;
		}


	    temp_max = -1.0e30;
	    integer imax, imin;
		temp_min = +1.0e30;
	    for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
		    if (temp_potent[i_1] > temp_max) {
			    temp_max = temp_potent[i_1];
			    imax = i_1;
		    }
		    if (temp_potent[i_1] < temp_min) {
			    temp_min = temp_potent[i_1];
			    imin = i_1;
		    }
	    }
	    told_iter = temp_max;
	    printf("temperature: min = %e, max=%e\n", temp_min, temp_max);

	    // Обязательно после сборки граничных условий.
		for (integer j_1 = 0; j_1 < maxelm_global; j_1++) {
			temp_potent_old[j_1] = temp_potent[j_1];
		}
		
		// инициализация.
		for (integer i_1 = 0; i_1 < 8; i_1++) {
			for (integer j_1 = 0; j_1 < 8; j_1++) {
				Kmatrix_local[i_1][j_1] = 0.0;
			}
		}

		IMatrix sparseS; // Разреженная матрица в формате IMatrix.
		//if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == DIRECT_SECOND_T_SOLVER) {
		initIMatrix(&sparseS, maxelm_global);
		//	}
		SIMPLESPARSE sparseM; // разреженная матрица
		initsimplesparse(sparseM, maxelm_global);


		//printf("t.maxelm=%lld ncell_shadow_gl=%lld ", t.maxelm, ncell_shadow_gl); getchar();

		// Сборка СЛАУ.
		for (integer ie = 0; ie < ncell_shadow_gl; ie++) {
			// инициализация.
			for (integer i_1 = 0; i_1 < 8; i_1++) {
				for (integer j_1 = 0; j_1 < 8; j_1++) {
					Kmatrix_local[i_1][j_1] = 0.0;
				}
			}
			// Сборка локальной матрицы жёсткости.
			// Термоупругость сборка матрицы Жёсткости для шестигранной призмы. 19.05.2018.
			
			//if (lu == 0) {
				// структурированная сетка.
				Thermal_ALICE_assemble(ie, nvtx_global,
					pa_global, prop_global, Kmatrix_local,
					t.ptr, Ux_arr, Uy_arr, Uz_arr, mut_arr);
			//}
			/*else {
				Thermal_ALICE_assemble_old(ie, nvtx_global,
				pa_global, prop_global, Kmatrix_local);
			}*/

			/*
					for (integer i_1 = 0; i_1 < 8; i_1++) {
						for (integer j_1 = 0; j_1 < 8; j_1++) {
							printf("%e ",Kmatrix_local[i_1][j_1]);
						}
						printf("\n");
					}
					getchar();
					*/

			for (integer i_4 = 0; i_4 < 8; i_4++) {
				for (integer j_4 = 0; j_4 < 8; j_4++) {

					//Kmatrix_local[i_4][j_4] *= 4.0;// Компенсация 0.5 в матрице частных производных.
				}
			}





			// Сборка правой части 
			// TODO.

			// Собирать ли правую часть СЛАУ.
			bool bsecond_member_of_equation = true;
			// Добавление локальной матрицы жёсткости в глобальную.


			if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER) {
				// Будет использован прямой метод решения. 
			elembdSparse3(ie, sparseS, nvtx_global,
				constr, rthdsd,
				Kmatrix_local, temp_potent,
				bsecond_member_of_equation);
			}
			else {


			elembdSparse4(ie, sparseM, nvtx_global,
				constr, rthdsd,
				Kmatrix_local, temp_potent,
				bsecond_member_of_equation);
			}


		}

		

		if (!(iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER)) {
			for (integer i_check = 0; i_check < maxelm_global; i_check++) {
				if (sparseM.root[i_check] == nullptr) {
					printf("error: zero string %lld \n", i_check);
					system("pause");
				}
				else {
					NONZEROELEM* p;
					p = sparseM.root[i_check];
					if (p != nullptr) {
						NONZEROELEM* q = nullptr;
						//printf("%e %d %d\n", p->aij, i_check, p->key);
						//getchar();

						if (p->key == i_check) {
							if (!constr[i_check]) {
								// диагональный элемент.
								doublereal ap0 = 0.0;
								ap0 = rho_export[i_check] * Cp_export[i_check] * Vol_export[i_check] / timestep_seq_current_step;

								if (btimedep) {
									p->aij += ap0;
									rthdsd[i_check] += ap0 * toldtimestep[i_check];
								}
							}


						}

						q = p->next;
						//p->next = nullptr;

						while (q != nullptr) {
							p = q;
							if (p->key == i_check) {
								if (!constr[i_check]) {
									// диагональный элемент.
									doublereal ap0 = 0.0;
									ap0 = rho_export[i_check] * Cp_export[i_check] * Vol_export[i_check] / timestep_seq_current_step;
									if (btimedep) {
										p->aij += ap0;
										rthdsd[i_check] += ap0 * toldtimestep[i_check];
									}
								}

							}

							//printf("%e %d %d\n", p->aij, i_check, p->key);
							//getchar();
							if (fabs(p->aij) < 1.0e-15) {
								if (p->key == i_check) {
									printf("%e %lld %lld\n", p->aij, i_check, p->key);
									system("PAUSE");
								}
							}

							//printf(" Dirichlet p-aij=%d\n",p->aij);
							//getchar();
							q = p->next;
							//p->next = nullptr;
							//delete p;
							p = nullptr;
							//M.n--;
						}
						//delete M.root[i];
						//M.root[i] = nullptr;
						//M.n--;
					}
				}
			}
		}


		
		/*
		for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
			if ((fabs(t.Sc[i_1]) > 1.0e-24)&&(boundary[i_1])) {
				if (btimedep) {
					printf("ERROR!!! unsteady modeling not realysation\n.");
					system("PAUSE");
					exit(1);
				}
				else {
					rthdsd[i_1] = 0.0;// Эта ячейка в которой задаётся мощность тепловыделения. 
					// Она не может быть одновременно излучающей.
					boundary[i_1] = false;
					printf("incomming\n");
					getchar();
				}
			}
		}*/

		doublereal Pdiss_actual = 0.0;
		// Объёмные источники тепла.
		for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
			// вычисление размеров текущего контрольного объёма:
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
			volume3D(i_1, nvtx_global, pa_global, dx, dy, dz);

			if (fabs(t.Sc[i_1]) > 1.0e-30) {
				for (integer j = 0; j <= 7; j++) {
					if (!constr[nvtx_global[j][i_1] - 1]) {
						if (i_1 < t.maxelm) {
							rthdsd[nvtx_global[j][i_1] - 1] +=
								0.125*dx*dy*dz*t.Sc[i_1] * poweron_multiplier_sequence;
							Pdiss_actual += 0.125*dx*dy*dz*t.Sc[i_1];
							//if (t.Sc[i_1] > 1.0e-30) {
								//printf("1 %e \n", 0.125 * dx * dy * dz * t.Sc[i_1] * poweron_multiplier_sequence);
								//getchar();
							//}
						}
					}
					else {
						printf("ERROR!!! 23.08.2020.  zanichenie power Pdiss.\n");
						system("PAUSE");
					}
				}
			}

		}
		
		integer ic_now = t.maxelm;
		for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
			for (integer i_1 = 0; i_1 < my_union[iu_74].t.maxelm; i_1++) {
				// вычисление размеров текущего контрольного объёма:
				doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
				volume3D(i_1, my_union[iu_74].t.nvtx, my_union[iu_74].t.pa, dx, dy, dz);

				for (integer j = 0; j <= 7; j++) {
					if (!constr[hash_table_pa[iu_74][my_union[iu_74].t.nvtx[j][i_1] - 1]]) {
						if (i_1 < my_union[iu_74].t.maxelm) {
							rthdsd[hash_table_pa[iu_74][my_union[iu_74].t.nvtx[j][i_1] - 1]] +=
								0.125*dx*dy*dz*my_union[iu_74].t.Sc[i_1] * poweron_multiplier_sequence;
							Pdiss_actual += 0.125*dx*dy*dz*my_union[iu_74].t.Sc[i_1];
							//if (my_union[iu_74].t.Sc[i_1] > 1.0e-30) {
								//printf("2 %e \n", 0.125 * dx * dy * dz * my_union[iu_74].t.Sc[i_1] * poweron_multiplier_sequence);
								//getchar();
							//}
						}
					}
				}
			}
		}
		printf("Pdiss actual value is equal=%e\n", Pdiss_actual);
		/*
			doublereal* square = new doublereal[3 * t.maxnod];
			for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1++) {
				square[i_1] = 0.0;
			}
			// Вычисляем площадь.
			for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
				doublereal hx = 1.0, hy = 1.0, hz = 1.0;
				volume3D(i_1, t.nvtx, t.pa, hx, hy, hz);
				for (integer j = 0; j <= 7; j++) {
					integer j_1 = t.nvtx[j][i_1] - 1;
					//X
					square[3 * j_1] += 0.25*hy*hz;
					//Y
					square[3 * j_1 + 1] += 0.25*hx*hz;
					//Z
					square[3 * j_1 + 2] += 0.25*hx*hy;
				}
			}
			*/



		for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
			//rthdsd[i_1] *= 1.0e-6;
		}




		// Умножаем силу на площадь.
		// Сила линейного теплового расширения E*vol*betaT*gradT или
		// E*Square_ortho*betaT*DeltaT, DeltaT=l*gradT.
		for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
			//rthdsd[i_1] *= square[i_1]; // Newton*m!2.
		}
		//delete[] square;


		printf("matrix is assemble.\n");
		//getchar();
		// Решение СЛАУ TODO.

		if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER) {
			// Прямой метод решения.
			calculateSPARSEgaussArray(&sparseS, temp_potent, rthdsd);
		}
		bool bprintmessage = true;
		integer maxiter = 20000; // !!!
								 //ICCG(TOTALDEFORMATIONVAR, sparseM, rthdsd, deformation, 3 * t.maxnod, bprintmessage, false, maxiter); //->//
								 //doublereal *val = nullptr;
								 //integer *col_ind = nullptr, *row_ptr = nullptr;
								 //simplesparsetoCRS(sparseM, val, col_ind, row_ptr, (3 * t.maxnod)); // преобразование матрицы из одного формата хранения в другой.
								 //simplesparsefree(sparseM, 3 * t.maxnod);

								 // Разрешающих свойств данного метода без предобуславливателя явно недостаточно.
								 //Bi_CGStabCRS((3 * t.maxnod), val, col_ind, row_ptr, rthdsd, deformation, maxiter);//->//

								 // BiCGStab + ILU6 сходимость есть.
		if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER) {
			// BiCGStab + ILU(lfil), lfil=1..6.
			Bi_CGStab_internal4(sparseM, (maxelm_global), rthdsd, temp_potent, maxiter, bprintmessage, m,w,lw,boundary,TEMP);			
		}
		// amg1r5 нет сходимости на задачи напряженно-деформированного состояния.
		//amg_loc_memory_Stress(sparseM, (3*t.maxnod), rthdsd, deformation, m);
		if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::CAMG_RUMBA_v0_14_SECOND_T_SOLVER) {
			my_agr_amg_loc_memory_Stress(sparseM, maxelm_global, rthdsd, temp_potent, m,b,lb,w,lw, boundary, TEMP);
		}
		if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::AMG1R5_SECOND_T_SOLVER) {
			// if (NONE_only_amg1r5==stabilization_amg1r5_algorithm) -> amg1r5 Руге и Штубена.
			// if (BiCGStab_plus_amg1r5==stabilization_amg1r5_algorithm) -> BiCGStab + amg1r5 Хенк Ван Дер Ворст + Руге и Штубен.
			// if (FGMRes_plus_amg1r5==stabilization_amg1r5_algorithm) -> FGMres + amg1r5 Ю. Саад и Мартин Шульц + Руге и Штубен. // 16.10.2018
			amg_loc_memory_for_Matrix_assemble2(sparseM, (maxelm_global), rthdsd, temp_potent, maxiter, bprintmessage, m, w, lw, boundary);//13.10.2018
		}
		if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::AMGCL_SECONT_T_SOLVER)
		{
#if AMGCL_INCLUDE_IN_MY_PROJECT == 1
			// AMGCL Denis Demidov
			// 20.11.2019
			amgcl_secondT_solver(sparseM, (maxelm_global), rthdsd, temp_potent, bprintmessage,w,lw,boundary,false);
#endif
		}
		// Нужна специальная версия BicgStab+ILU2.

		

		printf("SLAU is solve.\n");
		//getchar();

		// Освобождение оперативной памяти.
		if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER) {
			freeIMatrix(&sparseS);
		}
		//simplesparsefree(sparseM, 3 * t.maxnod);
		
		doublereal alpha = 1.0; 

	    // Это не специальная нелинейная версия кода amgcl CAMG.
		for (integer k = 0; k < lw; k++) {
			if ((w[k].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) ||
				(w[k].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY)) {
				alpha = 0.2; // Для того чтобы решение СЛАУ сходилось.
				
			}
		}

		if ((adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) ||
			(adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC)) {
			alpha = 0.2; // Для того чтобы решение СЛАУ сходилось.
		}


		for (integer j_1 = 0; j_1 < maxelm_global; j_1++) {
			// Нижняя релаксация.
			if (temp_potent[j_1] < -272.15) {
				temp_potent[j_1] = -272.15;				
			}
			else {
				temp_potent[j_1] = temp_potent_old[j_1] + /*0.25; 1.0*/alpha *(temp_potent[j_1] - temp_potent_old[j_1]);
			}
			if (temp_potent[j_1] <= -272.15) {
				temp_potent[j_1] = -272.15;
			}
			
			if (btimedep) {
				// Результат вычисления.
				tnewtimestep[j_1]= temp_potent[j_1];
			}
		}

		temp_max = -1.0e30;
		temp_min = +1.0e30;
		for (integer i_1 = 0; i_1 < maxelm_global; i_1++) {
			if (temp_potent[i_1] > temp_max) temp_max = temp_potent[i_1];
			if (temp_potent[i_1] < temp_min) temp_min = temp_potent[i_1];
		}
		printf("temperature: min = %e, max=%e\n", temp_min, temp_max);

		

		simplesparsefree(sparseM, maxelm_global); // Очистка памяти из под матрицы sparseM.
		freeIMatrix(&sparseS);

		printf("prohod=%lld %e %e\n", iprohod, told_iter, temp_max);

	} while ((iprohod<=2)||(fabs(told_iter - temp_max) > 0.05*fabs(temp_max-temp_min)));
	
	


	delete[] temp_potent_old;
	
	delete[] boundary;

	delete[] bwall_active;
	delete[] icount_number_visit;
	delete[] dsquare_default_non_linear;
	delete[] dsquare_default_non_linear1;
	delete[] emissivity_default_non_linear;

	if (Kmatrix_local != nullptr) {
		for (integer i_1 = 0; i_1 < 8; i_1++) {
			if (Kmatrix_local[i_1] != nullptr) {
				delete[] Kmatrix_local[i_1];
				Kmatrix_local[i_1] = nullptr;
			}
		}
		delete[] Kmatrix_local;
		Kmatrix_local = nullptr;
	}

	/*
	// Запись результата для визуализации.
	if (t.total_deformation == nullptr) {
		t.total_deformation = new doublereal*[SIZE_DEFORMATION_ARRAY];
		for (integer j_6 = 0; j_6 < SIZE_DEFORMATION_ARRAY; j_6++) {
			t.total_deformation[j_6] = nullptr;
			if (t.total_deformation[j_6] == nullptr) {
				t.total_deformation[j_6] = new doublereal[t.maxelm + t.maxbound];
			}
		}

	}
	else {
		for (integer j_6 = 0; j_6 < SIZE_DEFORMATION_ARRAY; j_6++) {
			delete[] t.total_deformation[j_6];
			t.total_deformation[j_6] = nullptr;
		}
		delete[] t.total_deformation;
		t.total_deformation = nullptr;

		t.total_deformation = new doublereal*[SIZE_DEFORMATION_ARRAY];
		for (integer j_6 = 0; j_6 < SIZE_DEFORMATION_ARRAY; j_6++) {
			t.total_deformation[j_6] = nullptr;
			if (t.total_deformation[j_6] == nullptr) {
				t.total_deformation[j_6] = new doublereal[t.maxelm + t.maxbound];
			}
		}
	}
	for (integer j_6 = 0; j_6 < SIZE_DEFORMATION_ARRAY; j_6++) {
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			t.total_deformation[j_6][i_1] = 0.0;
		}
	}
	*/
	/*	
	if (1) {

		// Метод линейного порядка.
		doublereal min_x = 1e60;
		doublereal min_y = 1e60;
		doublereal min_z = 1e60;
		doublereal max_x = -1e60;
		doublereal max_y = -1e60;
		doublereal max_z = -1e60;

		for (integer i = 0; i < t.maxnod; i++) {
			if (t.pa[i].x < min_x) {
				min_x = t.pa[i].x;
			}
			if (t.pa[i].y < min_y) {
				min_y = t.pa[i].y;
			}
			if (t.pa[i].z < min_z) {
				min_z = t.pa[i].z;
			}
			if (t.pa[i].x > max_x) {
				max_x = t.pa[i].x;
			}
			if (t.pa[i].y > max_y) {
				max_y = t.pa[i].y;
			}
			if (t.pa[i].z > max_z) {
				max_z = t.pa[i].z;
			}
		}

		//min_x *= 1.2;
		//min_y *= 1.2;
		//min_z *= 1.2;



		min_x = 1.05*fabs(max_x - min_x);
		if (min_x < 1.0e-30) {
			min_x = 1.05*fabs(max_x);
		}
		min_y = 1.05*fabs(max_y - min_y);
		if (min_y < 1.0e-30) {
			min_y = 1.05*fabs(max_y);
		}
		min_z = 1.05*fabs(max_z - min_z);
		if (min_z < 1.0e-30) {
			min_z = 1.05*fabs(max_z);
		}

		*/
		/*
		if (min_x < 1.0e-30) {
		printf("error!!! negative min_x MNK!\n");
		printf("min_x=%e max_x=%e\n",min_x,max_x);
		}
		if (min_y < 1.0e-30) {
		printf("error!!! negative min_y MNK!\n");
		printf("min_y=%e max_y=%e\n", min_y, max_y);
		}
		if (min_z < 1.0e-30) {
		printf("error!!! negative min_z MNK!\n");
		printf("min_z=%e max_z=%e\n", min_z, max_z);
		}
		*/
	   /*
		TOCHKA** pointerlist = new TOCHKA*[t.maxelm];
		doublereal** rthdsd_Gauss = new doublereal*[t.maxelm];
		for (integer i_47 = 0; i_47 < t.maxelm; i_47++) {
			pointerlist[i_47] = new TOCHKA[8];
			rthdsd_Gauss[i_47] = new doublereal[8];
		}

		for (integer j_6 = 0; j_6 < 4; j_6++) {

			for (integer i = 0; i < t.maxelm; i++) {
				//doublereal xc47, yc47, zc47;

				TOCHKA p;
				center_cord3D(i, t.nvtx, t.pa, p, 100);
				//xc47 = p.x;
				//yc47 = p.y;
				//zc47 = p.z;


				p.x = p.x + min_x;
				p.y = p.y + min_y;
				p.z = p.z + min_z;

				for (integer j = 0; j <= 7; j++) {
					TOCHKA p1;
					p1.x = t.pa[t.nvtx[j][i] - 1].x;
					p1.y = t.pa[t.nvtx[j][i] - 1].y;
					p1.z = t.pa[t.nvtx[j][i] - 1].z;
					p1.x = p1.x + min_x;
					p1.y = p1.y + min_y;
					p1.z = p1.z + min_z;

					pointerlist[i][j] = p1;
					if (fabs(p1.x) < 1.0e-40) {
						printf("problem x=%e\n", p1.x);
						getchar();
					}
					if (fabs(p1.y) < 1.0e-40) {
						printf("problem y=%e\n", p1.y);
						getchar();
					}
					if (fabs(p1.z) < 1.0e-40) {
						printf("problem z=%e\n", p1.z);
						getchar();
					}
					integer j_1 = t.nvtx[j][i] - 1;
					switch (j_6) {
					case 0:// TOTAL DEFORMATION
						rthdsd_Gauss[i][j] = temp_potent[j_1];

						break;
					case 1: // X deformation
						rthdsd_Gauss[i][j] = temp_potent[j_1];

						break;
					case 2: // Y deformation
						rthdsd_Gauss[i][j] = temp_potent[j_1];

						break;
					case 3: // Z deformation
						rthdsd_Gauss[i][j] = temp_potent[j_1];

						break;
					}

				}


				doublereal** Xmatr = new doublereal*[4];
				for (integer j = 0; j <= 3; j++) {
					Xmatr[j] = new doublereal[4];
				}


				doublereal* bmatr = new doublereal[4];
				doublereal* koefmatr = new doublereal[4];

				for (integer j1 = 0; j1 <= 3; j1++) {
					for (integer j2 = 0; j2 <= 3; j2++) {
						Xmatr[j1][j2] = 0.0;
					}
					bmatr[j1] = 0.0;
					koefmatr[j1] = 0.0;
				}




				for (integer j = 0; j < 8; j++) {

					Xmatr[0][0] += 1.0;
					Xmatr[0][1] += pointerlist[i][j].x;
					Xmatr[0][2] += pointerlist[i][j].y;
					Xmatr[0][3] += pointerlist[i][j].z;

					Xmatr[1][0] += pointerlist[i][j].x;
					Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
					Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
					Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;

					Xmatr[2][0] += pointerlist[i][j].y;
					Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
					Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
					Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;

					Xmatr[3][0] += pointerlist[i][j].z;
					Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
					Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
					Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;

					bmatr[0] += rthdsd_Gauss[i][j];
					bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
					bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
					bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];
				}


				for (integer j1 = 0; j1 <= 100; j1++) {
					koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
					koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
					koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
					koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
				}
				//t.total_deformation[j_6][i] = koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z);
				// Температура.
				t.potent[i]= 2.0*(koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));


				for (integer j = 0; j <= 3; j++) {
					delete[] Xmatr[j];
				}
				delete[] Xmatr;
				delete[] bmatr;
				delete[] koefmatr;

			}
		} // j_6

		for (integer i = 0; i < t.maxelm; i++) {
			delete[] pointerlist[i];
			delete[] rthdsd_Gauss[i];
		}
		delete[] pointerlist;
		delete[] rthdsd_Gauss;

	}
	*/
	// Сохранение деформации.
	// TODO.
    delete[] t_for_Mechanical;
    t_for_Mechanical = new doublereal[t.maxnod + 2];
	for (integer i = 0; i < t.maxnod; i++) {
		// Сохранение температуры в узлах для механики.
		t_for_Mechanical[i] = temp_potent[i];
	}

    // С асемблесами репорт не работает корректно.

    // Без учёта асемблесов 25.07.2020.
	for (integer i = 0; i < t.maxelm; i++) {
		t.potent[i] = 0.0;
		for (integer j = 0; j < 8; j++) {
			t.potent[i]+=temp_potent[nvtx_global[j][i] - 1];
		}
		t.potent[i] /= 8.0;
	}
	for (integer i = 0; i < t.maxbound; i++) {
		if ((t.border_neighbor[i].MCB < (ls + lw)) &&
			(t.border_neighbor[i].MCB >= ls) &&
			(w[t.border_neighbor[i].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY)) {
			t.potent[i + t.maxelm] = w[t.border_neighbor[i].MCB - ls].Tamb; // На стенке точное значение.
		}
		else {
			integer iI = t.border_neighbor[i].iI;
			doublereal dc = 0.0;
			//t.potent[i + t.maxelm] = t.potent[iI]; // из ближайшего внутреннего узла.
			switch (t.border_neighbor[i].Norm) { // Внутренняя нормаль.
			case E_SIDE: t.potent[i + t.maxelm] = 0.0;
				for (integer j = 0; j < 8; j++) {
					if (fabs(t.pa[nvtx_global[j][iI] - 1].x - t.pa[nvtx_global[0][iI] - 1].x) < 1.0e-24) {
						t.potent[i + t.maxelm] += temp_potent[nvtx_global[j][iI] - 1];
						dc += 1.0;
					}
				}
				if (fabs(dc) < 0.5) {
					printf("E_SIDE division by zero in solve_Thermal function...\n");
					system("PAUSE");
					exit(1);
				}
				t.potent[i + t.maxelm] /= dc;
				break;
			case W_SIDE: t.potent[i + t.maxelm] = 0.0;
				for (integer j = 0; j < 8; j++) {
					if (fabs(t.pa[nvtx_global[j][iI] - 1].x - t.pa[nvtx_global[1][iI] - 1].x) < 1.0e-24) {
						t.potent[i + t.maxelm] += temp_potent[nvtx_global[j][iI] - 1];
						dc += 1.0;
					}
				}
				if (fabs(dc) < 0.5) {
					printf("W_SIDE division by zero in solve_Thermal function...\n");
					system("PAUSE");
					exit(1);
				}
				t.potent[i + t.maxelm] /= dc;
				break;
			case N_SIDE: t.potent[i + t.maxelm] = 0.0;
				for (integer j = 0; j < 8; j++) {
					if (fabs(t.pa[nvtx_global[j][iI] - 1].y - t.pa[nvtx_global[0][iI] - 1].y) < 1.0e-24) {
						t.potent[i + t.maxelm] += temp_potent[nvtx_global[j][iI] - 1];
						dc += 1.0;
					}
				}
				if (fabs(dc) < 0.5) {
					printf("N_SIDE division by zero in solve_Thermal function...\n");
					system("PAUSE");
					exit(1);
				}
				t.potent[i + t.maxelm] /= dc;
				break;
			case S_SIDE: t.potent[i + t.maxelm] = 0.0;
				for (integer j = 0; j < 8; j++) {
					if (fabs(t.pa[nvtx_global[j][iI] - 1].y - t.pa[nvtx_global[2][iI] - 1].y) < 1.0e-24) {
						t.potent[i + t.maxelm] += temp_potent[nvtx_global[j][iI] - 1];
						dc += 1.0;
					}
				}
				if (fabs(dc) < 0.5) {
					printf("S_SIDE division by zero in solve_Thermal function...\n");
					system("PAUSE");
					exit(1);
				}
				t.potent[i + t.maxelm] /= dc;
				break;
			case T_SIDE: t.potent[i + t.maxelm] = 0.0;
				for (integer j = 0; j < 8; j++) {
					if (fabs(t.pa[nvtx_global[j][iI] - 1].z - t.pa[nvtx_global[0][iI] - 1].z) < 1.0e-24) {
						t.potent[i + t.maxelm] += temp_potent[nvtx_global[j][iI] - 1];
						dc += 1.0;
					}
				}
				if (fabs(dc) < 0.5) {
					printf("T_SIDE division by zero in solve_Thermal function...\n");
					system("PAUSE");
					exit(1);
				}
				t.potent[i + t.maxelm] /= dc;
				break;
			case B_SIDE: t.potent[i + t.maxelm] = 0.0;
				for (integer j = 0; j < 8; j++) {
					if (fabs(t.pa[nvtx_global[j][iI] - 1].z - t.pa[nvtx_global[4][iI] - 1].z) < 1.0e-24) {
						t.potent[i + t.maxelm] += temp_potent[nvtx_global[j][iI] - 1];
						dc += 1.0;
					}
				}
				if (fabs(dc) < 0.5) {
					printf("B_SIDE division by zero in solve_Thermal function...\n");
					system("PAUSE");
					exit(1);
				}
				t.potent[i + t.maxelm] /= dc;
				break;
			}
		}
	}



	if (!btimedep) {

		// для преобразования из центра в края.
		//lam_for_export;
		//sum_vol; // суммарный объем.
		for (integer i_72 = 0; i_72 < maxelm_global; i_72++) {
			// инициализация.
			lam_for_export[i_72] = 0.0;
			sum_vol[i_72] = 0.0;
		}
		for (integer i_72 = 0; i_72 < ncell_shadow_gl; i_72++) {
			doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
			volume3D(i_72, nvtx_global, pa_global, hx, hy, hz);
			for (integer i_73 = 0; i_73 < 8; i_73++) {
				lam_for_export[nvtx_global[i_73][i_72] - 1] += hx * hy*hz*lam_export[i_72];
				sum_vol[nvtx_global[i_73][i_72] - 1] += hx * hy*hz;
			}
		}
		for (integer i_72 = 0; i_72 < maxelm_global; i_72++) {
			// инициализация.
			lam_for_export[i_72] = lam_for_export[i_72] / sum_vol[i_72];
		}

		// Вычисление тепловых потоков.
		doublereal* Tx = new doublereal[ncell_shadow_gl];
		doublereal* Ty = new doublereal[ncell_shadow_gl];
		doublereal* Tz = new doublereal[ncell_shadow_gl];
		doublereal* HeatFluxMag = new doublereal[ncell_shadow_gl];

		for (integer i_72 = 0; i_72 < ncell_shadow_gl; i_72++) {
			doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
			volume3D(i_72, nvtx_global, pa_global, hx, hy, hz);
			Tx[i_72] = -lam_export[i_72] * (temp_potent[nvtx_global[1][i_72] - 1] - temp_potent[nvtx_global[0][i_72] - 1]) / hx;
			Ty[i_72] = -lam_export[i_72] * (temp_potent[nvtx_global[3][i_72] - 1] - temp_potent[nvtx_global[0][i_72] - 1]) / hy;
			Tz[i_72] = -lam_export[i_72] * (temp_potent[nvtx_global[4][i_72] - 1] - temp_potent[nvtx_global[0][i_72] - 1]) / hz;
			HeatFluxMag[i_72] = sqrt(Tx[i_72] * Tx[i_72] + Ty[i_72] * Ty[i_72] + Tz[i_72] * Tz[i_72]);
		}

		doublereal* Txgl = new doublereal[maxelm_global];
		doublereal* Tygl = new doublereal[maxelm_global];
		doublereal* Tzgl = new doublereal[maxelm_global];
		doublereal* HeatFluxMaggl = new doublereal[maxelm_global];

		// Если есть асемблесы или сетка адаптивная локально измельченная.
		if ((lu > 0)/*||(b_on_adaptive_local_refinement_mesh)*/) {

			for (integer i_72 = 0; i_72 < maxelm_global; i_72++) {
				// инициализация.
				Txgl[i_72] = 0.0;
				Tygl[i_72] = 0.0;
				Tzgl[i_72] = 0.0;
				HeatFluxMaggl[i_72] = 0.0;
				sum_vol[i_72] = 0.0;
			}
			for (integer i_72 = 0; i_72 < ncell_shadow_gl; i_72++) {
				doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
				volume3D(i_72, nvtx_global, pa_global, hx, hy, hz);
				for (integer i_73 = 0; i_73 < 8; i_73++) {
					Txgl[nvtx_global[i_73][i_72] - 1] += hx * hy*hz*Tx[i_72];
					Tygl[nvtx_global[i_73][i_72] - 1] += hx * hy*hz*Ty[i_72];
					Tzgl[nvtx_global[i_73][i_72] - 1] += hx * hy*hz*Tz[i_72];
					HeatFluxMaggl[nvtx_global[i_73][i_72] - 1] += hx * hy*hz*HeatFluxMag[i_72];
					sum_vol[nvtx_global[i_73][i_72] - 1] += hx * hy*hz;
				}
			}
			for (integer i_72 = 0; i_72 < maxelm_global; i_72++) {
				// инициализация.
				Txgl[i_72] = Txgl[i_72] / sum_vol[i_72];
				Tygl[i_72] = Tygl[i_72] / sum_vol[i_72];
				Tzgl[i_72] = Tzgl[i_72] / sum_vol[i_72];
				HeatFluxMaggl[i_72] = HeatFluxMaggl[i_72] / sum_vol[i_72];
			}

		}
		else {
			if (lu == 0) {
				// Только для структурированной сетки.
				// Принятая нумерация в nvtx.
				// 2 3 | 6 7
				// 0 1 | 4 5
				integer nstop = t.maxelm;// ncell_shadow_gl
				if (t.neighbors_for_the_internal_node != nullptr) {
					for (integer i_72 = 0; i_72 < nstop; i_72++) {
						if ((t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE1 > -1) && (t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE1 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[1][i_72] - 1];
							hxplus = fabs(pp.x - pb.x);
							// почему то возникает сбой. Инструкция обратилась не по адресу TODO 20.05.2018
							hxminus = fabs(pp.x - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE1] - 1].x);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE1] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE1]);
							Txgl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[1][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE2 > -1) && (t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE2 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[1][i_72] - 1];
							hxplus = fabs(pp.x - pb.x);
							// почему то возникает сбой. Инструкция обратилась не по адресу TODO 20.05.2018
							hxminus = fabs(pp.x - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE2] - 1].x);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE2] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE2]);
							Txgl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[1][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE3 > -1) && (t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE3 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[1][i_72] - 1];
							hxplus = fabs(pp.x - pb.x);
							// почему то возникает сбой. Инструкция обратилась не по адресу TODO 20.05.2018
							hxminus = fabs(pp.x - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE3] - 1].x);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE3] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE3]);
							Txgl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[1][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE4 > -1) && (t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE4 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[1][i_72] - 1];
							hxplus = fabs(pp.x - pb.x);
							// почему то возникает сбой. Инструкция обратилась не по адресу TODO 20.05.2018
							hxminus = fabs(pp.x - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE4] - 1].x);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE4] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE4]);
							Txgl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[1][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[W_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(i_72, nvtx_global, pa_global, hx, hy, hz);
							Txgl[nvtx_global[0][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[1][i_72] - 1] - temp_potent[nvtx_global[0][i_72] - 1]) / hx;
							Txgl[nvtx_global[2][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[3][i_72] - 1] - temp_potent[nvtx_global[2][i_72] - 1]) / hx;
							Txgl[nvtx_global[4][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[5][i_72] - 1] - temp_potent[nvtx_global[4][i_72] - 1]) / hx;
							Txgl[nvtx_global[6][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[7][i_72] - 1] - temp_potent[nvtx_global[6][i_72] - 1]) / hx;
						}
						if ((t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE1 > -1) && (t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE1 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[1][i_72] - 1];
							hxminus = fabs(pp.x - pb.x);
							hxplus = fabs(pb.x - pa_global[nvtx_global[1][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE1] - 1].x);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE1] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE1]);
							Txgl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE2 > -1) && (t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE2 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[1][i_72] - 1];
							hxminus = fabs(pp.x - pb.x);
							hxplus = fabs(pb.x - pa_global[nvtx_global[1][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE2] - 1].x);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE2] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE2]);
							Txgl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE3 > -1) && (t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE3 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[1][i_72] - 1];
							hxminus = fabs(pp.x - pb.x);
							hxplus = fabs(pb.x - pa_global[nvtx_global[1][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE3] - 1].x);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE3] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE3]);
							Txgl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE4 > -1) && (t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE4 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[1][i_72] - 1];
							hxminus = fabs(pp.x - pb.x);
							hxplus = fabs(pb.x - pa_global[nvtx_global[1][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE4] - 1].x);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE4] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE4]);
							Txgl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.
							Txgl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[E_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(i_72, nvtx_global, pa_global, hx, hy, hz);
							Txgl[nvtx_global[1][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[1][i_72] - 1] - temp_potent[nvtx_global[0][i_72] - 1]) / hx;
							Txgl[nvtx_global[3][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[3][i_72] - 1] - temp_potent[nvtx_global[2][i_72] - 1]) / hx;
							Txgl[nvtx_global[5][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[5][i_72] - 1] - temp_potent[nvtx_global[4][i_72] - 1]) / hx;
							Txgl[nvtx_global[7][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[7][i_72] - 1] - temp_potent[nvtx_global[6][i_72] - 1]) / hx;
						}

						if ((t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE1 > -1) && (t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE1 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[2][i_72] - 1];
							hxplus = fabs(pp.y - pb.y);
							hxminus = fabs(pp.y - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE1] - 1].y);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE1] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE1]);
							Tygl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[2][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE2 > -1) && (t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE2 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[2][i_72] - 1];
							hxplus = fabs(pp.y - pb.y);
							hxminus = fabs(pp.y - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE2] - 1].y);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE2] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE2]);
							Tygl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[2][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE3 > -1) && (t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE3 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[2][i_72] - 1];
							hxplus = fabs(pp.y - pb.y);
							hxminus = fabs(pp.y - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE3] - 1].y);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE3] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE3]);
							Tygl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[2][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE4 > -1) && (t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE4 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[2][i_72] - 1];
							hxplus = fabs(pp.y - pb.y);
							hxminus = fabs(pp.y - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE4] - 1].y);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE4] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE4]);
							Tygl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[2][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[S_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(i_72, nvtx_global, pa_global, hx, hy, hz);
							Tygl[nvtx_global[0][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[2][i_72] - 1] - temp_potent[nvtx_global[0][i_72] - 1]) / hy;
							Tygl[nvtx_global[1][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[3][i_72] - 1] - temp_potent[nvtx_global[1][i_72] - 1]) / hy;
							Tygl[nvtx_global[4][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[6][i_72] - 1] - temp_potent[nvtx_global[4][i_72] - 1]) / hy;
							Tygl[nvtx_global[5][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[7][i_72] - 1] - temp_potent[nvtx_global[5][i_72] - 1]) / hy;
						}

						if ((t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE1 > -1) && (t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE1 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[2][i_72] - 1];
							hxminus = fabs(pp.y - pb.y);
							hxplus = fabs(pb.y - pa_global[nvtx_global[2][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE1] - 1].y);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE1] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE1]);
							Tygl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE2 > -1) && (t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE2 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[2][i_72] - 1];
							hxminus = fabs(pp.y - pb.y);
							hxplus = fabs(pb.y - pa_global[nvtx_global[2][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE2] - 1].y);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE2] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE2]);
							Tygl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE3 > -1) && (t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE3 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[2][i_72] - 1];
							hxminus = fabs(pp.y - pb.y);
							hxplus = fabs(pb.y - pa_global[nvtx_global[2][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE3] - 1].y);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE3] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE3]);
							Tygl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE4 > -1) && (t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE4 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[2][i_72] - 1];
							hxminus = fabs(pp.y - pb.y);
							hxplus = fabs(pb.y - pa_global[nvtx_global[2][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE4] - 1].y);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE4] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE4]);
							Tygl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.
							Tygl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[N_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(i_72, nvtx_global, pa_global, hx, hy, hz);
							Tygl[nvtx_global[2][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[2][i_72] - 1] - temp_potent[nvtx_global[0][i_72] - 1]) / hy;
							Tygl[nvtx_global[3][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[3][i_72] - 1] - temp_potent[nvtx_global[1][i_72] - 1]) / hy;
							Tygl[nvtx_global[6][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[6][i_72] - 1] - temp_potent[nvtx_global[4][i_72] - 1]) / hy;
							Tygl[nvtx_global[7][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[7][i_72] - 1] - temp_potent[nvtx_global[5][i_72] - 1]) / hy;
						}

						if ((t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE1 > -1) && (t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE1 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[4][i_72] - 1];
							hxplus = fabs(pp.z - pb.z);
							hxminus = fabs(pp.z - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE1] - 1].z);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE1] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE1]);
							Tzgl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[4][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE1] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE2 > -1) && (t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE2 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[4][i_72] - 1];
							hxplus = fabs(pp.z - pb.z);
							hxminus = fabs(pp.z - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE2] - 1].z);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE2] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE2]);
							Tzgl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[4][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE2] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE3 > -1) && (t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE3 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[4][i_72] - 1];
							hxplus = fabs(pp.z - pb.z);
							hxminus = fabs(pp.z - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE3] - 1].z);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE3] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE3]);
							Tzgl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[4][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE3] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE4 > -1) && (t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE4 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[4][i_72] - 1];
							hxplus = fabs(pp.z - pb.z);
							hxminus = fabs(pp.z - pa_global[nvtx_global[0][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE4] - 1].z);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE4] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE4]);
							Tzgl[nvtx_global[0][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[4][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[1][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[2][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[3][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[3][t.neighbors_for_the_internal_node[B_SIDE][i_72].iNODE4] - 1], temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(i_72, nvtx_global, pa_global, hx, hy, hz);
							Tzgl[nvtx_global[0][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[4][i_72] - 1] - temp_potent[nvtx_global[0][i_72] - 1]) / hz;
							Tzgl[nvtx_global[1][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[5][i_72] - 1] - temp_potent[nvtx_global[1][i_72] - 1]) / hz;
							Tzgl[nvtx_global[2][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[6][i_72] - 1] - temp_potent[nvtx_global[2][i_72] - 1]) / hz;
							Tzgl[nvtx_global[3][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[7][i_72] - 1] - temp_potent[nvtx_global[3][i_72] - 1]) / hz;
						}
						if ((t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE1 > -1) && (t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE1 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[4][i_72] - 1];
							hxminus = fabs(pp.z - pb.z);
							hxplus = fabs(pb.z - pa_global[nvtx_global[4][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE1] - 1].z);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE1] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE1]);
							Tzgl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE1] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE2 > -1) && (t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE2 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[4][i_72] - 1];
							hxminus = fabs(pp.z - pb.z);
							hxplus = fabs(pb.z - pa_global[nvtx_global[4][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE2] - 1].z);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE2] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE2]);
							Tzgl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE2] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE3 > -1) && (t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE3 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[4][i_72] - 1];
							hxminus = fabs(pp.z - pb.z);
							hxplus = fabs(pb.z - pa_global[nvtx_global[4][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE3] - 1].z);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE3] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE3]);
							Tzgl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE3] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else if ((t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE4 > -1) && (t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE4 < nstop)) {
							TOCHKA pp, pb;
							doublereal hxminus, hxplus;
							pp = pa_global[nvtx_global[0][i_72] - 1];
							pb = pa_global[nvtx_global[4][i_72] - 1];
							hxminus = fabs(pp.z - pb.z);
							hxplus = fabs(pb.z - pa_global[nvtx_global[4][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE4] - 1].z);
							doublereal lam_mix = 2.0*lam_export[i_72] * lam_export[t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE4] / (lam_export[i_72] + lam_export[t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE4]);
							Tzgl[nvtx_global[4][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[0][i_72] - 1], temp_potent[nvtx_global[4][i_72] - 1], temp_potent[nvtx_global[4][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[5][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[1][i_72] - 1], temp_potent[nvtx_global[5][i_72] - 1], temp_potent[nvtx_global[5][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[6][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[2][i_72] - 1], temp_potent[nvtx_global[6][i_72] - 1], temp_potent[nvtx_global[6][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.
							Tzgl[nvtx_global[7][i_72] - 1] = -lam_mix * rgradF(temp_potent[nvtx_global[3][i_72] - 1], temp_potent[nvtx_global[7][i_72] - 1], temp_potent[nvtx_global[7][t.neighbors_for_the_internal_node[T_SIDE][i_72].iNODE4] - 1], hxminus, hxplus); // второй порядок точности.

						}
						else {
							doublereal hx = 1.0, hy = 1.0, hz = 1.0; // размеры кубика
							volume3D(i_72, nvtx_global, pa_global, hx, hy, hz);
							Tzgl[nvtx_global[4][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[4][i_72] - 1] - temp_potent[nvtx_global[0][i_72] - 1]) / hz;
							Tzgl[nvtx_global[5][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[5][i_72] - 1] - temp_potent[nvtx_global[1][i_72] - 1]) / hz;
							Tzgl[nvtx_global[6][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[6][i_72] - 1] - temp_potent[nvtx_global[2][i_72] - 1]) / hz;
							Tzgl[nvtx_global[7][i_72] - 1] = -lam_export[i_72] * (temp_potent[nvtx_global[7][i_72] - 1] - temp_potent[nvtx_global[3][i_72] - 1]) / hz;
						}

					}
				}
				else {
					printf("t.neighbors_for_the_internal_node==nullptr in function solve_Thermal in module mysolverv0_03.c\n");
					system("PAUSE");
					exit(1);
				}
			}

			for (integer i_72 = 0; i_72 < maxelm_global; i_72++) {
				HeatFluxMaggl[i_72] = sqrt(Txgl[i_72] * Txgl[i_72] + Tygl[i_72] * Tygl[i_72] + Tzgl[i_72] * Tzgl[i_72]);
			}
		}

		// Запись для визуализации.
		// ncell_shadow_gl - количество ячеек - конечных элементов (прямых прямоугольных призм).
		// maxelm_global - количество узлов сетки.
		export_tecplot_temperature_ass(nvtx_global, bcheck_visible, pa_global, temp_potent, lam_for_export, Txgl, Tygl, Tzgl, HeatFluxMaggl, maxelm_global, ncell_shadow_gl);

		printf("temperature is writing.\n");
		//getchar();


		delete[] Txgl;
		delete[] Tygl;
		delete[] Tzgl;
		delete[] HeatFluxMaggl;
		delete[] Tx;
		delete[] Ty;
		delete[] Tz;
		delete[] HeatFluxMag;
	}

	// Освобождение оперативной памяти.
	delete[] Ux_arr;
	delete[] Uy_arr;
	delete[] Uz_arr;
	delete[] mut_arr;
	delete[] bcheck_visible;	
	delete[] lam_for_export;
	delete[] rho_export;
	delete[] Cp_export;
	delete[] Vol_export;
	delete[] sum_vol;
	delete[] lam_export;
	delete[] rthdsd;
	delete[] temp_potent;
	delete[] constr;

	
	for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
		delete[] hash_table_pa[iu_74];
	}
	delete[] hash_table_pa;

	// Освобождение памяти из под новейшей сборки учитывающей асемблесы.
	delete[] pa_global;
	for (integer i = 0; i<8; i++) { // -8N
		if (nvtx_global != nullptr) {
			if (nvtx_global[i] != nullptr) {
				delete[] nvtx_global[i];
				nvtx_global[i] = nullptr;
			}
		}
		
	}
	if (nvtx_global != nullptr) {
		delete[] nvtx_global;
		nvtx_global = nullptr;
	}
	if (prop_global != nullptr) {
		for (integer i = 0; i<8; i++) {
			if (prop_global[i] != nullptr) {
				delete[] prop_global[i]; // -3N
			}
		}
	}
	if (prop_global != nullptr) {
		delete[] prop_global;
		prop_global = nullptr;
	}
} // solve_Thermal

void Stress2Thermal_vector_translate(TEMPER &t, 
	doublereal* &input_vector_stress_3dim,
	integer id_translate,
	doublereal* &output_vector_Thermal_dim) {

	// t - вспомогательные данные: координаты узлов и пр. (pa, nvtx, maxelm, maxnod).
	// id_translate: 0-total, 1 - x def, 2 - y def, 3 - z def.


	if (1) {

		doublereal min_v = 1e60;
		doublereal max_v = -1e60;
		

		if ((id_translate == XDEFORMATION) || (id_translate == YDEFORMATION) || (id_translate == ZDEFORMATION)) {
			for (integer i_1 = 0; i_1 < t.maxnod; i_1++) {
				if (input_vector_stress_3dim[i_1] < min_v) {
					min_v = input_vector_stress_3dim[i_1];
				}
				if (input_vector_stress_3dim[i_1] > max_v) {
					max_v = input_vector_stress_3dim[i_1];
				}
			}
		}
		else if (id_translate == TOTALDEFORMATION) {
			for (integer j_1 = 0; j_1 < t.maxnod; j_1++) {
				doublereal td = sqrt(input_vector_stress_3dim[index_of(j_1, 'x')] * input_vector_stress_3dim[index_of(j_1, 'x')]
					+ input_vector_stress_3dim[index_of(j_1, 'y')] * input_vector_stress_3dim[index_of(j_1, 'y')] +
					input_vector_stress_3dim[index_of(j_1, 'z')] * input_vector_stress_3dim[index_of(j_1, 'z')]);
				if (td < min_v) {
					min_v = td;
				}
				if (td > max_v) {
					max_v = td;
				}
			}
		}

		// Метод линейного порядка.
		doublereal min_x = 1e60;
		doublereal min_y = 1e60;
		doublereal min_z = 1e60;
		doublereal max_x = -1e60;
		doublereal max_y = -1e60;
		doublereal max_z = -1e60;

		for (integer i = 0; i < t.maxnod; i++) {
			if (t.pa[i].x < min_x) {
				min_x = t.pa[i].x;
			}
			if (t.pa[i].y < min_y) {
				min_y = t.pa[i].y;
			}
			if (t.pa[i].z < min_z) {
				min_z = t.pa[i].z;
			}
			if (t.pa[i].x > max_x) {
				max_x = t.pa[i].x;
			}
			if (t.pa[i].y > max_y) {
				max_y = t.pa[i].y;
			}
			if (t.pa[i].z > max_z) {
				max_z = t.pa[i].z;
			}
		}

		//min_x *= 1.2;
		//min_y *= 1.2;
		//min_z *= 1.2;



		min_x = 1.05*fabs(max_x - min_x);
		if (min_x < 1.0e-30) {
			min_x = 1.05*fabs(max_x);
		}
		min_y = 1.05*fabs(max_y - min_y);
		if (min_y < 1.0e-30) {
			min_y = 1.05*fabs(max_y);
		}
		min_z = 1.05*fabs(max_z - min_z);
		if (min_z < 1.0e-30) {
			min_z = 1.05*fabs(max_z);
		}


		/*
		if (min_x < 1.0e-30) {
		printf("error!!! negative min_x MNK!\n");
		printf("min_x=%e max_x=%e\n",min_x,max_x);
		}
		if (min_y < 1.0e-30) {
		printf("error!!! negative min_y MNK!\n");
		printf("min_y=%e max_y=%e\n", min_y, max_y);
		}
		if (min_z < 1.0e-30) {
		printf("error!!! negative min_z MNK!\n");
		printf("min_z=%e max_z=%e\n", min_z, max_z);
		}
		*/

		TOCHKA** pointerlist = new TOCHKA*[t.maxelm];
		doublereal** rthdsd_Gauss = new doublereal*[t.maxelm];
		for (integer i_47 = 0; i_47 < t.maxelm; i_47++) {
			pointerlist[i_47] = new TOCHKA[8];
			rthdsd_Gauss[i_47] = new doublereal[8];
		}

		doublereal min_v1 = 1e60;
		doublereal max_v1 = -1e60;

		doublereal** Xmatr = new doublereal * [4];
		for (integer j = 0; j <= 3; j++) {
			Xmatr[j] = new doublereal[4];
		}


		doublereal* bmatr = new doublereal[4];
		doublereal* koefmatr = new doublereal[4];

		for (integer i = 0; i < t.maxelm; i++) {
			//doublereal xc47, yc47, zc47;

			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			//xc47 = p.x;
			//yc47 = p.y;
			//zc47 = p.z;


			p.x = p.x + min_x;
			p.y = p.y + min_y;
			p.z = p.z + min_z;

			for (integer j = 0; j <= 7; j++) {
				TOCHKA p1;
				p1.x = t.pa[t.nvtx[j][i] - 1].x;
				p1.y = t.pa[t.nvtx[j][i] - 1].y;
				p1.z = t.pa[t.nvtx[j][i] - 1].z;
				p1.x = p1.x + min_x;
				p1.y = p1.y + min_y;
				p1.z = p1.z + min_z;

				pointerlist[i][j] = p1;
				if (fabs(p1.x) < 1.0e-40) {
					printf("problem x=%e\n", p1.x);
					system("PAUSE");
				}
				if (fabs(p1.y) < 1.0e-40) {
					printf("problem y=%e\n", p1.y);
					system("PAUSE");
				}
				if (fabs(p1.z) < 1.0e-40) {
					printf("problem z=%e\n", p1.z);
					system("PAUSE");
				}
				integer j_1 = t.nvtx[j][i] - 1;
				switch (id_translate) {
				case TOTALDEFORMATION:// TOTAL DEFORMATION
					rthdsd_Gauss[i][j] = sqrt(input_vector_stress_3dim[index_of(j_1, 'x')] * input_vector_stress_3dim[index_of(j_1, 'x')]
						+ input_vector_stress_3dim[index_of(j_1, 'y')] * input_vector_stress_3dim[index_of(j_1, 'y')] +
						input_vector_stress_3dim[index_of(j_1, 'z')] * input_vector_stress_3dim[index_of(j_1, 'z')]);
					break;
				case XDEFORMATION: // X deformation
					rthdsd_Gauss[i][j] = input_vector_stress_3dim[index_of(j_1, 'x')]; // rthdsd[3 * j_1];
					break;
					// Поменял местами Y и Z 19.04.2019
				case YDEFORMATION: // Y deformation
					rthdsd_Gauss[i][j] = input_vector_stress_3dim[index_of(j_1, 'y')]; // rthdsd[3 * j_1 + 1];
					break;
				case ZDEFORMATION: // Z deformation
					rthdsd_Gauss[i][j] = input_vector_stress_3dim[index_of(j_1, 'z')]; // rthdsd[3 * j_1 + 2];
					break;
				default:
					printf("ERROR in Stress2Thermal_vector_translate in module mysolverv0_03.c!!!\n");
					printf("UNKNOWN id_translate variable == %lld\n", id_translate);
					system("PAUSE");
					exit(1);
					break;
				}

			}


			

			for (integer j1 = 0; j1 <= 3; j1++) {
				for (integer j2 = 0; j2 <= 3; j2++) {
					Xmatr[j1][j2] = 0.0;
				}
				bmatr[j1] = 0.0;
				koefmatr[j1] = 0.0;
			}




			for (integer j = 0; j < 8; j++) {

				Xmatr[0][0] += 1.0;
				Xmatr[0][1] += pointerlist[i][j].x;
				Xmatr[0][2] += pointerlist[i][j].y;
				Xmatr[0][3] += pointerlist[i][j].z;

				Xmatr[1][0] += pointerlist[i][j].x;
				Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
				Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
				Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;

				Xmatr[2][0] += pointerlist[i][j].y;
				Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
				Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
				Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;

				Xmatr[3][0] += pointerlist[i][j].z;
				Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
				Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
				Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;

				bmatr[0] += rthdsd_Gauss[i][j];
				bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
				bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
				bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];
			}


			for (integer j1 = 0; j1 <= 100; j1++) {
				koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
				koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
				koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
				koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
			}
			// Не забываем коэффициент 2.0
			//22,02,2019 Не нужен такой коэффициент. Проверено.
			output_vector_Thermal_dim[i] = 1.0*(koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));
			//if (j_6>0&&t.total_deformation[j_6][i] > 0) {
			//printf("%e\n", t.total_deformation[j_6][i]);
			//getchar();
			//}

			if (output_vector_Thermal_dim[i] < min_v1) {
				min_v1 = output_vector_Thermal_dim[i];
			}
			if (output_vector_Thermal_dim[i] > max_v1) {
				max_v1 = output_vector_Thermal_dim[i];
			}
			

		}

		for (integer j = 0; j <= 3; j++) {
			if (Xmatr[j] != nullptr) {
				delete[] Xmatr[j];
				Xmatr[j] = nullptr;
			}
		}
		delete[] Xmatr;
		Xmatr = nullptr;
		delete[] bmatr;
		bmatr = nullptr;
		delete[] koefmatr;
		koefmatr = nullptr;

		for (integer i = 0; i < t.maxelm; i++) {
			// Преобразование из вершин в центры ячеек с сохранением модуля величины.
			output_vector_Thermal_dim[i] *= ((max_v - min_v)/(max_v1-min_v1));
		}

		for (integer i = 0; i < t.maxelm; i++) {
			if (pointerlist[i] != nullptr) {
				delete[] pointerlist[i];
				pointerlist[i] = nullptr;
			}
			if (rthdsd_Gauss[i] != nullptr) {
				delete[] rthdsd_Gauss[i];
				rthdsd_Gauss[i] = nullptr;
			}
		}
		delete[] pointerlist;
		pointerlist = nullptr;
		delete[] rthdsd_Gauss;
		rthdsd_Gauss = nullptr;

	}

} //Stress2Thermal_vector_translate

void init_total_deformation(TEMPER &t) {
	
	// allocation memomory
	if (t.total_deformation == nullptr) {
		t.total_deformation = new doublereal*[SIZE_DEFORMATION_ARRAY];
		for (integer j_6 = 0; j_6 < SIZE_DEFORMATION_ARRAY; j_6++) {
			t.total_deformation[j_6] = nullptr;
			if (t.total_deformation[j_6] == nullptr) {
				t.total_deformation[j_6] = new doublereal[t.maxelm + t.maxbound];
			}
		}

	}
	else {
		for (integer j_6 = 0; j_6 < SIZE_DEFORMATION_ARRAY; j_6++) {
			delete[] t.total_deformation[j_6];
			t.total_deformation[j_6] = nullptr;
		}
		delete[] t.total_deformation;
		t.total_deformation = nullptr;

		t.total_deformation = new doublereal*[SIZE_DEFORMATION_ARRAY];
		for (integer j_6 = 0; j_6 < SIZE_DEFORMATION_ARRAY; j_6++) {
			t.total_deformation[j_6] = nullptr;
			if (t.total_deformation[j_6] == nullptr) {
				t.total_deformation[j_6] = new doublereal[t.maxelm + t.maxbound];
			}
		}
	}
	// init zero (0.0)
	for (integer j_6 = 0; j_6 < SIZE_DEFORMATION_ARRAY; j_6++) {
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {
			t.total_deformation[j_6][i_1] = 0.0;
		}
	}
} // init_total_deformation

void Stress2Thermal_vector_translate(TEMPER &t,
	doublereal* &input_vector_stress_1dim,
	doublereal* &output_vector_Thermal_dim) {

	// t - вспомогательные данные: координаты узлов и пр. (pa, nvtx, maxelm, maxnod).
	// id_translate: 0-total, 1 - x def, 2 - y def, 3 - z def.


	if (1) {

		// Метод линейного порядка.
		doublereal min_x = 1e60;
		doublereal min_y = 1e60;
		doublereal min_z = 1e60;
		doublereal max_x = -1e60;
		doublereal max_y = -1e60;
		doublereal max_z = -1e60;

		for (integer i = 0; i < t.maxnod; i++) {
			if (t.pa[i].x < min_x) {
				min_x = t.pa[i].x;
			}
			if (t.pa[i].y < min_y) {
				min_y = t.pa[i].y;
			}
			if (t.pa[i].z < min_z) {
				min_z = t.pa[i].z;
			}
			if (t.pa[i].x > max_x) {
				max_x = t.pa[i].x;
			}
			if (t.pa[i].y > max_y) {
				max_y = t.pa[i].y;
			}
			if (t.pa[i].z > max_z) {
				max_z = t.pa[i].z;
			}
		}

		//min_x *= 1.2;
		//min_y *= 1.2;
		//min_z *= 1.2;



		min_x = 1.05*fabs(max_x - min_x);
		if (min_x < 1.0e-30) {
			min_x = 1.05*fabs(max_x);
		}
		min_y = 1.05*fabs(max_y - min_y);
		if (min_y < 1.0e-30) {
			min_y = 1.05*fabs(max_y);
		}
		min_z = 1.05*fabs(max_z - min_z);
		if (min_z < 1.0e-30) {
			min_z = 1.05*fabs(max_z);
		}


		/*
		if (min_x < 1.0e-30) {
		printf("error!!! negative min_x MNK!\n");
		printf("min_x=%e max_x=%e\n",min_x,max_x);
		}
		if (min_y < 1.0e-30) {
		printf("error!!! negative min_y MNK!\n");
		printf("min_y=%e max_y=%e\n", min_y, max_y);
		}
		if (min_z < 1.0e-30) {
		printf("error!!! negative min_z MNK!\n");
		printf("min_z=%e max_z=%e\n", min_z, max_z);
		}
		*/

		TOCHKA** pointerlist = new TOCHKA*[t.maxelm];
		doublereal** rthdsd_Gauss = new doublereal*[t.maxelm];
		for (integer i_47 = 0; i_47 < t.maxelm; i_47++) {
			pointerlist[i_47] = new TOCHKA[8];
			rthdsd_Gauss[i_47] = new doublereal[8];
		}


		doublereal** Xmatr = new doublereal * [4];
		for (integer j = 0; j <= 3; j++) {
			Xmatr[j] = new doublereal[4];
		}


		doublereal* bmatr = new doublereal[4];
		doublereal* koefmatr = new doublereal[4];

		for (integer i = 0; i < t.maxelm; i++) {
			//doublereal xc47, yc47, zc47;

			TOCHKA p;
			center_cord3D(i, t.nvtx, t.pa, p, 100);
			//xc47 = p.x;
			//yc47 = p.y;
			//zc47 = p.z;


			p.x = p.x + min_x;
			p.y = p.y + min_y;
			p.z = p.z + min_z;

			for (integer j = 0; j <= 7; j++) {
				TOCHKA p1;
				p1.x = t.pa[t.nvtx[j][i] - 1].x;
				p1.y = t.pa[t.nvtx[j][i] - 1].y;
				p1.z = t.pa[t.nvtx[j][i] - 1].z;
				p1.x = p1.x + min_x;
				p1.y = p1.y + min_y;
				p1.z = p1.z + min_z;

				pointerlist[i][j] = p1;
				if (fabs(p1.x) < 1.0e-40) {
					printf("problem x=%e\n", p1.x);
					system("PAUSE");
				}
				if (fabs(p1.y) < 1.0e-40) {
					printf("problem y=%e\n", p1.y);
					system("PAUSE");
				}
				if (fabs(p1.z) < 1.0e-40) {
					printf("problem z=%e\n", p1.z);
					system("PAUSE");
				}
				integer j_1 = t.nvtx[j][i] - 1;
				rthdsd_Gauss[i][j] = input_vector_stress_1dim[j_1]; // rthdsd[j_1];
				
			}


			

			for (integer j1 = 0; j1 <= 3; j1++) {
				for (integer j2 = 0; j2 <= 3; j2++) {
					Xmatr[j1][j2] = 0.0;
				}
				bmatr[j1] = 0.0;
				koefmatr[j1] = 0.0;
			}




			for (integer j = 0; j < 8; j++) {

				Xmatr[0][0] += 1.0;
				Xmatr[0][1] += pointerlist[i][j].x;
				Xmatr[0][2] += pointerlist[i][j].y;
				Xmatr[0][3] += pointerlist[i][j].z;

				Xmatr[1][0] += pointerlist[i][j].x;
				Xmatr[1][1] += pointerlist[i][j].x*pointerlist[i][j].x;
				Xmatr[1][2] += pointerlist[i][j].x*pointerlist[i][j].y;
				Xmatr[1][3] += pointerlist[i][j].x*pointerlist[i][j].z;

				Xmatr[2][0] += pointerlist[i][j].y;
				Xmatr[2][1] += pointerlist[i][j].y*pointerlist[i][j].x;
				Xmatr[2][2] += pointerlist[i][j].y*pointerlist[i][j].y;
				Xmatr[2][3] += pointerlist[i][j].y*pointerlist[i][j].z;

				Xmatr[3][0] += pointerlist[i][j].z;
				Xmatr[3][1] += pointerlist[i][j].z*pointerlist[i][j].x;
				Xmatr[3][2] += pointerlist[i][j].z*pointerlist[i][j].y;
				Xmatr[3][3] += pointerlist[i][j].z*pointerlist[i][j].z;

				bmatr[0] += rthdsd_Gauss[i][j];
				bmatr[1] += pointerlist[i][j].x*rthdsd_Gauss[i][j];
				bmatr[2] += pointerlist[i][j].y*rthdsd_Gauss[i][j];
				bmatr[3] += pointerlist[i][j].z*rthdsd_Gauss[i][j];
			}

			// 22.02.2019 Достаточно 100 итераций.
			for (integer j1 = 0; j1 <= 100; j1++) {
				koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
				koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
				koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
				koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
			}
			// Не забываем коэффициент 2.0
			//22,02,2019 Не нужен такой коэффициент. Проверено.
			output_vector_Thermal_dim[i] = 1.0*(koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));
			//if (j_6>0&&t.total_deformation[j_6][i] > 0) {
			//printf("%e\n", t.total_deformation[j_6][i]);
			//getchar();
			//}			

		}

		for (integer j = 0; j <= 3; j++) {
			delete[] Xmatr[j];
		}
		delete[] Xmatr;
		delete[] bmatr;
		delete[] koefmatr;

		for (integer i = 0; i < t.maxelm; i++) {
			delete[] pointerlist[i];
			delete[] rthdsd_Gauss[i];
		}
		delete[] pointerlist;
		delete[] rthdsd_Gauss;

	}

} //Stress2Thermal_vector_translate

// Возвращает квадрат числа.
doublereal  sqr(doublereal x) {
	return x * x;
}

doublereal epsilon(doublereal current, doublereal min, doublereal max, doublereal center) {
	// Для фиксированного координатного напрпвления.
	// min,max - границы объекта по оси Ох.
	// center - координата центра масс.
	// current - координата текущей точки.
	if (current - center > 0.0) {
		return ((current - center) / (max-center));
	}
	else {
		return ((current - center) / (center-min));
	}
} // для силы линейного теплового расширения.


  


// Решение прочностной задачи в 3D.
// 6 августа 2017. август 2020.
void solve_Structural(TEMPER &t, WALL* &w, integer lw, 
	bool bThermalStress, doublereal operatingtemperature, 
	BLOCK* &b, integer &lb, integer &lu,
	bool btimedep, doublereal timestep_sizenow,
	doublereal* &uoldtimestep, doublereal* &uolddoubletimestep,
	doublereal poweron_multiplyer_sequence, TPROP* &matlist,
	doublereal* &t_for_Mechanical) {

	// btimedep==true - нестационарное моделирование,
	// btimedep==false - стационарная задача механики.
	// timestep_sizenow - размер шага по времени,
	// uoldtimestep - перемещения на один шаг назад,
	// uolddoubletimestep - перемещения на два шага назад.
	// uoldtimestep, uolddoubletimestep - память выделена заранее в вызывающем внешнем коде.
	// poweron_multiplyer_sequence==0.0 вектор силы выключен,
	// poweron_multiplyer_sequence==1.0 вектор силы полностью включён,
	// 0 < poweron_multiplyer_sequence < 1 - вектор силы частично активен.


	printf("Stress n=%lld\n", 3 * t.maxnod);

	doublereal* rthdsd = new doublereal[3*t.maxnod]; // Правая часть.
	doublereal* deformation = new doublereal[3*t.maxnod]; // Деформация.
	bool* constr = new bool[3 * t.maxnod]; // Фиксированное смещение.
	CylindricalSupport* cylsup = new CylindricalSupport[3 * t.maxnod];

	// инициализация.
	for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1++) {
		rthdsd[i_1] = 0.0;
		deformation[i_1] = 0.0;
		constr[i_1] = false; // По умолчанию все узлы свободны.
		cylsup[i_1].bactive = false; // Не используется.
	}

	// Вычисление допусков.
	doublereal epsx = 1.0e+30, epsy = 1.0e+30, epsz = 1.0e+30;
	for (integer ie = 0; ie < t.maxelm; ie++) {
		doublereal hx = 0.0, hy = 0.0, hz = 0.0;
		volume3D(ie, t.nvtx, t.pa, hx, hy, hz);
		if (0.3*hx < epsx) epsx = 0.3*hx;
		if (0.3*hy < epsy) epsy = 0.3*hy;
		if (0.3*hz < epsz) epsz = 0.3*hz;
	}

	for (integer iu_74 = 0; iu_74 < lu; iu_74++) {
		for (integer ie = 0; ie < my_union[iu_74].t.maxelm; ie++) {
			doublereal hx = 0.0, hy = 0.0, hz = 0.0;
			volume3D(ie, my_union[iu_74].t.nvtx, my_union[iu_74].t.pa, hx, hy, hz);
			if (0.3*hx < epsx) epsx = 0.3*hx;
			if (0.3*hy < epsy) epsy = 0.3*hy;
			if (0.3*hz < epsz) epsz = 0.3*hz;
		}
	}
	// Для асемблесов не работает.


	// В каждом узле хранит список nvtx которые имеют вершиной этот узел.
	integer *nvtx_link_count = new integer[t.maxnod];
	integer **nvtx_link = new integer*[t.maxnod];
	for (integer i_4 = 0; i_4 < t.maxnod; i_4++) {
		nvtx_link_count[i_4] = 0;
		nvtx_link[i_4] = new integer[8];
		for (integer i_5 = 0; i_5 < 8; i_5++) nvtx_link[i_4][i_5] = -1;
	}
	for (integer j_1 = 0; j_1 < t.maxelm; j_1++) {
		for (integer j_75 = 0; j_75 < 8; j_75++) {
			nvtx_link[t.nvtx[j_75][j_1] - 1][nvtx_link_count[t.nvtx[j_75][j_1] - 1]] = j_1;
			nvtx_link_count[t.nvtx[j_75][j_1] - 1]++;
		}
	}

	doublereal hxl = 0.0, hyl = 0.0, hzl = 0.0, dSl = 0.0, dln = 0.0, vol_l = 0.0;
	integer k_1l = 0;

	//integer ie = 0;
	//for (integer j = 0; j < 8; j++) {
	//	printf("%e %e %e\n", t.pa[t.nvtx[j][ie] - 1].x, t.pa[t.nvtx[j][ie] - 1].y, t.pa[t.nvtx[j][ie] - 1].z);
	//}
    //	getchar();
	// Учёт граничных условий.
	// На зафиксированных участках границы мы выставляем флаг true.
	for (integer i_1 = 0; i_1 < lw; i_1++) {
		
			//const doublereal eps1 = 1.0e-30;
			// pa нумеруется  нуля.
			for (integer j_1 = 0; j_1 < t.maxnod; j_1++) {
				bool bfound = false;
				switch (w[i_1].iPlane) {
				case XY_PLANE: 
					if ((fabs(t.pa[j_1].z-w[i_1].g.zS) < epsz)&&
					         (t.pa[j_1].x<w[i_1].g.xE + epsx)&&
					         (t.pa[j_1].x>w[i_1].g.xS - epsx)&&
					         (t.pa[j_1].y>w[i_1].g.yS - epsy)&&
					         (t.pa[j_1].y<w[i_1].g.yE + epsy)) 
				{
					//printf("found: plane XY wall[%lld]\n",i_1);
					bfound = true;
				}
					break;
				case YZ_PLANE:
					if ((fabs(t.pa[j_1].x - w[i_1].g.xS)<epsx)
						&& (t.pa[j_1].z<w[i_1].g.zE + epsz)
						&& (t.pa[j_1].z>w[i_1].g.zS - epsz) 
						&& (t.pa[j_1].y>w[i_1].g.yS - epsy)
						&& (t.pa[j_1].y<w[i_1].g.yE + epsy)) 
					{
						//printf("found: plane YZ wall[%lld]\n", i_1);
						bfound = true;
					}
					break;
				case XZ_PLANE:
					if ((fabs(t.pa[j_1].y - w[i_1].g.yS)<epsy)
						&& (t.pa[j_1].z<w[i_1].g.zE + epsz) 
						&& (t.pa[j_1].z>w[i_1].g.zS - epsz)
						&& (t.pa[j_1].x>w[i_1].g.xS - epsx)
						&& (t.pa[j_1].x<w[i_1].g.xE + epsx))
					{
						//printf("found: plane XZ wall[%lld]\n", i_1);
						bfound = true;
					}
					break;
				}
				if (bfound) {
					// Фиксированный потенциал.
					// Thermal-Stress boundary condition
					// 0 - free,
					// 1 - x fixit,
					// 2 - y fixit,
					// 3 - z fixit,
					// 4 - xy fixit,
					// 5 - xz fixit,
					// 6 - yz fixit,
					// 7 - fixit all,
					// 8 - x Force,
					// 9 - y Force,
					// 10 - z Force.
					switch (w[i_1].ithermal_Stress_boundary_condition) {
					case THERMAL_STRESS_BOUNDARY_CONDITION::FREE: //FREE all
						// Граница не фиксирована.
						break;
					case THERMAL_STRESS_BOUNDARY_CONDITION::X_FIXIT:  constr[index_of(j_1, 'x')] = true; // X
						break;
					case THERMAL_STRESS_BOUNDARY_CONDITION::Y_FIXIT:  constr[index_of(j_1, 'y')] = true; // Y
						break;
					case THERMAL_STRESS_BOUNDARY_CONDITION::Z_FIXIT:  constr[index_of(j_1, 'z')] = true; // Z
						break;
					case THERMAL_STRESS_BOUNDARY_CONDITION::XY_FIXIT: 
							constr[index_of(j_1, 'x')] = true; // X
						     constr[index_of(j_1, 'y')] = true; // Y
						break;
					case THERMAL_STRESS_BOUNDARY_CONDITION::XZ_FIXIT:  
							constr[index_of(j_1, 'x')] = true; //X
						     constr[index_of(j_1, 'z')] = true; //Z
						break;
					case THERMAL_STRESS_BOUNDARY_CONDITION::YZ_FIXIT:
							constr[index_of(j_1, 'y')] = true; // Y 
						    constr[index_of(j_1, 'z')] = true; // Z
						break;
					case THERMAL_STRESS_BOUNDARY_CONDITION::ALL_FIXIT:
						//printf("ok");
						//getchar();//ok
						// Граница полностью фиксирована
						// по всем трём координатам.
						constr[index_of(j_1, 'x')] = true;//X
						constr[index_of(j_1, 'y')] = true;//Y
						constr[index_of(j_1, 'z')] = true;//Z
						//printf("ALL Fixit Ok\n");
						break;
					case  THERMAL_STRESS_BOUNDARY_CONDITION::X_FORCE:
					case  THERMAL_STRESS_BOUNDARY_CONDITION::Y_FORCE:
					case  THERMAL_STRESS_BOUNDARY_CONDITION::Z_FORCE:
						// Здесь обязательно нужно умножить на площадь.
						// По размерности в правой части стоит именно сила в Ньютонах.
						// Площадь нужна если задано давление.
						switch (w[i_1].iPlane) {
						case XY_PLANE:
							dSl = 0.0;
							for (k_1l = 0; k_1l < 8; k_1l++) {
								if (nvtx_link[j_1][k_1l] > -1) {
									volume3D(nvtx_link[j_1][k_1l], t.nvtx, t.pa, hxl, hyl, hzl);
									dSl += 0.25*hyl*hxl;
								}
							}
							dSl = 1.0;
							if (btimedep) {
								rthdsd[index_of(j_1, 'z')] = poweron_multiplyer_sequence*w[i_1].zForce; // Normal component.
							}
							else {
								rthdsd[index_of(j_1, 'z')] = w[i_1].zForce; // Normal component.
							}
							//rthdsd[index_of(j_1, 'z')] *= hzl * hzl;
							// На границе где приложена нормальная сила 
							//разрешаем лишь нормальные деформации.
							//constr[3 * j_1] = true;//X
							//constr[3 * j_1 + 1] = true;//Y
							//constr[3 * j_1 + 2] = true;//Z

							dSl = 0;
							for (k_1l = 0; k_1l < 8; k_1l++) {
								if (nvtx_link[j_1][k_1l] > -1) {
									volume3D(nvtx_link[j_1][k_1l], t.nvtx, t.pa, hxl, hyl, hzl);
									dSl += 0.25*hyl*hzl;
								}
							}
							if (nvtx_link_count[j_1] == 4) dSl *= 0.5;
							dSl = 1.0;
							if (btimedep) {
								rthdsd[index_of(j_1, 'x')] = poweron_multiplyer_sequence*w[i_1].xForce;
							}
							else {
								rthdsd[index_of(j_1, 'x')] = w[i_1].xForce;
							}
							//rthdsd[index_of(j_1, 'x')] *= hxl * hxl;
							dSl = 0;
							for (k_1l = 0; k_1l < 8; k_1l++) {
								if (nvtx_link[j_1][k_1l] > -1) {
									volume3D(nvtx_link[j_1][k_1l], t.nvtx, t.pa, hxl, hyl, hzl);
									dSl += 0.25*hxl*hzl;
								}
							}
							if (nvtx_link_count[j_1] == 4) dSl *= 0.5;
							dSl = 1.0;
							if (btimedep) {
								rthdsd[index_of(j_1, 'y')] = poweron_multiplyer_sequence*w[i_1].yForce;
							}
							else {
								rthdsd[index_of(j_1, 'y')] = w[i_1].yForce;
							}
							//rthdsd[index_of(j_1, 'y')] *= hyl * hyl;
							break;
						case YZ_PLANE:
							dSl = 0.0; dln = 0.0; vol_l = 0.0;
							for (k_1l = 0; k_1l < 8; k_1l++) {
								if (nvtx_link[j_1][k_1l] > -1) {
									volume3D(nvtx_link[j_1][k_1l], t.nvtx, t.pa, hxl, hyl , hzl);
									dSl += 0.25*hyl*hzl;
									dln = 0.5*hxl;
									vol_l += 0.125*hxl*hyl*hzl;
								}
							}
							//dSl = 1.0;
							if (btimedep) {
								rthdsd[index_of(j_1, 'x')] = poweron_multiplyer_sequence*w[i_1].xForce;// Normal component.
							}
							else {
								rthdsd[index_of(j_1, 'x')] = w[i_1].xForce;// Normal component.
							}
							//rthdsd[index_of(j_1, 'x')] *= hxl * hxl;
							// На границе где приложена нормальная сила 
							//разрешаем лишь нормальные деформации.
							//constr[3 * j_1] = true;//X
								//constr[3 * j_1 + 1] = true;//Y
								//constr[3 * j_1 + 2] = true;//Z

							dSl = 0;
							for (k_1l = 0; k_1l < 8; k_1l++) {
								if (nvtx_link[j_1][k_1l] > -1) {
									volume3D(nvtx_link[j_1][k_1l], t.nvtx, t.pa, hxl, hyl, hzl);
									dSl += 0.25*hxl*hzl;
								}
							}
							if (nvtx_link_count[j_1] == 4) dSl *= 0.5; dSl = 1.0;
							if (btimedep) {
								rthdsd[index_of(j_1, 'y')] = poweron_multiplyer_sequence*w[i_1].yForce;
							}
							else {
								rthdsd[index_of(j_1, 'y')] = w[i_1].yForce;
							}
							//rthdsd[index_of(j_1, 'y')] *= hyl * hyl;
							dSl = 0;
							for (k_1l = 0; k_1l < 8; k_1l++) {
								if (nvtx_link[j_1][k_1l] > -1) {
									volume3D(nvtx_link[j_1][k_1l], t.nvtx, t.pa, hxl, hyl, hzl);
									dSl += 0.25*hxl*hyl;
								}
							}
							if (nvtx_link_count[j_1] == 4) dSl *= 0.5; dSl = 1.0;
							if (btimedep) {
								rthdsd[index_of(j_1, 'z')] = poweron_multiplyer_sequence*w[i_1].zForce;
							}
							else {
								rthdsd[index_of(j_1, 'z')] = w[i_1].zForce;
							}
							//rthdsd[index_of(j_1, 'z')] *= hzl * hzl;
							break;
						case XZ_PLANE:
							dSl = 0.0;
							for (k_1l = 0; k_1l < 8; k_1l++) {
								if (nvtx_link[j_1][k_1l] > -1) {
									volume3D(nvtx_link[j_1][k_1l], t.nvtx, t.pa, hxl, hyl, hzl);
									dSl += 0.25*hxl*hzl;
								}
							}
							dSl = 1.0;
							if (btimedep) {
								rthdsd[index_of(j_1, 'y')] = poweron_multiplyer_sequence*w[i_1].yForce;
							}
							else {
								rthdsd[index_of(j_1, 'y')] = w[i_1].yForce;// Normal component.
							}
							//rthdsd[index_of(j_1, 'y')] *= hyl * hyl;
							//printf("w[i_1].yForce=%e\n", w[i_1].yForce);
							//getchar();
							// На границе где приложена нормальная сила 
							//разрешаем лишь нормальные деформации.
								//constr[3 * j_1] = true;//X
							//constr[3 * j_1 + 1] = true;//Y
								//constr[3 * j_1 + 2] = true;//Z

							dSl = 0;
							for (k_1l = 0; k_1l < 8; k_1l++) {
								if (nvtx_link[j_1][k_1l] > -1) {
									volume3D(nvtx_link[j_1][k_1l], t.nvtx, t.pa, hxl, hyl, hzl);
									dSl += 0.25*hyl*hzl;
								}
							}
							if (nvtx_link_count[j_1] == 4) dSl *= 0.5; dSl = 1.0;
							if (btimedep) {
								rthdsd[index_of(j_1, 'x')] = poweron_multiplyer_sequence*w[i_1].xForce;
							}
							else {
								rthdsd[index_of(j_1, 'x')] = w[i_1].xForce;
							}
							///rthdsd[index_of(j_1, 'x')] *= hxl * hxl;
							dSl = 0;
							for (k_1l = 0; k_1l < 8; k_1l++) {
								if (nvtx_link[j_1][k_1l] > -1) {
									volume3D(nvtx_link[j_1][k_1l], t.nvtx, t.pa, hxl, hyl, hzl);
									dSl += 0.25*hxl*hyl;
								}
							}
							if (nvtx_link_count[j_1] == 4) dSl *= 0.5; dSl = 1.0;
							if (btimedep) {
								rthdsd[index_of(j_1, 'z')] = poweron_multiplyer_sequence*w[i_1].zForce;
							}
							else {
								rthdsd[index_of(j_1, 'z')] = w[i_1].zForce;
							}
							//rthdsd[index_of(j_1, 'z')] *= hzl * hzl;
							break;
						}
						//rthdsd[3 * j_1] = w[i_1].xForce;
						//rthdsd[3 * j_1 + 1] = w[i_1].yForce;
						//rthdsd[3 * j_1 + 2] = w[i_1].zForce;
						//printf("Fotce X =%e %e %e\n", w[i_1].xForce, w[i_1].yForce, w[i_1].zForce);
						//getchar();//ok
						break;
					}					
					
				}
			}		
	}

	// Учёт граничных условий закрепления.
	for (integer j_11 = 0; j_11 < t.maxelm; j_11++) {
		for (integer i_11 = 0; i_11 < 8; i_11++) {
			//for (integer j_1 = 0; j_1 < t.maxnod; j_1++) {
			integer j_1 = t.nvtx[i_11][j_11] - 1;
			/*
			if (constr[3 * j_1+1] && fabs(rthdsd[3 * j_1]) > 0.0) {
				for (integer i_111 = 0; i_111 < 8; i_111++) {
					integer j_111 = t.nvtx[i_111][j_11] - 1;
					constr[3 * j_111] = false;
					rthdsd[3 * j_111] = 0.0;
					constr[3 * j_111 + 2] = 0.0;//Z -fix
				}
			}
			*/
			if (constr[index_of(j_1, 'x')]) rthdsd[index_of(j_1, 'x')] = 0.0; // X
			if (constr[index_of(j_1, 'y')]) rthdsd[index_of(j_1, 'y')] = 0.0; // Y 
			if (constr[index_of(j_1, 'z')]) rthdsd[index_of(j_1, 'z')] = 0.0; // Z 
		}
	}

	delete[] nvtx_link_count;
	nvtx_link_count = nullptr;
	for (integer i_4 = 0; i_4 < t.maxnod; i_4++) {
		delete[] nvtx_link[i_4];
	}
	delete[] nvtx_link;
	nvtx_link = nullptr;


	if (0) {
		// Включить если нужно обязательно зафиксировать смещения на боковых стенках цилиндров.

		for (integer i_1 = 0; i_1 < lb; i_1++) {
			// Фиксация на боковых стенках цилиндров.
			if (b[i_1].g.itypegeom == CYLINDER) {
				// определяет принадлежность точки Цилиндру.
				for (integer j_1 = 0; j_1 < t.maxnod; j_1++) {

					switch (b[i_1].g.iPlane) {
					case XY_PLANE: if (sqrt((t.pa[j_1].x - b[i_1].g.xC)*(t.pa[j_1].x - b[i_1].g.xC) + (t.pa[j_1].y - b[i_1].g.yC)*(t.pa[j_1].y - b[i_1].g.yC)) < b[i_1].g.R_out_cyl + sqrt(9 * epsx*epsx + 9 * epsy*epsy)) {
						if ((t.pa[j_1].z > b[i_1].g.zC - 3 * epsz) && (t.pa[j_1].z < b[i_1].g.zC + b[i_1].g.Hcyl + 3 * epsz)) {
							constr[3 * j_1] = true;//X
							constr[3 * j_1 + 1] = true;//Y
							constr[3 * j_1 + 2] = true;//Z
						}
					}
							 break;
					case XZ_PLANE:
						if (sqrt((t.pa[j_1].x - b[i_1].g.xC)*(t.pa[j_1].x - b[i_1].g.xC) + (t.pa[j_1].z - b[i_1].g.zC)*(t.pa[j_1].z - b[i_1].g.zC)) < b[i_1].g.R_out_cyl + sqrt(9 * epsx*epsx + 9 * epsz*epsz)) {
							if ((t.pa[j_1].y > b[i_1].g.yC - 3 * epsy) && (t.pa[j_1].y < b[i_1].g.yC + b[i_1].g.Hcyl + 3 * epsy)) {
								constr[3 * j_1] = true;//X
								constr[3 * j_1 + 1] = true;//Y
								constr[3 * j_1 + 2] = true;//Z
							}
						}
						break;
					case YZ_PLANE:
						if (sqrt((t.pa[j_1].y - b[i_1].g.yC)*(t.pa[j_1].y - b[i_1].g.yC) + (t.pa[j_1].z - b[i_1].g.zC)*(t.pa[j_1].z - b[i_1].g.zC)) < b[i_1].g.R_out_cyl + sqrt(9 * epsz*epsz + 9 * epsy*epsy)) {
							if ((t.pa[j_1].x > b[i_1].g.xC - 3 * epsx) && (t.pa[j_1].x < b[i_1].g.xC + b[i_1].g.Hcyl + 3 * epsx)) {
								constr[3 * j_1] = true;//X
								constr[3 * j_1 + 1] = true;//Y
								constr[3 * j_1 + 2] = true;//Z
							}
						}
						break;
					}

				}
			}
		}
	}

	if (1) {
		// Включить если нужно обязательно зафиксировать смещения на боковых стенках цилиндров.
		// Более точный допуск eps.

		for (integer i_1 = 0; i_1 < lb; i_1++) {
			// Фиксация на боковых стенках цилиндров.
			// b[i_1].CylinderFixed - только если пользователь отметил в интерфейсе галочку для i_1 блока.
			if ((b[i_1].g.itypegeom == CYLINDER)&&(b[i_1].CylinderFixed)) {
				for (integer k_1 = 0; k_1 < t.maxelm; k_1++) {
					doublereal hx = 0.0, hy = 0.0, hz = 0.0;
					volume3D(k_1, t.nvtx, t.pa, hx, hy, hz);

					// определяет принадлежность точки Цилиндру.
					for (integer j_11 = 0; j_11 < 8; j_11++) {
						integer j_1 = t.nvtx[j_11][k_1] - 1;

						doublereal epsx1=0.3*hx;
						doublereal epsy1=0.3*hy;
						doublereal epsz1=0.3*hz;

						switch (b[i_1].g.iPlane) {
						case XY_PLANE: if (sqrt((t.pa[j_1].x - b[i_1].g.xC)*(t.pa[j_1].x - b[i_1].g.xC) + (t.pa[j_1].y - b[i_1].g.yC)*(t.pa[j_1].y - b[i_1].g.yC)) < b[i_1].g.R_out_cyl + sqrt(epsx1*epsx1 + epsy1*epsy1)) {
							if ((t.pa[j_1].z > b[i_1].g.zC - 3 * epsz1) && (t.pa[j_1].z < b[i_1].g.zC + b[i_1].g.Hcyl + 3 * epsz1)) {
								constr[index_of(j_1, 'x')] = true; //X
								constr[index_of(j_1, 'y')] = true; //Y
								constr[index_of(j_1, 'z')] = true; //Z
							}
						}
								 break;
						case XZ_PLANE:
							if (sqrt((t.pa[j_1].x - b[i_1].g.xC)*(t.pa[j_1].x - b[i_1].g.xC) + (t.pa[j_1].z - b[i_1].g.zC)*(t.pa[j_1].z - b[i_1].g.zC)) < b[i_1].g.R_out_cyl + sqrt(epsx1*epsx1 + epsz1*epsz1)) {
								if ((t.pa[j_1].y > b[i_1].g.yC - 3 * epsy1) && (t.pa[j_1].y < b[i_1].g.yC + b[i_1].g.Hcyl + 3 * epsy1)) {
									constr[index_of(j_1, 'x')] = true; //X
									constr[index_of(j_1, 'y')] = true; //Y
									constr[index_of(j_1, 'z')] = true; //Z
								}
							}
							break;
						case YZ_PLANE:
							if (sqrt((t.pa[j_1].y - b[i_1].g.yC)*(t.pa[j_1].y - b[i_1].g.yC) + (t.pa[j_1].z - b[i_1].g.zC)*(t.pa[j_1].z - b[i_1].g.zC)) < b[i_1].g.R_out_cyl + sqrt(epsz1*epsz1 + epsy1*epsy1)) {
								if ((t.pa[j_1].x > b[i_1].g.xC - 3 * epsx1) && (t.pa[j_1].x < b[i_1].g.xC + b[i_1].g.Hcyl + 3 * epsx1)) {
									constr[index_of(j_1, 'x')] = true; //X
									constr[index_of(j_1, 'y')] = true; //Y
									constr[index_of(j_1, 'z')] = true; //Z
								}
							}
							break;
						}

					}
				}
			}
		}
	}

	if (0) {
		// Включить если нужно обязательно разрешить вращение вокруг боковых стенок цилиндров.

		for (integer i_1 = 0; i_1 < lb; i_1++) {
			// Фиксация на боковых стенках цилиндров.
			if (b[i_1].g.itypegeom == CYLINDER) {
				for (integer k_1 = 0; k_1 < t.maxelm; k_1++) {
					doublereal hx = 0.0, hy = 0.0, hz = 0.0;
					volume3D(k_1, t.nvtx, t.pa, hx, hy, hz);

					// определяет принадлежность точки Цилиндру.
					for (integer j_11 = 0; j_11 < 8; j_11++) {
						integer j_1 = t.nvtx[j_11][k_1] - 1;

						doublereal epsx1 = 0.3*hx;
						doublereal epsy1 = 0.3*hy;
						doublereal epsz1 = 0.3*hz;

						switch (b[i_1].g.iPlane) {
						case XY_PLANE: if (sqrt((t.pa[j_1].x - b[i_1].g.xC)*(t.pa[j_1].x - b[i_1].g.xC) + (t.pa[j_1].y - b[i_1].g.yC)*(t.pa[j_1].y - b[i_1].g.yC)) < b[i_1].g.R_out_cyl + sqrt(epsx1*epsx1 + epsy1*epsy1)) {
							if ((t.pa[j_1].z > b[i_1].g.zC - 3 * epsz1) && (t.pa[j_1].z < b[i_1].g.zC + b[i_1].g.Hcyl + 3 * epsz1)) {
								constr[3 * j_1] = true;//X
								constr[3 * j_1 + 1] = true;//Y
								cylsup[3 * j_1].bactive = true;
								cylsup[3 * j_1].iPlane = b[i_1].g.iPlane;
								cylsup[3 * j_1].Radius = b[i_1].g.R_out_cyl;
								cylsup[3 * j_1].xC = b[i_1].g.xC;
								cylsup[3 * j_1].yC = b[i_1].g.yC;
								cylsup[3 * j_1].zC = b[i_1].g.zC;
								cylsup[3 * j_1].x1 = t.pa[j_1].x;
								cylsup[3 * j_1].y1 = t.pa[j_1].y;
								cylsup[3 * j_1].z1 = t.pa[j_1].z;
								cylsup[3 * j_1 + 1].bactive = true;
								cylsup[3 * j_1 + 1].iPlane = b[i_1].g.iPlane;
								cylsup[3 * j_1 + 1].Radius = b[i_1].g.R_out_cyl;
								cylsup[3 * j_1 + 1].xC = b[i_1].g.xC;
								cylsup[3 * j_1 + 1].yC = b[i_1].g.yC;
								cylsup[3 * j_1 + 1].zC = b[i_1].g.zC;
								cylsup[3 * j_1 + 1].x1 = t.pa[j_1].x;
								cylsup[3 * j_1 + 1].y1 = t.pa[j_1].y;
								cylsup[3 * j_1 + 1].z1 = t.pa[j_1].z;
								constr[3 * j_1 + 2] = true;//Z
							}
						}
								 break;
						case XZ_PLANE:
							if (sqrt((t.pa[j_1].x - b[i_1].g.xC)*(t.pa[j_1].x - b[i_1].g.xC) + (t.pa[j_1].z - b[i_1].g.zC)*(t.pa[j_1].z - b[i_1].g.zC)) < b[i_1].g.R_out_cyl + sqrt(epsx1*epsx1 + epsz1*epsz1)) {
								if ((t.pa[j_1].y > b[i_1].g.yC - 3 * epsy1) && (t.pa[j_1].y < b[i_1].g.yC + b[i_1].g.Hcyl + 3 * epsy1)) {
									constr[3 * j_1] = true;//X
									cylsup[3 * j_1].bactive = true;
									cylsup[3 * j_1].iPlane = b[i_1].g.iPlane;
									cylsup[3 * j_1].Radius = b[i_1].g.R_out_cyl;
									cylsup[3 * j_1].xC = b[i_1].g.xC;
									cylsup[3 * j_1].yC = b[i_1].g.yC;
									cylsup[3 * j_1].zC = b[i_1].g.zC;
									cylsup[3 * j_1].x1 = t.pa[j_1].x;
									cylsup[3 * j_1].y1 = t.pa[j_1].y;
									cylsup[3 * j_1].z1 = t.pa[j_1].z;
									constr[3 * j_1 + 1] = true;//Y
									constr[3 * j_1 + 2] = true;//Z
									cylsup[3 * j_1 + 2].bactive = true;
									cylsup[3 * j_1 + 2].iPlane = b[i_1].g.iPlane;
									cylsup[3 * j_1 + 2].Radius = b[i_1].g.R_out_cyl;
									cylsup[3 * j_1 + 2].xC = b[i_1].g.xC;
									cylsup[3 * j_1 + 2].yC = b[i_1].g.yC;
									cylsup[3 * j_1 + 2].zC = b[i_1].g.zC;
									cylsup[3 * j_1 + 2].x1 = t.pa[j_1].x;
									cylsup[3 * j_1 + 2].y1 = t.pa[j_1].y;
									cylsup[3 * j_1 + 2].z1 = t.pa[j_1].z;
								}
							}
							break;
						case YZ_PLANE:
							if (sqrt((t.pa[j_1].y - b[i_1].g.yC)*(t.pa[j_1].y - b[i_1].g.yC) + (t.pa[j_1].z - b[i_1].g.zC)*(t.pa[j_1].z - b[i_1].g.zC)) < b[i_1].g.R_out_cyl + sqrt(epsz1*epsz1 + epsy1*epsy1)) {
								if ((t.pa[j_1].x > b[i_1].g.xC - 3 * epsx1) && (t.pa[j_1].x < b[i_1].g.xC + b[i_1].g.Hcyl + 3 * epsx1)) {
									constr[3 * j_1] = true;//X
									constr[3 * j_1 + 1] = true;//Y
									cylsup[3 * j_1 + 1].bactive = true;
									cylsup[3 * j_1 + 1].iPlane = b[i_1].g.iPlane;
									cylsup[3 * j_1 + 1].Radius = b[i_1].g.R_out_cyl;
									cylsup[3 * j_1 + 1].xC = b[i_1].g.xC;
									cylsup[3 * j_1 + 1].yC = b[i_1].g.yC;
									cylsup[3 * j_1 + 1].zC = b[i_1].g.zC;
									cylsup[3 * j_1 + 1].x1 = t.pa[j_1].x;
									cylsup[3 * j_1 + 1].y1 = t.pa[j_1].y;
									cylsup[3 * j_1 + 1].z1 = t.pa[j_1].z;
									constr[3 * j_1 + 2] = true;//Z
									cylsup[3 * j_1 + 2].bactive = true;
									cylsup[3 * j_1 + 2].iPlane = b[i_1].g.iPlane;
									cylsup[3 * j_1 + 2].Radius = b[i_1].g.R_out_cyl;
									cylsup[3 * j_1 + 2].xC = b[i_1].g.xC;
									cylsup[3 * j_1 + 2].yC = b[i_1].g.yC;
									cylsup[3 * j_1 + 2].zC = b[i_1].g.zC;
									cylsup[3 * j_1 + 2].x1 = t.pa[j_1].x;
									cylsup[3 * j_1 + 2].y1 = t.pa[j_1].y;
									cylsup[3 * j_1 + 2].z1 = t.pa[j_1].z;
								}
							}
							break;
						}

					}
				}
			}
		}
	}

	doublereal* square = new doublereal[3 * t.maxnod];

#pragma omp parallel for
	for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1++) {
		square[i_1] = 0.0;
	}
	bool** zashita_ot_dublirovaniq = new bool*[3 * t.maxnod];
	for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1++) {
		zashita_ot_dublirovaniq[i_1] = new bool[8];
		for (integer i_7 = 0; i_7 < 8; i_7++) zashita_ot_dublirovaniq[i_1][i_7] = true;
	}
	// Вычисляем площадь.
	// Т.к. для строго внутренних узлов площадь завышена в двое  (nvtx проходится в четвертинке дважды),
	// а для граничных узлов такого завышения скорее всего нет то здесь мы избавились от дублирования
	// и площадь теперь вычисляется верно.
	for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
		doublereal hx = 1.0, hy = 1.0, hz = 1.0;
		volume3D(i_1, t.nvtx, t.pa, hx, hy, hz);
		for (integer j = 0; j <= 7; j++) {
			integer j_1 = t.nvtx[j][i_1] - 1;
			//X
			//if (j == 0 || j == 3 || j == 4 || j == 7) {
			if ((j==0||(j==1))&&zashita_ot_dublirovaniq[index_of(j_1, 'x')][0]&& zashita_ot_dublirovaniq[index_of(j_1, 'x')][1]) {
				square[index_of(j_1, 'x')] += 0.25*hy*hz;
				zashita_ot_dublirovaniq[index_of(j_1, 'x')][j] = false;
			}
			if ((j == 2 || (j == 3)) && zashita_ot_dublirovaniq[index_of(j_1, 'x')][2] && zashita_ot_dublirovaniq[index_of(j_1, 'x')][3]) {
				square[index_of(j_1, 'x')] += 0.25*hy*hz;
				zashita_ot_dublirovaniq[index_of(j_1, 'x')][j] = false;
			}
			if ((j == 4 || (j == 5)) && zashita_ot_dublirovaniq[index_of(j_1, 'x')][4] && zashita_ot_dublirovaniq[index_of(j_1, 'x')][5]) {
				square[index_of(j_1, 'x')] += 0.25*hy*hz;
				zashita_ot_dublirovaniq[index_of(j_1, 'x')][j] = false;
			}
			if ((j == 6 || (j == 7)) && zashita_ot_dublirovaniq[index_of(j_1, 'x')][6] && zashita_ot_dublirovaniq[index_of(j_1, 'x')][7]) {
				square[index_of(j_1, 'x')] += 0.25*hy*hz;
				zashita_ot_dublirovaniq[index_of(j_1, 'x')][j] = false;
			}

			//Y
			integer indY = index_of(j_1, 'y');
			if ((j == 0 || (j == 3)) && zashita_ot_dublirovaniq[indY][0] && zashita_ot_dublirovaniq[indY][3]) {
				square[indY] += 0.25*hx*hz;
				zashita_ot_dublirovaniq[indY][j] = false;
			}
			if ((j == 1 || (j == 2)) && zashita_ot_dublirovaniq[indY][1] && zashita_ot_dublirovaniq[indY][2]) {
				square[indY] += 0.25*hx*hz;
				zashita_ot_dublirovaniq[indY][j] = false;
			}
			if ((j == 4 || (j == 7)) && zashita_ot_dublirovaniq[indY][4] && zashita_ot_dublirovaniq[indY][7]) {
				square[indY] += 0.25*hx*hz;
				zashita_ot_dublirovaniq[indY][j] = false;
			}
			if ((j == 5 || (j == 6)) && zashita_ot_dublirovaniq[indY][5] && zashita_ot_dublirovaniq[indY][6]) {
				square[indY] += 0.25*hx*hz;
				zashita_ot_dublirovaniq[indY][j] = false;
			}
			//Z
			integer indZ = index_of(j_1, 'z');
			if ((j == 0 || (j == 4)) && zashita_ot_dublirovaniq[indZ][0] && zashita_ot_dublirovaniq[indZ][4]) {
				square[indZ] += 0.25*hx*hy;
				zashita_ot_dublirovaniq[indZ][j] = false;
			}
			if ((j == 1 || (j == 5)) && zashita_ot_dublirovaniq[indZ][1] && zashita_ot_dublirovaniq[indZ][5]) {
				square[indZ] += 0.25*hx*hy;
				zashita_ot_dublirovaniq[indZ][j] = false;
			}
			if ((j == 7 || (j == 3)) && zashita_ot_dublirovaniq[indZ][7] && zashita_ot_dublirovaniq[indZ][3]) {
				square[indZ] += 0.25*hx*hy;
				zashita_ot_dublirovaniq[indZ][j] = false;
			}
			if ((j == 6 || (j == 2)) && zashita_ot_dublirovaniq[indZ][6] && zashita_ot_dublirovaniq[indZ][2]) {
				square[indZ] += 0.25*hx*hy;
				zashita_ot_dublirovaniq[indZ][j] = false;
			}
		}
	}

	for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1++) {
		delete[] zashita_ot_dublirovaniq[i_1];
		zashita_ot_dublirovaniq[i_1] = nullptr;
	}
	delete[] zashita_ot_dublirovaniq;

	

	// Деформации вызванные тепловым расширением.
	if (bThermalStress) {

		TOCHKA center_source;
		doublereal dcount = 0.0;
		center_source.x = 0.0;
		center_source.y = 0.0;
		center_source.z = 0.0;
		for (integer i_43 = 0; i_43 < t.maxelm; i_43++) {
			/*if (fabs(t.Sc[i_43]) > 1.0e-30) {
				TOCHKA p;
				center_cord3D(i_43, t.nvtx, t.pa, p, 100);
				dcount += 1.0;
				center_source.x += p.x;
				center_source.y += p.y;
				center_source.z += p.z;
			}*/
			// Центр масс
			doublereal hx = 1.0, hy = 1.0, hz = 1.0;
			volume3D(i_43, t.nvtx, t.pa, hx, hy, hz);

			TOCHKA p;
			center_cord3D(i_43, t.nvtx, t.pa, p, 100);
			dcount += t.prop[RHO][i_43] * hx*hy*hz;
			center_source.x += p.x*t.prop[RHO][i_43]*hx*hy*hz;
			center_source.y += p.y*t.prop[RHO][i_43] * hx*hy*hz;
			center_source.z += p.z*t.prop[RHO][i_43] * hx*hy*hz;

		}
		//if (dcount > 0.5) 
		integer indf = -1;
		{
			center_source.x /= dcount;
			center_source.y /= dcount;
			center_source.z /= dcount;
			printf("center mass: x=%e y=%e z=%e\n", center_source.x, center_source.y, center_source.z);
			doublereal dist_m = 1.0e30;
			
			for (integer i_43 = 0; i_43 < t.maxelm; i_43++) {
				TOCHKA p;
				center_cord3D(i_43, t.nvtx, t.pa, p, 100);
				if (sqrt((p.x - center_source.x)*(p.x - center_source.x) + (p.y - center_source.y)*(p.y - center_source.y) + (p.z - center_source.z)*(p.z - center_source.z)) < dist_m) {
					dist_m = sqrt((p.x - center_source.x)*(p.x - center_source.x) + (p.y - center_source.y)*(p.y - center_source.y) + (p.z - center_source.z)*(p.z - center_source.z));
					indf = i_43;
				}
			}
		}

		doublereal temp_cm0 = t.potent[indf];
		printf("temperature center mass=%e\n", temp_cm0);

		//else {
			//printf("center power model not be found. Power source not be found...\n");
			//system("PAUSE");
			//exit(1);
		//}

		doublereal* Tx_transform = new doublereal[t.maxnod];
		doublereal* Ty_transform = new doublereal[t.maxnod];
		doublereal* Tz_transform = new doublereal[t.maxnod];
		doublereal* T_transform = new doublereal[t.maxnod];

		doublereal* volume = new doublereal[t.maxnod];


		// Вычисление градиентов температуры.		
		
		doublereal *gradTx1 = nullptr;
		doublereal *gradTy1 = nullptr;
		doublereal *gradTz1 = nullptr;
		gradTx1 = new doublereal[t.maxelm + t.maxbound];
		gradTy1 = new doublereal[t.maxelm + t.maxbound];
		gradTz1 = new doublereal[t.maxelm + t.maxbound];

		doublereal* potent2 = nullptr;
		doublereal *gradTx2 = nullptr;
		doublereal *gradTy2 = nullptr;
		doublereal *gradTz2 = nullptr;
		potent2 = new doublereal[t.maxelm + t.maxbound];
		gradTx2 = new doublereal[t.maxelm + t.maxbound];
		gradTy2 = new doublereal[t.maxelm + t.maxbound];
		gradTz2 = new doublereal[t.maxelm + t.maxbound];

		

		// инициализация нулём.
#pragma omp parallel for
		for (integer i = 0; i<t.maxelm + t.maxbound; i++) {
			

			gradTx1[i] = 0.0;
			gradTy1[i] = 0.0;
			gradTz1[i] = 0.0;

			potent2[i] = 0.0;
			gradTx2[i] = 0.0;
			gradTy2[i] = 0.0;
			gradTz2[i] = 0.0;

		}



		// нахождение градиентов.
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, gradTx1, gradTy1, gradTz1, t.ilevel_alice);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, gradTx1, gradTy1, gradTz1, t.ilevel_alice);
		}
		
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm + t.maxbound; i++) {
			potent2[i] = t.potent[i] - operatingtemperature;
		}
		// нахождение градиентов.
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gaussTemperature(i, potent2, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, gradTx2, gradTy2, gradTz2, t.ilevel_alice);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gaussTemperature(i, potent2, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, gradTx2, gradTy2, gradTz2, t.ilevel_alice);
		}


		// Метод линейного порядка.
		doublereal min_x = 1e60;
		doublereal min_y = 1e60;
		doublereal min_z = 1e60;
		doublereal max_x = -1e60;
		doublereal max_y = -1e60;
		doublereal max_z = -1e60;

		for (integer i = 0; i < t.maxnod; i++) {
			if (t.pa[i].x < min_x) {
				min_x = t.pa[i].x;
			}
			if (t.pa[i].y < min_y) {
				min_y = t.pa[i].y;
			}
			if (t.pa[i].z < min_z) {
				min_z = t.pa[i].z;
			}
			if (t.pa[i].x > max_x) {
				max_x = t.pa[i].x;
			}
			if (t.pa[i].y > max_y) {
				max_y = t.pa[i].y;
			}
			if (t.pa[i].z > max_z) {
				max_z = t.pa[i].z;
			}
		}

		//min_x *= 1.2;
		//min_y *= 1.2;
		//min_z *= 1.2;
		doublereal minx1 = min_x, miny1 = min_y, minz1 = min_z;
		doublereal maxx1 = max_x, maxy1 = max_y, maxz1 = max_z;
			

		// Силы вызванные тепловыми деформациями.
	
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxnod; i_1++) {
			volume[i_1] = 0.0; // inicialization.
		}
		doublereal* YoungModule = new doublereal[t.maxnod];
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxnod; i_1++) {
			YoungModule[i_1] = 0.0; // inicialization.
		}	
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxnod; i_1++) {
			Tx_transform[i_1] = 0.0; // inicialization.
		}
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxnod; i_1++) {
			Ty_transform[i_1] = 0.0; // inicialization.
		}	
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxnod; i_1++) {
			Tz_transform[i_1] = 0.0; // inicialization.
		}
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxnod; i_1++) {
			T_transform[i_1] = 0.0; // inicialization.
		}

		if (0) {
			/* // 28.08.2020 Устаревший код.
			for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
				doublereal hx = 1.0, hy = 1.0, hz = 1.0;
				volume3D(i_1, t.nvtx, t.pa, hx, hy, hz);

				if (hx <= 0.0) {
					printf("negative hx\n");
					system("PAUSE");
				}
				if (hy <= 0.0) {
					printf("negative hy\n");
					system("PAUSE");
				}
				if (hz <= 0.0) {
					printf("negative hz\n");
					system("PAUSE");
				}

				//doublereal mu, lambda; // Коэффициенты Ламе.
				doublereal beta_t_solid; //  коэффициент линейного теплового расширения.

				//mu = t.prop[MU_LAME][i_1];
				//lambda = t.prop[LAMBDA_LAME][i_1];
				beta_t_solid = t.prop[BETA_T_MECHANICAL][i_1];// Коэффициент линейного теплового расширения
				//printf("beta_t_solid=%e\n", beta_t_solid); getchar(); // debug Ok.

				for (integer j_1 = 0; j_1 <= 7; j_1++) {
					volume[t.nvtx[j_1][i_1] - 1] += fabs(0.125 * hx * hy * hz);

					YoungModule[t.nvtx[j_1][i_1] - 1] += 0.125 * hx * hy * hz * ((mu * (3 * lambda + 2 * mu)) / (lambda + mu));

					//Tx_transform[t.nvtx[j_1][i_1] - 1] += 0.125*hx*hy*hz*Tx[i_1];

					//Ty_transform[t.nvtx[j_1][i_1] - 1] += 0.125*hx*hy*hz*Ty[i_1];

					//Tz_transform[t.nvtx[j_1][i_1] - 1] += 0.125*hx*hy*hz*Tz[i_1];

					// умножаем на объём, на модуль Юнга, на градиент, на коэффициент линейного теплового расширения
					// Минус градиент !!! чтобы расширялось в направлении от источника тепла.
					//Tx_transform[t.nvtx[j_1][i_1] - 1] -= beta_t_solid * 0.125 * hx * hy * hz * 0.125 * hx * hy * hz * ((mu * (3 * lambda + 2 * mu)) / (lambda + mu)) * gradTx[i_1];

					//Ty_transform[t.nvtx[j_1][i_1] - 1] -= beta_t_solid * 0.125 * hx * hy * hz * 0.125 * hx * hy * hz * ((mu * (3 * lambda + 2 * mu)) / (lambda + mu)) * gradTy[i_1];

					//Tz_transform[t.nvtx[j_1][i_1] - 1] -= beta_t_solid * 0.125 * hx * hy * hz * 0.125 * hx * hy * hz * ((mu * (3 * lambda + 2 * mu)) / (lambda + mu)) * gradTz[i_1];

					//T_transform[t.nvtx[j_1][i_1] - 1] += 0.125 * hx * hy * hz * (t.potent[i_1]);

				}

			}

			for (integer i_1 = 0; i_1 < t.maxnod; i_1++) {
				YoungModule[i_1] = YoungModule[i_1] / volume[i_1];
				Tx_transform[i_1] = Tx_transform[i_1] / volume[i_1];
				Ty_transform[i_1] = Ty_transform[i_1] / volume[i_1];
				Tz_transform[i_1] = Tz_transform[i_1] / volume[i_1];
				T_transform[i_1] = T_transform[i_1] / volume[i_1];
			}
			// E*vol*beta_t_solid*gradT
			*/

		}
		else if (0) {// На основе температуры хранящейся в центре ячейки (Метод Контрольного Объёма).
			// 22.08.2020 
			// Такое представление не подходит по видимому из за очень низкой
			// точности процедуры переноса из центра ячейки в вершины. 
			// Полностью переходим на хранение в вершинах ячейки для задач Механики.

			// Метод линейного порядка.
			doublereal min_x = 1e60;
			doublereal min_y = 1e60;
			doublereal min_z = 1e60;
			doublereal max_x = -1e60;
			doublereal max_y = -1e60;
			doublereal max_z = -1e60;

			for (integer i = 0; i < t.maxnod; i++) {
				if (t.pa[i].x < min_x) {
					min_x = t.pa[i].x;
				}
				if (t.pa[i].y < min_y) {
					min_y = t.pa[i].y;
				}
				if (t.pa[i].z < min_z) {
					min_z = t.pa[i].z;
				}
				if (t.pa[i].x > max_x) {
					max_x = t.pa[i].x;
				}
				if (t.pa[i].y > max_y) {
					max_y = t.pa[i].y;
				}
				if (t.pa[i].z > max_z) {
					max_z = t.pa[i].z;
				}
			}

			//min_x *= 1.2;
			//min_y *= 1.2;
			//min_z *= 1.2;
			doublereal minx1 = min_x, miny1=min_y, minz1=min_z;
			doublereal maxx1 = max_x, maxy1 = max_y, maxz1 = max_z;



			min_x = 1.05 * fabs(max_x - min_x);
			if (min_x < 1.0e-30) {
				min_x = 1.05 * fabs(max_x);
			}
			min_y = 1.05 * fabs(max_y - min_y);
			if (min_y < 1.0e-30) {
				min_y = 1.05 * fabs(max_y);
			}
			min_z = 1.05 * fabs(max_z - min_z);
			if (min_z < 1.0e-30) {
				min_z = 1.05 * fabs(max_z);
			}
			doublereal eps_mashine = 1.0e-308; // double

			doublereal* vol = new doublereal[t.maxnod];
			doublereal* potent_loc=new doublereal[t.maxelm];
			for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
				doublereal hx = 1.0, hy = 1.0, hz = 1.0;
				volume3D(i_1, t.nvtx, t.pa, hx, hy, hz);

				TOCHKA p;
				center_cord3D(i_1, t.nvtx, t.pa, p, 100);

				doublereal dist = sqrt((p.x - center_source.x)*(p.x - center_source.x) + (p.y - center_source.y)*(p.y - center_source.y) + (p.z - center_source.z)*(p.z - center_source.z));


				
				doublereal directional_forcex = 1.0;
				if (gradTx1[i_1] > 0.0) {
					//directional_forcex = -1.0;
				}
				doublereal directional_forcey = 1.0;
				if (gradTy1[i_1] > 0.0) {
					//directional_forcey = -1.0;
				}
				doublereal directional_forcez = 1.0;
				if (gradTz1[i_1] > 0.0) {
					//directional_forcez = -1.0;
				}
				if (p.x < center_source.x) {
					directional_forcex = -1.0;
				}
				if (p.y < center_source.y) {
					directional_forcey = -1.0;
				}
				if (p.z < center_source.z) {
					directional_forcez = -1.0;
				}
				

				if (hx <= 0.0) {
					printf("negative hx\n");
					system("PAUSE");
					exit(1);
				}
				if (hy <= 0.0) {
					printf("negative hy\n");
					system("PAUSE");
					exit(1);
				}
				if (hz <= 0.0) {
					printf("negative hz\n");
					system("PAUSE");
					exit(1);
				}

				doublereal  beta_t_solid; // Коэффициенты Ламе, коэффициент линейного теплового расширения.

				//mu = t.prop[MU_LAME][i_1];
				//lambda = t.prop[LAMBDA_LAME][i_1];
				beta_t_solid = t.prop[BETA_T_MECHANICAL][i_1];
				doublereal beta_t_solid_x = t.prop[MULT_BETA_T_MECHANICAL_X][i_1] * t.prop[BETA_T_MECHANICAL][i_1];// Коэффициент линейного теплового расширения 1/K.
				doublereal beta_t_solid_y = t.prop[MULT_BETA_T_MECHANICAL_Y][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
				doublereal beta_t_solid_z = t.prop[MULT_BETA_T_MECHANICAL_Z][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
				doublereal Ex = t.prop[MULT_YOUNG_MODULE_X][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
				doublereal Ey = t.prop[MULT_YOUNG_MODULE_Y][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
				doublereal Ez = t.prop[MULT_YOUNG_MODULE_Z][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
				doublereal E =  t.prop[YOUNG_MODULE][i_1];
				doublereal nu_x = t.prop[MULT_POISSON_RATIO_YZ][i_1] * t.prop[POISSON_RATIO][i_1];
				doublereal nu_y = t.prop[MULT_POISSON_RATIO_XZ][i_1] * t.prop[POISSON_RATIO][i_1];
				doublereal nu_z = t.prop[MULT_POISSON_RATIO_XY][i_1] * t.prop[POISSON_RATIO][i_1];
				doublereal nu= t.prop[POISSON_RATIO][i_1];
				doublereal mu_x = (1 - 2.0*nu_x) / (2.0*(1.0 - nu_x));
				doublereal mu_y = (1 - 2.0*nu_y) / (2.0*(1.0 - nu_y));
				doublereal mu_z = (1 - 2.0*nu_z) / (2.0*(1.0 - nu_z));
				doublereal mu = (1 - 2.0*nu) / (2.0*(1.0 - nu));
				doublereal lambda_x = nu_x / (1.0 - nu_x);
				doublereal lambda_y = nu_y / (1.0 - nu_y);
				doublereal lambda_z = nu_z / (1.0 - nu_z);
				doublereal lambda = nu / (1.0 - nu);
				doublereal Em = E * (1.0 - nu) / ((1.0+nu)*(1.0-2.0*nu));
				doublereal Em_x = Ex * (1.0 - nu_x) / ((1.0 + nu_x)*(1.0 - 2.0*nu_x));
				doublereal Em_y = Ey * (1.0 - nu_y) / ((1.0 + nu_y)*(1.0 - 2.0*nu_y));
				doublereal Em_z = Ez * (1.0 - nu_z) / ((1.0 + nu_z)*(1.0 - 2.0*nu_z));

				// E==((mu * (3 * lambda + 2 * mu)) / (lambda + mu));

				// Ненужно домножать на объём !!!
				//potent_loc[i_1] = -(hx * hy * hz) * ( beta_t_solid * E * gradTx[i_1]);
				// Сила это тензор. Без учёта ортотропности материала.
				//potent_loc[i_1] = -(hx * hy * hz) * (beta_t_solid * Em * gradTx[i_1]+
					//beta_t_solid * Em * lambda * gradTy[i_1] +
					//beta_t_solid * Em * lambda * gradTz[i_1]);
				// Сила это тензор. С учётом ортотропности материала.
				//potent_loc[i_1] = -(hx * hy * hz) * (beta_t_solid_x * Em_x * gradTx[i_1] +
					//beta_t_solid_y * Em_y * lambda_y * gradTy[i_1] +
					//beta_t_solid_z * Em_z * lambda_z * gradTz[i_1]);
				
				//potent_loc[i_1] = E * hy*hz*beta_t_solid_x*(t.potent[i_1] - operatingtemperature)*(p.x - center_source.x) / dist;
				//epsilon(doublereal current, doublereal min, doublereal max, doublereal center)
				//potent_loc[i_1] = E * hy*hz*beta_t_solid_x*(t.potent[i_1] - operatingtemperature)*epsilon(p.x,minx1,maxx1, center_source.x);
				//temp_cm0
			    //potent_loc[i_1] = E * hy*hz*beta_t_solid_x*fabs(t.potent[i_1] - temp_cm0)*epsilon(p.x, minx1, maxx1, center_source.x);
				//potent_loc[i_1] =  (hx * hy * hz) *Em * (beta_t_solid_x + lambda * (beta_t_solid_y + beta_t_solid_z))*gradTx1[i_1];

				// рабочий вариант.
				potent_loc[i_1] = (hx * hy * hz) *Em * (beta_t_solid_x + lambda * (beta_t_solid_y + beta_t_solid_z))*gradTx2[i_1];
				

				// Вычитал в статье https://cyberleninka.ru/article/v/vliyanie-koeffitsienta-teplovogo-rasshireniya-na-termouprugie-napryazheniya-v-keramicheskoy-probke
				//potent_loc[i_1] = (hx*hy*hz)*((beta_t_solid * E * (t.potent[i_1]-operatingtemperature))/(hx));
			}
			//ZERO_ORDER_RECONSTRUCT(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Tx_transform, eps_mashine,
				//potent_loc, nullptr, nullptr, true);
			SECOND_ORDER_QUADRATIC_RECONSTRUCTA(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Tx_transform, min_x, min_y, min_z, potent_loc, t, eps_mashine, false);
			for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
				doublereal hx = 1.0, hy = 1.0, hz = 1.0;
				volume3D(i_1, t.nvtx, t.pa, hx, hy, hz);

				TOCHKA p;
				center_cord3D(i_1, t.nvtx, t.pa, p, 100);

				doublereal dist = sqrt((p.x - center_source.x)*(p.x - center_source.x) + (p.y - center_source.y)*(p.y - center_source.y) + (p.z - center_source.z)*(p.z - center_source.z));

				
				doublereal directional_forcex = 1.0;
				if (gradTx1[i_1] > 0.0) {
					//directional_forcex = -1.0;
				}
				doublereal directional_forcey = 1.0;
				if (gradTy1[i_1] > 0.0) {
					//directional_forcey = -1.0;
				}
				doublereal directional_forcez = 1.0;
				if (gradTz1[i_1] > 0.0) {
					//directional_forcez = -1.0;
				}
				if (p.x < center_source.x) {
					directional_forcex = -1.0;
				}
				if (p.y < center_source.y) {
					directional_forcey = -1.0;
				}
				if (p.z < center_source.z) {
					directional_forcez = -1.0;
				}

				if (hx <= 0.0) {
					printf("negative hx\n");
					system("PAUSE");
				}
				if (hy <= 0.0) {
					printf("negative hy\n");
					system("PAUSE");
				}
				if (hz <= 0.0) {
					printf("negative hz\n");
					system("PAUSE");
				}

				doublereal /*mu, lambda,*/ beta_t_solid; // Коэффициенты Ламе, коэффициент линейного теплового расширения.

				//mu = t.prop[MU_LAME][i_1];
				//lambda = t.prop[LAMBDA_LAME][i_1];
				beta_t_solid = t.prop[BETA_T_MECHANICAL][i_1];
				doublereal beta_t_solid_x = t.prop[MULT_BETA_T_MECHANICAL_X][i_1] * t.prop[BETA_T_MECHANICAL][i_1];// Коэффициент линейного теплового расширения 1/K.
				doublereal beta_t_solid_y = t.prop[MULT_BETA_T_MECHANICAL_Y][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
				doublereal beta_t_solid_z = t.prop[MULT_BETA_T_MECHANICAL_Z][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
				doublereal Ex = t.prop[MULT_YOUNG_MODULE_X][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
				doublereal Ey = t.prop[MULT_YOUNG_MODULE_Y][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
				doublereal Ez = t.prop[MULT_YOUNG_MODULE_Z][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
				doublereal E = t.prop[YOUNG_MODULE][i_1];
				doublereal nu_x = t.prop[MULT_POISSON_RATIO_YZ][i_1] * t.prop[POISSON_RATIO][i_1];
				doublereal nu_y = t.prop[MULT_POISSON_RATIO_XZ][i_1] * t.prop[POISSON_RATIO][i_1];
				doublereal nu_z = t.prop[MULT_POISSON_RATIO_XY][i_1] * t.prop[POISSON_RATIO][i_1];
				doublereal nu = t.prop[POISSON_RATIO][i_1];
				doublereal mu_x = (1 - 2.0*nu_x) / (2.0*(1.0 - nu_x));
				doublereal mu_y = (1 - 2.0*nu_y) / (2.0*(1.0 - nu_y));
				doublereal mu_z = (1 - 2.0*nu_z) / (2.0*(1.0 - nu_z));
				doublereal mu = (1 - 2.0*nu) / (2.0*(1.0 - nu));
				doublereal lambda_x = nu_x / (1.0 - nu_x);
				doublereal lambda_y = nu_y / (1.0 - nu_y);
				doublereal lambda_z = nu_z / (1.0 - nu_z);
				doublereal lambda = nu / (1.0 - nu);
				doublereal Em = E * (1.0 - nu) / ((1.0 + nu)*(1.0 - 2.0*nu));
				doublereal Em_x = Ex * (1.0 - nu_x) / ((1.0 + nu_x)*(1.0 - 2.0*nu_x));
				doublereal Em_y = Ey * (1.0 - nu_y) / ((1.0 + nu_y)*(1.0 - 2.0*nu_y));
				doublereal Em_z = Ez * (1.0 - nu_z) / ((1.0 + nu_z)*(1.0 - 2.0*nu_z));
				// E==((mu * (3 * lambda + 2 * mu)) / (lambda + mu));

				// Ненужно домножать на объём !!!
				//potent_loc[i_1] = -(hx * hy * hz) * (beta_t_solid * E * gradTy[i_1]);
				// Сила это тензор. Без учёта ортотропности материала.
				//potent_loc[i_1] = -(hx * hy * hz) * (beta_t_solid * Em * lambda * gradTx[i_1]+
					//beta_t_solid * Em * gradTy[i_1] + 
					//beta_t_solid * Em * lambda * gradTz[i_1]);
				// Сила это тензор. С учётом ортотропности материала.
				//potent_loc[i_1] = -(hx * hy * hz) * (beta_t_solid_x * Em_x * lambda_x * gradTx[i_1] +
				//	beta_t_solid_y * Em_y * gradTy[i_1] +
				//	beta_t_solid_z * Em_z * lambda_z * gradTz[i_1]);
				

				//potent_loc[i_1] = E * hx*hz*beta_t_solid_y*(t.potent[i_1] - operatingtemperature)*(p.y - center_source.y) / dist;
				//potent_loc[i_1] = E * hx*hz*beta_t_solid_y*(t.potent[i_1] - operatingtemperature)*epsilon(p.y, miny1, maxy1, center_source.y);
				//temp_cm0
				//potent_loc[i_1] = E * hx*hz*beta_t_solid_y*fabs(t.potent[i_1] - temp_cm0)*epsilon(p.y, miny1, maxy1, center_source.y);
				
				// рабочий вариант.
				potent_loc[i_1] = (hx * hy * hz) *Em * (beta_t_solid_y + lambda * (beta_t_solid_x + beta_t_solid_z))*gradTy2[i_1];
				

				// Вычитал в статье https://cyberleninka.ru/article/v/vliyanie-koeffitsienta-teplovogo-rasshireniya-na-termouprugie-napryazheniya-v-keramicheskoy-probke
				//potent_loc[i_1] = (hx * hy * hz) * ((beta_t_solid * E * (t.potent[i_1] - operatingtemperature)) / (hy));
			}
			//ZERO_ORDER_RECONSTRUCT(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Ty_transform, eps_mashine,
				//potent_loc, nullptr, nullptr, true);
			SECOND_ORDER_QUADRATIC_RECONSTRUCTA(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Ty_transform, min_x, min_y, min_z, potent_loc, t, eps_mashine, false);
			for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
				doublereal hx = 1.0, hy = 1.0, hz = 1.0;
				volume3D(i_1, t.nvtx, t.pa, hx, hy, hz);

				
				TOCHKA p;
				center_cord3D(i_1, t.nvtx, t.pa, p, 100);

				doublereal dist = sqrt((p.x - center_source.x)*(p.x - center_source.x) + (p.y - center_source.y)*(p.y - center_source.y) + (p.z - center_source.z)*(p.z - center_source.z));


				doublereal directional_forcex = 1.0;
				if (gradTx1[i_1] > 0.0) {
					//directional_forcex = -1.0;
				}
				doublereal directional_forcey = 1.0;
				if (gradTy1[i_1] > 0.0) {
					//directional_forcey = -1.0;
				}
				doublereal directional_forcez = 1.0;
				if (gradTz1[i_1] > 0.0) {
					//directional_forcez = -1.0;
				}
				if (p.x < center_source.x) {
					directional_forcex = -1.0;
				}
				if (p.y < center_source.y) {
					directional_forcey = -1.0;
				}
				if (p.z < center_source.z) {
					directional_forcez = -1.0;
				}

				if (hx <= 0.0) {
					printf("negative hx\n");
					system("PAUSE");
				}
				if (hy <= 0.0) {
					printf("negative hy\n");
					system("PAUSE");
				}
				if (hz <= 0.0) {
					printf("negative hz\n");
					system("PAUSE");
				}

				doublereal /*mu, lambda,*/ beta_t_solid; // Коэффициенты Ламе, коэффициент линейного теплового расширения.

				//mu = t.prop[MU_LAME][i_1];
				//lambda = t.prop[LAMBDA_LAME][i_1];
				beta_t_solid = t.prop[BETA_T_MECHANICAL][i_1];
				doublereal beta_t_solid_x = t.prop[MULT_BETA_T_MECHANICAL_X][i_1] * t.prop[BETA_T_MECHANICAL][i_1];// Коэффициент линейного теплового расширения 1/K.
				doublereal beta_t_solid_y = t.prop[MULT_BETA_T_MECHANICAL_Y][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
				doublereal beta_t_solid_z = t.prop[MULT_BETA_T_MECHANICAL_Z][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
				doublereal Ex = t.prop[MULT_YOUNG_MODULE_X][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
				doublereal Ey = t.prop[MULT_YOUNG_MODULE_Y][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
				doublereal Ez = t.prop[MULT_YOUNG_MODULE_Z][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
				doublereal E = t.prop[YOUNG_MODULE][i_1];
				doublereal nu_x = t.prop[MULT_POISSON_RATIO_YZ][i_1] * t.prop[POISSON_RATIO][i_1];
				doublereal nu_y = t.prop[MULT_POISSON_RATIO_XZ][i_1] * t.prop[POISSON_RATIO][i_1];
				doublereal nu_z = t.prop[MULT_POISSON_RATIO_XY][i_1] * t.prop[POISSON_RATIO][i_1];
				doublereal nu = t.prop[POISSON_RATIO][i_1];
				doublereal mu_x = (1 - 2.0*nu_x) / (2.0*(1.0 - nu_x));
				doublereal mu_y = (1 - 2.0*nu_y) / (2.0*(1.0 - nu_y));
				doublereal mu_z = (1 - 2.0*nu_z) / (2.0*(1.0 - nu_z));
				doublereal mu = (1 - 2.0*nu) / (2.0*(1.0 - nu));
				doublereal lambda_x = nu_x / (1.0 - nu_x);
				doublereal lambda_y = nu_y / (1.0 - nu_y);
				doublereal lambda_z = nu_z / (1.0 - nu_z);
				doublereal lambda = nu / (1.0 - nu);
				doublereal Em = E * (1.0 - nu) / ((1.0 + nu)*(1.0 - 2.0*nu));
				doublereal Em_x = Ex * (1.0 - nu_x) / ((1.0 + nu_x)*(1.0 - 2.0*nu_x));
				doublereal Em_y = Ey * (1.0 - nu_y) / ((1.0 + nu_y)*(1.0 - 2.0*nu_y));
				doublereal Em_z = Ez * (1.0 - nu_z) / ((1.0 + nu_z)*(1.0 - 2.0*nu_z));
				// E==((mu * (3 * lambda + 2 * mu)) / (lambda + mu));


				// Ненужно домножать на объём !!!
				//potent_loc[i_1] = -(hx * hy * hz) * (beta_t_solid * E * gradTz[i_1]);
				// Сила это тензор. Без учёта ортотропности материала.
				//potent_loc[i_1] = -(hx * hy * hz) * (beta_t_solid * Em * lambda * gradTx[i_1]+
					//beta_t_solid * Em * lambda * gradTy[i_1] + 
					//beta_t_solid * Em * gradTz[i_1]);
				// Сила это тензор. С учётом ортотропности материала.
				//potent_loc[i_1] = -(hx * hy * hz) * (beta_t_solid_x * Em_x * lambda_x * gradTx[i_1] +
					//beta_t_solid_y * Em_y * lambda_y * gradTy[i_1] +
					//beta_t_solid_z * Em_z * gradTz[i_1]);
				

				//potent_loc[i_1] = E * hx*hy*beta_t_solid_z*(t.potent[i_1] - operatingtemperature)*(p.z - center_source.z) / dist;
				//potent_loc[i_1] = E * hx*hy*beta_t_solid_z*(t.potent[i_1] - operatingtemperature)*epsilon(p.z, minz1, maxz1, center_source.z);
				//temp_cm0
				//potent_loc[i_1] = E * hx*hy*beta_t_solid_z*fabs(t.potent[i_1] - temp_cm0)*epsilon(p.z, minz1, maxz1, center_source.z);
				
				// рабочий вариант.
				potent_loc[i_1] = (hx * hy * hz) *Em * (beta_t_solid_z + lambda * (beta_t_solid_x + beta_t_solid_y))*gradTz2[i_1];
			
				// Вычитал в статье https://cyberleninka.ru/article/v/vliyanie-koeffitsienta-teplovogo-rasshireniya-na-termouprugie-napryazheniya-v-keramicheskoy-probke
				//potent_loc[i_1] = (hx * hy * hz) * ((beta_t_solid * E * (t.potent[i_1] - operatingtemperature)) / (hz));
			}
			//ZERO_ORDER_RECONSTRUCT(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Tz_transform, eps_mashine,
				//potent_loc, nullptr, nullptr, true);
			SECOND_ORDER_QUADRATIC_RECONSTRUCTA(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Tz_transform, min_x, min_y, min_z, potent_loc, t, eps_mashine, false);
			delete[] vol;
			delete[] potent_loc;

			//T_transform не заполнен
		}
		else {
		   
            // Абсолютно верно 22.08.2020 Но требует второго температурного солвера.

			for (integer j_1 = 0; j_1 < t.maxelm; j_1++) {

					doublereal hx = 1.0, hy = 1.0, hz = 1.0;
					volume3D(j_1, t.nvtx, t.pa, hx, hy, hz);

				    
					doublereal beta_t_solid = t.prop[BETA_T_MECHANICAL][j_1]; //  Коэффициент линейного теплового расширения.
					doublereal beta_t_solid_x = t.prop[MULT_BETA_T_MECHANICAL_X][j_1] * t.prop[BETA_T_MECHANICAL][j_1];// Коэффициент линейного теплового расширения 1/K.
					doublereal beta_t_solid_y = t.prop[MULT_BETA_T_MECHANICAL_Y][j_1] * t.prop[BETA_T_MECHANICAL][j_1];
					doublereal beta_t_solid_z = t.prop[MULT_BETA_T_MECHANICAL_Z][j_1] * t.prop[BETA_T_MECHANICAL][j_1];
					doublereal Ex = t.prop[MULT_YOUNG_MODULE_X][j_1] * t.prop[YOUNG_MODULE][j_1]; // Модуль Юнга Па.
					doublereal Ey = t.prop[MULT_YOUNG_MODULE_Y][j_1] * t.prop[YOUNG_MODULE][j_1]; // Модуль Юнга Па.
					doublereal Ez = t.prop[MULT_YOUNG_MODULE_Z][j_1] * t.prop[YOUNG_MODULE][j_1]; // Модуль Юнга Па.
					doublereal E = t.prop[YOUNG_MODULE][j_1];
					doublereal nuyz = t.prop[MULT_POISSON_RATIO_YZ][j_1] * t.prop[POISSON_RATIO][j_1];
					doublereal nuxz = t.prop[MULT_POISSON_RATIO_XZ][j_1] * t.prop[POISSON_RATIO][j_1];
					doublereal nuxy = t.prop[MULT_POISSON_RATIO_XY][j_1] * t.prop[POISSON_RATIO][j_1];
					doublereal nuzy = t.prop[MULT_POISSON_RATIO_ZY][j_1] * t.prop[POISSON_RATIO][j_1];
					doublereal nuzx = t.prop[MULT_POISSON_RATIO_ZX][j_1] * t.prop[POISSON_RATIO][j_1];
					doublereal nuyx = t.prop[MULT_POISSON_RATIO_YX][j_1] * t.prop[POISSON_RATIO][j_1];
					doublereal nu = t.prop[POISSON_RATIO][j_1];
					
					doublereal Gxy, Gyz, Gxz;
					if (!t.bActiveShearModule[j_1]) {
						Gxy = Gyz = Gxz = Ex / (2.0 * (1.0 + nuxy));
					}
					else {
						Gyz = t.prop[SHEAR_MODULE_YZ][j_1];
						Gxz = t.prop[SHEAR_MODULE_XZ][j_1];
						Gxy = t.prop[SHEAR_MODULE_XY][j_1];
					}

					// 0.125 09.09.2020
					doublereal Am = 0.125; // 1.0;

					Tx_transform[t.nvtx[0][j_1] - 1] = (Am *hx*hy*hz)*(-1.0 / hx * (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[0][j_1] - 1] - 1.0*
							operatingtemperature) + 1 / hx * (nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
								* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[0][j_1] - 1] - 1.0*operatingtemperature) +
						1 / hx * (nuxy*nuyz + nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz *
							nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[0][j_1] - 1] - 1.0*operatingtemperature));
					Ty_transform[t.nvtx[0][j_1] - 1] = (Am*hx*hy*hz)*(1 / hy * (nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
						* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[0][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hy * (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz
							* nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[0][j_1] - 1] - 1.0*operatingtemperature) + 1 / hy * (nuxz*nuyx + nuyz)
						/ (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*
						beta_t_solid_z*(t_for_Mechanical[t.nvtx[0][j_1] - 1] - 1.0*operatingtemperature));
					Tz_transform[t.nvtx[0][j_1] - 1] = (Am*hx*hy*hz)*(1 / hz * (nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
						* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[0][j_1] - 1] - 1.0*operatingtemperature) +
						1 / hz * (nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz *
							nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[0][j_1] - 1] - 1.0*operatingtemperature) - 1.0 / hz * (nuxy*nuyx - 1.0)
						/ (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*
						beta_t_solid_z*(t_for_Mechanical[t.nvtx[0][j_1] - 1] - 1.0*operatingtemperature));
					Tx_transform[t.nvtx[1][j_1] - 1] = (Am*hx*hy*hz)*(1 / hx * (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
						nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[1][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hx * (nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx +
							nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[1][j_1] - 1] - 1.0*operatingtemperature) - 1.0 / hx * (nuxy*nuyz
								+ nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*
						beta_t_solid_z*(t_for_Mechanical[t.nvtx[1][j_1] - 1] - 1.0*operatingtemperature));
					Ty_transform[t.nvtx[1][j_1] - 1] = (Am*hx*hy*hz)*(1 / hy * (nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
						* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[1][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hy * (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz
							* nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[1][j_1] - 1] - 1.0*operatingtemperature) + 1 / hy * (nuxz*nuyx + nuyz)
						/ (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*
						beta_t_solid_z*(t_for_Mechanical[t.nvtx[1][j_1] - 1] - 1.0*operatingtemperature));
					Tz_transform[t.nvtx[1][j_1] - 1] = (Am*hx*hy*hz)*(1 / hz * (nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
						* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[1][j_1] - 1] - 1.0*operatingtemperature) +
						1 / hz * (nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz *
							nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[1][j_1] - 1] - 1.0*operatingtemperature) - 1.0 / hz * (nuxy*nuyx - 1.0)
						/ (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*
						beta_t_solid_z*(t_for_Mechanical[t.nvtx[1][j_1] - 1] - 1.0*operatingtemperature));
					Tx_transform[t.nvtx[2][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hx * (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[2][j_1] - 1] - 1.0*
							operatingtemperature) + 1 / hx * (nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
								* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[2][j_1] - 1] - 1.0*operatingtemperature) +
						1 / hx * (nuxy*nuyz + nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz *
							nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[2][j_1] - 1] - 1.0*operatingtemperature));
					Ty_transform[t.nvtx[2][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hy * (nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[2][j_1] - 1] - 1.0*
							operatingtemperature) + 1 / hy * (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
								nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[2][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hy * (nuxz*nuyx + nuyz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx +
							nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[2][j_1] - 1] - 1.0*operatingtemperature));
					Tz_transform[t.nvtx[2][j_1] - 1] = (Am*hx*hy*hz)*(1 / hz * (nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
						* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[2][j_1] - 1] - 1.0*operatingtemperature) +
						1 / hz * (nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz *
							nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[2][j_1] - 1] - 1.0*operatingtemperature) - 1.0 / hz * (nuxy*nuyx - 1.0)
						/ (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*
						beta_t_solid_z*(t_for_Mechanical[t.nvtx[2][j_1] - 1] - 1.0*operatingtemperature));
					Tx_transform[t.nvtx[3][j_1] - 1] = (Am*hx*hy*hz)*(1 / hx * (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
						nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[3][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hx * (nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx +
							nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[3][j_1] - 1] - 1.0*operatingtemperature) - 1.0 / hx * (nuxy*nuyz
								+ nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*
						beta_t_solid_z*(t_for_Mechanical[t.nvtx[3][j_1] - 1] - 1.0*operatingtemperature));
					Ty_transform[t.nvtx[3][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hy * (nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[3][j_1] - 1] - 1.0*
							operatingtemperature) + 1 / hy * (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
								nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[3][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hy * (nuxz*nuyx + nuyz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx +
							nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[3][j_1] - 1] - 1.0*operatingtemperature));
					Tz_transform[t.nvtx[3][j_1] - 1] = (Am*hx*hy*hz)*(1 / hz * (nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[3][j_1] - 1] - 1.0*
							operatingtemperature) + 1 / hz * (nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
								* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[3][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hz * (nuxy*nuyx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz
							* nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[3][j_1] - 1] - 1.0*operatingtemperature));
					Tx_transform[t.nvtx[4][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hx * (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[4][j_1] - 1] - 1.0*
							operatingtemperature) + 1 / hx * (nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
								* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[4][j_1] - 1] - 1.0*operatingtemperature) +
						1 / hx * (nuxy*nuyz + nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz *
							nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[4][j_1] - 1] - 1.0*operatingtemperature));
					Ty_transform[t.nvtx[4][j_1] - 1] = (Am*hx*hy*hz)*(1 / hy * (nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[4][j_1] - 1] - 1.0*
							operatingtemperature) - 1.0 / hy * (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
								nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[4][j_1] - 1] - 1.0*
									operatingtemperature) + 1 / hy * (nuxz*nuyx + nuyz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
										* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[4][j_1] - 1] - 1.0*operatingtemperature));
					Tz_transform[t.nvtx[4][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hz * (nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[4][j_1] - 1] - 1.0*
							operatingtemperature) - 1.0 / hz * (nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
								nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[4][j_1] - 1] - 1.0*
									operatingtemperature) + 1 / hz * (nuxy*nuyx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
										nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[4][j_1] - 1] - 1.0*operatingtemperature));
					Tx_transform[t.nvtx[5][j_1] - 1] = (Am*hx*hy*hz)*(1 / hx * (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
						* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[5][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hx * (nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx +
							nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[5][j_1] - 1] - 1.0*operatingtemperature) - 1.0 / hx * (nuxy*nuyz
								+ nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*
						beta_t_solid_z*(t_for_Mechanical[t.nvtx[5][j_1] - 1] - 1.0*operatingtemperature));
					Ty_transform[t.nvtx[5][j_1] - 1] = (Am*hx*hy*hz)*(1 / hy * (nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[5][j_1] - 1] - 1.0*
							operatingtemperature) - 1.0 / hy * (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
								nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[5][j_1] - 1] - 1.0*
									operatingtemperature) + 1 / hy * (nuxz*nuyx + nuyz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
										* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[5][j_1] - 1] - 1.0*operatingtemperature));
					Tz_transform[t.nvtx[5][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hz * (nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[5][j_1] - 1] - 1.0*
							operatingtemperature) - 1.0 / hz * (nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
								nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[5][j_1] - 1] - 1.0*
									operatingtemperature) + 1 / hz * (nuxy*nuyx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
										nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[5][j_1] - 1] - 1.0*operatingtemperature));
					Tx_transform[t.nvtx[6][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hx * (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[6][j_1] - 1] - 1.0*
							operatingtemperature) + 1 / hx * (nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
								* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[6][j_1] - 1] - 1.0*operatingtemperature) +
						1 / hx * (nuxy*nuyz + nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz *
							nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[6][j_1] - 1] - 1.0*operatingtemperature));
					Ty_transform[t.nvtx[6][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hy * (nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[6][j_1] - 1] - 1.0*
							operatingtemperature) + 1 / hy * (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
								nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[6][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hy * (nuxz*nuyx + nuyz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx +
							nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[6][j_1] - 1] - 1.0*operatingtemperature));
					Tz_transform[t.nvtx[6][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hz * (nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[6][j_1] - 1] - 1.0*
							operatingtemperature) - 1.0 / hz * (nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
								nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[6][j_1] - 1] - 1.0*
									operatingtemperature) + 1 / hz * (nuxy*nuyx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
										nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[6][j_1] - 1] - 1.0*operatingtemperature));
					Tx_transform[t.nvtx[7][j_1] - 1] = (Am*hx*hy*hz)*(1 / hx * (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy
						* nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[7][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hx * (nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx +
							nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[7][j_1] - 1] - 1.0*operatingtemperature) - 1.0 / hx * (nuxy*nuyz
								+ nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*
						beta_t_solid_z*(t_for_Mechanical[t.nvtx[7][j_1] - 1] - 1.0*operatingtemperature));
					Ty_transform[t.nvtx[7][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hy * (nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[7][j_1] - 1] - 1.0*
							operatingtemperature) + 1 / hy * (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
								nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[7][j_1] - 1] - 1.0*operatingtemperature)
						- 1.0 / hy * (nuxz*nuyx + nuyz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz * nuzx +
							nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[7][j_1] - 1] - 1.0*operatingtemperature));
					Tz_transform[t.nvtx[7][j_1] - 1] = (Am*hx*hy*hz)*(-1.0 / hz * (nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
						nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ex*beta_t_solid_x*(t_for_Mechanical[t.nvtx[7][j_1] - 1] - 1.0*
							operatingtemperature) - 1.0 / hz * (nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy +
								nuxy * nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ey*beta_t_solid_y*(t_for_Mechanical[t.nvtx[7][j_1] - 1] - 1.0*
									operatingtemperature) + 1 / hz * (nuxy*nuyx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy *
										nuyx + nuxz * nuzx + nuyz * nuzy - 1.0)*Ez*beta_t_solid_z*(t_for_Mechanical[t.nvtx[7][j_1] - 1] - 1.0*operatingtemperature));

									   

			}

        }

		integer ibconstrX = 0;
		integer ibconstrY = 0;
		integer ibconstrZ = 0;
		//for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1++) {
		for (integer i_10 = 0; i_10 < t.maxnod; i_10++) {
			integer i_1 = index_of(i_10, 'x'); // X
			if ((!constr[i_1]) && (!cylsup[i_1].bactive)) {

				//beta_t_solid = t.prop[BETA_T_MECHANICAL][i_1]; см. использование выше.
				// Если узел не зафиксирован.
				doublereal gradT = 0.0;
				gradT = Tx_transform[i_10];

				// Именно добавляем, т.к. изначально могла быть приложена сосредоточенная сила в Ньютонах.
				// Здесь добавляется сила вызванная линейным тепловым расширением.
				// В правой части должна стоять сила в Ньютонах.По размерности.				
				//rthdsd[i_1] += YoungModule[inode]* volume[inode]* betaT*gradT;
				//rthdsd[i_1] += YoungModule[inode] * square[i_1] * betaT*(T_transform[inode] - operatingtemperature);
				//21.02.2019 Уже все учтено, см. выше.

				rthdsd[i_1] += gradT;		

			}
			else {
				rthdsd[i_1] = 0.0;// FIXIT
				ibconstrX++;
			}

			i_1 = index_of(i_10, 'y'); // Y
			if ((!constr[i_1]) && (!cylsup[i_1].bactive)) {

				//beta_t_solid = t.prop[BETA_T_MECHANICAL][i_1]; см. использование выше.
				// Если узел не зафиксирован.
				doublereal gradT = 0.0;
				gradT = Ty_transform[i_10];

				// Именно добавляем, т.к. изначально могла быть приложена сосредоточенная сила в Ньютонах.
				// Здесь добавляется сила вызванная линейным тепловым расширением.
				// В правой части должна стоять сила в Ньютонах.По размерности.				
				//rthdsd[i_1] += YoungModule[inode]* volume[inode]* betaT*gradT;
				//rthdsd[i_1] += YoungModule[inode] * square[i_1] * betaT*(T_transform[inode] - operatingtemperature);
				//21.02.2019 Уже все учтено, см. выше.

				rthdsd[i_1] += gradT;

			}
			else {
				rthdsd[i_1] = 0.0;// FIXIT
				ibconstrY++;
			}

			i_1 = index_of(i_10, 'z'); // Z
			if ((!constr[i_1]) && (!cylsup[i_1].bactive)) {

				//beta_t_solid = t.prop[BETA_T_MECHANICAL][i_1]; см. использование выше.
				// Если узел не зафиксирован.
				doublereal gradT = 0.0;
				gradT = Tz_transform[i_10];

				// Именно добавляем, т.к. изначально могла быть приложена сосредоточенная сила в Ньютонах.
				// Здесь добавляется сила вызванная линейным тепловым расширением.
				// В правой части должна стоять сила в Ньютонах.По размерности.				
				//rthdsd[i_1] += YoungModule[inode]* volume[inode]* betaT*gradT;
				//rthdsd[i_1] += YoungModule[inode] * square[i_1] * betaT*(T_transform[inode] - operatingtemperature);
				//21.02.2019 Уже все учтено, см. выше.

				rthdsd[i_1] += gradT;

			}
			else {
				rthdsd[i_1] = 0.0;// FIXIT
				ibconstrZ++;
			}
		}


		if (btimedep) {
			{//3226

// Метод линейного порядка.
				doublereal min_x = 1e60;
				doublereal min_y = 1e60;
				doublereal min_z = 1e60;
				doublereal max_x = -1e60;
				doublereal max_y = -1e60;
				doublereal max_z = -1e60;

				for (integer i = 0; i < t.maxnod; i++) {
					if (t.pa[i].x < min_x) {
						min_x = t.pa[i].x;
					}
					if (t.pa[i].y < min_y) {
						min_y = t.pa[i].y;
					}
					if (t.pa[i].z < min_z) {
						min_z = t.pa[i].z;
					}
					if (t.pa[i].x > max_x) {
						max_x = t.pa[i].x;
					}
					if (t.pa[i].y > max_y) {
						max_y = t.pa[i].y;
					}
					if (t.pa[i].z > max_z) {
						max_z = t.pa[i].z;
					}
				}

				//min_x *= 1.2;
				//min_y *= 1.2;
				//min_z *= 1.2;



				min_x = 1.05 * fabs(max_x - min_x);
				if (min_x < 1.0e-30) {
					min_x = 1.05 * fabs(max_x);
				}
				min_y = 1.05 * fabs(max_y - min_y);
				if (min_y < 1.0e-30) {
					min_y = 1.05 * fabs(max_y);
				}
				min_z = 1.05 * fabs(max_z - min_z);
				if (min_z < 1.0e-30) {
					min_z = 1.05 * fabs(max_z);
				}
				doublereal eps_mashine = 1.0e-308; // double

				doublereal* vol = new doublereal[t.maxnod];
				doublereal* potent_loc = new doublereal[t.maxelm];
				for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
					doublereal hx = 1.0, hy = 1.0, hz = 1.0;
					volume3D(i_1, t.nvtx, t.pa, hx, hy, hz);

					if (hx <= 0.0) {
						printf("negative hx\n");
						system("PAUSE");
					}
					if (hy <= 0.0) {
						printf("negative hy\n");
						system("PAUSE");
					}
					if (hz <= 0.0) {
						printf("negative hz\n");
						system("PAUSE");
					}

					doublereal rho; // плотность в ячейке кг/м!3.

					rho = t.prop[RHO][i_1];

					// Восстанавливаем значение добавки к правой части
					// из вершин кубика в центр ячейки.
					doublereal tmp = 0.0;
					for (integer k_1 = 0; k_1 < 8; k_1++) {
						integer i_10 = t.nvtx[k_1][i_1] - 1;
						integer i_1a = index_of(i_10, 'x'); // X
						tmp += 0.125*(2.0*uoldtimestep[i_1a] - uolddoubletimestep[i_1a]);
					}
					
					potent_loc[i_1] = (rho*(hx * hy * hz) / (timestep_sizenow*timestep_sizenow))*(tmp);
					
				}
				//ZERO_ORDER_RECONSTRUCT(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Tx_transform, eps_mashine,
				//potent_loc, nullptr, nullptr, true);
				SECOND_ORDER_QUADRATIC_RECONSTRUCTA(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Tx_transform, min_x, min_y, min_z, potent_loc, t, eps_mashine, false);
				for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
					doublereal hx = 1.0, hy = 1.0, hz = 1.0;
					volume3D(i_1, t.nvtx, t.pa, hx, hy, hz);

					if (hx <= 0.0) {
						printf("negative hx\n");
						system("PAUSE");
					}
					if (hy <= 0.0) {
						printf("negative hy\n");
						system("PAUSE");
					}
					if (hz <= 0.0) {
						printf("negative hz\n");
						system("PAUSE");
					}

					doublereal rho; // плотность в ячейке кг/м!3.

					rho = t.prop[RHO][i_1];

					// Восстанавливаем значение добавки к правой части
					// из вершин кубика в центр ячейки.
					doublereal tmp = 0.0;
					for (integer k_1 = 0; k_1 < 8; k_1++) {
						integer i_10 = t.nvtx[k_1][i_1] - 1;
						integer i_1a = index_of(i_10, 'y'); // Y
						tmp += 0.125*(2.0*uoldtimestep[i_1a] - uolddoubletimestep[i_1a]);
					}

					potent_loc[i_1] = (rho*(hx * hy * hz) / (timestep_sizenow*timestep_sizenow))*(tmp);

				}
				//ZERO_ORDER_RECONSTRUCT(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Ty_transform, eps_mashine,
				//potent_loc, nullptr, nullptr, true);
				SECOND_ORDER_QUADRATIC_RECONSTRUCTA(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Ty_transform, min_x, min_y, min_z, potent_loc, t, eps_mashine, false);
				for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {
					doublereal hx = 1.0, hy = 1.0, hz = 1.0;
					volume3D(i_1, t.nvtx, t.pa, hx, hy, hz);

					if (hx <= 0.0) {
						printf("negative hx\n");
						system("PAUSE");
					}
					if (hy <= 0.0) {
						printf("negative hy\n");
						system("PAUSE");
					}
					if (hz <= 0.0) {
						printf("negative hz\n");
						system("PAUSE");
					}

					doublereal rho; // плотность в ячейке кг/м!3.

					rho = t.prop[RHO][i_1];

					// Восстанавливаем значение добавки к правой части
					// из вершин кубика в центр ячейки.
					doublereal tmp = 0.0;
					for (integer k_1 = 0; k_1 < 8; k_1++) {
						integer i_10 = t.nvtx[k_1][i_1] - 1;
						integer i_1a = index_of(i_10, 'z'); // Z
						tmp += 0.125*(2.0*uoldtimestep[i_1a] - uolddoubletimestep[i_1a]);
					}

					potent_loc[i_1] = (rho*(hx * hy * hz) / (timestep_sizenow*timestep_sizenow))*(tmp);

				}
				//ZERO_ORDER_RECONSTRUCT(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Tz_transform, eps_mashine,
				//potent_loc, nullptr, nullptr, true);
				SECOND_ORDER_QUADRATIC_RECONSTRUCTA(t.maxnod, t.maxelm, t.pa, t.nvtx, vol, Tz_transform, min_x, min_y, min_z, potent_loc, t, eps_mashine, false);
				delete[] vol;
				delete[] potent_loc;

				//T_transform не заполнен
			}

			integer ibconstrX = 0;
			integer ibconstrY = 0;
			integer ibconstrZ = 0;
			//for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1++) {
			for (integer i_10 = 0; i_10 < t.maxnod; i_10++) {
				integer i_1 = index_of(i_10, 'x'); // X
				if ((!constr[i_1]) && (!cylsup[i_1].bactive)) {

					// Если узел не зафиксирован.
					doublereal ForceAdditional = 0.0;
					ForceAdditional = Tx_transform[i_10];

					rthdsd[i_1] += ForceAdditional;

				}
				else {
					rthdsd[i_1] = 0.0;// FIXIT
					//ibconstrX++;
				}

				i_1 = index_of(i_10, 'y'); // Y
				if ((!constr[i_1]) && (!cylsup[i_1].bactive)) {

					
					// Если узел не зафиксирован.
					doublereal ForceAdditional = 0.0;
					ForceAdditional = Ty_transform[i_10];

					rthdsd[i_1] += ForceAdditional;

				}
				else {
					rthdsd[i_1] = 0.0;// FIXIT
					//ibconstrY++;
				}

				i_1 = index_of(i_10, 'z'); // Z
				if ((!constr[i_1]) && (!cylsup[i_1].bactive)) {

					// Если узел не зафиксирован.
					doublereal ForceAdditional = 0.0;
					ForceAdditional = Tz_transform[i_10];

					rthdsd[i_1] += ForceAdditional;

				}
				else {
					rthdsd[i_1] = 0.0;// FIXIT
					//ibconstrZ++;
				}
			}
		}

		// При обнаружении фиксации вершины
		// зануляем правую часть во всей ячейке содержащей 
		// данную вершину. 09.08.2020
		// Учёт граничных условий закрепления.
		for (integer j_11 = 0; j_11 < t.maxelm; j_11++) {
			for (integer i_11 = 0; i_11 < 8; i_11++) {
				//for (integer j_1 = 0; j_1 < t.maxnod; j_1++) {
				integer j_1 = t.nvtx[i_11][j_11] - 1;
				/*
				if (constr[3 * j_1+1] && fabs(rthdsd[3 * j_1]) > 0.0) {
					for (integer i_111 = 0; i_111 < 8; i_111++) {
						integer j_111 = t.nvtx[i_111][j_11] - 1;
						constr[3 * j_111] = false;
						rthdsd[3 * j_111] = 0.0;
						constr[3 * j_111 + 2] = 0.0;//Z -fix
					}
				}
				*/
				integer i_1 = index_of(j_1, 'x'); // X
				if (constr[i_1]||(cylsup[i_1].bactive)) {
					for (integer i_1a = 0; i_1a < 8; i_1a++) {
						integer j_1a = t.nvtx[i_1a][j_11] - 1;
						rthdsd[index_of(j_1a, 'x')] = 0.0; //X
					}
				}
				i_1 = index_of(j_1, 'y');
				if (constr[i_1]||(cylsup[i_1].bactive)) {
					for (integer i_1a = 0; i_1a < 8; i_1a++) {
						integer j_1a = t.nvtx[i_1a][j_11] - 1;
						rthdsd[index_of(j_1a, 'y')] = 0.0; //Y 
					}
				}
				i_1 = index_of(j_1, 'z');
				if (constr[i_1]|| (cylsup[i_1].bactive)) {
					for (integer i_1a = 0; i_1a < 8; i_1a++) {
						integer j_1a = t.nvtx[i_1a][j_11] - 1;
						rthdsd[index_of(j_1a, 'z')] = 0.0; //Z 
					}
				}
			}
		}

		printf("FIXIT: X=%lld, Y=%lld, Z=%lld\n", ibconstrX, ibconstrY, ibconstrZ);

		// Освобождение оперативной памяти из под градиентов температуры.
		

		if (gradTx1 != nullptr) {
			delete[] gradTx1;
			gradTx1 = nullptr;
		}
		if (gradTy1 != nullptr) {
			delete[] gradTy1;
			gradTy1 = nullptr;
		}
		if (gradTz1 != nullptr) {
			delete[] gradTz1;
			gradTz1 = nullptr;
		}

		if (potent2 != nullptr) {
			delete[] potent2;
			potent2 = nullptr;
		}

		if (gradTx2 != nullptr) {
			delete[] gradTx2;
			gradTx2 = nullptr;
		}
		if (gradTy2 != nullptr) {
			delete[] gradTy2;
			gradTy2 = nullptr;
		}
		if (gradTz2 != nullptr) {
			delete[] gradTz2;
			gradTz2 = nullptr;
		}	

		if (YoungModule != nullptr) {
			delete[] YoungModule;
			YoungModule = nullptr;
		}
		
		if (volume != nullptr) {
			delete[] volume;
			volume = nullptr;
		}

		if (Tx_transform != nullptr) {
			delete[] Tx_transform;
			Tx_transform = nullptr;
		}
		if (Ty_transform != nullptr) {
			delete[] Ty_transform;
			Ty_transform = nullptr;
		}
		if (Tz_transform != nullptr) {
			delete[] Tz_transform;
			Tz_transform = nullptr;
		}
		if (T_transform != nullptr) {
			delete[] T_transform;
			T_transform = nullptr;
		}
		
	}

	for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1++) {
		//rthdsd[i_1] *= 1.0e-6;
		//if ((!constr[i_1]) && (!cylsup[i_1].bactive)) {
			//integer inode = (integer)(i_1 / 3);
			//rthdsd[i_1] *= (1e-15 / volume[inode]);
		//}
	}

	
	// Умножаем силу на площадь.
	// Сила линейного теплового расширения E*vol*betaT*gradT или
	// E*Square_ortho*betaT*DeltaT, DeltaT=l*gradT.
	for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1++) {
		//rthdsd[i_1] *= square[i_1]; // Newton*m!2.
	}
	delete[] square;




	doublereal** Kmatrix_local = nullptr;
	Kmatrix_local = new doublereal*[24];
	for (integer i_1 = 0; i_1 < 24; i_1++) {
		Kmatrix_local[i_1]= new doublereal[24];
	}
	// инициализация.
	for (integer i_1 = 0; i_1 < 24; i_1++) {
		for (integer j_1 = 0; j_1 < 24; j_1++) {
			Kmatrix_local[i_1][j_1] = 0.0;
		}
	}

	IMatrix sparseS; // Разреженная матрица в формате IMatrix.
	if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER) {
		initIMatrix(&sparseS, 3 * t.maxnod);
	}
	SIMPLESPARSE sparseM; // разреженная матрица
	initsimplesparse(sparseM, 3 * t.maxnod);

	
	// Сборка СЛАУ.
	for (integer ie = 0; ie < t.maxelm; ie++) {
		// инициализация.
		for (integer i_1 = 0; i_1 < 24; i_1++) {
			for (integer j_1 = 0; j_1 < 24; j_1++) {
				Kmatrix_local[i_1][j_1] = 0.0;
			}
		}
		// Сборка локальной матрицы жёсткости.
		// Термоупругость сборка матрицы Жёсткости
		// для шестигранной призмы. 4.08.2017.
		//Thermal_Structural_assemble(ie, t.nvtx,
			//t.pa, t.prop, Kmatrix_local);
		//Thermal_Structural_assemble_Volk2(ie, t.nvtx,
			//t.pa, t.prop, Kmatrix_local, btimedep, timestep_sizenow);
		Thermal_Structural_assemble_Volk3(ie, t.nvtx,
			t.pa, t.prop, Kmatrix_local, btimedep, timestep_sizenow,t.bActiveShearModule[ie]);

		for (integer i_4 = 0; i_4 < 24; i_4++) {
			for (integer j_4 = 0; j_4 < 24; j_4++) {
				//Kmatrix_local[i_4][j_4] *= 1.0e-6;
				//Kmatrix_local[i_4][j_4] *= 4.0;// Компенсация 0.5 в матрице частных производных.
			}
		}

		//if (ie == 3) {
			for (integer i_4 = 0; i_4 < 24; i_4++) {
				integer ipositive = 0;
				for (integer j_4 = 0; j_4 < 24; j_4++) {
					//printf("%1.2f ", Kmatrix_local[i_4][j_4]);
					//printf("%e ", Kmatrix_local[i_4][j_4]);
					if (fabs(Kmatrix_local[i_4][j_4] - Kmatrix_local[j_4][i_4]) > 1.0e-3) {
						// Нарушена симметричность.
						printf("i=%lld j=%lld %e %e", i_4 + 1, j_4 + 1, Kmatrix_local[i_4][j_4], Kmatrix_local[j_4][i_4]);
					}
					//printf("%d \n%d \n%1.9f\n", i_4 + 1, j_4 + 1, Kmatrix_local[i_4][j_4]);
					/*
					if (Kmatrix_local[i_4][j_4] > 0.0) {
						printf("%d \n%d \n1\n",i_4+1,j_4+1);
						ipositive++;
					}
					else {
						printf("%d \n%d \n0\n",i_4+1,j_4+1);
					}
					*/
				}
				//printf("i=%d %d\n",i_4, ipositive);
				//printf("\n");
			}
			//getchar();
		//}
		

		// Сборка правой части 
		// TODO.

		// Собирать ли правую часть СЛАУ.
		bool bsecond_member_of_equation = true;
		bsecond_member_of_equation = false;//
		// Добавление локальной матрицы жёсткости в глобальную.
		if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER)
		{
			
			elembdSparse2(ie, sparseS, t.nvtx,
				constr, rthdsd,
				Kmatrix_local, deformation,
				bsecond_member_of_equation);
			
		}
		else {
			
			//elembdSparse(ie, sparseM, t.nvtx,
				//constr, rthdsd,
				//Kmatrix_local, deformation,
				//bsecond_member_of_equation, cylsup, epsx, epsy, epsz);
				
			elembdSparse_noCylindricalSupport(ie, sparseM, t.nvtx,
				constr, rthdsd,
				Kmatrix_local, deformation,
				bsecond_member_of_equation, cylsup, epsx, epsy, epsz, t.maxnod);

		}
	}

	if (!(iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER)) {
		for (integer i_check = 0; i_check < 3 * t.maxnod; i_check++) {
			if (sparseM.root[i_check] == nullptr) {
				printf("error: zero string %lld \n", i_check);
				system("pause");
			}
			else {
				NONZEROELEM* p;
				p = sparseM.root[i_check];
				if (p != nullptr) {
					NONZEROELEM* q = nullptr;

					/*
					if (i_check == 11) {
						q = p;
						while (q != nullptr) {
							printf("val=%e col_ind=%d row_ind=%d\n", q->aij, q->key, i_check);
							q = q->next;
						}
						getchar();
					}
					*/
					q = nullptr;
					q = p->next;
					//p->next = nullptr;

					while (q != nullptr) {
						p = q;
						if (fabs(p->aij) < 1.0e-300) {
							if (p->key == i_check) {
								printf("%e %lld %lld\n", p->aij, i_check, p->key);
								system("PAUSE");
							}
						}

						//printf(" Dirichlet p-aij=%d\n",p->aij);
						//getchar();
						q = p->next;
						//p->next = nullptr;
						//delete p;
						p = nullptr;
						//M.n--;
					}
					//delete M.root[i];
					//M.root[i] = nullptr;
					//M.n--;
				}
			}
		}
	}
	

	printf("matrix is assemble.\n");
	//getchar();
	// Решение СЛАУ TODO.

	if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER) {
		calculateSPARSEgaussArray(&sparseS, deformation, rthdsd);
	}
	bool bprintmessage = true;
	integer maxiter = 20000; // !!!
	//ICCG(TOTALDEFORMATIONVAR, sparseM, rthdsd, deformation, 3 * t.maxnod, bprintmessage, false, maxiter); //->//
	//doublereal *val = nullptr;
	//integer *col_ind = nullptr, *row_ptr = nullptr;
	//simplesparsetoCRS(sparseM, val, col_ind, row_ptr, (3 * t.maxnod)); // преобразование матрицы из одного формата хранения в другой.
	//simplesparsefree(sparseM, 3 * t.maxnod);

	QuickMemVorst m;
	m.ballocCRSt = false; // Выделять память
	m.bsignalfreeCRSt = true; // и сразу освобождать.

							  // инициализация указателей.
	m.tval = nullptr;
	m.tcol_ind = nullptr;
	m.trow_ptr = nullptr;
	m.tri = nullptr;
	m.troc = nullptr;
	m.ts = nullptr;
	m.tt = nullptr;
	m.tvi = nullptr;
	m.tpi = nullptr;
	m.tdx = nullptr;
	m.tdax = nullptr;
	m.ty = nullptr;
	m.tz = nullptr;
	m.ta = nullptr;
	m.tja = nullptr;
	m.tia = nullptr;
	m.talu = nullptr;
	m.tjlu = nullptr;
	m.tju = nullptr;
	m.tiw = nullptr;
	m.tlevs = nullptr;
	m.tw = nullptr;
	m.tjw = nullptr;
	m.icount_vel = 100000; // очень большое число.

	// Разрешающих свойств данного метода без предобуславливателя явно недостаточно.
	//Bi_CGStabCRS((3 * t.maxnod), val, col_ind, row_ptr, rthdsd, deformation, maxiter);//->//
	
	// BiCGStab + ILU6 сходимость есть.
	if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::BICGSTAB_PLUS_ILU6_SECOND_T_SOLVER) {
		 // BiCGStab +ILU(lfil), lfil=1..6.
		 bool* boundary = nullptr;
	     Bi_CGStab_internal4(sparseM, (3 * t.maxnod), rthdsd, deformation, maxiter, bprintmessage, m, w,lw,boundary, TOTALDEFORMATION);
	}
	// amg1r5 нет сходимости на задачи напряженно-деформированного состояния.
	//amg_loc_memory_Stress(sparseM, (3*t.maxnod), rthdsd, deformation, m);
	if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::CAMG_RUMBA_v0_14_SECOND_T_SOLVER) {
		bool* boundary = nullptr;
		my_agr_amg_loc_memory_Stress(sparseM, (3 * t.maxnod), rthdsd, deformation, m,b,lb,w,lw, boundary, TOTALDEFORMATION);
	}
	if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::AMG1R5_SECOND_T_SOLVER)  {
		bool* boundary = nullptr;
	    //  if (NONE_only_amg1r5==stabilization_amg1r5_algorithm)  -> amg1r5 Руге и Штубена.
		//  if (BiCGStab_plus_amg1r5==stabilization_amg1r5_algorithm)  -> BiCGStab + amg1r5 Хенк Ван дер Ворст + Руге и Штубен.
		//  if (FGMRes_plus_amg1r5==stabilization_amg1r5_algorithm)  -> FGMres + amg1r5 Ю. Саад и Мартин Шульц + Руге и Штубен. // 16.10.2018.
		amg_loc_memory_for_Matrix_assemble2(sparseM, (3 * t.maxnod), rthdsd, deformation, maxiter, bprintmessage, m, w, lw, boundary); // 13.10.2018
	}
	if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::AMGCL_SECONT_T_SOLVER)
	{
#if AMGCL_INCLUDE_IN_MY_PROJECT == 1
		// Denis Demidov AMGCL
		bool* boundary = nullptr;
		amgcl_secondT_solver(sparseM, (3 * t.maxnod),
			rthdsd, deformation, bprintmessage,w,lw,boundary,true);
#else
		std::cout << "ERROR!!! Library AMGCL ddemidov not connected." << std::endl;
		system("PAUSE");
#endif
	}

	// Нужна специальная версия BicgStab+ILU2.

	doublereal deform_max_total = -1.0e30;
	doublereal deform_min_total = +1.0e30;
	doublereal deform_max_x = -1.0e30;
	doublereal deform_min_x = +1.0e30;
	doublereal deform_max_y = -1.0e30;
	doublereal deform_min_y = +1.0e30;
	doublereal deform_max_z = -1.0e30;
	doublereal deform_min_z = +1.0e30;
	for (integer i_1 = 0; i_1 < 3 * t.maxnod; i_1=i_1+3) {
		if (deformation[i_1] > deform_max_x) deform_max_x = deformation[i_1];
		if (deformation[i_1] < deform_max_x) deform_min_x = deformation[i_1];

		if (deformation[i_1+1] > deform_max_y) deform_max_y = deformation[i_1+1];
		if (deformation[i_1+1] < deform_max_y) deform_min_y = deformation[i_1+1];

		if (deformation[i_1+2] > deform_max_z) deform_max_z = deformation[i_1+2];
		if (deformation[i_1+2] < deform_max_z) deform_min_z = deformation[i_1+2];

		doublereal td = sqrt(deformation[i_1]* deformation[i_1]+
			deformation[i_1+1] * deformation[i_1+1] + 
			deformation[i_1+2] * deformation[i_1+2]);
		if (td > deform_max_total) deform_max_total = td;
		if (td < deform_max_total) deform_min_total = td;

		if (0) {
			if (rthdsd[i_1] > 0.0) {
				rthdsd[i_1] = log10(rthdsd[i_1]);
			}
			else if (rthdsd[i_1] < 0.0) {
				rthdsd[i_1] = -log10(fabs(rthdsd[i_1]));
			}
		}
	}
	printf("deformation x directional: min = %e, max=%e\n", deform_min_x, deform_max_x);
	printf("deformation y directional: min = %e, max=%e\n", deform_min_y, deform_max_y);
	printf("deformation z directional: min = %e, max=%e\n", deform_min_z, deform_max_z);
	printf("total deformation directional: min = %e, max=%e\n", deform_min_total, deform_max_total);

	printf("SLAU is solve.\n");
	//getchar();

	if (btimedep) {
		// Подготовка к следующему шагу по времени:
		for (integer i = 0; i < 3 * t.maxnod; i++) {			
			uolddoubletimestep[i] = uoldtimestep[i];
			uoldtimestep[i] = deformation[i];
		}
	}

	// Запись результата для визуализации.
	// деформация: выделение памяти и инициализация нулём.
	init_total_deformation(t);

	// Преобразование векторной величины из хранения в вершинах кубикой в
	// хранение в центре кубиков.
	for (integer j_6 = 0; j_6 < 4; j_6++) {
		Stress2Thermal_vector_translate(t,
			/*rthdsd,*/ deformation, // input
			j_6,
			t.total_deformation[j_6]); // output
	}
	
	// Сохранение деформации.
	// TODO.
	
	
	

	printf("deformation writing.\n");
	//getchar();

	// Освобождение оперативной памяти.
	if (iswitchsolveramg_vs_BiCGstab_plus_ILU6 == SECOND_T_SOLVER_ID_SWITCH::DIRECT_SECOND_T_SOLVER) {
	    freeIMatrix(&sparseS);
	}
	//simplesparsefree(sparseM, 3 * t.maxnod);
	delete[] rthdsd;
	delete[] deformation;
	delete[] constr;
	delete[] cylsup;
	for (integer i_1 = 0; i_1 < 24; i_1++) {
		delete[] Kmatrix_local[i_1];
	}
	delete[] Kmatrix_local;
} // solve_Structural



  // Решение прочностной задачи в 3D.
  // 6 августа 2017. август 2020.
// Природа силы непонятна, непонятно на поверхности она или в объеме.
// Более того непонятно если ее учесть даст ли это совпадение картинок с ансис.
void solve_Structural_foundation_stiffness(TEMPER &t, WALL* &w, integer lw,
	bool bThermalStress, doublereal operatingtemperature,
	BLOCK* &b, integer &lb, integer &lu,
	bool btimedep, doublereal timestep_sizenow,
	doublereal* &uoldtimestep, doublereal* &uolddoubletimestep,
	doublereal poweron_multiplyer_sequence, TPROP* &matlist,
	doublereal* &t_for_Mechanical) {

	// btimedep==true - нестационарное моделирование,
	// btimedep==false - стационарная задача механики.
	// timestep_sizenow - размер шага по времени,
	// uoldtimestep - перемещения на один шаг назад,
	// uolddoubletimestep - перемещения на два шага назад.
	// uoldtimestep, uolddoubletimestep - память выделена заранее в вызывающем внешнем коде.
	// poweron_multiplyer_sequence==0.0 вектор силы выключен,
	// poweron_multiplyer_sequence==1.0 вектор силы полностью включён,
	// 0 < poweron_multiplyer_sequence < 1 - вектор силы частично активен.

	doublereal** deformation_old_iteration = new doublereal*[4];
	for (integer i_1 = 0; i_1 < 4; i_1++) {
		deformation_old_iteration[i_1] = new doublereal[t.maxelm + t.maxbound];
	}
	// initialization
	for (integer i_1 = 0; i_1 < 4; i_1++) {
#pragma omp parallel for
		for (integer j_1 = 0; j_1 < t.maxelm + t.maxbound; j_1++) {
			deformation_old_iteration[i_1][j_1] = 0.0;
		}
	}

	for (integer iter = 0; iter < 4; iter++) {

		solve_Structural(t, w, lw,
			bThermalStress, operatingtemperature,
			b, lb, lu,
			btimedep, timestep_sizenow,
			uoldtimestep, uolddoubletimestep,
			poweron_multiplyer_sequence, matlist,
			t_for_Mechanical);


		doublereal dmr = 0.0;
		for (integer i_1 = 0; i_1 < 4; i_1++) {
#pragma omp parallel for
			for (integer j_1 = 0; j_1 < t.maxelm + t.maxbound; j_1++) {
				deformation_old_iteration[i_1][j_1] = t.total_deformation[i_1][j_1];
			}
		}


		for (integer i_1 = 1; i_1 < 4; i_1++) {
#pragma omp parallel for
			for (integer j_1 = 0; j_1 < t.maxelm + t.maxbound; j_1++) {
				dmr += (deformation_old_iteration[i_1][j_1] - t.total_deformation[i_1][j_1])*(deformation_old_iteration[i_1][j_1] - t.total_deformation[i_1][j_1]);
			}
		}
		dmr /= 3.0*(t.maxelm + t.maxbound);
		std::cout << "residual foundation Mechanical=" << dmr << std::endl;

		// gamma_xy
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 2, 4, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 2, 4, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 1, 5, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 1, 5, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {

			t.total_deformation[STRAIN_XY][i_1] = t.total_deformation[4][i_1] + t.total_deformation[5][i_1];

		}

		// gamma_yz
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 3, 4, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 3, 4, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 2, 5, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 2, 5, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {

			t.total_deformation[STRAIN_YZ][i_1] = t.total_deformation[4][i_1] + t.total_deformation[5][i_1];

		}


		// gamma_zx
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 1, 4, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 1, 4, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, 3, 5, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, 3, 5, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {

			t.total_deformation[STRAIN_ZX][i_1] = t.total_deformation[4][i_1] + t.total_deformation[5][i_1];

		}

		// epsilon_x
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, XDEFORMATION, STRAIN_X, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, XDEFORMATION, STRAIN_X, LINE_DIRECTIONAL::X_LINE_DIRECTIONAL);
		}


		// epsilon_y
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, YDEFORMATION, STRAIN_Y, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, YDEFORMATION, STRAIN_Y, LINE_DIRECTIONAL::Y_LINE_DIRECTIONAL);
		}

		// epsilon_z
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только внутренние узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, false,
				t.border_neighbor, t.ilevel_alice, ZDEFORMATION, STRAIN_Z, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

#pragma omp parallel for
		for (integer i = 0; i < t.maxelm; i++) {
			// Только граничные узлы.
			green_gauss_Stress(i, t.total_deformation, t.nvtx, t.pa,
				t.neighbors_for_the_internal_node, t.maxelm, true,
				t.border_neighbor, t.ilevel_alice, ZDEFORMATION, STRAIN_Z, LINE_DIRECTIONAL::Z_LINE_DIRECTIONAL);
		}

		double **Dirichlet = new doublereal*[6];
		for (integer i_11 = 0; i_11 < 6; i_11++) {
			Dirichlet[i_11] = new doublereal[6];
		}
		for (integer i_1 = 0; i_1 < t.maxelm; i_1++) {


			doublereal beta_t_solid = t.prop[BETA_T_MECHANICAL][i_1]; // Коэффициенты Ламе, коэффициент линейного теплового расширения.
			doublereal beta_t_solid_x = t.prop[MULT_BETA_T_MECHANICAL_X][i_1] * t.prop[BETA_T_MECHANICAL][i_1];// Коэффициент линейного теплового расширения 1/K.
			doublereal beta_t_solid_y = t.prop[MULT_BETA_T_MECHANICAL_Y][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
			doublereal beta_t_solid_z = t.prop[MULT_BETA_T_MECHANICAL_Z][i_1] * t.prop[BETA_T_MECHANICAL][i_1];
			doublereal Ex = t.prop[MULT_YOUNG_MODULE_X][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
			doublereal Ey = t.prop[MULT_YOUNG_MODULE_Y][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
			doublereal Ez = t.prop[MULT_YOUNG_MODULE_Z][i_1] * t.prop[YOUNG_MODULE][i_1]; // Модуль Юнга Па.
			doublereal E = t.prop[YOUNG_MODULE][i_1];
			doublereal nuyz = t.prop[MULT_POISSON_RATIO_YZ][i_1] * t.prop[POISSON_RATIO][i_1];
			doublereal nuxz = t.prop[MULT_POISSON_RATIO_XZ][i_1] * t.prop[POISSON_RATIO][i_1];
			doublereal nuxy = t.prop[MULT_POISSON_RATIO_XY][i_1] * t.prop[POISSON_RATIO][i_1];
			doublereal nuzy = t.prop[MULT_POISSON_RATIO_ZY][i_1] * t.prop[POISSON_RATIO][i_1];
			doublereal nuzx = t.prop[MULT_POISSON_RATIO_ZX][i_1] * t.prop[POISSON_RATIO][i_1];
			doublereal nuyx = t.prop[MULT_POISSON_RATIO_YX][i_1] * t.prop[POISSON_RATIO][i_1];

			doublereal nu = t.prop[POISSON_RATIO][i_1];


			doublereal Gxy, Gyz, Gxz;
			if (!t.bActiveShearModule[i_1]) {
				Gxy = Gyz = Gxz = Ex / (2.0 * (1.0 + nuxy));
			}
			else {
				Gyz = t.prop[SHEAR_MODULE_YZ][i_1];
				Gxz = t.prop[SHEAR_MODULE_XZ][i_1];
				Gxy = t.prop[SHEAR_MODULE_XY][i_1];
			}


			Dirichlet[0][0] = (nuyz*nuzy - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz
				* nuzx + nuyz * nuzy - 1.0)*Ex;
			Dirichlet[0][1] = -(nuxz*nuzy + nuxy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ey;
			Dirichlet[0][2] = -(nuxy*nuyz + nuxz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ez;
			Dirichlet[0][3] = 0.0;
			Dirichlet[0][4] = 0.0;
			Dirichlet[0][5] = 0.0;
			Dirichlet[1][0] = -(nuyz*nuzx + nuyx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ex;
			Dirichlet[1][1] = (nuxz*nuzx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz
				* nuzx + nuyz * nuzy - 1.0)*Ey;
			Dirichlet[1][2] = -(nuxz*nuyx + nuyz) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ez;
			Dirichlet[1][3] = 0.0;
			Dirichlet[1][4] = 0.0;
			Dirichlet[1][5] = 0.0;
			Dirichlet[2][0] = -(nuyx*nuzy + nuzx) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ex;
			Dirichlet[2][1] = -(nuxy*nuzx + nuzy) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx +
				nuxz * nuzx + nuyz * nuzy - 1.0)*Ey;
			Dirichlet[2][2] = (nuxy*nuyx - 1.0) / (nuxy*nuyz*nuzx + nuxz * nuyx*nuzy + nuxy * nuyx + nuxz
				* nuzx + nuyz * nuzy - 1.0)*Ez;
			Dirichlet[2][3] = 0.0;
			Dirichlet[2][4] = 0.0;
			Dirichlet[2][5] = 0.0;
			Dirichlet[3][0] = 0.0;
			Dirichlet[3][1] = 0.0;
			Dirichlet[3][2] = 0.0;
			Dirichlet[3][3] = Gxy;
			Dirichlet[3][4] = 0.0;
			Dirichlet[3][5] = 0.0;
			Dirichlet[4][0] = 0.0;
			Dirichlet[4][1] = 0.0;
			Dirichlet[4][2] = 0.0;
			Dirichlet[4][3] = 0.0;
			Dirichlet[4][4] = Gyz;
			Dirichlet[4][5] = 0.0;
			Dirichlet[5][0] = 0.0;
			Dirichlet[5][1] = 0.0;
			Dirichlet[5][2] = 0.0;
			Dirichlet[5][3] = 0.0;
			Dirichlet[5][4] = 0.0;
			Dirichlet[5][5] = Gxz;


			t.total_deformation[STRESS_X][i_1] = 0.0;
			for (integer i_11 = 0; i_11 < 6; i_11++) {
				if (((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))) {
					if ((i_11 == 0)) {
						t.total_deformation[STRESS_X][i_1] += Dirichlet[0][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] -
							beta_t_solid_x * (t.potent[i_1] - t.operatingtemperature_copy));
					}
					if ((i_11 == 1)) {
						t.total_deformation[STRESS_X][i_1] += Dirichlet[0][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] -
							beta_t_solid_y * (t.potent[i_1] - t.operatingtemperature_copy));
					}
					if ((i_11 == 2)) {
						t.total_deformation[STRESS_X][i_1] += Dirichlet[0][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] -
							beta_t_solid_z * (t.potent[i_1] - t.operatingtemperature_copy));
					}
				}
				else {
					t.total_deformation[STRESS_X][i_1] += Dirichlet[0][i_11] * t.total_deformation[i_11 + STRAIN_X][i_1];
				}
			}
			t.total_deformation[STRESS_Y][i_1] = 0.0;
			for (integer i_11 = 0; i_11 < 6; i_11++) {
				if (((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))) {
					if ((i_11 == 0)) {
						t.total_deformation[STRESS_Y][i_1] += Dirichlet[1][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] -
							beta_t_solid_x * (t.potent[i_1] - t.operatingtemperature_copy));
					}
					if ((i_11 == 1)) {
						t.total_deformation[STRESS_Y][i_1] += Dirichlet[1][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] -
							beta_t_solid_y * (t.potent[i_1] - t.operatingtemperature_copy));
					}
					if ((i_11 == 2)) {
						t.total_deformation[STRESS_Y][i_1] += Dirichlet[1][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] -
							beta_t_solid_z * (t.potent[i_1] - t.operatingtemperature_copy));
					}
				}
				else {
					t.total_deformation[STRESS_Y][i_1] += Dirichlet[1][i_11] * t.total_deformation[i_11 + STRAIN_X][i_1];
				}
			}
			t.total_deformation[STRESS_Z][i_1] = 0.0;
			for (integer i_11 = 0; i_11 < 6; i_11++) {
				if (((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
					(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))) {
					if ((i_11 == 0)) {
						t.total_deformation[STRESS_Z][i_1] += Dirichlet[2][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] -
							beta_t_solid_x * (t.potent[i_1] - t.operatingtemperature_copy));
					}
					if ((i_11 == 1)) {
						t.total_deformation[STRESS_Z][i_1] += Dirichlet[2][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] -
							beta_t_solid_y * (t.potent[i_1] - t.operatingtemperature_copy));
					}
					if ((i_11 == 2)) {
						t.total_deformation[STRESS_Z][i_1] += Dirichlet[2][i_11] * (t.total_deformation[i_11 + STRAIN_X][i_1] -
							beta_t_solid_z * (t.potent[i_1] - t.operatingtemperature_copy));
					}
				}
				else {
					t.total_deformation[STRESS_Z][i_1] += Dirichlet[2][i_11] * t.total_deformation[i_11 + STRAIN_X][i_1];
				}
			}

			t.total_deformation[STRESS_XY][i_1] = Dirichlet[3][3] * t.total_deformation[STRAIN_XY][i_1];
			t.total_deformation[STRESS_YZ][i_1] = Dirichlet[4][4] * t.total_deformation[STRAIN_YZ][i_1];
			t.total_deformation[STRESS_ZX][i_1] = Dirichlet[5][5] * t.total_deformation[STRAIN_ZX][i_1];
		}
		for (integer i_11 = 0; i_11 < 6; i_11++) {
			delete[] Dirichlet[i_11];
		}
		delete[] Dirichlet;

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxbound; i_1++) {
			t.total_deformation[STRESS_X][i_1] = t.total_deformation[STRESS_X][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_Y][i_1] = t.total_deformation[STRESS_Y][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_Z][i_1] = t.total_deformation[STRESS_Z][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_XY][i_1] = t.total_deformation[STRESS_XY][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_YZ][i_1] = t.total_deformation[STRESS_YZ][t.border_neighbor[i_1].iI];
			t.total_deformation[STRESS_ZX][i_1] = t.total_deformation[STRESS_ZX][t.border_neighbor[i_1].iI];
		}

		// epsilon (STRAIN) von Mizes
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < t.maxelm + t.maxbound; i_1++) {

			t.total_deformation[STRAIN_VON_MIZES][i_1] = sqrt(0.5*((t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Y][i_1])*
				(t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Y][i_1]) + (t.total_deformation[STRAIN_Y][i_1] - t.total_deformation[STRAIN_Z][i_1]) *
				(t.total_deformation[STRAIN_Y][i_1] - t.total_deformation[STRAIN_Z][i_1]) + (t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Z][i_1]) *
				(t.total_deformation[STRAIN_X][i_1] - t.total_deformation[STRAIN_Z][i_1])));

			t.total_deformation[LOG10_STRAIN_VON_MIZES][i_1] = log10(t.total_deformation[STRAIN_VON_MIZES][i_1]);


			// STRESS

			t.total_deformation[STRESS_VON_MIZES][i_1] = sqrt(0.5*((t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Y][i_1])*
				(t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Y][i_1]) + (t.total_deformation[STRESS_Y][i_1] - t.total_deformation[STRESS_Z][i_1]) *
				(t.total_deformation[STRESS_Y][i_1] - t.total_deformation[STRESS_Z][i_1]) + (t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Z][i_1]) *
				(t.total_deformation[STRESS_X][i_1] - t.total_deformation[STRESS_Z][i_1])));

			t.total_deformation[LOG10_STRESS_VON_MIZES][i_1] = log10(t.total_deformation[STRESS_VON_MIZES][i_1]);

		}

	}

}


// Предупреждает в случае нарушения физики (консервативности).
// Посыл: положительная мощность приводит только к росту температуры,
// иначе нарушена консервативность.
void debug_signal(TEMPER& t, doublereal operating_temperature);

// решает одно уравнение
// например уравнение теплопроводности.
// res - возвращаемая невязка.
void solve(integer iVar, doublereal &res, FLOW &f,
	FLOW* &fglobal, TEMPER &t,
	doublereal** &rhie_chow,
	SOURCE* &s, WALL* &w, BLOCK* &b,
	integer ls, integer lw, integer lb,
	doublereal dbeta,
	integer flow_interior,
	bool bconvective,
	bool bfirst_start,
	doublereal* toldtimestep, // рассчитываемое поле температур с предыдущего шага по времени,
	doublereal* told_iter, // Температура с предыдущей итерации.
	doublereal** speedoldtimestep, // рассчитываемое поле скоростей с предыдущего шага по времени,
	doublereal** mfoldtimestep, // конвективный поток через грани контрольного объёма с предыдущего временного слоя, 
	doublereal tauparam,  // размер шага по времени
	bool btimedep, // стационарный или нестационарный солвер.
	doublereal dgx, doublereal dgy, doublereal dgz,
	TPROP* &matlist,
	integer inumiter,// номер итерации SIMPLE алгоритма
	bool bprintmessage, doublereal RCh,
	bool bVERYStable,
	doublereal** tau,
	doublereal** &sumanb,
	bool bmyhighorder, bool bdeltapfinish,
	doublereal poweron_multiplier_sequence,
	QuickMemVorst& m, doublereal* &rthdsd,
	doublereal &rfluentresval, 
	integer lu, UNION* &my_union,
	integer* &color, integer dist_max)
{

	

	// QuickMemVorst& m - дополнительная память для алгоритма Хенка ван дер Ворста, чтобы избежать частых выделений и уничтожений памяти.

	// btimedep передаётся по значению его можно менять внутри данной функции и это никак не отразится на этом параметре за пределами
	// данной функции.
	// если bfirst_start==true то мы имеем дело с первой итерацией алгоритма SIMPLE и этот случай должен быть обработан по особому.

	// сборка матрицы СЛАУ
	//integer i = 0; // счётчик цикла for
	//doublereal RCh=0.1;//0.1; 1.0;

	bool bRhieChowiPAM = true;
	/*
	// Оказалось что данная вставка не влияет на сходимость
	// поэтому от неё и отказались.
	// Отключение или включение поправки Рхи-Чоу на первых итерациях
	// не влияет на процесс сходимости.
	if (inumiter<2) {
		 bRhieChowiPAM=false;
	}
	*/
	integer i75 = 0;

	bool breversedflow = false;
	integer icell = 0; // количество ячеек с возвратно-циркуляционным течением.
	integer imyscheme = UDS;

#ifdef _OPENMP 
	// Узнаёт количество ядер в системе.
	// 15млн неизвестных время параллельного кода 21мин 39с.
	// Время однопоточного кода 27мин 48с.
	unsigned int nthreads = number_cores();
	omp_set_num_threads(nthreads); // установка числа потоков
#endif

	bool brthdsd_ON1 = true;
	if (inumiter < 7) brthdsd_ON1 = false;

	bool brthdsd_ON = true;
	if (inumiter < 7) brthdsd_ON = false;

	switch (iVar) {
	case PAM:

		// Симметризация СЛАУ выполняется внутри
		// сборки матрицы:
		
#ifdef _OPENMP

		if (bparallelizm_old) {
			printf("error bparallelizm_old\n");
			system("pause");

			if (inumcore == 1) {
				// serial
				// Граничные условия Дирихле обязательно 
				// должны собираться в первую очередь
				for (integer  i = 0; i < f.maxbound; i++) {

					breversedflow = false;
					// Заполнение граничных условий.
					// В уравнении для поправки давления.
					// Сначала собираем только условия Дирихле
					// последний параметр bool bDirichlet равен true.
					/*
					my_elmatr_quad_PAm_bon( f.slau_bon,
					f.slau, i,
					f.maxelm,
					f.maxbound,
					f.border_neighbor,
					f.nvtx,
					f.bPressureFix,
					dbeta, f.pa,
					f.potent,
					f.prop, f.prop_b,
					f.alpha,
					ls, lw, w,
					true, f.neighbors_for_the_internal_node,
					f.diag_coef, RCh,
					breversedflow);
					//*/

					// Граничные условия на основе сглаженного коэффициент tau.
					// Сначала собираем только условия Дирихле
					// последний параметр bool bDirichlet равен true.
					//*
					my_elmatr_quad_PAm_bon3(f.slau_bon,
						f.slau, i,
						f.maxelm,
						f.maxbound,
						f.border_neighbor,
						f.nvtx,
						f.bPressureFix,
						dbeta, f.pa,
						f.potent,
						f.prop, f.prop_b,
						f.alpha,
						ls, lw, w,
						true,
						f.neighbors_for_the_internal_node, RCh,
						breversedflow,
						tau);//*/

					if (breversedflow) icell++;
				}

			}

			if (inumcore == 2) {
				if (nd.b0.active) {

					// первый поток
					for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
						integer iPloc = f.ifrontregulationgl[iscan_par];
						if (iPloc >= f.maxelm) {
							breversedflow = false;
							// Заполнение граничных условий.
							// В уравнении для поправки давления.
							// Сначала собираем только условия Дирихле
							// последний параметр bool bDirichlet равен true.
							/*
							my_elmatr_quad_PAm_bon( f.slau_bon,
											 f.slau, iPloc-f.maxelm,
											 f.maxelm,
											 f.maxbound,
											 f.border_neighbor,
											 f.nvtx,
											 f.bPressureFix,
											 dbeta, f.pa,
											 f.potent,
											 f.prop, f.prop_b,
											 f.alpha,
											 ls, lw, w,
											 true, f.neighbors_for_the_internal_node,
											 f.diag_coef, RCh,
											 breversedflow);
											 //*/

											 // Граничные условия на основе сглаженного коэффициент tau.
											 // Сначала собираем только условия Дирихле
											 // последний параметр bool bDirichlet равен true.
											 //*
							my_elmatr_quad_PAm_bon3(f.slau_bon,
								f.slau, iPloc - f.maxelm,
								f.maxelm,
								f.maxbound,
								f.border_neighbor,
								f.nvtx,
								f.bPressureFix,
								dbeta, f.pa,
								f.potent,
								f.prop, f.prop_b,
								f.alpha,
								ls, lw, w,
								true,
								f.neighbors_for_the_internal_node, RCh,
								breversedflow,
								tau);//*/

							if (breversedflow) icell++;
						}
					}
					// второй поток
					for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
						integer iPloc = f.ifrontregulationgl[iscan_par];
						if (iPloc >= f.maxelm) {
							breversedflow = false;
							// Заполнение граничных условий.
							// В уравнении для поправки давления.
							// Сначала собираем только условия Дирихле
							// последний параметр bool bDirichlet равен true.
							/*
							my_elmatr_quad_PAm_bon( f.slau_bon,
											f.slau, iPloc-f.maxelm,
											f.maxelm,
											f.maxbound,
											f.border_neighbor,
											f.nvtx,
											f.bPressureFix,
											dbeta, f.pa,
											f.potent,
											f.prop, f.prop_b,
											f.alpha,
											ls, lw, w,
											true, f.neighbors_for_the_internal_node,
											f.diag_coef, RCh,
											breversedflow);
											//*/

											// Граничные условия на основе сглаженного коэффициент tau.
											// Сначала собираем только условия Дирихле
											// последний параметр bool bDirichlet равен true.
											//*
							my_elmatr_quad_PAm_bon3(f.slau_bon,
								f.slau, iPloc - f.maxelm,
								f.maxelm,
								f.maxbound,
								f.border_neighbor,
								f.nvtx,
								f.bPressureFix,
								dbeta, f.pa,
								f.potent,
								f.prop, f.prop_b,
								f.alpha,
								ls, lw, w,
								true,
								f.neighbors_for_the_internal_node, RCh,
								breversedflow,
								tau);//*/

							if (breversedflow) icell++;
						}
					}
					// серийный смыкающий кусок
					for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc = f.ifrontregulationgl[iscan_par];
						if (iPloc >= f.maxelm) {
							breversedflow = false;
							// Заполнение граничных условий.
							// В уравнении для поправки давления.
							// Сначала собираем только условия Дирихле
							// последний параметр bool bDirichlet равен true.
							/*
							my_elmatr_quad_PAm_bon( f.slau_bon,
											 f.slau, iPloc-f.maxelm,
											 f.maxelm,
											 f.maxbound,
											 f.border_neighbor,
											 f.nvtx,
											 f.bPressureFix,
											 dbeta, f.pa,
											 f.potent,
											 f.prop, f.prop_b,
											 f.alpha,
											 ls, lw, w,
											 true, f.neighbors_for_the_internal_node,
											 f.diag_coef, RCh,
											 breversedflow);
											 //*/

											 // Граничные условия на основе сглаженного коэффициент tau.
											 // Сначала собираем только условия Дирихле
											 // последний параметр bool bDirichlet равен true.
											 //*
							my_elmatr_quad_PAm_bon3(f.slau_bon,
								f.slau, iPloc - f.maxelm,
								f.maxelm,
								f.maxbound,
								f.border_neighbor,
								f.nvtx,
								f.bPressureFix,
								dbeta, f.pa,
								f.potent,
								f.prop, f.prop_b,
								f.alpha,
								ls, lw, w,
								true,
								f.neighbors_for_the_internal_node, RCh,
								breversedflow,
								tau);//*/

							if (breversedflow) icell++;
						}
					}


				}
			}
		}
		else {
			// Граничные условия Дирихле обязательно 
			// должны собираться в первую очередь
			integer icell_loc = 0;

#pragma omp parallel for reduction(+:icell_loc)
			for (integer  i = 0; i < f.maxbound; i++) {

				breversedflow = false;
				// Заполнение граничных условий.
				// В уравнении для поправки давления.
				// Сначала собираем только условия Дирихле
				// последний параметр bool bDirichlet равен true.
				/*
				my_elmatr_quad_PAm_bon( f.slau_bon,
				f.slau, i,
				f.maxelm,
				f.maxbound,
				f.border_neighbor,
				f.nvtx,
				f.bPressureFix,
				dbeta, f.pa,
				f.potent,
				f.prop, f.prop_b,
				f.alpha,
				ls, lw, w,
				true, f.neighbors_for_the_internal_node,
				f.diag_coef, RCh,
				breversedflow);
				//*/

				// Граничные условия на основе сглаженного коэффициент tau.
				// Сначала собираем только условия Дирихле
				// последний параметр bool bDirichlet равен true.
				//*
				my_elmatr_quad_PAm_bon3(f.slau_bon,
					f.slau, i,
					f.maxelm,
					f.maxbound,
					f.border_neighbor,
					f.nvtx,
					f.bPressureFix,
					dbeta, f.pa,
					f.potent,
					f.prop, f.prop_b,
					f.alpha,
					ls, lw, w,
					true,
					f.neighbors_for_the_internal_node, RCh,
					breversedflow,
					tau);//*/

				if (breversedflow) icell_loc++;
			}
		
			icell += icell_loc;
}

#else

				  // Граничные условия Дирихле обязательно 
				  // должны собираться в первую очередь
		for (integer  i = 0; i < f.maxbound; i++) {

			breversedflow = false;
			// Заполнение граничных условий.
			// В уравнении для поправки давления.
			// Сначала собираем только условия Дирихле
			// последний параметр bool bDirichlet равен true.
			/*
			my_elmatr_quad_PAm_bon( f.slau_bon,
									f.slau, i,
									f.maxelm,
									f.maxbound,
									f.border_neighbor,
									f.nvtx,
									f.bPressureFix,
									dbeta, f.pa,
									f.potent,
									f.prop, f.prop_b,
									f.alpha,
									ls, lw, w,
									true, f.neighbors_for_the_internal_node,
									f.diag_coef, RCh,
									breversedflow);
									//*/

									// Граничные условия на основе сглаженного коэффициент tau.
									 // Сначала собираем только условия Дирихле
									// последний параметр bool bDirichlet равен true.
									//*
			my_elmatr_quad_PAm_bon3(f.slau_bon,
				f.slau, i,
				f.maxelm,
				f.maxbound,
				f.border_neighbor,
				f.nvtx,
				f.bPressureFix,
				dbeta, f.pa,
				f.potent,
				f.prop, f.prop_b,
				f.alpha,
				ls, lw, w,
				true,
				f.neighbors_for_the_internal_node, RCh,
				breversedflow,
				tau);//*/

			if (breversedflow) icell++;
		}

#endif

		// printf("step Dirichlet...\n");
		// getchar();


#ifdef _OPENMP

		if (bparallelizm_old) {
			printf("error bparallelizm_old\n");
			system("pause");

			if (inumcore == 1) {
				// serial
				for (integer  i = 0; i < f.maxbound; i++) {


					breversedflow = false;
					// Заполнение граничных условий.
					// В уравнении для поправки давления.
					// Метод отложенной коррекции.
					// Собираем однородные условия Неймана.
					// последний параметр bool bDirichlet равен false.
					/*
					my_elmatr_quad_PAm_bon( f.slau_bon,
					f.slau, i,
					f.maxelm,
					f.maxbound,
					f.border_neighbor,
					f.nvtx,
					f.bPressureFix,
					dbeta, f.pa,
					f.potent, f.prop,
					f.prop_b, f.alpha,
					ls, lw, w, false,
					f.neighbors_for_the_internal_node,
					f.diag_coef, RCh,
					breversedflow);
					//					   */

					// Граничные условия на основе сглаженного коэффициент tau.
					// Метод отложенной коррекции.
					// Собираем однородные условия Неймана.
					// последний параметр bool bDirichlet равен false.
					//*
					my_elmatr_quad_PAm_bon3(f.slau_bon,
						f.slau, i,
						f.maxelm,
						f.maxbound,
						f.border_neighbor,
						f.nvtx,
						f.bPressureFix,
						dbeta, f.pa,
						f.potent,
						f.prop, f.prop_b,
						f.alpha,
						ls, lw, w,
						false,
						f.neighbors_for_the_internal_node, RCh,
						breversedflow,
						tau);

					//						  */
					if (breversedflow) icell++;


				}
			}


			if (inumcore == 2) {
				if (nd.b0.active) {

					// первый поток
					for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
						integer iPloc = f.ifrontregulationgl[iscan_par];
						if (iPloc >= f.maxelm) {
							breversedflow = false;
							// Заполнение граничных условий.
							// В уравнении для поправки давления.
							// Метод отложенной коррекции.
							// Собираем однородные условия Неймана.
							// последний параметр bool bDirichlet равен false.
							/*
							my_elmatr_quad_PAm_bon( f.slau_bon,
											  f.slau, iPloc-f.maxelm,
											  f.maxelm,
											  f.maxbound,
											  f.border_neighbor,
											  f.nvtx,
											  f.bPressureFix,
											  dbeta, f.pa,
											  f.potent, f.prop,
											  f.prop_b, f.alpha,
											  ls, lw, w, false,
											  f.neighbors_for_the_internal_node,
											  f.diag_coef, RCh,
											  breversedflow);
							 //					   */

							 // Граничные условия на основе сглаженного коэффициент tau.
							 // Метод отложенной коррекции.
							 // Собираем однородные условия Неймана.
							  // последний параметр bool bDirichlet равен false.
							  //*
							my_elmatr_quad_PAm_bon3(f.slau_bon,
								f.slau, iPloc - f.maxelm,
								f.maxelm,
								f.maxbound,
								f.border_neighbor,
								f.nvtx,
								f.bPressureFix,
								dbeta, f.pa,
								f.potent,
								f.prop, f.prop_b,
								f.alpha,
								ls, lw, w,
								false,
								f.neighbors_for_the_internal_node, RCh,
								breversedflow,
								tau);

							//						  */
							if (breversedflow) icell++;


						}

					}
					// второй поток
					for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
						integer iPloc = f.ifrontregulationgl[iscan_par];
						if (iPloc >= f.maxelm) {
							breversedflow = false;
							// Заполнение граничных условий.
							// В уравнении для поправки давления.
							// Метод отложенной коррекции.
							// Собираем однородные условия Неймана.
							// последний параметр bool bDirichlet равен false.
							/*
							my_elmatr_quad_PAm_bon( f.slau_bon,
											  f.slau, iPloc-f.maxelm,
											  f.maxelm,
											  f.maxbound,
											  f.border_neighbor,
											  f.nvtx,
											  f.bPressureFix,
											  dbeta, f.pa,
											  f.potent, f.prop,
											  f.prop_b, f.alpha,
											  ls, lw, w, false,
											  f.neighbors_for_the_internal_node,
											  f.diag_coef, RCh,
											  breversedflow);
							//					   */

							// Граничные условия на основе сглаженного коэффициент tau.
							// Метод отложенной коррекции.
							// Собираем однородные условия Неймана.
							// последний параметр bool bDirichlet равен false.
							//*
							my_elmatr_quad_PAm_bon3(f.slau_bon,
								f.slau, iPloc - f.maxelm,
								f.maxelm,
								f.maxbound,
								f.border_neighbor,
								f.nvtx,
								f.bPressureFix,
								dbeta, f.pa,
								f.potent,
								f.prop, f.prop_b,
								f.alpha,
								ls, lw, w,
								false,
								f.neighbors_for_the_internal_node, RCh,
								breversedflow,
								tau);

							//						  */
							if (breversedflow) icell++;


						}
					}
					// серийный смыкающий кусок
					for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc = f.ifrontregulationgl[iscan_par];
						if (iPloc >= f.maxelm) {
							breversedflow = false;
							// Заполнение граничных условий.
							// В уравнении для поправки давления.
							// Метод отложенной коррекции.
							// Собираем однородные условия Неймана.
							// последний параметр bool bDirichlet равен false.
							/*
							my_elmatr_quad_PAm_bon( f.slau_bon,
											  f.slau, iPloc-f.maxelm,
											  f.maxelm,
											  f.maxbound,
											  f.border_neighbor,
											  f.nvtx,
											  f.bPressureFix,
											  dbeta, f.pa,
											  f.potent, f.prop,
											  f.prop_b, f.alpha,
											  ls, lw, w, false,
											  f.neighbors_for_the_internal_node,
											  f.diag_coef, RCh,
											  breversedflow);
							//					   */

							// Граничные условия на основе сглаженного коэффициент tau.
							// Метод отложенной коррекции.
							// Собираем однородные условия Неймана.
							// последний параметр bool bDirichlet равен false.
							//*
							my_elmatr_quad_PAm_bon3(f.slau_bon,
								f.slau, iPloc - f.maxelm,
								f.maxelm,
								f.maxbound,
								f.border_neighbor,
								f.nvtx,
								f.bPressureFix,
								dbeta, f.pa,
								f.potent,
								f.prop, f.prop_b,
								f.alpha,
								ls, lw, w,
								false,
								f.neighbors_for_the_internal_node, RCh,
								breversedflow,
								tau);

							//						  */
							if (breversedflow) icell++;


						}
					}


				}
	}
}
else {

integer icell_loc = 0;

#pragma omp parallel for reduction(+:icell_loc)
	for (integer  i = 0; i<f.maxbound; i++) {


		breversedflow = false;
		// Заполнение граничных условий.
		// В уравнении для поправки давления.
		// Метод отложенной коррекции.
		// Собираем однородные условия Неймана.
		// последний параметр bool bDirichlet равен false.
		/*
		my_elmatr_quad_PAm_bon( f.slau_bon,
		f.slau, i,
		f.maxelm,
		f.maxbound,
		f.border_neighbor,
		f.nvtx,
		f.bPressureFix,
		dbeta, f.pa,
		f.potent, f.prop,
		f.prop_b, f.alpha,
		ls, lw, w, false,
		f.neighbors_for_the_internal_node,
		f.diag_coef, RCh,
		breversedflow);
		//					   */

		// Граничные условия на основе сглаженного коэффициент tau.
		// Метод отложенной коррекции.
		// Собираем однородные условия Неймана.
		// последний параметр bool bDirichlet равен false.
		//*
		my_elmatr_quad_PAm_bon3(f.slau_bon,
			f.slau, i,
			f.maxelm,
			f.maxbound,
			f.border_neighbor,
			f.nvtx,
			f.bPressureFix,
			dbeta, f.pa,
			f.potent,
			f.prop, f.prop_b,
			f.alpha,
			ls, lw, w,
			false,
			f.neighbors_for_the_internal_node, RCh,
			breversedflow,
			tau);

		//						  */
		if (breversedflow) icell_loc++;


				  }


	icell += icell_loc;
}

#else

			      for (integer  i=0; i<f.maxbound; i++) {


					  breversedflow=false;
                       // Заполнение граничных условий.
                       // В уравнении для поправки давления.
                       // Метод отложенной коррекции.
					   // Собираем однородные условия Неймана.
                       // последний параметр bool bDirichlet равен false.
					  /*
					   my_elmatr_quad_PAm_bon( f.slau_bon,
						                       f.slau, i,
											   f.maxelm, 
											   f.maxbound,
											   f.border_neighbor,
											   f.nvtx,
											   f.bPressureFix,
											   dbeta, f.pa,
											   f.potent, f.prop,
											   f.prop_b, f.alpha,
											   ls, lw, w, false, 
											   f.neighbors_for_the_internal_node, 
											   f.diag_coef, RCh,
											   breversedflow);
						//					   */

					   // Граничные условия на основе сглаженного коэффициент tau.
					  // Метод отложенной коррекции.
					   // Собираем однородные условия Неймана.
                       // последний параметр bool bDirichlet равен false.
					//*
					my_elmatr_quad_PAm_bon3(f.slau_bon,
						                      f.slau, i, 
											  f.maxelm,
											  f.maxbound, 
							                  f.border_neighbor,
											  f.nvtx,
											  f.bPressureFix,
							                  dbeta, f.pa,
											  f.potent,
							                  f.prop, f.prop_b,
											  f.alpha,
							                  ls, lw, w, 
											  false, 
							                  f.neighbors_for_the_internal_node, RCh,
							                  breversedflow, 
											  tau);

					//						  */
					   if (breversedflow) icell++;

			      
				  }
#endif

				 // printf("step Neiman...\n");
				 // getchar();
#pragma omp parallel for
				  for (integer  i = 0; i < f.maxbound; i++) {
					  if (f.slau_bon[PAM][i].aw < 0.0) {
						  printf("maxbound=%lld i=%lld problem PAM aw=%e\n", f.maxbound, i, f.slau_bon[PAM][i].aw);
						  system("pause");
						  //f.slau_bon[iVar][i].aw
					  }
				  }

				  /*
				  // debug
				   if (inumiter==1) {
					   for ( i=t.maxbound-2; i<t.maxbound; i++) {
					   #if doubleintprecision == 1
							 printf("id=%lld aw=%e ai=%e b=%e\n", i+f.maxelm, f.slau_bon[PAM][i].aw, f.slau_bon[PAM][i].ai, f.slau_bon[PAM][i].b);
					   #else
							printf("id=%d aw=%e ai=%e b=%e\n", i+f.maxelm, f.slau_bon[PAM][i].aw, f.slau_bon[PAM][i].ai, f.slau_bon[PAM][i].b);
					   #endif
		         	getchar();
					}
		          }
				  */

				  if (icell>0) {
					  // уменьшение icell по идее должно говорить о сходимости,
					  // увеличение icell должно говорить о расходимости.
					  // стабильное значение icell>0 говорит о сошедшемся решении через 
					  // выходную границу которого проходит возвратно циркуляционное течение.
#if doubleintprecision == 1
					  printf("reversed flow in %lld cell pressure outlet...\n", icell);
#else
					  printf("reversed flow in %d cell pressure outlet...\n", icell);
#endif

					  
				  }
				  

				#ifdef _OPENMP

				  if (bparallelizm_old) {
					  printf("error bparallelizm_old\n");
					  system("pause");
					  if (inumcore == 1) {

						  for (integer  i = 0; i < f.maxelm; i++) {
							  rhie_chow[0][i] = 0.0;
							  rhie_chow[1][i] = 0.0;
							  rhie_chow[2][i] = 0.0;

							  if ((mfoldtimestep == nullptr) && (speedoldtimestep == nullptr)) {

								  // стационарный солвер.
								  // будем использовать отрицательный шаг по времени,
								  // это будет сигнал использовать среднее арифметическое за опорное время.
								  /*
								  my_elmatr_quad_PAm(i, f.slau,
								  f.slau_bon,
								  f.potent,
								  f.pa, f.prop,
								  f.prop_b, f.nvtx,
								  f.neighbors_for_the_internal_node, f.maxelm,
								  f.diag_coef,
								  f.alpha, dbeta,
								  rhie_chow, RCh,
								  btimedep,
								  tauparam,
								  nullptr,
								  f.mf[i],
								  nullptr);
								  */

								  // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
								  // см. опыт Гаврилова Анрея sigma cfd.
								  //*
								  my_elmatr_quad_PAm3(i, f.slau,
									  f.slau_bon,
									  f.potent,
									  f.pa,
									  f.prop,
									  f.prop_b,
									  f.nvtx,
									  f.neighbors_for_the_internal_node,
									  f.maxelm,
									  f.alpha,
									  dbeta,
									  rhie_chow, RCh,
									  btimedep,
									  -0.1, // tauparam
									  nullptr,
									  f.mf[i],
									  nullptr,
									  tau,
									  bmyhighorder, bdeltapfinish,
									  bRhieChowiPAM, false,
									  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
									  sumanb);//*/

							  }
							  else
							  {
								  // поправка давления
								  /*
								  my_elmatr_quad_PAm(i, f.slau,
								  f.slau_bon,
								  f.potent,
								  f.pa, f.prop,
								  f.prop_b, f.nvtx,
								  f.neighbors_for_the_internal_node, f.maxelm,
								  f.diag_coef,
								  f.alpha, dbeta,
								  rhie_chow, RCh,
								  btimedep,
								  tauparam,
								  mfoldtimestep[i],
								  f.mf[i],
								  speedoldtimestep);
								  */

								  // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
								  // см. опыт Гаврилова Анрея sigma cfd.
								  ///*
								  my_elmatr_quad_PAm3(i, f.slau,
									  f.slau_bon,
									  f.potent,
									  f.pa,
									  f.prop,
									  f.prop_b,
									  f.nvtx,
									  f.neighbors_for_the_internal_node,
									  f.maxelm,
									  f.alpha, dbeta,
									  rhie_chow, RCh,
									  btimedep,
									  tauparam,
									  mfoldtimestep[i],
									  f.mf[i],
									  speedoldtimestep,
									  tau,
									  bmyhighorder, bdeltapfinish,
									  bRhieChowiPAM, false, 
									  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
									  sumanb);
								  //*/
							  }
						  }
					  }

					  if (inumcore == 2) {
						  if ((mfoldtimestep == nullptr) && (speedoldtimestep == nullptr)) {

							  if (nd.b0.active) {

#pragma omp parallel
								  {
#pragma omp sections 
									  {
#pragma omp section
										  {
											  // первый поток
											  for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
												  integer iPloc = f.ifrontregulationgl[iscan_par];
												  if (iPloc < f.maxelm) {
													  rhie_chow[0][iPloc] = 0.0;
													  rhie_chow[1][iPloc] = 0.0;
													  rhie_chow[2][iPloc] = 0.0;



													  // стационарный солвер.
													  // будем использовать отрицательный шаг по времени,
													  // это будет сигнал использовать среднее арифметическое за опорное время.
													 /*
													  my_elmatr_quad_PAm(iPloc, f.slau,
																	 f.slau_bon,
																	 f.potent,
																	 f.pa, f.prop,
																	 f.prop_b, f.nvtx,
																	 f.neighbors_for_the_internal_node, f.maxelm,
																	 f.diag_coef,
																	 f.alpha, dbeta,
																	 rhie_chow, RCh,
																	 btimedep,
																	 tauparam,
																	 nullptr,
																	 f.mf[iPloc],
																	 nullptr);
													 */

													 // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
													 // см. опыт Гаврилова Анрея sigma cfd.
													 //*
													  my_elmatr_quad_PAm3(iPloc, f.slau,
														  f.slau_bon,
														  f.potent,
														  f.pa,
														  f.prop,
														  f.prop_b,
														  f.nvtx,
														  f.neighbors_for_the_internal_node,
														  f.maxelm,
														  f.alpha,
														  dbeta,
														  rhie_chow, RCh,
														  btimedep,
														  -0.1, // tauparam
														  nullptr,
														  f.mf[iPloc],
														  nullptr,
														  tau,
														  bmyhighorder, bdeltapfinish,
														  bRhieChowiPAM, false, 
														  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
														  sumanb);//*/


												  }

											  }
										  }
#pragma omp section
										  {
											  // второй поток
											  for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
												  integer iPloc = f.ifrontregulationgl[iscan_par];
												  if (iPloc < f.maxelm) {
													  rhie_chow[0][iPloc] = 0.0;
													  rhie_chow[1][iPloc] = 0.0;
													  rhie_chow[2][iPloc] = 0.0;



													  // стационарный солвер.
													  // будем использовать отрицательный шаг по времени,
													  // это будет сигнал использовать среднее арифметическое за опорное время.
													 /*
													  my_elmatr_quad_PAm(iPloc, f.slau,
																	 f.slau_bon,
																	 f.potent,
																	 f.pa, f.prop,
																	 f.prop_b, f.nvtx,
																	 f.neighbors_for_the_internal_node, f.maxelm,
																	 f.diag_coef,
																	 f.alpha, dbeta,
																	 rhie_chow, RCh,
																	 btimedep,
																	 tauparam,
																	 nullptr,
																	 f.mf[iPloc],
																	 nullptr);
													 */

													 // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
													 // см. опыт Гаврилова Анрея sigma cfd.
													 //*
													  my_elmatr_quad_PAm3(iPloc, f.slau,
														  f.slau_bon,
														  f.potent,
														  f.pa,
														  f.prop,
														  f.prop_b,
														  f.nvtx,
														  f.neighbors_for_the_internal_node,
														  f.maxelm,
														  f.alpha,
														  dbeta,
														  rhie_chow, RCh,
														  btimedep,
														  -0.1, // tauparam
														  nullptr,
														  f.mf[iPloc],
														  nullptr,
														  tau,
														  bmyhighorder, bdeltapfinish,
														  bRhieChowiPAM, false, 
														  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
														  sumanb);//*/


												  }
											  }
										  }
									  }
								  }
								  // серийный смыкающий кусок
								  for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
									  integer iPloc = f.ifrontregulationgl[iscan_par];
									  if (iPloc < f.maxelm) {
										  rhie_chow[0][iPloc] = 0.0;
										  rhie_chow[1][iPloc] = 0.0;
										  rhie_chow[2][iPloc] = 0.0;



										  // стационарный солвер.
										  // будем использовать отрицательный шаг по времени,
										  // это будет сигнал использовать среднее арифметическое за опорное время.
										 /*
										  my_elmatr_quad_PAm(iPloc, f.slau,
														 f.slau_bon,
														 f.potent,
														 f.pa, f.prop,
														 f.prop_b, f.nvtx,
														 f.neighbors_for_the_internal_node, f.maxelm,
														 f.diag_coef,
														 f.alpha, dbeta,
														 rhie_chow, RCh,
														 btimedep,
														 tauparam,
														 nullptr,
														 f.mf[iPloc],
														 nullptr);
										 */

										 // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
										 // см. опыт Гаврилова Анрея sigma cfd.
										 //*
										  my_elmatr_quad_PAm3(iPloc, f.slau,
											  f.slau_bon,
											  f.potent,
											  f.pa,
											  f.prop,
											  f.prop_b,
											  f.nvtx,
											  f.neighbors_for_the_internal_node,
											  f.maxelm,
											  f.alpha,
											  dbeta,
											  rhie_chow, RCh,
											  btimedep,
											  -0.1, // tauparam
											  nullptr,
											  f.mf[iPloc],
											  nullptr,
											  tau,
											  bmyhighorder, bdeltapfinish,
											  bRhieChowiPAM, false, 
											  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
											  sumanb);//*/


									  }
								  }


							  }

						  }
						  else {

							  if (nd.b0.active) {

#pragma omp parallel
								  {
#pragma omp sections 
									  {
#pragma omp section
										  {
											  // первый поток
											  for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
												  integer iPloc = f.ifrontregulationgl[iscan_par];
												  if (iPloc < f.maxelm) {
													  rhie_chow[0][iPloc] = 0.0;
													  rhie_chow[1][iPloc] = 0.0;
													  rhie_chow[2][iPloc] = 0.0;


													  // поправка давления
													 /*
													  my_elmatr_quad_PAm(iPloc, f.slau,
																		 f.slau_bon,
																		 f.potent,
																		 f.pa, f.prop,
																		 f.prop_b, f.nvtx,
																		 f.neighbors_for_the_internal_node, f.maxelm,
																		 f.diag_coef,
																		 f.alpha, dbeta,
																		 rhie_chow, RCh,
																		 btimedep,
																		 tauparam,
																		 mfoldtimestep[iPloc],
																		 f.mf[iPloc],
																		 speedoldtimestep);
													 */

													 // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
													 // см. опыт Гаврилова Анрея sigma cfd.
													 ///*
													  my_elmatr_quad_PAm3(iPloc, f.slau,
														  f.slau_bon,
														  f.potent,
														  f.pa,
														  f.prop,
														  f.prop_b,
														  f.nvtx,
														  f.neighbors_for_the_internal_node,
														  f.maxelm,
														  f.alpha, dbeta,
														  rhie_chow, RCh,
														  btimedep,
														  tauparam,
														  mfoldtimestep[iPloc],
														  f.mf[iPloc],
														  speedoldtimestep,
														  tau,
														  bmyhighorder, bdeltapfinish,
														  bRhieChowiPAM, false,
														  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
														  sumanb);
													  //*/

												  }

											  }
										  }
#pragma omp section
										  {
											  // второй поток
											  for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
												  integer iPloc = f.ifrontregulationgl[iscan_par];
												  if (iPloc < f.maxelm) {
													  rhie_chow[0][iPloc] = 0.0;
													  rhie_chow[1][iPloc] = 0.0;
													  rhie_chow[2][iPloc] = 0.0;


													  // поправка давления
													 /*
													  my_elmatr_quad_PAm(iPloc, f.slau,
																		 f.slau_bon,
																		 f.potent,
																		 f.pa, f.prop,
																		 f.prop_b, f.nvtx,
																		 f.neighbors_for_the_internal_node, f.maxelm,
																		 f.diag_coef,
																		 f.alpha, dbeta,
																		 rhie_chow, RCh,
																		 btimedep,
																		 tauparam,
																		 mfoldtimestep[iPloc],
																		 f.mf[iPloc],
																		 speedoldtimestep);
													 */

													 // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
													 // см. опыт Гаврилова Анрея sigma cfd.
													 ///*
													  my_elmatr_quad_PAm3(iPloc, f.slau,
														  f.slau_bon,
														  f.potent,
														  f.pa,
														  f.prop,
														  f.prop_b,
														  f.nvtx,
														  f.neighbors_for_the_internal_node,
														  f.maxelm,
														  f.alpha, dbeta,
														  rhie_chow, RCh,
														  btimedep,
														  tauparam,
														  mfoldtimestep[iPloc],
														  f.mf[iPloc],
														  speedoldtimestep,
														  tau,
														  bmyhighorder, bdeltapfinish,
														  bRhieChowiPAM, false,
														  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
														  sumanb);
													  //*/

												  }
											  }
										  }
									  }
								  }
								  // серийный смыкающий кусок
								  for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
									  integer iPloc = f.ifrontregulationgl[iscan_par];
									  if (iPloc < f.maxelm) {
										  rhie_chow[0][iPloc] = 0.0;
										  rhie_chow[1][iPloc] = 0.0;
										  rhie_chow[2][iPloc] = 0.0;


										  // поправка давления
										 /*
										  my_elmatr_quad_PAm(iPloc, f.slau,
															 f.slau_bon,
															 f.potent,
															 f.pa, f.prop,
															 f.prop_b, f.nvtx,
															 f.neighbors_for_the_internal_node, f.maxelm,
															 f.diag_coef,
															 f.alpha, dbeta,
															 rhie_chow, RCh,
															 btimedep,
															 tauparam,
															 mfoldtimestep[iPloc],
															 f.mf[iPloc],
															 speedoldtimestep);
										 */

										 // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
										 // см. опыт Гаврилова Анрея sigma cfd.
										 ///*
										  my_elmatr_quad_PAm3(iPloc, f.slau,
											  f.slau_bon,
											  f.potent,
											  f.pa,
											  f.prop,
											  f.prop_b,
											  f.nvtx,
											  f.neighbors_for_the_internal_node,
											  f.maxelm,
											  f.alpha, dbeta,
											  rhie_chow, RCh,
											  btimedep,
											  tauparam,
											  mfoldtimestep[iPloc],
											  f.mf[iPloc],
											  speedoldtimestep,
											  tau,
											  bmyhighorder, bdeltapfinish,
											  bRhieChowiPAM, false,
											  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
											  sumanb);
										  //*/

									  }
								  }


							  }

						  }
					  }
				  }
				  else {

#pragma omp parallel for
					  for (integer  i = 0; i<f.maxelm; i++) {
						  rhie_chow[0][i] = 0.0;
						  rhie_chow[1][i] = 0.0;
						  rhie_chow[2][i] = 0.0;

						  if ((mfoldtimestep == nullptr) && (speedoldtimestep == nullptr)) {

							  // стационарный солвер.
							  // будем использовать отрицательный шаг по времени,
							  // это будет сигнал использовать среднее арифметическое за опорное время.
							  /*
							  my_elmatr_quad_PAm(i, f.slau,
							  f.slau_bon,
							  f.potent,
							  f.pa, f.prop,
							  f.prop_b, f.nvtx,
							  f.neighbors_for_the_internal_node, f.maxelm,
							  f.diag_coef,
							  f.alpha, dbeta,
							  rhie_chow, RCh,
							  btimedep,
							  tauparam,
							  nullptr,
							  f.mf[i],
							  nullptr);
							  */

							  // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
							  // см. опыт Гаврилова Анрея sigma cfd.
							  //*
							  my_elmatr_quad_PAm3(i, f.slau,
								  f.slau_bon,
								  f.potent,
								  f.pa,
								  f.prop,
								  f.prop_b,
								  f.nvtx,
								  f.neighbors_for_the_internal_node,
								  f.maxelm,
								  f.alpha,
								  dbeta,
								  rhie_chow, RCh,
								  btimedep,
								  -0.1, // tauparam
								  nullptr,
								  f.mf[i],
								  nullptr,
								  tau,
								  bmyhighorder, bdeltapfinish,
								  bRhieChowiPAM, false,
								  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
								  sumanb);//*/

						  }
						  else
						  {
							  // поправка давления
							  /*
							  my_elmatr_quad_PAm(i, f.slau,
							  f.slau_bon,
							  f.potent,
							  f.pa, f.prop,
							  f.prop_b, f.nvtx,
							  f.neighbors_for_the_internal_node, f.maxelm,
							  f.diag_coef,
							  f.alpha, dbeta,
							  rhie_chow, RCh,
							  btimedep,
							  tauparam,
							  mfoldtimestep[i],
							  f.mf[i],
							  speedoldtimestep);
							  */

							  // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
							  // см. опыт Гаврилова Анрея sigma cfd.
							  ///*
							  if (mfoldtimestep != nullptr) {
								  my_elmatr_quad_PAm3(i, f.slau,
									  f.slau_bon,
									  f.potent,
									  f.pa,
									  f.prop,
									  f.prop_b,
									  f.nvtx,
									  f.neighbors_for_the_internal_node,
									  f.maxelm,
									  f.alpha, dbeta,
									  rhie_chow, RCh,
									  btimedep,
									  tauparam,
									  mfoldtimestep[i],
									  f.mf[i],
									  speedoldtimestep,
									  tau,
									  bmyhighorder, bdeltapfinish,
									  bRhieChowiPAM, false,
									  f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
									  sumanb);
							  }
							  else {
								  printf("Fatal error!!!\n");
								  printf("mfoldtimestep==nullptr error in call  my_elmatr_quad_PAm3.\n");
								  system("PAUSE");
								  exit(1);
							  }
							  //*/
						  }
				   }
				  }

			

#else  

			      for (integer  i=0; i<f.maxelm; i++) {
					   rhie_chow[0][i]=0.0;
                       rhie_chow[1][i]=0.0;
					   rhie_chow[2][i]=0.0;

					   if ((mfoldtimestep==nullptr)&&(speedoldtimestep==nullptr)) {

						   // стационарный солвер.
						   // будем использовать отрицательный шаг по времени,
						   // это будет сигнал использовать среднее арифметическое за опорное время.
						  /*
						   my_elmatr_quad_PAm(i, f.slau, 
						                  f.slau_bon,
										  f.potent, 
										  f.pa, f.prop,
										  f.prop_b, f.nvtx,
										  f.neighbors_for_the_internal_node, f.maxelm,
										  f.diag_coef,
										  f.alpha, dbeta,
										  rhie_chow, RCh,
										  btimedep,
										  tauparam,
										  nullptr,
										  f.mf[i],
										  nullptr);
						  */

						   // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
						   // см. опыт Гаврилова Анрея sigma cfd.
						   //*
                           my_elmatr_quad_PAm3(i, f.slau,
							                   f.slau_bon,  
						                       f.potent,
											   f.pa,
											   f.prop,
											   f.prop_b,
					                           f.nvtx, 
											   f.neighbors_for_the_internal_node,
											   f.maxelm, 
						                       f.alpha,
											   dbeta,
											   rhie_chow, RCh,
						                       btimedep, 
											   -0.1, // tauparam
											   nullptr,
						                       f.mf[i],
											   nullptr,
											   tau,
											   bmyhighorder, bdeltapfinish,
											   bRhieChowiPAM,false,
							                   f.border_neighbor, t.ilevel_alice, f.ptr, f.maxbound,
							                   sumanb);//*/

					   }
					   else
					   {
					   // поправка давления
					  /*
					   my_elmatr_quad_PAm(i, f.slau, 
						                  f.slau_bon,
										  f.potent, 
										  f.pa, f.prop,
										  f.prop_b, f.nvtx,
										  f.neighbors_for_the_internal_node, f.maxelm,
										  f.diag_coef,
										  f.alpha, dbeta,
										  rhie_chow, RCh,
										  btimedep,
										  tauparam,
										  mfoldtimestep[i],
										  f.mf[i],
										  speedoldtimestep);
					  */

					   // Рекомендуется собирать матрицу на основе сглаженного псевдовремени.
					   // см. опыт Гаврилова Анрея sigma cfd.
					   ///*
						   if (mfoldtimestep != nullptr) {
							   my_elmatr_quad_PAm3(i, f.slau,
								   f.slau_bon,
								   f.potent,
								   f.pa,
								   f.prop,
								   f.prop_b,
								   f.nvtx,
								   f.neighbors_for_the_internal_node,
								   f.maxelm,
								   f.alpha, dbeta,
								   rhie_chow, RCh,
								   btimedep,
								   tauparam,
								   mfoldtimestep[i],
								   f.mf[i],
								   speedoldtimestep,
								   tau,
								   bmyhighorder, bdeltapfinish,
								   bRhieChowiPAM, false, 
								   f.border_neighbor, t.ilevel_alice, f.ptr,f.maxbound,
								   sumanb);
						   }
						   else {
							   printf("Fatal error!!!\n");
							   printf("mfoldtimestep==nullptr error in call  my_elmatr_quad_PAm3.\n");
							   system("PAUSE");
							   exit(1);
						   }
										   //*/
					   }
				   }

#endif

				  // printf("step internal...\n");
				  // getchar();
				   
				   #ifdef _OPENMP
				   if (bparallelizm_old) {
					   printf("error bparallelizm_old\n");
					   system("pause");

				  if (inumcore == 1) {
					  // serial
					  // Нормировка:
					  for (integer  i = 0; i<f.maxelm; i++) {
						  f.slau[PAM][i].ae /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].aw /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].an /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].as /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].at /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].ab /= f.slau[PAM][i].ap;
						  // АЛИС сетка.
						  f.slau[PAM][i].ae2 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].aw2 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].an2 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].as2 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].at2 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].ab2 /= f.slau[PAM][i].ap;

						  f.slau[PAM][i].ae3 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].aw3 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].an3 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].as3 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].at3 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].ab3 /= f.slau[PAM][i].ap;

						  f.slau[PAM][i].ae4 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].aw4 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].an4 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].as4 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].at4 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].ab4 /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].b /= f.slau[PAM][i].ap;
						  f.slau[PAM][i].ap = 1.0;
					  }
					  for (integer i = 0; i<f.maxbound; i++) {
						  f.slau_bon[PAM][i].ai /= f.slau_bon[PAM][i].aw;
						  f.slau_bon[PAM][i].b /= f.slau_bon[PAM][i].aw;
						  f.slau_bon[PAM][i].aw = 1.0;
					  }
				  }

			if (inumcore==2) {
				if (nd.b0.active) {

					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=f.ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
                               f.slau[PAM][iPloc].ae/=f.slau[PAM][iPloc].ap;
						       f.slau[PAM][iPloc].aw/=f.slau[PAM][iPloc].ap;
						       f.slau[PAM][iPloc].an/=f.slau[PAM][iPloc].ap;
						       f.slau[PAM][iPloc].as/=f.slau[PAM][iPloc].ap;
				        	   f.slau[PAM][iPloc].at/=f.slau[PAM][iPloc].ap;
						       f.slau[PAM][iPloc].ab/=f.slau[PAM][iPloc].ap;

							   // АЛИС сетка.
							   f.slau[PAM][iPloc].ae2 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].aw2 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].an2 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].as2 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].at2 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].ab2 /= f.slau[PAM][iPloc].ap;

							   f.slau[PAM][iPloc].ae3 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].aw3 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].an3 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].as3 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].at3 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].ab3 /= f.slau[PAM][iPloc].ap;

							   f.slau[PAM][iPloc].ae4 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].aw4 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].an4 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].as4 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].at4 /= f.slau[PAM][iPloc].ap;
							   f.slau[PAM][iPloc].ab4 /= f.slau[PAM][iPloc].ap;
						       f.slau[PAM][iPloc].b/=f.slau[PAM][iPloc].ap;
						       f.slau[PAM][iPloc].ap=1.0;
						}

					}
					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=f.ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							 f.slau[PAM][iPloc].ae/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].aw/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].an/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].as/=f.slau[PAM][iPloc].ap;
				        	 f.slau[PAM][iPloc].at/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].ab/=f.slau[PAM][iPloc].ap;
							 // АЛИС сетка.
							 f.slau[PAM][iPloc].ae2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].aw2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].an2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].as2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].at2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].ab2 /= f.slau[PAM][iPloc].ap;

							 f.slau[PAM][iPloc].ae3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].aw3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].an3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].as3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].at3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].ab3 /= f.slau[PAM][iPloc].ap;

							 f.slau[PAM][iPloc].ae4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].aw4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].an4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].as4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].at4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].ab4 /= f.slau[PAM][iPloc].ap;

						     f.slau[PAM][iPloc].b/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].ap=1.0;
						}

					}
					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=f.ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							 f.slau[PAM][iPloc].ae/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].aw/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].an/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].as/=f.slau[PAM][iPloc].ap;
				        	 f.slau[PAM][iPloc].at/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].ab/=f.slau[PAM][iPloc].ap;
							 // АЛИС сетка.
							 f.slau[PAM][iPloc].ae2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].aw2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].an2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].as2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].at2 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].ab2 /= f.slau[PAM][iPloc].ap;

							 f.slau[PAM][iPloc].ae3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].aw3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].an3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].as3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].at3 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].ab3 /= f.slau[PAM][iPloc].ap;

							 f.slau[PAM][iPloc].ae4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].aw4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].an4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].as4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].at4 /= f.slau[PAM][iPloc].ap;
							 f.slau[PAM][iPloc].ab4 /= f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].b/=f.slau[PAM][iPloc].ap;
						     f.slau[PAM][iPloc].ap=1.0;
						}

					}


				}
			}

			if (inumcore==2) {
				if (nd.b0.active) {

					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=f.ifrontregulationgl[iscan_par];
						if (iPloc>=f.maxelm) {
							iPloc-=f.maxelm;
							f.slau_bon[PAM][iPloc].ai/=f.slau_bon[PAM][iPloc].aw;
					    	f.slau_bon[PAM][iPloc].b/=f.slau_bon[PAM][iPloc].aw;
					    	f.slau_bon[PAM][iPloc].aw=1.0;
						}

					}
					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=f.ifrontregulationgl[iscan_par];
						if (iPloc>=f.maxelm) {
                            iPloc-=f.maxelm;
							f.slau_bon[PAM][iPloc].ai/=f.slau_bon[PAM][iPloc].aw;
					    	f.slau_bon[PAM][iPloc].b/=f.slau_bon[PAM][iPloc].aw;
					    	f.slau_bon[PAM][iPloc].aw=1.0;
						}

					}
					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=f.ifrontregulationgl[iscan_par];
						if (iPloc>=f.maxelm) {
							iPloc-=f.maxelm;
							f.slau_bon[PAM][iPloc].ai/=f.slau_bon[PAM][iPloc].aw;
					    	f.slau_bon[PAM][iPloc].b/=f.slau_bon[PAM][iPloc].aw;
					    	f.slau_bon[PAM][iPloc].aw=1.0;
						}

					}


				}
			}
			}
				   else {
					   // Нормировка:
					   // Дополнено 19.03.2019
#pragma omp parallel for
					   for (integer  i = 0; i<f.maxelm; i++) {
						   f.slau[PAM][i].ae /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].aw /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].an /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].as /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].at /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].ab /= f.slau[PAM][i].ap;
						   // АЛИС сетка.
						   f.slau[PAM][i].ae2 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].aw2 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].an2 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].as2 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].at2 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].ab2 /= f.slau[PAM][i].ap;

						   f.slau[PAM][i].ae3 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].aw3 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].an3 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].as3 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].at3 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].ab3 /= f.slau[PAM][i].ap;

						   f.slau[PAM][i].ae4 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].aw4 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].an4 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].as4 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].at4 /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].ab4 /= f.slau[PAM][i].ap;

						   f.slau[PAM][i].b /= f.slau[PAM][i].ap;
						   f.slau[PAM][i].ap = 1.0;
					   }
#pragma omp parallel for
					   for (integer  i = 0; i<f.maxbound; i++) {
						   f.slau_bon[PAM][i].ai /= f.slau_bon[PAM][i].aw;
						   f.slau_bon[PAM][i].b /= f.slau_bon[PAM][i].aw;
						   f.slau_bon[PAM][i].aw = 1.0;
					   }
				   }

#else
				   
				   // Нормировка:
				   // Дополнено 19.03.2019
				    for (integer  i=0; i<f.maxelm; i++) {
						f.slau[PAM][i].ae/=f.slau[PAM][i].ap;
						f.slau[PAM][i].aw/=f.slau[PAM][i].ap;
						f.slau[PAM][i].an/=f.slau[PAM][i].ap;
						f.slau[PAM][i].as/=f.slau[PAM][i].ap;
						f.slau[PAM][i].at/=f.slau[PAM][i].ap;
						f.slau[PAM][i].ab/=f.slau[PAM][i].ap;
						// АЛИС сетка.
						f.slau[PAM][i].ae2 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].aw2 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].an2 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].as2 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].at2 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].ab2 /= f.slau[PAM][i].ap;

						f.slau[PAM][i].ae3 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].aw3 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].an3 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].as3 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].at3 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].ab3 /= f.slau[PAM][i].ap;

						f.slau[PAM][i].ae4 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].aw4 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].an4 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].as4 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].at4 /= f.slau[PAM][i].ap;
						f.slau[PAM][i].ab4 /= f.slau[PAM][i].ap;

						f.slau[PAM][i].b/=f.slau[PAM][i].ap;
						f.slau[PAM][i].ap=1.0;
					}
					for (integer  i=0; i<f.maxbound; i++) {
						f.slau_bon[PAM][i].ai/=f.slau_bon[PAM][i].aw;
						f.slau_bon[PAM][i].b/=f.slau_bon[PAM][i].aw;
						f.slau_bon[PAM][i].aw=1.0;
					}
					
					
#endif

#if doubleintprecision == 1
					// debug
					//if (inumiter==2) {
					/*
						for (integer  i=0; i<f.maxelm; i++) {
							printf("id=%lld ap=%e ae=%e aw=%e an=%e as=%e at=%e ab=%e b=%e\n",i, f.slau[PAM][i].ap, f.slau[PAM][i].ae, f.slau[PAM][i].aw, f.slau[PAM][i].an, f.slau[PAM][i].as, f.slau[PAM][i].at, f.slau[PAM][i].ab, f.slau[PAM][i].b);
							getchar();
						}*/
					//}

					/*
						for (integer  i=0; i<f.maxbound; i++) {
							if (fabs(f.slau_bon[PAM][i].ai)<1.0e-20) {
								printf("id=%lld aw=%e ai=%e b=%e %lld %lld\n",i,f.slau_bon[PAM][i].aw,f.slau_bon[PAM][i].ai,f.slau_bon[PAM][i].b, f.slau_bon[PAM][i].iW, f.slau_bon[PAM][i].iI);
								getchar();
							}
						}
					*/
#else
					// debug
					//if (inumiter==2) {
						/*
						for (integer  i=0; i<f.maxelm; i++) {
							printf("id=%d ap=%e ae=%e aw=%e an=%e as=%e at=%e ab=%e b=%e\n",i, f.slau[PAM][i].ap, f.slau[PAM][i].ae, f.slau[PAM][i].aw, f.slau[PAM][i].an, f.slau[PAM][i].as, f.slau[PAM][i].at, f.slau[PAM][i].ab, f.slau[PAM][i].b);
							getchar();
						}*/
					//}

					/*
					for (integer i=0; i<f.maxbound; i++) {
						if (fabs(f.slau_bon[PAM][i].ai)<1.0e-20) {
							printf("id=%d aw=%e ai=%e b=%e %d %d\n",i,f.slau_bon[PAM][i].aw,f.slau_bon[PAM][i].ai,f.slau_bon[PAM][i].b, f.slau_bon[PAM][i].iW, f.slau_bon[PAM][i].iI);
							getchar();
						}
					}
					*/
#endif
#pragma omp parallel for
					for (integer  i = 0; i < f.maxelm; i++) {
						if (f.slau[PAM][i].ap < 0.0) {
							printf("maxelm=%lld i=%lld problem PAM\n",f.maxelm, i);
							system("pause");
							//f.slau_bon[iVar][i].aw
						}
					}
				  

			       //getchar(); // debug
				   // В res записывается норма нескомпенсированных источников массы:
				   // только для внутренних КО.
				   // По этому показателю можно судить о 
				  // сходимости всей системы уравнений Навье-Стокса.
				  if (0) {
				  for (integer  i=0; i<f.maxelm+f.maxbound; i++) {
					  if (i<f.maxelm) {
						  f.potent[FBUF][i]=f.slau[PAM][i].b;
					  }
					  else {
						  f.potent[FBUF][i]=f.slau_bon[PAM][i-f.maxelm].b;
					  }
				  }
				  } // debug
				  if (0) {
#if doubleintprecision == 1
					  for (integer  i = 0; i<f.maxelm; i++) {
						  if (f.slau[PAM][i].b >= 1e-20) {
							  printf("zero internal elem %lld", i);
							  //getchar();
							  system("pause");
					  }

				  }
					  for (integer  i = 0; i<f.maxbound; i++) {
						  if (f.slau_bon[PAM][i].b >= 1e-20) {
							  printf("zero boundary elem %lld", i);
							  //getchar();
							  system("pause");
						  }
					  }
				  }
#else
					  for (integer  i = 0; i<f.maxelm; i++) {
						  if (f.slau[PAM][i].b >= 1e-20) {
							  printf("zero internal elem %d", i);
							  //getchar();
							  system("pause");
						  }

					  }
					  for (integer  i = 0; i<f.maxbound; i++) {
						  if (f.slau_bon[PAM][i].b >= 1e-20) {
							  printf("zero boundary elem %d", i);
							  //getchar();
							  system("pause");
						  }
					  }
				  }
#endif
				  

				  if (0) {
		              xyplot( fglobal, flow_interior, t);
		              printf("xy plot continity calc. OK.\n");
	                 // getchar(); // debug avtosave
					  system("pause");
	              }
				  // нельзя этого делать.
				  /*
				  if (bHORF) {
					  doublereal fHORF = 0.25;
					  if (btimedep) {
						  fHORF = 0.75;
					  }
					  for (i = 0; i < f.maxelm; i++) {
						  f.slau[PAM][i].b = bPamendment_source_old[i] + fHORF*(f.slau[PAM][i].b - bPamendment_source_old[i]);
						  bPamendment_source_old[i] = f.slau[PAM][i].b;
					  }
					  for (i = 0; i < f.maxbound; i++) {
						  f.slau_bon[PAM][i].b = bPamendment_source_old[i + f.maxelm] + fHORF*(f.slau_bon[PAM][i].b - bPamendment_source_old[i + f.maxelm]);
						  bPamendment_source_old[i + f.maxelm] = f.slau_bon[PAM][i].b;
					  }
				  }
				  */
			       res=0.0;
			       for (integer  i=0; i<f.maxelm; i++) {
					   // Здесь выбрана норма Чебышева т.к.
					   // она наибольшая из всех норм. (худший случай).
					   //if (inumiter==2) printf("b=%e, res=%e\n",f.slau[PAM][i].b,res);
					   //if (i%10==0) getchar(); // debug
				       res=fmax(res,fabs(f.slau[PAM][i].b));
				       //res+=fabs(f.slau[PAM][i].b);
				       //res+=f.slau[PAM][i].b*f.slau[PAM][i].b;
			       }
			       //res=sqrt((doublereal)(res/f.maxelm));

			       break;

		case TEMP:   
			        
			      

			      // 3 мая 2016. ошибка с 2D источником тепла.
			      t.alpha = 1.0; // Это очень важно иначе одна сплошная недоэтерированность.
				  if ((fabs(dgx) > 1.0e-20) || (fabs(dgy) > 1.0e-20) || (fabs(dgz) > 1.0e-20)) {
					  if (bSIMPLErun_now_for_temperature) {
						  // Естественная конвекция.
						  // Для задач с естественной конвекцией (CFD + Temp)
						  // Обязательно нужно вводить нижнюю релаксацию на температуру.
						  // 06.08.2019.
						  // При t.alpha =1.0 сходимость была очень плохой. Расходимость, переполнения.
						  // t.alpha = 0.9; при таком alpha температуры очень сильно занижены.
						  // Управляемый коэффициент нижней релаксации.
						  //t.alpha = my_amg_manager.F_to_F_Stress;
						  // Компромисс температуры немного занижены но сходимость нормальная,
						  // переполнений нет.
						  t.alpha = 0.99999;// именно пять девяток. Переполнений нет. Сходимость.
						  // Точность значения температуры очень сильно зависит от качества расчетной сетки.
					  }
				  }

				  if (sourse2Dproblem != nullptr) {
					  delete[] sourse2Dproblem;
					  sourse2Dproblem = nullptr;
					  sourse2Dproblem = new bool[t.maxbound];
				  }

				  if (conductivity2Dinsource != nullptr) {
					  delete[] conductivity2Dinsource;
					  conductivity2Dinsource = nullptr;
					  conductivity2Dinsource = new doublereal[t.maxbound];
				  }

#pragma omp parallel for
			      for (integer i74 = 0; i74 < t.maxbound; i74++) {
				        sourse2Dproblem[i74] = true;
						conductivity2Dinsource[i74] = -2.0; // инициализация любым отрицательным значением.
			      }
				  // Вычисляет теплопроводность на грани внутреннего источника тепла
				  // используется в граничном условии на источнике тепла.
				  for (integer  i = 0; i < t.maxelm; i++) {
					  conduct2Dsourceconstruct(i, t.slau, t.neighbors_for_the_internal_node, t.maxelm,
						  t.border_neighbor, ls, t.prop);
				  }

			       // теперь внутрь функции передаётся параметр который показывает
			       // стационарную или нестационарную задачу мы решаем
			       // btimedep=false; // стационарный
			      
				 
                   // симметризация СЛАУ
			       // Теперь симметризация матрицы производится внутри
				   // сборки для каждого контрольного объёма отдельно.
				   // Это позволило значительно ускорить вычисления, 
				   // путём избавления от двойного цикла 
				   // который был реализован раньше и сильно замедлял вычисления.

				   // учёт граничных условий для 
                   // уравнения теплопроводности
#pragma omp parallel for
				   for (integer  i=0; i<t.maxbound; i++) {
					   // Условия Дирихле:
					   my_elmatr_quad_T3D_bound(i, t.maxbound,
						                        t.maxelm,
												t.border_neighbor,
												ls, lw, w, s, 
												t.slau_bon, 
												true, dbeta,
												t.nvtx, t.pa, 
												t.prop,
												t.prop_b, 
												t.potent,
												told_iter,
												t.ptr,
												fglobal,
												poweron_multiplier_sequence);
				   }

				 
#pragma omp parallel for
                   for (integer  i=0; i<t.maxbound; i++) {
					   // Условия Неймана однородные и неоднородные:
					   my_elmatr_quad_T3D_bound(i, t.maxbound,
						                        t.maxelm,
												t.border_neighbor,
												ls, lw, w, s,
												t.slau_bon,
												false, dbeta,
												t.nvtx, t.pa,
												t.prop, 
												t.prop_b,
												t.potent,
												told_iter,
												t.ptr,
												fglobal,
												poweron_multiplier_sequence);
				   }
				 

				   imyscheme = UDS;
				   if (bVERYStable) {
					   // стабильная противопоточная схема первого порядка
					   imyscheme = UDS; // UDS BULG
				   }
				   else {
					   imyscheme = UNEVEN_WACEB;
					   //imyscheme=UDS; 
				   }
				   imyscheme = UDS; //UNEVEN_WACEB;
				   imyscheme = UNEVEN_WACEB; // 30 07 2015.
											 //imyscheme=UNEVENQUICK;// не работает на 30 07 2015.
											 // Выбор делается из графического интерфейса пользователя !
				   imyscheme = iTEMPScheme; //  5 02 2017 // 31 07 2015


				   // Сборка строк матрицы для внутренних КО.

#pragma omp parallel for
			       for (integer  i=0; i<t.maxelm; i++) {
					   if (toldtimestep==nullptr) {
						   my_elmatr_quad_T3D(i, t.slau,
						                  t.slau_bon,
										  btimedep, // стационарный или нестационарный
										  tauparam, // размер шага по времени.
							              imyscheme, // UDS BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
										  t.nvtx,
										  t.potent, 
										  t.pa, t.prop,
										  t.prop_b,
										  t.neighbors_for_the_internal_node, 
										  t.alpha, 
										  dbeta, 
										  bconvective,
										  t.Sc[i], 
                                          t.ipower_time_depend[i],
										  t.maxelm,
										  flow_interior,
										  t.ptr,
										  t.border_neighbor,
										  ls,fglobal,
										  t.binternalsource,
										  nullptr, // поле температур с предыдущего временного слоя
										  b,lb,
										  matlist,
										  inumiter, s, 
                                          poweron_multiplier_sequence,
										  t.ilevel_alice,
							              t.whot_is_block);
				
					   }
					   else
					   {
						   //if (b[t.whot_is_block[i]].ipower_time_depend > 0) {
							   //printf("%d %d\n", t.ipower_time_depend[i], b[t.whot_is_block[i]].ipower_time_depend);
							   //getchar();
						   //}
					   my_elmatr_quad_T3D(i, t.slau,
						                  t.slau_bon,
										  btimedep, // стационарный или нестационарный
										  tauparam, // размер шага по времени.
						                  imyscheme, // UDS BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
										  t.nvtx,
										  t.potent, 
										  t.pa, t.prop,
										  t.prop_b,
										  t.neighbors_for_the_internal_node, 
										  t.alpha, 
										  dbeta, 
										  bconvective,
										  t.Sc[i], 
                                          t.ipower_time_depend[i],
										  t.maxelm,
										  flow_interior,
										  t.ptr,
										  t.border_neighbor,
										  ls,fglobal,
										  t.binternalsource,
										  toldtimestep, // поле температур с предыдущего временного слоя
										  b,lb,
										  matlist,
										  inumiter, s, 
										  poweron_multiplier_sequence,
										  t.ilevel_alice,
						                  t.whot_is_block);
					   }
				   }
                   

				   
				  


				   // debug
				   if ((0)&&(inumiter>=66)) {
					   // печать матрицы в файл.
					   print_temp_slau(t);
					   printf("temperature matrix print...\n");
					   //getchar();
					   system("pause");
				   }

                   /*
                   if (t.free_temper_level1) {
					  free_level1_temp(t);
                      printf("memory: 43N+15N-42N=16N\n");
				      getchar();
				   }
				   */
                   res=0.0;
				   
				   for (integer i74 = 0; i74 < t.maxbound; i74++) {
					   if (sourse2Dproblem[i74] == false) {
						  // printf("found non internal source.\n");
						   i75++;
					   }
				   }
#if doubleintprecision == 1
				   // printf("found non internal source %lld \n ",i75);
#else
				   // printf("found non internal source %d \n ",i75);
#endif
				  

				   // для SIMPLE алгоритма.
				   rfluentresval = fluent_residual_for_x(t.slau, t.slau_bon, t.potent, t.maxelm, t.maxbound, TEMP); // невязка по формуле fluent.

				   break;

				   case NUSHA:
					   // Модифицированная кинематическая турбулентная вязкость.

					   // Граничные условия Дирихле обязательно 
				       // должны собираться в первую очередь
					   for (integer  i = 0; i < f.maxbound; i++) {
						   // условия Дирихле третий параметр равен true
						   my_elmatr_quad_SpallartAllmares3D_bound(i, f.maxelm,
							   true, f.border_neighbor, ls, lw, w,
							   f.slau_bon[NUSHA_SL],
							   f.pa, f.nvtx, f.prop_b, f.prop,
							   f.potent);
					   }

					   // Во вторую очередь собираются однородные условия Неймана.
					   for (integer i = 0; i < f.maxbound; i++) {
						   // однородные условия Неймана третий параметр равен false
						   my_elmatr_quad_SpallartAllmares3D_bound(i, f.maxelm,
							   false, f.border_neighbor, ls, lw, w,
							   f.slau_bon[NUSHA_SL],
							   f.pa, f.nvtx, f.prop_b, f.prop,
							   f.potent);
					   }

					   // Схема для аппроксимации конвективного члена в модели 
					   // Спаларта Аллмареса такая же как и для аппроксимации компонент 
					   // вектора скорости.

					   // Выбор делается из графического интерфейса пользователя !
					   imyscheme = iFLOWScheme; // 31 07 2015

					   // Сборка строк матрицы для внутренних КО.
					   for (integer  i = 0; i < f.maxelm; i++) {
						   my_elmatr_quad_SpallartAllmares3D(
							   i,
							   f.border_neighbor,
							   lw, ls,
							   f.slau,
							   f.slau_bon,
							   //doublereal** diag_coef,
							   //integer iVar,
							   //bool btimedep,
							   //doublereal tauparam,
							   f.ptr,
							   f.nvtx,
							   f.potent,
							   //doublereal* potent_temper,
							   f.pa,
							   f.prop,
							   f.prop_b,
							   f.maxelm,
							   f.neighbors_for_the_internal_node,
							   //doublereal* alpha,
							   //doublereal dgx,
							   //doublereal dgy,
							   //doublereal dgz,
							   //doublereal dbeta, 
							   imyscheme,
							   //bool bBussineskApproach,
							   //doublereal temp_ref,
							   //bool bfirst_start,
							   //doublereal RCh,
							   //integer iflowregime,
							   //doublereal* speedoldtimestep,
							   //doublereal* toldtimestep,
							   b, lb, matlist,
							   f.mf,
							   //bool bVERYStable,
							   //doublereal &sumanb,
							   t.ilevel_alice,
							   f.rdistWall,
							   f.SInvariantStrainRateTensor
						   );
					   }


					break;

				   case TURBULENT_KINETIK_ENERGY:
					   // Кинетическая энергия турбулентных пульсаций в ммодели SST Ментера.

					   // init для проверки
#pragma omp parallel for
					   for (integer i = 0; i < f.maxbound; i++) {
						   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL][i].ai = 0.0;
						   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL][i].aw = -100.0;
						   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL][i].b = 0.0;
					   }

					   // Граничные условия Дирихле обязательно 
					   // должны собираться в первую очередь
#pragma omp parallel for
					   for (integer  i = 0; i < f.maxbound; i++) {
						   // условия Дирихле третий параметр равен true
						   my_elmatr_quad_kinetik_turbulence_energy_3D_bound(i, f.maxelm,
							   true, f.border_neighbor, ls, lw, w,
							   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL],
							   f.pa, f.nvtx, f.prop_b, f.prop,
							   f.potent);
					   }

					   // Во вторую очередь собираются однородные условия Неймана.
#pragma omp parallel for
					   for (integer  i = 0; i < f.maxbound; i++) {
						   // однородные условия Неймана третий параметр равен false
						   my_elmatr_quad_kinetik_turbulence_energy_3D_bound(i, f.maxelm,
							   false, f.border_neighbor, ls, lw, w,
							   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL],
							   f.pa, f.nvtx, f.prop_b, f.prop,
							   f.potent);
					   }

					   // Проверка.
#pragma omp parallel for
					   for (integer  i = 0; i < f.maxbound; i++) {
						   if (f.slau_bon[TURBULENT_KINETIK_ENERGY_SL][i].aw <0.0) {
							   printf("maxbound=%lld i=%lld problem KE\n",f.maxbound, i);
							   system("pause");
							   //f.slau_bon[TURBULENT_KINETIK_ENERGY_SL][i].aw
						   }
					   }

					   // Схема для аппроксимации конвективного члена в модели 
					   // SST Ментера такая же как и для аппроксимации компонент 
					   // вектора скорости.

					   // Выбор делается из графического интерфейса пользователя !
					   imyscheme = iFLOWScheme; // 31 07 2015


					   

					   // Сборка строк матрицы для внутренних КО.
#pragma omp parallel for
					   for (integer  i = 0; i < f.maxelm; i++) {

						   

						   my_elmatr_quad_turbulent_kinetik_energy_MenterSST_3D(
							   i,
							   f.border_neighbor,
							   lw, ls,
							   f.slau,
							   f.slau_bon,
							   //doublereal** diag_coef,
							   //integer iVar,
							   //bool btimedep,
							   //doublereal tauparam,
							   f.ptr,
							   f.nvtx,
							   f.potent,
							   //doublereal* potent_temper,
							   f.pa,
							   f.prop,
							   f.prop_b,
							   f.maxelm,
							   f.neighbors_for_the_internal_node,
							   //doublereal* alpha,
							   //doublereal dgx,
							   //doublereal dgy,
							   //doublereal dgz,
							   //doublereal dbeta, 
							   imyscheme,
							   //bool bBussineskApproach,
							   //doublereal temp_ref,
							   //bool bfirst_start,
							   //doublereal RCh,
							   //integer iflowregime,
							   //doublereal* speedoldtimestep,
							   //doublereal* toldtimestep,
							   b, lb, matlist,
							   f.mf,
							   //bool bVERYStable,
							   //doublereal &sumanb,
							   t.ilevel_alice,
							   f.rdistWall,
							   f.SInvariantStrainRateTensor,
							   brthdsd_ON
						   );
					   }


#pragma omp parallel for
					   for (integer  i = 0; i < f.maxelm; i++) {
						   if (f.slau[TURBULENT_KINETIK_ENERGY_SL][i].ap < 0.0) {
							   printf("maxelm=%lld i=%lld problem TURBULENT_KINETIK_ENERGY_SL ap=%e\n", f.maxelm, i, f.slau[TURBULENT_KINETIK_ENERGY_SL][i].ap);
							   system("pause");
							   //f.slau_bon[iVar][i].aw
						   }
					   }

					   break;

				   case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
					   // Удельная скорость диссипации кинетической энергии
					   // турбулентных пульсаций.

						// init для проверки
#pragma omp parallel for
					   for (integer  i = 0; i < f.maxbound; i++) {
						   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].ai = 0.0;
						   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].aw = -100.0;
						   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].b = 0.0;
					   }

					   // Граничные условия Дирихле обязательно 
					   // должны собираться в первую очередь
#pragma omp parallel for
					   for (integer  i = 0; i < f.maxbound; i++) {
						   // условия Дирихле третий параметр равен true
						   my_elmatr_quad_OmegaSSTMenter3D_bound(i, f.maxelm,
							   true, f.border_neighbor, ls, lw, w,
							   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
							   f.pa, f.nvtx, f.prop_b, f.prop,
							   f.potent);
					   }

					   // Во вторую очередь собираются однородные условия Неймана.
#pragma omp parallel for
					   for (integer  i = 0; i < f.maxbound; i++) {
						   // однородные условия Неймана третий параметр равен false
						   my_elmatr_quad_OmegaSSTMenter3D_bound(i, f.maxelm,
							   false, f.border_neighbor, ls, lw, w,
							   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
							   f.pa, f.nvtx, f.prop_b, f.prop,
							   f.potent);
					   }

					   // Проверка.
#pragma omp parallel for
					   for (integer  i = 0; i < f.maxbound; i++) {
						   if (f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].aw < 0.0) {
							   printf("maxbound=%lld i=%lld problem OMEGA\n", f.maxbound, i);
							   system("pause");
							   //f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].aw
						   }
					   }

					   // Схема для аппроксимации конвективного члена в модели 
					   // SST Ментера такая же как и для аппроксимации компонент 
					   // вектора скорости.

					   // Выбор делается из графического интерфейса пользователя !
					   imyscheme = iFLOWScheme; // 31 07 2015

					   

					   // Сборка строк матрицы для внутренних КО.
#pragma omp parallel for
					   for (integer  i = 0; i < f.maxelm; i++) {					   

						   my_elmatr_quad_specific_dissipation_rate_omega_MenterSST3D(
							   i,
							   f.border_neighbor,
							   lw, ls,
							   f.slau,
							   f.slau_bon,
							   //doublereal** diag_coef,
							   //integer iVar,
							   //bool btimedep,
							   //doublereal tauparam,
							   f.ptr,
							   f.nvtx,
							   f.potent,
							   //doublereal* potent_temper,
							   f.pa,
							   f.prop,
							   f.prop_b,
							   f.maxelm,
							   f.neighbors_for_the_internal_node,
							   //doublereal* alpha,
							   //doublereal dgx,
							   //doublereal dgy,
							   //doublereal dgz,
							   //doublereal dbeta, 
							   imyscheme,
							   //bool bBussineskApproach,
							   //doublereal temp_ref,
							   //bool bfirst_start,
							   //doublereal RCh,
							   //integer iflowregime,
							   //doublereal* speedoldtimestep,
							   //doublereal* toldtimestep,
							   b, lb, matlist,
							   f.mf,
							   //bool bVERYStable,
							   //doublereal &sumanb,
							   t.ilevel_alice,
							   f.rdistWall,
							   f.SInvariantStrainRateTensor,
							   brthdsd_ON1
						   );
					   }

#pragma omp parallel for
					   for (integer  i = 0; i < f.maxelm; i++) {
						   if (f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].ap < 0.0) {
							   printf("maxelm=%lld i=%lld problem TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL ap=%e\n", f.maxelm, i, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].ap);
							   system("pause");
							   //f.slau_bon[iVar][i].aw
						   }
					   }

					   break;

					case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
							// Модифицированная кинематическая турбулентная вязкость.

							// Граничные условия Дирихле обязательно 
							// должны собираться в первую очередь
							for (integer  i = 0; i < f.maxbound; i++) {
								// условия Дирихле третий параметр равен true
								my_elmatr_quad_kinetik_turbulence_energy_3D_bound_standart_k_epsilon(i, f.maxelm,
									true, f.border_neighbor, ls, lw, w,
									f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
									f.pa, f.nvtx, f.prop_b, f.prop,
									f.potent);
							}

							// Во вторую очередь собираются однородные условия Неймана.
							for (integer  i = 0; i < f.maxbound; i++) {
								// однородные условия Неймана третий параметр равен false
								my_elmatr_quad_kinetik_turbulence_energy_3D_bound_standart_k_epsilon(i, f.maxelm,
									false, f.border_neighbor, ls, lw, w,
									f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
									f.pa, f.nvtx, f.prop_b, f.prop,
									f.potent);
							}

							// Схема для аппроксимации конвективного члена в модели 
							// Спаларта Аллмареса такая же как и для аппроксимации компонент 
							// вектора скорости.

							// Выбор делается из графического интерфейса пользователя !
							imyscheme = iFLOWScheme; // 31 07 2015

													 // Сборка строк матрицы для внутренних КО.
							for (integer  i = 0; i < f.maxelm; i++) {
								my_elmatr_quad_turbulent_kinetik_energy_Standart_KE_3D(
									i,
									f.border_neighbor,
									lw, ls,
									f.slau,
									f.slau_bon,
									//doublereal** diag_coef,
									//integer iVar,
									//bool btimedep,
									//doublereal tauparam,
									f.ptr,
									f.nvtx,
									f.potent,
									//doublereal* potent_temper,
									f.pa,
									f.prop,
									f.prop_b,
									f.maxelm,
									f.neighbors_for_the_internal_node,
									//doublereal* alpha,
									//doublereal dgx,
									//doublereal dgy,
									//doublereal dgz,
									//doublereal dbeta, 
									imyscheme,
									//bool bBussineskApproach,
									//doublereal temp_ref,
									//bool bfirst_start,
									//doublereal RCh,
									//integer iflowregime,
									//doublereal* speedoldtimestep,
									//doublereal* toldtimestep,
									b, lb, matlist,
									f.mf,
									//bool bVERYStable,
									//doublereal &sumanb,
									t.ilevel_alice,
									f.rdistWall,
									f.SInvariantStrainRateTensor
								);
							}


							break;

					case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
						// Модифицированная кинематическая турбулентная вязкость.

						// Граничные условия Дирихле обязательно 
						// должны собираться в первую очередь
						for (integer  i = 0; i < f.maxbound; i++) {
							// условия Дирихле третий параметр равен true
							my_elmatr_quad_dissipation_rate_epsilon_3D_bound_standart_k_epsilon(i, f.maxelm,
								true, f.border_neighbor, ls, lw, w,
								f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
								f.pa, f.nvtx, f.prop_b, f.prop,
								f.potent, f.rdistWall, f.rdistWallmax);
						}

						// Во вторую очередь собираются однородные условия Неймана.
						for (integer  i = 0; i < f.maxbound; i++) {
							// однородные условия Неймана третий параметр равен false
							my_elmatr_quad_dissipation_rate_epsilon_3D_bound_standart_k_epsilon(i, f.maxelm,
								false, f.border_neighbor, ls, lw, w,
								f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
								f.pa, f.nvtx, f.prop_b, f.prop,
								f.potent, f.rdistWall, f.rdistWallmax);
						}

						// Схема для аппроксимации конвективного члена в модели 
						// Спаларта Аллмареса такая же как и для аппроксимации компонент 
						// вектора скорости.

						// Выбор делается из графического интерфейса пользователя !
						imyscheme = iFLOWScheme; // 31 07 2015

												 // Сборка строк матрицы для внутренних КО.
						for (integer i = 0; i < f.maxelm; i++) {
							my_elmatr_quad_turbulent_dissipation_rate_epsilon_Standart_KE_3D(
								i,
								f.border_neighbor,
								lw, ls,
								f.slau,
								f.slau_bon,
								//doublereal** diag_coef,
								//integer iVar,
								//bool btimedep,
								//doublereal tauparam,
								f.ptr,
								f.nvtx,
								f.potent,
								//doublereal* potent_temper,
								f.pa,
								f.prop,
								f.prop_b,
								f.maxelm,
								f.neighbors_for_the_internal_node,
								//doublereal* alpha,
								//doublereal dgx,
								//doublereal dgy,
								//doublereal dgz,
								//doublereal dbeta, 
								imyscheme,
								//bool bBussineskApproach,
								//doublereal temp_ref,
								//bool bfirst_start,
								//doublereal RCh,
								//integer iflowregime,
								//doublereal* speedoldtimestep,
								//doublereal* toldtimestep,
								b, lb, matlist,
								f.mf,
								//bool bVERYStable,
								//doublereal &sumanb,
								t.ilevel_alice,
								f.rdistWall,
								f.SInvariantStrainRateTensor
							);
						}


						break;


		default: if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT)) {

			// симметризация СЛАУ
			// Теперь симметризация матрицы производится внутри
			// сборки для каждого контрольного объёма отдельно.
			// Это позволило значительно ускорить вычисления, 
			// путём избавления от двойного цикла 
			// который был реализован раньше и сильно замедлял вычисления.


		   // теперь реализован стационарный и нестационарный солвер,
		   // в зависимости от передаваемого параметра btimedep.
		   // Значение tauparam содержит значение шага по времени.

			// Симметризация СЛАУ выполняется внутри
			// сборки матрицы:


#ifdef _OPENMP

			if (bparallelizm_old) {
				printf("error bparallelizm_old\n");
				system("pause");
			if (inumcore == 1) {
				// Граничные условия Дирихле обязательно 
				// должны собираться в первую очередь
				for (integer i = 0; i < f.maxbound; i++) {
					// условия Дирихле третий параметр равен true
					my_elmatr_quad_F3D_bound(i, f.maxelm,
						true,
						f.border_neighbor,
						ls, lw, w,
						iVar,
						f.slau_bon[iVar],
						dbeta, f.pa,
						f.nvtx,
						f.prop_b,
						f.prop,
						f.potent,
						f.iflowregime);
				}

			}
			
			if (inumcore == 2) {
				if (nd.b0.active) {

#pragma omp parallel
					{
#pragma omp sections
						{
#pragma omp section 
							{
								// первый поток
								for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
									integer iPloc = f.ifrontregulationgl[iscan_par];
									if (iPloc >= f.maxelm) {
										// Граничные условия Дирихле обязательно 
										// должны собираться в первую очередь
										// условия Дирихле третий параметр равен true
										my_elmatr_quad_F3D_bound(iPloc - f.maxelm, f.maxelm,
											true,
											f.border_neighbor,
											ls, lw, w,
											iVar,
											f.slau_bon[iVar],
											dbeta, f.pa,
											f.nvtx,
											f.prop_b,
											f.prop,
											f.potent,
											f.iflowregime);
									} // if

								} // for
							} // omp section
#pragma omp section
							{
								// второй поток
								for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
									integer iPloc = f.ifrontregulationgl[iscan_par];
									if (iPloc >= f.maxelm) {
										// Граничные условия Дирихле обязательно 
										// должны собираться в первую очередь
										// условия Дирихле третий параметр равен true
										my_elmatr_quad_F3D_bound(iPloc - f.maxelm, f.maxelm,
											true,
											f.border_neighbor,
											ls, lw, w,
											iVar,
											f.slau_bon[iVar],
											dbeta, f.pa,
											f.nvtx,
											f.prop_b,
											f.prop,
											f.potent,
											f.iflowregime);
									} // if
								} // for
							} // omp section
						}// omp sections
					} // omp parallel
					// серийный смыкающий кусок
					for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc = f.ifrontregulationgl[iscan_par];
						if (iPloc >= f.maxelm) {
							// Граничные условия Дирихле обязательно 
							// должны собираться в первую очередь
							// условия Дирихле третий параметр равен true
							my_elmatr_quad_F3D_bound(iPloc - f.maxelm, f.maxelm,
								true,
								f.border_neighbor,
								ls, lw, w,
								iVar,
								f.slau_bon[iVar],
								dbeta, f.pa,
								f.nvtx,
								f.prop_b,
								f.prop,
								f.potent,
								f.iflowregime);
						}
					}


				}
			}
		}
			else {
				// Свободный openmp без задания предварительно inumcore
				// Граничные условия Дирихле обязательно 
				// должны собираться в первую очередь
#pragma omp parallel for
				for (integer i = 0; i<f.maxbound; i++) {
					// условия Дирихле третий параметр равен true
					my_elmatr_quad_F3D_bound(i, f.maxelm,
						true,
						f.border_neighbor,
						ls, lw, w,
						iVar,
						f.slau_bon[iVar],
						dbeta, f.pa,
						f.nvtx,
						f.prop_b,
						f.prop,
						f.potent,
						f.iflowregime);
				}

			}

#else

			          // Граничные условия Дирихле обязательно 
			          // должны собираться в первую очередь
                      for (integer i=0; i<f.maxbound; i++) {
						  // условия Дирихле третий параметр равен true
						  my_elmatr_quad_F3D_bound(i, f.maxelm,
							                       true,
												   f.border_neighbor, 
												   ls, lw, w,
												   iVar, 
												   f.slau_bon[iVar],
												   dbeta, f.pa, 
												   f.nvtx,
												   f.prop_b, 
												   f.prop,
												   f.potent,
												   f.iflowregime);
					  }

#endif


					 

					  #ifdef _OPENMP

					  if (bparallelizm_old) {
						  printf("error bparallelizm_old\n");
						  system("pause");

						  if (inumcore == 1) {
							  // Во вторую очередь собираются однородные условия Неймана.
							  for (integer i = 0; i < f.maxbound; i++) {
								  // однородные условия Неймана третий параметр равен false
								  my_elmatr_quad_F3D_bound(i, f.maxelm,
									  false,
									  f.border_neighbor,
									  ls, lw, w, iVar,
									  f.slau_bon[iVar],
									  dbeta, f.pa,
									  f.nvtx,
									  f.prop_b,
									  f.prop,
									  f.potent,
									  f.iflowregime);
							  }
						  }

						  if (inumcore == 2) {
							  if (nd.b0.active) {


#pragma omp parallel
								  {
#pragma omp sections
									  {
#pragma omp section 
										  {
											  // первый поток
											  for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
												  integer iPloc = f.ifrontregulationgl[iscan_par];
												  if (iPloc >= f.maxelm) {
													  // Во вторую очередь собираются однородные условия Неймана.
													  // однородные условия Неймана третий параметр равен false
													  my_elmatr_quad_F3D_bound(iPloc - f.maxelm, f.maxelm,
														  false,
														  f.border_neighbor,
														  ls, lw, w, iVar,
														  f.slau_bon[iVar],
														  dbeta, f.pa,
														  f.nvtx,
														  f.prop_b,
														  f.prop,
														  f.potent,
														  f.iflowregime);
												  } // if

											  } // for
										  } // omp section
#pragma omp section 
										  {
											  // второй поток
											  for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
												  integer iPloc = f.ifrontregulationgl[iscan_par];
												  if (iPloc >= f.maxelm) {
													  // Во вторую очередь собираются однородные условия Неймана.
													  // однородные условия Неймана третий параметр равен false
													  my_elmatr_quad_F3D_bound(iPloc - f.maxelm, f.maxelm,
														  false,
														  f.border_neighbor,
														  ls, lw, w, iVar,
														  f.slau_bon[iVar],
														  dbeta, f.pa,
														  f.nvtx,
														  f.prop_b,
														  f.prop,
														  f.potent,
														  f.iflowregime);
												  } // if


											  } // for
										  } // omp section
									  } // omp sections
								  } // omp parallel
								  // серийный смыкающий кусок
								  for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
									  integer iPloc = f.ifrontregulationgl[iscan_par];
									  if (iPloc >= f.maxelm) {
										  // Во вторую очередь собираются однородные условия Неймана.
										  // однородные условия Неймана третий параметр равен false
										  my_elmatr_quad_F3D_bound(iPloc - f.maxelm, f.maxelm,
											  false,
											  f.border_neighbor,
											  ls, lw, w, iVar,
											  f.slau_bon[iVar],
											  dbeta, f.pa,
											  f.nvtx,
											  f.prop_b,
											  f.prop,
											  f.potent,
											  f.iflowregime);
									  }


								  }


							  }
						  }
					  }
					  else {
						  // Во вторую очередь собираются однородные условия Неймана.
#pragma omp parallel for
						  for (integer i = 0; i<f.maxbound; i++) {
							  // однородные условия Неймана третий параметр равен false
							  my_elmatr_quad_F3D_bound(i, f.maxelm,
								  false,
								  f.border_neighbor,
								  ls, lw, w, iVar,
								  f.slau_bon[iVar],
								  dbeta, f.pa,
								  f.nvtx,
								  f.prop_b,
								  f.prop,
								  f.potent,
								  f.iflowregime);
						  }
					  }

#else

					  // Во вторую очередь собираются однородные условия Неймана.
                      for (integer i=0; i<f.maxbound; i++) {
						  // однородные условия Неймана третий параметр равен false
						  my_elmatr_quad_F3D_bound(i, f.maxelm, 
							                       false, 
												   f.border_neighbor,
												   ls, lw, w, iVar,
												   f.slau_bon[iVar],
												   dbeta, f.pa,
												   f.nvtx,
												   f.prop_b, 
												   f.prop, 
												   f.potent,
												   f.iflowregime);
					  }

					  

#endif
#pragma omp parallel for
					  for (integer i = 0; i < f.maxbound; i++) {
						  if (f.slau_bon[iVar][i].aw < 0.0) {
							  printf("maxbound=%lld i=%lld problem iVar==%lld\n",f.maxbound,i, iVar);
							  system("pause");
							  //f.slau_bon[iVar][i].aw
						  }
					  }

					  imyscheme=UDS;
					  if (bVERYStable) {
						  // стабильная противопоточная схема первого порядка
						  imyscheme=UDS; // UDS BULG
					  }
					  else {
						  imyscheme=UNEVEN_WACEB;
						  //imyscheme=UDS; 
					  }
					  imyscheme=UDS; //UNEVEN_WACEB;
					  imyscheme=UNEVEN_WACEB; // 30 07 2015.
					  //imyscheme=UNEVENQUICK;// не работает на 30 07 2015.
					  // Выбор делается из графического интерфейса пользователя !
					  imyscheme=iFLOWScheme; // 31 07 2015

#ifdef _OPENMP

					  if (bparallelizm_old) {
						  printf("error bparallelizm_old\n");
						  system("pause");
					  if (inumcore == 1) {
						  doublereal temp_ref = 0.0;
						  temp_ref = f.OpTemp; // Operating Temperature
											   // Определяем минимальное значение температуры в расчётной области:
						  doublereal tavg = 0.0; // Минимальная температура в расчётной области.
						  for (integer i_1 = 0; i_1 < f.maxelm; i_1++) {
							  //if (potent_temper[f.ptr[i_1]] < tmin) {
							  tavg += t.potent[f.ptr[i_1]];
							  //}
						  }
						  tavg = tavg / f.maxelm;
						  // Архимедова сила встлытия на основе tmin всегода противоположно направлена силе тяжести.
						  temp_ref = tavg;

						  // Сборка строк матрицы для внутренних КО.
						  for (integer i = 0; i<f.maxelm; i++) {
							  TOCHKA p;
							  center_cord3D(i, f.nvtx, f.pa, p,100);
							  integer ib; // номер искомого блока
							  in_model_flow(p, ib, b, lb);
							  
							  bool bBussineskApproach = false;
							  bBussineskApproach = matlist[b[ib].imatid].bBussineskApproach;
							  
							 

							  if ((speedoldtimestep == nullptr) && (toldtimestep == nullptr)) {
								  // сборка матрицы для компонент скорости
								  my_elmatr_quad_F3D(i, f.border_neighbor, lw, ls, f.slau,
									  f.slau_bon,
									  f.diag_coef,
									  iVar,
									  btimedep, // стационарный или нестационарный
									  tauparam, // размер шага по времени
									  f.ptr,
									  f.nvtx,
									  f.potent,
									  t.potent,
									  f.pa, f.prop,
									  f.prop_b,
									  f.maxelm,
									  f.neighbors_for_the_internal_node,
									  f.alpha,
									  dgx, dgy, dgz,
									  dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
									  bBussineskApproach, temp_ref, // Bussinesk Approach
									  bfirst_start,
									  RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
									  f.iflowregime, // режим течения Ламинарный или Турбулентный
									  nullptr, // компонента скорости с предыдущего временного слоя
									  nullptr, // температура с предыдущего слоя по времени.
									  b, lb, matlist,
									  f.mf, // массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
									  bVERYStable,
									  sumanb[iVar][i],
									  t.ilevel_alice);
							  }
							  else
							  {
								  

								  // сборка матрицы для компонент скорости
								  my_elmatr_quad_F3D(i, f.border_neighbor, lw, ls, f.slau,
									  f.slau_bon,
									  f.diag_coef,
									  iVar,
									  btimedep, // стационарный или нестационарный
									  tauparam, // размер шага по времени
									  f.ptr,
									  f.nvtx,
									  f.potent,
									  t.potent,
									  f.pa, f.prop,
									  f.prop_b,
									  f.maxelm,
									  f.neighbors_for_the_internal_node,
									  f.alpha,
									  dgx, dgy, dgz,
									  dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
									  bBussineskApproach, temp_ref, // Bussinesk Approach
									  bfirst_start,
									  RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
									  f.iflowregime, // режим течения Ламинарный или Турбулентный
									  speedoldtimestep[iVar], // компонента скорости с предыдущего временного слоя
									  toldtimestep, // температура с предыдущего слоя по времени.
									  b, lb, matlist,
									  f.mf, bVERYStable,// массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
									  sumanb[iVar][i],
									  t.ilevel_alice);
							  }
						  }
					  
						  for (integer i = 0; i < f.maxbound; i++) {
							  sumanb[iVar][f.maxelm + i] = f.slau_bon[iVar][i].aw;
						  }

						  // Вычисление действия градиента давления.
						  // Это необходимо сделать строго после того как все
						  // sumanb вычислены.
						  for (integer  i = 0; i < f.maxelm; i++) {
							  pterm(i, f.slau,
								  f.potent, iVar,
								  f.maxelm, f.nvtx, f.pa, 
								  sumanb);
						  }
					  }

			if (inumcore==2) {
				if (nd.b0.active) {


					doublereal temp_ref = 0.0;
					temp_ref = f.OpTemp; // Operating Temperature
										 // Определяем минимальное значение температуры в расчётной области:
					//doublereal tavg = 0.0; // Минимальная температура в расчётной области.
					//for (integer i_1 = 0; i_1 < f.maxelm; i_1++) {
						//if (potent_temper[f.ptr[i_1]] < tmin) {
						//tavg += t.potent[f.ptr[i_1]];
						//}
					//}
					//tavg = tavg / f.maxelm;
					// Архимедова сила встлытия на основе tmin всегода противоположно направлена силе тяжести.
					//temp_ref = tavg;



					#pragma omp parallel
					{
#pragma omp sections
						{
#pragma omp section 
							{
					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=f.ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							TOCHKA p;
						  center_cord3D(iPloc,f.nvtx,f.pa,p,100);
						  integer ib; // номер искомого блока
						  in_model_flow(p,ib,b,lb);
						  
						  bool bBussineskApproach=false;
						  bBussineskApproach=matlist[b[ib].imatid].bBussineskApproach;
						

						  if ((speedoldtimestep==nullptr)&&(toldtimestep==nullptr)) {
							  // сборка матрицы для компонент скорости
						  my_elmatr_quad_F3D(iPloc, f.border_neighbor, lw,ls, f.slau,
							                 f.slau_bon,
											 f.diag_coef, 
											 iVar, 
											 btimedep, // стационарный или нестационарный
											 tauparam, // размер шага по времени
											 f.ptr,
											 f.nvtx,
											 f.potent,
											 t.potent,
											 f.pa, f.prop,
											 f.prop_b,
											 f.maxelm,
											 f.neighbors_for_the_internal_node, 
											 f.alpha,
											 dgx, dgy, dgz,
											 dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
											 bBussineskApproach, temp_ref, // Bussinesk Approach
											 bfirst_start,
											 RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
											 f.iflowregime, // режим течения Ламинарный или Турбулентный
											 nullptr, // компонента скорости с предыдущего временного слоя
											 nullptr, // температура с предыдущего слоя по времени.
											 b,lb, matlist,
											 f.mf, // массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
											 bVERYStable,
											 sumanb[iVar][iPloc],
							                 t.ilevel_alice);
						  }
						  else
						  {

						  // сборка матрицы для компонент скорости
						  my_elmatr_quad_F3D(iPloc, f.border_neighbor, lw,ls, f.slau,
							                 f.slau_bon,
											 f.diag_coef, 
											 iVar, 
											 btimedep, // стационарный или нестационарный
											 tauparam, // размер шага по времени
											 f.ptr,
											 f.nvtx,
											 f.potent,
											 t.potent,
											 f.pa, f.prop,
											 f.prop_b,
											 f.maxelm,
											 f.neighbors_for_the_internal_node, 
											 f.alpha,
											 dgx, dgy, dgz,
											 dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
											 bBussineskApproach, temp_ref, // Bussinesk Approach
											 bfirst_start,
											 RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
											 f.iflowregime, // режим течения Ламинарный или Турбулентный
											 speedoldtimestep[iVar], // компонента скорости с предыдущего временного слоя
											 toldtimestep, // температура с предыдущего слоя по времени.
											 b,lb, matlist,
											 f.mf,bVERYStable,// массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
											 sumanb[iVar][iPloc],
							                 t.ilevel_alice);
						  }
				      }
						}
							}
#pragma omp section
							{
					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=f.ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							TOCHKA p;
						  center_cord3D(iPloc,f.nvtx,f.pa,p,100);
						  integer ib; // номер искомого блока
						  in_model_flow(p,ib,b,lb);
						  
						  bool bBussineskApproach=false;
						  bBussineskApproach=matlist[b[ib].imatid].bBussineskApproach;
						 

						  if ((speedoldtimestep==nullptr)&&(toldtimestep==nullptr)) {
							  // сборка матрицы для компонент скорости
						  my_elmatr_quad_F3D(iPloc, f.border_neighbor, lw,ls, f.slau,
							                 f.slau_bon,
											 f.diag_coef, 
											 iVar, 
											 btimedep, // стационарный или нестационарный
											 tauparam, // размер шага по времени
											 f.ptr,
											 f.nvtx,
											 f.potent,
											 t.potent,
											 f.pa, f.prop,
											 f.prop_b,
											 f.maxelm,
											 f.neighbors_for_the_internal_node, 
											 f.alpha,
											 dgx, dgy, dgz,
											 dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
											 bBussineskApproach, temp_ref, // Bussinesk Approach
											 bfirst_start,
											 RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
											 f.iflowregime, // режим течения Ламинарный или Турбулентный
											 nullptr, // компонента скорости с предыдущего временного слоя
											 nullptr, // температура с предыдущего слоя по времени.
											 b,lb, matlist,
											 f.mf, // массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
											 bVERYStable,
											 sumanb[iVar][iPloc],
							                 t.ilevel_alice);
						  }
						  else
						  {

						  // сборка матрицы для компонент скорости
						  my_elmatr_quad_F3D(iPloc, f.border_neighbor, lw,ls, f.slau,
							                 f.slau_bon,
											 f.diag_coef, 
											 iVar, 
											 btimedep, // стационарный или нестационарный
											 tauparam, // размер шага по времени
											 f.ptr,
											 f.nvtx,
											 f.potent,
											 t.potent,
											 f.pa, f.prop,
											 f.prop_b,
											 f.maxelm,
											 f.neighbors_for_the_internal_node, 
											 f.alpha,
											 dgx, dgy, dgz,
											 dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
											 bBussineskApproach, temp_ref, // Bussinesk Approach
											 bfirst_start,
											 RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
											 f.iflowregime, // режим течения Ламинарный или Турбулентный
											 speedoldtimestep[iVar], // компонента скорости с предыдущего временного слоя
											 toldtimestep, // температура с предыдущего слоя по времени.
											 b,lb, matlist,
											 f.mf,bVERYStable,// массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
											 sumanb[iVar][iPloc],
							                 t.ilevel_alice);
						  }
				      }
						}

							} // omp section
							} // omp sections
							} // omp parallel

					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=f.ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
TOCHKA p;
						  center_cord3D(iPloc,f.nvtx,f.pa,p,100);
						  integer ib; // номер искомого блока
						  in_model_flow(p,ib,b,lb);
						  
						  bool bBussineskApproach=false;
						  bBussineskApproach=matlist[b[ib].imatid].bBussineskApproach;
						 

						  if ((speedoldtimestep==nullptr)&&(toldtimestep==nullptr)) {
							  // сборка матрицы для компонент скорости
						  my_elmatr_quad_F3D(iPloc, f.border_neighbor, lw,ls, f.slau,
							                 f.slau_bon,
											 f.diag_coef, 
											 iVar, 
											 btimedep, // стационарный или нестационарный
											 tauparam, // размер шага по времени
											 f.ptr,
											 f.nvtx,
											 f.potent,
											 t.potent,
											 f.pa, f.prop,
											 f.prop_b,
											 f.maxelm,
											 f.neighbors_for_the_internal_node, 
											 f.alpha,
											 dgx, dgy, dgz,
											 dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
											 bBussineskApproach, temp_ref, // Bussinesk Approach
											 bfirst_start,
											 RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
											 f.iflowregime, // режим течения Ламинарный или Турбулентный
											 nullptr, // компонента скорости с предыдущего временного слоя
											 nullptr, // температура с предыдущего слоя по времени.
											 b,lb, matlist,
											 f.mf, // массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
											 bVERYStable,
											 sumanb[iVar][iPloc],
							                 t.ilevel_alice);
						  }
						  else
						  {

						  // сборка матрицы для компонент скорости
						  my_elmatr_quad_F3D(iPloc, f.border_neighbor, lw,ls, f.slau,
							                 f.slau_bon,
											 f.diag_coef, 
											 iVar, 
											 btimedep, // стационарный или нестационарный
											 tauparam, // размер шага по времени
											 f.ptr,
											 f.nvtx,
											 f.potent,
											 t.potent,
											 f.pa, f.prop,
											 f.prop_b,
											 f.maxelm,
											 f.neighbors_for_the_internal_node, 
											 f.alpha,
											 dgx, dgy, dgz,
											 dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
											 bBussineskApproach, temp_ref, // Bussinesk Approach
											 bfirst_start,
											 RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
											 f.iflowregime, // режим течения Ламинарный или Турбулентный
											 speedoldtimestep[iVar], // компонента скорости с предыдущего временного слоя
											 toldtimestep, // температура с предыдущего слоя по времени.
											 b,lb, matlist,
											 f.mf,bVERYStable,// массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
											 sumanb[iVar][iPloc],
							                 t.ilevel_alice);
						  }
				      }
						}

					}


					for (integer i = 0; i < f.maxbound; i++) {
						sumanb[iVar][f.maxelm + i] = f.slau_bon[iVar][i].aw;
					}

					// Вычисление действия градиента давления.
				    // Это необходимо сделать строго после того как все
					// sumanb вычислены.
					for (integer i = 0; i < f.maxelm; i++) {
						pterm(i, f.slau,
							f.potent, iVar,
							f.maxelm, f.nvtx, f.pa,
							sumanb);
					}
				}

				}
					  else {
						  doublereal temp_ref = 0.0;
						  temp_ref = f.OpTemp; // Operating Temperature
											   // Определяем минимальное значение температуры в расчётной области:
						  //doublereal tavg = 0.0; // Минимальная температура в расчётной области.
						  //for (integer i_1 = 0; i_1 < f.maxelm; i_1++) {
							  //if (potent_temper[f.ptr[i_1]] < tmin) {
							//  tavg += t.potent[f.ptr[i_1]];
							  //}
						  //}
						  //tavg = tavg / f.maxelm;
						  // Архимедова сила встлытия на основе tmin всегода противоположно направлена силе тяжести.
						 // temp_ref = tavg;

						  // Сборка строк матрицы для внутренних КО.
#pragma omp parallel for
						  for (integer i = 0; i<f.maxelm; i++) {
							  TOCHKA p;
							  center_cord3D(i, f.nvtx, f.pa, p, 100);
							  integer ib; // номер искомого блока
							  in_model_flow(p, ib, b, lb);

							  bool bBussineskApproach = false;
							  bBussineskApproach = matlist[b[ib].imatid].bBussineskApproach;


							  if ((speedoldtimestep == nullptr) && (toldtimestep == nullptr)) {
								  // сборка матрицы для компонент скорости
								  my_elmatr_quad_F3D(i, f.border_neighbor, lw, ls, f.slau,
									  f.slau_bon,
									  f.diag_coef,
									  iVar,
									  btimedep, // стационарный или нестационарный
									  tauparam, // размер шага по времени
									  f.ptr,
									  f.nvtx,
									  f.potent,
									  t.potent,
									  f.pa, f.prop,
									  f.prop_b,
									  f.maxelm,
									  f.neighbors_for_the_internal_node,
									  f.alpha,
									  dgx, dgy, dgz,
									  dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
									  bBussineskApproach, temp_ref, // Bussinesk Approach
									  bfirst_start,
									  RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
									  f.iflowregime, // режим течения Ламинарный или Турбулентный
									  nullptr, // компонента скорости с предыдущего временного слоя
									  nullptr, // температура с предыдущего слоя по времени.
									  b, lb, matlist,
									  f.mf, // массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
									  bVERYStable,
									  sumanb[iVar][i],
									  t.ilevel_alice);
							  }
							  else
							  {

								  // сборка матрицы для компонент скорости
								  my_elmatr_quad_F3D(i, f.border_neighbor, lw, ls, f.slau,
									  f.slau_bon,
									  f.diag_coef,
									  iVar,
									  btimedep, // стационарный или нестационарный
									  tauparam, // размер шага по времени
									  f.ptr,
									  f.nvtx,
									  f.potent,
									  t.potent,
									  f.pa, f.prop,
									  f.prop_b,
									  f.maxelm,
									  f.neighbors_for_the_internal_node,
									  f.alpha,
									  dgx, dgy, dgz,
									  dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
									  bBussineskApproach, temp_ref, // Bussinesk Approach
									  bfirst_start,
									  RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
									  f.iflowregime, // режим течения Ламинарный или Турбулентный
									  speedoldtimestep[iVar], // компонента скорости с предыдущего временного слоя
									  toldtimestep, // температура с предыдущего слоя по времени.
									  b, lb, matlist,
									  f.mf, bVERYStable,// массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
									  sumanb[iVar][i],
									  t.ilevel_alice);
							  }
						  }
					  
#pragma omp parallel for
						  for (integer i = 0; i < f.maxbound; i++) {
							  sumanb[iVar][f.maxelm + i] = f.slau_bon[iVar][i].aw;
						  }

						  // Вычисление действия градиента давления.
						// Это необходимо сделать строго после того как все
						   // sumanb вычислены.
#pragma omp parallel for
						  for (integer i = 0; i < f.maxelm; i++) {
							  pterm(i, f.slau,
								  f.potent, iVar,
								  f.maxelm, f.nvtx, f.pa,
								  sumanb);
						  }
}

#else
//+

                      doublereal temp_ref = 0.0;
                      temp_ref = f.OpTemp; // Operating Temperature
					  // Определяем минимальное значение температуры в расчётной области:
                      //doublereal tavg = 0.0; // Минимальная температура в расчётной области.
                      //for (integer i_1 = 0; i_1 < f.maxelm; i_1++) {
	                       //if (potent_temper[f.ptr[i_1]] < tmin) {
	                    //  tavg += t.potent[f.ptr[i_1]];
	                      //}
                      //}
                      //tavg = tavg / f.maxelm;
                      // Архимедова сила встлытия на основе tmin всегода противоположно направлена силе тяжести.
                    //  temp_ref = tavg;

                      // Сборка строк матрицы для внутренних КО.
				      for (integer i=0; i<f.maxelm; i++) { 
						  TOCHKA p;
						  center_cord3D(i,f.nvtx,f.pa,p,100);
						  integer ib; // номер искомого блока
						  in_model_flow(p,ib,b,lb);
						  
						  bool bBussineskApproach=false;
						  bBussineskApproach=matlist[b[ib].imatid].bBussineskApproach;
						  

						  if ((speedoldtimestep==nullptr)&&(toldtimestep==nullptr)) {
							  // сборка матрицы для компонент скорости
						  my_elmatr_quad_F3D(i, f.border_neighbor, lw,ls, f.slau,
							                 f.slau_bon,
											 f.diag_coef, 
											 iVar, 
											 btimedep, // стационарный или нестационарный
											 tauparam, // размер шага по времени
											 f.ptr,
											 f.nvtx,
											 f.potent,
											 t.potent,
											 f.pa, f.prop,
											 f.prop_b,
											 f.maxelm,
											 f.neighbors_for_the_internal_node, 
											 f.alpha,
											 dgx, dgy, dgz,
											 dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
											 bBussineskApproach, temp_ref, // Bussinesk Approach
											 bfirst_start,
											 RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
											 f.iflowregime, // режим течения Ламинарный или Турбулентный
											 nullptr, // компонента скорости с предыдущего временного слоя
											 nullptr, // температура с предыдущего слоя по времени.
											 b,lb, matlist,
											 f.mf, // массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
											 bVERYStable,
											 sumanb[iVar][i],
							                 t.ilevel_alice);
						  }
						  else
						  {

						  // сборка матрицы для компонент скорости
						  my_elmatr_quad_F3D(i, f.border_neighbor, lw,ls, f.slau,
							                 f.slau_bon,
											 f.diag_coef, 
											 iVar, 
											 btimedep, // стационарный или нестационарный
											 tauparam, // размер шага по времени
											 f.ptr,
											 f.nvtx,
											 f.potent,
											 t.potent,
											 f.pa, f.prop,
											 f.prop_b,
											 f.maxelm,
											 f.neighbors_for_the_internal_node, 
											 f.alpha,
											 dgx, dgy, dgz,
											 dbeta, imyscheme, // BULG EXP2 QUICK UNEVENQUICK UNEVEN_WACEB
											 bBussineskApproach, temp_ref, // Bussinesk Approach
											 bfirst_start,
											 RCh, // 0.0 RCh // на примере задачи Блазиуса подстановка RCh отличной от нуля приводит к расходимости.
											 f.iflowregime, // режим течения Ламинарный или Турбулентный
											 speedoldtimestep[iVar], // компонента скорости с предыдущего временного слоя
											 toldtimestep, // температура с предыдущего слоя по времени.
											 b,lb, matlist,
											 f.mf,bVERYStable,// массовый поток через грани контрольного объёма с учётом монотонизирующей поправки Рхи-Чоу.
											 sumanb[iVar][i],
							                 t.ilevel_alice);
						  }
				      }

					  for (integer i = 0; i < f.maxbound; i++) {
						  sumanb[iVar][f.maxelm + i] = f.slau_bon[iVar][i].aw;
					  }

					  // Вычисление действия градиента давления.
					  // Это необходимо сделать строго после того как все
					  // sumanb вычислены.
					  for (integer i = 0; i < f.maxelm; i++) {
						  pterm(i, f.slau,
							  f.potent, iVar,
							  f.maxelm, f.nvtx, f.pa,
							  sumanb);
					  }

#endif
#pragma omp parallel for
					  for (integer  i = 0; i < f.maxelm; i++) {
						  if (f.slau[iVar][i].ap < 0.0) {
							  printf("maxelm=%lld i=%lld problem iVar==%lld ap=%e\n", f.maxelm, i, iVar, f.slau[iVar][i].ap);
							  system("pause");
							  //f.slau_bon[iVar][i].aw
						  }
					  }
			          res=0.0;
					  //printf("assemblage finish...\n"); // diagnostic message
					  //getchar();
				  }

#if doubleintprecision == 1
				  /* // отладочная информация.
				  // несоответствующие физическому смыслу осцилляции вблизи выходной границы. виновато deltaF при сборке матрицы.
				  // следуя Патанкару от deltaF вообще следует отказаться.
				  if (iVar==VY) {
						integer id6[6]={17325,17347,17369,17391,17413,35495}; // последовательность узлов в линию вблизи выходной границы
						// debug
						if (inumiter==2) {
							for ( i=0; i<5; i++) {
								printf("id=%lld ap=%e ae=%e aw=%e an=%e as=%e at=%e ab=%e b=%e\n",id6[i], f.slau[VY][id6[i]].ap, f.slau[VY][id6[i]].ae, f.slau[VY][id6[i]].aw, f.slau[VY][id6[i]].an, f.slau[VY][id6[i]].as, f.slau[VY][id6[i]].at, f.slau[VY][id6[i]].ab, f.slau[VY][id6[i]].b);
								printf("mf[SSIDE]=%e, mf[NSIDE]=%e\n",f.mf[id6[i]][SSIDE],f.mf[id6[i]][NSIDE]);

							}
							printf("id=%lld aw=%e ai=%e b=%e\n",id6[5], f.slau_bon[VY][id6[i]-f.maxelm].aw, f.slau_bon[VY][id6[i]-f.maxelm].ai,f.slau_bon[VY][id6[i]-f.maxelm].b);
							getchar();
					}

				  }
				  */
#else
				  /* // отладочная информация.
				  // несоответствующие физическому смыслу осцилляции вблизи выходной границы. виновато deltaF при сборке матрицы.
				  // следуя Патанкару от deltaF вообще следует отказаться.
				  if (iVar==VY) {
						integer id6[6]={17325,17347,17369,17391,17413,35495}; // последовательность узлов в линию вблизи выходной границы
						// debug
						if (inumiter==2) {
							for ( i=0; i<5; i++) {
								printf("id=%lld ap=%e ae=%e aw=%e an=%e as=%e at=%e ab=%e b=%e\n",id6[i], f.slau[VY][id6[i]].ap, f.slau[VY][id6[i]].ae, f.slau[VY][id6[i]].aw, f.slau[VY][id6[i]].an, f.slau[VY][id6[i]].as, f.slau[VY][id6[i]].at, f.slau[VY][id6[i]].ab, f.slau[VY][id6[i]].b);
								printf("mf[SSIDE]=%e, mf[NSIDE]=%e\n",f.mf[id6[i]][SSIDE],f.mf[id6[i]][NSIDE]);

							}
							printf("id=%lld aw=%e ai=%e b=%e\n",id6[5], f.slau_bon[VY][id6[i]-f.maxelm].aw, f.slau_bon[VY][id6[i]-f.maxelm].ai,f.slau_bon[VY][id6[i]-f.maxelm].b);
							getchar();
						}

				  }
				  */
#endif
				  
			break;
	}

	
	/*
	// Для параллельного программирования очень подходит этот паттерн.
	#ifdef _OPENMP

			if (inumcore==2) {
				if (nd.b0.active) {

					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];

					}
					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];

					}
					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];

					}


				}
			}

#else
	*/

#ifdef _OPENMP 
	omp_set_num_threads(1); // установка числа потоков
#endif

	/*
    for (i=0; i<maxelm; i++) {
		printf("%.3f %.3f %.3f %.3f %.3f = %.3f \n", slau[iVar][i].ap, slau[iVar][i].ae, slau[iVar][i].an, slau[iVar][i].as, slau[iVar][i].aw, slau[iVar][i].b);
	}
	printf("\n");
	getchar();
	//*/
	if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar==NUSHA) || (iVar==PAM)||
		(iVar == TURBULENT_KINETIK_ENERGY)||(iVar== TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA)||
		(iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS)||(iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)) {
		for (integer i_1 = 0; i_1 < f.maxelm; i_1++) {
			integer iVar_in = iVar;
			if (iVar == NUSHA) iVar_in = NUSHA_SL;
			if (iVar == TURBULENT_KINETIK_ENERGY) iVar_in = TURBULENT_KINETIK_ENERGY_SL;
			if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) iVar_in = TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL;
			if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) iVar_in = TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL;
			if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) iVar_in = TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL;
			if (f.slau[iVar_in][i_1].b != f.slau[iVar_in][i_1].b) {
				switch (iVar) {
				case VELOCITY_X_COMPONENT: printf("VX problem\n");
					break;
				case VELOCITY_Y_COMPONENT: printf("VY problem\n");
					break;
				case VELOCITY_Z_COMPONENT: printf("VZ problem\n");
					break;
				case PAM: printf("PAM problem\n");
					break;
				case NUSHA: printf("NUSHA problem\n");
					break;
				case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY problem\n");
					break;
				case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA problem\n");
					break;
				case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_Standart K-EPSILON problem\n");
					break;
				case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_Standart K-EPSILON problem\n");
					break;
				case TEMP: printf("TEMP problem\n");
					break;
				}
				printf("POST ASSEMBLE CONTROL b part.\n");
				printf("NAN or INF in premeshin.txt file. Power in control volume= %lld is undefined...\n", i_1);
				printf("ispolzuite poslednuu versiu Mesh generator AliceMesh. 23.09.2018.\n");
				system("PAUSE");
				exit(1);
			}
		}
	}

	/*doublereal *rthdsd=nullptr; // правая часть системы уравнений
	if (iVar!=TEMP) {
		rthdsd = new doublereal[f.maxelm+f.maxbound];
	}
	else {
		rthdsd = new doublereal[t.maxelm+t.maxbound];
	}*/
	if (rthdsd==nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	SIMPLESPARSE sparseM; // разреженная матрица
	IMatrix sparseS; // разреженная матрица в формате IMatrix
	

	// Если эта переменная равна истине то мы не выделяем здесь никакой памяти,
	// Всё выделяется внутри BiCGStab или вне функции solve.
	bool bBiCGStabSaad=true;
    
    // разреженная матрица в формате CSIR
    //doublereal *adiag, *altr;
    //integer *jptr, *iptr;

	//switch (iVar) {
	//case VX: printf("VX \n"); break;
	//case VY: printf("VY \n"); break;
	//case VZ: printf("VZ \n"); break;
	//case NUSHA: printf("NU \n"); break;
	//case TURBULENT_KINETIK_ENERGY: printf("TURBULENT_KINETIK_ENERGY \n"); break;
	//case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA \n"); break;
	//case TURBULENT_KINETIK_ENERGY_STD_K_EPS: printf("TURBULENT_KINETIK_ENERGY_Standart K-EPSILON problem\n"); break;
	//case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: printf("TURBULENT_DISSIPATION_RATE_EPSILON_Standart K-EPSILON problem\n"); break;
	//case PAM: printf("PAM \n"); break;
	//case TEMP: printf("TEMP \n"); break;
	//}

	if (iVar!=TEMP) {

#pragma omp parallel for
		for (integer  i = 0; i < f.maxelm + f.maxbound; i++) {
			// инициализация.
			rthdsd[i] = 0.0;
		}

		// выделение памяти и инициализация для 
	    // простейшей разреженной матрицы.
		if (!bBiCGStabSaad) {
		    initsimplesparse(sparseM, f.maxelm + f.maxbound );
		    initIMatrix(&sparseS, f.maxelm + f.maxbound);
		}
		else {
			/*Lr1sk_up
			if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2) {
				// LR1sk
				initsimplesparse(sparseM, f.maxelm + f.maxbound);
				initIMatrix(&sparseS, f.maxelm + f.maxbound);
			}
			*/
		}

		// Без диагонального предобуславливателя. TODO

		// Диагональный предобуславливатель:

		    // Нижеследующего делать ни в коем случае нельзя по непонятной причине. TODO
		    // после его применения все элементы нужно сравнивать с единицей. 
			// 1.0 характерное значение.

			// Внутренние узлы:
            //for (i=0; i<f.maxelm; i++) {
                //f.slau[iVar][i].ae/=f.slau[iVar][i].ap;
                //f.slau[iVar][i].an/=f.slau[iVar][i].ap;
                //f.slau[iVar][i].at/=f.slau[iVar][i].ap;
                //f.slau[iVar][i].as/=f.slau[iVar][i].ap;
			    //f.slau[iVar][i].aw/=f.slau[iVar][i].ap;
			    //f.slau[iVar][i].ab/=f.slau[iVar][i].ap;
                //f.slau[iVar][i].b/=f.slau[iVar][i].ap;
                //f.slau[iVar][i].ap=1.0;
		    //}

			// Граничные узлы:
			//for (i=0; i<f.maxbound; i++) {
				//f.slau_bon[iVar][i].ai/=f.slau_bon[iVar][i].aw;
				//f.slau_bon[iVar][i].b/=f.slau_bon[iVar][i].aw;
				//f.slau_bon[iVar][i].aw=1.0;
			//}
			
		// В принципе можно было бы ввести верхнюю релаксацию для поправки давления,
		// но непонятно к какому значению нужно релаксировать. Можно релаксировать к текущему
		// состоянию но это надо пересобирать матрицу для поправки давления после того как сделано несколько итераций.
		// Здесь правда можно ввести многоуровневую систему - организовать аппроксимацию высокого порядка
		// реализовав её постепенно через метод отложенной коррекции - это должно улучшить точность с которой решается уравнение для давления.

		// Замечание о том как правильно ввести нижнюю релаксацию.
		// Есть два способа и оба они модифицируют собранную матрицу СЛАУ. Причём можно использовать их только по отдельности.
		// Первый способ введение псевдо времени tau ~ alpha/ap/(1.0-alpha) что даёт apzero=(1-alpha)*ap/alpha. Видно что при alpha=1.0 что соответствует отсутствию релаксации apzero==0.0.
		// Второй способ модификация матрицы СЛАУ после её сборки: ap->ap/alpha; b->b+ap*(1.0-alpha)*Vpolditer/alpha.
		// Информация получена в письме Гаврилова Андрея.

		

		// Для внутренних узлов расчётной сетки:
#pragma omp parallel for
		for (integer i = 0; i < f.maxelm; i++) {
			switch (iVar) {
			case VELOCITY_X_COMPONENT: if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau[iVar][i].ap / f.alpha[iVar], f.slau[iVar][i].iP, f.slau[iVar][i].iP, true); }
					 if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iP, f.slau[iVar][i].ap / f.alpha[iVar]); }
#if doubleintprecision == 1
					 // printf("iP=%lld mm=%lld i=%lld iVar=%lld", f.slau[iVar][i].iP, f.maxelm,i,iVar);
#else
					 // printf("iP=%d mm=%d i=%d iVar=%d", f.slau[iVar][i].iP, f.maxelm,i,iVar);
#endif
					 
					 rthdsd[f.slau[iVar][i].iP] = f.slau[iVar][i].b + (1 - f.alpha[iVar])*f.slau[iVar][i].ap*f.potent[VXCOR][f.slau[iVar][i].iP] / f.alpha[iVar];
					 break;
			case VELOCITY_Y_COMPONENT: if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau[iVar][i].ap / f.alpha[iVar], f.slau[iVar][i].iP, f.slau[iVar][i].iP, true); }
					 if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iP, f.slau[iVar][i].ap / f.alpha[iVar]); }
					 rthdsd[f.slau[iVar][i].iP] = f.slau[iVar][i].b + (1 - f.alpha[iVar])*f.slau[iVar][i].ap*f.potent[VYCOR][f.slau[iVar][i].iP] / f.alpha[iVar];
					 break;
			case VELOCITY_Z_COMPONENT: if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau[iVar][i].ap / f.alpha[iVar], f.slau[iVar][i].iP, f.slau[iVar][i].iP, true); }
					 if ((!bBiCGStabSaad)/* || (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iP, f.slau[iVar][i].ap / f.alpha[iVar]); }
                         rthdsd[f.slau[iVar][i].iP]=f.slau[iVar][i].b+(1-f.alpha[iVar])*f.slau[iVar][i].ap*f.potent[VZCOR][f.slau[iVar][i].iP] / f.alpha[iVar];
				         break;
			case NUSHA:
				rthdsd[f.slau[NUSHA_SL][i].iP] = f.slau[NUSHA_SL][i].b + (1 - f.alpha[NUSHA_SL])*
					f.slau[NUSHA_SL][i].ap*f.potent[NUSHA][f.slau[NUSHA_SL][i].iP] / f.alpha[NUSHA_SL];
				break;
			case TURBULENT_KINETIK_ENERGY:
				rthdsd[f.slau[TURBULENT_KINETIK_ENERGY_SL][i].iP] = f.slau[TURBULENT_KINETIK_ENERGY_SL][i].b +
					(1 - f.alpha[TURBULENT_KINETIK_ENERGY_SL])*f.slau[TURBULENT_KINETIK_ENERGY_SL][i].ap*
					f.potent[TURBULENT_KINETIK_ENERGY][f.slau[TURBULENT_KINETIK_ENERGY_SL][i].iP] / f.alpha[TURBULENT_KINETIK_ENERGY_SL];
				break;
			case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
				rthdsd[f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].iP] = f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].b +
					(1 - f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL])*f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].ap*
					f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].iP] 
					/ f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL];
				break;
			case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
				rthdsd[f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL][i].iP] = f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL][i].b +
					(1 - f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL])*f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL][i].ap*
					f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL][i].iP] / f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL];
				break;
			case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
				rthdsd[f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL][i].iP] = f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL][i].b +
					(1 - f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL])*f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL][i].ap*
					f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL][i].iP]
					/ f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL];
				break;
			   case PAM: // PRESSURE:
				   if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau[iVar][i].ap, f.slau[iVar][i].iP, f.slau[iVar][i].iP, true); }
				   if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iP, f.slau[iVar][i].ap); }

				   rthdsd[f.slau[iVar][i].iP] = f.slau[iVar][i].b;
				 break;
		     } // switch

			 
		
			if ((!bBiCGStabSaad)/* || (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) {
		         const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

			     if ((f.slau[iVar][i].iE>-1) && (fabs(f.slau[iVar][i].ae) > nonzeroEPS)) {
				    addelmsimplesparse(sparseM, -f.slau[iVar][i].ae, f.slau[iVar][i].iP, f.slau[iVar][i].iE, true);
				    setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iE, -f.slau[iVar][i].ae);
			     }
			     if ((f.slau[iVar][i].iN>-1) && (fabs(f.slau[iVar][i].an) > nonzeroEPS)) {
				    addelmsimplesparse(sparseM, -f.slau[iVar][i].an, f.slau[iVar][i].iP, f.slau[iVar][i].iN, true);
				    setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iN, -f.slau[iVar][i].an);
			     }
			     if ((f.slau[iVar][i].iT>-1) && (fabs(f.slau[iVar][i].at) > nonzeroEPS)) {
				    addelmsimplesparse(sparseM, -f.slau[iVar][i].at, f.slau[iVar][i].iP, f.slau[iVar][i].iT, true);
                    setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iT, -f.slau[iVar][i].at);
			     }
			     if ((f.slau[iVar][i].iS>-1) && (fabs(f.slau[iVar][i].as) > nonzeroEPS)) {
				    addelmsimplesparse(sparseM, -f.slau[iVar][i].as, f.slau[iVar][i].iP, f.slau[iVar][i].iS, true);
                    setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iS, -f.slau[iVar][i].as);
			     }
			     if ((f.slau[iVar][i].iW>-1) && (fabs(f.slau[iVar][i].aw) > nonzeroEPS)) {
				    addelmsimplesparse(sparseM, -f.slau[iVar][i].aw, f.slau[iVar][i].iP, f.slau[iVar][i].iW, true);
                    setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iW, -f.slau[iVar][i].aw);
			     }
			     if ((f.slau[iVar][i].iB>-1) && (fabs(f.slau[iVar][i].ab) > nonzeroEPS)) {
				    addelmsimplesparse(sparseM, -f.slau[iVar][i].ab, f.slau[iVar][i].iP, f.slau[iVar][i].iB, true);
				    setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iB, -f.slau[iVar][i].ab);
			     }

				 if (b_on_adaptive_local_refinement_mesh) {

					 if ((f.slau[iVar][i].iE2 > -1) && (fabs(f.slau[iVar][i].ae2) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].ae2, f.slau[iVar][i].iP, f.slau[iVar][i].iE2, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iE2, -f.slau[iVar][i].ae2);
					 }
					 if ((f.slau[iVar][i].iN2 > -1) && (fabs(f.slau[iVar][i].an2) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].an2, f.slau[iVar][i].iP, f.slau[iVar][i].iN2, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iN2, -f.slau[iVar][i].an2);
					 }
					 if ((f.slau[iVar][i].iT2 > -1) && (fabs(f.slau[iVar][i].at2) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].at2, f.slau[iVar][i].iP, f.slau[iVar][i].iT2, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iT2, -f.slau[iVar][i].at2);
					 }
					 if ((f.slau[iVar][i].iS2 > -1) && (fabs(f.slau[iVar][i].as2) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].as2, f.slau[iVar][i].iP, f.slau[iVar][i].iS2, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iS2, -f.slau[iVar][i].as2);
					 }
					 if ((f.slau[iVar][i].iW2 > -1) && (fabs(f.slau[iVar][i].aw2) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].aw2, f.slau[iVar][i].iP, f.slau[iVar][i].iW2, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iW2, -f.slau[iVar][i].aw2);
					 }
					 if ((f.slau[iVar][i].iB2 > -1) && (fabs(f.slau[iVar][i].ab2) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].ab2, f.slau[iVar][i].iP, f.slau[iVar][i].iB2, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iB2, -f.slau[iVar][i].ab2);
					 }


					 if ((f.slau[iVar][i].iE3 > -1) && (fabs(f.slau[iVar][i].ae3) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].ae3, f.slau[iVar][i].iP, f.slau[iVar][i].iE3, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iE3, -f.slau[iVar][i].ae3);
					 }
					 if ((f.slau[iVar][i].iN3 > -1) && (fabs(f.slau[iVar][i].an3) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].an3, f.slau[iVar][i].iP, f.slau[iVar][i].iN3, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iN3, -f.slau[iVar][i].an3);
					 }
					 if ((f.slau[iVar][i].iT3 > -1) && (fabs(f.slau[iVar][i].at3) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].at3, f.slau[iVar][i].iP, f.slau[iVar][i].iT3, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iT3, -f.slau[iVar][i].at3);
					 }
					 if ((f.slau[iVar][i].iS3 > -1) && (fabs(f.slau[iVar][i].as3) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].as3, f.slau[iVar][i].iP, f.slau[iVar][i].iS3, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iS3, -f.slau[iVar][i].as3);
					 }
					 if ((f.slau[iVar][i].iW3 > -1) && (fabs(f.slau[iVar][i].aw3) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].aw3, f.slau[iVar][i].iP, f.slau[iVar][i].iW3, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iW3, -f.slau[iVar][i].aw3);
					 }
					 if ((f.slau[iVar][i].iB3 > -1) && (fabs(f.slau[iVar][i].ab3) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].ab3, f.slau[iVar][i].iP, f.slau[iVar][i].iB3, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iB3, -f.slau[iVar][i].ab3);
					 }


					 if ((f.slau[iVar][i].iE4 > -1) && (fabs(f.slau[iVar][i].ae4) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].ae4, f.slau[iVar][i].iP, f.slau[iVar][i].iE4, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iE4, -f.slau[iVar][i].ae4);
					 }
					 if ((f.slau[iVar][i].iN4 > -1) && (fabs(f.slau[iVar][i].an4) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].an4, f.slau[iVar][i].iP, f.slau[iVar][i].iN4, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iN4, -f.slau[iVar][i].an4);
					 }
					 if ((f.slau[iVar][i].iT4 > -1) && (fabs(f.slau[iVar][i].at4) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].at4, f.slau[iVar][i].iP, f.slau[iVar][i].iT4, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iT4, -f.slau[iVar][i].at4);
					 }
					 if ((f.slau[iVar][i].iS4 > -1) && (fabs(f.slau[iVar][i].as4) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].as4, f.slau[iVar][i].iP, f.slau[iVar][i].iS4, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iS4, -f.slau[iVar][i].as4);
					 }
					 if ((f.slau[iVar][i].iW4 > -1) && (fabs(f.slau[iVar][i].aw4) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].aw4, f.slau[iVar][i].iP, f.slau[iVar][i].iW4, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iW4, -f.slau[iVar][i].aw4);
					 }
					 if ((f.slau[iVar][i].iB4 > -1) && (fabs(f.slau[iVar][i].ab4) > nonzeroEPS)) {
						 addelmsimplesparse(sparseM, -f.slau[iVar][i].ab4, f.slau[iVar][i].iP, f.slau[iVar][i].iB4, true);
						 setValueIMatrix(&sparseS, f.slau[iVar][i].iP, f.slau[iVar][i].iB4, -f.slau[iVar][i].ab4);
					 }

				 }

			 }

			
	   } // for

	   if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar == NUSHA) || (iVar == PAM) ||
		   (iVar == TURBULENT_KINETIK_ENERGY) || (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		   (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)) {
#pragma omp parallel for
		   for (integer i_1 = 0; i_1 < f.maxelm + f.maxbound; i_1++) {
			   if (rthdsd[i_1] != rthdsd[i_1]) {
				   switch (iVar) {
				   case VELOCITY_X_COMPONENT: printf("VX problem\n");
					   break;
				   case VELOCITY_Y_COMPONENT: printf("VY problem\n");
					   break;
				   case VELOCITY_Z_COMPONENT: printf("VZ problem\n");
					   break;
				   case PAM: printf("PAM problem\n");
					   break;
				   case NUSHA:
					   printf("NUSHA problem\n");
					   break;
				   case TURBULENT_KINETIK_ENERGY:
					   printf("TURBULENT_KINETIK_ENERGY problem\n");
					   break;
				   case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
					   printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA problem\n");
					   break;
				   case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
					   printf("TURBULENT_KINETIK_ENERGY_Standart K-EPSILON problem\n");
					   break;
				   case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
					   printf("TURBULENT_DISSIPATION_RATE_EPSILON_Standart K-EPSILON problem\n");
					   break;
				   case TEMP: printf("TEMP problem\n");
					   break;
				   }
				   printf("pre bound ASSEMBLE CONTROL rthdsd part.\n");
				   printf("NAN or INF in premeshin.txt file. Power in control volume= %lld is undefined...\n", i_1);
				   printf("ispolzuite poslednuu versiu Mesh generator AliceMesh. 23.09.2018.\n");
				   system("PAUSE");
				   exit(1);
			   }
		   }
	   }

       // Для граничных узлов расчётной сетки:
#pragma omp parallel for
	   for (integer i=0; i<f.maxbound; i++) {

           // К граничным условиям Дирихле релаксации применять ненужно ? Это спорный вопрос.
           // Пока пробный запуск будет сделан в случае когда к граничным условиям Дирихле применяется нижняя релаксация.
		   integer iVar_in = iVar;
		   if (iVar == NUSHA) iVar_in = NUSHA_SL;
		   if (iVar == TURBULENT_KINETIK_ENERGY) iVar_in = TURBULENT_KINETIK_ENERGY_SL;
		   if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) iVar_in = TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL;
		   if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) iVar_in = TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL;
		   if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) iVar_in = TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL;

		   if (f.slau_bon[iVar_in][i].iI>-1) {
			   // Если условие Неймана то нижняя релаксация:

			   /*
		       switch (iVar) {
			   case VX: if (!bBiCGStabSaad) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw/f.alpha[iVar], f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true); }
						 if (!bBiCGStabSaad) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw/f.alpha[iVar]);  }
				         rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b+(1-f.alpha[iVar])*f.slau_bon[iVar][i].aw*f.potent[VXCOR][f.slau_bon[iVar][i].iW]/f.alpha[iVar];
				         break;
			   case VY: if (!bBiCGStabSaad) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw/f.alpha[iVar], f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true);}
						 if (!bBiCGStabSaad) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw/f.alpha[iVar]);}
                         rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b+(1-f.alpha[iVar])*f.slau_bon[iVar][i].aw*f.potent[VYCOR][f.slau_bon[iVar][i].iW]/f.alpha[iVar];
				         break;
			   case VZ: if (!bBiCGStabSaad) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw/f.alpha[iVar], f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true);}
						 if (!bBiCGStabSaad) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw/f.alpha[iVar]);}
                         rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b+(1-f.alpha[iVar])*f.slau_bon[iVar][i].aw*f.potent[VZCOR][f.slau_bon[iVar][i].iW]/f.alpha[iVar];
				         break;
			       case PAM: // PRESSURE: 
					   if (!bBiCGStabSaad) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true);}
					   if (!bBiCGStabSaad) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw); }
			             rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b;
				         break;
		       } // switch
			   */

			   //*
			   // Граничные значения не должны подвергаться нижней релаксации.
			   // Был проведён вычислительныйэксперимент и было выяснено, что невязки опускаются ниже в случае
			   // если в граничном условии отсутствует нижняя релаксация. Функция equation3DtoCRS исправлена согласованным образом !.
			   // Но вообще в принципе особо страшного не было и в первом варианте !.
			   switch (iVar) {
			   case VELOCITY_X_COMPONENT: if ((!bBiCGStabSaad)/* || (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true); }
						if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw); }
				         rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b;
				         break;
			   case VELOCITY_Y_COMPONENT: if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true); }
						if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw); }
                         rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b;
				         break;
			   case VELOCITY_Z_COMPONENT: if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true); }
						if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw); }
                         rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b;
				         break;
			   case NUSHA: // Модифицированная кинематическая турбулентная вязкость.
				        rthdsd[f.slau_bon[NUSHA_SL][i].iW] = f.slau_bon[NUSHA_SL][i].b;
				        break;
			   case TURBULENT_KINETIK_ENERGY: // Кинетическая энергия турбулентных пульсаций.
				   rthdsd[f.slau_bon[TURBULENT_KINETIK_ENERGY_SL][i].iW] = f.slau_bon[TURBULENT_KINETIK_ENERGY_SL][i].b;
				   break;
			   case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: // Удельная скорость диссипации кинетической энергии турбулентных пульсаций.
				   rthdsd[f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].iW] = f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].b;
				   break;
			   case TURBULENT_KINETIK_ENERGY_STD_K_EPS: // Кинетическая энергия турбулентных пульсаций..
				   rthdsd[f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL][i].iW] = f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL][i].b;
				   break;
			   case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: // Скорость диссипации кинетической энергии турбулентных пульсаций.
				   rthdsd[f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL][i].iW] = f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL][i].b;
				   break;
			   case PAM: // PRESSURE: 
					   if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true); }
					   if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw); }
			             rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b;
				         break;
		       } // switch
			   //*/
		   }
		   else {
			   // Если стоит условие Дирихле:

               switch (iVar) {
			   case VELOCITY_X_COMPONENT: case VELOCITY_Y_COMPONENT:case VELOCITY_Z_COMPONENT:
				   if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true); }
						if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw); }
				         rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b;
#if doubleintprecision == 1
						 //printf("%lld, eto ob=%lld \n",i, f.slau_bon[iVar][i].iW); getchar(); // debug
#else
						 //printf("%d, eto ob=%d \n",i, f.slau_bon[iVar][i].iW); getchar(); // debug
#endif
						 
                         /* // Граничные узлы не подвергаются нижней релаксации 
						 if (f.slau_bon[iVar][i].iW==13671) {
							 printf("in=%e, nominal=%e\n",f.potent[VXCOR][13671],f.potent[VXCOR][f.slau_bon[iVar][i].iW]);
	                         getchar();
						 }*/
				   break;			                              
			   case NUSHA: // Модифицированная кинематическая турбулентная вязкость.
				   rthdsd[f.slau_bon[NUSHA_SL][i].iW] = f.slau_bon[NUSHA_SL][i].b;
				   break;
			   case TURBULENT_KINETIK_ENERGY: // Кинетическая энергия турбулентных пульсаций.
				   rthdsd[f.slau_bon[TURBULENT_KINETIK_ENERGY_SL][i].iW] = f.slau_bon[TURBULENT_KINETIK_ENERGY_SL][i].b;
				   break;
			   case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA: // Удельная скорость диссипации кинетической энергии турбулентных пульсаций.
				   rthdsd[f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].iW] = f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL][i].b;
				   break;
			   case TURBULENT_KINETIK_ENERGY_STD_K_EPS: // Кинетическая энергия турбулентных пульсаций.
				   rthdsd[f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL][i].iW] = f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL][i].b;
				   break;
			   case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS: // Скорость диссипации кинетической энергии турбулентных пульсаций.
				   rthdsd[f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL][i].iW] = f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL][i].b;
				   break;
			     case PAM: // PRESSURE: 
					 if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { addelmsimplesparse(sparseM, f.slau_bon[iVar][i].aw, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, true); }
					 if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) { setValueIMatrix(&sparseS, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].iW, f.slau_bon[iVar][i].aw); }
			             rthdsd[f.slau_bon[iVar][i].iW]=f.slau_bon[iVar][i].b;
				         break;
		       } // switch

		   }

		   if ((!bBiCGStabSaad) /*|| (bBiCGStabSaad && (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2))*/) {
               const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

			   if ((f.slau_bon[iVar_in][i].iI>-1) && (fabs(f.slau_bon[iVar_in][i].ai) > nonzeroEPS)) {
				   addelmsimplesparse(sparseM, -f.slau_bon[iVar_in][i].ai, f.slau_bon[iVar_in][i].iW, f.slau_bon[iVar_in][i].iI, true);
                   setValueIMatrix(&sparseS, f.slau_bon[iVar_in][i].iW, f.slau_bon[iVar_in][i].iI, -f.slau_bon[iVar_in][i].ai);
				   // if (iVar==PAM) printf("neiman...\n"); // debug
			   }
		   }

	   } // for

	   
	   if ((iVar == VELOCITY_X_COMPONENT) || (iVar == VELOCITY_Y_COMPONENT) || (iVar == VELOCITY_Z_COMPONENT) || (iVar==NUSHA) || (iVar==PAM)||
		   (iVar== TURBULENT_KINETIK_ENERGY)||(iVar== TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) ||
		   (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) || (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS)) {
#pragma omp parallel for
		   for (integer i_1 = 0; i_1 < f.maxelm+f.maxbound; i_1++) {
			   if (rthdsd[i_1] != rthdsd[i_1]) {
				   switch (iVar) {
				   case VELOCITY_X_COMPONENT: printf("VX problem\n");
					   break;
				   case VELOCITY_Y_COMPONENT: printf("VY problem\n");
					   break;
				   case VELOCITY_Z_COMPONENT: printf("VZ problem\n");
					   break;
				   case PAM: printf("PAM problem\n");
					   break;
				   case NUSHA:
					   printf("NUSHA problem\n");
					   break;
				   case TURBULENT_KINETIK_ENERGY:
					   printf("TURBULENT_KINETIK_ENERGY problem\n");
					   break;
				   case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA:
					   printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA problem\n");
					   break;
				   case TURBULENT_KINETIK_ENERGY_STD_K_EPS:
					   printf("TURBULENT_KINETIK_ENERGY_Standart K-EPSILON problem\n");
					   break;
				   case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS:
					   printf("TURBULENT_DISSIPATION_RATE_EPSILON_Standart K-EPSILON problem\n");
					   break;
				   case TEMP: printf("TEMP problem\n");
					   break;
				   }
				   printf("POST ASSEMBLE CONTROL rthdsd part.\n");
				   printf("NAN or INF in premeshin.txt file. Power in control volume= %lld is undefined...\n", i_1);
				   printf("ispolzuite poslednuu versiu Mesh generator AliceMesh. 23.09.2018.\n");
				   system("PAUSE");
				   exit(1);
			   }
		   }
	   }

	}
	else { 
		    // Температура
		
		   if (!bBiCGStabSaad) {
        	   // выделение памяти и инициализация для 
	           // простейшей разреженной матрицы.
	           initsimplesparse(sparseM, t.maxelm + t.maxbound);
		       initIMatrix(&sparseS, t.maxelm + t.maxbound);
		   }

		  

            // Диагональный предобуславливатель:
			// после его применения все элементы нужно сравнивать с единицей. 
			// 1.0 характерное значение.

			/*
			// По видимому диагональный предобуславливатель 
			// влияент на LR вобщем в если его раскоментировать то будет ошибка
			// Пусть лучше он будет в закомментированном состоянии.
			// Внутренние узлы:
            for (i=0; i<t.maxelm; i++) {
                t.slau[i].ae/=t.slau[i].ap;
                t.slau[i].an/=t.slau[i].ap;
                t.slau[i].at/=t.slau[i].ap;
                t.slau[i].as/=t.slau[i].ap;
			    t.slau[i].aw/=t.slau[i].ap;
			    t.slau[i].ab/=t.slau[i].ap;
                t.slau[i].b/=t.slau[i].ap;
                t.slau[i].ap=1.0;
		    }

			// Граничные узлы:
			for (i=0; i<t.maxbound; i++) {
				t.slau_bon[i].ai/=t.slau_bon[i].aw;
				t.slau_bon[i].b/=t.slau_bon[i].aw;
				t.slau_bon[i].aw=1.0;
			}
			*/
#pragma omp parallel for
			for (integer i=0; i<t.maxelm+t.maxbound; i++) {
				rthdsd[i]=0.0; // обнуление
			}
			
			// 7 августа 2016 введена нижняя релаксация для температуры в матрицу СЛАУ.

		    // запись уравнений для внутренних узлов в матрицу.
#pragma omp parallel for
		    for (integer i=0; i<t.maxelm; i++) {
#if doubleintprecision == 1
				//printf("%lld %e %e %e %e %e %e %e %e\n",i, t.slau[i].ap,  t.slau[i].ae,  t.slau[i].aw,  t.slau[i].an,  t.slau[i].as,  t.slau[i].at,  t.slau[i].ab);
#else
				//printf("%d %e %e %e %e %e %e %e %e\n",i, t.slau[i].ap,  t.slau[i].ae,  t.slau[i].aw,  t.slau[i].an,  t.slau[i].as,  t.slau[i].at,  t.slau[i].ab);
#endif
				//getchar(); // debug

				if (!bBiCGStabSaad) {
                   addelmsimplesparse(sparseM, t.slau[i].ap, t.slau[i].iP, t.slau[i].iP, true);
				  // addelmsimplesparse(sparseM, t.slau[i].ap/t.alpha, t.slau[i].iP, t.slau[i].iP, true);
				}

				 /*if (!(t.slau[i].ap==t.slau[i].ap)) {
					 printf("error matrix construct node=%i\n",i);
					 getchar();
				 }
				 */

				if (!bBiCGStabSaad) {
					setValueIMatrix(&sparseS,t.slau[i].iP,t.slau[i].iP,t.slau[i].ap);
					//setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iP, t.slau[i].ap/t.alpha);
				}
		       //rthdsd[t.slau[i].iP]=t.slau[i].b;
				rthdsd[t.slau[i].iP] = t.slau[i].b + (1 - t.alpha)*t.slau[i].ap*t.potent[t.slau[i].iP] / t.alpha;
				//rthdsd[t.slau[i].iP] = t.slau[i].b + rthdsd[t.slau[i].iP] + ((1 - t.alpha)*t.slau[i].ap*t.potent[t.slau[i].iP] )/ t.alpha;
				
				//told_iter
            
				if (!bBiCGStabSaad) {
			         const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

			
			         if ((t.slau[i].iE>-1) && (fabs(t.slau[i].ae) > nonzeroEPS)){
				         addelmsimplesparse(sparseM, -t.slau[i].ae, t.slau[i].iP, t.slau[i].iE, true);
                         setValueIMatrix(&sparseS,t.slau[i].iP,t.slau[i].iE,-t.slau[i].ae);
			         }
			         if ((t.slau[i].iN>-1) && (fabs(t.slau[i].an) > nonzeroEPS)) {
				        addelmsimplesparse(sparseM, -t.slau[i].an, t.slau[i].iP, t.slau[i].iN, true);
				        setValueIMatrix(&sparseS,t.slau[i].iP,t.slau[i].iN,-t.slau[i].an);
			         }
			         if ((t.slau[i].iT>-1) && (fabs(t.slau[i].at) > nonzeroEPS)) {
				        addelmsimplesparse(sparseM, -t.slau[i].at, t.slau[i].iP, t.slau[i].iT, true);
                        setValueIMatrix(&sparseS,t.slau[i].iP,t.slau[i].iT,-t.slau[i].at);
			         }
			         if ((t.slau[i].iS>-1) && (fabs(t.slau[i].as) > nonzeroEPS)) {
				        addelmsimplesparse(sparseM, -t.slau[i].as, t.slau[i].iP, t.slau[i].iS, true);
                       setValueIMatrix(&sparseS,t.slau[i].iP,t.slau[i].iS,-t.slau[i].as);
			         }
			         if ((t.slau[i].iW>-1) && (fabs(t.slau[i].aw) > nonzeroEPS)) {
				        addelmsimplesparse(sparseM, -t.slau[i].aw, t.slau[i].iP, t.slau[i].iW, true);
                        setValueIMatrix(&sparseS,t.slau[i].iP,t.slau[i].iW,-t.slau[i].aw);
			         }
			         if ((t.slau[i].iB>-1) && (fabs(t.slau[i].ab) > nonzeroEPS)) {
				        addelmsimplesparse(sparseM, -t.slau[i].ab, t.slau[i].iP, t.slau[i].iB, true);
                        setValueIMatrix(&sparseS,t.slau[i].iP,t.slau[i].iB,-t.slau[i].ab);
		             }

					 if (b_on_adaptive_local_refinement_mesh) {
						 if ((t.slau[i].iE2>-1) && (fabs(t.slau[i].ae2) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].ae2, t.slau[i].iP, t.slau[i].iE2, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iE2, -t.slau[i].ae2);
						 }
						 if ((t.slau[i].iN2>-1) && (fabs(t.slau[i].an2) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].an2, t.slau[i].iP, t.slau[i].iN2, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iN2, -t.slau[i].an2);
						 }
						 if ((t.slau[i].iT2>-1) && (fabs(t.slau[i].at2) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].at2, t.slau[i].iP, t.slau[i].iT2, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iT2, -t.slau[i].at2);
						 }
						 if ((t.slau[i].iS2>-1) && (fabs(t.slau[i].as2) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].as2, t.slau[i].iP, t.slau[i].iS2, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iS2, -t.slau[i].as2);
						 }
						 if ((t.slau[i].iW2>-1) && (fabs(t.slau[i].aw2) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].aw2, t.slau[i].iP, t.slau[i].iW2, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iW2, -t.slau[i].aw2);
						 }
						 if ((t.slau[i].iB2>-1) && (fabs(t.slau[i].ab2) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].ab2, t.slau[i].iP, t.slau[i].iB2, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iB2, -t.slau[i].ab2);
						 }

						 if ((t.slau[i].iE3>-1) && (fabs(t.slau[i].ae3) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].ae3, t.slau[i].iP, t.slau[i].iE3, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iE3, -t.slau[i].ae3);
						 }
						 if ((t.slau[i].iN3>-1) && (fabs(t.slau[i].an3) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].an3, t.slau[i].iP, t.slau[i].iN3, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iN3, -t.slau[i].an3);
						 }
						 if ((t.slau[i].iT3>-1) && (fabs(t.slau[i].at3) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].at3, t.slau[i].iP, t.slau[i].iT3, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iT3, -t.slau[i].at3);
						 }
						 if ((t.slau[i].iS3>-1) && (fabs(t.slau[i].as3) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].as3, t.slau[i].iP, t.slau[i].iS3, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iS3, -t.slau[i].as3);
						 }
						 if ((t.slau[i].iW3>-1) && (fabs(t.slau[i].aw3) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].aw3, t.slau[i].iP, t.slau[i].iW3, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iW3, -t.slau[i].aw3);
						 }
						 if ((t.slau[i].iB3>-1) && (fabs(t.slau[i].ab3) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].ab3, t.slau[i].iP, t.slau[i].iB3, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iB3, -t.slau[i].ab3);
						 }

						 if ((t.slau[i].iE4>-1) && (fabs(t.slau[i].ae4) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].ae4, t.slau[i].iP, t.slau[i].iE4, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iE4, -t.slau[i].ae4);
						 }
						 if ((t.slau[i].iN4>-1) && (fabs(t.slau[i].an4) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].an4, t.slau[i].iP, t.slau[i].iN4, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iN4, -t.slau[i].an4);
						 }
						 if ((t.slau[i].iT4>-1) && (fabs(t.slau[i].at4) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].at4, t.slau[i].iP, t.slau[i].iT4, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iT4, -t.slau[i].at4);
						 }
						 if ((t.slau[i].iS4>-1) && (fabs(t.slau[i].as4) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].as4, t.slau[i].iP, t.slau[i].iS4, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iS4, -t.slau[i].as4);
						 }
						 if ((t.slau[i].iW4>-1) && (fabs(t.slau[i].aw4) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].aw4, t.slau[i].iP, t.slau[i].iW4, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iW4, -t.slau[i].aw4);
						 }
						 if ((t.slau[i].iB4>-1) && (fabs(t.slau[i].ab4) > nonzeroEPS)) {
							 addelmsimplesparse(sparseM, -t.slau[i].ab4, t.slau[i].iP, t.slau[i].iB4, true);
							 setValueIMatrix(&sparseS, t.slau[i].iP, t.slau[i].iB4, -t.slau[i].ab4);
						 }

					 }
			
                 // Построчное заполнение матрицы для внутренних КО.
				}
		    }

			

			if (0) {
				// debug
			    for (integer i=0; i<t.maxelm; i++) {
#if doubleintprecision == 1
					printf("%lld %e\n", i, rthdsd[i]);
#else
					printf("%d %e\n", i, rthdsd[i]);
#endif
			         
					 if (i % 10 == 0) {
						// getchar();
						 system("pause");
					 }
			    }
			}

            // Запись уравнений для граничных узлов в матрицу:
#pragma omp parallel for
            for (integer i=0; i<t.maxbound; i++) {
                
				if (!bBiCGStabSaad) {
					
				    addelmsimplesparse(sparseM, t.slau_bon[i].aw, t.slau_bon[i].iW, t.slau_bon[i].iW, true);
				    setValueIMatrix(&sparseS,t.slau_bon[i].iW ,t.slau_bon[i].iW, t.slau_bon[i].aw);
					//addelmsimplesparse(sparseM, t.slau_bon[i].aw/t.alpha, t.slau_bon[i].iW, t.slau_bon[i].iW, true);
					//setValueIMatrix(&sparseS, t.slau_bon[i].iW, t.slau_bon[i].iW, t.slau_bon[i].aw/t.alpha);
				}
		         rthdsd[t.slau_bon[i].iW]=t.slau_bon[i].b;				
			    //rthdsd[t.slau_bon[i].iW] = t.slau_bon[i].b + (1 - t.alpha)*t.slau_bon[i].aw*t.potent[t.slau_bon[i].iW] / t.alpha;

				 if (!bBiCGStabSaad) {
				     const doublereal nonzeroEPS=1e-37; // для отделения вещественного нуля

				     if ((t.slau_bon[i].iI>-1) && (fabs(t.slau_bon[i].ai) > nonzeroEPS)) {
				        addelmsimplesparse(sparseM, -t.slau_bon[i].ai, t.slau_bon[i].iW, t.slau_bon[i].iI, true);
                        setValueIMatrix(&sparseS,t.slau_bon[i].iW,t.slau_bon[i].iI,-t.slau_bon[i].ai);
		             }
				 }
            }

			if (rthdsd_no_radiosity_patch == nullptr) {
				rthdsd_no_radiosity_patch = new doublereal[t.maxelm+t.maxbound + 1];
				if (rthdsd_no_radiosity_patch == nullptr) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment for rthdsd_no_radiosity_patch my_solverv0_03 my_agregat_amg.cpp...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
			}

			if (rthdsd_no_radiosity_patch != nullptr) {
#pragma omp parallel for
				for (integer i23 = 0; i23 < t.maxelm + t.maxbound; i23++) {
					rthdsd_no_radiosity_patch[i23] = rthdsd[i23];
				}
			}

			radiosity_patch_for_vacuum_Prism_Object_(rthdsd, b, lb, t.maxelm);

            //printM_and_CSIR(sparseM, t.maxelm + t.maxbound);
			

			//printf("matrix printing...\n");
			//getchar();

		    /*
		    printf("add elem in matrix ready\n");
		    // освобождение памяти из под t.slau
		    if (t.free_temper_level2) {
			    printf("max memory control point1...\n");
			    getchar();
			    free_level2_temp(t);
			    printf("t.slau free memory 15N...\n");
			    getchar();
		    }
		    */

			// Учёт модели излучения в Prism Object:
			// Она модифицирует источниковый член в СЛАУ.
			//radiosity_patch_for_vacuum_Prism_Object(t.slau, t.slau_bon, b, lb, t.maxelm);
			
			
		 
	}

	


	/*if (iVar==PAm) {
		printM_and_CSIR(sparseM, maxelm); 
	    getchar();
	    printf("\n");
	}// */ // debug 


    // Блок решения СЛАУ:
	// Следует отметить статью П.Н. Вабищевича который применял в ней ICCG солвер для поправки давления,
	// в ряде других статей также отдаётся предпочтение методам Крыловского типа.

    //1. Прямые методы решения СЛАУ:

    // 1.1. Для плотных матриц.

	// решает СЛАУ методом разложения Холесского
	// для симметричной положительно определённой матрицы s.
	// без учёта разреженности матрицы. В два раза более 
	// эффективна чем метод прямого исключения Гаусса.
	//eqsolv_simple_holesskii(s, nodes, rthdsd, potent);
	// решает СЛАУ методом исключения Гаусса
	// без выбора главного элемента и без 
	// учёта разреженности матрицы.
	//eqsolve_simple_gauss(s, nodes, rthdsd, potent[Temp]);

    // 1.2. Для разреженных матриц.

	// Метод К.Ф. Гаусса для разреженной возможно несимметричной матрицы
    //calculateSPARSEgaussArray(&sparseS, potent, rthdsd);


	//2. Итерационные методы решения СЛАУ:

	// 2.1. Для плотных матриц.

	// Метод Гаусса-Зейделя-Ричардсона-Либмана SOR
	// для неразреженной матрицы СЛАУ.
	// метод характеризуется медленной сходимостью.
	//Seidel(s, rthdsd, potent, nodes, 1e-5, 1.855);
	// Решает неразреженную СЛАУ методом
	// сопряжённых градиентов.
	// для симметричной положительно определённой матрицы s.
	// Метод является уточняющим, т.е. если известен вектор
	// начального приближения то его нужно поставить вместо nullptr.
	//potent=SoprGrad(s, rthdsd, nullptr, nodes);
		
    // 2.2. Для разреженных матриц.

	// 2.2.1 МСГ для CRS (CSIR, CSIR_ITL) формата симметричный случай.
    ///*
    //simplesparsetoCRS(sparseM, val, col_ind, row_ptr, nodes);
	//printM_and_CRS(sparseM,val,col_ind,row_ptr,nodes);
    //potent=SoprGradCRS(val, col_ind, row_ptr, rthdsd, nullptr, nodes);//*/

	

	// МСГ (Congruate Gradients) для CSIR формата хранения
	// симметричный положительно определённый случай.
    //simplesparsetoCSIR(sparseM, adiag, altr, jptr, iptr, nodes);
	//printM_and_CSIR(sparseM, adiag, altr, jptr, iptr,  nodes);
	//integer inz=(int64_t)((sparseM.n-nodes)/2.0);
	//potent=SoprGradCSIR(adiag, altr, jptr, iptr, rthdsd, nullptr, nodes, inz);
	//potent=SoloveichikAlgCSIR_SPD(nodes, adiag, altr, jptr, iptr, rthdsd, nullptr, true);
	//potent=SoloveichikAlgCSIR_SPDgood(nodes, inz, adiag, altr, jptr, iptr, rthdsd, nullptr, true);
	if (iVar==PAM) {

		
         //for (integer i=0; i<230000; i++) { // debug
		// Убрал инициализацию нулём 9 августа 2016 года.
		// Обязательно нужна инициализация нулем и никакой нижней релаксации быть недолжно.
		 for (integer i=0; i<(f.maxelm+f.maxbound); i++) f.potent[PAM][i]=0.0;
		 integer maxiter=2600; //2000 120
		 doublereal *val=nullptr;
         integer *col_ind=nullptr, *row_ptr=nullptr;

		 if (!bBiCGStabSaad) {
			
		    freeIMatrix(&sparseS); // освобождение оперативной памяти.
		    simplesparsetoCRS(sparseM, val, col_ind, row_ptr, (f.maxelm+f.maxbound)); // преобразование матрицы из одного формата хранения в другой.
			simplesparsefree(sparseM, (f.maxelm + f.maxbound));
		 }
		 else {
			 if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2) {
				 /*Lr1sK_up
				 freeIMatrix(&sparseS); // освобождение оперативной памяти.
				 simplesparsetoCRS(sparseM, val, col_ind, row_ptr, (f.maxelm + f.maxbound)); // преобразование матрицы из одного формата хранения в другой.
				 simplesparsefree(sparseM, (f.maxelm + f.maxbound));
				 if (val == nullptr) {
					 // отладка.
					 printf("nullptr pointer val\n");
					 getchar();
				 }
				 */
			 }
		 }

		 
		
		 // Пользователю желательно предоставить несколько солверов на выбор:
		 // Рекомендации по выбору того или иного решающего устройства можно найти в документации.
		 // Если коротко, то никогда не используйте SOR3DS для задач большой размерности. 
		 // Для задач большой размерности обычно подходит MICCGS - сопряжённые градиенты с предобуславливанием 
		 // неполным разложением Холецкого. Правда, наверное (это требует проверки), для задач ещё большей размерности
		 // видимо лучше LR1sk решатель. LR1sk - сходится за меньшее число итераций чем MICCG однако каждая его итерация более ресурсоёмка.
		 // В настоящее время предпочтение отдаётся LR1sk Solver. LR1sk - солвер даёт расходимость на реальных задачах.
		 const integer SOR3DS=0; // SOR3D Solver - метод последовательной верхней релаксации в 3D пространстве. // простой нерабочий для задач большой размерности.
		 const integer MICCGS=1; // Modify Incomplete Cholesky Congruate Gradient SOLVER // рекомендуемый старый
		 const integer LR1SKS=2; // LR1sk - А.А. Фомин, Л.Н. Фомина // рекомедуемый новый
		 const integer LRN=3; // Данный солвер позволяет на всей границе гидродинамической подобласти ставить условие Неймана, что очень важно для сходимости. К тому же данный солвер рекомендован С.Патанкаром.
		 const integer BICGSTAB=4;
		 const integer BICGSTABP=5; // предобусловленный алгоритм Ван дер Ворста. Предобуславливатели на основе SPARSKIT2.
		 integer iSOLVER=MICCGS;//LR1SKS; //MICCGS; // выбор решающего алгоритма.

		 if (bBiCGStabSaad) {
			 iSOLVER=BICGSTABP;
		 }

		 // В методичке по CFX рекомендуют уровень среднеквадратичных невязок равный 1.0e-4.
		 // На уравнении теплопередачи в твёрдом теле было установлено что уровень среднеквадратичных невязок должен быть гдето 1.e-5 
		 dterminatedTResudual=3.0e-6; //f.resICCG; // O(h!2) 1e-40

		 //printf("Ok\n");
		 //getchar();
		 
		 //SOR3D(f.slau[PAM], f.slau_bon[PAM], f.potent[PAM], f.maxelm, f.maxbound, PAM);
		 if (f.bLR1free)
		 {

			 if (0) {
		         // Условие однозначности состоит в том чтобы
		         // среднее значение поправки давления было равно нулю.

		         //printf("regularization condition for Pamendment...\n");
		         //printf("please, press any key to continue...\n");
		         //getchar();

		         // Следуя П.Н. Вабищевичу, А.Н.Павлову, А.Г.Чурбанову
		         // будем использовать условие однозначности для 
		         // определения поля поправки давления в том случае
		         // если на всей границе расчётной области стоит 
		         // однородное условие Неймана.
	             // Условие однозначности:
	             doublereal pamscal1=0.0;
			     doublereal Vol=0.0;
	         
			     // условие однозначной разрешимости.
		         for (integer iP1=0; iP1<(f.maxelm); iP1++) {
		             // вычисление размеров текущего контрольного объёма:
	                 doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	                 volume3D(iP1, f.nvtx, f.pa, dx, dy, dz);
				     Vol+=dx*dy*dz;

			         pamscal1+=rthdsd[iP1]*dx*dy*dz;     
		         }
	
			     pamscal1=(doublereal)(pamscal1/(f.maxelm*Vol));
			     for (integer iP1=0; iP1<(f.maxelm); iP1++) {
				     rthdsd[iP1]-=pamscal1;
			     }

			     // Когда на всей границе стоит условие Неймана
			     // то давление необходимо фиксировать в одной точке
			     // для того чтобы сходились такие методы как ICCG(0) и LR1sk.
			 }
			
			 // метод Гаусса Зейделя быстрее сходится когда давление не фиксировано ни в одной точке.
            //SOR3D(f.slau[PAM], f.slau_bon[PAM], f.potent[PAM], f.maxelm, f.maxbound, PAM);

			 //printf("bLR1free Ok\n");
		     //getchar();

			 // плавающий полилинейный солвер показывает более быструю сходимость.
			 //solveLRn(f.potent[PAM], rthdsd, (f.maxelm+f.maxbound), PAM, 20, true, true, f.neighbors_for_the_internal_node, f.maxelm, f.slau, f.slau_bon, f.iN, f.id, f.iWE, f.iSN, f.iBT, f.alpha, f.maxbound);
			 // максимальное количество итераций полилинейного метода оценивается равным
			 // 4*Nequations, где Nequations - число линейных уравнений. 

			 // ICCG solver работает и при условиях Неймана на всей границе расчётной области ???
			 
			 bool bexporttecplot=false;
			 bool bprintPAM=false;
			 bool bnorelaxPAM=true;
			 //iSOLVER=SOR3DS; // debug TODO

			 switch (iSOLVER) {
			 case MICCGS: // метод сопряжённых градиентов с предобуславливанием неполным разложением Холесского.
				           delete val; delete col_ind; delete row_ptr;
				           //printf("Ok!\n"); // debug
				           //getchar();
			               ICCG(iVar, sparseM, rthdsd, f.potent[PAM], f.maxelm + f.maxbound ,bprintmessage,false,maxiter); //->//
		                   // при использовании ICCG память из под sparseM освобождается автоматом
				          break;
			 case LR1SKS: // Алгоритм Фомин и Фомина (BiCGStab усиленный полилинейным методом). 
				
				           simplesparsefree(sparseM,(f.maxelm+f.maxbound));
		                   //->//->//solveLRn(f, f.potent[PAM], rthdsd, (f.maxelm+f.maxbound), PAM, maxiter, true, true); // полилинейный солвер для поправки давления

						   // А.А. Фомин, Л.Н. Фомина
		                   // Ускорение полилинейного рекуррентного метода в подпространствах крылова.
                           // Вестник томского государственного университета. Математика и механика №2(14) 2011год.
		                   LR1sK(f, f.slau[PAM], f.slau_bon[PAM], val, col_ind, row_ptr, f.maxelm, f.maxbound, PAM, rthdsd, f.potent[PAM], maxiter, bprintPAM, bexporttecplot);
			               delete val; delete col_ind; delete row_ptr;

				           if (bexporttecplot) {
					           // Экспорт в программу tecplot360 
					           // в случае обнаружения расходимости.
							   const int ianimate = 0;
					           exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior, ianimate,false,0,b,lb);
				           }
				          break;
			 case LRN: // полилинейный метод рекомендованный Патанкаром.
				        simplesparsefree(sparseM, (f.maxelm+f.maxbound));
						delete val; delete col_ind; delete row_ptr;
				        maxiter=10; // В CFX говорят что 10 прогонов должно хватать.
						
						solveLRn(f.potent[PAM], rthdsd, (f.maxelm+f.maxbound), PAM, maxiter,
							     bprintPAM, bnorelaxPAM, f.neighbors_for_the_internal_node, f.maxelm, f.slau, f.slau_bon,
								 f.iN, f.id, f.iWE, f.iSN, f.iBT, f.alpha, f.maxbound); // полилинейный солвер для поправки давления				        
				       break;
		     case BICGSTAB: // Би сопряжённые градиенты со стабилизацией. Просто алгоритм без какого-либо предобуславливания.
			                // В основе лежит CRS формат хранения.
							simplesparsefree(sparseM,(f.maxelm+f.maxbound));
				            Bi_CGStabCRS((f.maxelm+f.maxbound), val, col_ind, row_ptr, rthdsd, f.potent[PAM], maxiter);//->//
							
							delete val; delete col_ind; delete row_ptr;
				          break;
			 case SOR3DS: // Данный решатель используется как отладочный для проверки корректности сборки матрицы.
						   if (!bBiCGStabSaad) {
				               simplesparsefree(sparseM,(f.maxelm+f.maxbound));
						       delete val; delete col_ind; delete row_ptr;
						   }
				           SOR3D(f.slau[PAM], f.slau_bon[PAM], f.potent[PAM], nullptr, f.maxelm, f.maxbound, PAM, 1.0);
						 break;

			 case BICGSTABP: // Предобусловленый алгоритм Ван Дер Ворста.
				          if (!bBiCGStabSaad) {
				              delete val; delete col_ind; delete row_ptr;
						      simplesparsefree(sparseM,(f.maxelm+f.maxbound));
						  }
						 Bi_CGStab(&sparseS, f.slau[PAM], f.slau_bon[PAM], f.maxelm, f.maxbound, rthdsd, f.potent[PAM], maxiter, 1.0,PAM, m, f.bLR1free,b,lb,f.ifrontregulationgl,f.ibackregulationgl, dgx, dgy, dgz,s,ls, inumiter,color,dist_max,w,lw);
				 break;
			 }
			
		 }
		 else {

			 bool bexporttecplot=false;
			 bool bprintPAM=false;
			 bool bnorelaxPAM=true;

			 switch (iSOLVER) {
			 case MICCGS: // метод сопряжённых градиентов с предобуславливанием неполным разложением Холесского.
				           delete val; delete col_ind; delete row_ptr;
				           //printf("Ok!\n"); // debug
				           //getchar();
			               ICCG(iVar, sparseM, rthdsd, f.potent[PAM], f.maxelm + f.maxbound ,bprintmessage,false,maxiter); //->//
		                   // при использовании ICCG память из под sparseM освобождается автоматом
				          break;
			 case LR1SKS: // Алгоритм Фомин и Фомина (BiCGStab усиленный полилинейным методом). 
				
				           simplesparsefree(sparseM,(f.maxelm+f.maxbound));
		                   //->//->//solveLRn(f, f.potent[PAM], rthdsd, (f.maxelm+f.maxbound), PAM, maxiter, true, true); // полилинейный солвер для поправки давления
						  
		                   // А.А. Фомин, Л.Н. Фомина
		                   // Ускорение полилинейного рекуррентного метода в подпространствах крылова.
                           // Вестник томского государственного университета. Математика и механика №2(14) 2011год.
		                   LR1sK(f, f.slau[PAM], f.slau_bon[PAM], val, col_ind, row_ptr, f.maxelm, f.maxbound, PAM, rthdsd, f.potent[PAM], maxiter, bprintmessage, bexporttecplot);
			               delete val; delete col_ind; delete row_ptr;

				           if (bexporttecplot) {
					           // Экспорт в программу tecplot360 
					           // в случае обнаружения расходимости.
							   const int ianimate = 0;
					           exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior, ianimate,false,0, b, lb);
				           }
				          break;
			 case LRN: // полилинейный метод рекомендованный Патанкаром.
				        simplesparsefree(sparseM, (f.maxelm+f.maxbound));
						delete val; delete col_ind; delete row_ptr;
				        maxiter=20; // В CFX говорят что 10 прогонов должно хватать.
						
						solveLRn(f.potent[PAM], rthdsd, (f.maxelm+f.maxbound), PAM, maxiter,
							     bprintPAM, bnorelaxPAM, f.neighbors_for_the_internal_node, f.maxelm, f.slau, f.slau_bon,
								 f.iN, f.id, f.iWE, f.iSN, f.iBT, f.alpha, f.maxbound); // полилинейный солвер для поправки давления				        
				       break;

			 case BICGSTABP: 
				       if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2) {
					       // LR1SK
					       iVar = PAM;
						   Lr1sk_up(f, t,  f.slau[PAM], f.slau_bon[PAM], f.maxelm, f.maxbound, rthdsd, f.potent[PAM], maxiter, 1.0, iVar, false);//->//

						   
					       //LR1sK(f, f.slau[PAM], f.slau_bon[PAM], val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, f.potent[PAM], maxiter, bprintmessage, bexporttecplot);//->//
					       //freeIMatrix(&sparseS);
					       if (val != nullptr) {
						       delete val;
						       val = nullptr;
					       }
					       if (col_ind != nullptr) {
						       delete col_ind;
						       col_ind = nullptr;
					        }
					       if (row_ptr != nullptr) {
						      delete row_ptr;
						      row_ptr = nullptr;
					       }
						   
				        }
				        else {
					        // Предобусловленый алгоритм Ван Дер Ворста.
					        if (!bBiCGStabSaad) {
						       delete val; delete col_ind; delete row_ptr;
						       simplesparsefree(sparseM, (f.maxelm + f.maxbound));
					        }
					        Bi_CGStab(&sparseS, f.slau[PAM], f.slau_bon[PAM], f.maxelm, f.maxbound, rthdsd, f.potent[PAM], maxiter, 1.0, PAM, m, f.bLR1free, b,lb, f.ifrontregulationgl, f.ibackregulationgl, dgx, dgy, dgz,s,ls, inumiter, color, dist_max, w, lw);
				        }
						break;
			 
			 }

			 
		 }
		   

		 // Нижняя релаксация:

         
		 //for (i=0; i<(f.maxelm+f.maxbound); i++) f.potent[PAM][i]=0.0;
			 //SOR3D(f.slau[PAM], f.slau_bon[PAM], f.potent[PAM], nullptr, f.maxelm, f.maxbound, PAM, 1.0);
			//BTrules(f.slau[PAM], f.slau_bon[PAM], f.potent[PAM], nullptr, f.maxelm, f.maxbound, PAM, 1.0); // BT процесс с управлением.

		 // Пост сглаживание. Нужно проверить влияет оно или нет.
		 // метод Гаусса Зейделя быстрее сходится когда давление не фиксировано ни в одной точке.
         //SOR3D(f.slau[PAM], f.slau_bon[PAM], f.potent[PAM], f.maxelm, f.maxbound, PAM);

		 // getchar(); // debug pressure solver...

         // освобождение памяти
	    // delete rthdsd;
		//}
	}
	//*/

	// 2.2.2 Методы для CRS (CSIR, CSIR_ITL) формата возможно несимметричный случай.
	if (iVar!=PAM) {
		// разреженная матрица в формате CRS
        doublereal *val=nullptr;
        integer *col_ind=nullptr, *row_ptr=nullptr;

		if (iVar != TEMP) {

			

		    if (!bBiCGStabSaad) {
			    simplesparsetoCRS(sparseM, val, col_ind, row_ptr, (f.maxelm+f.maxbound));
                simplesparsefree(sparseM,(f.maxelm+f.maxbound));
			}
			else {
				if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2) {
					/*LR1sK_up
					freeIMatrix(&sparseS);
					simplesparsetoCRS(sparseM, val, col_ind, row_ptr, (f.maxelm + f.maxbound));
					simplesparsefree(sparseM, (f.maxelm + f.maxbound));
					if (val == nullptr) {
						printf("Null pointer val\n");
						getchar();
					}
					*/
				}
			}

			

			// Решение производится итерационным методом.
			// для метода LR1SK который рекомендован к использованию максимальное количество итераций
			// ограничено совсем небольшим значением. Т.к. данный метод обладает хорошей сходимостью, но
			// тратит на каждую итерацию значительно большее количество арифметических операций чем более
			// простые методы, например, алгоритм Соловейчика без предобуславливания.
	        integer maxiter=320; //120
			// limititer - максимальное количество итераций для простого метода (без предобуславливателя). простой
			// метод используется как правило без предобуславливания и на каждую итерацию расходует мало арефмитических
			// операций поэтому обладает низкой скоростью сходимости и для достижения данной сходимости требует больше итераций.
			const integer limititer=1000; // 120 в определённых итерационных методах используется увеличенное количество итераций.

	        
			// Перечень реализованных решателей СЛАУ:
			const integer SOLOVEICHIK=0; // алгоритм Соловейчика.
			const integer SOR3DALG=1; // метод релаксации основанный на методе Гаусса-Зейделя.
			const integer SOLVELRN=2; // полилинейный метод.
			const integer BICG=3; // бисопряжённые градиенты.
			const integer BICGSTAB=4; // алгоритм Ван-Дер-Ворста
			const integer LR1SK=5; // Фомин-Фомина данный алгоритм рекомендован к использованию.
			const integer BICGSTABP=6;// Предобусловленный алгоритм Ван Дер Ворста. Предобуславливатели из библиотеки SPARSKIT2.

			integer iSOLVER=LR1SK; // LR1SK выбор решающего устройства.

			if (bBiCGStabSaad) {
				iSOLVER=BICGSTABP;
			}

			// выбор предобуславливателя:
			const integer NOPRECOND=0; // без предобуславливания.
			const integer ILU0ITL=1; // неполное LU разложение уровня 0 реализованное в библиотеке ITL.
			const integer ILU0SAAD=2; // неполное LU разложение уровня 0 описанное в книге Й.Саада.

			integer iPRECOND=ILU0SAAD; // выбор предобуславливателя.

			 /*switch (iSOLVER) {
		        case BICGCRS: potent[iVar]=BiSoprGradCRS(val, col_ind, row_ptr, rthdsd, potent[iVar], maxelm, maxiter); break;
		        case SOLOVALGCRS: SoloveichikAlgCRS(maxelm, val, col_ind, row_ptr, rthdsd, potent[iVar], true, maxiter); break;
		        case BICGSTABCRS: potent[iVar]=Bi_CGStab(maxelm, val, col_ind, row_ptr, rthdsd, potent[iVar], maxiter); break;
		      default: potent[iVar]=Bi_CGStab(maxelm, val, col_ind, row_ptr, rthdsd, potent[iVar], maxiter); break;
	        }*/

			//dterminatedTResudual=1e-40*f.resLR1sk; // O(h!3) // 1e-16 обязательный эмпирический запас (для него надо найти какую-нибудь теоретическую оценку).
			// В одной методичке по CFX сказано, что уровня среднеквадратической невязки 1.e-4 должно хватать. Нами было установлено экспериментально
			// на примере уравнения теплопередачи что достаточен уровень 1.0e-5.
			dterminatedTResudual=1.0e-5;

			// в случае обнаружения проблем со сходимостью будет произведён экспорт
			// в программу tecplot360.
		    bool bexporttecplot=false;
			integer iVarCor=VXCOR; // идентификатор скорректированной скорости к которой будет осуществляться нижняя релаксация.

			

			switch (iSOLVER) {
			case SOLOVEICHIK: // Алгоритм Ю.Г.Соловейчика 1993 года.
				         switch(iPRECOND) {
						 case NOPRECOND: // Алгоритм Ю.Г.Соловейчика 1993 года. без предобуславливателя.
							           maxiter=limititer;
                                       SoloveichikAlgCRS((f.maxelm+f.maxbound), val, col_ind, row_ptr, rthdsd, f.potent[iVar], true, maxiter); // без предобуславливателя.  
									    freeIMatrix(&sparseS);
							           break;
						 case ILU0ITL: // Алгоритм Ю.Г.Соловейчика 1993 года. Встроен ILU0 предобуславливатель.
							           // Данный алгоритм сам освобождает память из под sparseS.
							           maxiter=limititer;
									   if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
										   SoloveichikAlg(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
											   f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], true, false, maxiter,
											   f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], b, lb, s, ls);
									   }
									   else if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
										   SoloveichikAlg(&sparseS, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS],
											   f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS], true, false,
											   maxiter, f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], b, lb, s, ls);
									   }
									   else
									   if (iVar == TURBULENT_KINETIK_ENERGY) {
										   SoloveichikAlg(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_SL], 
											   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[TURBULENT_KINETIK_ENERGY], true, false, maxiter,
											   f.alpha[TURBULENT_KINETIK_ENERGY_SL], b, lb, s, ls);
									   }
									   else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
										   SoloveichikAlg(&sparseS, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
											   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], true, false, maxiter,
											   f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], b, lb, s, ls);
									   }
									   else if (iVar == NUSHA) {
										   SoloveichikAlg(&sparseS, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[NUSHA], true, false, maxiter, f.alpha[NUSHA_SL], b, lb, s, ls);
									   }
									   else {
										   SoloveichikAlg(&sparseS, f.slau[iVar], f.slau_bon[iVar], f.maxelm, f.maxbound,
											   rthdsd, f.potent[iVar], true, false, maxiter, f.alpha[iVar], b, lb, s, ls);
									   }
									   break;
						 case ILU0SAAD: // Алгоритм Ю.Г.Соловейчика 1993 года. Встроен ILU0 предобуславливатель. Данный алгоритм сам освобождает память из под sparseS.
							           maxiter=limititer;
									   if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
										   SoloveichikAlg(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
											   f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], true, true, maxiter,
											   f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], b, lb, s, ls);
									   }
									   else  if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
										   SoloveichikAlg(&sparseS, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
											   f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS], true, true, maxiter,
											   f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], b, lb, s, ls);
									   }
									   else  if (iVar == TURBULENT_KINETIK_ENERGY) {
										   SoloveichikAlg(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_SL], 
											   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[TURBULENT_KINETIK_ENERGY], true, true, maxiter,
											   f.alpha[TURBULENT_KINETIK_ENERGY_SL], b, lb, s, ls);
									   }
									   else  if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
										   SoloveichikAlg(&sparseS, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
											   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], true, true, maxiter,
											   f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], b, lb, s, ls);
									   }
									   else if (iVar == NUSHA) {
										   SoloveichikAlg(&sparseS, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], f.maxelm, f.maxbound,
											   rthdsd, f.potent[NUSHA], true, true, maxiter, f.alpha[NUSHA_SL], b, lb, s, ls);
									   }
									   else {
										   SoloveichikAlg(&sparseS, f.slau[iVar], f.slau_bon[iVar], f.maxelm, f.maxbound, rthdsd, 
											   f.potent[iVar], true, true, maxiter, f.alpha[iVar], b, lb, s, ls);
									   }
									   break;
						 default: // Алгоритм Ю.Г.Соловейчика 1993 года. Встроен ILU0 предобуславливатель. 
							       // Данный алгоритм сам освобождает память из под sparseS.
							       // ilu0 из книги Й. Саада - Итерационные методы для больших разреженных СЛАУ.
							       maxiter=limititer;
								   if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
									   SoloveichikAlg(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
										   f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], f.maxelm, f.maxbound,
										   rthdsd, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], true, true, maxiter,
										   f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], b, lb, s, ls);
								   }
								   else if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
									   SoloveichikAlg(&sparseS, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
										   f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], f.maxelm, f.maxbound,
										   rthdsd, f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS], true, true, maxiter,
										   f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], b, lb, s, ls);
								   }
								   else
								   if (iVar == TURBULENT_KINETIK_ENERGY) {
									   SoloveichikAlg(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_SL],
										   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], f.maxelm, f.maxbound,
										   rthdsd, f.potent[TURBULENT_KINETIK_ENERGY], true, true, maxiter,
										   f.alpha[TURBULENT_KINETIK_ENERGY_SL], b, lb, s, ls);
								   }
								   else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
									   SoloveichikAlg(&sparseS, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
										   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], f.maxelm, f.maxbound, 
										   rthdsd, f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], true, true, maxiter,
										   f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], b, lb, s, ls);
								   }
								   else if (iVar == NUSHA) {
									   SoloveichikAlg(&sparseS, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], f.maxelm, f.maxbound, 
										   rthdsd, f.potent[NUSHA], true, true, maxiter, f.alpha[NUSHA_SL], b, lb, s, ls);
								   }
								   else {
									   SoloveichikAlg(&sparseS, f.slau[iVar], f.slau_bon[iVar], f.maxelm, f.maxbound, 
										   rthdsd, f.potent[iVar], true, true, maxiter, f.alpha[iVar], b, lb, s, ls);
								   }
								   break;
						 }
				         break;
			case SOR3DALG: 
				
				         
				         switch (iVar) {
				            case VELOCITY_X_COMPONENT: iVarCor=VXCOR; 
					                  break;
				            case VELOCITY_Y_COMPONENT: iVarCor=VYCOR;
					                  break;
				            case VELOCITY_Z_COMPONENT: iVarCor=VZCOR;
					                  break;
				         }
				         // метод Гаусса-Зейделя используется главным образом для проверки корректности сборки мавтрицы СЛАУ.
				         // данный метод обладает совсем слабой скоростью сходимости и не рекомендован к использованию.
				         // Одна из возможных целей его применения проверка корректности сборки матрицы СЛАУ.
						 if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
							 SOR3D(f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], 
								 f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], 
								 f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], 
								 f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], f.maxelm, f.maxbound,
								 iVar, f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL]);
						 }
						 else if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
							 SOR3D(f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
								 f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
								 f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS], 
								 f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS], 
								 f.maxelm, f.maxbound, iVar, f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL]);
						 }
						 else
						 if (iVar == TURBULENT_KINETIK_ENERGY) {
							 SOR3D(f.slau[TURBULENT_KINETIK_ENERGY_SL],
								 f.slau_bon[TURBULENT_KINETIK_ENERGY_SL],
								 f.potent[TURBULENT_KINETIK_ENERGY], 
								 f.potent[TURBULENT_KINETIK_ENERGY], f.maxelm, f.maxbound,
								 iVar, f.alpha[TURBULENT_KINETIK_ENERGY_SL]);
						 }
						 else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
							 SOR3D(f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
								 f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], 
								 f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA],
								 f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], f.maxelm, f.maxbound, 
								 iVar, f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL]);
						 }
						 else if (iVar == NUSHA) {
							 SOR3D(f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], 
								 f.potent[NUSHA], f.potent[NUSHA], f.maxelm, f.maxbound,
								 iVar, f.alpha[NUSHA_SL]);
						 }
						 else {
							 SOR3D(f.slau[iVar], f.slau_bon[iVar], 
								 f.potent[iVar], f.potent[iVarCor], f.maxelm, f.maxbound,
								 iVar, f.alpha[iVar]);
						 }
						 freeIMatrix(&sparseS);
				         break; 
			case SOLVELRN: // полилинейный метод
				          // Данный метод требует небольшое число итераций, он рекомендован к использованию С. Патанкаром. Но метод
				          // LR1sk обладает гораздо лучшей сходимостью.
				          solveLRn(f.potent[iVar], rthdsd, (f.maxelm+f.maxbound), iVar, maxiter, bprintmessage, false, f.neighbors_for_the_internal_node, f.maxelm, f.slau, f.slau_bon, f.iN, f.id, f.iWE, f.iSN, f.iBT, f.alpha, f.maxbound);
						  freeIMatrix(&sparseS);
				         break;
			case BICG:  // Алгоритм Бисопряжённых градиентов.
				         // последний параметр определяет способ ilu0 разложения
			             // Если true то используется алгоритм из книги Й. Саада
			             // Внимание: ilu0 из книги Саада работает только для матриц с симметричным портретом.
			             // Особенность реализации данного метода в том что он сам освобождает память из под sparseS
				         switch(iPRECOND) {
						 case NOPRECOND: maxiter=limititer; 
							             BiSoprGradCRS(val, col_ind, row_ptr, rthdsd, f.potent[iVar], (f.maxelm+f.maxbound), maxiter); 
							             freeIMatrix(&sparseS);  
							             break; // bug ?
												
						   // с предобуславливателем:
                           case ILU0ITL: maxiter=limititer;
							   if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], 
									   f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], rthdsd, 
									   f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], f.maxelm, f.maxbound, 
									   false, f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], maxiter, b, lb, s, ls);
							   }
							   else if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], 
									   f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], rthdsd,
									   f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS], f.maxelm, f.maxbound,
									   false, f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], maxiter, b, lb, s, ls);
							   }
							   else
							   if (iVar == TURBULENT_KINETIK_ENERGY) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_SL],
									   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], rthdsd,
									   f.potent[TURBULENT_KINETIK_ENERGY], f.maxelm, f.maxbound, 
									   false, f.alpha[TURBULENT_KINETIK_ENERGY_SL], maxiter, b, lb, s, ls);
							   }
							   else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
									   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], rthdsd, 
									   f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], f.maxelm, f.maxbound,
									   false, f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], maxiter, b, lb, s, ls);
							   }
							   else if (iVar == NUSHA) {
								   BiSoprGrad(&sparseS, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], rthdsd, 
									   f.potent[NUSHA], f.maxelm, f.maxbound, false, f.alpha[NUSHA_SL], maxiter, b, lb, s, ls);
							   }
							   else {
								   BiSoprGrad(&sparseS, f.slau[iVar], f.slau_bon[iVar], rthdsd, f.potent[iVar],
									   f.maxelm, f.maxbound, false, f.alpha[iVar], maxiter, b, lb, s, ls);
							   }
							break;
						   case ILU0SAAD: maxiter=limititer;
							   if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
									   f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], rthdsd,
									   f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], f.maxelm, f.maxbound,
									   true, f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], maxiter, b, lb, s, ls);
							   }
							   else if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
									   f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], rthdsd,
									   f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS], f.maxelm, f.maxbound,
									   true, f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], maxiter, b, lb, s, ls);
							   }
							   else
							   if (iVar == TURBULENT_KINETIK_ENERGY) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_SL],
									   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], rthdsd, 
									   f.potent[TURBULENT_KINETIK_ENERGY], f.maxelm, f.maxbound, 
									   true, f.alpha[TURBULENT_KINETIK_ENERGY_SL], maxiter, b, lb, s, ls);
							   }
							   else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], 
									   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], rthdsd,
									   f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], f.maxelm, f.maxbound,
									   true, f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
									   maxiter, b, lb, s, ls);
							   }
							   else if (iVar == NUSHA) {
								   BiSoprGrad(&sparseS, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], 
									   rthdsd, f.potent[NUSHA], f.maxelm, f.maxbound, true, 
									   f.alpha[NUSHA_SL], maxiter, b, lb, s, ls);
							   }
							   else {
								   BiSoprGrad(&sparseS, f.slau[iVar], f.slau_bon[iVar], rthdsd,
									   f.potent[iVar], f.maxelm, f.maxbound, true, f.alpha[iVar],
									   maxiter, b, lb, s, ls);
							   }
							break;
						   default: maxiter=limititer;
							   if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
									   f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], rthdsd,
									   f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], f.maxelm, f.maxbound,
									   true, f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], maxiter, b, lb, s, ls);
							   }
							   else if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
									   f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], rthdsd, 
									   f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS], f.maxelm, f.maxbound, 
									   true, f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], maxiter, b, lb, s, ls);
							   }
							   else
							   if (iVar == TURBULENT_KINETIK_ENERGY) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_SL], 
									   f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], rthdsd, 
									   f.potent[TURBULENT_KINETIK_ENERGY], f.maxelm, f.maxbound,
									   true, f.alpha[TURBULENT_KINETIK_ENERGY_SL], maxiter, b, lb, s, ls);
							   }
							   else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
								   BiSoprGrad(&sparseS, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
									   f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], rthdsd,
									   f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], f.maxelm, f.maxbound,
									   true, f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], maxiter, b, lb, s, ls);
							   }
							   else if (iVar == NUSHA) {
								   BiSoprGrad(&sparseS, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], rthdsd, f.potent[NUSHA], 
									   f.maxelm, f.maxbound, true, f.alpha[NUSHA_SL], maxiter, b, lb, s, ls);
							   }
							   else {
								   BiSoprGrad(&sparseS, f.slau[iVar], f.slau_bon[iVar], rthdsd, f.potent[iVar], 
									   f.maxelm, f.maxbound, true, f.alpha[iVar], maxiter, b, lb, s, ls);
							   }
								   break;
						 }
				         break;
			case BICGSTAB: // Би сопряжённые градиенты со стабилизацией. Просто алгоритм без какого-либо предобуславливания.
			                // В основе лежит CRS формат хранения.
				           
				              	maxiter = limititer;
				             	Bi_CGStabCRS((f.maxelm + f.maxbound), val, col_ind, row_ptr, rthdsd, f.potent[iVar], maxiter);//->//
					            freeIMatrix(&sparseS);
				            							
				          break;
			case LR1SK: 
				
				         // А.А. Фомин, Л.Н. Фомина
			             // Ускорение полилинейного рекуррентного метода в подпространствах Крылова.
                         // Вестник томского государственного университета. Математика и механика №2(14) 2011год.
				         if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
				              LR1sK(f, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], 
								  f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], val, col_ind, row_ptr,
								  f.maxelm, f.maxbound, iVar, rthdsd, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS],
								  maxiter, bprintmessage, bexporttecplot);//->//
				         }
				         else if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
					         LR1sK(f, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], 
								 f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], val, col_ind, row_ptr,
								 f.maxelm, f.maxbound, iVar, rthdsd, f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS],
								 maxiter, bprintmessage, bexporttecplot);//->//
				         }
				         else
				         if (iVar == TURBULENT_KINETIK_ENERGY) {
					         LR1sK(f, f.slau[TURBULENT_KINETIK_ENERGY_SL],
								 f.slau_bon[TURBULENT_KINETIK_ENERGY_SL],
								 val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd,
								 f.potent[TURBULENT_KINETIK_ENERGY], maxiter, bprintmessage,
								 bexporttecplot);//->//
				         }
				         else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
					         LR1sK(f, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
								 f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
								 val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd,
								 f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA],
								 maxiter, bprintmessage, bexporttecplot);//->//
				         }
				         else if (iVar == NUSHA) {
					        LR1sK(f, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, f.potent[NUSHA], maxiter, bprintmessage, bexporttecplot);//->//
				         }
				         else {
					        LR1sK(f, f.slau[iVar], f.slau_bon[iVar], val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, f.potent[iVar], maxiter, bprintmessage, bexporttecplot);//->//
				         }

						 freeIMatrix(&sparseS);

						 if (bexporttecplot) {
							 // Экспорт в программу tecplot360 
					         // в случае обнаружения расходимости.
							 const int ianimate = 0;
					         exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior, ianimate,false,0, b, lb);
						 }

				      break;
			case BICGSTABP:
				      if (iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2) {
						  
					      // LR1SK
						  if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
							  // LR1sK(f, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
							  //val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
							  //maxiter, bprintmessage, bexporttecplot);//->//
							  Lr1sk_up(f, t, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], 
								  f.maxelm, f.maxbound, rthdsd, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], maxiter, 
								  f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], iVar, false);//->//

						  }
						  else if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
							  // LR1sK(f, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
							  //val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], maxiter, 
							  //bprintmessage, bexporttecplot);//->//
							  Lr1sk_up(f, t, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], 
								  f.maxelm, f.maxbound, rthdsd, f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS], maxiter,
								  f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], iVar, false);//->//

						  }
						  else
						  if (iVar == TURBULENT_KINETIK_ENERGY) {
							  // LR1sK(f, f.slau[TURBULENT_KINETIK_ENERGY_SL], f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, f.potent[TURBULENT_KINETIK_ENERGY_SL], maxiter, bprintmessage, bexporttecplot);//->//
							  Lr1sk_up(f, t, f.slau[TURBULENT_KINETIK_ENERGY_SL], f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], f.maxelm, f.maxbound, rthdsd, f.potent[TURBULENT_KINETIK_ENERGY], maxiter, f.alpha[TURBULENT_KINETIK_ENERGY_SL], iVar, false);//->//

						  }
						  else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
							  // LR1sK(f, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], maxiter, bprintmessage, bexporttecplot);//->//
							  Lr1sk_up(f, t, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], f.maxelm, f.maxbound, rthdsd, f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], maxiter, f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], iVar, false);//->//

						  }
						  else if (iVar == NUSHA) {
							  // LR1sK(f, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, f.potent[NUSHA_SL], maxiter, bprintmessage, bexporttecplot);//->//
							  Lr1sk_up(f, t, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], f.maxelm, f.maxbound, rthdsd, f.potent[NUSHA], maxiter, f.alpha[NUSHA_SL], iVar, false);//->//

						  }
						  else {
							  // LR1sK(f, f.slau[iVar], f.slau_bon[iVar], val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, f.potent[iVar], maxiter, bprintmessage, bexporttecplot);//->//
							  Lr1sk_up(f, t, f.slau[iVar], f.slau_bon[iVar], f.maxelm, f.maxbound, rthdsd, f.potent[iVar], maxiter, f.alpha[iVar], iVar, false);//->//
						  }
						  
						  //freeIMatrix(&sparseS);
						  if (val != nullptr) {
							  delete val;
							  val = nullptr;
						  }
						  if (col_ind != nullptr) {
							  delete col_ind;
							  col_ind = nullptr;
						  }
						  if (row_ptr != nullptr) {
							  delete row_ptr;
							  row_ptr = nullptr;
						  }
				      }
				      else {

				      	if (!bBiCGStabSaad) {
						    delete val; delete col_ind; delete row_ptr;
						    val = nullptr; col_ind = nullptr; row_ptr = nullptr;
					    }
				     	// Освобождение памяти из под sparseS происходит внутри метода Bi_CGStab.
					    //dterminatedTResudual=1e-6;
						if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
							Bi_CGStab(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], 
								f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], f.maxelm, f.maxbound,
								rthdsd, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], maxiter,
								f.alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], iVar, m, false, b, lb,
								f.ifrontregulationgl, f.ibackregulationgl, dgx, dgy, dgz, s, ls, inumiter, color, dist_max, w, lw);
						}
						else if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
							Bi_CGStab(&sparseS, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
								f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
								f.maxelm, f.maxbound, rthdsd, f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS],
								maxiter, f.alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL], 
								iVar, m, false, b, lb, f.ifrontregulationgl, f.ibackregulationgl, dgx, dgy, dgz, s, ls, inumiter, color, dist_max, w, lw);
						}
						else
						if (iVar == TURBULENT_KINETIK_ENERGY) {
							Bi_CGStab(&sparseS, f.slau[TURBULENT_KINETIK_ENERGY_SL], 
								f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], f.maxelm, f.maxbound, rthdsd,
								f.potent[TURBULENT_KINETIK_ENERGY], maxiter, f.alpha[TURBULENT_KINETIK_ENERGY_SL],
								iVar, m, false, b, lb, f.ifrontregulationgl, f.ibackregulationgl, dgx, dgy, dgz, s, ls, inumiter, color, dist_max, w, lw);
						}
						else if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
							Bi_CGStab(&sparseS, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], 
								f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], f.maxelm, f.maxbound, rthdsd,
								f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], maxiter, 
								f.alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], iVar, m, false, b, lb, 
								f.ifrontregulationgl, f.ibackregulationgl, dgx, dgy, dgz, s, ls, inumiter, color, dist_max, w, lw);
						}
						else if (iVar==NUSHA) {
							Bi_CGStab(&sparseS, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], f.maxelm, f.maxbound, 
								rthdsd, f.potent[NUSHA], maxiter, f.alpha[NUSHA_SL], iVar, m, false, b, lb, 
								f.ifrontregulationgl, f.ibackregulationgl, dgx, dgy, dgz, s, ls, inumiter, color, dist_max, w, lw);
						}
						else {
							Bi_CGStab(&sparseS, f.slau[iVar], f.slau_bon[iVar], f.maxelm, f.maxbound, 
								rthdsd, f.potent[iVar], maxiter, f.alpha[iVar], iVar, m, false, b, lb, 
								f.ifrontregulationgl, f.ibackregulationgl, dgx, dgy, dgz, s, ls, inumiter, color, dist_max, w, lw);
						}
					}
				      break;
			default: // алгоритм рекомендованный по умолчанию.
				      // А.А. Фомин, Л.Н. Фомина
			          // Ускорение полилинейного рекуррентного метода в подпространствах крылова.
                      // Вестник томского государственного университета. Математика и механика №2(14) 2011год.
				      if (iVar == TURBULENT_KINETIK_ENERGY_STD_K_EPS) {
					        LR1sK(f, f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL], 
								f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
								val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd,
								f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS], 
								maxiter, bprintmessage, bexporttecplot);//->//
				      }
				      else  if (iVar == TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS) {
				        	LR1sK(f, f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
								f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
								val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd, 
								f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS],
								maxiter, bprintmessage, bexporttecplot);//->//
				      }
					  else  if (iVar == TURBULENT_KINETIK_ENERGY) {
					     LR1sK(f, f.slau[TURBULENT_KINETIK_ENERGY_SL],
							 f.slau_bon[TURBULENT_KINETIK_ENERGY_SL],
							 val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd,
							 f.potent[TURBULENT_KINETIK_ENERGY], maxiter, bprintmessage, bexporttecplot);//->//
				      }
				      else  if (iVar == TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
					     LR1sK(f, f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
							 f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL],
							 val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar, rthdsd,
							 f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA],
							 maxiter, bprintmessage, bexporttecplot);//->//
				      }
				      else if (iVar == NUSHA) {
					     LR1sK(f, f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL],
							 val, col_ind, row_ptr, f.maxelm, f.maxbound, iVar,
							 rthdsd, f.potent[NUSHA], maxiter, bprintmessage, bexporttecplot);//->//
				      }
				      else {
					     LR1sK(f, f.slau[iVar], f.slau_bon[iVar], val, col_ind, row_ptr,
							 f.maxelm, f.maxbound, iVar, rthdsd, f.potent[iVar],
							 maxiter, bprintmessage, bexporttecplot);//->//
				      }
					  freeIMatrix(&sparseS);
					  if (bexporttecplot) {
							 // Экспорт в программу tecplot360 
					         // в случае обнаружения расходимости.
						     const int ianimate = 0;
					         exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior, ianimate,false,0, b, lb);
					 }

				      break;
			}

			if (!bBiCGStabSaad) {
			   if (val!=nullptr) {
			      delete val; 
			   }
			   if (col_ind!=nullptr) {
			      delete col_ind;
			   }
			   if (row_ptr!=nullptr) {
				   delete row_ptr;
			   }
			}
			// освобождение памяти
	       // delete rthdsd;
			//getchar(); // debug
			
		}
		else {

		    

			// Существенная экономия оперативной памяти за счёт уничтожения значительной части объектов
			// с последующим их перепостроением с нуля. 
			// Работает только для структурированной сетки.

			bool bmemory_saving = true; // true

			if (!bonly_solid_calculation) {
                // Только в случае решения задачи чистой теплопередачи. 
				bmemory_saving = false;
			}

			if (b_on_adaptive_local_refinement_mesh) {
				// Только для структурированной сетки.
				bmemory_saving = false; 
			}

			if (ireconstruction_free_construct_alloc == 0) {
				// 8 september 2017.
				// Если пользователь хочет это выключить.
				// При использовании полигонов в модели повторное построение 
				// расчётной сетки приводит как правло к новой сетке пока по 
				// непонятным причинам, размер вектора температуры следовательно 
				// меняется но он не меняется т.к. память под поле температур выделялась единожды.
				// Нужна переинтерполяция с сетки на сетку. // TODO переинтерполяция.
				bmemory_saving = false;
			}

			


			// iVar==TEMP
			if (!bBiCGStabSaad) {
			    simplesparsetoCRS(sparseM, val, col_ind, row_ptr, (t.maxelm+t.maxbound));
                simplesparsefree(sparseM,(t.maxelm+t.maxbound));
			}
            
			integer maxiter=4000; // 10 4000 Для возможно несимметричных матриц.

			// Перечень реализованных решателей СЛАУ:
			const integer NOSOLVE=7; // не решать уравнение теплопроводности.
			const integer SOLOVEICHIK=0; // алгоритм Соловейчика.
			const integer SOR3DALG=1; // метод релаксации основанный на методе Гаусса-Зейделя.
			const integer SOLVELRN=2; // полинейный метод.
			const integer BICG=3; // бисопряжённые градиенты.
			const integer BICGSTAB=4; // алгоритм Ван-Дер-Ворста
			const integer LR1SK=5; // Фомин-Фомина данный алгоритм рекомендован к использованию.
			const integer BICGSTABP=6; // Предобусловленный алгоритм Ван дер Ворста. Предобуславливатели из SPARSKIT2.
 
			integer iSOLVER=NOSOLVE; // LR1SK выбор решающего устройства.
			if (inumiter>100) {
				// Здесь мы прследуем следующую цель.
				// Если сразу включить расчёт теплопередачи на первых итерациях SIMPL`е
				// алгоритма то это усилит нелинейность (через свойства материалов) и вероятно может привести
				// к расходимости. Хотя более позднее включение дополнительного уравнения тоже приводит
				// к "удару по невязкам" вероятно это более предпочтительно.
				 maxiter=1000; // 120
                 iSOLVER=LR1SK;
			}
			iSOLVER=LR1SK;
			if (bBiCGStabSaad) {
			   iSOLVER=BICGSTABP;//LR1SK;//BICGSTAB;
			}

			if ((iSOLVER == BICGSTABP) && (iVar == TEMP)&&(iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 2)) {
				// В случае LR1sk солвера освобождать память ненадо.
				// LR1sk использует в результате своей работы t.neighbors_for_the_internal_node, t.slau_bon, t.slau.
				// Поэтому освобождать из под них память нельзя.
				bmemory_saving = false;
			}

			if (bmemory_saving) {
				// Освобождение тех структур данных которые не используются в момент решения СЛАУ.
				free_level1_temp(t);
				free_level1_flow(fglobal, flow_interior);
			}

			// Для сошедшейся задачи подходит уровень среднеквадратических невязок 1.0e-4
            // Если норма Чебышёва то невязка сошедшаяся равна 1.0e-3.
            // Источник опять же мануал по CFX.
			// Однако для задачи трёхмерной чистой теплопроводности установлено,
			// что необходимый уровень среднеквадратиченой невязки скорее 1.0e-5,
			// а уровень 1.0e-4 недостаточен. (это около 287 итераций).
			dterminatedTResudual=1e-7; //t.resLR1sk; // O(h!3)
			 doublereal alpharelaxmethod=1.0; // параметр релаксации.

			// В случае возникновения проблем со сходимостью
			// будет произведён экспорт в программу tecplot360.
			bool bexporttecplot=false;

			if (bconvective) {
				
				switch (iSOLVER) {
				case SOLVELRN: freeIMatrix(&sparseS);
					            delete val; delete col_ind; delete row_ptr;
					            solveLRn_temp(t, t.potent, rthdsd, (t.maxelm + t.maxbound), maxiter, bprintmessage); // полилинейный метод
					           break;
				case BICGSTAB: freeIMatrix(&sparseS);
					            Bi_CGStabCRS((t.maxelm+t.maxbound), val, col_ind, row_ptr, rthdsd, t.potent, maxiter);
                                delete val; delete col_ind; delete row_ptr;
					            break; 
				case BICG:  // последний параметр определяет способ ilu0 разложения
			                 // Если true то используется алгоритм из книги Й. Саада
			                 // Внимание: ilu0 из книги Саада работает только для матриц с симметричным портретом.
			                 // Особенность реализации данного метода в том что он сам освобождает память из под sparseS
			                 // bug ?//
				             // 9-ый параметр коэффициент релаксации == 1.0.
			                 // bug ?//
					         delete val; delete col_ind; delete row_ptr;
					         BiSoprGrad(&sparseS, t.slau, t.slau_bon, rthdsd, t.potent, t.maxelm, t.maxbound, false, 1.0, maxiter, b, lb, s, ls);
							 freeIMatrix(&sparseS);
					        break;
				case LR1SK:  freeIMatrix(&sparseS);				

				              // А.А. Фомин, Л.Н. Фомина
			                  // Ускорение полилинейного рекуррентного метода в подпространствах крылова.
                              // Вестник томского государственного университета. Математика и механика №2(14) 2011год.
			                  LR1sK_temp(t, t.slau, t.slau_bon, val, col_ind, row_ptr, t.maxelm, t.maxbound, rthdsd, t.potent, maxiter, inumiter,bprintmessage,bexporttecplot);
							  delete val; delete col_ind; delete row_ptr;
					        break;	
				case BICGSTABP: if (!bBiCGStabSaad) {
					                delete val; delete col_ind; delete row_ptr;
								 }
					            // Освобождение памяти из под sparseS происходит внутри метода Bi_CGStab.
				                dterminatedTResudual=1e-6;
								alpharelaxmethod = t.alpha; // 1.0
								if (maxiter<m.icount_vel) {
								    m.icount_vel=maxiter;
								}
								if (iswitchsolveramg_vs_BiCGstab_plus_ILU2==2) {
									dterminatedTResudual=1e-8;
									// Специальная версия решателя для уравнения теплопередачи.
									Lr1sk_up(f, t,  t.slau, t.slau_bon,t.maxelm, t.maxbound,   rthdsd, t.potent, maxiter, alpharelaxmethod, TEMP, false);
								}
								else {
									

									dterminatedTResudual=1e-6;
				                   Bi_CGStab(&sparseS, t.slau, t.slau_bon, t.maxelm, t.maxbound, rthdsd, t.potent, maxiter, alpharelaxmethod,TEMP, m, false, b,lb, t.ifrontregulationgl, t.ibackregulationgl, dgx, dgy, dgz,s,ls, inumiter, color, dist_max, w, lw);
								  // doublereal tmax = 0.0;
								  // for (integer i1 = 0; i1<t.maxelm + t.maxbound; i1++) tmax = fmax(tmax, fabs(t.potent[i1]));
								   //printf("apost solver: maximum temperature in default interior is %1.4e\n", tmax);
								   //getchar();
								}
					         break;
				case NOSOLVE: freeIMatrix(&sparseS);
					           delete val; delete col_ind; delete row_ptr; 
							   // не будем решать уравнение теплопроводности так как оно может разойтись если нет теплоотвода.
							   // жидкость течёт на первых итерациях очень слабо и не успевает унести тепло.
					         break;
				}
				
			}
			else 
			{
				// В действительности остаётся только вариант 0 с предобуславливателями Саада.
				const integer iwork_variant=0; // 0 or 1
				//if (bBiCGStabSaad) {
					//   iwork_variant=0; 
				//}

				/*
				// Проверка корректности матрицы СЛАУ:
				for (integer i=0; i<t.maxelm+t.maxbound; i++) {
					if (i<t.maxelm) {
					#if doubleintprecision == 1
						printf("numberCV=%lld ap=%1.4e ae=%1.4e aw=%1.4e an=%1.4e as=%1.4e at=%1.4e ab=%1.4e b=%1.4e\n",i,t.slau[i].ap,t.slau[i].ae,t.slau[i].aw,t.slau[i].an,t.slau[i].as,t.slau[i].at,t.slau[i].ab,t.slau[i].b);
					#else
						printf("numberCV=%d ap=%1.4e ae=%1.4e aw=%1.4e an=%1.4e as=%1.4e at=%1.4e ab=%1.4e b=%1.4e\n",i,t.slau[i].ap,t.slau[i].ae,t.slau[i].aw,t.slau[i].an,t.slau[i].as,t.slau[i].at,t.slau[i].ab,t.slau[i].b);
					#endif
					}
					if (i % 1==0) getchar();
				}
				*/
			
			//SOR3D(t.slau, t.slau_bon, t.potent, t.maxelm, t.maxbound);
			//ICCG(sparseM, rthdsd, t.potent, t.maxelm);
			
			//simplesparsetoCRS(sparseM, val, col_ind, row_ptr, t.maxelm);
			//BiSoprGradCRS(val, col_ind, row_ptr, rthdsd, t.potent, t.maxelm, maxiter);
			
			//equation3DtoCRS(t.slau, t.slau_bon, val, col_ind, row_ptr, t.maxelm, t.maxbound);
			//printf("max memory precalculation ...\n");
			//getchar();
			//free_level2_temp(t);
	        
            //SoloveichikAlgCRS(t.maxelm+t.maxbound, val, col_ind, row_ptr, rthdsd, t.potent, true, maxiter);
				// рабочий код !
				/*
				freeIMatrix(&sparseS); // Внимание !!! нужно обязательно освободить память из под sparseS, иначе будет утечка памяти
			    Bi_CGStabCRS(t.maxelm+t.maxbound, val, col_ind, row_ptr, rthdsd, t.potent, maxiter);
			    delete val; delete col_ind; delete row_ptr;
			//*/
			//


			  // последний параметр определяет способ ilu0 разложения
			  // Если true то используется алгоритм из книги Й. Саада
			     //->//BiSoprGrad(&sparseS, t.slau, t.slau_bon, rthdsd, t.potent, t.maxelm, t.maxbound, true);
			  //SoloveichikAlg(&sparseS, t.slau, t.slau_bon, t.maxelm, t.maxbound, rthdsd, t.potent, true, maxiter);
#if (iwork_variant==0) 
				{
				   //*
				   // Данный код почему то не сходится или очень медленно сходится, гораздо медленнее чем чистый Bi_CGStabCRS.
				   // Данное замечание относится к устаревшей версии кода.
					if (!bBiCGStabSaad) {
				       delete val; delete col_ind; delete row_ptr;
					}
				   doublereal alpharelaxmethod=t.alpha; // 1.0 параметр релаксации.
				   // Освобождение памяти из под sparseS происходит внутри метода Bi_CGStab.
				  // dterminatedTResudual=1e-6;
                  dterminatedTResudual=1.0e-24; // это соответствует евклидовой норме 1.0Е-12.
				   
				   if (iswitchsolveramg_vs_BiCGstab_plus_ILU2==2) {
					    dterminatedTResudual=1e-8;
						Lr1sk_up(f, t,  t.slau, t.slau_bon,t.maxelm, t.maxbound,   rthdsd, t.potent, maxiter, alpharelaxmethod, TEMP, false);
				   }
				   else {

					   //doublereal tmax = -1.0e30;
					   //for (integer i23 = 0; i23 <t.maxelm + t.maxbound; i23++) {

						  // if (t.potent[i23] > tmax) {
							//   tmax = t.potent[i23];
						   //}
					   //}
					   //printf("solve tmax=%e\n", tmax); getchar();

					    dterminatedTResudual=1e-6;
				        Bi_CGStab(&sparseS, t.slau, t.slau_bon, t.maxelm, t.maxbound, rthdsd, t.potent, maxiter, alpharelaxmethod,TEMP, m, false,b,lb, t.ifrontregulationgl, t.ibackregulationgl, dgx, dgy, dgz,s,ls, inumiter, color, dist_max, w, lw);
				   }
				   //*/
				}
#endif
			  //free_level2_temp(t); // delete t.slau;
			  //freeIMatrix(&sparseS);
			     //-->//solveLRn_temp(t, t.potent, rthdsd, (t.maxelm + t.maxbound), maxiter, true); // полилинейный метод

				

#if (iwork_variant==1) 
				{
					// Внимание ! Было проведено тестирование и было выяснено, что данный метод малопригоден для реальных задач.
					// Двигаться надо в сторону SPARSKIT2.
			       //*
			       // А.А. Фомин, Л.Н. Фомина
			       // Ускорение полилинейного рекуррентного метода в подпространствах крылова.
                   // Вестник томского государственного университета. Математика и механика №2(14) 2011год.
				   // рабочий код !
				   freeIMatrix(&sparseS); // Внимание !!! нужно обязательно освободить память из под sparseS, иначе будет утечка памяти
			       LR1sK_temp(t, t.slau, t.slau_bon, val, col_ind, row_ptr, t.maxelm, t.maxbound, rthdsd, t.potent, maxiter, inumiter,bprintmessage,bexporttecplot);
                   delete val; delete col_ind; delete row_ptr; 
			       //*/ 
				}
#endif
			}

			if (bmemory_saving) {
				// Восстановление тех структур данных которые не используются в момент решения СЛАУ.
				bool bextendedprint = false; // печать на граничных узлах рассчитанных полей.
				// надо запомнить в t массивы координат inx, iny, inz.
				integer iCabinetMarker = 0;
				load_TEMPER_and_FLOW(t, fglobal, t.inx_copy, t.iny_copy, t.inz_copy,
					t.xpos_copy, t.ypos_copy, t.zpos_copy, flow_interior,
					b, lb, lw, w, s, ls, lu, my_union, t.operatingtemperature_copy, matlist,
					bextendedprint, dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, true, iCabinetMarker);

				for (integer iu = 0; iu < lu; iu++) {
					integer iup1 = iu + 1;
					load_TEMPER_and_FLOW(my_union[iu].t, my_union[iu].f,
						my_union[iu].inx, my_union[iu].iny, my_union[iu].inz,
						my_union[iu].xpos, my_union[iu].ypos, my_union[iu].zpos,
						my_union[iu].flow_interior,
						b, lb, lw, w, s, ls, lu, my_union, my_union[iu].t.operatingtemperature, matlist, bextendedprint,
						dgx, dgy, dgz, b_on_adaptive_local_refinement_mesh, true, iup1);
				}

			}

			//debug_signal(t, t.operatingtemperature_copy);

		    // освобождение памяти
			
	        // delete rthdsd;

			if (bexporttecplot) {
				// Экспортируем в tecplot только в том случае если мы на обычной а не на АЛИС сетке.
				if (!b_on_adaptive_local_refinement_mesh) {
					const int ianimate = 0;
					exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, ianimate, false,0, b, lb);
				}
			}
	        
		}
         
	}
	
	

	// методы с интегрированным ILU предобуславливателем
	// На симметричной матрице наблюдается резкое ухудшене
	// скорости сходимости и даже расходимость.
	//potent=BiSoprGrad(&sparseS, sparseM,  rthdsd, nullptr, nodes);
	// Для алгоритма Ю.Г. Соловейчика предобуславливание замедляет сходимость ?
	//potent[Temp]=SoloveichikAlg( &sparseS, sparseM, rthdsd, nullptr, true);

} // solve

// В случае мощности зависящей от температуры производит 
// обновление рабочего значения мощности.
void update_power_temperature_depend(SOURCE* s, integer ls, TEMPER &t, BOUND* border_neighbor,
	                                 TEMP_DEP_POWER* gtdps, integer ltdp, doublereal* toldtimestep,
	bool bavgpower, bool bmultipowersourse) {

	// Таблица может быть одна, а источников тепла использующих данную таблицу может быть несколько.
	// Поскольку разные источники тепла имеют разную температуру (мы предполагаем что источник тепла
	// характерезуется всего одним значением температуры, максимальной температурой на поверхности данного источника тепла)
	// то разные источники тепла характеризуются разной мощностью (зависящей от температуры).
	// Это предположение верно, если переменная bmultipowersourse==true; Если false то для всех источников тепла
	// использующих данную таблицу будет присвоена одна и таже мощность рассеяния, соответствующая максимальной температуре
	// среди всех источников тепла использующих данную таблицу.
	//bool bmultipowersourse=true; // true каждый источник со своей мощностью, false - у группы источников одинаковая мощность.

	// следующие величины используются для улучшенного прогноза мощности рассеяния которая будет действовать на текущем шаге по времени.
	// toldtimestep - температура на предыдущем шаге по времени.
	// На основе двух температур toldtimestep и текущей t.potent можно рассчитать более точное
	// значение мощности действующей на текущем временном шаге, чем мощность основывающаяся только на значении t.potent.
	// можно рассчитать как мощность от среднего арифметического температур:
	// toldtimestep && t.potent. (можно также использовать среднее гармоническое значение).
	const integer AREFM=0; // среднее арифметическое.
	const integer GARMONIC=1; // среднее гармоническое.
	integer itypepowerdepend=GARMONIC;
	//bool bavgpower=true; // если true то мощность считается на основе среднего арифметического или среднего гармонического двух температур. 

	doublereal Tmaxtable; // максимальная температура на поверхности источника тепла среди источников тепла использующих данную таблицу.
	doublereal* Tmaxtablei = nullptr;
	Tmaxtablei=new doublereal[ls];
	doublereal Tmaxtableold; // максимальная температура на поверхности источника тепла среди источников тепла использующих данную таблицу с предыдущего временного слоя.
	doublereal* Tmaxtableiold = nullptr;
	Tmaxtableiold=new doublereal[ls];

	// пройдёмся по всем таблицам.
	for (integer itable=0; itable<ltdp; itable++) {
		
		if (bmultipowersourse==false) {
		    Tmaxtable=-1.0e+30;
			Tmaxtableold=-1.0e+30;
		}
		else if (bmultipowersourse==true) {
			for (integer isor=0; isor<ls; isor++) {
				Tmaxtablei[isor]=-1.0e+30;
				Tmaxtableiold[isor]=-1.0e+30;
			}
		}
		
		// пройдёмся по всем граничным контрольным объёмам.
		for (integer inumber=0; inumber<t.maxbound; inumber++) {
            // inumber - номер граничного узла
            if (border_neighbor[inumber].MCB<ls) {
				// источник тепла
				if ((s[border_neighbor[inumber].MCB].bgarber_depend)&&(s[border_neighbor[inumber].MCB].igarber_depend==itable)) {
					// мощность рассеяния задана таблично от температуры и смещения стока,
					// и таблица имеет идентификатор равный itable.

					if (bmultipowersourse==false) {
					    // Вычисление максимальной температуры соответствующей 
					    // данной таблице. (короче для всех источников тепла использующих данную таблицу).
						if (bavgpower==false) {
					        Tmaxtable=fmax(Tmaxtable,t.potent[t.maxelm+inumber]); 
						}
						else if (bavgpower==true) {
							Tmaxtableold=fmax(Tmaxtableold,toldtimestep[t.maxelm+inumber]);
						}
					}
					else if (bmultipowersourse==true) {
						if (bavgpower==false) {
						    Tmaxtablei[border_neighbor[inumber].MCB]=fmax(Tmaxtablei[border_neighbor[inumber].MCB],t.potent[t.maxelm+inumber]);
						}
						else if (bavgpower==true) {
							Tmaxtableiold[border_neighbor[inumber].MCB]=fmax(Tmaxtableiold[border_neighbor[inumber].MCB],toldtimestep[t.maxelm+inumber]);
						}
					}
				}
			}
		}

		// пройдёмся по всем источникам тепла
		for (integer isource=0; isource<ls; isource++) {
			if ((s[isource].bgarber_depend)&&(s[isource].igarber_depend==itable)) {

				if (bmultipowersourse==false) {

					doublereal Tmaxf=Tmaxtable;

					if (bavgpower==true) {
						switch (itypepowerdepend) {
						case AREFM: Tmaxf=0.5*(Tmaxtable+Tmaxtableold);
							         break; 
						case GARMONIC: Tmaxf=2.0*Tmaxtable*Tmaxtableold/(Tmaxtable+Tmaxtableold);
							         break;
						default: Tmaxf=2.0*Tmaxtable*Tmaxtableold/(Tmaxtable+Tmaxtableold);
							break;
						}
					}

				    // если для данного источника мощность зависит от температуры и смещения стока,
				    // а также для задания мощности используется таблица с идентификатором itable, 
				    // то мощность определяется максимальной температурой среди группы источников тепла которым приписана данная таблица.
				    s[isource].power=my_splain_interpol_power_table(gtdps[s[isource].igarber_depend].intemp, 
					                                      gtdps[s[isource].igarber_depend].inoffset_drain, 
														  gtdps[s[isource].igarber_depend].rtemp,
														  gtdps[s[isource].igarber_depend].roffset_drain,
														  gtdps[s[isource].igarber_depend].rpower_table, 
						                                  Tmaxf,
														  s[isource].roperation_offset_drain);
				    s[isource].power*=s[isource].power_multiplyer; // домножение на корректирующий множитель.
				}
				else if (bmultipowersourse==true) {

					doublereal Tmaxf=Tmaxtablei[isource];

					if (bavgpower==true) {
						switch (itypepowerdepend) {
						case AREFM: Tmaxf=0.5*(Tmaxtablei[isource]+Tmaxtableiold[isource]);
							         break; 
						case GARMONIC: Tmaxf=2.0*Tmaxtablei[isource]*Tmaxtableiold[isource]/(Tmaxtablei[isource]+Tmaxtableiold[isource]);
							         break;
						default: Tmaxf=2.0*Tmaxtablei[isource]*Tmaxtableiold[isource]/(Tmaxtablei[isource]+Tmaxtableiold[isource]);
							break;
						}
					}

					// если для данного источника мощность зависит от температуры и смещения стока,
				    // а также для задания мощности используется таблица с идентификатором itable, 
				    // то мощность определяется максимальной температурой на поверхности данного источника тепла.
				    s[isource].power=my_splain_interpol_power_table(gtdps[s[isource].igarber_depend].intemp, 
					                                      gtdps[s[isource].igarber_depend].inoffset_drain, 
														  gtdps[s[isource].igarber_depend].rtemp,
														  gtdps[s[isource].igarber_depend].roffset_drain,
														  gtdps[s[isource].igarber_depend].rpower_table, 
						                                  Tmaxf,
														  s[isource].roperation_offset_drain);
				    s[isource].power*=s[isource].power_multiplyer; // домножение на корректирующий множитель.

				}

			}
		}
	}

	// Освобождение оперативной памяти.
	
	if (Tmaxtablei != nullptr) {
		delete[] Tmaxtablei;
		Tmaxtablei = nullptr;
	}
	if (Tmaxtableiold != nullptr) {
		delete[] Tmaxtableiold;
		Tmaxtableiold = nullptr;
	}

} // update_power_temperature_depend

/*
// Вычисление средних температур в К на границах вакуумных промежитков:
// Это нужно для surface - 2 - surface модели излучения внутри Prism Object.
void update_avg_temperatures(doublereal* &potent, BLOCK &b) {
	if (b.radiation.binternalRadiation) {
		// Если внутри блока работает модель излучения.

		if ((b.radiation.nodelistW != nullptr) && (b.radiation.nodelistE != nullptr) && (b.radiation.nodelistS != nullptr) && (b.radiation.nodelistN != nullptr) && (b.radiation.nodelistB != nullptr) && (b.radiation.nodelistT != nullptr))
		{
			// Сделанор в полном соответствии с Theory Guide ANSYS Fluent.

			// potent - поле темератур.
			doublereal Tavg = 0.0;
			doublereal dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistWsize; i++) {
				doublereal temp= 0.5*((273.15 + potent[b.radiation.nodelistW[i].node1]) + (273.15 + potent[b.radiation.nodelistW[i].node2]));
				Tavg += b.radiation.nodelistW[i].dS*temp*temp*temp*temp;
				dS += b.radiation.nodelistW[i].dS;
			}
			Tavg = pow(Tavg / dS,0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempW is negative...");
			}
			b.radiation.TempW = Tavg;

			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistEsize; i++) {
				doublereal temp= 0.5*((273.15 + potent[b.radiation.nodelistE[i].node1]) + (273.15 + potent[b.radiation.nodelistE[i].node2]));
				Tavg += b.radiation.nodelistE[i].dS*temp*temp*temp*temp;
				dS += b.radiation.nodelistE[i].dS;
			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempE is negative...");
			}
			b.radiation.TempE = Tavg;

			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistSsize; i++) {
				doublereal temp = 0.5*((273.15 + potent[b.radiation.nodelistS[i].node1]) + (273.15 + potent[b.radiation.nodelistS[i].node2]));
				Tavg += b.radiation.nodelistS[i].dS*temp*temp*temp*temp;
				dS += b.radiation.nodelistS[i].dS;
			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempS is negative...");
			}
			b.radiation.TempS = Tavg;

			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistNsize; i++) {
				doublereal temp= 0.5*((273.15 + potent[b.radiation.nodelistN[i].node1]) + (273.15 + potent[b.radiation.nodelistN[i].node2]));
				Tavg += b.radiation.nodelistN[i].dS*temp*temp*temp*temp;
				dS += b.radiation.nodelistN[i].dS;
			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempN is negative...");
			}
			b.radiation.TempN = Tavg;


			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistBsize; i++) {
				doublereal temp = 0.5*((273.15 + potent[b.radiation.nodelistB[i].node1]) + (273.15 + potent[b.radiation.nodelistB[i].node2]));
				Tavg += b.radiation.nodelistB[i].dS*temp*temp*temp*temp;
				dS += b.radiation.nodelistB[i].dS;
			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempB is negative...");
			}
			b.radiation.TempB = Tavg;


			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistTsize; i++) {
				doublereal temp = 0.5*((273.15 + potent[b.radiation.nodelistT[i].node1]) + (273.15 + potent[b.radiation.nodelistT[i].node2]));
				Tavg += b.radiation.nodelistT[i].dS*temp*temp*temp*temp;
				dS += b.radiation.nodelistT[i].dS;
			}
			Tavg = pow(Tavg / dS,0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempT is negative...");
			}
			b.radiation.TempT = Tavg;
		}
		else {
			printf("error memory: can not calculate avg Temperature in bounds vacuum Prism Object. Memory no allocate.\n");
			getchar();
			exit(1);
		}
	}

} // update_avg_temperatures
*/
// Вычисление средних температур в К на границах вакуумных промежутков:
// Это нужно для surface - 2 - surface модели излучения внутри Prism Object.
// Модифицировано с учётом АЛИС сетки 20 сентября 2016.
void update_avg_temperatures(doublereal* &potent, BLOCK &b) {
	if (b.radiation.binternalRadiation) {
		// Если внутри блока работает модель излучения.
		// Мы полностью перешли на Кельвины.

		if ((b.radiation.nodelistW != nullptr) && (b.radiation.nodelistE != nullptr) &&
			(b.radiation.nodelistS != nullptr) && (b.radiation.nodelistN != nullptr) &&
			(b.radiation.nodelistB != nullptr) && (b.radiation.nodelistT != nullptr))
		{
			// Сделано в полном соответствии с Theory Guide ANSYS Fluent.

			// potent - поле темератур.
			doublereal Tavg = 0.0;
			doublereal dS = 0.0;

			// potent в ГРАДУСАХ Цельсия !!!!

			for (integer i = 0; i < b.radiation.nodelistWsize; i++) {
				doublereal temp = 0.0;
				if (b.radiation.nodelistW[i].node21>-1) {
					temp = 0.5*((potent[b.radiation.nodelistW[i].node1]) + (potent[b.radiation.nodelistW[i].node21]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistW[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistW[i].dS1;
				}
				if (b.radiation.nodelistW[i].node22>-1) {
					temp = 0.5*((potent[b.radiation.nodelistW[i].node1]) + (potent[b.radiation.nodelistW[i].node22]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistW[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistW[i].dS2;
				}
				if (b.radiation.nodelistW[i].node23>-1) {
					temp = 0.5*((potent[b.radiation.nodelistW[i].node1]) + (potent[b.radiation.nodelistW[i].node23]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistW[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistW[i].dS3;
				}
				if (b.radiation.nodelistW[i].node24>-1) {
					temp = 0.5*((potent[b.radiation.nodelistW[i].node1]) + (potent[b.radiation.nodelistW[i].node24]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistW[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistW[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempW is negative...");
			}
			b.radiation.TempW = Tavg;

			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistEsize; i++) {

				if (b.radiation.nodelistE[i].node21>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistE[i].node1]) + (potent[b.radiation.nodelistE[i].node21]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistE[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistE[i].dS1;
				}
				if (b.radiation.nodelistE[i].node22>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistE[i].node1]) + (potent[b.radiation.nodelistE[i].node22]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistE[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistE[i].dS2;
				}
				if (b.radiation.nodelistE[i].node23>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistE[i].node1]) + (potent[b.radiation.nodelistE[i].node23]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistE[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistE[i].dS3;
				}
				if (b.radiation.nodelistE[i].node24>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistE[i].node1]) + (potent[b.radiation.nodelistE[i].node24]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistE[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistE[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempE is negative...");
			}
			b.radiation.TempE = Tavg;

			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistSsize; i++) {

				if (b.radiation.nodelistS[i].node21>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistS[i].node1]) + (potent[b.radiation.nodelistS[i].node21]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistS[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistS[i].dS1;
				}
				if (b.radiation.nodelistS[i].node22>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistS[i].node1]) + (potent[b.radiation.nodelistS[i].node22]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistS[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistS[i].dS2;
				}
				if (b.radiation.nodelistS[i].node23>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistS[i].node1]) + (potent[b.radiation.nodelistS[i].node23]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistS[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistS[i].dS3;
				}
				if (b.radiation.nodelistS[i].node24>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistS[i].node1]) + (potent[b.radiation.nodelistS[i].node24]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistS[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistS[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempS is negative...");
			}
			b.radiation.TempS = Tavg;

			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistNsize; i++) {

				if (b.radiation.nodelistN[i].node21>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistN[i].node1]) + (potent[b.radiation.nodelistN[i].node21]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistN[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistN[i].dS1;
				}
				if (b.radiation.nodelistN[i].node22>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistN[i].node1]) + (potent[b.radiation.nodelistN[i].node22]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistN[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistN[i].dS2;
				}
				if (b.radiation.nodelistN[i].node23>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistN[i].node1]) + (potent[b.radiation.nodelistN[i].node23]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistN[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistN[i].dS3;
				}
				if (b.radiation.nodelistN[i].node24>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistN[i].node1]) + (potent[b.radiation.nodelistN[i].node24]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistN[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistN[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempN is negative...");
			}
			b.radiation.TempN = Tavg;


			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistBsize; i++) {

				if (b.radiation.nodelistB[i].node21>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistB[i].node1]) + (potent[b.radiation.nodelistB[i].node21]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistB[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistB[i].dS1;
				}
				if (b.radiation.nodelistB[i].node22>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistB[i].node1]) + (potent[b.radiation.nodelistB[i].node22]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistB[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistB[i].dS2;
				}
				if (b.radiation.nodelistB[i].node23>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistB[i].node1]) + (potent[b.radiation.nodelistB[i].node23]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistB[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistB[i].dS3;
				}
				if (b.radiation.nodelistB[i].node24>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistB[i].node1]) + (potent[b.radiation.nodelistB[i].node24]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistB[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistB[i].dS4;
				}


			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempB is negative...");
			}
			b.radiation.TempB = Tavg;


			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistTsize; i++) {

				if (b.radiation.nodelistT[i].node21>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistT[i].node1]) + (potent[b.radiation.nodelistT[i].node21]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistT[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistT[i].dS1;
				}
				if (b.radiation.nodelistT[i].node22>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistT[i].node1]) + (potent[b.radiation.nodelistT[i].node22]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistT[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistT[i].dS2;
				}
				if (b.radiation.nodelistT[i].node23>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistT[i].node1]) + (potent[b.radiation.nodelistT[i].node23]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistT[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistT[i].dS3;
				}
				if (b.radiation.nodelistT[i].node24>-1) {
					doublereal temp = 0.5*((potent[b.radiation.nodelistT[i].node1]) + (potent[b.radiation.nodelistT[i].node24]));
					temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistT[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistT[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempT is negative...");
			}
			b.radiation.TempT = Tavg;
		}
		else {
			printf("error memory: can not calculate avg Temperature in bounds vacuum Prism Object. Memory no allocate.\n");
			//getchar();
			system("PAUSE");
			exit(1);
		}
	}

} // update_avg_temperatures_old_05.12.2019


  // Вычисление средних температур в К на границах вакуумных промежутков:
  // Это нужно для surface - 2 - surface модели излучения внутри Prism Object.
  // Модифицировано с учётом АЛИС сетки 20 сентября 2016.
void update_avg_temperatures_new(doublereal* &potent, BLOCK &b) {
	if (b.radiation.binternalRadiation) {
		// Если внутри блока работает модель излучения.
		// Мы полностью перешли на Кельвины.

		if ((b.radiation.nodelistW != nullptr) && (b.radiation.nodelistE != nullptr) &&
			(b.radiation.nodelistS != nullptr) && (b.radiation.nodelistN != nullptr) &&
			(b.radiation.nodelistB != nullptr) && (b.radiation.nodelistT != nullptr))
		{
			// Сделано в полном соответствии с Theory Guide ANSYS Fluent.

			// potent - поле темератур.
			doublereal Tavg = 0.0;
			doublereal dS = 0.0;

			// potent в ГРАДУСАХ Цельсия !!!!

			for (integer i = 0; i < b.radiation.nodelistWsize; i++) {
				doublereal temp = 0.0;
				if (b.radiation.nodelistW[i].node21>-1) {
					temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistW[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistW[i].node21])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistW[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistW[i].dS1;
				}
				if (b.radiation.nodelistW[i].node22>-1) {
					temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistW[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistW[i].node22])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistW[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistW[i].dS2;
				}
				if (b.radiation.nodelistW[i].node23>-1) {
					temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistW[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistW[i].node23])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistW[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistW[i].dS3;
				}
				if (b.radiation.nodelistW[i].node24>-1) {
					temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistW[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistW[i].node24])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistW[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistW[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempW is negative...");
			}
			b.radiation.TempW = Tavg;

			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistEsize; i++) {

				if (b.radiation.nodelistE[i].node21>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistE[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistE[i].node21])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistE[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistE[i].dS1;
				}
				if (b.radiation.nodelistE[i].node22>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistE[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistE[i].node22])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistE[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistE[i].dS2;
				}
				if (b.radiation.nodelistE[i].node23>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistE[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistE[i].node23])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistE[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistE[i].dS3;
				}
				if (b.radiation.nodelistE[i].node24>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistE[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistE[i].node24])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistE[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistE[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempE is negative...");
			}
			b.radiation.TempE = Tavg;

			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistSsize; i++) {

				if (b.radiation.nodelistS[i].node21>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistS[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistS[i].node21])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistS[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistS[i].dS1;
				}
				if (b.radiation.nodelistS[i].node22>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistS[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistS[i].node22])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistS[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistS[i].dS2;
				}
				if (b.radiation.nodelistS[i].node23>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistS[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistS[i].node23])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistS[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistS[i].dS3;
				}
				if (b.radiation.nodelistS[i].node24>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistS[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistS[i].node24])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistS[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistS[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempS is negative...");
			}
			b.radiation.TempS = Tavg;

			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistNsize; i++) {

				if (b.radiation.nodelistN[i].node21>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistN[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistN[i].node21])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistN[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistN[i].dS1;
				}
				if (b.radiation.nodelistN[i].node22>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistN[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistN[i].node22])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistN[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistN[i].dS2;
				}
				if (b.radiation.nodelistN[i].node23>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistN[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistN[i].node23])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistN[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistN[i].dS3;
				}
				if (b.radiation.nodelistN[i].node24>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistN[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistN[i].node24])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistN[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistN[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempN is negative...");
			}
			b.radiation.TempN = Tavg;


			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistBsize; i++) {

				if (b.radiation.nodelistB[i].node21>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistB[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistB[i].node21])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistB[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistB[i].dS1;
				}
				if (b.radiation.nodelistB[i].node22>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistB[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistB[i].node22])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistB[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistB[i].dS2;
				}
				if (b.radiation.nodelistB[i].node23>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistB[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistB[i].node23])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistB[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistB[i].dS3;
				}
				if (b.radiation.nodelistB[i].node24>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistB[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistB[i].node24])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistB[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistB[i].dS4;
				}


			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempB is negative...");
			}
			b.radiation.TempB = Tavg;


			Tavg = 0.0;
			dS = 0.0;

			for (integer i = 0; i < b.radiation.nodelistTsize; i++) {

				if (b.radiation.nodelistT[i].node21>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistT[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistT[i].node21])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistT[i].dS1*temp*temp*temp*temp;
					dS += b.radiation.nodelistT[i].dS1;
				}
				if (b.radiation.nodelistT[i].node22>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistT[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistT[i].node22])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistT[i].dS2*temp*temp*temp*temp;
					dS += b.radiation.nodelistT[i].dS2;
				}
				if (b.radiation.nodelistT[i].node23>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistT[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistT[i].node23])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistT[i].dS3*temp*temp*temp*temp;
					dS += b.radiation.nodelistT[i].dS3;
				}
				if (b.radiation.nodelistT[i].node24>-1) {
					doublereal temp = 1.0/((1.0/(273.15+potent[b.radiation.nodelistT[i].node1])) + (1.0/(273.15+potent[b.radiation.nodelistT[i].node24])));
					//temp += 273.15; // перевод в Кельвины.
					Tavg += b.radiation.nodelistT[i].dS4*temp*temp*temp*temp;
					dS += b.radiation.nodelistT[i].dS4;
				}

			}
			Tavg = pow(Tavg / dS, 0.25);
			if (Tavg < 0.0) {
				Tavg = 0.0; // Нарушение физического смысла.
				printf("error: TempT is negative...");
			}
			b.radiation.TempT = Tavg;
		}
		else {
			printf("error memory: can not calculate avg Temperature in bounds vacuum Prism Object. Memory no allocate.\n");
			//getchar();
			system("PAUSE");
			exit(1);
		}
	}

} // update_avg_temperatures



void update_Stefan_Bolcman_condition_double_vacuum_PRISM(WALL* w, integer lw, integer ls, integer maxbound, BOUND* border_neighbor,
	doublereal* potent, doublereal** prop_b, TOCHKA* pa, integer** nvtx, integer maxelm)
{
	for (integer inumber = 0; inumber < maxbound; inumber++) {
		if ((bBlockStefanBolcman && (
			(((border_neighbor[inumber].MCB < (ls + lw)) &&
			(border_neighbor[inumber].MCB >= ls) && 
				(w[border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY))))))
		{
		    

			doublereal dl=1.0,  dS=1.0; // геометрические параметры
			

			// внутренняя нормаль
			switch (border_neighbor[inumber].Norm) {
			case E_SIDE:
				dl = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
				dS = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[1][border_neighbor[inumber].iI] - 1].y;
				dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

								
				break;

			case N_SIDE:
				dl = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y;
				dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
				dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

				
				break;

			case T_SIDE:
				dl = pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z;
				dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
				dS *= (pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y); // площадь грани


				break;

			case W_SIDE:
				dl = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
				dS = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[1][border_neighbor[inumber].iI] - 1].y;
				dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

			
				break;

			case S_SIDE:
				dl = pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y;
				dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
				dS *= (pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z); // площадь грани

				break;

			case B_SIDE:
				dl = pa[nvtx[4][border_neighbor[inumber].iI] - 1].z - pa[nvtx[0][border_neighbor[inumber].iI] - 1].z;
				dS = pa[nvtx[1][border_neighbor[inumber].iI] - 1].x - pa[nvtx[0][border_neighbor[inumber].iI] - 1].x;
				dS *= (pa[nvtx[2][border_neighbor[inumber].iI] - 1].y - pa[nvtx[0][border_neighbor[inumber].iI] - 1].y); // площадь грани
			
			
				break;

			} // end switch

			doublereal lamB = prop_b[LAM][border_neighbor[inumber].iB - maxelm];

			// Модуль плотности кондукционного теплового потока на границе.
			//doublereal density_heat_flux = lamB*fabs(potent[border_neighbor[inumber].iB]-potent[border_neighbor[inumber].iI])/dl;
			doublereal density_heat_flux = lamB * (potent[border_neighbor[inumber].iB] - potent[border_neighbor[inumber].iI]) / dl;
			doublereal sigma = STEFAN_BOLCMAN_CONST;
			doublereal r1 = density_heat_flux / (sigma*w[border_neighbor[inumber].MCB - ls].emissivity);
			doublereal r2 = (273.15+w[border_neighbor[inumber].MCB - ls].Tamb);
			r1 += r2*r2*r2*r2;
			w[border_neighbor[inumber].MCB - ls].Tamb = sqrt(sqrt(r1))-273.15;
			// Мы обновили температуру к которой надо стремиться.
		}

	}
}

// возвращает максимальный элемент в массиве.
template<typename doublerealT>
doublerealT get_max_array_elm(doublerealT *& a, const integer n) {

#ifdef _OPENMP
	// многопоточная версия
	doublerealT dmax = -1.0e27;

#pragma omp parallel
	{
		doublerealT dmax_loc = -1.0e27;

#pragma omp for 
		for (integer i=0; i < n; i++) {
			doublerealT z = a[i];
			if (z > dmax_loc) dmax_loc = z;
        }

#pragma omp critical
		{
			if (dmax_loc > dmax) {
				dmax = dmax_loc;
			}
		}
	}

	return dmax;
#else
	// однопточная версия
	doublerealT dmax = -1.0e27;
	for (integer i=0; i < n; i++) {
		if (a[i] > dmax) dmax = a[i];
	}
	return dmax;
#endif

}

// возвращает максимальный элемент в массиве.
template<typename doublerealT>
doublerealT get_max_array_elm(doublerealT *& a1, doublerealT *& a2, const integer n) {

#ifdef _OPENMP
	// многопоточная версия
	doublerealT dmax = -1.0e27;

#pragma omp parallel
	{
		doublerealT dmax_loc = -1.0e27;

#pragma omp for 
		for (integer i = 0; i < n; i++) {
			doublerealT z = fabs(a1[i]-a2[i]);
			if (z > dmax_loc) dmax_loc = z;
		}

#pragma omp critical
		{
			if (dmax_loc > dmax) {
				dmax = dmax_loc;
			}
		}
	}

	return dmax;
#else
	// однопточная версия
	doublerealT dmax = -1.0e27;
	for (integer i = 0; i < n; i++) {
		doublerealT z = fabs(a1[i]-a2[i]);
		if (z > dmax) dmax = z;
	}
	return dmax;
#endif

}

// возвращает максимальный элемент в массиве.
template<typename doublerealT>
doublerealT get_min_array_elm(doublerealT *& a, const integer n) {

#ifdef _OPENMP
	// многопоточная версия
	doublerealT dmin = 1.0e37;

#pragma omp parallel
	{
		doublerealT dmin_loc = 1.0e37;

#pragma omp for 
		for (integer i=0; i < n; i++) {
			doublerealT z = a[i];
			if (z < dmin_loc) dmin_loc = z;
		}

#pragma omp critical
		{
			if (dmin_loc < dmin) {
				dmin = dmin_loc;
			}
		}
	}

	return dmin;
#else
	// однопточная версия
	doublerealT dmin = 1.0e37;
	for (integer i=0; i < n; i++) {
		if (a[i] < dmin) dmin = a[i];
	}
	return dmin;
#endif

}

// Если уравнение теплопроводности нелинейно, т.е.
// коэффициенты теплоёмкости и тепловодности входящие в него 
// зависят от температуры то следует применять данный решатель, 
// представленный здесь ниже по тексту:
// Примечание: более того, даже если свойства материалов постоянны то
// при применении высокого порядка аппроксимации теплового потока на границе
// порядка выше первого (второй BETA_PRECISION==1.33333, третий BETA_PRECISION==1.2) нужно применять
// данный нелинейный солвер т.к. в случае аппроксимации высокого порядка коэффициенты
// дискретного аналога зависят от текущей температуры и нужны итерации чтобы справиться с нелинейностью.
void solve_nonlinear_temp(FLOW &f, FLOW* &fglobal, TEMPER &t, doublereal** &rhie_chow,
	BLOCK* b, integer lb, SOURCE* s, integer ls, WALL* w, integer lw,
	doublereal dbeta, integer flow_interior, bool bconvective,
	doublereal* toldtimestep, doublereal tauparam, doublereal tauparamold, bool btimedep,
	TPROP* matlist, integer inumiter, bool bprintmessage,
	TEMP_DEP_POWER* gtdps, integer ltdp, doublereal  poweron_multiplier_sequence,
	QuickMemVorst& m, doublereal** speedoldtimestep, doublereal** mfoldtimestep,
	integer lu, UNION* &my_union, integer* &color, integer dist_max) {



	//if (bconvective) {
		//printf("bconvective. Ok.\n");
	//}
	//else {
		//printf("SOLID STATE only.\n");
	//}
	//getchar();

	// Таблица может быть одна, а источников тепла использующих данную таблицу может быть несколько.
	// Поскольку разные источники тепла имеют разную температуру (мы предполагаем что источник тепла
	// характерезуется всего одним значением температуры, максимальной температурой на поверхности данного источника тепла)
	// то разные источники тепла характеризуются разной мощностью (зависящей от температуры).
	// Это предположение верно, если переменная bmultipowersourse==true; Если false то для всех источников тепла
	// использующих данную таблицу будет присвоена одна и таже мощность рассеяния, соответствующая максимальной температуре
	// среди всех источников тепла использующих данную таблицу.
	// для функции update_power_temperature_depend().
	bool bmultipowersourse=true; // true каждый источник со своей мощностью, false - у группы источников одинаковая мощность.

	doublereal power_diss_message_06_10_2018 = 0.0;

	// Проверка есть ли файл с распределением начальной скорости.	
	
	FILE* fp_inicialization_data = nullptr;
#ifdef MINGW_COMPILLER
	int err_inicialization_data = 0;// используется далее по коду. Код нужен.
	fp_inicialization_data = fopen64("load.txt", "r");
	if (fp_inicialization_data == nullptr) {
		// Ошибка открытия файла.
		err_inicialization_data = 1;
    }
#else
	errno_t err_inicialization_data = 0;// используется далее по коду. Код нужен.
	err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif
	if (0 == err_inicialization_data) {
		fclose(fp_inicialization_data);
	}

	t.alpha = 1.0;

	doublereal balancet = 0.0;
	if (bglobal_first_start_radiation) {
		doublereal pdiss = 0.0;

#pragma omp parallel for reduction(+:pdiss)
		for (integer i = 0; i < ls; i++) {
			if (s[i].power < 0.0) {
				//printf("warning source [%lld] is negative power = %e\n",i, s[i].power);
				std::cout << "warning source [" << i << "] is negative power = " << s[i].power << std::endl;
			}
			pdiss += s[i].power;
		}
		//for (integer i = 0; i < lb; i++) {
			//pdiss += b[i].Sc*(fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS));
		//}
		// 19 november 2016.
		// Обновление мощности тепловыделения во всех внутренних узлах.

#pragma omp parallel for reduction(+:pdiss)
		for (integer i47 = 0; i47 < t.maxelm; i47++) {
			// Скорость в том что значение не вычисляется как раньше а просто хранится.
			integer ib = t.whot_is_block[i47];
			
			t.Sc[i47] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, t.potent[i47]);
			// вычисление размеров текущего контрольного объёма:
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
			volume3D(i47, t.nvtx, t.pa, dx, dy, dz);
			if (t.Sc[i47] * dx*dy*dz < 0.0) {
				//printf("ERROR!!!  control volume [%lld] is negative power = %e\n", i47, t.Sc[i47] * dx*dy*dz);
				std::cout << "ERROR!!!  control volume [" << i47 << "] is negative power =" << (t.Sc[i47] * dx*dy*dz) << std::endl;
				//system("PAUSE");
			}
			pdiss += t.Sc[i47] * dx*dy*dz;
			/*
			if ((ib >= 114)&&(ib<=120)) {
				// debug
				printf("ib=%lld i47=%lld t.Sc=%e dx=%e dy=%e dz=%e\n", ib, i47, t.Sc[i47],dx,dy,dz);
				printf("n_Sc=%lld TSc=%e Sc=%e T=%e\n",b[ib].n_Sc, b[ib].temp_Sc[0], b[ib].arr_Sc[0], t.potent[i47]);
				system("pause");
			}
			*/
		}
		//printf("power generation is equal=%e\n",pdiss);
		std::cout << "power generation is equal=" << pdiss << std::endl;
		if (fabs(d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL - pdiss)
			> 0.01 * d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL) {

			// Корректируем значение тепловой мощности так чтобы она в 
			// точности была равна заданной пользователем.
			doublereal m_power = 1.0;
			if (fabs(pdiss) > 1.0e-30) {
				m_power = d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL / pdiss;
				printf("Correction power: m_power=%e\n", m_power);
				// m_power - множитель на который нужно домнодить ненулевую внутрипрограммную
				// тепловую мощность, чтобы в точности получить заданную пользователем мощность.
				pdiss = 0.0;

#pragma omp parallel for reduction(+:pdiss)
				for (integer i = 0; i < ls; i++) {
					if (s[i].power < 0.0) {
						//printf("warning source [%lld] is negative power = %e\n",i, s[i].power);
						std::cout << "warning source [" << i << "] is negative power = " << s[i].power << std::endl;
					}
					s[i].power *= m_power;
					pdiss += s[i].power;
				}
				//for (integer i = 0; i < lb; i++) {
				//b[i].Sc*=m_power;
				//pdiss += b[i].Sc*(fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS));
				//}
				// 19 november 2016.
				// Обновление мощности тепловыделения во всех внутренних узлах.

#pragma omp parallel for 
				for (int ib = 1; ib < lb; ib++) {
					for (int i23 = 0; i23 < b[ib].n_Sc; i23++) {
						if (fabs(b[ib].arr_Sc[i23]) > 1.0e-30) {
							//printf("%e ", b[ib].arr_Sc[i23]);
							b[ib].arr_Sc[i23] *= m_power;
							//printf("%e \n", b[ib].arr_Sc[i23]);
							//getchar();
						}

					}
				}

#pragma omp parallel for reduction(+:pdiss)
				for (integer i47 = 0; i47 < t.maxelm; i47++) {
					// Скорость в том что значение не вычисляется как раньше а просто хранится.
					integer ib = t.whot_is_block[i47];

					

					/*if (t.Sc[i47] > 0.0) {
						printf("%e ", t.Sc[i47]);
						doublereal t12= get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, t.potent[i47]);
						printf("%e ", t12);
						getchar();
					}*/
					t.Sc[i47] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, t.potent[i47]);
					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
					volume3D(i47, t.nvtx, t.pa, dx, dy, dz);
					if (t.Sc[i47] * dx*dy*dz < 0.0) {
						//printf("ERROR!!!  control volume [%lld] is negative power = %e\n", i47, t.Sc[i47] * dx*dy*dz);
						std::cout << "ERROR!!!  control volume [" << i47 << "] is negative power =" << (t.Sc[i47] * dx*dy*dz) << std::endl;
						//system("PAUSE");
					}
					pdiss += t.Sc[i47] * dx*dy*dz;
					/*
					if ((ib >= 114)&&(ib<=120)) {
					// debug
					printf("ib=%lld i47=%lld t.Sc=%e dx=%e dy=%e dz=%e\n", ib, i47, t.Sc[i47],dx,dy,dz);
					printf("n_Sc=%lld TSc=%e Sc=%e T=%e\n",b[ib].n_Sc, b[ib].temp_Sc[0], b[ib].arr_Sc[0], t.potent[i47]);
					system("pause");
					}
					*/
				}
			}

			if (fabs(d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL - pdiss)
		          > 0.2 * d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL) {
				// Проблемы при построении модели. Возможна сильная разномасштабность геометрии.
				//printf("Apriory Pdiss=%e\n", d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL);
				std::cout << "Apriory Pdiss=" << d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL << std::endl;
				printf("FATAL ERROR!!! Your model is incorrect. Power leak.\n");
				printf("Please send you message on kirill7785@mail.ru\n");
				system("pause");
				exit(1);

			}
			else {
				// Мощности отличаются менее чем на 10% прощаем, но пишем диагностическое предупреждение.
				std::cout << "Apriory Pdiss=" << d_GLOBAL_POWER_HEAT_GENERATION_IN_CURRENT_MODEL << std::endl;
				printf("WARNING!!! Your model is incorrect. Power leak <16%%.\n");
				printf("Please send you message(model) on kirill7785@mail.ru\n");
				printf("May be your model is incorrect. WARNING!!!\n");
			}
			
		}
		if (pdiss>0.0) {
			doublereal square_bolc = 0.0;
			doublereal emissivity = 1.0;

#pragma omp parallel for reduction(+:square_bolc) 
			for (integer i = 0; i < lw; i++) {
				if (w[i].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) {
					switch (w[i].iPlane) {
					  case XY_PLANE: square_bolc += fabs(w[i].g.xE - w[i].g.xS)*fabs(w[i].g.yE - w[i].g.yS); break;
					  case XZ_PLANE: square_bolc += fabs(w[i].g.xE - w[i].g.xS)*fabs(w[i].g.zE - w[i].g.zS); break;
					  case YZ_PLANE: square_bolc += fabs(w[i].g.yE - w[i].g.yS)*fabs(w[i].g.zE - w[i].g.zS); break;
					}
					// Здесь мы предполагаем что на всех излучающих поверхностях излучающая способность одна и таже.
					// Если это не так то возникнет ошибка.
				}
			}
			for (integer i = 0; i < lw; i++) {
				if (w[i].ifamily == WALL_BOUNDARY_CONDITION::STEFAN_BOLCMAN_FAMILY) {
					emissivity = w[i].emissivity;
				}
			}

			if (fabs(square_bolc)>1e-23) {
				//printf("Pdiss=%e, S=%e\n",pdiss, square_bolc);
				std::cout << "Pdiss=" << pdiss << ", S=" << square_bolc << std::endl;
				balancet = sqrt(sqrt((pdiss / (square_bolc * STEFAN_BOLCMAN_CONST*emissivity)))) - 273.15;
				//printf("balance temperature =%f\n", balancet);
				std::cout << "balance temperature =" << balancet << std::endl;
				//t.alpha = 0.5;
				//balancet = 50.0;
			}
			else {
				printf("square_bolc is zero!!!\n");
			}
		}
		else {
			// идентифицирована отрицательная мощность тепловыделения.
			bool bDirichlet = false;
#pragma omp parallel for reduction(|| : bDirichlet)
			for (integer i = 0; i < lw; i++) {
				if (w[i].ifamily == WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY) {
					// условие идеального теплоотвода обнаружено.
					
					bDirichlet = bDirichlet || true;
					
				}
			}
			if (!bDirichlet) {
				printf("negative power and the lack of Dirichlet conditions \n");
				//getchar();
				balancet = -272.15;
			}
		}
		power_diss_message_06_10_2018 = pdiss;
	}
	else {
		doublereal pdiss = 0.0;
		
#pragma omp parallel for reduction(+:pdiss)
		for (integer i = 0; i < ls; i++) {
			if (s[i].power < 0.0) {
				//printf("warning source [%lld] is negative power = %e\n", i, s[i].power);
				std::cout << "warning source [" << i << "] is negative power = " << s[i].power << std::endl;
			}
			pdiss += s[i].power;
		}
		//for (integer i = 0; i < lb; i++) {
			//pdiss += b[i].Sc*(fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS));
		//}
		// 19 november 2016.
		// Обновление мощности тепловыделения во всех внутренних узлах.
#pragma omp parallel for reduction(+:pdiss)
		for (integer i47 = 0; i47 < t.maxelm; i47++) {
			// Скорость в том что значение не вычисляется как раньше а просто хранится.
			integer ib = t.whot_is_block[i47];
			t.Sc[i47] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, t.potent[i47]);
			// вычисление размеров текущего контрольного объёма:
			doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
			volume3D(i47, t.nvtx, t.pa, dx, dy, dz);
			if (t.Sc[i47] * dx*dy*dz < 0.0) {
				//printf("ERROR!!!  control volume [%lld] is negative power = %e\n", i47, t.Sc[i47] * dx*dy*dz);
				std::cout << "ERROR!!!  control volume [" << i47
					<< "] is negative power = " << (t.Sc[i47] * dx*dy*dz) << std::endl;
				system("PAUSE");
			}
			pdiss += t.Sc[i47] * dx*dy*dz;
		}
		power_diss_message_06_10_2018 = pdiss;
		printf("no Stefan - Bolcman boundary condition...\n");
	}
	//system("pause");
	
	power_diss_message_06_10_2018*=poweron_multiplier_sequence;  // Включено или выключено.

	// К этому значению источникового члена мы будем релаксировать.
#pragma omp parallel for 
	for (integer i_init = 0; i_init < t.maxelm; i_init++) {
		bsource_term_radiation_for_relax[i_init] = 0.0;
	}

	

	if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {
		//printf("film coefficient=%e, operating_temperature=%f\n", film_coefficient, operating_temperature_for_film_coeff);
		std::cout << "film coefficient=" << film_coefficient 
		<< ", operating_temperature=" << operating_temperature_for_film_coeff 
		<< std::endl;
		//t.alpha = 0.8; // по видимому нужна нижняя релаксация.
		//printf("temperature relax factor is equal %e\n", t.alpha);
		std::cout << "temperature relax factor is equal " << t.alpha << std::endl;
		// system("pause");
	}

	if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::STEFAN_BOLCMAN_BC) {
		//printf(" operating_temperature=%f\n", operating_temperature_for_film_coeff);
		std::cout << " operating_temperature=" << operating_temperature_for_film_coeff << std::endl;
		//t.alpha = 0.8; // по видимому нужна нижняя релаксация.
		//printf("temperature relax factor is equal %e\n", t.alpha);
		std::cout << "temperature relax factor is equal " << t.alpha << std::endl;
		// system("pause");
	}

	if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::MIX_CONDITION_BC) {
		//printf("film coefficient=%e, operating_temperature=%f\n", film_coefficient, operating_temperature_for_film_coeff);
		std::cout << "film coefficient=" << film_coefficient << ", operating_temperature=" << operating_temperature_for_film_coeff << std::endl;
		//t.alpha = 0.8; // по видимому нужна нижняя релаксация.
		//printf("temperature relax factor is equal %e\n", t.alpha);
		std::cout << "temperature relax factor is equal " << t.alpha << std::endl;
		// system("pause");
	}

	doublereal RCh = 1.0;

	//if (bfirst_start_nonlinear_process) {
	// Во избежании расходимости по начальному условию в двойном 
	// вакуумном промежутке.
	//bfirst_start_nonlinear_process = false;
	// В диапазоне от 1 до 25 градусов.
	// for (integer i_init = 0; i_init < t.maxelm; i_init++) t.potent[i_init] = rand() % 2 + 1;
	//}

	

	// начальное поле температур определяется константой в ноль градусов Цельсия во всей расчётной области.
	// Это определено в файле constr_struct.c.
	bool bupdateproperties = true; // начальные свойства при температуре ноль (0) градусов Цельсия.
	doublereal deltat = 100.0; // начальная разница между итерациями 100 градусов Цельсия.
	doublereal *told = nullptr;
	doublereal *told_nonlinear=nullptr;
	told = new doublereal[t.maxelm + t.maxbound];
	told_nonlinear = new doublereal[t.maxelm + t.maxbound];
	integer i = 0;
	doublereal res = 0.0;

	bool bRichman = false;

#pragma omp parallel for reduction(|| : bRichman)
	for (integer inumber = 0; inumber < t.maxbound; inumber++) {
			
			//std::cout << t.border_neighbor[inumber].MCB << " " << ls + lw << " " << w[t.border_neighbor[inumber].MCB - ls].ifamily << std::endl;
			if (adiabatic_vs_heat_transfer_coeff == DEFAULT_CABINET_BOUNDARY_CONDITION::NEWTON_RICHMAN_BC) {

				bRichman = bRichman || true;
				
				//break;
			}
			
			if ((t.border_neighbor[inumber].MCB < (ls + lw)) &&
				(t.border_neighbor[inumber].MCB >= ls) &&
				(w[t.border_neighbor[inumber].MCB - ls].ifamily == WALL_BOUNDARY_CONDITION::NEWTON_RICHMAN_FAMILY)) {
				// на твёрдой стенке задано условие Ньютона Рихмана.

					bRichman = bRichman || true;
					
					// break;
				}		
	}

#pragma omp parallel for shared (t,told) private (i) firstprivate(bglobal_first_start_radiation, bRichman) schedule (guided)
	for (i = 0; i < t.maxelm + t.maxbound; i++) {
		//t.potent[i] = -269.0;
		told[i] = t.potent[i]; // копирование
		if (bglobal_first_start_radiation) {
			//printf("bglobal_first_start_radiation %e\n", balancet); getchar();
			told[i] = balancet;
			t.potent[i] = balancet;

			if (bRichman) {
				//std::cout << "incomming" << std::endl;
				told[i] = operating_temperature_for_film_coeff;
				t.potent[i] = operating_temperature_for_film_coeff;
			}

			//if (bdouble_vacuum_PRISM) {
				//if (i < t.maxelm) {
					//printf("bdouble vacuum prism\n");
					//system("pause");
					//told[i] = 0.1*(rand() % 10) + balancet;
					//t.potent[i] = 0.1*(rand() % 10) + balancet;
				//}
			//}
		}
		told_nonlinear[i]=told[i];
	}


	//{
		//doublereal tmax = -1.0e30;
		//for (integer i23 = 0; i23 <t.maxelm + t.maxbound; i23++) {

			//if (t.potent[i23] > tmax) {
				//tmax = t.potent[i23];
			//}
		//}
		//printf("solve1 tmax=%e\n", tmax); getchar();
	//}

	integer ic = 1;
	doublereal tmax = 0.0;

	bool bfreeflag = m.bsignalfreeCRSt;
	if (bfreeflag) {
		// Сразу память не освобождаем а немного погодя.
		m.bsignalfreeCRSt = false;
	}


	doublereal* rthdsdt = nullptr;
	rthdsdt = new doublereal[t.maxelm + t.maxbound];

	integer iprohod = 1;
	if (bdouble_vacuum_PRISM) iprohod = 1; // На первом проходе условия Стефана - Больцмана блокируются и заменяются на термостабильнуую стенку.
	for (integer iprohodtek = 0; iprohodtek < iprohod; iprohodtek++) {

		deltat = 100.0; // для того чтобы снова запустить цикл.

		//if ((iprohodtek == 0) && (iprohod == 2)) {
		if (bdouble_vacuum_PRISM) {
			// Устанавливаем Блокировку Стефана - Больцмана.
			//bBlockStefanBolcman = true;
		}
		//}
		//else {
		// Снимаем блокировку Стефана - Больцмана.
		//	bBlockStefanBolcman = false;
		//}
		// Уравнение теплопроводности в данном случае нелинейно.
		// Нелинейность вызвана во первых свойствами материала зависящими от температуры,
		// во вторых мощностью зависящей от температуры. Это приводит к тому что внутри шага по времени
		// требуется организовать итерационный процесс.
		//integer i87 = 0;
		integer ibreak_counter_25_07_2017 = 0;
		doublereal fporogmax = -1.0e30;


		doublereal procent_porog=0.05;
		//if ((1 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) &&
		  //    (Non_Linear_amg1r5 == stabilization_amg1r5_algorithm)) {
				 // procent_porog=0.01;
		//}

		// не менее 10 итераций.
		while (deltat > procent_porog*fporogmax) {
			ibreak_counter_25_07_2017++;
			// Досрочный выход из итерационного процесса.
			//if ((err_inicialization_data == 0)&&(ibreak_counter_25_07_2017 > 19)) break;
			// На некоторых задачах разница температур между прогонами не может стать менее 7С поэтому нужен досрочный выход.
			if ((0 == err_inicialization_data) && (ibreak_counter_25_07_2017 > 6)&&
				((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_TEMPERATURE))) break;
			// 25.07.2017 Схема SuperC почему-то плохо сходится,
			// температура основания корпуса не устанавливается, нестановится меньше 0.5 град С.
			// Решалось уравнение конвекции-диффузии для модуля ВУМ на радиаторе водяного 
			// охлаждения 3л/мин  10М узлов с предварительно рассчитанным в отдельной 
			// задаче полем скорости.
			// Схема Super C обеспечивает точность нахождения поля температуры, т.к. схема UDS даёт заниженные на 25% температуры из-за 
			// сильной искусственной диффузии.

					/*if (i87 == 0) {
						bBlockStefanBolcman = true;
					}
					else if (i87==1) {
						bBlockStefanBolcman = false;
					}
					else*/
					//if (i87%2 == 0) {
						// На старте мы блокируем Стефана Больцмана дав сойтись лучистым потокам.
						//bBlockStefanBolcman = true;
					////}
					//else {
						bBlockStefanBolcman = false;
					//}
					
						

					// пока разница в температурах между итерациями больше одного градуса, то:
					if (bupdateproperties) {
						// В случае мощности зависящей от температуры производит 
						// обновление рабочего значения мощности.
						const bool bavgpower = true;
						update_power_temperature_depend(s, ls, t, t.border_neighbor, gtdps, ltdp, toldtimestep, bavgpower, bmultipowersourse);
						update_temp_properties(t, fglobal, b, lb, matlist); // обновляем свойства твёрдых материалов
						// и обновляем свойства жидких материалов т.к. они связаны через плотность в поправке Рхи-Чоу
						// плотность, динамическую вязкость, коэффициент линейного температурного расширения.
						update_flow_properties(t, fglobal, b, lb, flow_interior, matlist, false);

						

						/*if (i87 == 0) {
							// На старте мы блокируем Стефана Больцмана дав сойтись лучистым потокам.
							// Вычисление осреднённых температур в К на границах вакуумных промежутков:
							for (integer i23 = 0; i23 < lb; i23++) {
								update_avg_temperatures(t.potent, b[i23]);
							}
							// Вычисление плотностей радиационных тепловых потоков:
							for (integer i23 = 0; i23 < lb; i23++) {
								calculation_density_radiation_heat_flux(b[i23]);
							}
						}
						else if (i87 == 1) {
							// ничего не делаем.
						}
						else*/// if (i87 % 2 == 0) {
						
							// На старте мы блокируем Стефана Больцмана дав сойтись лучистым потокам.
							// Вычисление осреднённых температур в К на границах вакуумных промежутков:
							for (integer i23 = 0; i23 < lb; i23++) {
								update_avg_temperatures(t.potent, b[i23]);
							}
							
							// Вычисление плотностей радиационных тепловых потоков:
							for (integer i23 = 0; i23 < lb; i23++) {
								calculation_density_radiation_heat_flux(b[i23]);
							}
							
						//}

							
					}
					res = 0.0;

					//printf("%e\n",s[0].power); // debug
					//getchar();

					// Обновление мощности тепловыделения во всех внутренних узлах.
#pragma omp parallel for 
					for (integer i47 = 0; i47 < t.maxelm; i47++) {
						// Скорость в том что значение не вычисляется как раньше а просто хранится.
						integer ib = t.whot_is_block[i47];
						t.Sc[i47]= get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, t.potent[i47]);
						if (t.Sc[i47]  < 0.0) {
							//printf("ERROR!!! control volume [%lld] is negative t.Sc = %e\n", i47, t.Sc[i47] );
							std::cout << "ERROR!!! control volume [" << i47 << "] is negative t.Sc = " << t.Sc[i47] << std::endl;
							//system("PAUSE");
						}
					}


					
					if (t.ptr != nullptr) {
						// У нас есть жидкие ячейки в которых задано поле скорости.

						if ((0 == err_inicialization_data) || 
							(starting_speed_Vx*starting_speed_Vx +
								starting_speed_Vy * starting_speed_Vy + 
								starting_speed_Vz * starting_speed_Vz > 1.0e-30)) {
							// Диагностика ошибки: конвекцию надо учитывать но она не учитывается.



							//if (fglobal[0].maxelm == 0) {
								//printf("Speed not constr struct.\n");
								//system("PAUSE");
							//}

							for (integer iP = 0; iP < fglobal[0].maxelm; iP++) {

								// вычисляем скорректированный массовый поток через грани КО.
								// Массовый поток вычисляется по обычным формулам но в данном
								// случае без монотонизирующей поправки Рхи-Чоу. При его вычислении используются
								// простая линейная интерполяция скорости на грань КО.



								return_calc_correct_mass_flux_only_interpolation(iP,
									fglobal[0].potent,
									fglobal[0].pa,
									fglobal[0].prop,
									fglobal[0].prop_b,
									fglobal[0].nvtx,
									fglobal[0].neighbors_for_the_internal_node,
									fglobal[0].maxelm,
									fglobal[0].mf[iP], // возвращаемое значение массового потока
									fglobal[0].border_neighbor,
									ls, lw,
									t.ilevel_alice, fglobal[0].ptr);


							}


							// Обязательная проверка корректности выполнения условий прилипания.
							iscorrectmf(fglobal[0].mf, fglobal[0].maxelm, fglobal[0].neighbors_for_the_internal_node, fglobal[0].border_neighbor, ls, lw, w);

						}

					}

					//doublereal tmax = -1.0e30;
					//for (integer i23 = 0; i23 <t.maxelm + t.maxbound; i23++) {

						//if (t.potent[i23] > tmax) {
							//tmax = t.potent[i23];
						//}
					//}
					//printf("solve1 tmax=%e\n", tmax); getchar();

						// решаем номинально линейную СЛАУ:
						// параметр dbeta отвечает за точность аппроксимации граничных условий:
						// dbeta==1.0 первый порядок, dbeta==1.33333333333 второй порядок, dbeta=1.2 третий порядок аппроксимации.
						doublereal** rsumanbstuff = nullptr;
						doublereal resfluent_temp = 0.0;
						solve(TEMP, // идентификатор уравнения (это уравнение теплопередачи).
							res, // невязка
							f,
							fglobal,
							t,
							rhie_chow,
							s, w, b, ls, lw, lb, // объекты (источники, стенки, блоки).
							dbeta,
							flow_interior,
							bconvective,
							false,
							toldtimestep, // поле температур с предыдущего временного слоя
							told,
							speedoldtimestep, // скорость с предыдущего временного слоя
							mfoldtimestep, // конвективный поток через грани КО с предыдущего временного слоя
							tauparam, // размер шага по времени
							btimedep, // стационарный или нестационарный солвер
							0.0, 0.0, 0.0, // ускорение свободного падения
							matlist, // параметры материалов
							inumiter, bprintmessage, RCh, false,
							nullptr, rsumanbstuff, false, false, poweron_multiplier_sequence,
							m, rthdsdt, resfluent_temp,
							lu, my_union, color, dist_max); // номер глобальной итерации

						if ((bglobal_restart_06_10_2018)) {
							// сетка была перестроена глобально.
							iprohodtek = 100000;
							break;
						}
						else {
							
#pragma omp parallel for private(i)
							for (i = 0; i < t.maxelm + t.maxbound; i++) {
								if (t.potent[i] < -272.15) {
									t.potent[i] = -272.15; // Идентифицируем абсолютный ноль.
								}
							}

							bglobal_first_start_radiation = false;
							t.alpha = 1.0;

							bupdateproperties = true;
							tmax = 0.0;

							// Это требуется для отладки
							//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, 0, false, 0);

							const integer ISIZE = t.maxelm + t.maxbound;

							// вычисленние максимума разности температур между итерациями:
							//for (i=0; i<t.maxelm+t.maxbound; i++) tmax=fmax(tmax,fabs(t.potent[i]-told[i])); 
							// 23 декабря 2015
							// На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
							// физичнее не смущать людей и приводить температуру только во внутренних КО.

							tmax=get_max_array_elm(t.potent, told, ISIZE);

							doublereal maxdomain = get_max_array_elm(t.potent, ISIZE); // Наименьшая температура в градусах Цельсия.

							

							if ((1 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) &&
								(AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5 == stabilization_amg1r5_algorithm)) {
									// Нелинейная версия amg1r5 алгоритма.
									printf("incomming\n");
									tmax= get_max_array_elm(t.potent, told_nonlinear, ISIZE);

							}
							deltat = tmax;
							//if (deltat > 13.0) {
								//t.alpha = 0.9;
							//}
							//else if (deltat<0.3) {
								//t.alpha = 1.0;
						//	}
							// В случае нелинейного граничного условия применяем нижнюю релаксацию.
							//if ((0 == err_inicialization_data) || (adiabatic_vs_heat_transfer_coeff > ADIABATIC_WALL_BC) || (breakRUMBAcalc_for_nonlinear_boundary_condition)) {
							if ((adiabatic_vs_heat_transfer_coeff > DEFAULT_CABINET_BOUNDARY_CONDITION::ADIABATIC_WALL_BC) ||
								(breakRUMBAcalc_for_nonlinear_boundary_condition)) {


								// High Order Term Relaxation 23.3.1.10 Fluent Solver Theory
								// default relaxation factor on steady-state cases is 0.25.
								doublereal fHORF = 0.25; // for steady state problem.
								if (btimedep) { // unsteady problems.
									fHORF = 0.75; // ANSYS Fluent Theory Guide.
								}
								if ((((iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 3) ||
									(iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 7))) && 
									(my_amg_manager.istabilizationTemp == 3)) {
									// Румба CAMG нелинейная версия решателя.
									fHORF = 1.0;
								}
								if ((1 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) &&
									(AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5 == stabilization_amg1r5_algorithm)) {
									// специальная нелинейная версия amg1r5 алгоритма.
									fHORF = 1.0;
								}
								//fHORF = 0.00625;
								if (bvacuumPrism) {
									//printf("vacuum determinate 23_11_2016.\n");
									//system("pause");
									//немного понизим усилим нижнюю релаксацию.
									// Ни в коем случае не делать fHORF больше чем 0.02.
									if ((((iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 3) ||
										(iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 7))) && 
										(my_amg_manager.istabilizationTemp != 3)) {
										// Это не специальная нелинейная версия кода Румба CAMG.
										fHORF = 0.01; // 0.01!!!
									}

									if ((1 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) &&
										(AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5 != stabilization_amg1r5_algorithm)) {
										// Это не специальная нелинейная версия кода amg1r5 CAMG.
										fHORF = 0.01; // 0.01!!!
									}									

								}
								//printf("fHORF=%e\n", fHORF); getchar();
								// Глобальное сокращение итераций по устранению нелинейности в системе.
								//fHORF = 1.0;

#pragma omp parallel for private(i)
								for (i = 0; i < t.maxelm + t.maxbound; i++) {
									if (t.potent[i] > -272.15) {
										t.potent[i] = told[i] + fHORF * (t.potent[i] - told[i]);
									}
									else {
										t.potent[i] = -272.15;
									}
								}
							}

							// Проверка на физическую корректность, не может опуститься ниже абсоютного нуля.
#pragma omp parallel for private(i)
							for (i = 0; i < t.maxelm + t.maxbound; i++) {
								if (t.potent[i] < -272.15) {
									// В холодных режимах АППАРАТА прописанных Сидорчуком достижимы математически
									// температуры ниже абсолютного нуля.
									t.potent[i] = -272.15;
									//printf("fatal error: temperature is < -273.15C\n");
								  //getchar();
								  //	exit(1);
								}
							}



#pragma omp parallel for  private (i) schedule (guided)
							for (i = 0; i < t.maxelm + t.maxbound; i++) {
								told[i] = t.potent[i]; // копирование
							}
							// экспорт результата вычисления в программу tecplot360:
							//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior,0);
							if (bprintmessage) {
								//printf("exports in techplot successfully completed. \n");
#if doubleintprecision == 1
								printf("temperature nonlinear solver. Global iteration number %lld.\n", ic);
#else
								printf("temperature nonlinear solver. Global iteration number %d.\n", ic);
#endif

								//printf("temperature difference between iterations %3.2f  oC. %e\n", deltat, maxdomain);
								std::cout << "temperature difference between iterations " << deltat << "  oC. " << maxdomain << std::endl;
								doublereal tmaxloc = get_max_array_elm(t.potent, t.maxelm);

								

								if (bprintmessage) {
									printf("Intermediate maximum temperature in default interior\n");
									//printf("is equal %e  oC.\n", tmaxloc);
									std::cout << "is equal " << tmaxloc << "  oC." << std::endl;
								}
								doublereal tminloc = get_min_array_elm(t.potent, t.maxelm); // возвращает минимальный элемент в массиве.

								if (bprintmessage) {
									printf("Intermediate minimum temperature in default interior\n");
									//printf("is equal %e  oC.\n", tminloc);
									std::cout << "is equal " << tminloc << "  oC." << std::endl;
								}

								if ((((iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 3) ||
										(iswitchsolveramg_vs_BiCGstab_plus_ILU2 == 7))) && 
										(my_amg_manager.istabilizationTemp == 3)) {
										// Это  специальная нелинейная версия кода Румба CAMG.
										// Печатает отдаваемый (снимаемый) во внешнюю среду тепловой поток в Вт,
										// проходящий через выходную границу потока. 28.10.2019
										report_out_boundary(fglobal[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, fglobal[0].OpTemp);
									}

								if ((1 == iswitchsolveramg_vs_BiCGstab_plus_ILU2) &&
										(AMG1R5_OUT_ITERATOR::Non_Linear_amg1r5 == stabilization_amg1r5_algorithm)) {
									// Печатает отдаваемый (снимаемый) во внешнюю среду тепловой поток в Вт,
									// проходящий через выходную границу потока. 28.10.2019
									report_out_boundary(fglobal[0], my_global_temperature_struct, ls, lw, w, b, lb, matlist, fglobal[0].OpTemp);
								}


								//fporogmax = fmax(fabs(tmaxloc), fabs(tminloc));
								fporogmax = fabs(tmaxloc- tminloc);
								integer ic62 = 0;

#pragma omp parallel for private(i) reduction(+: ic62)
								for (i = 0; i < t.maxelm; i++) {
									if (t.potent[i] < t.operatingtemperature) {
										ic62++;
										//printf("anomal control volume %d\n", i);
									}
								}
								if (ic62 > 0) {
									std::cout << "maxelm="<< t.maxelm <<" maxbound= "<< t.maxbound<<" anomal internal temperature control volume=" << ic62 << std::endl;
									//getchar();
								}

								//getchar();
							}
							// Вычисление значения максимальной температуры внутри расчётной области и на её границах:
							tmax = get_max_array_elm(t.potent, t.maxelm + t.maxbound); 
							

							diagnostic_critical_temperature( tmax, fglobal, t, b,lb); // предупреждение при превышении максимально допустимой температуры.
							//getchar(); // debug режим отладки
							ic++;
							//if (ic > 20) break; // не сошлось за 20 пассов
							//break; // не будем делать циклов по ниленейности

							// true при первом запуске должно избавить от бесконечной температуры, нефизичности постановки задачи.
							blocker_Newton_Richman = false;
							//system("pause");
						}
				}

				

				if (bdouble_vacuum_PRISM) {
					//	update_Stefan_Bolcman_condition_double_vacuum_PRISM(w, lw, ls, t.maxbound, t.border_neighbor,
					//	t.potent, t.prop_b, t.pa, t.nvtx, t.maxelm);
					//	getchar();
				}
			
				printf("deltat=%e > 0.05*fporogmax=%e\n", deltat, 0.05*fporogmax);
				//getchar();

		}
		
	
	if (!((bglobal_restart_06_10_2018))) {
		// Вычисление значения максимальной температуры внутри расчётной области и на её границах:
		//for (i=0; i<t.maxelm+t.maxbound; i++) tmax=fmax(tmax,fabs(t.potent[i]));
		// 23 декабря 2015
		// На граничных гранях источников тепла мы имеем нефизично высокую температуру, поэтому
		// физичнее не смущать людей и приводить температуру только во внутренних КО. 
		//for (i = 0; i<t.maxelm; i++) tmax = fmax(tmax, fabs(t.potent[i]));
		tmax = get_max_array_elm(t.potent, t.maxelm); 


		if (bprintmessage) {
			printf("Finally maximum temperature in default interior\n");
			//printf("is equal %3.2f  oC. Power %e, W\n", tmax, power_diss_message_06_10_2018);
			std::cout << "is equal " << tmax << "  oC. Power " 
			<< power_diss_message_06_10_2018 << ", W" << std::endl;
		}


		if (bfreeflag) {

			// пока разница в температурах между итерациями больше одного градуса, то:
			if (bupdateproperties) {
				// В случае мощности зависящей от температуры производит 
				// обновление рабочего значения мощности.
				const bool bavgpower = true;
				update_power_temperature_depend(s, ls, t, t.border_neighbor, gtdps, ltdp, toldtimestep, bavgpower, bmultipowersourse);
				update_temp_properties(t, fglobal, b, lb, matlist); // обновляем свойства твёрдых материалов
																	// и обновляем свойства жидких материалов т.к. они связаны через плотность в поправке Рхи-Чоу
																	// плотность, динамическую вязкость, коэффициент линейного температурного расширения.
				update_flow_properties(t, fglobal, b, lb, flow_interior, matlist, false);

				/*if (i87 == 0) {
				// На старте мы блокируем Стефана Больцмана дав сойтись лучистым потокам.
				// Вычисление осреднённых температур в К на границах вакуумных промежутков:
				for (integer i23 = 0; i23 < lb; i23++) {
				update_avg_temperatures(t.potent, b[i23]);
				}
				// Вычисление плотностей радиационных тепловых потоков:
				for (integer i23 = 0; i23 < lb; i23++) {
				calculation_density_radiation_heat_flux(b[i23]);
				}
				}
				else if (i87 == 1) {
				// ничего не делаем.
				}
				else*/// if (i87 % 2 == 0) {

				// На старте мы блокируем Стефана Больцмана дав сойтись лучистым потокам.
				// Вычисление осреднённых температур в К на границах вакуумных промежутков:
				for (integer i23 = 0; i23 < lb; i23++) {
					update_avg_temperatures(t.potent, b[i23]);
				}
				// Вычисление плотностей радиационных тепловых потоков:
				for (integer i23 = 0; i23 < lb; i23++) {
					calculation_density_radiation_heat_flux(b[i23]);
				}

				//}

				
			}

			// Обновление мощности тепловыделения во всех внутренних узлах.
#pragma omp parallel for 
			for (integer i47 = 0; i47 < t.maxelm; i47++) {
				// Скорость в том что значение не вычисляется как раньше а просто хранится.
				integer ib = t.whot_is_block[i47];
				t.Sc[i47] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, t.potent[i47]);
				if (t.Sc[i47] < 0.0) {
					//printf("ERROR!!! control volume [%lld] is negative power t.Sc = %e\n", i47, t.Sc[i47]);
					std::cout << "ERROR!!! control volume [" << i47
					<< "] is negative power t.Sc = " << t.Sc[i47] << std::endl;
					//system("PAUSE");
				}
			}

			if (t.ptr != nullptr) {
				if ((0 == err_inicialization_data) ||
					(starting_speed_Vx*starting_speed_Vx + 
						starting_speed_Vy * starting_speed_Vy +
						starting_speed_Vz * starting_speed_Vz > 1.0e-30)) {
					// Диагностика ошибки: конвекцию надо учитывать но она не учитывается.



					//if (fglobal[0].maxelm == 0) {
						//printf("Speed not constr struct.\n");
						//system("PAUSE");
					//}

					for (integer iP = 0; iP < fglobal[0].maxelm; iP++) {

						// вычисляем скорректированный массовый поток через грани КО.
						// Массовый поток вычисляется по обычным формулам но в данном
						// случае без монотонизирующей поправки Рхи-Чоу. При его вычислении используются
						// простая линейная интерполяция скорости на грань КО.



						return_calc_correct_mass_flux_only_interpolation(iP,
							fglobal[0].potent,
							fglobal[0].pa,
							fglobal[0].prop,
							fglobal[0].prop_b,
							fglobal[0].nvtx,
							fglobal[0].neighbors_for_the_internal_node,
							fglobal[0].maxelm,
							fglobal[0].mf[iP], // возвращаемое значение массового потока
							fglobal[0].border_neighbor,
							ls, lw,
							t.ilevel_alice, fglobal[0].ptr);


					}



					// Обязательная проверка корректности выполнения условий прилипания.
					iscorrectmf(fglobal[0].mf, fglobal[0].maxelm, fglobal[0].neighbors_for_the_internal_node, fglobal[0].border_neighbor, ls, lw, w);

				}
			}

			// Освобождаем память.
			m.bsignalfreeCRSt = true;
			doublereal** rsumanbstuff = nullptr;
			doublereal resfluent_temp = 0.0;
			// Ещё один вызов решателя не должен повлиять на решение, т.к. оно и так уже получено.
			solve(TEMP, // идентификатор уравнения (это уравнение теплопередачи).
				res, // невязка
				f,
				fglobal,
				t,
				rhie_chow,
				s, w, b, ls, lw, lb, // объекты (источники, стенки, блоки).
				dbeta,
				flow_interior,
				bconvective,
				false,
				toldtimestep, // поле температур с предыдущего временного слоя
				told,
				speedoldtimestep, // скорость с предыдущего временного слоя
				mfoldtimestep, // конвективный поток через грани КО с предыдущего временного слоя
				tauparam, // размер шага по времени
				btimedep, // стационарный или нестационарный солвер
				0.0, 0.0, 0.0, // ускорение свободного падения
				matlist, // параметры материалов
				inumiter, bprintmessage, RCh, false,
				nullptr, rsumanbstuff, false, false, poweron_multiplier_sequence,
				m, rthdsdt, resfluent_temp, lu, my_union, color, dist_max); // номер глобальной итерации


#pragma omp parallel for private(i)
			for (i = 0; i < t.maxelm + t.maxbound; i++) {
				if (t.potent[i] < -273.15) {
					t.potent[i] = -273.15; // Идентифицируем абсолютный ноль.
				}
			}



			tmax = get_max_array_elm(t.potent, t.maxelm); 


			if (bprintmessage) {
				printf("Finally maximum temperature in default interior\n");
				//printf("is equal %3.2f  oC. Pdiss= %e, W\n", tmax, power_diss_message_06_10_2018);
				std::cout << "is equal " << tmax << "  oC. Pdiss= " 
				<< power_diss_message_06_10_2018 << ", W" << std::endl;
			}
			doublereal tmin = get_min_array_elm(t.potent, t.maxelm);


			if (bprintmessage) {
				printf("Finally minimum temperature in default interior\n");
				//printf("is equal %3.2f  oC.\n", tmin);
				std::cout << "is equal " << tmin << "  oC." << std::endl;
			}
			integer ic62 = 0;

#pragma omp parallel for private(i) reduction(+: ic62)
			for (i = 0; i < t.maxelm; i++) {
				if (t.potent[i] < t.operatingtemperature) {
					ic62++;
					//printf("anomal control volume %d\n",i);
				}
			}
			if (ic62 > 0) {
				printf("maxelm=%lld maxbound=%lld anomal internal control volume is negative power=%lld \n", t.maxelm, t.maxbound, ic62);
				//getchar();
			}

		}		
	}
	else {
		printf("Global Grid rebuild 06.10.2018\n.");
	}

	{
		// Внимание. Обязательный финишный update свойств материалов, иначе
		// при экспорте в tecplot материалы будут неверны. Вместо температурно зависимых матриалов
		// из программной библиотеки материалов будет Дюраль Д16Т.
		// 29.11.2018
		// В случае мощности зависящей от температуры производит 
		// обновление рабочего значения мощности.
		const bool bavgpower = true;
		update_power_temperature_depend(s, ls, t, t.border_neighbor, gtdps, ltdp, toldtimestep, bavgpower, bmultipowersourse);
		update_temp_properties(t, fglobal, b, lb, matlist); // обновляем свойства твёрдых материалов
															// и обновляем свойства жидких материалов т.к. они связаны через плотность в поправке Рхи-Чоу
															// плотность, динамическую вязкость, коэффициент линейного температурного расширения.
		update_flow_properties(t, fglobal, b, lb, flow_interior, matlist, false);

	}
	// Освобождение оперативной памяти:
	if (told != nullptr) {
		delete[] told;
		told = nullptr;
	}
	delete[] told_nonlinear;

	if (rthdsdt != nullptr) {
		delete[] rthdsdt;
		rthdsdt = nullptr;
	}
} // solve_nonlinear_temp

// Графовый метод решения уравнений кондукционной теплопередачи.
// begin 24.06.2020
// end **.**.**
#include "network_T_solver.cpp"

  
  // гиперболический косинус
doublereal mycosh(doublereal x) {
	doublereal r=0.0;
	try
	{
	  // наивная реализация, на самом деле нужна проверка на переполнение.
	  r=0.5*(exp(x)+exp(-x));
	}
	catch(integer iex)
	{
#if doubleintprecision == 1
		printf("exeption identifier=%lld in function cosh(x) for Smagorinsky Model/n", iex);
#else
		printf("exeption identifier=%d in function cosh(x) for Smagorinsky Model/n", iex);
#endif
		
		printf("please, press any key to exit is program/n");
		//getchar();
		system("pause");
		exit(0);
	}
	return r;
} // mycosh

// возвращает квадрат длины по Смагоринскому на неравномерной сетке.
doublereal smagorinsky_length2(doublereal Cs, doublereal dx, doublereal dy, doublereal dz, doublereal dist,
	                     doublereal roughness, integer ipowerroughness,
	                     bool bfdelta, bool bSmagorinsky_Lilly, 
						 bool bsurface_roughness) {
	// Cs - константа Смагоринского.
	// В данном коде она может быть как положительной так и отрицательной. 
	// если bfdelta==true то включается корректирующая поправка на неравномерной сетке.
	// если bSmagorinsky_Lilly==true то включается модель Smagorinsky_Lilly
	// Если bsurface_roughness==true значит нужно учитывать шероховатость поверхности,
	// roghness значение шероховатости поверхности в СИ [м], 
	// ipowerroughness - показатель степени в формуле для шероховатости.
	// Все эти опции не конфликтуют и могут использоваться в любом сочетании.

	doublereal dV=dx*dy*dz; // объём контрольного объёма.
	//doublereal delta=pow(dV,1.0/3.0); // кубический корень из объёма контрольного объёма.
	doublereal delta=exp((1.0/3.0)*log(dV)); 
	doublereal length=Cs*delta; // длина перемешивания.
	// Для неравномерной сетки рекомендуется ввести следующий поправочный коэффициент.
	doublereal fdelta=1.0; // по умолчанию корректирующая поправка выключена.
	if (bfdelta) {
	   doublereal rK=2.0/(3.0*sqrt(3.0));
	   doublereal *xi=new doublereal[3];
	   xi[0]=dx; xi[1]=dy; xi[2]=dz;
	   // отсортировать массив по возрастанию.
 	   bool bflag=true;
	   while (bflag) {
		   bflag=false;
		   for (integer j=0; j<2; j++) {
			   if (xi[j]>xi[j+1]) {
				   // swap
				   doublereal xbuf=xi[j];
				   xi[j]=xi[j+1];
				   xi[j+1]=xbuf;
				   bflag=true;
			   }
		   }
	   }
	   doublereal rloga1=log(xi[0]/xi[2]); // натуральный логарифм.
	   doublereal rloga2=log(xi[1]/xi[2]);
	   fdelta=mycosh(rK*sqrt(rloga1*rloga1-rloga1*rloga2+rloga2*rloga2));
	   if (xi != nullptr) {
		   delete[] xi;
	   }
	}
	length*=fdelta;
	if (bsurface_roughness) {
		// Дополнительный учёт шероховатости твёрдой стенки.
		if (ipowerroughness==2) {
			doublereal fr=1.0/sqrt(1.0+(length/(0.419*(dist+roughness)))*(length/(0.419*(dist+roughness))));
			length*=fr; // применение демпфирующей функции.
		}
		if (ipowerroughness==1) {
			length=(length*0.419*(dist+roughness))/(0.419*(dist+roughness)+length);
		}
	}
	if (bSmagorinsky_Lilly) {
		// модель Смагоринского-Лиллу
		// для длины перемешивания.
		length=fmin(length,0.419*dist);
	}
	doublereal length2;
	length2=length*length;
	if (length<0.0) length2*=-1.0;	
	return length2;
}

// интерполяция турбулентной динамической вязкости на границу Расчётной области.
void my_boundary_musgs_LES(integer inumber, integer maxelm, 
							  BOUND* border_neighbor, integer ls, integer lw,
							  WALL* w, TOCHKA* pa, integer** nvtx, doublereal* &potent
							  ) 
{

	 // potent - вектор со значениями динамической вязкости.
	 // На твёрдой неподвижной стенке значение турбулентной динамической вязкости равно нулю, а
	 // на остальных границах значение турбулентной динамической вязкости сносится (копируется) из 
	 // ближайшей внутренней точки на границу расчётной области.

     // inumber - номер граничного КО.
	 // inumber изменяется от 0..maxbound-1

     bool bzeromusgs=false; // нулевое значение турбулентной динамической вязкости. 


     // Сначала запишем граничные условия Дирихле
	 if ((border_neighbor[inumber].MCB<(ls + lw)) && (border_neighbor[inumber].MCB >= ls) && (!w[border_neighbor[inumber].MCB - ls].bsymmetry) && (!w[border_neighbor[inumber].MCB - ls].bpressure) && (!w[border_neighbor[inumber].MCB - ls].bopening)) {
		// Задана скорость на границе
        // Это не граница симметрии и не выходная граница.

		 doublereal vel_mag=sqrt(w[border_neighbor[inumber].MCB-ls].Vx*w[border_neighbor[inumber].MCB-ls].Vx+
			          w[border_neighbor[inumber].MCB-ls].Vy*w[border_neighbor[inumber].MCB-ls].Vy+
					  w[border_neighbor[inumber].MCB-ls].Vz*w[border_neighbor[inumber].MCB-ls].Vz);

		 doublereal epsilon0=1e-32;

		 if (fabs(vel_mag)<epsilon0) {
			 // скорость на стенке равна нулю.
			 // значит это твёрдая неподвижная стенка на 
			 // которой турбулентная динамическая вязкость должна быть равна нулю.
			 bzeromusgs=true;
		 } else bzeromusgs=false;
        
	}
	else if (( (border_neighbor[inumber].MCB==(ls+lw)) ||(border_neighbor[inumber].MCB<ls)) ) { // 
		// источник тоже является твёрдой неподвижной стенкой.
		// граница по умолчанию также является твёрдой неподвижной стенкой.

        // нулевое значение турбулентной вязкости.
		bzeromusgs=true;

	}
	else  {
		// Все остальные границы: выходная, с условием симметрии, входная.
		// На них должно стоять однородное условие Неймана.
        bzeromusgs=false;
	}

	 if (bzeromusgs) {
			 // на твёрдой неподвижной стенке значение турбулентной
		     // динамической вязкости равно нулю.
		     potent[border_neighbor[inumber].iB]=0.0; 
	} else  {
		// однородное условие Неймана.
		// иными словами мы скопируем значение турбулентной вязкости изнутри расчётной области на границу.
		potent[border_neighbor[inumber].iB]=potent[border_neighbor[inumber].iI];
	}
} // my_boundary_musgs_LES

  // коррекция массового расхода после решения уравнений на компоненты скорости или
  // после коррекции скоростей (в данном случае выбрано после коррекции скоростей).
void mass_balance(FLOW &f, integer lw, integer ls, WALL* w)
{
	// 6 августа 2016 Теперь с разделением на несвязанные гидрдинамические подобласти.
	bool bOk = true;
	integer imarker_fluid_color = 0; // последовательно для всех цветов начиная с единицы.

	bool nobcont = true;
	for (integer i_1 = 0; i_1 < lw; i_1++) {
		if ((!w[i_1].bopening) && (sqrt((w[i_1].Vx)*(w[i_1].Vx) + (w[i_1].Vy)*(w[i_1].Vy) + (w[i_1].Vz)*(w[i_1].Vz)) > 1.0e-20)) {
			nobcont = false;
		}
		if (w[i_1].bopening) nobcont = false;
	}

	// Если у нас есть opening то мы делаем коррекцию. 
	// Если же у нас все стенки с условием прилипания то мы никакой коррекции не делаем.
	if (nobcont) bOk = false;

	imarker_fluid_color = 0;

	while (bOk) {
		bOk = false;
		imarker_fluid_color++;


		doublereal Gset = 0.0; // заданный расход.
		doublereal rhosquare = 0.0;
		doublereal mysquare = 0.0; // только площадь выходной границы (или сумма площадей выходных границ).
		doublereal deltavel = 0.0; // аддитивная поправка к скорости на выходной границе.
								   // Будем считать, что если на границе задана нормальная компонента скорости, то
								   // то жидкость может только втекать внутрь расчётной области. Позже это ограничение
								   // можно снять.
		for (integer inumber = 0; inumber < f.maxbound; inumber++) {
			if ((f.border_neighbor[inumber].MCB < (ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) &&
				(!w[f.border_neighbor[inumber].MCB - ls].bsymmetry) && (!w[f.border_neighbor[inumber].MCB - ls].bpressure &&
				(!w[f.border_neighbor[inumber].MCB - ls].bopening))) {
				if ((f.border_neighbor[inumber].iI<0)||(f.border_neighbor[inumber].iI>=f.maxelm)) {
					printf("error 1 iI=%lld maxelm=%lld\n", f.border_neighbor[inumber].iI, f.maxelm);
					printf("fuinction mass balance. file mysolverv0_03.c.\n");
					system("PAUSE");
				}
				if (f.icolor_different_fluid_domain[f.border_neighbor[inumber].iI] == imarker_fluid_color) {
					bOk = true;

					doublereal dS, rho; // площадь грани.
					rho = f.prop_b[RHO][f.border_neighbor[inumber].iB - f.maxelm];
					// сканируем внутреннюю нормаль.
					switch (f.border_neighbor[inumber].Norm) {
					case E_SIDE:  dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						Gset += rho*dS*fabs(w[f.border_neighbor[inumber].MCB - ls].Vx);
						mysquare += dS;
						break;
					case W_SIDE:  dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						Gset += rho*dS*fabs(w[f.border_neighbor[inumber].MCB - ls].Vx);
						mysquare += dS;
						break;
					case N_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						Gset += rho*dS*fabs(w[f.border_neighbor[inumber].MCB - ls].Vy);
						mysquare += dS;
						break;
					case S_SIDE: dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						Gset += rho*dS*fabs(w[f.border_neighbor[inumber].MCB - ls].Vy);
						mysquare += dS;
						break;
					case T_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
						Gset += rho*dS*fabs(w[f.border_neighbor[inumber].MCB - ls].Vz);
						mysquare += dS;
						break;
					case B_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
						Gset += rho*dS*fabs(w[f.border_neighbor[inumber].MCB - ls].Vz);
						mysquare += dS;
						break;
					}
				}

			}
		}
		//printf("Gset=%e\n",Gset);
		//getchar();
		//system("pause");
		// Итак, мы вычислили Gset - заданный расход.
		// Поправка скорости на границе определяется следующей величиной.
		doublereal istinnjiRashod = 0.0; // расход полученный после решения уравнения на скорость.
		for (integer inumber = 0; inumber < f.maxbound; inumber++) {
			if ((f.border_neighbor[inumber].MCB < (ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) &&
				((w[f.border_neighbor[inumber].MCB - ls].bpressure) || w[f.border_neighbor[inumber].MCB - ls].bopening)) {
				
				if ((f.border_neighbor[inumber].iI < 0) || (f.border_neighbor[inumber].iI >= f.maxelm)) {
					printf("error 2 iI=%lld maxelm=%lld\n", f.border_neighbor[inumber].iI, f.maxelm);
					printf("fuinction mass balance. file mysolverv0_03.c.\n");
					system("PAUSE");
				}
				if (f.icolor_different_fluid_domain[f.border_neighbor[inumber].iI] == imarker_fluid_color) {

					doublereal dS, rho; // площадь грани.
					rho = f.prop_b[RHO][f.border_neighbor[inumber].iB - f.maxelm];
					// сканируем внутреннюю нормаль.
					switch (f.border_neighbor[inumber].Norm) {
					case E_SIDE:  dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						istinnjiRashod -= rho*dS*f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case W_SIDE:  dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						istinnjiRashod += rho*dS*f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case N_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						istinnjiRashod -= rho*dS*f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case S_SIDE: dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						istinnjiRashod += rho*dS*f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case T_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
						istinnjiRashod -= rho*dS*f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case B_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
						istinnjiRashod += rho*dS*f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					}
				}
			}
		}

		//printf("istinnjiRashod=%e\n",istinnjiRashod);
		//getchar();
		//system("pause");
		// Коррекция компоненты скорости и массового потока на выходной границе
		// в соотвествии с заданным расходом. См. Гаврилов Андрей опыт CFD. sigma flow.

		//deltavel=(Gset-istinnjiRashod)/rhosquare;
		deltavel = (Gset - istinnjiRashod) / mysquare; // См. Гаврилов Андрей.
													   // Осуществляем корректирующую аддитивную поправку.
													   // Более того мы искусственно поддерживаем однородные условия Неймана.
		for (integer inumber = 0; inumber < f.maxbound; inumber++) {
			if ((f.border_neighbor[inumber].MCB < (ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && ((w[f.border_neighbor[inumber].MCB - ls].bpressure) || (w[f.border_neighbor[inumber].MCB - ls].bopening))) {
				if (f.icolor_different_fluid_domain[f.border_neighbor[inumber].iI] == imarker_fluid_color) {

					doublereal dS, rho; // площадь грани.
					rho = f.prop_b[RHO][f.border_neighbor[inumber].iB - f.maxelm];

					bool bperpendicularzero = false; // перпендикулярная составляющая скорости к выходной границе равна нулю.

													 // сканируем внутреннюю нормаль.
													 // Будем задавать скорректированную скорость только в граничном узле не трогая ближайший внутренний узел.
					switch (f.border_neighbor[inumber].Norm) {
					case E_SIDE:  f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] += deltavel / rho;
						f.potent[VXCOR][f.border_neighbor[inumber].iB] += deltavel / rho;
						//f.potent[VX][f.border_neighbor[inumber].iI]=f.potent[VX][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						f.mf[f.border_neighbor[inumber].iI][W_SIDE] += deltavel*dS;
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VYCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VZCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case W_SIDE:  f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] += deltavel / rho;
						f.potent[VXCOR][f.border_neighbor[inumber].iB] += deltavel / rho;
						//f.potent[VX][f.border_neighbor[inumber].iI]=f.potent[VX][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						f.mf[f.border_neighbor[inumber].iI][E_SIDE] += deltavel*dS;
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VYCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VZCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case N_SIDE:  f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] += deltavel / rho;
						f.potent[VYCOR][f.border_neighbor[inumber].iB] += deltavel / rho;
						//f.potent[VY][f.border_neighbor[inumber].iI]=f.potent[VY][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						f.mf[f.border_neighbor[inumber].iI][S_SIDE] += deltavel*dS;
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VXCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VZCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case S_SIDE: f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] += deltavel / rho;
						f.potent[VYCOR][f.border_neighbor[inumber].iB] += deltavel / rho;
						//f.potent[VY][f.border_neighbor[inumber].iI]=f.potent[VY][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						f.mf[f.border_neighbor[inumber].iI][N_SIDE] += deltavel*dS;
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VXCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VZCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case T_SIDE:  f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] += deltavel / rho;
						f.potent[VZCOR][f.border_neighbor[inumber].iB] += deltavel / rho;
						//f.potent[VZ][f.border_neighbor[inumber].iI]=f.potent[VZ][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани 
						f.mf[f.border_neighbor[inumber].iI][B_SIDE] += deltavel*dS;
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VXCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VYCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case B_SIDE:  f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] += deltavel / rho;
						f.potent[VZCOR][f.border_neighbor[inumber].iB] += deltavel / rho;
						//f.potent[VZ][f.border_neighbor[inumber].iI]=f.potent[VZ][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
						f.mf[f.border_neighbor[inumber].iI][T_SIDE] += deltavel*dS;
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VXCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VYCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					}
				}
			}
		}
	} // Цикл по всем несвязанным гидродинамическим подобластям.
} // mass_balance


// коррекция массового расхода после решения уравнений на компоненты скорости или
// после коррекции скоростей (в данном случае выбрано после коррекции скоростей).
void mass_balance_09_5_2017(FLOW &f, integer lw, integer ls, WALL* w) 
{
	// 6 августа 2016 Теперь с разделением на несвязанные гидрдинамические подобласти.
	bool bOk = true;
	integer imarker_fluid_color = 0; // последовательно для всех цветов начиная с единицы.

	bool nobcont = true;
	for (integer i_1 = 0; i_1 < lw; i_1++) {
		if ((!w[i_1].bopening)&&(sqrt((w[i_1].Vx)*(w[i_1].Vx) + (w[i_1].Vy)*(w[i_1].Vy) + (w[i_1].Vz)*(w[i_1].Vz)) > 1.0e-20) ){
			nobcont = false;
	    }
		if (w[i_1].bopening) nobcont = false;
	}

	// Если у нас есть opening то мы делаем коррекцию. 
	// Если же у нас все стенки с условием прилипания то мы никакой коррекции не делаем.
	if (nobcont) bOk = false;

	imarker_fluid_color = 0;

	while (bOk) {
		bOk = false;
		imarker_fluid_color++;


		doublereal Gset = 0.0; // заданный расход.
		// Количество вытекаемой из расчётной области жидкости заданное пользователем в граничных условиях.
		doublereal rhosquare = 0.0;
		doublereal mysquare = 0.0; // только площадь выходной границы (или сумма площадей выходных границ).
		doublereal deltavel = 0.0; // аддитивная поправка к скорости на выходной границе.
		// Будем считать, что если на границе задана нормальная компонента скорости, то
		// то жидкость может только втекать внутрь расчётной области. Позже это ограничение
		// можно снять.
		for (integer inumber = 0; inumber < f.maxbound; inumber++) {
			if ((f.border_neighbor[inumber].MCB < (ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && (!w[f.border_neighbor[inumber].MCB - ls].bsymmetry) && (!w[f.border_neighbor[inumber].MCB - ls].bpressure && (!w[f.border_neighbor[inumber].MCB - ls].bopening))) {
				if (f.icolor_different_fluid_domain[f.border_neighbor[inumber].iI] == imarker_fluid_color) {
					bOk = true;

					if (w[f.border_neighbor[inumber].MCB - ls].Vx*w[f.border_neighbor[inumber].MCB - ls].Vx + w[f.border_neighbor[inumber].MCB - ls].Vy*w[f.border_neighbor[inumber].MCB - ls].Vy + w[f.border_neighbor[inumber].MCB - ls].Vz*w[f.border_neighbor[inumber].MCB - ls].Vz > 1.0e-20) {
						// Только те входные границы на которых действительно задана ненулевая скорость.


						doublereal dS, rho; // площадь грани.
						rho = f.prop_b[RHO][f.border_neighbor[inumber].iB - f.maxelm];
						// сканируем внутреннюю нормаль.
						switch (f.border_neighbor[inumber].Norm) {
						case E_SIDE:  dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
							dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
							Gset -= rho*dS*(w[f.border_neighbor[inumber].MCB - ls].Vx);
							// Только ненулевая нормальная компонента скорости.
							if (w[f.border_neighbor[inumber].MCB - ls].Vx*w[f.border_neighbor[inumber].MCB - ls].Vx > 1.0e-20) {
								mysquare += dS;
								rhosquare += rho*dS;
							}
							break;
						case W_SIDE:  dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
							dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
							Gset += rho*dS*(w[f.border_neighbor[inumber].MCB - ls].Vx);
							if (w[f.border_neighbor[inumber].MCB - ls].Vx*w[f.border_neighbor[inumber].MCB - ls].Vx > 1.0e-20) {
								mysquare += dS;
								rhosquare += rho*dS;
							}
							break;
						case N_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
							dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
							Gset -= rho*dS*(w[f.border_neighbor[inumber].MCB - ls].Vy);
							if (w[f.border_neighbor[inumber].MCB - ls].Vy*w[f.border_neighbor[inumber].MCB - ls].Vy > 1.0e-20) {
								mysquare += dS;
								rhosquare += rho*dS;
							}
							break;
						case S_SIDE: dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
							dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
							Gset += rho*dS*(w[f.border_neighbor[inumber].MCB - ls].Vy);
							if (w[f.border_neighbor[inumber].MCB - ls].Vy*w[f.border_neighbor[inumber].MCB - ls].Vy > 1.0e-20) {
								mysquare += dS;
								rhosquare += rho*dS;
							}
							break;
						case T_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
							dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
							Gset -= rho*dS*(w[f.border_neighbor[inumber].MCB - ls].Vz);
							if (w[f.border_neighbor[inumber].MCB - ls].Vz*w[f.border_neighbor[inumber].MCB - ls].Vz > 1.0e-20) {
								mysquare += dS;
								rhosquare += rho*dS;
							}
							break;
						case B_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
							dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
							Gset += rho*dS*(w[f.border_neighbor[inumber].MCB - ls].Vz);
							if (w[f.border_neighbor[inumber].MCB - ls].Vz*w[f.border_neighbor[inumber].MCB - ls].Vz > 1.0e-20) {
								mysquare += dS;
								rhosquare += rho*dS;
							}
							break;
						}
					}
				}

			}
		}
		//printf("Gset=%e\n",Gset);
		//getchar();
		//system("pause");
		// Итак, мы вычислили Gset - заданный расход. Это количество вытекаемой из расчётной области жидкости в полном соответствии с тем
		// что задал пользователь в граничных условиях.
		// Поправка скорости на границе определяется следующей величиной.
		doublereal istinnjiRashod = 0.0; // расход полученный после решения уравнения на скорость.
		// istinnjiRashod - количество вытекамой жидкости через pressure и opening границы.
		for (integer inumber = 0; inumber < f.maxbound; inumber++) {
			if ((f.border_neighbor[inumber].MCB < (ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && ((w[f.border_neighbor[inumber].MCB - ls].bpressure) || w[f.border_neighbor[inumber].MCB - ls].bopening)) {
				if (f.icolor_different_fluid_domain[f.border_neighbor[inumber].iI] == imarker_fluid_color) {

					doublereal dS, rho; // площадь грани.
					rho = f.prop_b[RHO][f.border_neighbor[inumber].iB - f.maxelm];
					// сканируем внутреннюю нормаль.
					switch (f.border_neighbor[inumber].Norm) {
					case E_SIDE:  dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						istinnjiRashod -= rho*dS*f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case W_SIDE:  dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						istinnjiRashod += rho*dS*f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case N_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						istinnjiRashod -= rho*dS*f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case S_SIDE: dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						istinnjiRashod += rho*dS*f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case T_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
						istinnjiRashod -= rho*dS*f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					case B_SIDE:  dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
						istinnjiRashod += rho*dS*f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB];
						rhosquare += rho*dS;
						mysquare += dS;
						break;
					}
				}
			}
		}

		//printf("istinnjiRashod=%e\n",istinnjiRashod);
		//getchar();
		//system("pause");
		// Коррекция компоненты скорости и массового потока на выходной границе
		// в соотвествии с заданным расходом. См. Гаврилов Андрей опыт CFD. sigma flow.

		doublereal alpharelax = 1.0;
		doublereal v_zv_zv;

		deltavel=-(Gset-istinnjiRashod)/rhosquare;
		//deltavel = -(Gset - istinnjiRashod) / mysquare; // См. Гаврилов Андрей.
		// deltavel - это плотность умноженая на скорость равную среднему значению лишней вытеченной из расчётной области скорости.
		// В идеале эта величина должна быть равна нулю.
		// лишняя вытеченная: где вытекает мы уменьшаем на эту величину (минус deltavel),
		// где втекает то мы увеличиваем втечение на эту величину (плюс deltavel)
		// Надо уменьшить на границах вытечения значение скорости вытекания на эту величину.
		// Осуществляем корректирующую аддитивную поправку.
		// Более того мы искусственно поддерживаем однородные условия Неймана.
		for (integer inumber = 0; inumber < f.maxbound; inumber++) {
			if ((f.border_neighbor[inumber].MCB < (ls + lw)) && (f.border_neighbor[inumber].MCB >= ls) && ((w[f.border_neighbor[inumber].MCB - ls].bpressure) || (w[f.border_neighbor[inumber].MCB - ls].bopening))) {
				if (f.icolor_different_fluid_domain[f.border_neighbor[inumber].iI] == imarker_fluid_color) {

					doublereal dS, rho; // площадь грани.
					rho = f.prop_b[RHO][f.border_neighbor[inumber].iB - f.maxelm];

					bool bperpendicularzero = false; // перпендикулярная составляющая скорости к выходной границе равна нулю.

					// сканируем внутреннюю нормаль.
					// Будем задавать скорректированную скорость только в граничном узле не трогая ближайший внутренний узел.
					switch (f.border_neighbor[inumber].Norm) {
					case E_SIDE: 
						v_zv_zv = f.potent[VXCOR][f.border_neighbor[inumber].iB] + deltavel;
						//f.potent[VX][f.border_neighbor[inumber].iB] += deltavel;
						//f.potent[VXCOR][f.border_neighbor[inumber].iB] += deltavel;
						f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv) + (1 - alpharelax)*f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB];
						f.potent[VXCOR][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VXCOR][f.border_neighbor[inumber].iB];
						//f.potent[VX][f.border_neighbor[inumber].iI]=f.potent[VX][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						//f.mf[f.border_neighbor[inumber].iI][WSIDE] += rho*deltavel*dS;
						f.mf[f.border_neighbor[inumber].iI][W_SIDE] = alpharelax*(v_zv_zv*rho*dS) + (1 - alpharelax)*f.mf[f.border_neighbor[inumber].iI][W_SIDE];
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VYCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VZCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case W_SIDE:
						v_zv_zv = f.potent[VXCOR][f.border_neighbor[inumber].iB] - deltavel;
						//f.potent[VX][f.border_neighbor[inumber].iB] -= deltavel;
						//f.potent[VXCOR][f.border_neighbor[inumber].iB] -= deltavel;
						f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB];
						f.potent[VXCOR][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VXCOR][f.border_neighbor[inumber].iB];
						//f.potent[VX][f.border_neighbor[inumber].iI]=f.potent[VX][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].y;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						//f.mf[f.border_neighbor[inumber].iI][ESIDE] -= rho*deltavel*dS;
						f.mf[f.border_neighbor[inumber].iI][E_SIDE] = alpharelax*(v_zv_zv*rho*dS) + (1 - alpharelax)*f.mf[f.border_neighbor[inumber].iI][E_SIDE];
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VYCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VZCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case N_SIDE:  
						v_zv_zv = f.potent[VYCOR][f.border_neighbor[inumber].iB] + deltavel;
						//f.potent[VY][f.border_neighbor[inumber].iB] += deltavel;
						//f.potent[VYCOR][f.border_neighbor[inumber].iB] += deltavel;
						f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB];
						f.potent[VYCOR][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VYCOR][f.border_neighbor[inumber].iB];
						//f.potent[VY][f.border_neighbor[inumber].iI]=f.potent[VY][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						//f.mf[f.border_neighbor[inumber].iI][SSIDE] += rho*deltavel*dS;
						f.mf[f.border_neighbor[inumber].iI][S_SIDE] = alpharelax*(v_zv_zv*rho*dS) + (1 - alpharelax)*f.mf[f.border_neighbor[inumber].iI][S_SIDE];
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VXCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VZCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case S_SIDE: 
						v_zv_zv = f.potent[VYCOR][f.border_neighbor[inumber].iB] - deltavel;
						// Вытекает и мы уменьшаем на лишнюю вытеченную.
						//f.potent[VY][f.border_neighbor[inumber].iB] -= deltavel;
						//f.potent[VYCOR][f.border_neighbor[inumber].iB] -= deltavel;
						f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB];
						f.potent[VYCOR][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VYCOR][f.border_neighbor[inumber].iB];
						//f.potent[VY][f.border_neighbor[inumber].iI]=f.potent[VY][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[4][f.border_neighbor[inumber].iI] - 1].z - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].z); // площадь грани
						//f.mf[f.border_neighbor[inumber].iI][NSIDE] -= rho*deltavel*dS;
						f.mf[f.border_neighbor[inumber].iI][N_SIDE] = alpharelax*(v_zv_zv*rho*dS) + (1 - alpharelax)*f.mf[f.border_neighbor[inumber].iI][N_SIDE];
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VXCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VZCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case T_SIDE: 
						v_zv_zv = f.potent[VZCOR][f.border_neighbor[inumber].iB] + deltavel;
						//f.potent[VZ][f.border_neighbor[inumber].iB] += deltavel;
						//f.potent[VZCOR][f.border_neighbor[inumber].iB] += deltavel;
						f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB];
						f.potent[VZCOR][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VZCOR][f.border_neighbor[inumber].iB];
						//f.potent[VZ][f.border_neighbor[inumber].iI]=f.potent[VZ][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани 
						//f.mf[f.border_neighbor[inumber].iI][BSIDE] += rho*deltavel*dS;
						f.mf[f.border_neighbor[inumber].iI][B_SIDE] = alpharelax*(v_zv_zv*rho*dS) + (1 - alpharelax)*f.mf[f.border_neighbor[inumber].iI][B_SIDE];
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VXCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VYCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					case B_SIDE: 
						v_zv_zv = f.potent[VZCOR][f.border_neighbor[inumber].iB] - deltavel;
						// Вытекает и мы уменьшаем на лишнюю вытеченную.
						//f.potent[VZ][f.border_neighbor[inumber].iB] -= deltavel;
						//f.potent[VZCOR][f.border_neighbor[inumber].iB] -= deltavel;
						f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VELOCITY_Z_COMPONENT][f.border_neighbor[inumber].iB];
						f.potent[VZCOR][f.border_neighbor[inumber].iB] = alpharelax*(v_zv_zv)+(1 - alpharelax)*f.potent[VZCOR][f.border_neighbor[inumber].iB];
						//f.potent[VZ][f.border_neighbor[inumber].iI]=f.potent[VZ][f.border_neighbor[inumber].iB];
						dS = f.pa[f.nvtx[1][f.border_neighbor[inumber].iI] - 1].x - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].x;
						dS *= (f.pa[f.nvtx[2][f.border_neighbor[inumber].iI] - 1].y - f.pa[f.nvtx[0][f.border_neighbor[inumber].iI] - 1].y); // площадь грани
						//f.mf[f.border_neighbor[inumber].iI][TSIDE] -= rho*deltavel*dS;
						f.mf[f.border_neighbor[inumber].iI][T_SIDE] = alpharelax*(v_zv_zv*rho*dS) + (1 - alpharelax)*f.mf[f.border_neighbor[inumber].iI][T_SIDE];
						// фиксируем перпендикулярную составляющую скорости равной нулю.
						if (bperpendicularzero) {
							f.potent[VELOCITY_X_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VXCOR][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VELOCITY_Y_COMPONENT][f.border_neighbor[inumber].iB] = 0.0;
							f.potent[VYCOR][f.border_neighbor[inumber].iB] = 0.0;
						}
						break;
					}
				}
			}
		}
	} // Цикл по всем несвязанным гидродинамическим подобластям.
} // mass_balance


// Полунеявная процедура для связывающих давление уравнений С.Патанкар, Б.Сполдинг 1972год.
void my_version_SIMPLE_Algorithm3D(doublereal &continity, integer inumiter, FLOW &f, FLOW* &fglobal, TEMPER &t, doublereal** &rhie_chow,
								   BLOCK* &b, integer lb, SOURCE* &s, integer ls, WALL* &w, integer lw, doublereal dbeta,
								   integer flow_interior, integer iflow, bool bfirst_start, doublereal dgx, doublereal dgy, doublereal dgz,
								   TPROP* &matlist, bool btimedep, doublereal dtimestep, doublereal dtimestepold, doublereal phisicaltime,
								   doublereal* toldtimestep, doublereal** speedoldtimestep, doublereal** mfoldtimestep, 
								   bool bprintmessage, TEMP_DEP_POWER* gtdps, integer ltdp, FLUENT_RESIDUAL &rfluentres,
								   doublereal &rfluentrestemp, doublereal* &smagconstolditer, doublereal** &mfold, integer itempersolve, 
								   QuickMemVorst& m, bool bextendedprint, doublereal** &SpeedCorOld, doublereal* &xb,
								   doublereal* &rthdsd, doublereal* &rthdsdt, integer lu, UNION* &my_union,
	                               integer* &color, integer dist_max_fluid, integer* &color_solid, integer dist_max_solid) 
{


	

	// Массовый расход на выходной границе удовлетворяющий уравнению неразрывности.
	mass_balance(f, lw, ls, w); // Это очень важно. См. Гаврилов Андрей.
	// отключил 6.05.2017
	
	// rfluentrestemp невязка в стиле Fluent для температуры.
	// если itempersolve == 0 то уравнение теплопроводности решать ненужно, 
	// а если itempersolve== 1 то требуется решить также и уравнение теплопроводности.

	bdontstartsolver = false;

	doublereal RCh=1.0;
	// RCh = 0.2;//08.03.2019 Похоже способна существенно улучшить начальные всплески на 1-6 итерации SIMPLE алгоритма.
	// Если bVERYStable==true то мы стремимся к стабильности вычислительного процесса:
	// 1. используем на границе аппроксимацию первого порядка dbeta==1.0;
	// 2. используем для аппроксимации конвективного члена схему первого порядка (противопоточную).
	// Будем считать что первые 100 итераций нужно применять схему bVERYStable==true. За 100 итераций решение 
	// "устаканится" и можно будет применить схему более высокого порядка точности за счёт метода отложенной коррекции.
	bool bVERYStable=true;
	if (inumiter>2000) bVERYStable=false; // переходим к аппроксимации высокого порядка.

	if (0) {
		xyplot( fglobal, flow_interior, t);
		printf("my version SIMPLE Algorith. OK.\n");
	    //getchar(); // debug avtosave
		system("pause");
	}

	// Если bfirst_start==true значит мы имеем дело с первой итерацией алгоритма SIMPLE.
	// smagconstolditer - константа Смагоринского с предыдущей итерации. (используется для нижней релаксации).
	const doublereal smagconstURF=0.001; // рекомендуемое в литературе значение.

	// inumiter - номер итерации.
	if (bprintmessage) {
		std::cout << "inumber iter SIMPLEC=" << inumiter << std::endl;
	}

	//printf("turb\n");
	//getchar();


	// начальное поле давления предполагается заданным.
	
	doublereal res=0.0;
	// находим промежуточное поле скоростей:
	// Граничные узлы не подвергаются нижней релаксации !!!.

	// UPDATE PROPERTIES 
	// обновляет свойства материалов зависящих от температуры и давления:
	// плотность, динамическую вязкость, коэффициент линейного температурного расширения.
    update_flow_properties(t, fglobal, b, lb, flow_interior, matlist, bfirst_start);

	// В случае Zero Equation Turbulence Model вычисление
	// коэффициента динамической турбулентной вязкости.
	if (f.iflowregime== VISCOSITY_MODEL::ZEROEQMOD) {
#pragma omp parallel for
		for (integer i=0; i<(f.maxelm+f.maxbound); i++) {

			doublereal rho=0.0; // Вычисление плотности.
			if (i<f.maxelm) rho=f.prop[RHO][i];
			else rho=f.prop_b[RHO][i-f.maxelm];
			doublereal PrandtlLength=fmin(0.419*f.rdistWall[i],0.09*f.rdistWallmax); // Формула Эскудиера (1966).
			//f.potent[MUT][i]=rho*PrandtlLength*PrandtlLength*f.SInvariantStrainRateTensor[i];
			f.potent[MUT][i] = rho*PrandtlLength*PrandtlLength*f.potent[CURL][i]/sqrt(2.0);
			
			if (0&&(inumiter == 40)) {
				printf("rho=%e dw=%e dwmax=%e plen=%e sigma=%e\n", rho, f.rdistWall[i], f.rdistWallmax, PrandtlLength, f.SInvariantStrainRateTensor[i]);
				system("pause");
			}
	    }
		// Принудительно сделаем чтобы граничные значения совпадали с ближайшими внутренними.
		// Это нужно прежде всего для корректности задания граничных условий.
		for (integer iP=0; iP<f.maxelm; iP++) {
            integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
			iE = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE1; iN = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE1; iT = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE1;
			iW = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE1; iS = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE1; iB = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE1;
			if (iE > -1) {
				if (iE >= f.maxelm) {
					f.potent[MUT][iE] = f.potent[MUT][iP];
				}
			}
			if (iW > -1) {
				if (iW >= f.maxelm) {
					f.potent[MUT][iW] = f.potent[MUT][iP];
				}
			}
			if (iN > -1) {
				if (iN >= f.maxelm) {
					f.potent[MUT][iN] = f.potent[MUT][iP];
				}
			}
			if (iS > -1) {
				if (iS >= f.maxelm) {
					f.potent[MUT][iS] = f.potent[MUT][iP];
				}
			}
			if (iT > -1) {
				if (iT >= f.maxelm) {
					f.potent[MUT][iT] = f.potent[MUT][iP];
				}
			}
			if (iB > -1) {
				if (iB >= f.maxelm) {
					f.potent[MUT][iB] = f.potent[MUT][iP];
				}
			}

			integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
			iE2 = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE2; iN2 = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE2; iT2 = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE2;
			iW2 = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE2; iS2 = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE2; iB2 = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE2;
			if (iE2 > -1) {
				if (iE2 >= f.maxelm) {
					f.potent[MUT][iE2] = f.potent[MUT][iP];
				}
			}
			if (iW2 > -1) {
				if (iW2 >= f.maxelm) {
					f.potent[MUT][iW2] = f.potent[MUT][iP];
				}
			}
			if (iN2 > -1) {
				if (iN2 >= f.maxelm) {
					f.potent[MUT][iN2] = f.potent[MUT][iP];
				}
			}
			if (iS2 > -1) {
				if (iS2 >= f.maxelm) {
					f.potent[MUT][iS2] = f.potent[MUT][iP];
				}
			}
			if (iT2 > -1) {
				if (iT2 >= f.maxelm) {
					f.potent[MUT][iT2] = f.potent[MUT][iP];
				}
			}
			if (iB2 > -1) {
				if (iB2 >= f.maxelm) {
					f.potent[MUT][iB2] = f.potent[MUT][iP];
				}
			}

			integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
			iE3 = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE3; iN3 = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE3; iT3 = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE3;
			iW3 = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE3; iS3 = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE3; iB3 = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE3;
			if (iE3 > -1) {
				if (iE3 >= f.maxelm) {
					f.potent[MUT][iE3] = f.potent[MUT][iP];
				}
			}
			if (iW3 > -1) {
				if (iW3 >= f.maxelm) {
					f.potent[MUT][iW3] = f.potent[MUT][iP];
				}
			}
			if (iN3 > -1) {
				if (iN3 >= f.maxelm) {
					f.potent[MUT][iN3] = f.potent[MUT][iP];
				}
			}
			if (iS3 > -1) {
				if (iS3 >= f.maxelm) {
					f.potent[MUT][iS3] = f.potent[MUT][iP];
				}
			}
			if (iT3 > -1) {
				if (iT3 >= f.maxelm) {
					f.potent[MUT][iT3] = f.potent[MUT][iP];
				}
			}
			if (iB3 > -1) {
				if (iB3 >= f.maxelm) {
					f.potent[MUT][iB3] = f.potent[MUT][iP];
				}
			}

			integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов
			iE4 = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE4; iN4 = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE4; iT4 = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE4;
			iW4 = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE4; iS4 = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE4; iB4 = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE4;
			if (iE4 > -1) {
				if (iE4 >= f.maxelm) {
					f.potent[MUT][iE4] = f.potent[MUT][iP];
				}
			}
			if (iW4 > -1) {
				if (iW4 >= f.maxelm) {
					f.potent[MUT][iW4] = f.potent[MUT][iP];
				}
			}
			if (iN4 > -1) {
				if (iN4 >= f.maxelm) {
					f.potent[MUT][iN4] = f.potent[MUT][iP];
				}
			}
			if (iS4 > -1) {
				if (iS4 >= f.maxelm) {
					f.potent[MUT][iS4] = f.potent[MUT][iP];
				}
			}
			if (iT4 > -1) {
				if (iT4 >= f.maxelm) {
					f.potent[MUT][iT4] = f.potent[MUT][iP];
				}
			}
			if (iB4 > -1) {
				if (iB4 >= f.maxelm) {
					f.potent[MUT][iB4] = f.potent[MUT][iP];
				}
			}
		}
	} // Zero Equation Model

	if (f.iflowregime== VISCOSITY_MODEL::SMAGORINSKY) {
		// Одна из разновидностей модели Смагоринского (LES моделирование).

		// Вычисление для Selective Smagorinsky
		bool* bfibeta=nullptr; // выделение памяти и инициализация происходят внутри функции calc_selective_smagorinsky.

		if (f.smaginfo.bSelectiveSmagorinsky) {
			// Используется модель Selective Smagorinsky
			calc_selective_smagorinsky(f, bfibeta, f.smaginfo.itypeFILTRSelectiveSmagorinsky, f.smaginfo.SSangle);
		}

		

		// динамическая модель турбулентности
		// Германо предложенная в 1991 году.
		// Dynamic Subgrid Scale Model Germano [1991].

		// квадрат постоянной Смагоринского, выделение памяти происходит внутри функции my_Germano_model.
		doublereal* Cs2=nullptr; // размер памяти 0..maxelm-1.

		if (f.smaginfo.bDynamic_Stress) {
			my_Germano_model(f,Cs2,f.smaginfo.itypeFiltrGermano); // определение квадрата константы Смагоринского.
			for (integer i=0; i<f.maxelm; i++) {
				doublereal myCs=f.smaginfo.Cs; // постоянная Смагоринского.
				// память под Cs2 выделена точлько от 0 до maxelm-1.
				myCs=sqrt(fabs(Cs2[i]));
				if (Cs2[i]<0.0) myCs*=-1.0;
				// произведём рекомендуемую в литературе нижнюю релаксацию.
				// Если так не делать то константа Смагоринского может принимать абсолютные значения +,- тысяча, а может быть и больше.
				// Вообще говоря релаксация не спасает от больших значений, она лишь замедляет скорость увеличения величины.
				myCs=(1.0-smagconstURF)*smagconstolditer[i]+smagconstURF*myCs;
				f.potent[FBUF][i]=myCs; // заносим вычисленное распределение константы Смагоринского в вектор для анализа в программе tecplot360.
		    }
		}
		

		for (integer i=0; i<f.maxelm; i++) {
			// Внутренние Контрольные объёмы.
            
			doublereal rho=0.0; // Вычисление плотности.
			rho=f.prop[RHO][i]; // значение плотности во внутреннем контрольном объёме.
			doublereal dx=0.0, dy=0.0, dz=0.0; // размеры текущего контрольного объёма.
			volume3D(i, f.nvtx, f.pa, dx, dy, dz);
			doublereal length2=0.0; // квадрат длины пути смешения
			doublereal myCs=f.smaginfo.Cs; // постоянная Смагоринского.
			if (f.smaginfo.bDynamic_Stress) {
				myCs=sqrt(fabs(Cs2[i]));
				if (Cs2[i]<0.0) myCs*=-1.0;
			}
			length2=smagorinsky_length2(myCs, dx, dy, dz, f.rdistWall[i],
				                        f.smaginfo.roughness, f.smaginfo.ipowerroughness,
										f.smaginfo.bfdelta, f.smaginfo.bSmagorinsky_Lilly, 
										f.smaginfo.bsurface_roughness);

			doublereal fRichardson=1.0; // поправка для течений с кривизной линий тока.
			if (f.smaginfo.bRichardsonCorrect) {
				doublereal SInv=f.SInvariantStrainRateTensor[i];
				if (fabs(SInv)<1e-30) {
					// течение скорее всего отсутствует, поэтому и поправка ненужна.
					fRichardson=1.0;
				}
				else {
					// Турбулентное число Ричардсона.
					doublereal Risgs=(f.potent[CURL][i]/SInv)*(f.potent[CURL][i]/SInv)+f.potent[CURL][i]/SInv;
					if ((1.0-f.smaginfo.rRichardsonMultiplyer*Risgs)>0.0) {
                       fRichardson=sqrt(1.0-f.smaginfo.rRichardsonMultiplyer*Risgs);
					}
					else {
						// Скорее всего инвариант тензора скоростей деформаций слишком мал,
						// что может означать отсутствие течения и значит поправка связанная
						// с крутизной линий тока не нужна => оставим значение 1.0.
						fRichardson=1.0;
					}
				}
			}

			f.potent[MUT][i]=rho*fRichardson*length2*f.SInvariantStrainRateTensor[i];

			if (f.smaginfo.bSelectiveSmagorinsky) {
				// Используется модель Selective Smagorinsky
				if (bfibeta[i]) {
				    // угол между вихрем и осреднённым вихрем достаточно
				    // велик поэтому требуется сгладить (размазать) градиенты включив
				    // дополнительную диссипацию в виде турбулентной вязкости.
					f.potent[MUT][i]=f.potent[MUT][i]*1.0;
			    }
			    else {
				    // угол между вихрем и осреднённым вихрем
				    // достаточно мал, меньше 15 градусов, 
				    // поэтому нужно выключить сглаживающее воздействие
				    // турбулентной вязкости.
				    f.potent[MUT][i]=0.0;
			    }
			}
		}

		if (f.smaginfo.bDynamic_Stress) {
			// Освобождение оперативной памяти.
		    delete[] Cs2;
		}

		if (f.smaginfo.bSelectiveSmagorinsky) {
			delete[] bfibeta; // Освобождение оперативной памяти.
		}

		// Значение турбулентной вязкости на твёрдой неподвижной стенке равно нулю.
		// Для всех остальных границ расчётной области значение турбулентной вязкости 
		// можно получить путём интерполяции (снесения) из ближайшего к границе выпуклого
		// контрольного объёма.
		for (integer i=0; i<f.maxbound; i++) {
			my_boundary_musgs_LES(i, f.maxelm, f.border_neighbor, ls, lw,
							  w, f.pa, f.nvtx, f.potent[MUT]);
		}
	}

	if (f.iflowregime== VISCOSITY_MODEL::RNG_LES) {
		// модель турбулентности базирующаяся на Renormalization Group Theory.
		// в рамках LES подхода. Данная реализация придерживается классического
		// описания на CFD-Wiki.

		for (integer i=0; i<f.maxelm; i++) {
			// Внутренние Контрольные объёмы.
            doublereal rho=0.0; // Вычисление плотности.
			rho=f.prop[RHO][i]; // значение плотности во внутреннем контрольном объёме.
			doublereal dx=0.0, dy=0.0, dz=0.0; // размеры текущего контрольного объёма.
			volume3D(i, f.nvtx, f.pa, dx, dy, dz);
			doublereal dV=dx*dy*dz; // объём контрольного объёма.
	        //doublereal delta=pow(dV,1.0/3.0); // кубический корень из объёма контрольного объёма.
	        doublereal delta=exp((1.0/3.0)*log(dV)); 
			doublereal Crng=0.157; // константа RNG_LES модели.
            doublereal length2=(Crng*delta)*(Crng*delta); // длина перемешивания.
			doublereal mu_sgs=rho*length2*fabs(f.SInvariantStrainRateTensor[i]);
			doublereal mu_lam=f.prop[MU_DYNAMIC_VISCOSITY][i]; 
			doublereal mu_eff=my_dixtomiq_RNG_LES(mu_lam, mu_sgs);
			f.potent[MUT][i]=fmax(mu_eff-mu_lam,0.0); // турбулентная динамическая вязкость.
		}
		// Значение турбулентной вязкости на твёрдой неподвижной стенке равно нулю.
		// Для всех остальных границ расчётной области значение турбулентной вязкости 
		// можно получить путём интерполяции (снесения) из ближайшего к границе выпуклого
		// контрольного объёма.
		for (integer i=0; i<f.maxbound; i++) {
			my_boundary_musgs_LES(i, f.maxelm, f.border_neighbor, ls, lw,
							  w, f.pa, f.nvtx, f.potent[MUT]);
		}
	}

	/*
	// Было решено не реализовывать отдельно модель Германо, а реализовать её в качестве опции
	// к модели Смагоринского и использовать её лишь для определения квадрата константы Смагоринского.
	if (f.iflowregime==GERMANO) {
		// динамическая модель турбулентности
		// Германо предложенная в 1991 году.
		// Dynamic Subgrid Scale Model Germamno [1991].

		// квадрат постоянной Смагоринского, выделение памяти происходит внутри функции my_Germano_model.
		doublereal* Cs2; // размер памяти 0..maxelm-1.
		my_Germano_model(f,Cs2,SIMPSON_FILTR); // определение евадрата константы Смагоринского.

		for (i=0; i<f.maxelm; i++) {
		    // Внутренние Контрольные объёмы.
            doublereal rho=0.0; // Вычисление плотности.
			rho=f.prop[RHO][i]; // значение плотности во внутреннем контрольном объёме.
			doublereal dx=0.0, dy=0.0, dz=0.0; // размеры текущего контрольного объёма.
			volume3D(i, f.nvtx, f.pa, dx, dy, dz);
			doublereal dV=dx*dy*dz; // объём контрольного объёма.
	        //doublereal delta=pow(dV,1.0/3.0); // кубический корень из объёма контрольного объёма.
	        doublereal delta2=exp((2.0/3.0)*log(dV)); 
			doublereal length2=Cs2[i]*delta2; // квадрат пути смешения Прандтля.
			f.potent[MUT][i]=rho*length2*f.SInvariantStrainRateTensor[i]; // турбулентная динамическая вязкость.
		}

		// Освобождение оперативной памяти.
		delete Cs2;

		// Значение турбулентной вязкости на твёрдой неподвижной стенке равно нулю.
		// Для всех остальных границ расчётной области значение турбулентной вязкости 
		// можно получить путём интерполяции (снесения) из ближайшего к границе выпуклого
		// контрольного объёма.
		for (i=0; i<f.maxbound; i++) {
			my_boundary_musgs_LES(i, f.maxelm, f.border_neighbor, ls, lw,
							  w, f.pa, f.nvtx, f.potent[MUT]);
		}

	}
	*/
	if (0&&(inumiter>=66)) {
		// debug
		integer i=0;
		while (i<f.maxelm+f.maxbound) {
			for (integer ic=0; (ic<10)&&(i+2*ic<f.maxelm+f.maxbound); ic++) {
#if doubleintprecision == 1
				printf("%lld %e %lld %e\n", i + ic, f.SInvariantStrainRateTensor[i + ic], i + 2 * ic, f.SInvariantStrainRateTensor[i + 2 * ic]);
#else
				printf("%d %e %d %e\n", i + ic, f.SInvariantStrainRateTensor[i + ic], i + 2 * ic, f.SInvariantStrainRateTensor[i + 2 * ic]);
#endif
				
			}
			i=i+20;
			//getchar();
			system("pause");
		}
	}
	

	

	// экспорт результата вычисления в программу tecplot360:
	if (0&&inumiter==60) {
		 exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior, iflow,false,0, b, lb);
	     printf("update flow properties. OK.\n");
	    // getchar(); // debug avtosave
		 system("pause");
	}

	
	doublereal bonbeta=dbeta;
	if (!bVERYStable) {
		// коээфициент обеспечивающий нужный порядок аппроксимации на границе.
		 bonbeta=1.0; // первый порядок
		//bonbeta=1.33333; // второй полрядок аппроксимации.
		// bonbeta=1.2; // третий порядок
	}

	// Псевдовремя определяется только по конвективно-диффузионной составляющей и 
	// обязательно до процесса симметризации СЛАУ. Поэтому запомним неполностью собранные диагональные элементы
	// соттветствующие лишь конвективно диффузионной составляющей.
	// Только для внутренних контрольных объёмов.
	doublereal** sumanb=new doublereal*[3];
	for (integer i=0; i<3; i++) sumanb[i]=new doublereal[f.maxelm+f.maxbound];
	// инициализация.
#pragma omp parallel for
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		sumanb[VELOCITY_X_COMPONENT][i]=0.0;
		sumanb[VELOCITY_Y_COMPONENT][i]=0.0;
		sumanb[VELOCITY_Z_COMPONENT][i]=0.0;
	}


	// Раньше здесь был цикл устранения нелинейности, НО	
	// дело в том что никакого цикла устранения нелинейности быть недолжно.
	// компоненты скорости фигурирующие при составлении дискретного аналога
	// для обобщённой переменной F на каждой итерации постоянны и равны VELCOR.
	// Если сделать коэффициенты зависящими от текущей скорости, которая меняется то
	// здесь по-видимому возможна расходимость вычислительного процесса.
	if (bprintmessage) {
		printf("VX \n");
	}

	// Если это последняя итерация SIMPLE алгоритма и нам нужно освободить память
	// в BICGSTAB_internal3 то мы временно снимем флаг освобождения памяти а установим его перед PAM.
	bool bflag_free_memory_cfd=m.bsignalfreeCRScfd;
	if (m.bsignalfreeCRScfd) {
        m.bsignalfreeCRScfd=false;
	}
#if doubleintprecision == 1
	//printf("inumiter=%lld\n",inumiter);
#else
	//printf("inumiter=%d\n",inumiter);
#endif
	
	//getchar();

	//printf("turb VX\n");
	//getchar();
	if (inumiter > 1) {
		bdontstartsolver = false;
	}
	else {
		bdontstartsolver = true;
	}

	bool bHORF_speed_on = true; // должно быть именно true.
	// если false то отсутствует сходимость на задачах с включёнными opening границами.
	// По крайней мере сходимость сильно деградирует при отключении этой опции.

		solve(VELOCITY_X_COMPONENT, res, f, fglobal, t, rhie_chow,
			s, w, b, ls, lw, lb, bonbeta,
			flow_interior, false,
			bfirst_start, toldtimestep, nullptr,
			speedoldtimestep, mfoldtimestep,
			dtimestep, btimedep, dgx, dgy, dgz,
			matlist, inumiter, bprintmessage,
			RCh, bVERYStable, nullptr, sumanb, false, false, 1.0, m,
			rthdsd, rfluentres.res_vx, lu, my_union, color, dist_max_fluid);

		// именно здесь верно. 04.05.2017
		rfluentres.res_vx = fluent_residual_for_x(f.slau[VELOCITY_X_COMPONENT], f.slau_bon[VELOCITY_X_COMPONENT], f.potent[VELOCITY_X_COMPONENT], f.maxelm, f.maxbound, VELOCITY_X_COMPONENT); // невязка по формуле fluent.

#pragma omp parallel for
		for (integer i = 0; i < f.maxbound; i++) {
			sumanb[VELOCITY_X_COMPONENT][f.maxelm + i] = f.slau_bon[VELOCITY_X_COMPONENT][i].aw;
		}

		if (bHORF_speed_on&&bHORF) {
#pragma omp parallel for
			for (integer i = 0; i < f.maxelm + f.maxbound; i++) {
				doublereal fHORF = 0.25;
				if (btimedep) {
					fHORF = 0.75;
				}
				f.potent[VELOCITY_X_COMPONENT][i] = f.potent[VXCOR][i] + fHORF*(f.potent[VELOCITY_X_COMPONENT][i] - f.potent[VXCOR][i]);
			}
		}
	

		if (0) {
			//for (integer i = 0; i < f.maxelm + f.maxbound; i++) {
				//f.potent[PRESS][i] = rthdsd[i];
			//}

			if (!b_on_adaptive_local_refinement_mesh) {
				exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, iflow, bextendedprint, 0, b, lb);
			}
			else {
				// Экспорт в программу tecplot температуры.
				//С АЛИС сетки.


				ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, t, fglobal, 0, b, lb);
			}
			// printf("temperature calculate begin now... OK.\n");
			 //getchar(); // debug avtosave
			//system("pause");//VX
		}


	//rfluentres.res_vx=fluent_residual_for_x(f.slau[VX], f.slau_bon[VX], f.potent[VX], f.maxelm, f.maxbound); // невязка по формуле fluent.
	//rfluentres.res_vx = fluent_residual_for_x_new(f.slau[VX], f.slau_bon[VX], f.potent[VX], f.maxelm, f.maxbound, rthdsd,f.alpha[VX]); // невязка по формуле fluent.
	bdontstartsolver = false;
	//interpolatevel(f, lw, ls,  w, VX); // ухудшает сходимость
	//getchar();
	
	// 28.07.2016
	//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, inumiter, bextendedprint);
	//getchar(); // debug
	
	
	if (bprintmessage) {
		printf("VY \n");
	}

	if (0) {
		xyplot( fglobal, flow_interior, t);
		printf("presolve. OK.\n");
	    //getchar(); // debug avtosave
		system("pause");
	}

	
	if (inumiter > 1) {
		bdontstartsolver = false;
	}
	else {
		bdontstartsolver = true;
	}

	

		//printf("VX VY\n");
		//getchar();
		solve(VELOCITY_Y_COMPONENT, res, f, fglobal, t, rhie_chow,
			s, w, b, ls, lw, lb, bonbeta,
			flow_interior, false,
			bfirst_start, toldtimestep, nullptr,
			speedoldtimestep, mfoldtimestep,
			dtimestep, btimedep, dgx, dgy, dgz,
			matlist, inumiter, bprintmessage,
			RCh, bVERYStable, nullptr, sumanb, false, false, 1.0, m,
			rthdsd, rfluentres.res_vy, lu, my_union, color, dist_max_fluid);
		

		// именно здесь верно. 04.05.2017
		rfluentres.res_vy = fluent_residual_for_x(f.slau[VELOCITY_Y_COMPONENT], f.slau_bon[VELOCITY_Y_COMPONENT], f.potent[VELOCITY_Y_COMPONENT], f.maxelm, f.maxbound, VELOCITY_Y_COMPONENT); // невязка по формуле fluent.

#pragma omp parallel for
		for (integer i = 0; i < f.maxbound; i++) {
			sumanb[VELOCITY_Y_COMPONENT][f.maxelm + i] = f.slau_bon[VELOCITY_Y_COMPONENT][i].aw;
		}

		if (bHORF_speed_on&&bHORF) {
#pragma omp parallel for
			for (integer i = 0; i < f.maxelm + f.maxbound; i++) {
				doublereal fHORF = 0.25;
				if (btimedep) {
					fHORF = 0.75;
				}
				f.potent[VELOCITY_Y_COMPONENT][i] = f.potent[VYCOR][i] + fHORF*(f.potent[VELOCITY_Y_COMPONENT][i] - f.potent[VYCOR][i]);
			}
		}
	
	//rfluentres.res_vy=fluent_residual_for_x(f.slau[VY], f.slau_bon[VY], f.potent[VY], f.maxelm, f.maxbound); // невязка по формуле fluent.
	//rfluentres.res_vy = fluent_residual_for_x_new(f.slau[VY], f.slau_bon[VY], f.potent[VY], f.maxelm, f.maxbound, rthdsd, f.alpha[VY]); // невязка по формуле fluent.
	
	bdontstartsolver = false;

	// 28.07.2016
	//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, inumiter, bextendedprint);
	//getchar(); // debug

	//interpolatevel(f, lw, ls,  w, VY); // ухудшает сходимость
	//getchar();
	if (bprintmessage) {
		printf("VZ \n");
	}
	//printf("VY VZ\n");
	//getchar();
	if (inumiter > 1) {
		bdontstartsolver = false;
	}
	else {
		bdontstartsolver = true;
	}

	

		solve(VELOCITY_Z_COMPONENT, res, f, fglobal, t, rhie_chow,
			s, w, b, ls, lw, lb, bonbeta,
			flow_interior, false,
			bfirst_start, toldtimestep, nullptr,
			speedoldtimestep, mfoldtimestep,
			dtimestep, btimedep, dgx, dgy, dgz,
			matlist, inumiter, bprintmessage,
			RCh, bVERYStable, nullptr, sumanb, false, false, 1.0, m,
			rthdsd, rfluentres.res_vz, lu, my_union, color, dist_max_fluid);
		

		// именно здесь верно. 04.05.2017
		rfluentres.res_vz = fluent_residual_for_x(f.slau[VELOCITY_Z_COMPONENT], f.slau_bon[VELOCITY_Z_COMPONENT], f.potent[VELOCITY_Z_COMPONENT], f.maxelm, f.maxbound, VELOCITY_Z_COMPONENT); // невязка по формуле fluent.

#pragma omp parallel for
		for (integer i = 0; i < f.maxbound; i++) {
			sumanb[VELOCITY_Z_COMPONENT][f.maxelm + i] = f.slau_bon[VELOCITY_Z_COMPONENT][i].aw;
		}

		

		if (bHORF_speed_on&&bHORF) {
#pragma omp parallel for
			for (integer i = 0; i < f.maxelm + f.maxbound; i++) {
				doublereal fHORF = 0.25;
				if (btimedep) {
					fHORF = 0.75;
				}
				f.potent[VELOCITY_Z_COMPONENT][i] = f.potent[VZCOR][i] + fHORF*(f.potent[VELOCITY_Z_COMPONENT][i] - f.potent[VZCOR][i]);
			}
		}
	
	//rfluentres.res_vz=fluent_residual_for_x(f.slau[VZ], f.slau_bon[VZ], f.potent[VZ], f.maxelm, f.maxbound); // невязка по формуле fluent.
	//rfluentres.res_vz = fluent_residual_for_x_new(f.slau[VZ], f.slau_bon[VZ], f.potent[VZ], f.maxelm, f.maxbound, rthdsd, f.alpha[VZ]);
	
	bdontstartsolver = false;

	// 28.07.2016
	//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, inumiter, bextendedprint);
	//getchar(); // debug

	//interpolatevel(f, lw, ls,  w, VZ); // ухудшает сходимость
	// printf("one iteration\n"); // debug
	//getchar();

	// debug:
	/*
	for (i=0; i<f.maxelm; i++) {
		if (fabs(f.potent[VY][i]+0.02089433689)<1e-6) {
		#if doubleintprecision == 1
			printf("iP==%lld, Vely=%e\n",i,f.potent[VY][i]);
		#else
			printf("iP==%d, Vely=%e\n",i,f.potent[VY][i]);
		#endif
			
			getchar();
		}
	}
	*/
	//printf("VZ end\n");
	//getchar();

	if (0) {
		xyplot( fglobal, flow_interior, t);
		printf("solve velocity. OK.\n");
	    //getchar(); // debug avtosave
		system("pause");
	}


	

	// Запоминает диагональные коэффициенты в матрице для каждой компоненты скорости:
	// диагональные коэффициенты понадобятся при вычислении поправки Рхи-Чоу.
	// В статье Винникова и Ревизникова из МАИ на стр.23 в пункте где рассматриваются
	// поправки скорости говорится, что скорость интерполированная на грань ячейки по 
	// способу Рхи-Чоу может быть не скорректирована, а вычислена заново на основе скорректированных
	// полей давления и скорости в центрах КО.
	// 10 августа 2012.
	// Теперь когда введены sumanb надо по-видимому использовать всюду их, а
	// f.diag_coef надо удалить за ненадобностью.
	// В sumanb повидимому понадобяться значения и на границе области это надо выяснить и если надо 
	// ввести в sumanb граничные значения.
	/* // закомментировано 28 февраля 2014.
	for (i=0; i<(f.maxelm+f.maxbound); i++) {
		if (i<f.maxelm) {
			f.diag_coef[VX][i]=f.slau[VX][i].ap; // VX
			f.diag_coef[VY][i]=f.slau[VY][i].ap; // VY
			f.diag_coef[VZ][i]=f.slau[VZ][i].ap; // VZ
		}
		else {
			f.diag_coef[VX][i]=f.slau_bon[VX][i-f.maxelm].aw; // VX
			f.diag_coef[VY][i]=f.slau_bon[VY][i-f.maxelm].aw; // VY 
			f.diag_coef[VZ][i]=f.slau_bon[VZ][i-f.maxelm].aw; // VZ
		}
	}
	*/
    
	// экспорт результата вычисления в программу tecplot360:
	if (0) {
		//if (inumiter>82) {
		   exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior,iflow,false,0, b, lb);
	       printf("solve momentum. OK.\n");
	       //getchar(); // debug avtosave
		   system("pause");
		//}
	}
	


	// решаем уравнение для поправки давления:
	// можно добавить инициализацию поправки давления нулём.
    if (bprintmessage) {
		printf("PAM\n"); // dbeta
	}

	

	// Массовый расход на выходной границе удовлетворяющий уравнению неразрывности.
	//mass_balance(f, lw, ls, w); // Это очень важно. См. Гаврилов Андрей.

	//doublereal* tau=new doublereal[f.maxelm+f.maxbound];

	doublereal** tau=new doublereal*[3];
	for (integer i=0; i<3; i++) tau[i]=new doublereal[f.maxelm+f.maxbound];

	bool bVERYstabtau=false;
	doublereal dusertimestep=dtimestep;
	// Отрицательный шаг по времени служит параметром сигналом означающим,
	// что нужно использовать постоянный шаг по времени равный просто среднему арифметическому
	// среди всех предварительно вычисленных псевдовременных шагов.
	if (!btimedep) {
	   dusertimestep=0.0;
	}
	// На первых 20 итерациях имеет смысл вообще использовать постоянный шаг по псевдовремени, т.к. это очень сильно восстанавливает стабильность.
	/*
	if (inumiter<20) {
		bVERYstabtau=true;	
	}
	*/
	bool boldscheme=true;
	if (0&&(inumiter==1)) {
		boldscheme=true;
	}
	// Вычисление сглаженного шага по псевдовремени.
	/*
	tau_calc(tau[0], f.maxelm, f.maxbound,
		     f.prop, f.prop_b, f.alpha, 
			 f.nvtx, f.pa, f.slau,
			 f.neighbors_for_the_internal_node, f.slau_bon,
			 btimedep, dusertimestep, 1.0,
			  inumiter,bVERYstabtau,boldscheme);*/


	doublereal CFL1=1.0; // 1.0
	// По видимому псевдовремя это не одно скалярное поле а целых
	// три скалярных поля в 3D случае, т.е. для каждой компоненты скорости
	// своё скалярное псевдовремя (как бы для каждой оси а в трёхмерном пространстве 3 оси Ox, Oy, Oz).
	tau_calc3(tau, f.maxelm, f.maxbound,
		     f.prop, f.prop_b, f.alpha, 
			 f.nvtx, f.pa, 
			 f.neighbors_for_the_internal_node, sumanb,
			 btimedep, dusertimestep, CFL1,
			 inumiter,bVERYstabtau,boldscheme);

	// Печатаем минимальное и максимальное значение tau по каждой из трех компонент скорости.
	doublereal min_tau_VX = 1.0e30; 
	doublereal min_tau_VY = 1.0e30;
	doublereal min_tau_VZ = 1.0e30;
	doublereal max_tau_VX = -1.0e30;
	doublereal max_tau_VY = -1.0e30;
	doublereal max_tau_VZ = -1.0e30;
	

	for (integer i = 0; i < f.maxelm; i++) {
		if (tau[VELOCITY_X_COMPONENT][i] < min_tau_VX) min_tau_VX = tau[VELOCITY_X_COMPONENT][i];
		if (tau[VELOCITY_Y_COMPONENT][i] < min_tau_VY) min_tau_VY = tau[VELOCITY_Y_COMPONENT][i];
		if (tau[VELOCITY_Z_COMPONENT][i] < min_tau_VZ) min_tau_VZ = tau[VELOCITY_Z_COMPONENT][i];

		if (tau[VELOCITY_X_COMPONENT][i] > max_tau_VX) max_tau_VX = tau[VELOCITY_X_COMPONENT][i];
		if (tau[VELOCITY_Y_COMPONENT][i] > max_tau_VY) max_tau_VY = tau[VELOCITY_Y_COMPONENT][i];
		if (tau[VELOCITY_Z_COMPONENT][i] > max_tau_VZ) max_tau_VZ = tau[VELOCITY_Z_COMPONENT][i];
	}

	doublereal avg_tau_VX = 0.0;
	doublereal avg_tau_VY = 0.0;
	doublereal avg_tau_VZ = 0.0;

#pragma omp parallel for reduction(+ : avg_tau_VX, avg_tau_VY, avg_tau_VZ)
	for (integer i = 0; i < f.maxelm; i++) {
		avg_tau_VX += tau[VELOCITY_X_COMPONENT][i] / f.maxelm;
		avg_tau_VY += tau[VELOCITY_Y_COMPONENT][i] / f.maxelm;
		avg_tau_VZ += tau[VELOCITY_Z_COMPONENT][i] / f.maxelm;

	}
	printf("\ntau statistics:\n");
	printf("	VX		VY		VZ\n");
	printf("min: %e %e %e\n", min_tau_VX, min_tau_VY, min_tau_VZ);
	printf("max: %e %e %e\n", max_tau_VX, max_tau_VY, max_tau_VZ);
	printf("avg: %e %e %e\n\n", avg_tau_VX, avg_tau_VY, avg_tau_VZ);
	if (1&&(inumiter < 20)) {
		doublereal m = 1.0; // (1.0 / my_amg_manager.theta_Stress);///debug //1.0

#pragma omp parallel for
		for (integer i = 0; i < f.maxelm; i++) {

			doublereal avg_tau_VX_loc = avg_tau_VX;
			if (!b_on_adaptive_local_refinement_mesh) {
				avg_tau_VX_loc = 0.0;
				avg_tau_VX_loc += 0.1666 * (tau[VELOCITY_X_COMPONENT][f.slau[VELOCITY_X_COMPONENT][i].iE] +
					tau[VELOCITY_X_COMPONENT][f.slau[VELOCITY_X_COMPONENT][i].iW] + tau[VELOCITY_X_COMPONENT][f.slau[VELOCITY_X_COMPONENT][i].iN] +
					tau[VELOCITY_X_COMPONENT][f.slau[VELOCITY_X_COMPONENT][i].iS] + tau[VELOCITY_X_COMPONENT][f.slau[VELOCITY_X_COMPONENT][i].iT] +
					tau[VELOCITY_X_COMPONENT][f.slau[VELOCITY_X_COMPONENT][i].iB]);
			}

			tau[VELOCITY_X_COMPONENT][i] = 1.0 / ((1.0 / (m*avg_tau_VX_loc)) + (1.0 / tau[VELOCITY_X_COMPONENT][i]));

			doublereal avg_tau_VY_loc = avg_tau_VY;
			if (!b_on_adaptive_local_refinement_mesh) {
				avg_tau_VY_loc = 0.0;
				avg_tau_VY_loc += 0.1666 * (tau[VELOCITY_Y_COMPONENT][f.slau[VELOCITY_Y_COMPONENT][i].iE] +
					tau[VELOCITY_Y_COMPONENT][f.slau[VELOCITY_Y_COMPONENT][i].iW] + tau[VELOCITY_Y_COMPONENT][f.slau[VELOCITY_Y_COMPONENT][i].iN] +
					tau[VELOCITY_Y_COMPONENT][f.slau[VELOCITY_Y_COMPONENT][i].iS] + tau[VELOCITY_Y_COMPONENT][f.slau[VELOCITY_Y_COMPONENT][i].iT] +
					tau[VELOCITY_Y_COMPONENT][f.slau[VELOCITY_Y_COMPONENT][i].iB]);
			}

			tau[VELOCITY_Y_COMPONENT][i] = 1.0 / ((1.0 / (m*avg_tau_VY_loc)) + (1.0 / tau[VELOCITY_Y_COMPONENT][i]));

			doublereal avg_tau_VZ_loc = avg_tau_VZ;
			if (!b_on_adaptive_local_refinement_mesh) {
				avg_tau_VZ_loc = 0.0;
				avg_tau_VZ_loc += 0.1666 * (tau[VELOCITY_Z_COMPONENT][f.slau[VELOCITY_Z_COMPONENT][i].iE] +
					tau[VELOCITY_Z_COMPONENT][f.slau[VELOCITY_Z_COMPONENT][i].iW] + tau[VELOCITY_Z_COMPONENT][f.slau[VELOCITY_Z_COMPONENT][i].iN] +
					tau[VELOCITY_Z_COMPONENT][f.slau[VELOCITY_Z_COMPONENT][i].iS] + tau[VELOCITY_Z_COMPONENT][f.slau[VELOCITY_Z_COMPONENT][i].iT] +
					tau[VELOCITY_Z_COMPONENT][f.slau[VELOCITY_Z_COMPONENT][i].iB]);
			}

			tau[VELOCITY_Z_COMPONENT][i] = 1.0 / ((1.0 / (m*avg_tau_VZ_loc)) + (1.0 / tau[VELOCITY_Z_COMPONENT][i]));
		}
	}

	/*
	if (inumiter==2) {
	for (i=0; i<f.maxelm; i++) {
	#if doubleintprecision == 1
		printf("%lld %e, %e %e\n",i, tau[VX][i],tau[VY][i],tau[VZ][i]);
	#else
		printf("%d %e, %e %e\n",i, tau[VX][i],tau[VY][i],tau[VZ][i]);
	#endif
		
		if (i%10==0) getchar();
	}
	}
	*/
	//Использование математического ожидания от псевдовремени (среднее арифметическое среди всех значений)
	// приводит к сильнейшей расходимости уже на второй итерации,
	// так как поле давления получается сильно завышенным.
	// Если же использовать минимальное среди всех возможных значений tau
	// (все значения псевдовремени одинаковы и равны минимальному значению среди исходных)
	// то сходимость продолжает иметь место.
	
	/*
	// Условие Куранта Фридрихса Леви.
	for (integer i=0; i<f.maxelm; i++) {
		// вычисление размеров текущего контрольного объёма:
	    doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	    volume3D(i, f.nvtx, f.pa, dx, dy, dz);
		doublereal CFL=f.alpha[VX]/(1.0-alpha[VX]); // Условие Куранта Фридрихса Леви.
		tau[VX][i]=fmin(fabs(CFL/(f.potent[VX][i])/dx),tau[VX][i]);
		if (f.neighbors_for_the_internal_node[ESIDE][i].iNODE1>=f.maxelm) tau[VX][f.neighbors_for_the_internal_node[ESIDE][i].iNODE1]=tau[VX][i];
		if (f.neighbors_for_the_internal_node[WSIDE][i].iNODE1>=f.maxelm) tau[VX][f.neighbors_for_the_internal_node[WSIDE][i].iNODE1]=tau[VX][i];
		tau[VY][i]=fmin(fabs(CFL/(f.potent[VY][i])/dy),tau[VY][i]);
		if (f.neighbors_for_the_internal_node[NSIDE][i].iNODE1>=f.maxelm) tau[VY][f.neighbors_for_the_internal_node[NSIDE][i].iNODE1]=tau[VY][i];
		if (f.neighbors_for_the_internal_node[SSIDE][i].iNODE1>=f.maxelm) tau[VY][f.neighbors_for_the_internal_node[SSIDE][i].iNODE1]=tau[VY][i];
		tau[VZ][i]=fmin(fabs(CFL/(f.potent[VZ][i])/dz),tau[VZ][i]);
		if (f.neighbors_for_the_internal_node[TSIDE][i].iNODE1>=f.maxelm) tau[VZ][f.neighbors_for_the_internal_node[TSIDE][i].iNODE1]=tau[VZ][i];
		if (f.neighbors_for_the_internal_node[BSIDE][i].iNODE1>=f.maxelm) tau[VZ][f.neighbors_for_the_internal_node[BSIDE][i].iNODE1]=tau[VZ][i];
	}

	// Условие на сеточное число Рейнольдса.
	for (integer i=0; i<f.maxelm; i++) {
		// вычисление размеров текущего контрольного объёма:
	    doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	    volume3D(i, f.nvtx, f.pa, dx, dy, dz);
		doublereal Rec=1.5; // сеточное число Рейнольдса
		tau[VX][i]=fmin((dx*dx)/(f.prop[RHO][i]*f.prop[MU][i]*Rec),tau[VX][i]);
		if (f.neighbors_for_the_internal_node[ESIDE][i].iNODE1>=f.maxelm) tau[VX][f.neighbors_for_the_internal_node[ESIDE][i].iNODE1]=tau[VX][i];
		if (f.neighbors_for_the_internal_node[WSIDE][i].iNODE1>=f.maxelm) tau[VX][f.neighbors_for_the_internal_node[WSIDE][i].iNODE1]=tau[VX][i];
		tau[VY][i]=fmin((dy*dy)/(f.prop[RHO][i]*f.prop[MU][i]*Rec),tau[VY][i]);
		if (f.neighbors_for_the_internal_node[NSIDE][i].iNODE1>=f.maxelm) tau[VY][f.neighbors_for_the_internal_node[NSIDE][i].iNODE1]=tau[VY][i];
		if (f.neighbors_for_the_internal_node[SSIDE][i].iNODE1>=f.maxelm) tau[VY][f.neighbors_for_the_internal_node[SSIDE][i].iNODE1]=tau[VY][i];
		tau[VZ][i]=fmin((dz*dz)/(f.prop[RHO][i]*f.prop[MU][i]*Rec),tau[VZ][i]);
		if (f.neighbors_for_the_internal_node[TSIDE][i].iNODE1>=f.maxelm) tau[VZ][f.neighbors_for_the_internal_node[TSIDE][i].iNODE1]=tau[VZ][i];
		if (f.neighbors_for_the_internal_node[BSIDE][i].iNODE1>=f.maxelm) tau[VZ][f.neighbors_for_the_internal_node[BSIDE][i].iNODE1]=tau[VZ][i];
	}

	// Условие на сеточное число Рейнольдса.
	for (integer i=0; i<f.maxelm; i++) {
		// вычисление размеров текущего контрольного объёма:
	    doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	    volume3D(i, f.nvtx, f.pa, dx, dy, dz);
		doublereal Rec=1.5; // сеточное число Рейнольдса
		tau[VX][i]=fmin((f.prop[RHO][i]*f.prop[MU][i]*Rec)/(fmax(1e-15,f.potent[VX][i]*f.potent[VX][i])),tau[VX][i]);
		if (f.neighbors_for_the_internal_node[ESIDE][i].iNODE1>=f.maxelm) tau[VX][f.neighbors_for_the_internal_node[ESIDE][i].iNODE1]=tau[VX][i];
		if (f.neighbors_for_the_internal_node[WSIDE][i].iNODE1>=f.maxelm) tau[VX][f.neighbors_for_the_internal_node[WSIDE][i].iNODE1]=tau[VX][i];
		tau[VY][i]=fmin((f.prop[RHO][i]*f.prop[MU][i]*Rec)/(fmax(1e-15,f.potent[VY][i]*f.potent[VY][i])),tau[VY][i]);
		if (f.neighbors_for_the_internal_node[NSIDE][i].iNODE1>=f.maxelm) tau[VY][f.neighbors_for_the_internal_node[NSIDE][i].iNODE1]=tau[VY][i];
		if (f.neighbors_for_the_internal_node[SSIDE][i].iNODE1>=f.maxelm) tau[VY][f.neighbors_for_the_internal_node[SSIDE][i].iNODE1]=tau[VY][i];
		tau[VZ][i]=fmin((f.prop[RHO][i]*f.prop[MU][i]*Rec)/(fmax(1e-15,f.potent[VZ][i]*f.potent[VZ][i])),tau[VZ][i]);
		if (f.neighbors_for_the_internal_node[TSIDE][i].iNODE1>=f.maxelm) tau[VZ][f.neighbors_for_the_internal_node[TSIDE][i].iNODE1]=tau[VZ][i];
		if (f.neighbors_for_the_internal_node[BSIDE][i].iNODE1>=f.maxelm) tau[VZ][f.neighbors_for_the_internal_node[BSIDE][i].iNODE1]=tau[VZ][i];
	}
	*/
	/*
	doublereal tauX=1.0e+30, tauY=1.0e+30, tauZ=1.0e+30;
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		tauX=fmin(tauX,tau[VX][i]);
		tauY=fmin(tauY,tau[VY][i]);
		tauZ=fmin(tauZ,tau[VZ][i]);
	}
	//doublereal rn=1.0*(f.maxelm+f.maxbound);
	//tauX/=rn;
	//tauY/=rn;
	//tauZ/=rn;
	for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		tau[VX][i]=0.1*tauX;
		tau[VY][i]=0.1*tauY;
		tau[VZ][i]=0.1*tauZ;
	}
	*/

	// Коррекция скорости и массового потока на выходной грнице
	// после коррекции скоростей см. Гаврилов Андрей.
	//mass_balance(f, lw, ls, w);

	
	//printf("PAM\n");
	//getchar();
	doublereal** rsumanbstuff=nullptr; // nullptr поинтер под сумму диагональных коэффициентов
	bool bhighorder_pressure=false;
	bool bdeltafinish=true; // если true то возвращаем значение массового потока на грани. Если false то невозвращаем. 
	// false может быть только в том случае если уравнение для поправки давления потом решается повторно с целью повысить точность
	// аппроксимации.
	// RANS Спалларт Аллмарес.
	if (!((f.iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES)||
		(f.iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)||
		(f.iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS))) {
		if (bflag_free_memory_cfd) {
			m.bsignalfreeCRScfd = true; // Освобождаем память.
		}
	}
	

	doublereal rfluentResPAM = 0.0;

	integer inum_iter_loc = 1000000;
	if ((fabs(dgx) > 1.0e-20) || (fabs(dgy) > 1.0e-20) || (fabs(dgz) > 1.0e-20)) {
		// Натуральная конвекция в cfd.
		bool b_natural_convection = true;
		for (int k21 = 0; k21 < lw; k21++) {
			if ((w[k21].Vx*w[k21].Vx + w[k21].Vy*w[k21].Vy + w[k21].Vz*w[k21].Vz) > 1.0e-20) {
				b_natural_convection = false; // Вынужденная конвекция.
			}
		}
		if (b_natural_convection) {
			// Ограничение работает только при естественной конвекции.
			inum_iter_loc = inumiter;
		}
	}


    solve(PAM,continity,f,fglobal,t,rhie_chow,
		  s, w, b, ls, lw, lb, 1.0, 
		  flow_interior, false, bfirst_start,
		  toldtimestep, nullptr, speedoldtimestep,
		  mfoldtimestep,dtimestep,btimedep,
		  dgx,dgy,dgz,matlist, inum_iter_loc,
		  bprintmessage, RCh,false,
		  tau, sumanb /*rsumanbstuff*/,bhighorder_pressure,
		  bdeltafinish, 1.0, m, rthdsd, rfluentResPAM,
		lu, my_union, color, dist_max_fluid);

	
	// 9 августа 2016 года. (нижняя релаксация для поправки давления).
	// Сильно портит сходимость.
	if (0&&bHORF) {

#pragma omp parallel for
		for (integer i = 0; i < f.maxelm + f.maxbound; i++) {
			doublereal fHORF = 0.25;
			if (btimedep) {
				fHORF = 0.75;
			}
			f.potent[PAM][i] = f.potent[PAMOLDITER][i] + fHORF*(f.potent[PAM][i] - f.potent[PAMOLDITER][i]);
			f.potent[PAMOLDITER][i] = f.potent[PAM][i];
		}
	}

	/*
	bhighorder=false;
	bdeltafinish=true;
	// здесь возможно нужен цикл в несколько проходов пока не будет достигнута необходимая точность аппроксимации.
	// dbeta 1.2 обеспечивает третий порядок на границе области.
	solve(PAM,continity,f,fglobal,t,rhie_chow,
		  s, w, b, ls, lw, lb, 1.0, 
		  flow_interior, false, bfirst_start,
		  toldtimestep,speedoldtimestep,
		  mfoldtimestep,dtimestep,btimedep,
		  dgx,dgy,dgz,matlist,inumiter,
		  bprintmessage, RCh,false, false, tau, rsumanbstuff,bhighorder,bdeltafinish);*/
		  

	// 28.07.2016
	//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, inumiter, bextendedprint);
	//getchar(); // debug

	//doublereal* temp = new doublereal[f.maxelm + f.maxbound];
	//for (integer i = 0; i < f.maxelm + f.maxbound; i++) {
		//temp[i] = f.potent[VZ][i];
		//f.potent[VZ][i] = rthdsd[i];
	//}


	// экспорт результата вычисления в программу tecplot360:
	if (0) {
		//for (integer i = 0; i < f.maxelm + f.maxbound; i++) {
			//f.potent[PRESS][i] = rthdsd[i];
		//}

		if (!b_on_adaptive_local_refinement_mesh) {
			exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, iflow, bextendedprint, 0, b, lb);
		}
		else {
			// Экспорт в программу tecplot температуры.
			//С АЛИС сетки.			

			ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, t, fglobal, 0, b, lb);
		}
	    // printf("temperature calculate begin now... OK.\n");
	     //getchar(); // debug avtosave
		 system("pause");//PAM
	}

	//for (integer i = 0; i < f.maxelm + f.maxbound; i++) {
		//f.potent[VZ][i]= temp[i];
	//}
	//delete[] temp;

	// 13 08 2015 вынесено в глобальную область видимости см.
	// my_unsteady_temperature.c файл.
	//doublereal* xb=new doublereal[f.maxelm+f.maxbound];

#pragma omp parallel for
	for (integer i=0; i<f.maxelm; i++) {
		xb[i]=f.slau[PAM][i].b;
	}
#pragma omp parallel for
	for (integer i=0; i<f.maxbound; i++) {
		xb[i+f.maxelm]=f.slau_bon[PAM][i].b;
	}
	rfluentres.res_no_balance=no_balance_mass_flux_fluent(xb, rfluentres.operating_value_b, f.maxelm+f.maxbound);
	//delete xb; // не забываем освобождать память.
	
	// Проверка получено ли решение СЛАУ.
	for (integer i = 0; i < f.maxelm + f.maxbound; i++) {
		if (f.potent[PAM][i] != f.potent[PAM][i]) {

			if (rfluentres.res_no_balance < 1.0e-10) {
				// Машинный ноль
				delete[] f.potent[PAM];
				f.potent[PAM] = new doublereal[f.maxelm + f.maxbound];
#pragma omp parallel for
				for (integer j = 0; j < f.maxelm + f.maxbound; j++) {
					f.potent[PAM][j] = 0.0;
				}
				break; // досрочный выход из цикла for.
			}
			else {
				printf("Solve(PAM) error NAN or INF f.potent[PAM][%lld]=%e\n", i, f.potent[PAM][i]);
				system("pause");
			}
		}
	}

	/*
	// Фильтрование не устраняет неровности.
	// Фильтрование действует следующим образом:
	// Пусть грубая сетка стыкуется с грубой сеткой через прослойку где сетка 
	// мелкая. Каждая из трёх сеток квазиравномерная. На стыке двух квазиравномерных сеток 
	// наблюдается разрыв (скачёк) первой производной. Этот скачёк изначально присутствует в самом
	// поле полученном после решения СЛАУ. Фильтрование не уберёт данный скачёк оно лишь усилит его:
	// Каждая из квазиравномерных областей становится ещё более равномерной, а место стыка двух 
	// квазиравномерных областей  как терпело разрыв так и продолжает его терпеть.
	// Теперь можно применить фильтрование к найденной поправке 
	// давления так чтобы она стала гладкой функцией.
	doublereal* potentPAMfiltr=new doublereal[f.maxelm+f.maxbound];
		
	doublereal *nullpointer=nullptr;
	// регулятор сглаженности поля градиента.
	integer iciclenumerical=300;  // Чем больше циклов тем больше сглаженность.
	for (integer icicle=0; icicle<iciclenumerical; icicle++) {
		double_average_potent(f.potent[PAM], potentPAMfiltr,
	                         f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, nullpointer,
							 SIMPSON_FILTR, f.border_neighbor,0); // VOLUME_AVERAGE_FILTR
		// копирование.
	    for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		    f.potent[PAM][i]=potentPAMfiltr[i];
	    }
	}

	// Освобождение оперативной памяти.
	delete potentPAMfiltr;
	*/

	//printf("Ok\n");
	//getchar();

	//if (0&&(f.bPressureFix)) {
	if (0&&(f.bLR1free)) {
	//if (1) {
		// Условие однозначности состоит в том чтобы
		// среднее значение поправки давления было равно нулю.

		//printf("regularization condition for Pamendment...\n");
		//printf("please, press any key to continue...\n");
		//getchar(); // отладочный контроль.

		// Следуя П.Н. Вабищевичу, А.Н.Павлову, А.Г.Чурбанову
		// будем использовать условие однозначности для 
		// определения поля поправки давления в том случае
		// если на всей границе расчётной области стоит 
		// однородное условие Неймана.
	    // Условие однозначности:
	    doublereal pamscal1=0.0;
	    // pamscal1=Scal(f.potent[PAM],1.0); должно равняться нулю
		/*
		// Этот кусок кода был задействован раньше.
		// Его следует раскоментировать при откате на прошлую версию.
	    for (i=0; i<(f.maxelm+f.maxbound); i++) pamscal1+=f.potent[PAM][i];
		pamscal1=(doublereal)(pamscal1/(f.maxelm+f.maxbound));
		for (i=0; i<(f.maxelm+f.maxbound); i++) f.potent[PAM][i]-=pamscal1;
		*/

		doublereal Vol=0.0;
		for (integer iP=0; iP<(f.maxelm); iP++) {
		    // вычисление размеров текущего контрольного объёма:
	        doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	        volume3D(iP, f.nvtx, f.pa, dx, dy, dz);
			Vol+=dx*dy*dz;

			pamscal1+=f.potent[PAM][iP]*dx*dy*dz;
		}
		pamscal1=(doublereal)(pamscal1/(f.maxelm*Vol));
		for (integer iP=0; iP<(f.maxelm+f.maxbound); iP++) {
			f.potent[PAM][iP]-=pamscal1;
		}

		// Возможно здесь лучше взять интеграл по объёму от поправки давления,
		// а потом разделить на объём рассматриваемой FLUID зоны.
		// Внесённое изменение должно быть проверено TODO.

		// Интерполяция поправки давления из внутренности расчётной 
	    // области на границу.
	    //free_pressure(f);
	}
	//else free_pressure(f);

	// Интерполяция поправки давления из внутренности расчётной 
	// области на границу.
	// под интерполяцией понимается линейная интерполяция.
	// Физический смысл в том, что однородные граничные условия Неймана для поправки давления
	// дают равенство поправки давления в граничном и ближайшем внутреннем узле в то время как
	// давление меняется во всей области включая приграничные узлы например линейно в случае отрезка трубы.
	// т.е. реально по физике должно быть давление в граничном узле и ближайшем внутреннем не равны, они должны быть
	// такими чтобы соответсвовать по непрерывности градиенту давления во внутренности расчётной области.
	// Вывод условия Неймана некорректны и надо задавать другие граничные условия, например неизменность градиента давления,
	// т.е. нужно например задать равенство второй производной нулю а это влет трёхточечную связь на границе что невозможно по идее.
	// (возможно если воспользоваться методом отложенной коррекции).
	// Ещё один способ добится равенства второй производной нулю на границе это продолжить давление по непрерывности изнутри расчётной области на границу.
	// Этот вариант здесь рассматривается.
	// при выключенном freepressure при прочих равных условиях в отличие от включённого freepressure амплитуда начальных колебаний невязки continity возрастает на 17%.
	// На конечный процесс сходимости наличие freepressure не влияет и его можно смело отключить. Никаких аномальных всплесков 
	// вблизи выходной границы при отключённом freepressure не наблюдается.
	//free_pressure(f); // Вообще говоря это искуственность и она может портить сходимость (откажемся от неё).


	//printf("PAM end\n");
	//getchar();

	
	
// Вычисление градиентов поправки давления:
	#ifdef _OPENMP

	// С параллельным кодом одна итерация SIMPLE алгоритма выполняется на 2.2% быстрее.
	if (bparallelizm_old) {
		printf("error bparallelizm_old\n");
		system("pause");
		if (inumcore == 2) {
			if (nd.b0.active) {

#pragma omp parallel 
				{
#pragma omp sections 
					{
#pragma omp section
						{
							// первый поток
							for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
								integer iPloc = f.ifrontregulationgl[iscan_par];
								if (iPloc < f.maxelm) {
									// градиенты поправки давления для внутренних КО.
									green_gaussPAM(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
								}

							}
						}
#pragma omp section
						{
							// второй поток
							for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
								integer iPloc = f.ifrontregulationgl[iscan_par];
								if (iPloc < f.maxelm) {
									// градиенты поправки давления для внутренних КО.
									green_gaussPAM(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
								}
							}
						} // section
					} // sections
				} // parallel

				// серийный смыкающий кусок
				for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
					integer iPloc = f.ifrontregulationgl[iscan_par];
					if (iPloc < f.maxelm) {
						// градиенты поправки давления для внутренних КО.
						green_gaussPAM(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
					}
				}


			}
		}

		if (inumcore == 2) {
			if (nd.b0.active) {

#pragma omp parallel 
				{
#pragma omp sections 
					{
#pragma omp section
						{

							// первый поток
							for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
								integer iPloc = f.ifrontregulationgl[iscan_par];
								if (iPloc < f.maxelm) {
									// градиенты поправки давления для внутренних КО.
									green_gaussPAM(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
								}

							}
						}
#pragma omp section
						{

							// второй поток
							for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
								integer iPloc = f.ifrontregulationgl[iscan_par];
								if (iPloc < f.maxelm) {
									// градиенты поправки давления для внутренних КО.
									green_gaussPAM(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
								}
							}
						} // section
					} // sections
				} // parallel

				// серийный смыкающий кусок
				for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
					integer iPloc = f.ifrontregulationgl[iscan_par];
					if (iPloc < f.maxelm) {
						// градиенты поправки давления для внутренних КО.
						green_gaussPAM(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
					}
				}


			}
		}
	}
	else {
		// Можно безболезненно распараллелить т.к. доступ между ячейками i независим.

		// Вычисление градиентов поправки давления:
#pragma omp parallel for
		for (integer i = 0; i < f.maxelm; i++) {
			// градиенты поправки давления для внутренних КО.
			green_gaussPAM(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
		}
#pragma omp parallel for
		for (integer i = 0; i < f.maxelm; i++) {
			// градиенты скоростей для граничных КО.
			green_gaussPAM(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
		}
	}
#else

	if (0) {
		// bLRfree==false
		if (f.bLR1free) {
			printf("bLRfree_on\n");
			system("pause");
		}
		else {
			printf("bLRfree_off\n");
			system("pause");
		}
	}

	// Вычисление градиентов поправки давления:
	for (integer i=0; i<f.maxelm; i++) {
		// градиенты поправки давления для внутренних КО.
		green_gaussPAM(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.border_neighbor,ls,lw,w,f.bLR1free,t.ilevel_alice,f.ptr);
	}
	for (integer i=0; i<f.maxelm; i++) {
		// градиенты скоростей для граничных КО.
	    green_gaussPAM(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.border_neighbor,ls,lw,w,f.bLR1free,t.ilevel_alice,f.ptr);
    }

#endif


	// 10 августа 2012 года отладочная печать.
	/*
	for (i=0; i<f.maxelm+f.maxbound; i++) {
	    f.potent[FBUF][i]=f.potent[GRADYPAM][i];
	}
	*/
	if (0) {
		xyplot( fglobal, flow_interior, t);
		printf("xy plot pam solve. OK.\n");
	   // getchar(); // debug avtosave
		system("pause");

		//exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior,iflow);
	    //printf("solve vel ok fbuf gradypam... OK.\n");
	    //getchar(); // debug avtosave
	}

	

	// Минимакстное сглаживание градиента поправки давления.
    //green_gaussPAMminmax(f.potent, f.neighbors_for_the_internal_node, f.maxelm, f.maxbound);

	/*
	// Диффузионный поток аппроксимирован с первым порядком на неравномерной сетке и
	// со вторым порядком на равномерной. Реальность заставляет нас использовать существенно 
	// неравномерные сетки. Рассмотрим проблему на примере уравнения для поправки давления.
	// Уравнение для порправки давления это чисто диффузионная задача с симметричной положительно
	// определённой матрицей. Итак мы имеем аппроксимацию первого порядка на неравномерной сетке.
	// Это значит что идеально точно найденное поле поправки давления после решения СЛАУ будет
	// непрерывно всюду, но при этом первая производная от поправки давления может быть и обязательно будет
	// РАЗРЫВНОЙ или ИСПЫТЫВАТЬ ЗНАЧИТЕЛЬНЫЕ КОЛЕБАНИЯ. Здесь предполагается, что причина этому неравномерность
	// расчётной сетки. Итак производная от найенной функции суть величина подверженная шуму. Посколько на основе 
	// этой подверженной шуму величине основаны другие величины - например скорректированная скорость и т.д. по цепочке.
	// То эти возмущения наследуются дальше усиливаются и мы имеем РАСХОДИМОСТЬ. Предполагается два способа возможной борьбы с 
	// расходимостью: 1. Сглаживание результата с помощью применения фильтра (многократно). Это простой рабочий вариант.
	// 2. увеличение точности аппроксимации дискреного аналога до второго так чтобы взятая по полученному решению первая производная
	// была гладкой функцией (непрерывной).
	// Вычисление и обязательное сглаживание градиента давления.
	// Под сглаживанием понимается усредняющее фильтрование.
	

	for (i=0; i<f.maxelm+f.maxbound; i++) {
		f.potent[FBUF][i]=f.potent[GRADYPAM][i];
	}
	if (1) {
		xyplot( fglobal, flow_interior, t);
		printf("grad PAM calc. OK.\n");
	    getchar(); // debug avtosave
	}

	// Теперь нужно применить фильтрование:
	doublereal** potentgradPAMfiltr=new doublereal*[3];
	for (i=0; i<3; i++) {
		potentgradPAMfiltr[i]=new doublereal[f.maxelm+f.maxbound];
	}
	
	doublereal *nullpointer=nullptr;
	// регулятор сглаженности поля градиента.
	integer iciclenumerical=3;  // Чем больше циклов тем больше сглаженность.
	for (integer icicle=0; icicle<iciclenumerical; icicle++) {
		double_average_potent(potentgradPAMfiltr[0], f.potent[GRADXPAM], 
	                         f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, nullpointer,
							 VOLUME_AVERAGE_FILTR, f.border_neighbor,0);
		double_average_potent(potentgradPAMfiltr[1], f.potent[GRADYPAM], 
	                         f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, nullpointer,
							 VOLUME_AVERAGE_FILTR, f.border_neighbor,0);
		double_average_potent(potentgradPAMfiltr[2], f.potent[GRADZPAM], 
	                         f.maxelm, f.maxbound, f.neighbors_for_the_internal_node, 
							 f.nvtx, f.pa, nullpointer,
							 VOLUME_AVERAGE_FILTR, f.border_neighbor,0);
		// копирование.
	    for (integer i=0; i<f.maxelm+f.maxbound; i++) {
		    f.potent[GRADXPAM][i]=potentgradPAMfiltr[0][i];
	    	f.potent[GRADYPAM][i]=potentgradPAMfiltr[1][i];
		    f.potent[GRADZPAM][i]=potentgradPAMfiltr[2][i];
	    }
	}

	// Освобождение оперативной памяти.
	for (i=0; i<3; i++) {
		delete potentgradPAMfiltr[i];
	}
	delete potentgradPAMfiltr;
	*/
	
	// Корректирующая поправка к поправке давления для устранения неправильного влияния
	// условий Неймана.
	// По ряду тестов включение следующих двух функций ухудшает сходимость.
    //correctpressureoutlet(f, lw, ls, w);
	//free_pressure(f);
	
	// экспорт результата вычисления в программу tecplot360:
	if (0) {
		if (1||inumiter>82) {
		   exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior,iflow,false,0, b, lb);
	       printf("solve pressure. OK.\n");
	       //getchar(); // debug avtosave
		   system("pause");
		}
	}
	if (0) {
		xyplot( fglobal, flow_interior, t);
		printf("solve pressure. OK.\n");
	    //getchar(); // debug avtosave
		system("pause");
	}

	// коррекция давления:
	doublereal ralphaP=1.0; // инициализация
	if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLEC_Van_Doormal_and_Raithby) {
		// SIMPLEC 
		ralphaP=1.0;
		//ralphaP=1.0-(f.alpha[VX]+f.alpha[VY]+f.alpha[VZ])/3.0;
	}
	if (iSIMPLE_alg== SIMPLE_CFD_ALGORITHM::SIMPLE_Carretto) {
		// SIMPLE
		ralphaP=f.alpha[PRESS];
	}
	// Только для внутренних КО.
#pragma omp parallel for
	for (integer i=0; i<f.maxelm; i++)  {
		//f.potent[PRESS][i]+=ralphaP*f.potent[PAM][i];
		// 06.05.2017 Со ссылкой на лекцию о SIMPLE алгоритме.
		// презентация: CFD - Solution Algorithms
		// SOE3213/4: CFD Lecture 3.
		// p_zv_zv - Presure звёздочка звёздочка.
		doublereal p_zv_zv = f.potent[PRESS][i] + f.potent[PAM][i];
		f.potent[PRESS][i] = ralphaP*p_zv_zv + (1.0 - ralphaP)*f.potent[PRESS][i];
	}

#pragma omp parallel for
	for (integer i=f.maxelm; i<(f.maxelm+f.maxbound); i++)  {
		// по-моему лучше не фиксировать давления, а позволить ему 
		// корректироваться в течении всего расчёта.
		// Это требует рассмотрения TODO.
		//if (f.slau_bon[PAM][i-f.maxelm].iI>-1) {

            // Старая версия (раскоментировать лишь одну строчку).
			// стоит однородное условие Неймана
           // f.potent[PRESS][i]+=ralphaP*f.potent[PAM][i];

		   // 06.05.2017 Со ссылкой на лекцию о SIMPLE алгоритме.
		   // презентация: CFD - Solution Algorithms
		   // SOE3213/4: CFD Lecture 3.
		   // p_zv_zv - Presure звёздочка звёздочка.
			doublereal p_zv_zv = f.potent[PRESS][i] + f.potent[PAM][i];
			f.potent[PRESS][i] = ralphaP*p_zv_zv + (1.0 - ralphaP)*f.potent[PRESS][i];

			//printf("corect boundary pressure...\n");  // debug
			//getchar();
		//}
/*
			integer inumber=i-f.maxelm;
			if ((f.border_neighbor[inumber].MCB>=ls) && (f.border_neighbor[inumber].MCB<(ls+lw)) && w[f.border_neighbor[inumber].MCB-ls].bpressure) {
				f.potent[PRESS][i]=w[f.border_neighbor[inumber].MCB-ls].P;
			}
			*/
	}

	
// Вычисление градиентов давления:
	// на основе скорректированного поля давления.
	// Градиенты давления понадобятся при вычислении поправки Рхи-Чоу.
	#ifdef _OPENMP

	// С параллельным кодом одна итерация SIMPLE алгоритма выполняется на 0.9% быстрее.
	if (bparallelizm_old) {
		printf("error bparallelizm_old\n");
		system("pause");
		if (inumcore == 2) {
			if (nd.b0.active) {

#pragma omp parallel 
				{
#pragma omp sections 
					{
#pragma omp section
						{
							// первый поток
							for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
								integer iPloc = f.ifrontregulationgl[iscan_par];
								if (iPloc < f.maxelm) {
									// градиенты поправки давления для внутренних КО.
									green_gaussPRESS(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
								}

							}
						}
#pragma omp section
						{
							// второй поток
							for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
								integer iPloc = f.ifrontregulationgl[iscan_par];
								if (iPloc < f.maxelm) {
									// градиенты поправки давления для внутренних КО.
									green_gaussPRESS(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
								}
							}
						} // section
					} // sections
				} // parallel

				// серийный смыкающий кусок
				for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
					integer iPloc = f.ifrontregulationgl[iscan_par];
					if (iPloc < f.maxelm) {
						// градиенты поправки давления для внутренних КО.
						green_gaussPRESS(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
					}
				}


			}
		}

		if (inumcore == 2) {
			if (nd.b0.active) {

#pragma omp parallel 
				{
#pragma omp sections 
					{
#pragma omp section
						{

							// первый поток
							for (integer iscan_par = nd.b0.ileft_start; iscan_par <= nd.b0.ileft_finish; iscan_par++) {
								integer iPloc = f.ifrontregulationgl[iscan_par];
								if (iPloc < f.maxelm) {
									// градиенты поправки давления для внутренних КО.
									green_gaussPRESS(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
								}

							}
						}
#pragma omp section
						{

							// второй поток
							for (integer iscan_par = nd.b0.iright_start; iscan_par <= nd.b0.iright_finish; iscan_par++) {
								integer iPloc = f.ifrontregulationgl[iscan_par];
								if (iPloc < f.maxelm) {
									// градиенты поправки давления для внутренних КО.
									green_gaussPRESS(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
								}
							}
						} // section
					} // sections
				} // parallel

				// серийный смыкающий кусок
				for (integer iscan_par = nd.b0.iseparate_start; iscan_par <= nd.b0.iseparate_finish; iscan_par++) {
					integer iPloc = f.ifrontregulationgl[iscan_par];
					if (iPloc < f.maxelm) {
						// градиенты поправки давления для внутренних КО.
						green_gaussPRESS(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
					}
				}


			}
		}
	}
	else {
		// Вычисление градиентов давления:
	    // на основе скорректированного поля давления.
	    // Градиенты давления понадобятся при вычислении поправки Рхи-Чоу.
#pragma omp parallel for
		for (integer i = 0; i < f.maxelm; i++) {
			// градиенты давления для внутренних КО.
			green_gaussPRESS(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
		}
#pragma omp parallel for
		for (integer i = 0; i < f.maxelm; i++) {
			// градиенты давления для граничных КО.
			green_gaussPRESS(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.border_neighbor, ls, lw, w, f.bLR1free, t.ilevel_alice, f.ptr);
		}
	}
#else

	// Вычисление градиентов давления:
	// на основе скорректированного поля давления.
	// Градиенты давления понадобятся при вычислении поправки Рхи-Чоу.
	for (integer i=0; i<f.maxelm; i++) {
		// градиенты давления для внутренних КО.
	    green_gaussPRESS(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false,f.border_neighbor,ls,lw,w,f.bLR1free,t.ilevel_alice,f.ptr);
	}
	for (integer i=0; i<f.maxelm; i++) {
		// градиенты давления для граничных КО.
	    green_gaussPRESS(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true,f.border_neighbor,ls,lw,w,f.bLR1free,t.ilevel_alice,f.ptr);
    }
	

	
#endif

	//printf("CORRECT start\n");
	//getchar();

	// Запоминаем скорректированную скорость с предыдущей итерации.
	// 13 08 2015 вынесено в глобальную область видимости см.
	// my_unsteady_temperature.c файл.
	//doublereal **SpeedCorOld=new doublereal*[3];
	//for (i=0; i<3; i++) {
	//	SpeedCorOld[i]=new doublereal[f.maxelm+f.maxbound];
	//}
	// На первой итерации и так должен быть ноль, не зачем 
	// выделять этот случай особо.
	//if (bfirst_start) {
		//for (i=VX; i<=VZ; i++) {
		  //  for (integer j=0; j<f.maxelm+f.maxbound; j++) {
			//    SpeedCorOld[i][j]=0.0;
		    //}
		//}
	//}
	//else {
#pragma omp parallel for
    for (integer j=0; j<f.maxelm+f.maxbound; j++) {
	    SpeedCorOld[VELOCITY_X_COMPONENT][j]=f.potent[VXCOR][j];
	    SpeedCorOld[VELOCITY_Y_COMPONENT][j]=f.potent[VYCOR][j];
	    SpeedCorOld[VELOCITY_Z_COMPONENT][j]=f.potent[VZCOR][j];
    }
	
	//}

	// коррекция скорости:
	// Только для внутренних КО.
	/*
	for (i=0; i<f.maxelm; i++) {
		correct_internal_volume(i, VX, f.slau, f.nvtx, f.potent, f.maxelm, f.alpha, f.pa, f.neighbors_for_the_internal_node, inumiter);
		correct_internal_volume(i, VY, f.slau, f.nvtx, f.potent, f.maxelm, f.alpha, f.pa, f.neighbors_for_the_internal_node, inumiter);
		correct_internal_volume(i, VZ, f.slau, f.nvtx, f.potent, f.maxelm, f.alpha, f.pa, f.neighbors_for_the_internal_node, inumiter);
        
		// ВНИМАНИЕ ВАЖНО !!!
		// Т.к. граничные узлы не подвергаются коррекции это может приводить к расходимости.
		// Граничный узел либо автономен - стоит условие Дирихле и коррекции не требуется, либо
		// на границе стоит однородное условие Неймана и скорость на границе равна скорости в 
		// ближайшем ненулевом КО.
		// Если будут проблемы со сходимостью имеет смысл снести скорректированную скорость, 
		// удовлетворяющую уравнению неразрывности из внутреннего узла в граничный.


	}
	*/
	
	
	// коррекция скорости на основе предварительно вычисленного по теореме Грина-Гаусса градиента поправки давления.
	/*
	for (i=0; i<f.maxelm; i++) {
		correct_internal_volume2(i, VX, f.slau, f.nvtx, f.potent, f.maxelm, f.alpha, f.pa, f.neighbors_for_the_internal_node, inumiter);
		correct_internal_volume2(i, VY, f.slau, f.nvtx, f.potent, f.maxelm, f.alpha, f.pa, f.neighbors_for_the_internal_node, inumiter);
		correct_internal_volume2(i, VZ, f.slau, f.nvtx, f.potent, f.maxelm, f.alpha, f.pa, f.neighbors_for_the_internal_node, inumiter);
	}
	*/

	// коррекция скорости на основе предварительно вычисленного по теореме Грина-Гаусса градиента поправки давления.
	//  с применением сглаженного псевдовремени !!!.
	// На существенно неравномерных сетках псевдовремя обладает сильной пространственной неоднородностью,
	// поэтому особенно на первых итерациях в целях подавления неустойчивостей имеет смысл взять всюду постоянный
	// в пространстве шаг по псевдоремени. Так рекомендует Гаврилов Андрей.
	/*
	for (i=0; i<f.maxelm; i++) {
	    correct_internal_volume3(i, VX, f.prop, f.potent,  tau[0]);
		correct_internal_volume3(i, VY, f.prop, f.potent,  tau[0]);
		correct_internal_volume3(i, VZ, f.prop, f.potent,  tau[0]);
	}*/

	// Основано на трёх скалярных сглаженных полях псевдовремени.
#pragma omp parallel for
	for (integer i=0; i<f.maxelm; i++) {
	    correct_internal_volume4(i, VELOCITY_X_COMPONENT, f.prop, f.potent,  tau);
		correct_internal_volume4(i, VELOCITY_Y_COMPONENT, f.prop, f.potent,  tau);
		correct_internal_volume4(i, VELOCITY_Z_COMPONENT, f.prop, f.potent,  tau);
	}
	


	// ВНИМАНИЕ ВАЖНО !!!
	// Т.к. граничные узлы не подвергаются коррекции это может приводить к расходимости.
	// Граничный узел либо автономен - стоит условие Дирихле и коррекции не требуется, либо
	// на границе стоит однородное условие Неймана и скорость на границе равна скорости в 
	// ближайшем ненулевом КО.
	// Если будут проблемы со сходимостью имеет смысл снести скорректированную скорость, 
	// удовлетворяющую уравнению неразрывности из внутреннего узла в граничный.

	// коррекция граничных скоростей только для границы на которой задано давление:
	// Также будем релаксировать к скорректированной скорости на предыдущей итерации 28.07.2016.
	const integer binterpol = 0; // см. correct_velocity.cpp
	correct_boundary_volume(VELOCITY_X_COMPONENT, f.potent, f.maxelm, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.border_neighbor, ls, lw, w, SpeedCorOld[VELOCITY_X_COMPONENT], binterpol);
	correct_boundary_volume(VELOCITY_Y_COMPONENT, f.potent, f.maxelm, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.border_neighbor, ls, lw, w, SpeedCorOld[VELOCITY_Y_COMPONENT], binterpol);
	correct_boundary_volume(VELOCITY_Z_COMPONENT, f.potent, f.maxelm, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.border_neighbor, ls, lw, w, SpeedCorOld[VELOCITY_Z_COMPONENT], binterpol);

#pragma omp parallel for
    for (integer i=0; i<(f.maxelm+f.maxbound); i++) {
		// запоминаем поле скорости 
		// удовлетворяющее уравнению
		// неразрывности, затем, чтобы
		// потом осуществить к нему нижнюю 
		// релаксацию.
		f.potent[VXCOR][i]=f.potent[VELOCITY_X_COMPONENT][i];
        f.potent[VYCOR][i]=f.potent[VELOCITY_Y_COMPONENT][i];
        f.potent[VZCOR][i]=f.potent[VELOCITY_Z_COMPONENT][i];
	}
	
	

	// Проверка корректности коррекции в граничных узлах при условии прилипания.
    iscorrectOk(f.potent, f.maxelm, f.neighbors_for_the_internal_node, f.border_neighbor, ls, lw, w);

	//printf("SINVARIANT end\n");
	//getchar();
	// Внимание его нужно вычислить обязательно до метода Грина Гаусса, т.к. в методе Гринна-Гаусса
	// используется интерполяция на грань, а значит требуется уже вычисленная поправка Рхи-Чоу.
	if (0) {

		// изменено 22 июня 2012 года.

		// не коррекция массового потока, а ВЫЧИСЛЕНИЕ на основе скорректированных
		// полей скорости и давления.
		// Равносильны ли такое ВЫЧИСЛЕНИЕ и коррекция масового потока остаётся открытым вопросом ?
		
		for (integer iP=0; iP<f.maxelm; iP++) {

			// вычисляем скорректированный массовый поток через грани КО.
			// Массовый поток вычисляется по обычным формулам с учётом
			// монотонизирующей поправки Рхи-Чоу но при его вычислении используются
			// скорректированные скорости и скорректированное давление.
			if (btimedep)
			{
            /*
			return_calc_correct_mass_flux(iP, 
									 f.potent,
									 f.pa,
									 f.prop,
									 f.prop_b,
					                 f.nvtx,
									 f.neighbors_for_the_internal_node,
									 f.maxelm,
									 f.diag_coef,
									 f.alpha,
						             RCh, // 1.0; 0.1;
									 btimedep,
						             dtimestep,
									 mfoldtimestep[iP],
						             f.mf[iP],  // возвращаемое значение массового потока
									 speedoldtimestep,false,
									 SpeedCorOld,
									 mfold[iP]);
									 */

		  // Данный вариант основывается на сглаженном псевдовремени &&
		  // рекомендован к использованию (см. опыт Гаврилова Андрея sigmaflow).
          /*return_calc_correct_mass_flux2(iP,
			                             f.potent,
										 f.pa, 
										 f.prop,
										 f.prop_b,
					                     f.nvtx,
										 f.neighbors_for_the_internal_node,
										 f.maxelm,
										 f.alpha,
										 RCh, // 1.0; 0.1;
						                 btimedep,
										 dtimestep, 
										 mfoldtimestep[iP],
						                 f.mf[iP], // возвращаемое значение массового потока
										 speedoldtimestep, false,
						                 SpeedCorOld,
										 mfold[iP],
										 tau);*/

			// учтена тройственность скалярного псевдовремени.
		   return_calc_correct_mass_flux3(iP,
			                             f.potent,
										 f.pa, 
										 f.prop,
										 f.prop_b,
					                     f.nvtx,
										 f.neighbors_for_the_internal_node,
										 f.maxelm,
										 f.alpha,
										 RCh, // 1.0; 0.1;
						                 btimedep,
										 dtimestep, 
										 mfoldtimestep[iP],
						                 f.mf[iP], // возвращаемое значение массового потока
										 speedoldtimestep, false,
						                 SpeedCorOld,
										 mfold[iP],
										 tau,
			                             f.border_neighbor, 
			                             t.ilevel_alice, 
			                             f.ptr);

			}
			else {
				// стационарный случай.
				/*
            return_calc_correct_mass_flux(iP, 
									 f.potent,
									 f.pa,
									 f.prop,
									 f.prop_b,
					                 f.nvtx,
									 f.neighbors_for_the_internal_node,
									 f.maxelm,
									 f.diag_coef,
									 f.alpha,
						             RCh, // 1.0; 0.1;
									 btimedep,
						             dtimestep,
									 nullptr,
						             f.mf[iP], // возвращаемое значение массового потока
									 nullptr,false,
									 SpeedCorOld,
									 mfold[iP]);
									 */

			// Данный вариант основывается на сглаженном псевдовремени &&
		    // рекомендован к использованию (см. опыт Гаврилова Андрея sigmaflow).
			// Так как мы пока не знаем какой шаг по времени выбрать в стационарном случае,
			// то зададим отрицательное значение, это будет сигналом использовать среднее арифметическое
			// среди предварительно вычисленных псевдовременных шагшов и сглаживание будет производится 
			// к этому среднему арифметическому значению.
           /* return_calc_correct_mass_flux2(iP,
			                             f.potent,
										 f.pa, 
										 f.prop,
										 f.prop_b,
					                     f.nvtx,
										 f.neighbors_for_the_internal_node,
										 f.maxelm,
										 f.alpha,
										 RCh, // 1.0; 0.1;
						                 btimedep,
										 -0.1,   //dtimestep, 
										 nullptr,
						                 f.mf[iP], // возвращаемое значение массового потока
										 nullptr, false,
						                 SpeedCorOld,
										 mfold[iP],
										 tau);*/

            // учтена тройственность скалярного псевдовремени.
			return_calc_correct_mass_flux3(iP,
			                             f.potent,
										 f.pa, 
										 f.prop,
										 f.prop_b,
					                     f.nvtx,
										 f.neighbors_for_the_internal_node,
										 f.maxelm,
										 f.alpha,
										 RCh, // 1.0; 0.1;
						                 btimedep,
										 -0.1,   //dtimestep, 
										 nullptr,
						                 f.mf[iP], // возвращаемое значение массового потока
										 nullptr, false,
						                 SpeedCorOld,
										 mfold[iP],
										 tau,
				                         f.border_neighbor,
				                         t.ilevel_alice,
				                         f.ptr);

			}
						         

		}
	}

	// коррекция массового потока по формуле схожей с формулой коррекции скорости.
	if (1) {
    	correct_mf(f.mf, f.potent,  tau,
	               f.pa, f.neighbors_for_the_internal_node, f.nvtx, f.maxelm,
				   f.border_neighbor, ls, lw, w, f.prop_b, t.ilevel_alice, f.ptr);
	}

	// Коррекция скорости и массового потока на выходной грнице
	// после коррекции скоростей см. Гаврилов Андрей.
	//mass_balance(f, lw, ls, w);

	
	// Освобождение оперативной памяти из кучи.
	// 13 августа 2015 вынесено в глобальную область видимости для ускорения быстродействия.
	//for (i=0; i<3; i++) {
		//delete SpeedCorOld[i];
	//}
	//delete SpeedCorOld; 
#pragma omp parallel for
	for (integer i=0; i<f.maxelm; i++) {
		for (integer j=0; j<6; j++) {
			// именно этот поток должен использоваться при нижней релаксации
			// поправки Рхи-Чоу предложенной в статье I.Sezai.
			mfold[i][j]=f.mf[i][j];
			if (mfold[i][j]!= mfold[i][j]) {
				printf("error nan or inf in function my_version_SIMPLE_Algorithm3D\n");
				printf("f.mf[%lld][%lld]=%e\n",i, j, f.mf[i][j]);
				system("pause");
			}
		}
	}

	// Обязательная проверка корректности выполнения условий прилипания.
	iscorrectmf(f.mf, f.maxelm, f.neighbors_for_the_internal_node, f.border_neighbor, ls, lw, w);
	
	//printf("MASS FLUX end\n");
	//getchar();
	/*
	for (i=0; i<f.maxelm; i++) {
		f.potent[FBUF][i]=tau[VY][i];
	}
	for (i=0; i<f.maxbound; i++) {
		f.potent[FBUF][i+f.maxelm]=tau[VY][i];
	}
	*/


	// экспорт результата вычисления в программу tecplot360:
	if (0) {
		//if (inumiter>82) {
		if (!b_on_adaptive_local_refinement_mesh) {
		   exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior,iflow,bextendedprint,0, b, lb);
	    }
	    else {
		   // Экспорт в программу tecplot температуры.
		   //С АЛИС сетки.
		   ANES_tecplot360_export_temperature(t.maxnod, t.pa, t.maxelm, t.nvtx, t.potent, t, fglobal, 0, b, lb);
	     }
	       printf("corect values. OK.\n");
	       //getchar(); // debug avtosave
		   system("pause");
		//}
	}
	if (0) {
		xyplot( fglobal, flow_interior, t);
		printf("xy plot corect values. OK.\n");
	    //getchar(); // debug avtosave
		system("pause");
	}

	//printf("CORRECT end\n");
	//getchar();

#ifdef D_OPENMP

	if (bparallelizm_old) {
			if (inumcore==2) {
				if (nd.b0.active) {

					#pragma omp parallel 
					{
#pragma omp sections 
						{
#pragma omp section
							{

					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							//green_gauss(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f);
							green_gaussO1(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.mf[iPloc], f.prop[RHO], f.prop_b[RHO], f.border_neighbor, t.ilevel_alice);
						}

					}
							}
#pragma omp section
							{
					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							//green_gauss(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f);
							green_gaussO1(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.mf[iPloc], f.prop[RHO], f.prop_b[RHO], f.border_neighbor, t.ilevel_alice);
						}

					}
							} // section
						} // sections
					} // parallel
					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							//green_gauss(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f);
							green_gaussO1(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f.mf[iPloc], f.prop[RHO], f.prop_b[RHO], f.border_neighbor, t.ilevel_alice);
						}
					}


				}
			}

			if (inumcore==2) {
				if (nd.b0.active) {

					#pragma omp parallel 
					{
#pragma omp sections 
						{
#pragma omp section
							{

					// первый поток
					for (integer iscan_par=nd.b0.ileft_start; iscan_par<=nd.b0.ileft_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							//green_gauss(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f);
						    green_gaussO1(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.mf[iPloc], f.prop[RHO], f.prop_b[RHO], f.border_neighbor, t.ilevel_alice);
						}

					}
							}
#pragma omp section
							{
					// второй поток
					for (integer iscan_par=nd.b0.iright_start; iscan_par<=nd.b0.iright_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							//green_gauss(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f);
					     	green_gaussO1(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.mf[iPloc], f.prop[RHO], f.prop_b[RHO], f.border_neighbor, t.ilevel_alice);
						}

					}
							} // section
						} // sections
					} // parallel
					// серийный смыкающий кусок
					for (integer iscan_par=nd.b0.iseparate_start; iscan_par<=nd.b0.iseparate_finish; iscan_par++) {
						integer iPloc=ifrontregulationgl[iscan_par];
						if (iPloc<f.maxelm) {
							//green_gauss(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f);
						    green_gaussO1(iPloc, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f.mf[iPloc], f.prop[RHO], f.prop_b[RHO], f.border_neighbor, t.ilevel_alice);
						}

					}


				}
			}
			}
#else

	
		// Вычисление градиентов скорости:
		// на основе поля скорости удовлетворяющего уравнению неразрывности.
#pragma omp parallel for
		for (integer i = 0; i < f.maxelm; i++) {
			// градиенты скоростей для внутренних КО.
			green_gauss(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, false, f, f.border_neighbor, t.ilevel_alice);
		}
#pragma omp parallel for
		for (integer i = 0; i < f.maxelm; i++) {
			// градиенты скоростей для граничных КО.
			green_gauss(i, f.potent, f.nvtx, f.pa, f.neighbors_for_the_internal_node, f.maxelm, true, f, f.border_neighbor, t.ilevel_alice);
		}
		
#endif

	
	
		// На твёрдой стенке турбулентная динамическая вязкость равна нулю.
		// Вычисление S инварианта тензора скоростей-деформаций для всех
		// внутренних и граничных контрольных объёмов расчётной области.
#pragma omp parallel for shared (f)  schedule (guided)
		for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
			// по поводу правильности формулы см. user_manual.
			doublereal sum = 0.0;
			sum += 2.0*f.potent[GRADXVX][i] * f.potent[GRADXVX][i];
			sum += 2.0*f.potent[GRADYVY][i] * f.potent[GRADYVY][i];
			sum += 2.0*f.potent[GRADZVZ][i] * f.potent[GRADZVZ][i];
			sum += (f.potent[GRADYVX][i] + f.potent[GRADXVY][i])*(f.potent[GRADYVX][i] + f.potent[GRADXVY][i]);
			sum += (f.potent[GRADZVX][i] + f.potent[GRADXVZ][i])*(f.potent[GRADZVX][i] + f.potent[GRADXVZ][i]);
			sum += (f.potent[GRADYVZ][i] + f.potent[GRADZVY][i])*(f.potent[GRADYVZ][i] + f.potent[GRADZVY][i]);
			// следующее слагаемое может сделать подкоренное выражение отрицательным.
			// вычитаем две трети квадрата дивергенции.
			sum -= (2.0 / 3.0)*(f.potent[GRADXVX][i] + f.potent[GRADYVY][i] + f.potent[GRADZVZ][i])*(f.potent[GRADXVX][i] + f.potent[GRADYVY][i] + f.potent[GRADZVZ][i]); // добавок связанный с несжимаемостью/сжимаемостью
			f.SInvariantStrainRateTensor[i] = sqrt(fmax(0.0, sum));
			if (f.SInvariantStrainRateTensor[i]!= f.SInvariantStrainRateTensor[i]) {
				printf("Error nan or inf in function my_version_SIMPLE_Algorithm3D\n");
				printf("f.SInvariantStrainRateTensor[%lld]=%e\n", i, f.SInvariantStrainRateTensor[i]);
				system("pause");
			}

			// Вихрь (модуль ротора скорости).
			sum = 0.0;
			sum += (f.potent[GRADYVZ][i] - f.potent[GRADZVY][i])*(f.potent[GRADYVZ][i] - f.potent[GRADZVY][i]); // проверено.
			sum += (f.potent[GRADZVX][i] - f.potent[GRADXVZ][i])*(f.potent[GRADZVX][i] - f.potent[GRADXVZ][i]);
			sum += (f.potent[GRADXVY][i] - f.potent[GRADYVX][i])*(f.potent[GRADXVY][i] - f.potent[GRADYVX][i]);
			f.potent[CURL][i] = sqrt(sum);
			if (f.potent[CURL][i]!= f.potent[CURL][i]) {
				printf("Error nan or inf in function my_version_SIMPLE_Algorithm3D\n");
				printf("f.potent[CURL][%lld]=%e\n", i, f.potent[CURL][i]);
				system("pause");
			}
		}
		
		

		// RANS Спалларт Аллмарес [1992].
		if (f.iflowregime == VISCOSITY_MODEL::RANS_SPALART_ALLMARES) {


			doublereal *nusha_aprior = new doublereal[f.maxelm+f.maxbound];

#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm+ f.maxbound); i++) {
				nusha_aprior[i] = f.potent[NUSHA][i];
			}

			if (bflag_free_memory_cfd) {
				m.bsignalfreeCRScfd = true; // Освобождаем память.
			}

			//printf("NUSHA \n");
			//getchar();
			solve(NUSHA, res, f, fglobal, t, rhie_chow,
				s, w, b, ls, lw, lb, bonbeta,
				flow_interior, false,
				bfirst_start, toldtimestep, nullptr,
				speedoldtimestep, mfoldtimestep,
				dtimestep, btimedep, dgx, dgy, dgz,
				matlist, inumiter, bprintmessage,
				RCh, bVERYStable, nullptr, sumanb, false, false, 1.0, m,
				rthdsd, rfluentres.res_nusha, lu, my_union, color, dist_max_fluid);
			

			doublereal fHORF = 0.25; // for steady state problem.
			if (btimedep) { // unsteady problems.
				fHORF = 0.75; // ANSYS Fluent Theory Guide.
			}
#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				f.potent[NUSHA][i] = nusha_aprior[i] +
					fHORF * (f.potent[NUSHA][i] - nusha_aprior[i]);
			}

			delete[] nusha_aprior; // Освобождение оперативной памяти.

			// именно здесь верно. 04.05.2017
			rfluentres.res_nusha = fluent_residual_for_x(f.slau[NUSHA_SL], f.slau_bon[NUSHA_SL], f.potent[NUSHA], f.maxelm, f.maxbound, NUSHA); // невязка по формуле fluent.


			// Вычисление градиента модифицированной кинетической турбулентой вязкости
#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_SpallartAllmares(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, false,
					f.border_neighbor, t.ilevel_alice);
			}

#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_SpallartAllmares(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, true,
					f.border_neighbor, t.ilevel_alice);
			}

			// Вычисление турбулентной вязкости.
			// 30.09.2019
			for (integer i = 0; i<(f.maxelm + f.maxbound); i++) {
				doublereal rho = 0.0, mu=0.0; // Вычисление плотности и динамической молекулярной вязкости.
				if (i < f.maxelm) {
					rho = f.prop[RHO][i];
					mu= f.prop[MU_DYNAMIC_VISCOSITY][i];
				}
				else {
					rho = f.prop_b[RHO][i - f.maxelm];
                    mu= f.prop_b[MU_DYNAMIC_VISCOSITY][i - f.maxelm];
				}
				//doublereal PrandtlLength = fmin(0.419*f.rdistWall[i], 0.09*f.rdistWallmax); // Формула Эскудиера (1966).
				//f.potent[MUT][i]=rho*PrandtlLength*PrandtlLength*f.SInvariantStrainRateTensor[i];
				
				doublereal kappa = rho*f.potent[NUSHA][i] / mu;
				if (kappa < 1.0e-4) kappa = 1.0e-4; // Емельянов, Волков, Дерюгин.
				doublereal f_nu1 = (kappa*kappa*kappa) / (kappa*kappa*kappa+eqin.fluidinfo[0].c_nu1*eqin.fluidinfo[0].c_nu1*eqin.fluidinfo[0].c_nu1);
				f.potent[MUT][i] = rho*f_nu1*f.potent[NUSHA][i];
				
				//if (0 && (inumiter == 40)) {
					//printf("rho=%e dw=%e dwmax=%e plen=%e sigma=%e\n", rho, f.rdistWall[i], f.rdistWallmax, PrandtlLength, f.SInvariantStrainRateTensor[i]);
					//system("pause");
				//}
			}

		}

		//if (inumiter == 8) system("PAUSE");

		// RANS SST Ментер [1993]
		if ((inumiter>1)&&(f.iflowregime == VISCOSITY_MODEL::RANS_MENTER_SST)) {

			// inumiter>1 - скорости должны немного устаканиться, сразу запускать расчёт турбулентности
			// на неустановившейся скорости нельзя, это приведет к расходимости итерационного процесса.

			doublereal *turb_kinetik_energy_aprior = new doublereal[f.maxelm + f.maxbound];

#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				turb_kinetik_energy_aprior[i] = f.potent[TURBULENT_KINETIK_ENERGY][i];
			}

			//printf("TURBULENT_KINETIK_ENERGY \n");
			//getchar();
			solve(TURBULENT_KINETIK_ENERGY, res, f, fglobal, t, rhie_chow,
				s, w, b, ls, lw, lb, bonbeta,
				flow_interior, false,
				bfirst_start, toldtimestep, nullptr,
				speedoldtimestep, mfoldtimestep,
				dtimestep, btimedep, dgx, dgy, dgz,
				matlist, inumiter, bprintmessage,
				RCh, bVERYStable, nullptr, sumanb, false, false, 1.0, m,
				rthdsd, rfluentres.res_turb_kinetik_energy, lu, my_union, color, dist_max_fluid);
			

			doublereal fHORF = 0.25; // for steady state problem.
			if (btimedep) { // unsteady problems.
				fHORF = 0.75; // ANSYS Fluent Theory Guide.
			}
#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				f.potent[TURBULENT_KINETIK_ENERGY][i] = fmax(K_limiter_min, turb_kinetik_energy_aprior[i] +
					fHORF * (f.potent[TURBULENT_KINETIK_ENERGY][i] - turb_kinetik_energy_aprior[i]));
			}

			delete[] turb_kinetik_energy_aprior; // Освобождение оперативной памяти.

			// именно здесь верно. 04.05.2017
			rfluentres.res_turb_kinetik_energy = fluent_residual_for_x(f.slau[TURBULENT_KINETIK_ENERGY_SL], f.slau_bon[TURBULENT_KINETIK_ENERGY_SL], f.potent[TURBULENT_KINETIK_ENERGY], f.maxelm, f.maxbound, TURBULENT_KINETIK_ENERGY); // невязка по формуле fluent.



			// Вычисление градиента кинетической энергии турбулентных пульсаций
#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_turbulent_kinetik_energy_MenterSST(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, false,
					f.border_neighbor, t.ilevel_alice);
			}

#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_turbulent_kinetik_energy_MenterSST(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, true,
					f.border_neighbor, t.ilevel_alice);
			}

			doublereal *omega_aprior = new doublereal[f.maxelm + f.maxbound];

#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				omega_aprior[i] = f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][i];
			}

			if (bflag_free_memory_cfd) {
				m.bsignalfreeCRScfd = true; // Освобождаем память.
			}

			//printf("TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA \n");
			//getchar();
			solve(TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA, res, f, fglobal, t, rhie_chow,
				s, w, b, ls, lw, lb, bonbeta,
				flow_interior, false,
				bfirst_start, toldtimestep, nullptr,
				speedoldtimestep, mfoldtimestep,
				dtimestep, btimedep, dgx, dgy, dgz,
				matlist, inumiter, bprintmessage,
				RCh, bVERYStable, nullptr, sumanb, false, false, 1.0, m,
				rthdsd, rfluentres.res_turb_omega, lu, my_union, color, dist_max_fluid);
			

#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][i] = fmax(Omega_limiter_min, omega_aprior[i] +
					fHORF * (f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][i] - omega_aprior[i]));
			}

			delete[] omega_aprior; // Освобождение оперативной памяти.

			// именно здесь верно. 04.05.2017
			rfluentres.res_turb_omega = fluent_residual_for_x(f.slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], f.slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL], f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA], f.maxelm, f.maxbound, TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA); // невязка по формуле fluent.


			// Вычисление градиента удельной скорости диссипации 
			// кинетической энергии турбулентных пульсаций.
#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_specific_dissipation_rate_omega_MenterSST(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, false,
					f.border_neighbor, t.ilevel_alice);
			}

#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_specific_dissipation_rate_omega_MenterSST(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, true,
					f.border_neighbor, t.ilevel_alice);
			}

			doublereal* Tx = nullptr;
			doublereal* Ty = nullptr;
			doublereal* Tz = nullptr;
			if (bSIMPLErun_now_for_natural_convection) 
			{
				Tx = new doublereal[t.maxelm + t.maxbound];
				Ty = new doublereal[t.maxelm + t.maxbound];
				Tz = new doublereal[t.maxelm + t.maxbound];

				// инициализация нулём.
				for (integer i = 0; i < t.maxelm + t.maxbound; i++) {
					Tx[i] = 0.0;
					Ty[i] = 0.0;
					Tz[i] = 0.0;
				}

				// нахождение градиентов.
#pragma omp parallel for
				for (integer i = 0; i < t.maxelm; i++) {
					// Только внутренние узлы.
					green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
						t.neighbors_for_the_internal_node, t.maxelm, false,
						t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
				}

#pragma omp parallel for
				for (integer i = 0; i < t.maxelm; i++) {
					// Только граничные узлы.
					green_gaussTemperature(i, t.potent, t.nvtx, t.pa,
						t.neighbors_for_the_internal_node, t.maxelm, true,
						t.border_neighbor, Tx, Ty, Tz, t.ilevel_alice);
				}
			}

			// Вычисление турбулентной вязкости.
			// 30.09.2019
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				doublereal rho = 0.0, mu = 0.0, beta=0.0; // Вычисление плотности и динамической молекулярной вязкости.
				if (i < f.maxelm) {
					rho = f.prop[RHO][i];
					mu = f.prop[MU_DYNAMIC_VISCOSITY][i];
					beta = t.prop[BETA_T][f.ptr[i]];
				}
				else {
					rho = f.prop_b[RHO][i - f.maxelm];
					mu = f.prop_b[MU_DYNAMIC_VISCOSITY][i - f.maxelm];
				}

				// Турбулентная динамическая вязкость.
				doublereal F2 = 1.0;
				
				doublereal arg2 = fmax(2.0*sqrt(fmax(K_limiter_min,f.potent[TURBULENT_KINETIK_ENERGY][i])) 
					/ (eqin.fluidinfo[0].beta_zvezda*fmax(Omega_limiter_min,f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][i]) *
						f.rdistWall[i]), (500.0*mu / rho) / (f.rdistWall[i] * f.rdistWall[i] *
							fmax(Omega_limiter_min,f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][i])));
				F2 = tanh(arg2*arg2);
				//doublereal m2 = f.potent[CURL][i] * F2;
				doublereal m2 = f.SInvariantStrainRateTensor[i] * F2;
				doublereal m1 = fmax(eqin.fluidinfo[0].menter_a1*fmax(Omega_limiter_min,f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][i]),m2);
				f.potent[MUT][i] = rho * eqin.fluidinfo[0].menter_a1*fmax(K_limiter_min,f.potent[TURBULENT_KINETIK_ENERGY][i]) / m1;
				f.potent[MUT][i] = fmax(0.0, f.potent[MUT][i]);

				// Коррекция турбулентной вязкости 
				// связанная с градиентным числом Ричардсона
				if (bSIMPLErun_now_for_natural_convection) {
					if (i < f.maxelm) {
						// Коррекция SST модели 2003 года.
						// Подход Роди-Лешцинера-Исаева
						// Прямая коррекция турбулентной вязкости.
						doublereal cr = 0.02 * beta * ((dgz * Tz[f.ptr[i]]) / (fmax(1.0e-30, (f.potent[GRADZVX][i]) * (f.potent[GRADZVX][i]) +
							(f.potent[GRADZVY][i]) * (f.potent[GRADZVY][i]) + (f.potent[GRADZVZ][i]) * (f.potent[GRADZVZ][i]))) +
							(dgy * Ty[f.ptr[i]]) / (fmax(1.0e-30, (f.potent[GRADYVX][i]) * (f.potent[GRADYVX][i]) +
							(f.potent[GRADYVY][i]) * (f.potent[GRADYVY][i]) + (f.potent[GRADYVZ][i]) * (f.potent[GRADYVZ][i]))) +
							(dgx * Tx[f.ptr[i]]) / (fmax(1.0e-30, (f.potent[GRADXVX][i]) * (f.potent[GRADXVX][i]) +
							(f.potent[GRADXVY][i]) * (f.potent[GRADXVY][i]) + (f.potent[GRADXVZ][i]) * (f.potent[GRADXVZ][i]))));
						cr /= 1.0 - beta * t.potent[f.ptr[i]];
						cr = fmax(-0.4, fmin(cr, 3.5));
						f.potent[MUT][i] *= 1.0 / (1.0 + cr);
					}
				}
				
			}

			if (bSIMPLErun_now_for_natural_convection) {
				delete[]  Tx;
				delete[]  Ty;
				delete[]  Tz;
			}

			// Принудительно сделаем чтобы граничные значения совпадали с ближайшими внутренними.
		// Это нужно прежде всего для корректности задания граничных условий.
			for (integer iP = 0; iP < f.maxelm; iP++) {
				integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
				iE = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE1; iN = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE1; iT = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE1;
				iW = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE1; iS = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE1; iB = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE1;
				if (iE > -1) {
					if (iE >= f.maxelm) {
						f.potent[MUT][iE] = f.potent[MUT][iP];
					}
				}
				if (iW > -1) {
					if (iW >= f.maxelm) {
						f.potent[MUT][iW] = f.potent[MUT][iP];
					}
				}
				if (iN > -1) {
					if (iN >= f.maxelm) {
						f.potent[MUT][iN] = f.potent[MUT][iP];
					}
				}
				if (iS > -1) {
					if (iS >= f.maxelm) {
						f.potent[MUT][iS] = f.potent[MUT][iP];
					}
				}
				if (iT > -1) {
					if (iT >= f.maxelm) {
						f.potent[MUT][iT] = f.potent[MUT][iP];
					}
				}
				if (iB > -1) {
					if (iB >= f.maxelm) {
						f.potent[MUT][iB] = f.potent[MUT][iP];
					}
				}

				integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
				iE2 = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE2; iN2 = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE2; iT2 = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE2;
				iW2 = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE2; iS2 = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE2; iB2 = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE2;
				if (iE2 > -1) {
					if (iE2 >= f.maxelm) {
						f.potent[MUT][iE2] = f.potent[MUT][iP];
					}
				}
				if (iW2 > -1) {
					if (iW2 >= f.maxelm) {
						f.potent[MUT][iW2] = f.potent[MUT][iP];
					}
				}
				if (iN2 > -1) {
					if (iN2 >= f.maxelm) {
						f.potent[MUT][iN2] = f.potent[MUT][iP];
					}
				}
				if (iS2 > -1) {
					if (iS2 >= f.maxelm) {
						f.potent[MUT][iS2] = f.potent[MUT][iP];
					}
				}
				if (iT2 > -1) {
					if (iT2 >= f.maxelm) {
						f.potent[MUT][iT2] = f.potent[MUT][iP];
					}
				}
				if (iB2 > -1) {
					if (iB2 >= f.maxelm) {
						f.potent[MUT][iB2] = f.potent[MUT][iP];
					}
				}

				integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
				iE3 = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE3; iN3 = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE3; iT3 = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE3;
				iW3 = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE3; iS3 = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE3; iB3 = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE3;
				if (iE3 > -1) {
					if (iE3 >= f.maxelm) {
						f.potent[MUT][iE3] = f.potent[MUT][iP];
					}
				}
				if (iW3 > -1) {
					if (iW3 >= f.maxelm) {
						f.potent[MUT][iW3] = f.potent[MUT][iP];
					}
				}
				if (iN3 > -1) {
					if (iN3 >= f.maxelm) {
						f.potent[MUT][iN3] = f.potent[MUT][iP];
					}
				}
				if (iS3 > -1) {
					if (iS3 >= f.maxelm) {
						f.potent[MUT][iS3] = f.potent[MUT][iP];
					}
				}
				if (iT3 > -1) {
					if (iT3 >= f.maxelm) {
						f.potent[MUT][iT3] = f.potent[MUT][iP];
					}
				}
				if (iB3 > -1) {
					if (iB3 >= f.maxelm) {
						f.potent[MUT][iB3] = f.potent[MUT][iP];
					}
				}

				integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов
				iE4 = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE4; iN4 = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE4; iT4 = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE4;
				iW4 = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE4; iS4 = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE4; iB4 = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE4;
				if (iE4 > -1) {
					if (iE4 >= f.maxelm) {
						f.potent[MUT][iE4] = f.potent[MUT][iP];
					}
				}
				if (iW4 > -1) {
					if (iW4 >= f.maxelm) {
						f.potent[MUT][iW4] = f.potent[MUT][iP];
					}
				}
				if (iN4 > -1) {
					if (iN4 >= f.maxelm) {
						f.potent[MUT][iN4] = f.potent[MUT][iP];
					}
				}
				if (iS4 > -1) {
					if (iS4 >= f.maxelm) {
						f.potent[MUT][iS4] = f.potent[MUT][iP];
					}
				}
				if (iT4 > -1) {
					if (iT4 >= f.maxelm) {
						f.potent[MUT][iT4] = f.potent[MUT][iP];
					}
				}
				if (iB4 > -1) {
					if (iB4 >= f.maxelm) {
						f.potent[MUT][iB4] = f.potent[MUT][iP];
					}
				}
			}


		} // Menter SST [1993]


		  // RANS двухслойная модель с двумя уравнениями [2001]
		if (f.iflowregime == VISCOSITY_MODEL::RANS_STANDART_K_EPS) {

			doublereal *turb_kinetik_energy_aprior = new doublereal[f.maxelm + f.maxbound];

#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				turb_kinetik_energy_aprior[i] = f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i];
			}

			//printf("TURBULENT_KINETIK_ENERGY_STD_K_EPS \n");
			//getchar();
			solve(TURBULENT_KINETIK_ENERGY_STD_K_EPS, res, f, fglobal, t, rhie_chow,
				s, w, b, ls, lw, lb, bonbeta,
				flow_interior, false,
				bfirst_start, toldtimestep, nullptr,
				speedoldtimestep, mfoldtimestep,
				dtimestep, btimedep, dgx, dgy, dgz,
				matlist, inumiter, bprintmessage,
				RCh, bVERYStable, nullptr, sumanb, false, false, 1.0, m,
				rthdsd, rfluentres.res_turb_kinetik_energy_std_ke, lu, my_union, color, dist_max_fluid);
			

			doublereal fHORF = 0.25; // for steady state problem.
			if (btimedep) { // unsteady problems.
				fHORF = 0.75; // ANSYS Fluent Theory Guide.
			}
#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i] = fmax(K_limiter_min, turb_kinetik_energy_aprior[i] +
					fHORF * (f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i] - turb_kinetik_energy_aprior[i]));
			}

			delete[] turb_kinetik_energy_aprior; // Освобождение оперативной памяти.

			// именно здесь верно. 04.05.2017
			rfluentres.res_turb_kinetik_energy_std_ke =
				fluent_residual_for_x(f.slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
					f.slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL],
					f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS],
					f.maxelm, f.maxbound, TURBULENT_KINETIK_ENERGY_STD_K_EPS); // невязка по формуле fluent.


			// Вычисление градиента кинетической энергии турбулентных пульсаций
#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_turbulent_kinetik_energy_standart_k_epsilon(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, false,
					f.border_neighbor, t.ilevel_alice);
			}

#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_turbulent_kinetik_energy_standart_k_epsilon(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, true,
					f.border_neighbor, t.ilevel_alice);
			}

			doublereal *epsilon_aprior = new doublereal[f.maxelm + f.maxbound];

#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				epsilon_aprior[i] = f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][i];
			}

			if (bflag_free_memory_cfd) {
				m.bsignalfreeCRScfd = true; // Освобождаем память.
			}

			//printf("TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS \n");
			//getchar();
			solve(TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS,
				res, f, fglobal, t, rhie_chow,
				s, w, b, ls, lw, lb, bonbeta,
				flow_interior, false,
				bfirst_start, toldtimestep, nullptr,
				speedoldtimestep, mfoldtimestep,
				dtimestep, btimedep, dgx, dgy, dgz,
				matlist, inumiter, bprintmessage,
				RCh, bVERYStable, nullptr, sumanb, false, false, 1.0, m,
				rthdsd, rfluentres.res_turb_epsilon, lu, my_union, color, dist_max_fluid);
			

#pragma omp parallel for
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][i] = fmax(Epsilon_limiter_min, epsilon_aprior[i] +
					fHORF * (f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][i] - epsilon_aprior[i]));
			}

			delete[] epsilon_aprior; // Освобождение оперативной памяти.

			// именно здесь верно. 04.05.2017
			rfluentres.res_turb_epsilon = fluent_residual_for_x(
				f.slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
				f.slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL],
				f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS],
				f.maxelm, f.maxbound, TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS); // невязка по формуле fluent.

			// Вычисление градиента скорости диссипации 
#pragma omp parallel for																																																																						 // кинетической энергии турбулентных пульсаций.
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_turbulent_dissipation_rate_epsilon_standart_k_epsilon(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, false,
					f.border_neighbor, t.ilevel_alice);
			}

#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_turbulent_dissipation_rate_epsilon_standart_k_epsilon(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, true,
					f.border_neighbor, t.ilevel_alice);
			}

			// Для cross-diffusion term.
			// omega=epsilon/(C_mu*k);
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				f.potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][i] = fmax(Omega_limiter_min,
					fmax(Epsilon_limiter_min, f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][i])/
					(eqin.fluidinfo[0].beta_zvezda*
						fmax(K_limiter_min, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i])));
			}

			// Для cross-diffusion term.
			// Вычисление градиента удельной скорости диссипации 
			// кинетической энергии турбулентных пульсаций.
#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_specific_dissipation_rate_omega_MenterSST(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, false,
					f.border_neighbor, t.ilevel_alice);
			}

#pragma omp parallel for
			for (integer i = 0; i < f.maxelm; i++) {
				green_gauss_specific_dissipation_rate_omega_MenterSST(i,
					f.potent, f.nvtx, f.pa,
					f.neighbors_for_the_internal_node, f.maxelm, true,
					f.border_neighbor, t.ilevel_alice);
			}


			// Вычисление турбулентной вязкости.
			// 30.09.2019
			for (integer i = 0; i < (f.maxelm + f.maxbound); i++) {
				doublereal rho = 0.0, mu = 0.0; // Вычисление плотности и динамической молекулярной вязкости.
				if (i < f.maxelm) {
					rho = f.prop[RHO][i];
					mu = f.prop[MU_DYNAMIC_VISCOSITY][i];
				}
				else {
					rho = f.prop_b[RHO][i - f.maxelm];
					mu = f.prop_b[MU_DYNAMIC_VISCOSITY][i - f.maxelm];
				}

				if (fabs(f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][i])>1.0e-20) {

					doublereal speed_or_sqrt_k = sqrt(fmax(K_limiter_min, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i]));
					//speed_or_sqrt_k = sqrt(f.potent[VX][i] * f.potent[VX][i] + f.potent[VY][i] * f.potent[VY][i] + f.potent[VZ][i] * f.potent[VZ][i]);
					doublereal Re_y = speed_or_sqrt_k*
						f.rdistWall[i] * rho / mu;
					doublereal lambda_switch = eqin.fluidinfo[0].lambda_switch(Re_y);
					doublereal lnu = eqin.fluidinfo[0].Cl_std_ke*f.rdistWall[i]*(1.0-
						exp(-Re_y/ eqin.fluidinfo[0].Anu_std_ke));


					// Последняя формула на странице 76 из статьи
					// А.В. Кузьминов, В.Н. Лапин, С.Г. Черный
					// Метод расчёта турбулентных течений несжимаемой жидкости на 
					// основе двухслойной (k-epsilon)-модели. Новосибирск 2001.

					f.potent[MUT][i] = rho*((lambda_switch*(eqin.fluidinfo[0].C_mu_std_ke*
						fmax(K_limiter_min, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i]*
							fmax(K_limiter_min, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i])
						) /(fmax(Epsilon_limiter_min,f.potent[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS][i])
							)))+
						(1.0- lambda_switch)*eqin.fluidinfo[0].C_mu_std_ke*lnu*
						sqrt(fmax(K_limiter_min, f.potent[TURBULENT_KINETIK_ENERGY_STD_K_EPS][i])));

					// Если не ограничивать сверху то сразу расходимость вычислительного процесса.
					// 
					const doublereal Kturm_C_k = 0.005;
					doublereal speed2 = f.potent[VELOCITY_X_COMPONENT][i] * f.potent[VELOCITY_X_COMPONENT][i] +
						f.potent[VELOCITY_Y_COMPONENT][i] * f.potent[VELOCITY_Y_COMPONENT][i] +
						f.potent[VELOCITY_Z_COMPONENT][i] * f.potent[VELOCITY_Z_COMPONENT][i];
					
					doublereal viscosity_ratio_max = 6000; // 101.2
					viscosity_ratio_max = 9.156*rho*eqin.fluidinfo[0].C_mu_std_ke * lnu * sqrt(Kturm_C_k * speed2) / mu;
					f.potent[MUT][i] = fmin(viscosity_ratio_max*mu, fmax(0.0, f.potent[MUT][i]));
					
					//f.potent[MUT][i] = fmax(0.0, f.potent[MUT][i]);
					
				}
			}


			// Принудительно сделаем чтобы граничные значения совпадали с ближайшими внутренними.
			// Это нужно прежде всего для корректности задания граничных условий.
			for (integer iP = 0; iP < f.maxelm; iP++) {
				integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
				iE = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE1; iN = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE1; iT = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE1;
				iW = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE1; iS = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE1; iB = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE1;
				if (iE > -1) {
					if (iE >= f.maxelm) {
						f.potent[MUT][iE] = f.potent[MUT][iP];
					}
				}
				if (iW > -1) {
					if (iW >= f.maxelm) {
						f.potent[MUT][iW] = f.potent[MUT][iP];
					}
				}
				if (iN > -1) {
					if (iN >= f.maxelm) {
						f.potent[MUT][iN] = f.potent[MUT][iP];
					}
				}
				if (iS > -1) {
					if (iS >= f.maxelm) {
						f.potent[MUT][iS] = f.potent[MUT][iP];
					}
				}
				if (iT > -1) {
					if (iT >= f.maxelm) {
						f.potent[MUT][iT] = f.potent[MUT][iP];
					}
				}
				if (iB > -1) {
					if (iB >= f.maxelm) {
						f.potent[MUT][iB] = f.potent[MUT][iP];
					}
				}

				integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
				iE2 = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE2; iN2 = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE2; iT2 = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE2;
				iW2 = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE2; iS2 = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE2; iB2 = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE2;
				if (iE2 > -1) {
					if (iE2 >= f.maxelm) {
						f.potent[MUT][iE2] = f.potent[MUT][iP];
					}
				}
				if (iW2 > -1) {
					if (iW2 >= f.maxelm) {
						f.potent[MUT][iW2] = f.potent[MUT][iP];
					}
				}
				if (iN2 > -1) {
					if (iN2 >= f.maxelm) {
						f.potent[MUT][iN2] = f.potent[MUT][iP];
					}
				}
				if (iS2 > -1) {
					if (iS2 >= f.maxelm) {
						f.potent[MUT][iS2] = f.potent[MUT][iP];
					}
				}
				if (iT2 > -1) {
					if (iT2 >= f.maxelm) {
						f.potent[MUT][iT2] = f.potent[MUT][iP];
					}
				}
				if (iB2 > -1) {
					if (iB2 >= f.maxelm) {
						f.potent[MUT][iB2] = f.potent[MUT][iP];
					}
				}

				integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
				iE3 = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE3; iN3 = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE3; iT3 = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE3;
				iW3 = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE3; iS3 = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE3; iB3 = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE3;
				if (iE3 > -1) {
					if (iE3 >= f.maxelm) {
						f.potent[MUT][iE3] = f.potent[MUT][iP];
					}
				}
				if (iW3 > -1) {
					if (iW3 >= f.maxelm) {
						f.potent[MUT][iW3] = f.potent[MUT][iP];
					}
				}
				if (iN3 > -1) {
					if (iN3 >= f.maxelm) {
						f.potent[MUT][iN3] = f.potent[MUT][iP];
					}
				}
				if (iS3 > -1) {
					if (iS3 >= f.maxelm) {
						f.potent[MUT][iS3] = f.potent[MUT][iP];
					}
				}
				if (iT3 > -1) {
					if (iT3 >= f.maxelm) {
						f.potent[MUT][iT3] = f.potent[MUT][iP];
					}
				}
				if (iB3 > -1) {
					if (iB3 >= f.maxelm) {
						f.potent[MUT][iB3] = f.potent[MUT][iP];
					}
				}

				integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов
				iE4 = f.neighbors_for_the_internal_node[E_SIDE][iP].iNODE4; iN4 = f.neighbors_for_the_internal_node[N_SIDE][iP].iNODE4; iT4 = f.neighbors_for_the_internal_node[T_SIDE][iP].iNODE4;
				iW4 = f.neighbors_for_the_internal_node[W_SIDE][iP].iNODE4; iS4 = f.neighbors_for_the_internal_node[S_SIDE][iP].iNODE4; iB4 = f.neighbors_for_the_internal_node[B_SIDE][iP].iNODE4;
				if (iE4 > -1) {
					if (iE4 >= f.maxelm) {
						f.potent[MUT][iE4] = f.potent[MUT][iP];
					}
				}
				if (iW4 > -1) {
					if (iW4 >= f.maxelm) {
						f.potent[MUT][iW4] = f.potent[MUT][iP];
					}
				}
				if (iN4 > -1) {
					if (iN4 >= f.maxelm) {
						f.potent[MUT][iN4] = f.potent[MUT][iP];
					}
				}
				if (iS4 > -1) {
					if (iS4 >= f.maxelm) {
						f.potent[MUT][iS4] = f.potent[MUT][iP];
					}
				}
				if (iT4 > -1) {
					if (iT4 >= f.maxelm) {
						f.potent[MUT][iT4] = f.potent[MUT][iP];
					}
				}
				if (iB4 > -1) {
					if (iB4 >= f.maxelm) {
						f.potent[MUT][iB4] = f.potent[MUT][iP];
					}
				}
			}

		}

    // экспорт результата вычисления в программу tecplot360:
	if (0) {
		 exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior,iflow,false,0, b, lb);
	     printf("temperature calculate begin now... OK.\n");
	     //getchar(); // debug avtosave
		 system("pause");
	}

	
	// btempdepend   если требуется совместное решение гидродинамики и теплопроводности.
	// совместное решение требуется в случае если параметры текучей среды зависят от температуры или
	// в случае если мы имеем дело с приближением Буссинеска.
	// btempdepend == false если уравнения гидродинамики и теплопроводности можно решить раздельно:
	// сначала гидродинамику, а потом теплопроводность.
	bool btempdepend=true;
	if (itempersolve==1) {
         btempdepend=true; // нужно решать уравнение теплопроводности.
	}
	else if (itempersolve==0) {
		btempdepend=false; // уравнение теплопроводности дополнительно решать ненужно.
	}
	if (1&&(btempdepend)) {
		// ВНИМАНИЕ !!! Решение уравнения теплопередачи пока выключено и если его необходимо включить то нужно
		// предусмотреть корректную обработку структуры m отвечающую за память в BICGStab_internal3. 
		// Пока всё просто выключено и работать не будет. 12 апреля 2013 года.

		// Обновляем свойства твёрдых материалов и свойства теплоносителя.
		update_temp_properties(t, fglobal, b, lb, matlist);

        // Т.к. подавляющее большинство свойств материалов зависит от температуры,
		// то необходимо решать уравнение теплопроводности на каждой итерации
		// SIMPLE алгоритма.
        //printf("Temp\n");
	    //solve(TEMP, res, f, fglobal, t, rhie_chow, s, w, b, ls, lw, lb, dbeta, flow_interior, true, false,nullptr,0.01,false,dgx,dgy,dgz); 
		doublereal dbeta_temp=1.0; // 1.0 - первый порядок, 1.2 третий порядок
		// Нелинейный солвер теплопроводности:
		if (bprintmessage) {
		    printf("TEMP\n");
	    }
		bool bOkSc = false;
		// Обновление мощности тепловыделения во всех внутренних узлах.
		for (integer i47 = 0; i47 < t.maxelm; i47++) {
			// Скорость в том что значение не вычисляется как раньше а просто хранится.
			integer ib = t.whot_is_block[i47];
			t.Sc[i47] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, t.potent[i47]);

			if (t.Sc[i47] > 1.0e-15) {
				bOkSc = true;
				//printf("myVersion SIMPLE algorithm... Sc TEMP >0\n");
				//getchar();
			}
		}
		if (!bOkSc) {
			// Должно быть по крайней мере две стенки с разной температурой 
			// тогда это предупреждение можно не писать 30.10.2019.
			bool bprint1 = true;
			if (lw > 2) {
				
				doublereal Tamb1 = w[0].Tamb;
					int j1 = -1;
					for (j1 = 0; j1 < lw; j1++) {
						//if (w[j1].ifamily == DIRICHLET_FAMILY) 
						if ((!w[j1].bpressure)&&(!w[j1].bsymmetry))
						{
							Tamb1 = w[j1].Tamb;
							break;
						}
					}
					if (j1 > -1) {
						for (int j2 = j1+1; j2 < lw; j2++) {
							//if (w[j2].ifamily == DIRICHLET_FAMILY) 
							if ((!w[j1].bpressure) && (!w[j1].bsymmetry))
							{
								if (Tamb1 != w[j2].Tamb) {
									// Найдено две стенки с разными температурами на них.
									bprint1 = false;
									break;
								}
							}
						}
					}
			}
			if (bprint1) {
				printf("WARNING !!! ZERO POWER in cfd / temperature calculation.\n");
				printf("function myVersion_SIMPLE_algorithm... Sc TEMP == 0.\n");
			}
			//system("PAUSE");
		}
		
		// Если мы вызываем solve_nonlinear_temp то он будет решать уравнение теплопередачи до упора,
		// в тоже время при моделировании естественной конвекции поле скорости еще не развито и мы получим колосальный перегрев и как 
		// вывод расходимость при решении задачи. В задачах с естественной конвекцией нужно сделать столько же итераций по температуре сколько и итераций 
		// для компонент скорости.
		
		/*		
		solve_nonlinear_temp(f, fglobal, t, rhie_chow,
			                 b, lb, s, ls, w, lw,
							 dbeta_temp, flow_interior,
							 true, toldtimestep, dtimestep, dtimestepold,
							 btimedep, matlist, inumiter,
							 bprintmessage, gtdps, ltdp, 1.0, m,
							 speedoldtimestep, mfoldtimestep);*/
		// последний параметр равный единице означает что мощность подаётся.

		// Освобождаем память.
            //m.bsignalfreeCRSt=true;
			doublereal** rsumanbstuff=nullptr;
			bool  bconvective=true;
			doublereal power_on=1.0;

			// Чтобы что-то показать значение невязки должно вычисляться до решения СЛАУ,
			// т.к. невязка после решения СЛАУ удоалетворяет уравнению тождественно.
			// Нельзя вычислять невязку т.к. матрица ещё несобрана.
			//rfluentrestemp = fluent_residual_for_x(t.slau, t.slau_bon, t.potent, t.maxelm, t.maxbound); // невязка по формуле fluent.
			
#pragma omp parallel for
			for (integer i = 0; i < t.maxelm + t.maxbound; i++) {
				told_temperature_global_for_HOrelax[i] = t.potent[i];
			}

			// Проверка на NAN до решения.
#pragma omp parallel for
			for (integer i = 0; i < t.maxelm + t.maxbound; i++) {
				if (t.potent[i] != t.potent[i]) {
					printf("function: my_version_SIMPLE_Algorithm3D apriory solve TEMP\n");
					printf("t.potent[%lld]=%e\n", i, t.potent[i]);
					system("pause");
				}
			}
			
			// Ещё один вызов решателя не должен повлиять на решение, т.к. оно и так уже получено.
			solve(TEMP, // идентификатор уравнения (это уравнение теплопередачи).
				   res, // невязка
				   f, 
				   fglobal,
				   t,
				   rhie_chow,
				   s, w, b, ls, lw, lb, // объекты (источники, стенки, блоки).
				   dbeta,
				   flow_interior,
				   bconvective,
				   false,
				   toldtimestep, // поле температур с предыдущего временного слоя
				   nullptr,
				   speedoldtimestep, // скорость с предыдущего временного слоя
				   mfoldtimestep, // конвективный поток через грани КО с предыдущего временного слоя
				   dtimestep, // размер шага по времени
				   btimedep, // стационарный или нестационарный солвер
				   dgx, dgy, dgz, // ускорение свободного падения (natural convection id)
				   matlist, // параметры материалов
				   inumiter,bprintmessage,RCh,false,
				   nullptr,rsumanbstuff,false,false, power_on, m,
				   rthdsdt, rfluentrestemp, lu, my_union, color_solid, dist_max_solid); // номер глобальной итерации
			//delete[] color_solid;

			// Проверка на NAN после решения.
#pragma omp parallel for
			for (integer i = 0; i < t.maxelm + t.maxbound; i++) {
				if (t.potent[i] != t.potent[i]) {
					printf("function: my_version_SIMPLE_Algorithm3D post solve TEMP\n");
					printf("t.potent[%lld]=%e\n",i, t.potent[i]);
					system("pause");
				}
			}

			//doublereal tmax = 0.0;
			//for (integer i1 = 0; i1<t.maxelm + t.maxbound; i1++) tmax = fmax(tmax, fabs(t.potent[i1]));
			//printf("apost SIMPLE internal: maximum temperature in default interior is %1.4e\n", tmax);
			//getchar();

		doublereal fHORF = 0.25; // for steady state problem.
		if ((fabs(dgx) > 1.0e-20) || (fabs(dgy) > 1.0e-20) || (fabs(dgz) > 1.0e-20)) {
			// Натуральная конвекция в cfd.
			// Оставим дефолтные 0.25.
			bool b_natural_convection = true;
			for (int k21 = 0; k21 < lw; k21++) {
				if ((w[k21].Vx*w[k21].Vx + w[k21].Vy*w[k21].Vy + w[k21].Vz*w[k21].Vz) > 1.0e-20) {
					b_natural_convection = false; // Вынужденная конвекция.
				}
			}
			if (b_natural_convection) {
				if (inumiter < 7) {
					fHORF = 0.01;
				}
			}
			//fHORF=my_amg_manager.F_to_F_Stress;
		}
		if (btimedep) { // unsteady problems.
			fHORF = 0.75; // ANSYS Fluent Theory Guide.
		}
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm + t.maxbound; i++) {
			t.potent[i] = told_temperature_global_for_HOrelax[i] + 
				fHORF * (t.potent[i] - told_temperature_global_for_HOrelax[i]);
		}
#pragma omp parallel for
		for (integer i = 0; i < t.maxelm + t.maxbound; i++) {
			told_temperature_global_for_HOrelax[i] = t.potent[i];
		}

		// Обновляем свойства твёрдых материалов и свойства теплоносителя.
		update_temp_properties(t, fglobal, b, lb, matlist);

	}

	//if (inumiter >= 96) getchar();

	// 28.07.2016
	//exporttecplotxy360T_3D_part2(t.maxelm, t.ncell, fglobal, t, flow_interior, inumiter, bextendedprint);
	//getchar(); // debug

	//printf("TEMP end\n");
	//getchar();

	//exporttecplotxy360_3D( f.maxelm, f.ncell, f.nvtx, f.nvtxcell, f.pa, f.potent, rhie_chow);
	//getchar();

	// экспорт результата вычисления в программу tecplot360:
	if (0) {
		 exporttecplotxy360T_3D_part2(t.maxelm,t.ncell, fglobal, t, flow_interior,iflow,false,0, b, lb);
	     printf("one iteration finish... OK.\n");
	     //getchar(); // debug avtosave
		 system("pause");
	}

	/*if (inumiter==30) {
		getchar();
	}*/

	// Освобождение занимаемой оперативной памяти 
	for (integer i=0; i<3; i++) delete[] tau[i];
	delete[] tau; // Освобождаем память из под сглаженного шага по псевдовремени.
	tau = nullptr;

	// Освобождение памяти из под недособранных диагональных коэффициентов СЛАУ.
	for (integer i=0; i<3; i++) delete[] sumanb[i];
	delete[] sumanb;
	sumanb = nullptr;

	

} // my_version_SIMPLE_Algorithm3D



#endif
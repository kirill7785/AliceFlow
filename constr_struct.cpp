// Файл constr_struct.cpp собирает структуры
// для решения задачи сопряжённого теплообмена.
// Данная версия работает на декартовой структурированной сетке с семиточечным шаблоном.
// Функционал данного модуля универсальный, он также предназначен также для вызова
// соответствующих функций из модуля constr_struct_alice.cpp на адаптивной локально измельченной сетке АЛИС.
// 26.09.2016 Написание данного модуля в основном завершено.

#pragma once
#ifndef _CONSTR_STRUCT_CPP_
#define _CONSTR_STRUCT_CPP_ 1

// тот самый стек.
typedef struct TSTEK {
	int i;
	int j;
	int k;

	TSTEK() {
		i = -1;
		j = -1;
		k = -1;
	}
} STEK;

// Число вершин конечного элемента.
// 23.09.2020
int NUMBER_OF_VERTEX_FINITE_ELEMENT() {
	// Конечный элемент в форме гексаэдра.
	return 8;
}

integer get_max_array_elm_integer1(integer*& a, const integer n) {

#ifdef _OPENMP
	// многопоточная версия
	integer dmax = INT_MIN;

#pragma omp parallel
	{
		integer dmax_loc = INT_MIN;

#pragma omp for 
		for (integer i = 1; i <= n; ++i) {
			integer z = a[i];
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
	integer dmax = INT_MIN;
	for (integer i = 1; i <= n; ++i) {
		if (a[i] > dmax) dmax = a[i];
	}
	return dmax;
#endif

}


integer get_max_array_elm_integer(integer *& a, const integer n) {

#ifdef _OPENMP
	// многопоточная версия
	integer dmax = INT_MIN;

#pragma omp parallel
	{
		integer dmax_loc = INT_MIN;

#pragma omp for 
		for (integer i = 0; i < n; ++i) {
			integer z = a[i];
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
	integer dmax = INT_MIN;
	for (integer i = 0; i < n; ++i) {
		if (a[i] > dmax) dmax = a[i];
	}
	return dmax;
#endif

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
		for (integer i = 0; i < n; ++i) {
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
	for (integer i = 0; i < n; ++i) {
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
		for (integer i = 0; i < n; ++i) {
			doublerealT z = fabs(a1[i] - a2[i]);
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
	for (integer i = 0; i < n; ++i) {
		doublerealT z = fabs(a1[i] - a2[i]);
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
		for (integer i = 0; i < n; ++i) {
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
	for (integer i = 0; i < n; ++i) {
		if (a[i] < dmin) dmin = a[i];
	}
	return dmin;
#endif

}

// Обёртка для сортировки (выбор метода сортировки).
// не использует дополнительной оперативной памяти.
// in - предполагается достаточно малым меньше 500,
// array[0...in]
template <typename myARRT>
void Sort_method(myARRT* &rb, integer in);

// Процедура использующаяся в функции allocation_memory_flow.
// Реализует два уровня octree дерева и используется при поиске точки на сетке.
// Пример вызова.
// CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(0, xc47,yc47,zc47,i_47,oct_load1,i_22,avgx,avgy,avgz);
void CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(integer identifikator,
	doublereal xc47, doublereal yc47, doublereal zc47, integer i_47,
	integer*** &oct_load1, integer** &i_22, 
	doublereal avgx[], doublereal avgy[], doublereal avgz[]) {

if (xc47 < avgx[identifikator]) {
	if (yc47 < avgy[identifikator]) {
		if (zc47 < avgz[identifikator]) {
			oct_load1[identifikator][0][i_22[identifikator][0]] = i_47;
			i_22[identifikator][0]++;
		}
		else {
			oct_load1[identifikator][1][i_22[identifikator][1]] = i_47;
			i_22[identifikator][1]++;
		}
	}
	else {
		if (zc47 < avgz[identifikator]) {
			oct_load1[identifikator][2][i_22[identifikator][2]] = i_47;
			i_22[identifikator][2]++;
		}
		else {
			oct_load1[identifikator][3][i_22[identifikator][3]] = i_47;
			i_22[identifikator][3]++;
		}
	}
}
else {
	if (yc47 < avgy[identifikator]) {
		if (zc47 < avgz[identifikator]) {
			oct_load1[identifikator][4][i_22[identifikator][4]] = i_47;
			i_22[identifikator][4]++;
		}
		else {
			oct_load1[identifikator][5][i_22[identifikator][5]] = i_47;
			i_22[identifikator][5]++;
		}
	}
	else {
		if (zc47 < avgz[identifikator]) {
			oct_load1[identifikator][6][i_22[identifikator][6]] = i_47;
			i_22[identifikator][6]++;
		}
		else {
			oct_load1[identifikator][7][i_22[identifikator][7]] = i_47;
			i_22[identifikator][7]++;
		}
	}
}
}// CONSTRUCT_SECONDARY_LEVEL_OCTREE


// Макрос №2 использующийся в функции allocation_memory_flow.
// Реализует два уровня octree дерева и используется при поиске точки на сетке.
// CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(0, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);
//
// Очень медленный вариант: перемудрил с параллельностью. Рекомендуется
// использовать существенно более быстрый вариант от 19.06.2021
void CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1_old(integer identifikator,
	doublereal xc47, doublereal yc47, doublereal zc47, 
	integer*** &oct_load1, integer** &i_22,
	doublereal avgx[], doublereal avgy[], doublereal avgz[],
	doublereal* &x47, doublereal* &y47, doublereal* &z47,
	int** &nvtx47, bool &bfound_gl, integer &ifound_gl) {

if (xc47 < avgx[identifikator]) {
	if (yc47 < avgy[identifikator]) {
		if (zc47 < avgz[identifikator]) {


#pragma omp parallel
			{

				bool bfound = false;
				integer ifound = -1;

#pragma omp for 
				for (integer i_471 = 0; i_471 < i_22[identifikator][0]; i_471++)
				{
					if (bfound) continue;

					integer i_47 = oct_load1[identifikator][0][i_471];
					if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
						(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
						(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
					{
						ifound = i_47;
						bfound = true;
						//break;
					}


#pragma omp critical 
					{
						if (bfound) {
							bfound_gl = true;
							ifound_gl = ifound;
						}
					}
				}

			}
		}
		else {



#pragma omp parallel
			{

				bool bfound = false;
				integer ifound = -1;

#pragma omp for 
				for (integer i_471 = 0; i_471 < i_22[identifikator][1]; i_471++)
				{
					if (bfound) continue;

					integer i_47 = oct_load1[identifikator][1][i_471];
					if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
						(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
						(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
					{
						ifound = i_47;
						bfound = true;
						//break;
					}


#pragma omp critical 
					{
						if (bfound) {
							bfound_gl = true;
							ifound_gl = ifound;
						}
					}
				}

			}
		}
	}
	else {


		if (zc47 < avgz[identifikator]) {


#pragma omp parallel
			{

				bool bfound = false;
				integer ifound = -1;

#pragma omp for 
				for (integer i_471 = 0; i_471 < i_22[identifikator][2]; i_471++)
				{
					if (bfound) continue;

					integer i_47 = oct_load1[identifikator][2][i_471];
					if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
						(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
						(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
					{
						ifound = i_47;
						bfound = true;
						//break;
					}



#pragma omp critical 
					{
						if (bfound) {
							bfound_gl = true;
							ifound_gl = ifound;
						}
					}
				}

			}
		}
		else {

#pragma omp parallel
			{

				bool bfound = false;
				integer ifound = -1;

#pragma omp for 
				for (integer i_471 = 0; i_471 < i_22[identifikator][3]; i_471++)
				{
					if (bfound) continue;

					integer i_47 = oct_load1[identifikator][3][i_471];
					if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
						(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
						(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
					{
						ifound = i_47;
						bfound = true;
						//break;
					}


#pragma omp critical 
					{
						if (bfound) {
							bfound_gl = true;
							ifound_gl = ifound;
						}
					}
				}

			}
		}
	}
}
else {
	if (yc47 < avgy[identifikator]) {
		if (zc47 < avgz[identifikator]) {


#pragma omp parallel
			{

				bool bfound = false;
				integer ifound = -1;

#pragma omp for 
				for (integer i_471 = 0; i_471 < i_22[identifikator][4]; i_471++)
				{
					if (bfound) continue;

					integer i_47 = oct_load1[identifikator][4][i_471];
					if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
						(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
						(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
					{
						ifound = i_47;
						bfound = true;
						//break;
					}


#pragma omp critical 
					{
						if (bfound) {
							bfound_gl = true;
							ifound_gl = ifound;
						}
					}
				}

			}

		}
		else {


#pragma omp parallel
			{

				bool bfound = false;
				integer ifound = -1;

#pragma omp for 
				for (integer i_471 = 0; i_471 < i_22[identifikator][5]; i_471++)
				{
					if (bfound) continue;

					integer i_47 = oct_load1[identifikator][5][i_471];
					if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
						(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
						(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
					{
						ifound = i_47;
						bfound = true;
						//break;
					}


#pragma omp critical 
					{
						if (bfound) {
							bfound_gl = true;
							ifound_gl = ifound;
						}
					}
				}
			}
		}
	}
	else {
		if (zc47 < avgz[identifikator]) {


#pragma omp parallel
			{

				bool bfound = false;
				integer ifound = -1;

#pragma omp for 
				for (integer i_471 = 0; i_471 < i_22[identifikator][6]; i_471++)
				{
					if (bfound) continue;

					integer i_47 = oct_load1[identifikator][6][i_471];
					if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
						(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
						(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
					{
						ifound = i_47;
						bfound = true;
						//break;
					}


#pragma omp critical 
					{
						if (bfound) {
							bfound_gl = true;
							ifound_gl = ifound;
						}
					}
				}
			}
		}
		else {

#pragma omp parallel
			{

				bool bfound = false;
				integer ifound = -1;

#pragma omp for 
				for (integer i_471 = 0; i_471 < i_22[identifikator][7]; i_471++)
				{
					if (bfound) continue;

					integer i_47 = oct_load1[identifikator][7][i_471];
					if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
						(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
						(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
					{
						ifound = i_47;
						bfound = true;
						//break;
					}


#pragma omp critical
					{
						if (bfound) {
							bfound_gl = true;
							ifound_gl = ifound;
						}
					}
				}
			}
		}
	}
}
}



// Макрос №2 использующийся в функции allocation_memory_flow.
// Реализует два уровня octree дерева и используется при поиске точки на сетке.
// CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(0, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);
//
// Более простой вариант распараллеливания. Работает более быстро в многопоточном режиме.
// Проверено, работает коректно. 19.06.2021.
void CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(integer identifikator,
	doublereal xc47, doublereal yc47, doublereal zc47,
	integer***& oct_load1, integer**& i_22,
	doublereal avgx[], doublereal avgy[], doublereal avgz[],
	doublereal*& x47, doublereal*& y47, doublereal*& z47,
	int**& nvtx47, bool& bfound_gl, integer& ifound_gl) {

	if (xc47 < avgx[identifikator]) {
		if (yc47 < avgy[identifikator]) {
			if (zc47 < avgz[identifikator]) {


//#pragma omp parallel
				{

					bool bfound = false;
					integer ifound = -1;

					const integer istop = i_22[identifikator][0];

#pragma omp parallel for 
					for (integer i_471 = 0; i_471 < istop; ++i_471)
					{
						//if (bfound) continue;

						integer i_47 = oct_load1[identifikator][0][i_471];
						if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
							(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
							(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
						{
							ifound = i_47;
							bfound = true;
							//break;
						}
					}

					if (bfound) {
						bfound_gl = true;
						ifound_gl = ifound;
					}


				}
			}
			else {




				{

					bool bfound = false;
					integer ifound = -1;

					const integer istop = i_22[identifikator][1];

#pragma omp parallel for 
					for (integer i_471 = 0; i_471 < istop; ++i_471)
					{
						//if (bfound) continue;

						integer i_47 = oct_load1[identifikator][1][i_471];
						if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
							(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
							(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
						{
							ifound = i_47;
							bfound = true;
							//break;
						}

					}

					if (bfound) {
						bfound_gl = true;
						ifound_gl = ifound;
					}

				}
			}
		}
		else {


			if (zc47 < avgz[identifikator]) {



				{

					bool bfound = false;
					integer ifound = -1;

					const integer istop = i_22[identifikator][2];

#pragma omp parallel for 
					for (integer i_471 = 0; i_471 < istop; ++i_471)
					{
						//if (bfound) continue;

						integer i_47 = oct_load1[identifikator][2][i_471];
						if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
							(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
							(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
						{
							ifound = i_47;
							bfound = true;
							//break;
						}
					}

					if (bfound) {
						bfound_gl = true;
						ifound_gl = ifound;
					}

				}
			}
			else {


				{

					bool bfound = false;
					integer ifound = -1;

					const integer istop = i_22[identifikator][3];

#pragma omp parallel for 
					for (integer i_471 = 0; i_471 < istop; ++i_471)
					{
						//if (bfound) continue;

						integer i_47 = oct_load1[identifikator][3][i_471];
						if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
							(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
							(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
						{
							ifound = i_47;
							bfound = true;
							//break;
						}
					}

					if (bfound) {
						bfound_gl = true;
						ifound_gl = ifound;
					}

				}
			}
		}
	}
	else {
		if (yc47 < avgy[identifikator]) {
			if (zc47 < avgz[identifikator]) {



				{

					bool bfound = false;
					integer ifound = -1;
					const integer istop = i_22[identifikator][4];


#pragma omp parallel for 
					for (integer i_471 = 0; i_471 < istop; ++i_471)
					{
						//if (bfound) continue;

						integer i_47 = oct_load1[identifikator][4][i_471];
						if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
							(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
							(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
						{
							ifound = i_47;
							bfound = true;
							//break;
						}
					}

					if (bfound) {
						bfound_gl = true;
						ifound_gl = ifound;
					}

				}

			}
			else {



				{

					bool bfound = false;
					integer ifound = -1;

					const integer istop = i_22[identifikator][5];

#pragma omp parallel for 
					for (integer i_471 = 0; i_471 < istop; ++i_471)
					{
						//if (bfound) continue;

						integer i_47 = oct_load1[identifikator][5][i_471];
						if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
							(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
							(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
						{
							ifound = i_47;
							bfound = true;
							//break;
						}
					}

					if (bfound) {
						bfound_gl = true;
						ifound_gl = ifound;
					}
				}
			}
		}
		else {
			if (zc47 < avgz[identifikator]) {



				{

					bool bfound = false;
					integer ifound = -1;

					const integer istop = i_22[identifikator][6];

#pragma omp parallel for 
					for (integer i_471 = 0; i_471 < istop; ++i_471)
					{
						//if (bfound) continue;

						integer i_47 = oct_load1[identifikator][6][i_471];
						if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
							(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
							(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
						{
							ifound = i_47;
							bfound = true;
							//break;
						}
					}

					if (bfound) {
						bfound_gl = true;
						ifound_gl = ifound;
					}
				}
			}
			else {


				{

					bool bfound = false;
					integer ifound = -1;

					const integer istop = i_22[identifikator][7];

#pragma omp parallel for 
					for (integer i_471 = 0; i_471 < istop; ++i_471)
					{
						//if (bfound) continue;

						integer i_47 = oct_load1[identifikator][7][i_471];
						if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
							(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
							(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
						{
							ifound = i_47;
							bfound = true;
							//break;
						}

					}

					if (bfound) {
						bfound_gl = true;
						ifound_gl = ifound;
					}
				}
			}
		}
	}
}

const unsigned char TEMPERATURE = 0;
const unsigned char HYDRODINAMIC = 1;


const unsigned char RHO = 0; // плотность
const unsigned char HEAT_CAPACITY = 1; // теплоёмкость при постоянном давлении
const unsigned char LAM = 2; // теплопроводность
const unsigned char MULT_LAM_X = 3;
const unsigned char MULT_LAM_Y = 4;
const unsigned char MULT_LAM_Z = 5;
//const unsigned char MU_LAME = 6;  // Коэффициент Ламе
//const unsigned char LAMBDA_LAME = 7; // Коэффициент Ламе.
const unsigned char BETA_T_MECHANICAL = 8; // коэффициент линейного теплового расширения для Механики.
const unsigned char YOUNG_MODULE = 9; // Молдуль Юнга.
const unsigned char POISSON_RATIO = 10; // Poisson ratio
// множитель к коэффициенту линейного теплового расширения.
const unsigned char MULT_BETA_T_MECHANICAL_X = 11;
const unsigned char MULT_BETA_T_MECHANICAL_Y = 12;
const unsigned char MULT_BETA_T_MECHANICAL_Z = 13;
// множитель к модулю Юнга.
const unsigned char MULT_YOUNG_MODULE_X = 14;
const unsigned char MULT_YOUNG_MODULE_Y = 15;
const unsigned char MULT_YOUNG_MODULE_Z = 16;
// Множитель к ортотропному коэффициенту Пуассона.
const unsigned char MULT_POISSON_RATIO_YZ = 17;
const unsigned char MULT_POISSON_RATIO_XZ = 18;
const unsigned char MULT_POISSON_RATIO_XY = 19;
const unsigned char MULT_POISSON_RATIO_ZY = 20;
const unsigned char MULT_POISSON_RATIO_ZX = 21;
const unsigned char MULT_POISSON_RATIO_YX = 22;
// Модуль сдвига
const unsigned char SHEAR_MODULE_YZ = 23;
const unsigned char SHEAR_MODULE_XZ = 24;
const unsigned char SHEAR_MODULE_XY = 25;

// Число параметров, задающих свойства материала при расчёте.
const unsigned char SIZE_PROPERTIES_ARRAY = 26;

const unsigned char MU_DYNAMIC_VISCOSITY = 1; // динамическая вязкость
const unsigned char BETA_T = 2; // коэффициент линейного температурного расширения

const unsigned char VELOCITY_X_COMPONENT = 0; // горизонтальная скорость (обязательно должен начинаться с нуля)
const unsigned char VELOCITY_Y_COMPONENT = 1; // вертикальная скорость
const unsigned char VELOCITY_Z_COMPONENT = 2;
const unsigned char PRESS = 3; // давление
const unsigned char PAM = 4; // поправка давления
const unsigned char NUSHA_SL = 5; // модифицированная кинематическая турбулентная вязкость. Спаларт-Аллмарес.
const unsigned char TURBULENT_KINETIK_ENERGY_SL = 6; // Кинетическая энергия турбулентных пульсаций.
const unsigned char TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL = 7; // Удельная скорость диссипации (omega).
const unsigned char TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL = 8; // Кинетическая энергия турбулентных пульсаций двухслойная модель на основе стандартной k-epsilon модели.
const unsigned char TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL = 9; // Скорость диссипации кинетической энергии турбулентности.
// СЛАУ для модели ламинарно турбулентного перехода Лангтрии - Ментера.
const unsigned char GAMMA_LANGTRY_MENTER_SL = 10;
const unsigned char RE_THETA_LANGTRY_MENTER_SL = 11;
const unsigned char VXCOR = 5; // скорости с предыдущей  
const unsigned char VYCOR = 6; // итерации удовлетворяющие
const unsigned char VZCOR = 7; // уравнению неразрывности
// Градиенты скоростей в центрах КО:
const unsigned char GRADXVX = 8;
const unsigned char GRADYVX = 9;
const unsigned char GRADZVX = 10;
const unsigned char GRADXVY = 11;
const unsigned char GRADYVY = 12;
const unsigned char GRADZVY = 13;
const unsigned char GRADXVZ = 14;
const unsigned char GRADYVZ = 15;
const unsigned char GRADZVZ = 16;
// Коэффициент динамической турбулентной вязкости
const unsigned char MUT = 17;
const unsigned char CURL = 18; // модуль ротора скорости (завихрённость).
const unsigned char FBUF = 19;// буфер в котором могут находится разные используемые внутри программы вектора: распределение константы Смагоринского в пространстве и т.п.
// Градиенты поправки давления.
const unsigned char GRADXPAM = 20;
const unsigned char GRADYPAM = 21;
const unsigned char GRADZPAM = 22;
// Градиенты давления используются для вычисления поправки Рхи-Чоу.
const unsigned char GRADXPRESS = 23;
const unsigned char GRADYPRESS = 24;
const unsigned char GRADZPRESS = 25;
const unsigned char PAMOLDITER = 26;
const unsigned char TOTALDEFORMATIONVAR = 27;// используется в сглаживателе Чебышева.
const unsigned char SPEED = 28; // Модуль скорости теплоносителя.
// Spalart-Allmares (RANS). 19.04.2019
const unsigned char NUSHA = 29; // Модифицированная кинематическая турбулентная вязкость.
// Градиенты модифицированной кинематической турбулентной вязкости.
const unsigned char GRADXNUSHA = 30;
const unsigned char GRADYNUSHA = 31;
const unsigned char GRADZNUSHA = 32;
// SST модель турбулентности Ментера (RANS) 03.10.2019.
const unsigned char TURBULENT_KINETIK_ENERGY = 33; // Кинетическая энергия турбулентных пульсаций.
const unsigned char TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA = 34; // Удельная скорость диссипации (omega).
// Градиенты  кинетической энергии турбулентных пульсаций.
const unsigned char GRADXTURBULENT_KINETIK_ENERGY = 35;
const unsigned char GRADYTURBULENT_KINETIK_ENERGY = 36;
const unsigned char GRADZTURBULENT_KINETIK_ENERGY = 37;
// Градиенты удельной скорости её диссипации
const unsigned char GRADXTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA = 38;
const unsigned char GRADYTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA = 39;
const unsigned char GRADZTURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA = 40;
// Стандартная k-epsilon модель на основе двухслойной модели.
const unsigned char TURBULENT_KINETIK_ENERGY_STD_K_EPS = 41; // Кинетическая энергия турбулентных пульсаций.
const unsigned char TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS = 42; // Скорость её диссипации (epsilon).
// Градиенты  кинетической энергии турбулентных пульсаций.
const unsigned char GRADXTURBULENT_KINETIK_ENERGY_STD_K_EPS = 43;
const unsigned char GRADYTURBULENT_KINETIK_ENERGY_STD_K_EPS = 44;
const unsigned char GRADZTURBULENT_KINETIK_ENERGY_STD_K_EPS = 45;
// Градиенты скорости её диссипации (градиенты epsilon)
const unsigned char GRADXTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS = 46;
const unsigned char GRADYTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS = 47;
const unsigned char GRADZTURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS = 48;
// Модель Лангри Ментера Ламнарно Турбулентного перехода.
const unsigned char GAMMA_LANGTRY_MENTER = 49;
const unsigned char RE_THETA_LANGTRY_MENTER = 50;
// Градиенты ненуждны gamma и Re_Thetta_T ненужны. !!!
// градиенты перемежаемости Gamma
//const unsigned char GRADX_GAMMA_LANGTRY_MENTER = 51;
//const unsigned char GRADY_GAMMA_LANGTRY_MENTER = 52;
//const unsigned char GRADZ_GAMMA_LANGTRY_MENTER = 53;
// градиенты ReTHETA
//const unsigned char GRADX_RE_THETA_LANGTRY_MENTER = 54;
//const unsigned char GRADY_RE_THETA_LANGTRY_MENTER = 55;
//const unsigned char GRADZ_RE_THETA_LANGTRY_MENTER = 56;

// Параметры турбулентности с предыдущего временного шага.
const unsigned char TURBULENT_KINETIK_ENERGY_MENTER_SST_OLD_TIME_STEP = 0;
const unsigned char TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_OLD_TIME_STEP = 1;
const unsigned char TURBULENT_NUSHA_OLD_TIME_STEP = 2;
const unsigned char GAMMA_LANGTRY_MENTER_OLD_TIME_STEP = 3;
const unsigned char RE_THETA_LANGTRY_MENTER_OLD_TIME_STEP = 4;
const unsigned char iNUMBER_FUNCTION_TURBULENT_OLD_TIME_STEP = 5; // Число функций для расчёта турбулентных течений на предыдущем временном шаге.

// Число вектор функции хранимых при расчёте cfd задачи.
const unsigned int SIZE_FLOW_POTENT_ARRAY = 51;

// Static Structural
// Деформации:
const unsigned char TOTALDEFORMATION = 0; // полная
const unsigned char XDEFORMATION = 1; // х компонента
const unsigned char YDEFORMATION = 2; // y компонента
const unsigned char ZDEFORMATION = 3; // z компонента
const unsigned char STRAIN_X = 4; // Деформация по оси ox
const unsigned char STRAIN_Y = 5; // Деформация по оси oy
const unsigned char STRAIN_Z = 6; // Деформация по оси oz
const unsigned char STRAIN_XY = 7;
const unsigned char STRAIN_YZ = 8;
const unsigned char STRAIN_ZX = 9;
// STRAIN von-Mises
const unsigned char STRAIN_VON_MIZES = 10;
const unsigned char LOG10_STRAIN_VON_MIZES = 11;
// 22.08.2020
const unsigned char STRESS_X = 12; // Напряжение по оси ox
const unsigned char STRESS_Y = 13; // Напряжение по оси oy
const unsigned char STRESS_Z = 14; // Напряжение по оси oz
const unsigned char STRESS_XY = 15;
const unsigned char STRESS_YZ = 16;
const unsigned char STRESS_ZX = 17;
// STRESS von-Mises
const unsigned char STRESS_VON_MIZES = 18;
const unsigned char LOG10_STRESS_VON_MIZES = 19;

// Число вектор функции хранимых при расчёте деформации.
const unsigned int SIZE_DEFORMATION_ARRAY = 20;


//const unsigned char SPEED = 99; // velocity magnitude (модуль скорости) константа используется при создании анимации.
const unsigned char TEMP = 100; // уравнение теплопроводности


const unsigned char ENUMERATECONTVOL = 0; // слой нумерации КО
const unsigned char MASKDOMAINFLUID = 1; // слой маска определяет зону FLUID




typedef struct TTOCKA_INT {
	integer i, j, k;
	TTOCKA_INT() : i(0), j(0), k(0) {
		// Конструктор.
		// i = 0, j = 0, k = 0;
	}
} TOCKA_INT;

typedef struct TTOCKA_SHORT_INT {
	int i, j, k;
	TTOCKA_SHORT_INT() : i(0), j(0), k(0) {
		// Конструктор.
		//i = 0, j = 0, k = 0;
	}
} TOCKA_SHORT_INT;

// Объявление функции код которой будет реализован ниже.
// проверка построенной сетки
// экспорт результата расчёта в программу tecplot360
// части 1 и 3.
void exporttecplotxy360T_3D_part1and3(TEMPER &t, integer maxelm, integer maxbound, bool bextendedprint, integer ncell, int** nvtx, int** nvtxcell,
	TOCHKA* pa, BOUND* border_neighbor, integer ivarexport, int** ptr_out);
// Объявление функции код которой будет реализован позже
void xyplot( FLOW* &fglobal, integer flow_interior, TEMPER &t);

// Объявление функции код которой будет реализован намного позже
// по тексту программы в модуле shortest_distance.cpp.
// Вычисление расстояния до ближайшей стенки путём решения СЛАУ
// с симметричной положительно определённой матрицей.
// Этот алгоритм реализован в ANSYS CFX.
void calcdistwallCFX(FLOW &f, integer ls, integer lw, WALL* w);



// освобожение памяти уровня 1
// для уравнения теплопроводности.
void free_level1_temp(TEMPER &t) {

	// По видимому операцию удаления лучше не распараллеливать т.к. 
	// каждое освобождение памяти это запрос к ядру ОС. 
	// По крайней мере приложение упало когда в этой функции стояло распараллеливание. 06.10.2020.

	if (print_log_message) {

		printf("delete temperature internally nodal neighbors\n");
	}
	// neighbors_for_the_internal_node

	delete[] t.bActiveShearModule;
	t.bActiveShearModule = nullptr;

	if (t.neighbors_for_the_internal_node != nullptr) {

		for (int i = 0; i < 12; ++i) {
			if (t.neighbors_for_the_internal_node[i] != nullptr) {
				for (int j = 0; j < 4; ++j) {
					if (t.neighbors_for_the_internal_node[i][j] != nullptr) {
						delete[] t.neighbors_for_the_internal_node[i][j]; 
						t.neighbors_for_the_internal_node[i][j] = nullptr;
					}
				}
			}
		}

		for (int i = 0; i<12; ++i) {
			if (t.neighbors_for_the_internal_node[i] != nullptr) {
				delete[] t.neighbors_for_the_internal_node[i]; // -12N
				t.neighbors_for_the_internal_node[i] = nullptr;
			}
		}
	}
	if (t.neighbors_for_the_internal_node != nullptr) {
		delete[] t.neighbors_for_the_internal_node;
		t.neighbors_for_the_internal_node = nullptr;
	}

//#pragma omp parallel sections
	{

		{
			if (t.whot_is_block != nullptr) {
				delete[] t.whot_is_block;
				t.whot_is_block = nullptr;
			}
		}

//#pragma omp section 
		{
			if (print_log_message) {
				printf("delete temperature Sc\n");
			}
			if (t.Sc != nullptr) {
				delete[] t.Sc; // -N
				t.Sc = nullptr;
			}
		}

//#pragma omp section 
		{
			if (print_log_message) {
				printf("delete temperature ipower_time_depend\n");
			}
			if (t.ipower_time_depend != nullptr) {
				delete[] t.ipower_time_depend; // -N
				t.ipower_time_depend = nullptr;
			}
		}


//#pragma omp section 
		{
			if (print_log_message) {
				printf("delete temperature pa\n");
			}
			if (t.pa != nullptr) {
				delete[] t.pa; // -3N
				t.pa = nullptr;
			}
		}

//#pragma omp section 
		{
			if (print_log_message) {
				printf("delete temperature neighbors in boundary nodes\n");
			}
			// border_neighbor
			if (t.border_neighbor != nullptr) {
				delete[] t.border_neighbor;
				t.border_neighbor = nullptr;
			}
		}

//#pragma omp section 
		{
			if (print_log_message) {
				printf("delete temperature binternal source\n");
			}
			if (t.binternalsource != nullptr) {
				delete[] t.binternalsource;
				t.binternalsource = nullptr;
			}
		}

//#pragma omp section 
		{
			if (print_log_message) {
				printf("delete temperature internal node properties\n");
			}
			// prop
			if (t.prop != nullptr) {
				for (int i = 0; i < SIZE_PROPERTIES_ARRAY; ++i) {
					if (t.prop[i] != nullptr) {
						delete[] t.prop[i]; // -3N
					}
				}
			}
			if (t.prop != nullptr) {
				delete[] t.prop;
				t.prop = nullptr;
			}
			
		}


//#pragma omp section 
		{
			if (print_log_message) {
				printf("delete temperature ptr\n");
			}
			if (t.ptr != nullptr) {
				for (int i = 0; i < 2; ++i) {
					if (t.ptr[i] != nullptr) {
						delete[] t.ptr[i];
						t.ptr[i] = nullptr;
					}
				}
				delete[] t.ptr;
				t.ptr = nullptr;
			}
		}

//#pragma omp section 
		{
			if (print_log_message) {
				printf("delete temperature boundary node properties\n");
			}
			// prop_b
			if (t.prop_b != nullptr) {
				for (int i = 0; i < 6; ++i) {
					if (t.prop_b[i] != nullptr) {
						delete[] t.prop_b[i]; // -3N
						t.prop_b[i] = nullptr;
					}
				}
			}
			if (t.prop_b != nullptr) {
				delete[] t.prop_b;
				t.prop_b = nullptr;
			}
		}

//#pragma omp section 
		{
			if (print_log_message) {
				printf("delete temperature nvtx and nvtxcell\n");
			}
			for (int i = 0; i < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) { // -NUMBER_OF_VERTEX_FINITE_ELEMENT()*N
				if (t.nvtx != nullptr) {
					if (t.nvtx[i] != nullptr) {
						delete[] t.nvtx[i];
						t.nvtx[i] = nullptr;
					}
				}

			}
			if (t.nvtx != nullptr) {
				delete[] t.nvtx;
				t.nvtx = nullptr;
			}
		}

//#pragma omp section 
		{
			for (int i = 0; i < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) { // -NUMBER_OF_VERTEX_FINITE_ELEMENT()*N

				if (t.nvtxcell != nullptr) {
					if (t.nvtxcell[i] != nullptr) {
						delete[] t.nvtxcell[i]; // может быть уже удалено
						t.nvtxcell[i] = nullptr;
					}
				}
			}
			if (t.nvtxcell != nullptr) {
				delete[] t.nvtxcell;
				t.nvtxcell = nullptr;
			}
		}

//#pragma omp section 
		{
			// 26_09_2016.
			if (print_log_message) {
				printf("delete ilevel_alice\n");
			}
			if (t.ilevel_alice != nullptr) {
				delete[] t.ilevel_alice;
				t.ilevel_alice = nullptr;
			}
		}
	}
	// итого -31N

} // free_level1_temp



// освобожение памяти уровня 2
// для уравнения теплопроводности.
void free_level2_temp(TEMPER &t) {

	// Освобождение оперативной памяти выделенной под матрицу СЛАУ.
	// -15N
	if (print_log_message) {
		printf("delete temperature slau\n");
	}
	if (t.slau != nullptr) {
		delete[] t.slau;
		t.slau = nullptr;
	}
	if (print_log_message) {
		printf("delete temperature slau_bon\n");
	}
	if (t.slau_bon != nullptr) {
		delete[] t.slau_bon;
		t.slau_bon = nullptr;
	}

	if (print_log_message) {
		printf("delete temperature potent\n");
	}
	if (t.potent != nullptr) {
		delete[] t.potent;
		t.potent = nullptr;
	}

	if (print_log_message) {
		printf("delete total deformation \n");
	}
	if (t.total_deformation != nullptr) {
		for (integer i_1 = 0; i_1 < SIZE_DEFORMATION_ARRAY; ++i_1) {
			if (t.total_deformation[i_1] != nullptr) {
				delete[] t.total_deformation[i_1];
				t.total_deformation[i_1] = nullptr;
			}
		}

		delete[] t.total_deformation;
		t.total_deformation = nullptr;
	}
} // free_level2_temp

// глобальная память для amg1r5
typedef struct TamgGlobalMemory {
	doublereal* a;
	integer* ia;
	integer* ja;
	doublereal* u;
	doublereal* f;
	integer* ig;
	integer nda;
	integer ndia;
	integer ndja;
	integer ndu;
	integer ndf;
	integer ndig;
} amgGlobalMemory;

// 17.04.2021
// Освобождение памяти в функции main_body().
void free_global_1(TPROP* &matlist, integer lmatmax, doublereal*& x_jacoby_buffer, bool bvery_big_memory, TEMPER& t,
	amgGlobalMemory &amgGM, doublereal*&xposadd, doublereal*& yposadd, doublereal*& zposadd) {

	// Освобождение оперативной памяти.
	if (t.xpos_copy != nullptr) {
		delete[] t.xpos_copy;
		t.xpos_copy = nullptr;
	}
	if (t.ypos_copy != nullptr) {
		delete[] t.ypos_copy;
		t.ypos_copy = nullptr;
	}
	if (t.zpos_copy != nullptr) {
		delete[] t.zpos_copy;
		t.zpos_copy = nullptr;
	}

	// Освобождение оперативной памяти.
	if (xposadd != nullptr) {
		delete[] xposadd;
		xposadd = nullptr;
	}
	if (yposadd != nullptr) {
		delete[] yposadd;
		yposadd = nullptr;
	}
	if (zposadd != nullptr) {
		delete[] zposadd;
		zposadd = nullptr;
	}


	// освобождение памяти из под amg1r5.
	if (amgGM.a != nullptr) {
		delete[] amgGM.a;
		amgGM.a = nullptr;
	}
	if (amgGM.ia != nullptr) {
		delete[] amgGM.ia;
		amgGM.ia = nullptr;
	}
	if (amgGM.ja != nullptr) {
		delete[] amgGM.ja;
		amgGM.ja = nullptr;
	}

	if (amgGM.u != nullptr) {
		delete[] amgGM.u;
		amgGM.u = nullptr;
	}
	if (amgGM.f != nullptr) {
		delete[] amgGM.f;
		amgGM.f = nullptr;
	}
	if (amgGM.ig != nullptr) {
		delete[] amgGM.ig;
		amgGM.ig = nullptr;
	}

	amgGM.nda = -1;
	amgGM.ndf = -1;
	amgGM.ndia = -1;
	amgGM.ndig = -1;
	amgGM.ndja = -1;
	amgGM.ndu = -1;

	for (integer i_7 = 0; i_7 < lmatmax; i_7++) {
		if (matlist[i_7].arr_cp != nullptr) {
			delete[] matlist[i_7].arr_cp;
			matlist[i_7].arr_cp = nullptr;
		}
		if (matlist[i_7].temp_cp != nullptr) {
			delete[] matlist[i_7].temp_cp;
			matlist[i_7].temp_cp = nullptr;
		}
		if (matlist[i_7].arr_lam != nullptr) {
			delete[] matlist[i_7].arr_lam;
			matlist[i_7].arr_lam = nullptr;
		}
		if (matlist[i_7].temp_lam != nullptr) {
			delete[] matlist[i_7].temp_lam;
			matlist[i_7].temp_lam = nullptr;
		}

		if (matlist[i_7].arr_Poisson_ratio != nullptr) {
			delete[] matlist[i_7].arr_Poisson_ratio;
			matlist[i_7].arr_Poisson_ratio = nullptr;
		}
		if (matlist[i_7].temp_Poisson_ratio != nullptr) {
			delete[] matlist[i_7].temp_Poisson_ratio;
			matlist[i_7].temp_Poisson_ratio = nullptr;
		}

		if (matlist[i_7].arr_Young_Module != nullptr) {
			delete[] matlist[i_7].arr_Young_Module;
			matlist[i_7].arr_Young_Module = nullptr;
		}
		if (matlist[i_7].temp_Young_Module != nullptr) {
			delete[] matlist[i_7].temp_Young_Module;
			matlist[i_7].temp_Young_Module = nullptr;
		}

		if (matlist[i_7].arr_beta_t_solid != nullptr) {
			delete[] matlist[i_7].arr_beta_t_solid;
			matlist[i_7].arr_beta_t_solid = nullptr;
		}
		if (matlist[i_7].temp_beta_t_solid != nullptr) {
			delete[] matlist[i_7].temp_beta_t_solid;
			matlist[i_7].temp_beta_t_solid = nullptr;
		}

	}
	delete[] matlist;
	matlist = nullptr;

	if (x_jacoby_buffer != nullptr) {
		// 30 октября 2016. 
		// В seidelsor2 сделан переключатель на метод нижней релаксации К.Г. Якоби.
		// Освобождение памяти из под Jacobi buffer.
		delete[] x_jacoby_buffer;
		x_jacoby_buffer = nullptr;
	}

	if (bvery_big_memory) {
		if (t.database.x != nullptr) {
			free(t.database.x);
			t.database.x = nullptr;
		}
		if (t.database.y != nullptr) {
			free(t.database.y);
			t.database.y = nullptr;
		}
		if (t.database.z != nullptr) {
			free(t.database.z);
			t.database.z = nullptr;
		}
		if (t.database.nvtxcell != nullptr) {
			for (integer i = 0; i <= 7; ++i) {
				delete[] t.database.nvtxcell[i];
			}
			delete[] t.database.nvtxcell;
			t.database.nvtxcell = nullptr;
		}
		if (t.database.ptr != nullptr) {
			if (t.database.ptr[0] != nullptr) {
				delete[] t.database.ptr[0];
			}
			if (t.database.ptr[1] != nullptr) {
				delete[] t.database.ptr[1];
			}
			delete[] t.database.ptr;
			t.database.ptr = nullptr;
		}
	}

	

}// free_global_1

// вычисление центральной точки центрального контрольного объёма
void center_cord3D_ray_tracing(integer iP, int** nvtx, TOCHKA* pa, TOCHKA& p, integer id) {
	// 0 1
	// координаты центра контрольного объёма
	const integer ilink0 = nvtx[0][iP] - 1;
	p.x = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[ilink0].x);
	p.y = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[ilink0].y);
	p.z = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[ilink0].z);
}

// вычисление центральной точки центрального контрольного объёма
void center_cord3D(integer iP, int** nvtx, TOCHKA* pa, TOCHKA &p, integer id) {
	// вычисление центральной точки центрального контрольного объёма
	if (iP < 0) {
		p.x = 0.0;
		p.y = 0.0;
		p.z = 0.0;
		printf("fatall error in center_cord3D !!! ");
		switch (id) {
		case W_SIDE: printf("W"); break;
		case E_SIDE: printf("E"); break;
		case N_SIDE: printf("N"); break;
		case S_SIDE: printf("S"); break;
		case B_SIDE: printf("B"); break;
		case T_SIDE: printf("T"); break;
		case WW_SIDE: printf("WW"); break;
		case EE_SIDE: printf("EE"); break;
		case NN_SIDE: printf("NN"); break;
		case SS_SIDE: printf("SS"); break;
		case BB_SIDE: printf("BB"); break;
		case TT_SIDE: printf("TT"); break;
		default: printf("P"); break;
		}
		system("PAUSE");
	}
	else {


		// 0 1
		// координаты центра контрольного объёма
		p.x = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		p.y = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		p.z = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);

		doublereal 	dx = 0.0, dy = 0.0, dz = 0.0;// размеры контрольного объёма
		dx = pa[nvtx[1][iP] - 1].x - pa[nvtx[0][iP] - 1].x;
		dy = pa[nvtx[2][iP] - 1].y - pa[nvtx[0][iP] - 1].y;
		dz = pa[nvtx[4][iP] - 1].z - pa[nvtx[0][iP] - 1].z;

		/*
		// Эта аномалия (другая нумерация) проявляется на всех АЛИС сетках.
		if (fabs(pa[nvtx[3][iP] - 1].x - pa[nvtx[0][iP] - 1].x) < 1.0e-20) {
			printf("iP=%lld\n",iP);
			printf("anomal: 2<->3\n");
			printf("mast be:\n");
			printf("2 3\n");
			printf("0 1\n");
			system("PAUSE");
		}
		*/

		//if (dz < shorter_length_for_simplificationZ(p.z)) {
		if (dz < 1.0e-36) {
			printf("ERROR Z: INCORRECT NUMERATION IN NVTX MAY BE...\n");
			//printf("error: center cord 3D slipanie po Z. delta < shorter_length_for_simplificationZ\n");
			printf("error: center cord 3D slipanie po Z. delta < 1.0e-36\n");
			std::cout << "iP=" << iP << " dx=" << dx << " dy=" << dy << " dz=" << dz << std::endl;
			std::cout << "nvtx[" << iP << "]: " << nvtx[0][iP] - 1 << " " << nvtx[1][iP] - 1 << " " << nvtx[2][iP] - 1 << " " << nvtx[3][iP] - 1 << " " << nvtx[4][iP] - 1 << " " << nvtx[5][iP] - 1 << " " << nvtx[6][iP] - 1 << " " << nvtx[7][iP] - 1 << std::endl;
			system("PAUSE");
		}
		//if (dx < shorter_length_for_simplificationX(p.x)) {
		if (dx < 1.0e-36) {
		    printf("ERROR X: INCORRECT NUMERATION IN NVTX MAY BE...\n");
			//printf("error: center cord 3D slipanie po X. delta < shorter_length_for_simplificationX\n");
			printf("error: center cord 3D slipanie po X. delta < 1.0e-36\n");
            std::cout << "iP=" << iP << " dx=" << dx << " dy=" << dy << " dz=" << dz << std::endl;
			std::cout << "nvtx[" << iP << "]: " << nvtx[0][iP] - 1 << " " << nvtx[1][iP] - 1 << " " << nvtx[2][iP] - 1 << " " << nvtx[3][iP] - 1 << " " << nvtx[4][iP] - 1 << " " << nvtx[5][iP] - 1 << " " << nvtx[6][iP] - 1 << " " << nvtx[7][iP] - 1 << std::endl;
			system("PAUSE");
		}
		//if (dy < shorter_length_for_simplificationY(p.y)) {
		if (dy < 1.0e-36) {
		    printf("ERROR Y: INCORRECT NUMERATION IN NVTX MAY BE...\n");
			//printf("error: center cord 3D slipanie po Y. delta < shorter_length_for_simplificationY\n");
			printf("error: center cord 3D slipanie po Y. delta <1.0e-36\n");
            std::cout << "iP=" << iP << " dx=" << dx << " dy=" << dy << " dz=" << dz << std::endl;
			std::cout << "nvtx[" << iP << "]: " << nvtx[0][iP] - 1 << " " << nvtx[1][iP] - 1 << " " << nvtx[2][iP] - 1 << " " << nvtx[3][iP] - 1 << " " << nvtx[4][iP] - 1 << " " << nvtx[5][iP] - 1 << " " << nvtx[6][iP] - 1 << " " << nvtx[7][iP] - 1 << std::endl;
			system("PAUSE");
		}
	}

	//printf("%e %e %e\n",dx,dy,dz); // debug GOOD
	//system("PAUSE");	
} // center_cord3D

// вычисляет размеры контрольного объёма
// 23.03.2019 Делаем более сильную проверку и на отрицательный КО тоже.
void volume3D(integer iP, int** nvtx, TOCHKA* pa, doublereal &dx, doublereal &dy, doublereal &dz) {
	// вычисление размеров текущего контрольного объёма:
	if (iP < 0) {
		dx = dy = dz = 0.0;
		printf("ERROR!!! iP out of cabinet: iP < 0.\n");
		system("PAUSE");
		exit(1);
	}
	else {

		TOCHKA p; // Для локальной привязки.
		// координаты центра контрольного объёма
		p.x = 0.5 * (pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
		p.y = 0.5 * (pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
		p.z = 0.5 * (pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);

		//dx = 0.0; dy = 0.0; dz = 0.0;// размеры контрольного объёма
		dx = pa[nvtx[1][iP] - 1].x - pa[nvtx[0][iP] - 1].x;
		dy = pa[nvtx[2][iP] - 1].y - pa[nvtx[0][iP] - 1].y;
		dz = pa[nvtx[4][iP] - 1].z - pa[nvtx[0][iP] - 1].z;
		//if (fabs(pa[nvtx[4][iP] - 1].z - pa[nvtx[0][iP] - 1].z) < 1.0e-36) {
		//if (dz < shorter_length_for_simplificationZ(p.z)) {
		if (dz < 1.0e-36) {
			printf("ERROR 23.03.2019 NVTX perenumeration is ERROR!!!\n");
			//printf("error: volume 3D slipanie po Z. delta < shorter_length_for_simplificationZ\n");
			printf("error: volume 3D slipanie po Z. delta < 1.0e-36\n");
            std::cout << "iP=" << iP << " dx=" << dx << " dy=" << dy << " dz=" << dz << std::endl;
			std::cout << "nvtx[" << iP << "]: " << nvtx[0][iP] - 1 << " " << nvtx[1][iP] - 1 << " " << nvtx[2][iP] - 1 << " " << nvtx[3][iP] - 1 << " " << nvtx[4][iP] - 1 << " " << nvtx[5][iP] - 1 << " " << nvtx[6][iP] - 1 << " " << nvtx[7][iP] - 1 << std::endl;
			system("PAUSE");
		}
		//if (fabs(pa[nvtx[1][iP] - 1].x - pa[nvtx[0][iP] - 1].x) < 1.0e-36) {
		//if (dx < shorter_length_for_simplificationX(p.x)) {
		if (dx < 1.0e-36) {
			printf("ERROR 23.03.2019 NVTX perenumeration is ERROR!!!\n");
			//printf("error: volume 3D slipanie po X. delta < shorter_length_for_simplificationX\n");
			printf("error: volume 3D slipanie po X. delta <  1.0e-36\n");
            std::cout << "iP=" << iP << " dx=" << dx << " dy=" << dy << " dz=" << dz << std::endl;
			std::cout << "nvtx[" << iP << "]: " << nvtx[0][iP] - 1 << " " << nvtx[1][iP] - 1 << " " << nvtx[2][iP] - 1 << " " << nvtx[3][iP] - 1 << " " << nvtx[4][iP] - 1 << " " << nvtx[5][iP] - 1 << " " << nvtx[6][iP] - 1 << " " << nvtx[7][iP] - 1 << std::endl;
			system("PAUSE");
		}
		//if (fabs(pa[nvtx[3][iP] - 1].y - pa[nvtx[0][iP] - 1].y) < 1.0e-36) {
	//if (dy < shorter_length_for_simplificationY(p.y)) {
		if (dy < 1.0e-36) {
			printf("ERROR 23.03.2019 NVTX perenumeration is ERROR!!!\n");
			//printf("error: volume 3D slipanie po Y. delta < shorter_length_for_simplificationY\n");
			printf("error: volume 3D slipanie po Y. delta < 1.0e-36\n");
            std::cout << "iP=" << iP << " dx=" << dx << " dy=" << dy << " dz=" << dz << std::endl;
			std::cout << "nvtx[" << iP << "]: " << nvtx[0][iP] - 1 << " " << nvtx[1][iP] - 1 << " " << nvtx[2][iP] - 1 << " " << nvtx[3][iP] - 1 << " " << nvtx[4][iP] - 1 << " " << nvtx[5][iP] - 1 << " " << nvtx[6][iP] - 1 << " " << nvtx[7][iP] - 1 << std::endl;
			system("PAUSE");
		}
		//printf("%e %e %e\n",dx,dy,dz); // debug GOOD
		//system("PAUSE");
	}
} // volume3D 

// для трассировки лучей.
doublereal volume3D_ray_tracing(integer iP, int** nvtx, TOCHKA* pa) {
	doublereal dx, dy, dz;
	const integer ilink0 = nvtx[0][iP] - 1;
	dx = pa[nvtx[1][iP] - 1].x - pa[ilink0].x;
	dy = pa[nvtx[2][iP] - 1].y - pa[ilink0].y;
	dz = pa[nvtx[4][iP] - 1].z - pa[ilink0].z;
	return (dx * dx + dy * dy + dz * dz);
}

// Сверхбыстрый volume3D 11.04.2021
void volume3D_q(integer iP, int** nvtx, TOCHKA* pa, doublereal& dx, doublereal& dy, doublereal& dz) {
	const integer ilink0 = nvtx[0][iP] - 1;
	dx = pa[nvtx[1][iP] - 1].x - pa[ilink0].x;
	dy = pa[nvtx[2][iP] - 1].y - pa[ilink0].y;
	dz = pa[nvtx[4][iP] - 1].z - pa[ilink0].z;
}

// 1 october 2016
// минимум из двух целых чисел
/*
integer min(integer i1, integer i2) {
	integer ir=i1;
	if (i2<i1) ir=i2;
	return ir;
} // min
*/



// 01.05.2015
typedef struct TQuickSearchBlockid {
	doublereal *x11, *y11, *z11;
	int ix11, iy11, iz11; // размерности в каждом координатном направлении.
	int *ijk_block_id;
	doublereal epsTolx, epsToly, epsTolz;
	int *b_non_prism;
	int lb_non_prism;

	TQuickSearchBlockid() : x11(nullptr), y11(nullptr), z11(nullptr),
		ix11(0), iy11(0), iz11(0), ijk_block_id(nullptr),
		epsTolx(1.0e36), epsToly(1.0e36), epsTolz(1.0e36),
		b_non_prism(nullptr), lb_non_prism(0)
	{
		//x11=nullptr; y11=nullptr; z11=nullptr;
		//ix11 = 0; iy11 = 0; iz11 = 0;
		//ijk_block_id = nullptr;
		//epsTolx = 1.0e36; epsToly = 1.0e36; epsTolz = 1.0e36;
		//b_non_prism=nullptr;
		//lb_non_prism=0;
	}
} QuickSearchBlockid;

QuickSearchBlockid QSBid; 

void init_QSBid(int lb, BLOCK* &b, WALL* &w, int &lw, SOURCE* &s, int &ls) {

#ifdef _OPENMP 
	omp_set_num_threads(number_cores()); // установка числа потоков
#endif


	// Вычисляем количество непризматических объектов QSBid.lb_non_prism и заносим их в отдельный список QSBid.b_non_prism.
	QSBid.lb_non_prism = 0;
	int isum = 0;
#pragma omp parallel for reduction(+ : isum)
	for (int i = 0; i < lb; ++i) {
		if (b[i].g.itypegeom != PRISM) {
			isum += 1;
		}
	}
	QSBid.lb_non_prism = isum;
	if (QSBid.lb_non_prism > 0) {
		if (QSBid.b_non_prism != nullptr) {
			delete[] QSBid.b_non_prism;
			QSBid.b_non_prism = nullptr;
		}
		QSBid.b_non_prism = new int[QSBid.lb_non_prism];
	}
	// Важно записывать их в прямом порядке, обрпатный порядок (приоритетов) учитывается позже по коду.
	int j = 0;
	for (int i = 0; i < lb; ++i) {
		if (b[i].g.itypegeom != PRISM) {
			QSBid.b_non_prism[j] = i;
			j++;
		}
	}

#pragma omp parallel sections 
	{
#pragma omp section
		{
			// X
			QSBid.ix11 = 1;
			QSBid.x11 = new doublereal[QSBid.ix11 + 1]; // число границ по оси x

			QSBid.x11[0] = b[0].g.xS; // начало области
			QSBid.x11[QSBid.ix11] = b[0].g.xE; // конец области

			for (integer i = 1; i < lb; ++i) {

				if (b[i].g.itypegeom == PRISM) {
					doublereal x_1 = b[i].g.xS;
					if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
						addboundary(QSBid.x11, QSBid.ix11, x_1, YZ_PLANE, b, lb, w, lw, s, ls);
					}
					x_1 = b[i].g.xE;
					if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
						addboundary(QSBid.x11, QSBid.ix11, x_1, YZ_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}

			Sort_method(QSBid.x11, QSBid.ix11);
		}
#pragma omp section
		{
			// Y
			QSBid.iy11 = 1;
			QSBid.y11 = new doublereal[QSBid.iy11 + 1]; // число границ по оси y

			QSBid.y11[0] = b[0].g.yS; // начало области
			QSBid.y11[QSBid.iy11] = b[0].g.yE; // конец области

			for (integer i = 1; i < lb; ++i) {

				if (b[i].g.itypegeom == PRISM) {
					doublereal y_1 = b[i].g.yS;
					if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
						addboundary(QSBid.y11, QSBid.iy11, y_1, XZ_PLANE, b, lb, w, lw, s, ls);
					}
					y_1 = b[i].g.yE;
					if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
						addboundary(QSBid.y11, QSBid.iy11, y_1, XZ_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}

			Sort_method(QSBid.y11, QSBid.iy11);
		}
#pragma omp section
		{
			// Z
			QSBid.iz11 = 1;
			QSBid.z11 = new doublereal[QSBid.iz11 + 1]; // число границ по оси z

			QSBid.z11[0] = b[0].g.zS; // начало области
			QSBid.z11[QSBid.iz11] = b[0].g.zE; // конец области

			for (integer i = 1; i < lb; ++i) {

				if (b[i].g.itypegeom == PRISM) {
					doublereal z_1 = b[i].g.zS;
					if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
						addboundary(QSBid.z11, QSBid.iz11, z_1, XY_PLANE, b, lb, w, lw, s, ls);
					}
					z_1 = b[i].g.zE;
					if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
						addboundary(QSBid.z11, QSBid.iz11, z_1, XY_PLANE, b, lb, w, lw, s, ls);
					}
				}
			}

			Sort_method(QSBid.z11, QSBid.iz11);
		}
	}

	delete[] QSBid.ijk_block_id;
	const integer isize_ijk_block_id = (QSBid.ix11 + 1) * (QSBid.iy11 + 1) * (QSBid.iz11 + 1);
	if (print_log_message) {
		printf("isize_ijk_block_id=%lld\n", isize_ijk_block_id);
	}
	QSBid.ijk_block_id = new int[isize_ijk_block_id];
	int* ilink_ijk_block_id = QSBid.ijk_block_id;
#pragma omp parallel for //shared(ilink_ijk_block_id)
	for (integer i_94 = 0; i_94 < isize_ijk_block_id; i_94++) {
		ilink_ijk_block_id[i_94] = lb;
	}
	ilink_ijk_block_id = nullptr;

	// Вычисление допуска.
#pragma omp parallel sections 
	{
#pragma omp section
		{
			QSBid.epsTolx = 1.0e36;
			for (integer i = 0; i < QSBid.ix11; ++i) {
				if (fabs(QSBid.x11[i + 1] - QSBid.x11[i]) < QSBid.epsTolx) {
					QSBid.epsTolx = 0.5*fabs(QSBid.x11[i + 1] - QSBid.x11[i]);
				}
			}
		}
#pragma omp section
		{
			QSBid.epsToly = 1.0e36;
			for (integer i = 0; i < QSBid.iy11; ++i) {
				if (fabs(QSBid.y11[i + 1] - QSBid.y11[i]) < QSBid.epsToly) {
					QSBid.epsToly = 0.5*fabs(QSBid.y11[i + 1] - QSBid.y11[i]);
				}
			}
		}
#pragma omp section
		{
			QSBid.epsTolz = 1.0e36;
			for (integer i = 0; i < QSBid.iz11; ++i) {
				if (fabs(QSBid.z11[i + 1] - QSBid.z11[i]) < QSBid.epsTolz) {
					QSBid.epsTolz = 0.5*fabs(QSBid.z11[i + 1] - QSBid.z11[i]);
				}
			}
		}
	}

	/*
	QSBid.epsToolx = -1.0e36;
	QSBid.epsTooly = -1.0e36;
	QSBid.epsToolz = -1.0e36;
	for (integer i = 0; i < inx; ++i) {
		if (fabs(xpos[i + 1] - xpos[i]) > QSBid.epsToolx) {
			QSBid.epsToolx = 0.5*fabs(xpos[i + 1] - xpos[i]);
		}
	}
	for (integer i = 0; i < iny; ++i) {
		if (fabs(ypos[i + 1] - ypos[i]) > QSBid.epsTooly) {
			QSBid.epsTooly = 0.5*fabs(ypos[i + 1] - ypos[i]);
		}
	}
	for (integer i = 0; i < inz; ++i) {
		if (fabs(zpos[i + 1] - zpos[i]) > QSBid.epsToolz) {
			QSBid.epsToolz = 0.5*fabs(zpos[i + 1] - zpos[i]);
		}
	}
	*/

	Block_indexes* block_indexes = new Block_indexes[lb];
	// Оператор new не требует проверки.
		//if (block_indexes == nullptr) {
			//printf("error in allocation memory for block_indexes in enumerate_volume_improved.\n");
			//system("pause");
			//exit(1);
		//}

	// initialization 30.07.2019
#pragma omp parallel for
	for (integer i_54 = 0; i_54 < lb; i_54++) {
		block_indexes[i_54].iL = -1;
		block_indexes[i_54].iR = -2;
		block_indexes[i_54].jL = -1;
		block_indexes[i_54].jR = -2;
		block_indexes[i_54].kL = -1;
		block_indexes[i_54].kR = -2;
	}

	integer  i=lb-1;
	//integer k;

	// Погрешность бывает абсолютная и относительная.
	// Вещественные числа в ЭВМ представляются с конечной точностью.
	// Лучше использовать относительную погрешность в 0.15%.
	//const doublereal otnositelnaq_tolerance_eps = 0.0015; // 0.15%

	integer iX_one_CELL_count_statistic = 0;
	integer iY_one_CELL_count_statistic = 0;
	integer iZ_one_CELL_count_statistic = 0;



	for (i = lb - 1; i >= 0; i--) {

		//if (b[i].g.itypegeom == 0) {

		// polygon (b[i].g.itypegeom == 2)
		// Значения прямоугольной призмы xS, xE, yS, yE, zS, zE - хранят окаймляющую
		// полигон прямоугольную призму, что позволит проверять принадлежность точки полигону
		// только для ячеек сетки находящихся внутри данной прямоугольной призмы, что сильно 
		// ускоряет обработку.
		//if ((b[i].g.itypegeom == 0) || (b[i].g.itypegeom == 1) || (b[i].g.itypegeom == 2))
		if ((b[i].g.itypegeom == PRISM))
		{

			doublereal x4 = b[i].g.xS;
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == XZ))) {
				//x4 = b[i].g.xC - b[i].g.R_out_cyl;
			//}
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == YZ))) {
				//if (b[i].g.Hcyl > 0.0) {
					//x4 = b[i].g.xC;
				//}
				//else {
					//x4 = b[i].g.xC + b[i].g.Hcyl;
				//}
			//}
			bool bfound = false;
			/*
			for (j = 0; j <= QSBid.ix11; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (QSBid.x11[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i].iL = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(QSBid.x11[j] - x4) < shorter_length_for_simplificationX(x4)) {
						block_indexes[i].iL = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (j = 0; j <= QSBid.ix11; ++j) {
					if ((!bfound) && (fabs(x4 - QSBid.x11[j]) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.

#pragma omp critical (CRITICAL_BLOCK_BODY1)
						{
							block_indexes[i].iL = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			if (!bfound) {
#pragma omp parallel for
				for (j = 0; j <= QSBid.ix11; ++j) {
					if ((!bfound) && (x4 > QSBid.x11[j])) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY2)
						{
							block_indexes[i].iL = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			x4 = b[i].g.xE;
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == XZ))) {
				//x4 = b[i].g.xC + b[i].g.R_out_cyl;
			//}
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == YZ))) {
				//if (b[i].g.Hcyl > 0.0) {
					//x4 = b[i].g.xC + b[i].g.Hcyl;
				//}
				//else {
					//x4 = b[i].g.xC;
				//}
			//}
			bfound = false;
			/*
			for (j = 0; j <= QSBid.ix11; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (QSBid.x11[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i].iR = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(QSBid.x11[j] - x4) < shorter_length_for_simplificationX(x4)) {
						block_indexes[i].iR = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (j = QSBid.ix11; j >= 0; j--) {
					if ((!bfound) && (fabs(x4 - QSBid.x11[j])<shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY3)
						{
							block_indexes[i].iR = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			if (!bfound) {
#pragma omp parallel for
				for (j = QSBid.ix11; j >= 0; j--) {
					if ((!bfound) && (x4 <= QSBid.x11[j])) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY4)
						{
							block_indexes[i].iR = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			x4 = b[i].g.yS;
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == YZ))) {
				//x4 = b[i].g.yC - b[i].g.R_out_cyl;
			//}
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XZ))) {
				//if (b[i].g.Hcyl > 0.0) {
					//x4 = b[i].g.yC;
				//}
				//else {
					//x4 = b[i].g.yC + b[i].g.Hcyl;
				//}
			//}
			bfound = false;
			/*
			for (j = 0; j <= QSBid.iy11; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (QSBid.y11[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i].jL = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(QSBid.y11[j] - x4) < shorter_length_for_simplificationY(x4)) {
						block_indexes[i].jL = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (j = 0; j <= QSBid.iy11; ++j) {
					if ((!bfound) && (fabs(x4 - QSBid.y11[j])<shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY5)
						{
							block_indexes[i].jL = j;
							bfound = true;
						//	break;
						}
					}
				}
			}
			if (!bfound) {
#pragma omp parallel for
				for (j = 0; j <= QSBid.iy11; ++j) {
					if ((!bfound) && (x4 > QSBid.y11[j])) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY6)
						{
							block_indexes[i].jL = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			x4 = b[i].g.yE;
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == YZ))) {
				//x4 = b[i].g.yC + b[i].g.R_out_cyl;
			//}
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XZ))) {
				//if (b[i].g.Hcyl > 0.0) {
					//x4 = b[i].g.yC + b[i].g.Hcyl;
				//}
				//else {
					//x4 = b[i].g.yC;
				//}
			//}
			bfound = false;
			/*
			for (j = 0; j <= QSBid.iy11; ++j) {

				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (QSBid.y11[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i].jR = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(QSBid.y11[j] - x4) < shorter_length_for_simplificationY(x4)) {
						block_indexes[i].jR = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (j = QSBid.iy11; j >= 0; j--) {
					if ((!bfound) && (fabs(x4 - QSBid.y11[j])<shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY7)
						{
							block_indexes[i].jR = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			if (!bfound) {
#pragma omp parallel for
				for (j = QSBid.iy11; j >= 0; j--) {
					if ((!bfound) && (x4 <= QSBid.y11[j])) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY8)
						{
							block_indexes[i].jR = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			x4 = b[i].g.zS;
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XZ) || (b[i].g.iPlane == YZ))) {
				//x4 = b[i].g.zC - b[i].g.R_out_cyl;
			//}
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XY))) {
				//if (b[i].g.Hcyl > 0.0) {
					//x4 = b[i].g.zC;
				//}
				//else {
					//x4 = b[i].g.zC + b[i].g.Hcyl;
				//}
			//}
			bfound = false;
			/*
			for (j = 0; j <= QSBid.iz11; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (QSBid.z11[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i].kL = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(QSBid.z11[j] - x4) < shorter_length_for_simplificationZ(x4)) {
						block_indexes[i].kL = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (j = 0; j <= QSBid.iz11; ++j) {
					if ((!bfound) && (fabs(x4 - QSBid.z11[j])<shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY9)
						{
							block_indexes[i].kL = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			if (!bfound) {
#pragma omp parallel for
				for (j = 0; j <= QSBid.iz11; ++j) {
					if ((!bfound) && (x4 > QSBid.z11[j])) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY10)
						{
							block_indexes[i].kL = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			x4 = b[i].g.zE;
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XZ) || (b[i].g.iPlane == YZ))) {
				//x4 = b[i].g.zC + b[i].g.R_out_cyl;
			//}
			//if ((b[i].g.itypegeom == 1) && ((b[i].g.iPlane == XY))) {
				//if (b[i].g.Hcyl > 0.0) {
					//x4 = b[i].g.zC + b[i].g.Hcyl;
				//}
				//else {
					//x4 = b[i].g.zC;
				//}
			//}
			bfound = false;
			/*
			for (j = 0; j <= QSBid.iz11; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (QSBid.z11[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i].kR = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(QSBid.z11[j] - x4) < shorter_length_for_simplificationZ(x4)) {
						block_indexes[i].kR = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (j = QSBid.iz11; j >= 0; j--) {
					if ((!bfound) && (fabs(x4 - QSBid.z11[j])<shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY11)
						{
							block_indexes[i].kR = j;
							bfound = true;
							//break;
						}
					}
				}
			}
			if (!bfound) {
#pragma omp parallel for
				for (j = QSBid.iz11; j >= 0; j--) {
					if ((!bfound) && (x4 <= QSBid.z11[j])) {
						// Нет точного совпаднения первая встреча.
#pragma omp critical (CRITICAL_BLOCK_BODY12)
						{
							block_indexes[i].kR = j;
							bfound = true;
							//break;
						}
					}
				}
			}

			if ((block_indexes[i].iL>=0)&&
				(block_indexes[i].iR>=0)&&
				(block_indexes[i].iL >= block_indexes[i].iR)) {
				printf("init_QSBid function\n");
				printf("violation of the order\n");
				std::cout << "i=" << i << " iL=" << block_indexes[i].iL << " iR=" << block_indexes[i].iR << " delta=" << shorter_length_for_simplificationX(QSBid.x11[block_indexes[i].iR], b, lb, w, lw, s, ls) << std::endl;
				system("pause");
			}

			if ((block_indexes[i].iL >= 0) &&
				(block_indexes[i].iR >= 0) &&
				(block_indexes[i].iL +1 == block_indexes[i].iR)) {
				// Всего одна клетка на блок.
				iX_one_CELL_count_statistic++;
			}

			if ((block_indexes[i].jL >= 0) &&
				(block_indexes[i].jR >= 0) &&
				(block_indexes[i].jL >= block_indexes[i].jR)) {
				printf("init_QSBid function\n");
				printf("violation of the order\n");
				std::cout << "i=" << i << " jL=" << block_indexes[i].jL << " jR=" <<  block_indexes[i].jR << " delta=" << shorter_length_for_simplificationY(QSBid.y11[block_indexes[i].jR], b, lb, w, lw, s, ls) << std::endl;
				system("pause");
			}


			if ((block_indexes[i].jL >= 0) &&
				(block_indexes[i].jR >= 0) &&
				(block_indexes[i].jL + 1 == block_indexes[i].jR)) {
				// Всего одна клетка на блок.
				iY_one_CELL_count_statistic++;
			}

			if ((block_indexes[i].kL >= 0) &&
				(block_indexes[i].kR >= 0) &&
				(block_indexes[i].kL >= block_indexes[i].kR)) {
				printf("init_QSBid function\n");
				printf("violation of the order\n");
				printf("i=%lld kL=%lld kR=%lld delta=%e\n", i, block_indexes[i].kL,
					block_indexes[i].kR, shorter_length_for_simplificationZ(QSBid.z11[block_indexes[i].kR], b, lb, w, lw, s, ls));
				system("pause");
			}


			if ((block_indexes[i].kL >= 0) &&
				(block_indexes[i].kR >= 0) &&
				(block_indexes[i].kL + 1 == block_indexes[i].kR)) {
				// Всего одна клетка на блок.
				iZ_one_CELL_count_statistic++;
			}

			if (0&&(i == 116)) {
				// debug
				printf("i=%lld iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
					block_indexes[i].iL,
					block_indexes[i].iR, block_indexes[i].jL,
					block_indexes[i].jR, block_indexes[i].kL,
					block_indexes[i].kR);
				printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[i].g.xS, b[i].g.xE,
					b[i].g.yS, b[i].g.yE, b[i].g.zS, b[i].g.zE);
				system("pause");
			}

			if ((block_indexes[i].iL < 0) ||
				(block_indexes[i].iR < 0) ||
				(block_indexes[i].jL < 0) ||
				(block_indexes[i].jR < 0) ||
				(block_indexes[i].kL < 0) ||
				(block_indexes[i].kR < 0)) {
				printf("init_QSBid function\n");
				printf("i=%lld iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
					block_indexes[i].iL,
					block_indexes[i].iR, block_indexes[i].jL,
					block_indexes[i].jR, block_indexes[i].kL,
					block_indexes[i].kR);
				printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[i].g.xS, b[i].g.xE,
					b[i].g.yS, b[i].g.yE, b[i].g.zS, b[i].g.zE);
				printf("cabinet: xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[0].g.xS, b[0].g.xE,
					b[0].g.yS, b[0].g.yE, b[0].g.zS, b[0].g.zE);
				printf("ERROR: may be your geometry out of cabinet...\n");
				system("pause");
			}


		}
		else {
			block_indexes[i].iL = -1;
			block_indexes[i].iR = -2;
			block_indexes[i].jL = -1;
			block_indexes[i].jR = -2;
			block_indexes[i].kL = -1;
			block_indexes[i].kR = -2;
		}
	}


	if ((iX_one_CELL_count_statistic > 0) || (iY_one_CELL_count_statistic > 0) || (iZ_one_CELL_count_statistic > 0)) {
		printf("WARNING ONE CELL ON BLOCK...\n");
		std::cout << "STATISTICS: X=" << iX_one_CELL_count_statistic << " Y=" << iY_one_CELL_count_statistic <<" Z=" << iZ_one_CELL_count_statistic <<"\n";
		// Работает и при такой ситуации. Проверено на одной из моделей.
		//system("pause");
	}


	// block_indexes подготовлен.

	bool* bvisit = nullptr;
	bvisit = new bool[isize_ijk_block_id];

	// Количество проходов существенно сократилось и в итоге это приводит к существенному
		// увеличению быстродействия.
	int m7 = lb - 1, m8;

#pragma omp parallel for
	for (integer iP = 0; iP < isize_ijk_block_id; ++iP) {
		bvisit[iP] = false;
	}



	// integer iP_k1, iP_j1, iP;
	for (m8 = lb - 1; m8 >= 0; m8--) {
		m7 = m8;
		if (b[m8].g.itypegeom == PRISM)
		{
			//#pragma omp parallel for
			for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++)
			{
				integer iP_k1 =k1 * (QSBid.ix11 + 1)*(QSBid.iy11 + 1);
				if ((k1 < 0) || (k1 > QSBid.iz11)) {
					// ERROR
					printf("ERROR PRISM\n");
					printf("inx=%d iny=%d inz=%d \n", QSBid.ix11, QSBid.iy11, QSBid.iz11);
					printf("k1=%lld \n",  k1);
					printf("iP_k1=%lld m8=%d", iP_k1, m8);
					printf("iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
					system("PAUSE");
				}
				for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++)
				{
					integer iP_j1=j1 * (QSBid.ix11 + 1) +iP_k1;
					if ((j1 < 0) || (j1 > QSBid.iy11)) {
							// ERROR
							printf("ERROR PRISM\n");
							printf("inx=%d iny=%d inz=%d \n", QSBid.ix11, QSBid.iy11, QSBid.iz11);
							printf("j1=%lld k1=%lld \n", j1, k1);
							printf("iP_j1=%lld m8=%d", iP_j1, m8);
							printf("iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
							system("PAUSE");
					}					
					for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++)
					{

						integer iP = i1 +  iP_j1;

						if ((i1 < 0) || (i1 > QSBid.ix11)) {
							// ERROR
							printf("ERROR PRISM\n");
							printf("inx=%d iny=%d inz=%d \n", QSBid.ix11, QSBid.iy11, QSBid.iz11);
							printf("i1=%lld j1=%lld k1=%lld \n", i1, j1, k1);
							printf("iP=%lld m8=%d", iP, m8);
							printf("iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
							system("PAUSE");
						}

						if (bvisit[iP] == false)
						{

							bvisit[iP] = true;

							QSBid.ijk_block_id[iP] = m8;
						}
					}
				}
			}
			//m7--;
		}
	}

	
	// Так делать нельзя т.к. мы заполняли только на призмах и не рассматривали цилиндры и полигоны.
	// Поэтому какие-то позиции из ijk_block_id[iP] останутся равными lb т.е. незаполненными.
	// Это нормально. НО!!! здесь мы добавили проверку if (bvisit[iP]  )  и поэтому гарантируется
	// Что блоки только призматические и только.
#pragma omp parallel for
	for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++)
	{
		integer iP_k1=k1 * (QSBid.ix11 + 1)*(QSBid.iy11 + 1);
		for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++)
		{
			integer iP_j1=j1 * (QSBid.ix11 + 1) +iP_k1;
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++)
			{
				integer iP = i1 +  iP_j1;
				if (bvisit[iP]  ) {
					// Если мы ранее посещали данную ячейку при обходе призматических объектов.
					if (QSBid.ijk_block_id[iP] == lb) {
						// Мы что-то пропустили и из-за этого возможен сбой в дальнейшем.
						// исправляем так чтобы сбоя не было 28.07.2019
						TOCHKA p;
						p.x = 0.5 * (QSBid.x11[i1] + QSBid.x11[i1 + 1]);
						p.y = 0.5 * (QSBid.y11[j1] + QSBid.y11[j1 + 1]);
						p.z = 0.5 * (QSBid.z11[k1] + QSBid.z11[k1 + 1]);
						// Лобовой надежный метод, правда очень медленный.
						// Чтобы работало быстро таких аномальных точек должен быть 
						// небольшой процент от общего числа.
						QSBid.ijk_block_id[iP] = myisblock_id_stab(lb, b, p);

					}
				}
			}
        }
    }
	

	delete[] bvisit;
	bvisit = nullptr;

	delete[] block_indexes;
	block_indexes = nullptr;

	

}

void free_QSBid() {
	QSBid.lb_non_prism = 0;

	if (QSBid.b_non_prism!=nullptr) {
	    delete[] QSBid.b_non_prism;
	    QSBid.b_non_prism = nullptr;
    }

	if (QSBid.x11 != nullptr) {
		delete[] QSBid.x11;
		QSBid.x11 = nullptr;
	}
	if (QSBid.y11 != nullptr) {
		delete[] QSBid.y11;
		QSBid.y11 = nullptr;
	}
	if (QSBid.z11 != nullptr) {
		delete[] QSBid.z11;
		QSBid.z11 = nullptr;
	}

	if (QSBid.ijk_block_id != nullptr) {
		delete[] QSBid.ijk_block_id;
		QSBid.ijk_block_id = nullptr;
	}

	QSBid.epsTolx = 1.0e36;
	QSBid.epsToly = 1.0e36;
	QSBid.epsTolz = 1.0e36;

	QSBid.ix11 = 0;
	QSBid.iy11 = 0;
	QSBid.iz11 = 0;
}

// Какому прямоугольному параллелепипеду принадлежит точка p.
int myisblock_id_PRISM_only(int lb, BLOCK* &b, TOCHKA p) {
	// return i_vacant + (inx + 1)*j_vacant + (inx + 1)*(iny + 1)*k_vacant;
	integer ikey = hash_key_alice33Q(QSBid.ix11, QSBid.iy11, QSBid.iz11, QSBid.x11, QSBid.y11, QSBid.z11, p, QSBid.epsTolx, QSBid.epsToly, QSBid.epsTolz);
	int ib= QSBid.ijk_block_id[ikey];
	if ((ib < 0) || (ib >= lb)) {
		printf("error in function myisblock_id_PRISM_only.\n");
		printf("p.x=%e p.y=%e p.z=%e ikey==%lld\n",p.x,p.y,p.z,ikey);
		printf("ikey=%lld size=%lld ib=%d\n",ikey, static_cast<integer>((QSBid.ix11 + 1)*(QSBid.iy11 + 1)*(QSBid.iz11 + 1)), ib);
		printf("lx=%d ly=%d iz=%d\n", QSBid.ix11, QSBid.iy11, QSBid.iz11);
		
		// Вызываем очень медленный наивный алгоритм.
		ib = lb;
		for (int i_1 = lb - 1; i_1 >= 0; i_1--) {
			if (b[i_1].g.itypegeom == PRISM) {
				if ((b[i_1].g.xS < p.x) && (b[i_1].g.xE > p.x) && (b[i_1].g.yS < p.y) && (b[i_1].g.yE > p.y) && (b[i_1].g.zS < p.z) && (b[i_1].g.zE > p.z)) {
					ib = i_1;
					// Если блок найден то сканирование сразу прекращается.
					break;
				}
			}
		}
		printf("ib==%d lb=%d\n",ib,lb);
		//system("PAUSE");
	}
	return ib;
} // myisblock_id_PRISM_only


// Хотел использовать для поиска соседей при асемблесах.
// Наивная реализация (очень медленная) функции whot_is_block.
// Не тестировал. 18,11,2021
int myisblock_id_naiv(int lb, BLOCK*& b, TOCHKA p) {

	int ib = 0;
	int i = 0, k = lb;

   // цикл по всем блокам
   for (i = lb - 1; i >= 0; i--) {
	   if (b[i].g.itypegeom == PRISM) {
		   // Мы уже выполнили ранее быстрое вычисление принадлежности точки призматическому объекту.

		   // Prism
		   // Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
		   //то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что
		   // приведет к неверному типу блока у данной ячейки.
		   if ((p.x >= b[i].g.xS) && (p.x <= b[i].g.xE) && (p.y >= b[i].g.yS) && (p.y <= b[i].g.yE) && (p.z >= b[i].g.zS) && (p.z <= b[i].g.zE)) {
			   k = i;
			   // Нашли и сразу завершили проверку в случае успеха.
			   goto OUTOF_IN_MODEL_TEMP51;
		   }
	   }
		if (b[i].g.itypegeom == CAD_STL) {
			  // 14.11.2020

			  int k_loc23 = -1;
			  if (b[i].g.in_CAD_STL_check(p, k_loc23, i)) {
				  k = i;
				  // Нашли и сразу завершили проверку в случае успеха.
				  goto OUTOF_IN_MODEL_TEMP51;
			  }
		  }
	   if (b[i].g.itypegeom == POLYGON) {

		   if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {

			   bool bfound = false;
			   // определяет принадлежность точки полигону.
			   bfound = in_polygon(p, b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2, k, i);
			   if (bfound) {
				   // Нашли и сразу завершили проверку в случае успеха.
				   goto OUTOF_IN_MODEL_TEMP51;
			   }
		   }

	   }
	   if (b[i].g.itypegeom == CYLINDER) {
		   // Cylinder
		   switch (b[i].g.iPlane) {
		   case XY_PLANE:
			   if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
				   if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
						   k = i;
						   // Нашли и сразу завершили проверку в случае успеха.
						   goto OUTOF_IN_MODEL_TEMP51;
					   }
				   }
			   }
			   else {
				   if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
						   if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) > b[i].g.R_in_cyl) {
							   k = i;
							   // Нашли и сразу завершили проверку в случае успеха.
							   goto OUTOF_IN_MODEL_TEMP51;
						   }
					   }
				   }
			   }
			   break;
			   case XZ_PLANE:
			   if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
				   if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
						   k = i;
						   // Нашли и сразу завершили проверку в случае успеха.
						   goto OUTOF_IN_MODEL_TEMP51;
					   }
				   }
			   }
			   else {
				   if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
						   if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
							   k = i;
							   // Нашли и сразу завершили проверку в случае успеха.
							   goto OUTOF_IN_MODEL_TEMP51;
						   }
					   }
				   }
			   }
			   break;
		   case YZ_PLANE:
			   if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
				   if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
						   k = i;
						   // Нашли и сразу завершили проверку в случае успеха.
						   goto OUTOF_IN_MODEL_TEMP51;
					   }
				   }
			   }
			   else {
				   if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
						   if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
							   k = i;
							   // Нашли и сразу завершили проверку в случае успеха.
							   goto OUTOF_IN_MODEL_TEMP51;
						   }
					   }
				   }
			   }
			   break;
		   }
	   }
   }
   
OUTOF_IN_MODEL_TEMP51:

  

   ib = k;

   if (ib == lb) {
	   printf("FATALL ERROR!!! ib==lb\n");
	   system("PAUSE");
   }

   return ib;
}


// Полнофункциональная медленная версия работающая и для цилиндров и для полигонов,
// а также для прямоугольных параллелепипедов.
// 01.05.2019
int myisblock_id(int lb, BLOCK* &b, TOCHKA p) {
	int ib = 0;
	//for (integer i_1 = lb - 1; i_1 >= 0; i_1--) {
		//if ((b[i_1].g.xS < p.x) && (b[i_1].g.xE > p.x) && (b[i_1].g.yS < p.y) && (b[i_1].g.yE > p.y) && (b[i_1].g.zS < p.z) && (b[i_1].g.zE > p.z)) {
			//ib = i_1;
			// Если блок найден то сканирование сразу прекращается.
			//break;
		//}
	//}

	int i = 0, k = 0;

	// Параллельный алгоритм существенно медленней однопоточного.
	// Однопоточный вариант полностью пересмотрен и использует двоичный поиск и хеш таблицы.
	// Откажемся от многопоточного варианта т.к. он неэффективен.
	//integer lb4 = lb / 4; // Это очень плохое деление пригодное только в случае,
	// если все блоки равноправны: содержат примерно одинаковое число ячеек.

/*
#ifdef _OPENMP
	integer k1 = 0, k2 = 0, k3 = 0, k4 = 0;
#pragma omp parallel shared(k1,k2,k3,k4) num_threads(4)
	{
#ifdef _OPENMP 
				int tid = omp_get_thread_num();
#else
				int tid = 0;
#endif
		if (tid == 0) {
			// цикл по всем блокам
			//for (integer i1 = 0; i1<lb4; i1++) {
			for (integer i1 = 0; i1 < lb; i1 += 4) {
				if (b[i1].g.itypegeom == PRISM) {
					// Prism
					// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
					//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
					// приведет к неверному типу блока у данной ячейки.
					if ((p.x >= b[i1].g.xS) && (p.x <= b[i1].g.xE) && (p.y >= b[i1].g.yS) && (p.y <= b[i1].g.yE) && (p.z >= b[i1].g.zS) && (p.z <= b[i1].g.zE)) {
						k1 = i1;
					}
				}
				if (b[i1].g.itypegeom == POLYGON) {

					if ((p.x > b[i1].g.xS) && (p.x < b[i1].g.xE) && (p.y > b[i1].g.yS) && (p.y < b[i1].g.yE) && (p.z > b[i1].g.zS) && (p.z < b[i1].g.zE)) {
						// определяет принадлежность точки полигону.
						in_polygon(p, b[i1].g.nsizei, b[i1].g.xi, b[i1].g.yi, b[i1].g.zi, b[i1].g.hi, b[i1].g.iPlane_obj2, k1, i1);
					}

				}

				if (b[i1].g.itypegeom == CAD_STL) {

					
					integer k_loc=-1;
					if (b[i1].g.in_CAD_STL_check(p,k_loc,i1)) {
					   k1=i1;
					}

				}

				if (b[i1].g.itypegeom == CYLINDER) {
					// Cylinder
					switch (b[i1].g.iPlane) {
					case XY:
						if (fabs(b[i1].g.R_in_cyl) < 1.0e-36) {
							if ((p.z > b[i1].g.zC) && (p.z < b[i1].g.zC + b[i1].g.Hcyl)) {
								if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.yC - p.y)*(b[i1].g.yC - p.y)) < b[i1].g.R_out_cyl) {
									k1 = i1;
								}
							}
						}
						else {
							if ((p.z > b[i1].g.zC) && (p.z < b[i1].g.zC + b[i1].g.Hcyl)) {
								if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.yC - p.y)*(b[i1].g.yC - p.y)) < b[i1].g.R_out_cyl) {
									if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.yC - p.y)*(b[i1].g.yC - p.y)) > b[i1].g.R_in_cyl) {
										k1 = i1;
									}
								}
							}
						}
						break;
					case XZ:
						if (fabs(b[i1].g.R_in_cyl) < 1.0e-36) {
							if ((p.y > b[i1].g.yC) && (p.y < b[i1].g.yC + b[i1].g.Hcyl)) {
								if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) < b[i1].g.R_out_cyl) {
									k1 = i1;
								}
							}
						}
						else {
							if ((p.y > b[i1].g.yC) && (p.y < b[i1].g.yC + b[i1].g.Hcyl)) {
								if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) < b[i1].g.R_out_cyl) {
									if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) > b[i1].g.R_in_cyl) {
										k1 = i1;
									}
								}
							}
						}
						break;
					case YZ:
						if (fabs(b[i1].g.R_in_cyl) < 1.0e-36) {
							if ((p.x > b[i1].g.xC) && (p.x < b[i1].g.xC + b[i1].g.Hcyl)) {
								if (sqrt((b[i1].g.yC - p.y)*(b[i1].g.yC - p.y) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) < b[i1].g.R_out_cyl) {
									k1 = i1;
								}
							}
						}
						else {
							if ((p.x > b[i1].g.xC) && (p.x < b[i1].g.xC + b[i1].g.Hcyl)) {
								if (sqrt((b[i1].g.yC - p.y)*(b[i1].g.yC - p.y) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) < b[i1].g.R_out_cyl) {
									if (sqrt((b[i1].g.yC - p.y)*(b[i1].g.yC - p.y) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) > b[i1].g.R_in_cyl) {
										k1 = i1;
									}
								}
							}
						}
						break;
					}
				}
			}
		}
		if (tid == 1) {
			// цикл по всем блокам
			//for (integer i2 = lb4; i2<2 * lb4; i2++) {
			for (integer i2 = 1; i2 < lb; i2 += 4) {
				if (b[i2].g.itypegeom == PRISM) {
					// Prism
					// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
					//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
					// приведет к неверному типу блока у данной ячейки.
					if ((p.x >= b[i2].g.xS) && (p.x <= b[i2].g.xE) && (p.y >= b[i2].g.yS) && (p.y <= b[i2].g.yE) && (p.z >= b[i2].g.zS) && (p.z <= b[i2].g.zE)) {
						k2 = i2;
					}
				}
				if (b[i2].g.itypegeom == POLYGON) {

					if ((p.x > b[i2].g.xS) && (p.x < b[i2].g.xE) && (p.y > b[i2].g.yS) && (p.y < b[i2].g.yE) && (p.z > b[i2].g.zS) && (p.z < b[i2].g.zE)) {
						// определяет принадлежность точки полигону.
						in_polygon(p, b[i2].g.nsizei, b[i2].g.xi, b[i2].g.yi, b[i2].g.zi, b[i2].g.hi, b[i2].g.iPlane_obj2, k2, i2);
					}

				}

				if (b[i2].g.itypegeom == CAD_STL) {


					integer k_loc=-1;
					if (b[i2].g.in_CAD_STL_check(p,k_loc,i2)) {
					   k2=i2;
					}

				}

				if (b[i2].g.itypegeom == CYLINDER) {
					// Cylinder
					switch (b[i2].g.iPlane) {
					case XY:
						if (fabs(b[i2].g.R_in_cyl) < 1.0e-36) {
							if ((p.z > b[i2].g.zC) && (p.z < b[i2].g.zC + b[i2].g.Hcyl)) {
								if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.yC - p.y)*(b[i2].g.yC - p.y)) < b[i2].g.R_out_cyl) {
									k2 = i2;
								}
							}
						}
						else {
							if ((p.z > b[i2].g.zC) && (p.z < b[i2].g.zC + b[i2].g.Hcyl)) {
								if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.yC - p.y)*(b[i2].g.yC - p.y)) < b[i2].g.R_out_cyl) {
									if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.yC - p.y)*(b[i2].g.yC - p.y)) > b[i2].g.R_in_cyl) {
										k2 = i2;
									}
								}
							}
						}
						break;
					case XZ:
						if (fabs(b[i2].g.R_in_cyl) < 1.0e-36) {
							if ((p.y > b[i2].g.yC) && (p.y < b[i2].g.yC + b[i2].g.Hcyl)) {
								if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) < b[i2].g.R_out_cyl) {
									k2 = i2;
								}
							}
						}
						else {
							if ((p.y > b[i2].g.yC) && (p.y < b[i2].g.yC + b[i2].g.Hcyl)) {
								if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) < b[i2].g.R_out_cyl) {
									if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) > b[i2].g.R_in_cyl) {
										k2 = i2;
									}
								}
							}
						}
						break;
					case YZ:
						if (fabs(b[i2].g.R_in_cyl) < 1.0e-36) {
							if ((p.x > b[i2].g.xC) && (p.x < b[i2].g.xC + b[i2].g.Hcyl)) {
								if (sqrt((b[i2].g.yC - p.y)*(b[i2].g.yC - p.y) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) < b[i2].g.R_out_cyl) {
									k2 = i2;
								}
							}
						}
						else {
							if ((p.x > b[i2].g.xC) && (p.x < b[i2].g.xC + b[i2].g.Hcyl)) {
								if (sqrt((b[i2].g.yC - p.y)*(b[i2].g.yC - p.y) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) < b[i2].g.R_out_cyl) {
									if (sqrt((b[i2].g.yC - p.y)*(b[i2].g.yC - p.y) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) > b[i2].g.R_in_cyl) {
										k2 = i2;
									}
								}
							}
						}
						break;
					}
				}
			}
		}
		if (tid == 2) {
			// цикл по всем блокам
			//for (integer i3 = 2 * lb4; i3<3 * lb4; i3++) {
			for (integer i3 = 2; i3 < lb; i3 += 4) {
				if (b[i3].g.itypegeom == PRISM) {
					// Prism
					// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
					//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
					// приведет к неверному типу блока у данной ячейки.
					if ((p.x >= b[i3].g.xS) && (p.x <= b[i3].g.xE) && (p.y >= b[i3].g.yS) && (p.y <= b[i3].g.yE) && (p.z >= b[i3].g.zS) && (p.z <= b[i3].g.zE)) {
						k3 = i3;
					}
				}
				if (b[i3].g.itypegeom == POLYGON) {

					if ((p.x > b[i3].g.xS) && (p.x < b[i3].g.xE) && (p.y > b[i3].g.yS) && (p.y < b[i3].g.yE) && (p.z > b[i3].g.zS) && (p.z < b[i3].g.zE)) {
						// определяет принадлежность точки полигону.
						in_polygon(p, b[i3].g.nsizei, b[i3].g.xi, b[i3].g.yi, b[i3].g.zi, b[i3].g.hi, b[i3].g.iPlane_obj2, k3, i3);
					}

				}

				if (b[i3].g.itypegeom == CAD_STL) {

					integer k_loc=-1;
					if (b[i3].g.in_CAD_STL_check(p,k_loc,i3)) {
					   k3=i3;
					}

				}

				if (b[i3].g.itypegeom == CYLINDER) {
					// Cylinder
					switch (b[i3].g.iPlane) {
					case XY:
						if (fabs(b[i3].g.R_in_cyl) < 1.0e-36) {
							if ((p.z > b[i3].g.zC) && (p.z < b[i3].g.zC + b[i3].g.Hcyl)) {
								if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.yC - p.y)*(b[i3].g.yC - p.y)) < b[i3].g.R_out_cyl) {
									k3 = i3;
								}
							}
						}
						else {
							if ((p.z > b[i3].g.zC) && (p.z < b[i3].g.zC + b[i3].g.Hcyl)) {
								if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.yC - p.y)*(b[i3].g.yC - p.y)) < b[i3].g.R_out_cyl) {
									if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.yC - p.y)*(b[i3].g.yC - p.y)) > b[i3].g.R_in_cyl) {
										k3 = i3;
									}
								}
							}
						}
						break;
					case XZ:
						if (fabs(b[i3].g.R_in_cyl) < 1.0e-36) {
							if ((p.y > b[i3].g.yC) && (p.y < b[i3].g.yC + b[i3].g.Hcyl)) {
								if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) < b[i3].g.R_out_cyl) {
									k3 = i3;
								}
							}
						}
						else {
							if ((p.y > b[i3].g.yC) && (p.y < b[i3].g.yC + b[i3].g.Hcyl)) {
								if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) < b[i3].g.R_out_cyl) {
									if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) > b[i3].g.R_in_cyl) {
										k3 = i3;
									}
								}
							}
						}
						break;
					case YZ:
						if (fabs(b[i3].g.R_in_cyl) < 1.0e-36) {
							if ((p.x > b[i3].g.xC) && (p.x < b[i3].g.xC + b[i3].g.Hcyl)) {
								if (sqrt((b[i3].g.yC - p.y)*(b[i3].g.yC - p.y) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) < b[i3].g.R_out_cyl) {
									k3 = i3;
								}
							}
						}
						else {
							if ((p.x > b[i3].g.xC) && (p.x < b[i3].g.xC + b[i3].g.Hcyl)) {
								if (sqrt((b[i3].g.yC - p.y)*(b[i3].g.yC - p.y) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) < b[i3].g.R_out_cyl) {
									if (sqrt((b[i3].g.yC - p.y)*(b[i3].g.yC - p.y) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) > b[i3].g.R_in_cyl) {
										k3 = i3;
									}
								}
							}
						}
						break;
					}
				}
			}
		}
		if (tid == 3) {
			// цикл по всем блокам
			//for (integer i4 = 3 * lb4; i4<lb; i4++) {
			for (integer i4 = 3; i4 < lb; i4 += 4) {
				if (b[i4].g.itypegeom == PRISM) {
					// Prism
					// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
					//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
					// приведет к неверному типу блока у данной ячейки.
					if ((p.x >= b[i4].g.xS) && (p.x <= b[i4].g.xE) && (p.y >= b[i4].g.yS) && (p.y <= b[i4].g.yE) && (p.z >= b[i4].g.zS) && (p.z <= b[i4].g.zE)) {
						k4 = i4;
					}
				}
				if (b[i4].g.itypegeom == POLYGON) {

					if ((p.x > b[i4].g.xS) && (p.x < b[i4].g.xE) && (p.y > b[i4].g.yS) && (p.y < b[i4].g.yE) && (p.z > b[i4].g.zS) && (p.z < b[i4].g.zE)) {
						// определяет принадлежность точки полигону.
						in_polygon(p, b[i4].g.nsizei, b[i4].g.xi, b[i4].g.yi, b[i4].g.zi, b[i4].g.hi, b[i4].g.iPlane_obj2, k4, i4);
					}

				}

				if (b[i4].g.itypegeom == CAD_STL) {


					integer k_loc=-1;
					if (b[i4].g.in_CAD_STL_check(p,k_loc,i4)) {
					   k4=i4;
					}

				}

				if (b[i4].g.itypegeom == CYLINDER) {
					// Cylinder
					switch (b[i4].g.iPlane) {
					case XY:
						if (fabs(b[i4].g.R_in_cyl) < 1.0e-36) {
							if ((p.z > b[i4].g.zC) && (p.z < b[i4].g.zC + b[i4].g.Hcyl)) {
								if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.yC - p.y)*(b[i4].g.yC - p.y)) < b[i4].g.R_out_cyl) {
									k4 = i4;
								}
							}
						}
						else {
							if ((p.z > b[i4].g.zC) && (p.z < b[i4].g.zC + b[i4].g.Hcyl)) {
								if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.yC - p.y)*(b[i4].g.yC - p.y)) < b[i4].g.R_out_cyl) {
									if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.yC - p.y)*(b[i4].g.yC - p.y)) > b[i4].g.R_in_cyl) {
										k4 = i4;
									}
								}
							}
						}
						break;
					case XZ:
						if (fabs(b[i4].g.R_in_cyl) < 1.0e-36) {
							if ((p.y > b[i4].g.yC) && (p.y < b[i4].g.yC + b[i4].g.Hcyl)) {
								if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) < b[i4].g.R_out_cyl) {
									k4 = i4;
								}
							}
						}
						else {
							if ((p.y > b[i4].g.yC) && (p.y < b[i4].g.yC + b[i4].g.Hcyl)) {
								if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) < b[i4].g.R_out_cyl) {
									if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) > b[i4].g.R_in_cyl) {
										k4 = i4;
									}
								}
							}
						}
						break;
					case YZ:
						if (fabs(b[i4].g.R_in_cyl) < 1.0e-36) {
							if ((p.x > b[i4].g.xC) && (p.x < b[i4].g.xC + b[i4].g.Hcyl)) {
								if (sqrt((b[i4].g.yC - p.y)*(b[i4].g.yC - p.y) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) < b[i4].g.R_out_cyl) {
									k4 = i4;
								}
							}
						}
						else {
							if ((p.x > b[i4].g.xC) && (p.x < b[i4].g.xC + b[i4].g.Hcyl)) {
								if (sqrt((b[i4].g.yC - p.y)*(b[i4].g.yC - p.y) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) < b[i4].g.R_out_cyl) {
									if (sqrt((b[i4].g.yC - p.y)*(b[i4].g.yC - p.y) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) > b[i4].g.R_in_cyl) {
										k4 = i4;
									}
								}
							}
						}
						break;
					}
				}
			}
		}

	}
	k = k1;
	if (k2 > k) k = k2;
	if (k3 > k) k = k3;
	if (k4 > k) k = k4;

#else
*/
   // 27_12_2017.
   // Блоки перечисляются согласно приоритетам.
   // Сначала идут блоки с низким приоритетом начиная с нуля.
   // Блок  большим номером перезаписывает свойства блока с меньшим номером.
   // Поэтому если сканировать блоки с конца (начиная с самомого большого номера lb и далее
   // в сторону уменьшения номеров, то первый найденный блок как раз и будет нужным и можно
   // досрочно прервать цикл for с помощью break и сэкономить существенную долю времени.

   // Быстрая проверка принадлежности для параллелипипедов
   // Она использует быстрый двоичный поиск.
   int kprism = myisblock_id_PRISM_only(lb,b,p); 

   k = lb; ib = lb;

   bool bfound_out = false;
   //integer *b_non_prism = nullptr;
   //integer lb_non_prism = 0;


// Внимание только однопоточный код, иначе неодназначность принадлежности блоку.
   for (int i_74 = QSBid.lb_non_prism - 1; i_74 >= 0; i_74--) {
	   if (!bfound_out) {
		   i = QSBid.b_non_prism[i_74]; // Номер непризматического объекта.
		   if (b[i].g.itypegeom == POLYGON) {

			   if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {

				   bool bfound = false;
				   int k_loc23 = -1;
				   // определяет принадлежность точки полигону.
				   bfound = in_polygon(p, b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2, k_loc23, i);
				   if (bfound) {
					   // Нашли и сразу завершили проверку в случае успеха.
					   {
						   	   k = i;
							   bfound_out = true;
							   continue;
					   }
				   }
			   }

		   }
		   else if (b[i].g.itypegeom == CAD_STL) {
			   // 14.11.2020

			   int k_loc23 = -1;
			   if (b[i].g.in_CAD_STL_check(p, k_loc23, i)) {

				   k = i;
				   bfound_out = true;
				   continue;
			   }
		   }
		   else if (b[i].g.itypegeom == CYLINDER) {
			   // Cylinder
			   switch (b[i].g.iPlane) {
			   case XY_PLANE:
				   if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.xC - p.x) * (b[i].g.xC - p.x) + (b[i].g.yC - p.y) * (b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {

						   
						   if ((fabs(b[i].g.R_in_cyl) < 1.0e-36) || (((sqrt((b[i].g.xC - p.x) * (b[i].g.xC - p.x) + (b[i].g.yC - p.y) * (b[i].g.yC - p.y)) > b[i].g.R_in_cyl)))) {

							   k = i;
							   bfound_out = true;
							   // Нашли и сразу завершили проверку в случае успеха.
							   continue;

						    }
					   }
				   }
				   break;
			   case XZ_PLANE:
				   if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.xC - p.x) * (b[i].g.xC - p.x) + (b[i].g.zC - p.z) * (b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {


						  if ((fabs(b[i].g.R_in_cyl) < 1.0e-36)||(((sqrt((b[i].g.xC - p.x) * (b[i].g.xC - p.x) + (b[i].g.zC - p.z) * (b[i].g.zC - p.z)) > b[i].g.R_in_cyl)))) {

						    	k = i;
							    bfound_out = true;
								// Нашли и сразу завершили проверку в случае успеха.
								continue;								
						  }						   
					   }
				   }
				   break;
			   case YZ_PLANE:
				   if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.yC - p.y) * (b[i].g.yC - p.y) + (b[i].g.zC - p.z) * (b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {

						  
							if ((fabs(b[i].g.R_in_cyl) < 1.0e-36)||(((sqrt((b[i].g.yC - p.y) * (b[i].g.yC - p.y) + (b[i].g.zC - p.z) * (b[i].g.zC - p.z)) > b[i].g.R_in_cyl)))) {

								  
								   k = i;
								   bfound_out = true;
								   // Нашли и сразу завершили проверку в случае успеха.
								   continue;	
							}						  
					   }
				   }
				   break;
			   }
		   }
	   }
   }
   
   /*
	// цикл по всем блокам
	for (i = lb - 1; i >= 0; i--) {
		if (b[i].g.itypegeom == PRISM) {
		    // Мы уже выполнили ранее быстрое вычисление принадлежности точки призматическому объекту.

			// Prism
			// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
			//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
			// приведет к неверному типу блока у данной ячейки.
			if ((p.x >= b[i].g.xS) && (p.x <= b[i].g.xE) && (p.y >= b[i].g.yS) && (p.y <= b[i].g.yE) && (p.z >= b[i].g.zS) && (p.z <= b[i].g.zE)) {
				k = i;
				// Нашли и сразу завершили проверку в случае успеха.
				goto OUTOF_IN_MODEL_TEMP1;
			}
		}
		 if (b[i].g.itypegeom == CAD_STL) {
			   // 14.11.2020

			   integer k_loc23 = -1;
			   if (b[i].g.in_CAD_STL_check(p, k_loc23, i)) {
				   k = i;
				   // Нашли и сразу завершили проверку в случае успеха.
				   goto OUTOF_IN_MODEL_TEMP1;
			   }
		   }
		if (b[i].g.itypegeom == POLYGON) {

			if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {

				bool bfound = false;
				// определяет принадлежность точки полигону.
				bfound = in_polygon(p, b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2, k, i);
				if (bfound) {
					// Нашли и сразу завершили проверку в случае успеха.
					goto OUTOF_IN_MODEL_TEMP1;
				}
			}

		}
		if (b[i].g.itypegeom == CYLINDER) {
			// Cylinder
			switch (b[i].g.iPlane) {
			case XY:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							k = i;
							// Нашли и сразу завершили проверку в случае успеха.
							goto OUTOF_IN_MODEL_TEMP1;
						}
					}
				}
				else {
					if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) > b[i].g.R_in_cyl) {
								k = i;
								// Нашли и сразу завершили проверку в случае успеха.
								goto OUTOF_IN_MODEL_TEMP1;
						    }
						}
					}
				}
				break;
		        case XZ:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Нашли и сразу завершили проверку в случае успеха.
							goto OUTOF_IN_MODEL_TEMP1;
						}
					}
				}
				else {
					if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Нашли и сразу завершили проверку в случае успеха.
								goto OUTOF_IN_MODEL_TEMP1;
							}
						}
					}
				}
				break;
			case YZ:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Нашли и сразу завершили проверку в случае успеха.
							goto OUTOF_IN_MODEL_TEMP1;
						}
					}
				}
				else {
					if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
					    if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Нашли и сразу завершили проверку в случае успеха.
								goto OUTOF_IN_MODEL_TEMP1;
							}
					    }
					}
				}
				break;
			}
		}
	}
	*/
	

//OUTOF_IN_MODEL_TEMP1:

//#endif

   if (kprism != lb) {
	   // Уже найден призматический блок с большим значением приоритета.
	   if (bfound_out) {
		   // Найдена принадлежность и цилиндрам и призмам.
		   // Выбираем сильнейшего с наивысшим приоритетом.
		   if (k < kprism) k = kprism;
	   }
	   else {
		   // Принадлежность к цилиндрам и полигонам не найдена.
		   //  Зато найдена принадлежность к призмам.
		   k = kprism;
	   }
	   //printf("kprism=%lld k=%lld\n", kprism,k);
	   //system("PAUSE");
   }

	ib = k;

	if (ib == lb) {
		printf("FATALL ERROR!!! ib==lb\n");
		system("PAUSE");
	}

	return ib;
} // myisblock_id 

// Интерфейс старого вызова.
integer myisblock_id(int lb, BLOCK* &b, doublereal x11, doublereal y11, doublereal z11) {
	TOCHKA p;
	p.x = x11;
	p.y = y11;
	p.z = z11;
	return (myisblock_id(lb, b, p));
} // Интерфейс старого вызова.


// проверяет принадлежит ли контрольный объём
// тепловой модели.
// Возвращает параметр ib равный номеру блока
// которому принадлежит контрольный объём.
bool in_model_temp(TOCHKA p, int &ib, BLOCK* b, int lb) {
	
	bool ret = true;// по умолчанию принадлежит модели
	ib = myisblock_id(lb, b, p);
	
	if (ib == lb) {
		ret = false;
	}
	else {
		if (b[ib].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) ret = false;
	}
	return ret;

} // in_model_temp

// проверяет принадлежит ли контрольный объём
// тепловой модели.
// Возвращает параметр ib равный номеру блока
// которому принадлежит контрольный объём.
bool in_model_temp_light(int ib, BLOCK* b, int lb) {

	//ib - передаётся внутрь функции

	bool ret = true; // по умолчанию принадлежит модели

	if (ib == lb) {
		ret = false;
	}
	else {
		if (b[ib].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) ret = false;
	}
	return ret;

} // in_model_temp

// проверяет принадлежит ли контрольный объём
// гидродинамической модели.
// Возвращает параметр ib равный номеру блока
// которому принадлежит контрольный объём.
bool in_model_flow(TOCHKA p, int &ib, BLOCK* b, int lb) {
	
	bool ret=true;// по умолчанию принадлежит модели
	ib = myisblock_id(lb, b, p);
	   	
	if (ib == lb) {
		ret = false;
	}
	else {
		if ((b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[ib].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) ret = false;
	}
	return ret;

} // in_model_flow

// проверяет принадлежит ли контрольный объём
// гидродинамической модели.
// Возвращает параметр ib равный номеру блока
// которому принадлежит контрольный объём.
bool in_model_flow_light(int ib, BLOCK* b, int lb) {

	//ib - передаётся внутрь функции

	bool ret = true;// по умолчанию принадлежит модели
	

	if (ib == lb) {
		ret = false;
	}
	else {
		if ((b[ib].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[ib].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) ret = false;
	}
	return ret;

} // in_model_flow

// проверяет принадлежит ли контрольный объём
// тепловой модели.
// Возвращает параметр ib равный номеру блока
// которому принадлежит контрольный объём.
bool in_model_temp_stab(TOCHKA p, int &ib, BLOCK* b, int lb) {

	int i = 0;
	int k = 0;
	bool ret = true;// по умолчанию принадлежит модели
	

	/*
	TOCHKA p_test;
	p_test.x = 2.100000e-03;
	p_test.y = 1.889375e-02;
	p_test.z = 2.330625e-02;
	 

	integer ib_loc;
	bool ret_loc;
	ret_loc = in_model_temp1(p_test, ib_loc, b, lb);

	if ((ib != ib_loc) || (ret_loc != ret)) {
		printf("p.x=%e p.y=%e p.z=%e\n", p.x, p.y, p.z);
		printf("lb=%lld ib=%lld ib_loc=%lld\n", lb, ib, ib_loc);
		system("PAUSE");
	}
	*/
	

	// параллельный код не является работоспособным 10.05.2019.
	// true - параллельно. false - однопоточно. 
	const bool bYES_PARALLELISM_ON = false;

#ifdef _OPENMP

	//integer lb4 = lb / 4; // Это очень плохое деление пригодное только в случае,
	// если все блоки равноправны: содержат примерно одинаковое число ячеек.

	if (bYES_PARALLELISM_ON) {
		omp_set_num_threads(4);
		int k1 = 0, k2 = 0, k3 = 0, k4 = 0;
#pragma omp parallel shared(k1,k2,k3,k4) num_threads(4)
		{
#ifdef _OPENMP 
			int tid = omp_get_thread_num();
#else
			int tid = 0;
#endif
			if (tid == 0) {
				// цикл по всем блокам
				//for (integer i1 = 0; i1<lb4; i1++) {
				for (int i1 = 0; i1 < lb; i1 += 4) {
					if (b[i1].g.itypegeom == PRISM) {
						// Prism
						// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
						//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
						// приведет к неверному типу блока у данной ячейки.
						if ((p.x >= b[i1].g.xS) && (p.x <= b[i1].g.xE) && (p.y >= b[i1].g.yS) && (p.y <= b[i1].g.yE) && (p.z >= b[i1].g.zS) && (p.z <= b[i1].g.zE)) {
							k1 = i1;
						}
					}
					if (b[i1].g.itypegeom == POLYGON) {

						if ((p.x > b[i1].g.xS) && (p.x < b[i1].g.xE) && (p.y > b[i1].g.yS) && (p.y < b[i1].g.yE) && (p.z > b[i1].g.zS) && (p.z < b[i1].g.zE)) {
							// определяет принадлежность точки полигону.
							in_polygon(p, b[i1].g.nsizei, b[i1].g.xi, b[i1].g.yi, b[i1].g.zi, b[i1].g.hi, b[i1].g.iPlane_obj2, k1, i1);
						}

					}
					if (b[i1].g.itypegeom == CAD_STL) {
						// 14.11.2020

						int k_loc23 = -1;
						if (b[i1].g.in_CAD_STL_check(p, k_loc23, i1)) {
							k1 = i1;
						}
					}

					if (b[i1].g.itypegeom == CYLINDER) {
						// Cylinder
						switch (b[i1].g.iPlane) {
						case XY_PLANE:
							if (fabs(b[i1].g.R_in_cyl) < 1.0e-36) {
								if ((p.z > b[i1].g.zC) && (p.z < b[i1].g.zC + b[i1].g.Hcyl)) {
									if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.yC - p.y)*(b[i1].g.yC - p.y)) < b[i1].g.R_out_cyl) {
										k1 = i1;
									}
								}
							}
							else {
								if ((p.z > b[i1].g.zC) && (p.z < b[i1].g.zC + b[i1].g.Hcyl)) {
									if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.yC - p.y)*(b[i1].g.yC - p.y)) < b[i1].g.R_out_cyl) {
										if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.yC - p.y)*(b[i1].g.yC - p.y)) > b[i1].g.R_in_cyl) {
											k1 = i1;
										}
									}
								}
							}
							break;
						case XZ_PLANE:
							if (fabs(b[i1].g.R_in_cyl) < 1.0e-36) {
								if ((p.y > b[i1].g.yC) && (p.y < b[i1].g.yC + b[i1].g.Hcyl)) {
									if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) < b[i1].g.R_out_cyl) {
										k1 = i1;
									}
								}
							}
							else {
								if ((p.y > b[i1].g.yC) && (p.y < b[i1].g.yC + b[i1].g.Hcyl)) {
									if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) < b[i1].g.R_out_cyl) {
										if (sqrt((b[i1].g.xC - p.x)*(b[i1].g.xC - p.x) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) > b[i1].g.R_in_cyl) {
											k1 = i1;
										}
									}
								}
							}
							break;
						case YZ_PLANE:
							if (fabs(b[i1].g.R_in_cyl) < 1.0e-36) {
								if ((p.x > b[i1].g.xC) && (p.x < b[i1].g.xC + b[i1].g.Hcyl)) {
									if (sqrt((b[i1].g.yC - p.y)*(b[i1].g.yC - p.y) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) < b[i1].g.R_out_cyl) {
										k1 = i1;
									}
								}
							}
							else {
								if ((p.x > b[i1].g.xC) && (p.x < b[i1].g.xC + b[i1].g.Hcyl)) {
									if (sqrt((b[i1].g.yC - p.y)*(b[i1].g.yC - p.y) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) < b[i1].g.R_out_cyl) {
										if (sqrt((b[i1].g.yC - p.y)*(b[i1].g.yC - p.y) + (b[i1].g.zC - p.z)*(b[i1].g.zC - p.z)) > b[i1].g.R_in_cyl) {
											k1 = i1;
										}
									}
								}
							}
							break;
						}
					}
				}
			}
			if (tid == 1) {
				// цикл по всем блокам
				//for (integer i2 = lb4; i2<2 * lb4; i2++) {
				for (int i2 = 1; i2 < lb; i2 += 4) {
					if (b[i2].g.itypegeom == PRISM) {
						// Prism
						// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
						//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
						// приведет к неверному типу блока у данной ячейки.
						if ((p.x >= b[i2].g.xS) && (p.x <= b[i2].g.xE) && (p.y >= b[i2].g.yS) && (p.y <= b[i2].g.yE) && (p.z >= b[i2].g.zS) && (p.z <= b[i2].g.zE)) {
							k2 = i2;
						}
					}
					if (b[i2].g.itypegeom == POLYGON) {

						if ((p.x > b[i2].g.xS) && (p.x < b[i2].g.xE) && (p.y > b[i2].g.yS) && (p.y < b[i2].g.yE) && (p.z > b[i2].g.zS) && (p.z < b[i2].g.zE)) {
							// определяет принадлежность точки полигону.
							in_polygon(p, b[i2].g.nsizei, b[i2].g.xi, b[i2].g.yi, b[i2].g.zi, b[i2].g.hi, b[i2].g.iPlane_obj2, k2, i2);
						}

					}

					if (b[i2].g.itypegeom == CAD_STL) {
						// 14.11.2020

						int k_loc23 = -1;
						if (b[i2].g.in_CAD_STL_check(p, k_loc23, i2)) {
							k2 = i2;
						}
					}

					if (b[i2].g.itypegeom == CYLINDER) {
						// Cylinder
						switch (b[i2].g.iPlane) {
						case XY_PLANE:
							if (fabs(b[i2].g.R_in_cyl) < 1.0e-36) {
								if ((p.z > b[i2].g.zC) && (p.z < b[i2].g.zC + b[i2].g.Hcyl)) {
									if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.yC - p.y)*(b[i2].g.yC - p.y)) < b[i2].g.R_out_cyl) {
										k2 = i2;
									}
								}
							}
							else {
								if ((p.z > b[i2].g.zC) && (p.z < b[i2].g.zC + b[i2].g.Hcyl)) {
									if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.yC - p.y)*(b[i2].g.yC - p.y)) < b[i2].g.R_out_cyl) {
										if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.yC - p.y)*(b[i2].g.yC - p.y)) > b[i2].g.R_in_cyl) {
											k2 = i2;
										}
									}
								}
							}
							break;
						case XZ_PLANE:
							if (fabs(b[i2].g.R_in_cyl) < 1.0e-36) {
								if ((p.y > b[i2].g.yC) && (p.y < b[i2].g.yC + b[i2].g.Hcyl)) {
									if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) < b[i2].g.R_out_cyl) {
										k2 = i2;
									}
								}
							}
							else {
								if ((p.y > b[i2].g.yC) && (p.y < b[i2].g.yC + b[i2].g.Hcyl)) {
									if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) < b[i2].g.R_out_cyl) {
										if (sqrt((b[i2].g.xC - p.x)*(b[i2].g.xC - p.x) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) > b[i2].g.R_in_cyl) {
											k2 = i2;
										}
									}
								}
							}
							break;
						case YZ_PLANE:
							if (fabs(b[i2].g.R_in_cyl) < 1.0e-36) {
								if ((p.x > b[i2].g.xC) && (p.x < b[i2].g.xC + b[i2].g.Hcyl)) {
									if (sqrt((b[i2].g.yC - p.y)*(b[i2].g.yC - p.y) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) < b[i2].g.R_out_cyl) {
										k2 = i2;
									}
								}
							}
							else {
								if ((p.x > b[i2].g.xC) && (p.x < b[i2].g.xC + b[i2].g.Hcyl)) {
									if (sqrt((b[i2].g.yC - p.y)*(b[i2].g.yC - p.y) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) < b[i2].g.R_out_cyl) {
										if (sqrt((b[i2].g.yC - p.y)*(b[i2].g.yC - p.y) + (b[i2].g.zC - p.z)*(b[i2].g.zC - p.z)) > b[i2].g.R_in_cyl) {
											k2 = i2;
										}
									}
								}
							}
							break;
						}
					}
				}
			}
			if (tid == 2) {
				// цикл по всем блокам
				//for (integer i3 = 2 * lb4; i3<3 * lb4; i3++) {
				for (int i3 = 2; i3 < lb; i3 += 4) {
					if (b[i3].g.itypegeom == PRISM) {
						// Prism
						// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
						//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
						// приведет к неверному типу блока у данной ячейки.
						if ((p.x >= b[i3].g.xS) && (p.x <= b[i3].g.xE) && (p.y >= b[i3].g.yS) && (p.y <= b[i3].g.yE) && (p.z >= b[i3].g.zS) && (p.z <= b[i3].g.zE)) {
							k3 = i3;
						}
					}
					if (b[i3].g.itypegeom == POLYGON) {

						if ((p.x > b[i3].g.xS) && (p.x < b[i3].g.xE) && (p.y > b[i3].g.yS) && (p.y < b[i3].g.yE) && (p.z > b[i3].g.zS) && (p.z < b[i3].g.zE)) {
							// определяет принадлежность точки полигону.
							in_polygon(p, b[i3].g.nsizei, b[i3].g.xi, b[i3].g.yi, b[i3].g.zi, b[i3].g.hi, b[i3].g.iPlane_obj2, k3, i3);
						}

					}

					if (b[i3].g.itypegeom == CAD_STL) {
						// 14.11.2020

						int k_loc23 = -1;
						if (b[i3].g.in_CAD_STL_check(p, k_loc23, i3)) {
							k3 = i3;
						}
					}

					if (b[i3].g.itypegeom == CYLINDER) {
						// Cylinder
						switch (b[i3].g.iPlane) {
						case XY_PLANE:
							if (fabs(b[i3].g.R_in_cyl) < 1.0e-36) {
								if ((p.z > b[i3].g.zC) && (p.z < b[i3].g.zC + b[i3].g.Hcyl)) {
									if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.yC - p.y)*(b[i3].g.yC - p.y)) < b[i3].g.R_out_cyl) {
										k3 = i3;
									}
								}
							}
							else {
								if ((p.z > b[i3].g.zC) && (p.z < b[i3].g.zC + b[i3].g.Hcyl)) {
									if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.yC - p.y)*(b[i3].g.yC - p.y)) < b[i3].g.R_out_cyl) {
										if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.yC - p.y)*(b[i3].g.yC - p.y)) > b[i3].g.R_in_cyl) {
											k3 = i3;
										}
									}
								}
							}
							break;
						case XZ_PLANE:
							if (fabs(b[i3].g.R_in_cyl) < 1.0e-36) {
								if ((p.y > b[i3].g.yC) && (p.y < b[i3].g.yC + b[i3].g.Hcyl)) {
									if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) < b[i3].g.R_out_cyl) {
										k3 = i3;
									}
								}
							}
							else {
								if ((p.y > b[i3].g.yC) && (p.y < b[i3].g.yC + b[i3].g.Hcyl)) {
									if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) < b[i3].g.R_out_cyl) {
										if (sqrt((b[i3].g.xC - p.x)*(b[i3].g.xC - p.x) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) > b[i3].g.R_in_cyl) {
											k3 = i3;
										}
									}
								}
							}
							break;
						case YZ_PLANE:
							if (fabs(b[i3].g.R_in_cyl) < 1.0e-36) {
								if ((p.x > b[i3].g.xC) && (p.x < b[i3].g.xC + b[i3].g.Hcyl)) {
									if (sqrt((b[i3].g.yC - p.y)*(b[i3].g.yC - p.y) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) < b[i3].g.R_out_cyl) {
										k3 = i3;
									}
								}
							}
							else {
								if ((p.x > b[i3].g.xC) && (p.x < b[i3].g.xC + b[i3].g.Hcyl)) {
									if (sqrt((b[i3].g.yC - p.y)*(b[i3].g.yC - p.y) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) < b[i3].g.R_out_cyl) {
										if (sqrt((b[i3].g.yC - p.y)*(b[i3].g.yC - p.y) + (b[i3].g.zC - p.z)*(b[i3].g.zC - p.z)) > b[i3].g.R_in_cyl) {
											k3 = i3;
										}
									}
								}
							}
							break;
						}
					}
				}
			}
			if (tid == 3) {
				// цикл по всем блокам
				//for (integer i4 = 3 * lb4; i4<lb; i4++) {
				for (int i4 = 3; i4 < lb; i4 += 4) {
					if (b[i4].g.itypegeom == PRISM) {
						// Prism
						// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
						//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
						// приведет к неверному типу блока у данной ячейки.
						if ((p.x >= b[i4].g.xS) && (p.x <= b[i4].g.xE) && (p.y >= b[i4].g.yS) && (p.y <= b[i4].g.yE) && (p.z >= b[i4].g.zS) && (p.z <= b[i4].g.zE)) {
							k4 = i4;
						}
					}
					if (b[i4].g.itypegeom == POLYGON) {

						if ((p.x > b[i4].g.xS) && (p.x < b[i4].g.xE) && (p.y > b[i4].g.yS) && (p.y < b[i4].g.yE) && (p.z > b[i4].g.zS) && (p.z < b[i4].g.zE)) {
							// определяет принадлежность точки полигону.
							in_polygon(p, b[i4].g.nsizei, b[i4].g.xi, b[i4].g.yi, b[i4].g.zi, b[i4].g.hi, b[i4].g.iPlane_obj2, k4, i4);
						}

					}

					if (b[i4].g.itypegeom == CAD_STL) {
						// 14.11.2020

						int k_loc23 = -1;
						if (b[i4].g.in_CAD_STL_check(p, k_loc23, i4)) {
							k4 = i4;
						}
					}

					if (b[i4].g.itypegeom == CYLINDER) {
						// Cylinder
						switch (b[i4].g.iPlane) {
						case XY_PLANE:
							if (fabs(b[i4].g.R_in_cyl) < 1.0e-36) {
								if ((p.z > b[i4].g.zC) && (p.z < b[i4].g.zC + b[i4].g.Hcyl)) {
									if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.yC - p.y)*(b[i4].g.yC - p.y)) < b[i4].g.R_out_cyl) {
										k4 = i4;
									}
								}
							}
							else {
								if ((p.z > b[i4].g.zC) && (p.z < b[i4].g.zC + b[i4].g.Hcyl)) {
									if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.yC - p.y)*(b[i4].g.yC - p.y)) < b[i4].g.R_out_cyl) {
										if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.yC - p.y)*(b[i4].g.yC - p.y)) > b[i4].g.R_in_cyl) {
											k4 = i4;
										}
									}
								}
							}
							break;
						case XZ_PLANE:
							if (fabs(b[i4].g.R_in_cyl) < 1.0e-36) {
								if ((p.y > b[i4].g.yC) && (p.y < b[i4].g.yC + b[i4].g.Hcyl)) {
									if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) < b[i4].g.R_out_cyl) {
										k4 = i4;
									}
								}
							}
							else {
								if ((p.y > b[i4].g.yC) && (p.y < b[i4].g.yC + b[i4].g.Hcyl)) {
									if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) < b[i4].g.R_out_cyl) {
										if (sqrt((b[i4].g.xC - p.x)*(b[i4].g.xC - p.x) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) > b[i4].g.R_in_cyl) {
											k4 = i4;
										}
									}
								}
							}
							break;
						case YZ_PLANE:
							if (fabs(b[i4].g.R_in_cyl) < 1.0e-36) {
								if ((p.x > b[i4].g.xC) && (p.x < b[i4].g.xC + b[i4].g.Hcyl)) {
									if (sqrt((b[i4].g.yC - p.y)*(b[i4].g.yC - p.y) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) < b[i4].g.R_out_cyl) {
										k4 = i4;
									}
								}
							}
							else {
								if ((p.x > b[i4].g.xC) && (p.x < b[i4].g.xC + b[i4].g.Hcyl)) {
									if (sqrt((b[i4].g.yC - p.y)*(b[i4].g.yC - p.y) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) < b[i4].g.R_out_cyl) {
										if (sqrt((b[i4].g.yC - p.y)*(b[i4].g.yC - p.y) + (b[i4].g.zC - p.z)*(b[i4].g.zC - p.z)) > b[i4].g.R_in_cyl) {
											k4 = i4;
										}
									}
								}
							}
							break;
						}
					}
				}
			}

		}
#ifdef _OPENMP
		omp_set_num_threads(1);
#endif
		k = k1;
		if (k2 > k) k = k2;
		if (k3 > k) k = k3;
		if (k4 > k) k = k4;

}
	else {
	// Однопоточная заглушка
	// цикл по всем блокам
	for (i = lb - 1; i >= 0; i--) {
		if (b[i].g.itypegeom == PRISM) {
			// Prism
			// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
			//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
			// приведет к неверному типу блока у данной ячейки.
			if ((p.x >= b[i].g.xS) && (p.x <= b[i].g.xE) && (p.y >= b[i].g.yS) && (p.y <= b[i].g.yE) && (p.z >= b[i].g.zS) && (p.z <= b[i].g.zE)) {
				k = i;
				// Нашли и сразу завершили проверку в случае успеха.
				goto OUTOF_IN_MODEL_TEMP1;
			}
		}
		if (b[i].g.itypegeom == POLYGON) {

			if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {

				bool bfound = false;
				// определяет принадлежность точки полигону.
				bfound = in_polygon(p, b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2, k, i);
				if (bfound) {
					// Нашли и сразу завершили проверку в случае успеха.
					goto OUTOF_IN_MODEL_TEMP1;
				}
			}

		}
		if (b[i].g.itypegeom == CAD_STL) {
			// 14.11.2020

			int k_loc23 = -1;
			if (b[i].g.in_CAD_STL_check(p, k_loc23, i)) {
				k = i;
				// Нашли и сразу завершили проверку в случае успеха.
				goto OUTOF_IN_MODEL_TEMP1;
			}
		}
		if (b[i].g.itypegeom == CYLINDER) {
			// Cylinder
			switch (b[i].g.iPlane) {
			case XY_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							k = i;
							// Нашли и сразу завершили проверку в случае успеха.
							goto OUTOF_IN_MODEL_TEMP1;
						}
					}
				}
				else {
					if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) > b[i].g.R_in_cyl) {
								k = i;
								// Нашли и сразу завершили проверку в случае успеха.
								goto OUTOF_IN_MODEL_TEMP1;
							}
						}
					}
				}
				break;
			case XZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Нашли и сразу завершили проверку в случае успеха.
							goto OUTOF_IN_MODEL_TEMP1;
						}
					}
				}
				else {
					if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Нашли и сразу завершили проверку в случае успеха.
								goto OUTOF_IN_MODEL_TEMP1;
							}
						}
					}
				}
				break;
			case YZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Нашли и сразу завершили проверку в случае успеха.
							goto OUTOF_IN_MODEL_TEMP1;
						}
					}
				}
				else {
					if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Нашли и сразу завершили проверку в случае успеха.
								goto OUTOF_IN_MODEL_TEMP1;
							}
						}
					}
				}
				break;
			}
		}
	}
	}

#else

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
			// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
			//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
			// приведет к неверному типу блока у данной ячейки.
			if ((p.x >= b[i].g.xS) && (p.x <= b[i].g.xE) && (p.y >= b[i].g.yS) && (p.y <= b[i].g.yE) && (p.z >= b[i].g.zS) && (p.z <= b[i].g.zE)) {
				k = i;
				// Нашли и сразу завершили проверку в случае успеха.
				goto OUTOF_IN_MODEL_TEMP1;
			}
		}
		if (b[i].g.itypegeom == POLYGON) {

			if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {

				bool bfound = false;
				// определяет принадлежность точки полигону.
				bfound = in_polygon(p, b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2, k, i);
				if (bfound) {
					// Нашли и сразу завершили проверку в случае успеха.
					goto OUTOF_IN_MODEL_TEMP1;
				}
			}

		}
		if (b[i].g.itypegeom == CAD_STL) {
			// 14.11.2020

			int k_loc23 = -1;
			if (b[i].g.in_CAD_STL_check(p, k_loc23, i)) {
				k = i;
				// Нашли и сразу завершили проверку в случае успеха.
				goto OUTOF_IN_MODEL_TEMP1;
			}
		}
		if (b[i].g.itypegeom == CYLINDER) {
			// Cylinder
			switch (b[i].g.iPlane) {
			case XY_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							k = i;
							// Нашли и сразу завершили проверку в случае успеха.
							goto OUTOF_IN_MODEL_TEMP1;
						}
					}
				}
				else {
					if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) > b[i].g.R_in_cyl) {
								k = i;
								// Нашли и сразу завершили проверку в случае успеха.
								goto OUTOF_IN_MODEL_TEMP1;
							}
						}
					}
				}
				break;
			case XZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Нашли и сразу завершили проверку в случае успеха.
							goto OUTOF_IN_MODEL_TEMP1;
						}
					}
				}
				else {
					if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Нашли и сразу завершили проверку в случае успеха.
								goto OUTOF_IN_MODEL_TEMP1;
							}
						}
					}
				}
				break;
			case YZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Нашли и сразу завершили проверку в случае успеха.
							goto OUTOF_IN_MODEL_TEMP1;
						}
					}
				}
				else {
					if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Нашли и сразу завершили проверку в случае успеха.
								goto OUTOF_IN_MODEL_TEMP1;
							}
						}
					}
				}
				break;
			}
		}
	}

#endif

OUTOF_IN_MODEL_TEMP1:

	if (b[k].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) ret = false;
	ib = k;

	/*
	integer ib_loc;
	bool ret_loc;
	ret_loc=in_model_temp1(p, ib_loc, b, lb);

	if ((ib != ib_loc) || (ret_loc != ret)) {
		printf("p.x=%e p.y=%e p.z=%e\n",p.x,p.y,p.z);
		printf("lb=%lld ib=%lld ib_loc=%lld\n",lb,ib,ib_loc);
		system("PAUSE");
	}
	*/

	return ret;

} // in_model_temp_stab

// проверяет принадлежит ли контрольный объём
// гидродинамической модели.
// Возвращает параметр ib равный номеру блока
// которому принадлежит контрольный объём.
bool in_model_flow_stab(TOCHKA p, int &ib, BLOCK* &b, int lb) {
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
			// Обязательно должно быть <= или >= а не просто < или >. Если будет просто < или >,
			//то центр некоторых ячеек обязательно попадет прямо на самую границу двух блоков, что 
			// приведет к неверному типу блока у данной ячейки.
			if ((p.x >= b[i].g.xS) && (p.x <= b[i].g.xE) && (p.y >= b[i].g.yS) && (p.y <= b[i].g.yE) && (p.z >= b[i].g.zS) && (p.z <= b[i].g.zE)) {
				k = i;
				// Только нашли и сразу закончили проверку.
				goto OUTOF_IN_MODEL_FLOW;
			}
		}
		if (b[i].g.itypegeom == POLYGON) {

			if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {

				bool found = false;
				// определяет принадлежность точки полигону.
				found = in_polygon(p, b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2, k, i);
				if (found) {
					// Только нашли и сразу закончили проверку.
					goto OUTOF_IN_MODEL_FLOW;
				}

			}

		}

		if (b[i].g.itypegeom == CAD_STL) {
			// 14.11.2020

			int k_loc23 = -1;
			if (b[i].g.in_CAD_STL_check(p, k_loc23, i)) {
				k = i;
				// Нашли и сразу завершили проверку в случае успеха.
				goto OUTOF_IN_MODEL_FLOW;
			}
		}

		if (b[i].g.itypegeom == CYLINDER) {
			// Cylinder
			switch (b[i].g.iPlane) {
			case XY_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							k = i;
							// Только нашли и сразу закончили проверку.
							goto OUTOF_IN_MODEL_FLOW;
						}
					}
				}
				else {
					if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) > b[i].g.R_in_cyl) {
								k = i;
								// Только нашли и сразу закончили проверку.
								goto OUTOF_IN_MODEL_FLOW;
							}
						}
					}
				}
				break;
			case XZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Только нашли и сразу закончили проверку.
							goto OUTOF_IN_MODEL_FLOW;
						}
					}
				}
				else {
					if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Только нашли и сразу закончили проверку.
								goto OUTOF_IN_MODEL_FLOW;
							}
						}
					}
				}
				break;
			case YZ_PLANE:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-36) {
					if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							k = i;
							// Только нашли и сразу закончили проверку.
							goto OUTOF_IN_MODEL_FLOW;
						}
					}
				}
				else {
					if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
						if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
							if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) > b[i].g.R_in_cyl) {
								k = i;
								// Только нашли и сразу закончили проверку.
								goto OUTOF_IN_MODEL_FLOW;
							}
						}
					}
				}
				break;
			}
		}
	}

OUTOF_IN_MODEL_FLOW:

	if ((b[k].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[k].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) ret = false;
	ib = k;

	return ret;

} // in_model_flow_stab


#include "constr_struct_alice.cpp" // Содержит функционал вызываемый в loadTemper_and_Flow при АЛИС сетке.

 

// только для enumerate_volume_improved
// См. заголовочный файл adaptive_local_refinement_mesh.cpp
//typedef struct TBlock_indexes {
	//integer iL, iR, jL, jR, kL, kR;
//} Block_indexes;

  // глобальная нумерация контрольных объёмов
  // для задач теплопроводности.
// 25.03.2017 improved версия более быстрая по скорости выполнения.
// данный метод ускорения быстродействия работает только для прямоугольных призм.
// 2.04.2017 распараллеленная версия.
void enumerate_volume_improved(int* &evt, int &maxelm, integer iflag,
	doublereal* &xpos, doublereal* &ypos, doublereal* &zpos, int* &whot_is_block,
	int inx, int iny, int inz, BLOCK* &b, int lb,
	WALL* &w, int &lw, SOURCE* &s, int &ls) {

#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_max_threads();
	omp_set_num_threads(8); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
#else
	int i_my_num_core_parallelesation = 1;
#endif
	

	evt = nullptr;
	evt = new int[inx*iny*inz];
	if (evt == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for evt constr struct...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}
	whot_is_block = nullptr;
	whot_is_block = new int[inx*iny*inz];
	if (whot_is_block == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for whot_is_block constr struct...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}
	

	if (print_log_message) {
		printf("enumerate_volume_improved start.\n");
	}
	Block_indexes* block_indexes = new Block_indexes[lb];
	// оператор new не требует проверки на null.
	//if (block_indexes==nullptr) {
		//printf("error in allocation memory for block_indexes in enumerate_volume_improved.\n");
		//system("pause");
		//exit(1);
	//}


	// 07.08.2021 Изменение допусков т.к. на большой сетке
	// структурированой для гидродинамики в 8млн узлов не распознались входная и выходная границы.
	// Для того чтобы входная и выходная границы распознавались введены допуски min_zazor_* вместо
	// допусков shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls).
	doublereal min_zazor_x = 1.0e36;
	doublereal min_zazor_y = 1.0e36;
	doublereal min_zazor_z = 1.0e36;

	for (int i = 0; i < inx; ++i) {
		if (fabs(xpos[i + 1] - xpos[i]) < min_zazor_x) {

			min_zazor_x = fabs(xpos[i + 1] - xpos[i]);
		}
	}
	min_zazor_x *= 0.25;

	for (int i = 0; i < iny; ++i) {
		if (fabs(ypos[i + 1] - ypos[i]) < min_zazor_y) {

			min_zazor_y = fabs(ypos[i + 1] - ypos[i]);
		}
	}
	min_zazor_y *= 0.25;

	for (int i = 0; i < inz; ++i) {
		if (fabs(zpos[i + 1] - zpos[i]) < min_zazor_z) {

			min_zazor_z = fabs(zpos[i + 1] - zpos[i]);
		}
	}
	min_zazor_z *= 0.25;


	//int i, j, k;

	// 08.04.2018
	for (int i = 0; i < lb; ++i) {
		// инициализация, на случай если блоки не будут распознаны.
		block_indexes[i].iL = -1;
		block_indexes[i].iR = -2;
		block_indexes[i].jL = -1;
		block_indexes[i].jR = -2;
		block_indexes[i].kL = -1;
		block_indexes[i].kR = -2;
	}

	for (int i = 0; i < lb; ++i) {
		int i_1 = i;
		doublereal x4 = b[i].g.xS;
		/*for (j = 0; j <= inx; ++j) {
			if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].iL = j;
				break;
			}
		}*/
		bool bfound = false;
		if (!bfound) {
			for (int j = 0; j <= inx; ++j) {
				if ((!bfound) && (fabs(x4 - xpos[j]) < min_zazor_x /*shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)*/))
				{
					// Точное совпадение !!!
					block_indexes[i_1].iL = j;
					bfound = true;
					break;
				}
			}
		}
		doublereal tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (int j = 0; j <= inx; ++j) {
				if (fabs(x4 - xpos[j]) < tzazor) {
					tzazor = fabs(x4 - xpos[j]);
					block_indexes[i_1].iL = j;
					bfound = true;
				}
			}
		}

		if (!bfound) {
			for (int j = 0; j <= inx; ++j) {
				if ((!bfound) && (x4 < xpos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].iL = j;
					bfound = true;
					break;
				}
			}
		}
		x4 = b[i].g.xE;
		/*for (j = 0; j <= inx; ++j) {
			if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].iR = j;
				break;
			}
		}*/
		bfound = false;
		if (!bfound) {
			for (int j = inx; j >= 0; j--) {
				if ((!bfound) && (fabs(x4 - xpos[j]) < min_zazor_x /*shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].iR = j;
					bfound = true;
					break;
				}
			}
		}
		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (int j = 0; j <= inx; ++j) {
				if (fabs(x4 - xpos[j]) < tzazor) {
					tzazor = fabs(x4 - xpos[j]);
					block_indexes[i_1].iR = j;
					bfound = true;
				}
			}
		}

		if (!bfound) {
			//std::cout << "Error incomming xE\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (int j = inx; j >= 0; j--) {
				if ((!bfound) && (x4 > xpos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].iR = j;
					bfound = true;
					break;
				}
			}
		}

		x4 = b[i].g.yS;
		/*for (j = 0; j <= iny; ++j) {
			if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].jL = j;
				break;
			}
		}*/
		bfound = false;
		if (!bfound) {
			for (int j = 0; j <= iny; ++j) {
				if ((!bfound) && (fabs(x4 - ypos[j]) < min_zazor_y /*shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].jL = j;
					bfound = true;
					break;
				}
			}
		}

		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (int j = 0; j <= iny; ++j) {
				if (fabs(x4 - ypos[j]) < tzazor) {
					tzazor = fabs(x4 - ypos[j]);
					block_indexes[i_1].jL = j;
					bfound = true;
				}
			}
		}
		if (!bfound) {
			//std::cout << "Error incomming yS\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (int j = 0; j <= iny; ++j) {
				if ((!bfound) && (x4 < ypos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].jL = j;
					bfound = true;
					break;
				}
			}
		}

		x4 = b[i].g.yE;
		/*for (j = 0; j <= iny; ++j) {
			if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].jR = j;
				break;
			}
		}*/
		bfound = false;
		if (!bfound) {
			for (int j = iny; j >= 0; j--) {
				if ((!bfound) && (fabs(x4 - ypos[j]) < min_zazor_y /*shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].jR = j;
					bfound = true;
					break;
				}
			}
		}
		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (int j = 0; j <= iny; ++j) {
				if (fabs(x4 - ypos[j]) < tzazor) {
					tzazor = fabs(x4 - ypos[j]);
					block_indexes[i_1].jR = j;
					bfound = true;
				}
			}
		}
		if (!bfound) {
			//std::cout << "Error incomming yE\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (int j = iny; j >= 0; j--) {
				if ((!bfound) && (x4 > ypos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].jR = j;
					bfound = true;
					break;
				}
			}
		}

		x4 = b[i].g.zS;
		/*for (j = 0; j <= inz; ++j) {
			if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].kL = j;
				break;
			}
		}*/
		bfound = false;
		if (!bfound) {
			for (int j = 0; j <= inz; ++j) {
				if ((!bfound) && (fabs(x4 - zpos[j]) < min_zazor_z /*shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].kL = j;
					bfound = true;
					break;
				}
			}
		}
		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (int j = 0; j <= inz; ++j) {
				if (fabs(x4 - zpos[j]) < tzazor) {
					tzazor = fabs(x4 - zpos[j]);
					block_indexes[i_1].kL = j;
					bfound = true;
				}
			}
		}
		if (!bfound) {

			//std::cout << "Error incomming zS\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (int j = 0; j <= inz; ++j) {
				if ((!bfound) && (x4 < zpos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].kL = j;
					bfound = true;
					break;
				}
			}
		}

		x4 = b[i].g.zE;
		/*for (j = 0; j <= inz; ++j) {
			if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].kR = j;
				break;
			}
		}*/
		bfound = false;
		// Следуюшая секция.
		if (!bfound) {
			for (int j = inz; j >= 0; j--) {
				if ((!bfound) && (fabs(x4 - zpos[j]) < min_zazor_z /* shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].kR = j;
					bfound = true;
					break;
				}
			}
		}
		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (int j = 0; j <= inz; ++j) {
				if (fabs(x4 - zpos[j]) < tzazor) {
					tzazor = fabs(x4 - zpos[j]);
					block_indexes[i_1].kR = j;
					bfound = true;
				}
			}
		}
		if (!bfound) {

			//std::cout << "Error incomming zE\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (int j = inz; j >= 0; j--) {
				if ((!bfound) && (x4 > zpos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].kR = j;
					bfound = true;
					break;
				}
			}
		}



		if ((block_indexes[i].iL >= 0) &&
			(block_indexes[i].iR >= 0) &&
			(block_indexes[i].iL >= block_indexes[i].iR)) {
			printf("enumerate_volume_improved function\n");
			printf("violation of the order block_indexes\n");
			std::cout << b[i].name << std::endl;
			std::cout << "xS=" << b[i].g.xS << " xE=" << b[i].g.xE << std::endl;
			printf("i=%d iL=%lld iR=%lld\n", i, block_indexes[i].iL,
				block_indexes[i].iR);
			system("pause");
			for (int i92 = 0; i92 <= inx; i92++) {
				std::cout << i92 << " " << xpos[i92] << " \n";
				system("pause");
			}
		}

		if ((block_indexes[i].iL >= 0) &&
			(block_indexes[i].iR >= 0) &&
			(block_indexes[i].iL + 1 == block_indexes[i].iR)) {
			// Всего одна клетка на блок.
			//iX_one_CELL_count_statistic++;
		}

		if ((block_indexes[i].jL >= 0) &&
			(block_indexes[i].jR >= 0) &&
			(block_indexes[i].jL >= block_indexes[i].jR)) {
			printf("enumerate_volume_improved function\n");
			printf("violation of the order block_indexes\n");
			std::cout << b[i].name << std::endl;
			std::cout << "yS=" << b[i].g.yS << " yE=" << b[i].g.yE << std::endl;
			printf("i=%d jL=%lld jR=%lld\n", i, block_indexes[i].jL,
				block_indexes[i].jR);
			system("pause");
			for (int i92 = 0; i92 <= iny; i92++) {
				std::cout << i92 << " " << ypos[i92] << " \n";
				system("pause");
			}
		}


		if ((block_indexes[i].jL >= 0) &&
			(block_indexes[i].jR >= 0) &&
			(block_indexes[i].jL + 1 == block_indexes[i].jR)) {
			// Всего одна клетка на блок.
			//iY_one_CELL_count_statistic++;
		}

		if ((block_indexes[i].kL >= 0) &&
			(block_indexes[i].kR >= 0) &&
			(block_indexes[i].kL >= block_indexes[i].kR)) {
			printf("enumerate_volume_improved function\n");
			printf("violation of the order block_indexes\n");
			std::cout << "zS=" << b[i].g.zS << " zE=" << b[i].g.zE << std::endl;
			std::cout << b[i].name << std::endl;
			printf("i=%d kL=%lld kR=%lld\n", i, block_indexes[i].kL,
				block_indexes[i].kR);
			system("pause");
			for (int i92 = 0; i92 <= inz; i92++) {
				std::cout << i92 << " " << zpos[i92] << " \n";
				system("pause");
			}
		}


		if ((block_indexes[i].kL >= 0) &&
			(block_indexes[i].kR >= 0) &&
			(block_indexes[i].kL + 1 == block_indexes[i].kR)) {
			// Всего одна клетка на блок.
			//iZ_one_CELL_count_statistic++;
		}

		if (0 && (i == 116)) {
			// debug
			printf("i=%d iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
				block_indexes[i].iL,
				block_indexes[i].iR, block_indexes[i].jL,
				block_indexes[i].jR, block_indexes[i].kL,
				block_indexes[i].kR);
			printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[i].g.xS, b[i].g.xE,
				b[i].g.yS, b[i].g.yE, b[i].g.zS, b[i].g.zE);
			system("pause");
		}

		if ((block_indexes[i].iL < 0) ||
			(block_indexes[i].iR < 0) ||
			(block_indexes[i].jL < 0) ||
			(block_indexes[i].jR < 0) ||
			(block_indexes[i].kL < 0) ||
			(block_indexes[i].kR < 0)) {
			printf("enumerate_volume_improved function\n");
			printf("i=%d iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
				block_indexes[i].iL,
				block_indexes[i].iR, block_indexes[i].jL,
				block_indexes[i].jR, block_indexes[i].kL,
				block_indexes[i].kR);
			printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[i].g.xS, b[i].g.xE,
				b[i].g.yS, b[i].g.yE, b[i].g.zS, b[i].g.zE);
			printf("cabinet: xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[0].g.xS, b[0].g.xE,
				b[0].g.yS, b[0].g.yE, b[0].g.zS, b[0].g.zE);
			printf("ERROR: may be your geometry out of cabinet...\n");
			system("pause");
		}


	}

	
	// Количество проходов существенно сократилось и в итоге это приводит к существенному
	// увеличению быстродействия.
	//int m7;

	/*
	integer ib_stub = -1;
	// Мы найдем самый большой по размеру Hollow block, иначе будет просто кабинет.
	ib_stub = 0;
	doublereal vol_stub = -1.0;
	for (i = 0; i < lb; ++i) {
		if (b[i].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
			if (fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS) > vol_stub) {
				ib_stub = i;
				vol_stub = fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS);
			}
		}
	}
	*/
const integer isize_xyz=inx*iny*inz;

#pragma omp parallel for
	for (integer iP=0; iP<isize_xyz; ++iP) {
	//for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
		//evt[i1 + j1*inx + k1*inx*iny] = -1;// -1
		evt[iP]=-1; // -1
	}

	// integer iP_k1, iP_j1, iP;

	for (int m7 = 0; m7 < lb; ++m7) {

#pragma omp parallel for
		for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; ++k1)
		{
			integer iP_k1=k1*inx*iny;
			for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; ++j1)
			{
				integer iP_j1=j1*inx +iP_k1;
				 for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; ++i1)
				{
                    integer iP=i1+iP_j1;
					switch (iflag) {
						case TEMPERATURE:
						if (b[m7].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
							evt[iP] = -1;
						}
						else {
							evt[iP] = m7;
						}
						break;
						case HYDRODINAMIC:
						if ((b[m7].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m7].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
							evt[iP] = -1;
						}
						else {
							evt[iP] = m7;
						}
						break;
					}
			
				}
			}
		}
	}
	if (print_log_message) {
		printf("enumerate_volume_improved 80 procent.\n");
	}
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	int l = 1, ib;
    //  integer iP_k, iP_j;
	for (int k = 0; k < inz; ++k)
	{
		integer iP_k = k * inx * iny;
		for (int j = 0; j < iny; ++j)
		{
			integer iP_j = j * inx + iP_k;
			for (int i = 0; i < inx; ++i)
			{
				integer iP = i + iP_j;
				if (evt[iP] > -1) {
					ib = evt[iP]; // номер блока был сохранён ранее.
					// Это очень нужно для записи репорта.
					whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
					evt[iP] = l;
					l++;
				}
				else
				{   // не принадлежит расчётной области
					evt[iP] = 0;
				}
			}
		}
	}
	maxelm = l - 1;
	
	delete[] block_indexes;

	if (print_log_message) {
		printf("enumerate_volume_improved end.\n");
	}

	/*
	bool inDomain = false;
	// нумерация контрольных объёмов начинается с единицы
	// если контрольный объём не принадлежит расчётной области то ставится 0.
	for (i = 0; i<inx; ++i) for (j = 0; j<iny; ++j) for (k = 0; k<inz; ++k) {
		p.x = 0.5*(xpos[i] + xpos[i + 1]);
		p.y = 0.5*(ypos[j] + ypos[j + 1]);
		p.z = 0.5*(zpos[k] + zpos[k + 1]);
		switch (iflag) {
		case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
		case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
		}
		if (inDomain) {
			// принадлежит расчётной области
			evt[i + j*inx + k*inx*iny] = l;
			// Это очень нужно для записи репорта.
			whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
			l++;
		}
		else
		{   // не принадлежит расчётной области
			evt[i + j*inx + k*inx*iny] = 0;
		}
	}

	// maxelm - число контрольных объёмов принадлежащих расчётной области
	maxelm = l - 1;
	*/

#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif

} // enumerate_volume_improved

// глобальная нумерация контрольных объёмов
// для задач теплопроводности.
// 25.03.2017 improved версия более быстрая по скорости выполнения.
// данный метод ускорения быстродействия работает только для прямоугольных призм.
// Цилиндры и полигоны обрабатываются обычным образом.
// 2.04.2017 распараллеленная версия.
void enumerate_volume_improved_obobshenie(int* &evt, int &maxelm, integer iflag,
	doublereal* &xpos, doublereal* &ypos, doublereal* &zpos, int* &whot_is_block,
	int inx, int iny, int inz, BLOCK* &b, int lb, TOCKA_SHORT_INT* &tck_int_list,
	WALL* &w, int &lw, SOURCE* &s, int &ls) {

#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_max_threads();
	omp_set_num_threads(8); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
#else
	int i_my_num_core_parallelesation = 1;
#endif

	// 07.08.2021 Изменение допусков т.к. на большой сетке
	// структурированой для гидродинамики в 8млн узлов не распознались входная и выходная границы.
	// Для того чтобы входная и выходная границы распознавались введены допуски min_zazor_* вместо
	// допусков shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls).
	doublereal min_zazor_x = 1.0e36;
	doublereal min_zazor_y = 1.0e36;
	doublereal min_zazor_z = 1.0e36;

	for (int i = 0; i < inx; ++i) {
		if (fabs(xpos[i + 1] - xpos[i]) < min_zazor_x) {

			min_zazor_x = fabs(xpos[i + 1] - xpos[i]);
		}
	}
	min_zazor_x *= 0.25;

	for (int i = 0; i < iny; ++i) {
		if (fabs(ypos[i + 1] - ypos[i]) < min_zazor_y) {

			min_zazor_y = fabs(ypos[i + 1] - ypos[i]);
		}
	}
	min_zazor_y *= 0.25;

	for (int i = 0; i < inz; ++i) {
		if (fabs(zpos[i + 1] - zpos[i]) < min_zazor_z) {

			min_zazor_z = fabs(zpos[i + 1] - zpos[i]);
		}
	}
	min_zazor_z *= 0.25;


	const integer isize_xyz=inx*iny*inz;

	tck_int_list = new TOCKA_SHORT_INT[isize_xyz];
	// оператор new не требует проверки на null.
	//if (tck_int_list == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for evt tck_int_list...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	bool* bvisit = nullptr;
	bvisit = new bool[isize_xyz];
	// оператор new не требует проверки на null.
	//if (bvisit == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	evt = nullptr;
	evt = new int[isize_xyz];
	// оператор new не требует проверки на null.
	//if (evt == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for evt constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}
	whot_is_block = nullptr;
	whot_is_block = new int[isize_xyz];
	// оператор new не требует проверки на null.
	//if (whot_is_block == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for whot_is_block constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	if (print_log_message) {
		printf("enumerate_volume_improved_obobshenie start.\n");
	}
	Block_indexes* block_indexes = new Block_indexes[lb];
	// оператор new не требует проверки на null.
	//if (block_indexes == nullptr) {
		//printf("error in allocation memory for block_indexes in enumerate_volume_improved.\n");
		//system("pause");
		//exit(1);
	//}
	int i=0, j=0, k=0, i_1=lb-1;//lb-1

	// 08.04.2018
#pragma omp parallel for
	for (i = 0; i < lb; ++i) {
		// инициализация, на случай если блоки не будут распознаны.
		block_indexes[i].iL = -1;
		block_indexes[i].iR = -2;
		block_indexes[i].jL = -1;
		block_indexes[i].jR = -2;
		block_indexes[i].kL = -1;
		block_indexes[i].kR = -2;	
	}

	// Погрешность бывает абсолютная и относительная.
	// Вещественные числа в ЭВМ представляются с конечной точностью.
	// Лучше использовать относительную погрешность в 0.15%.
	//const doublereal otnositelnaq_tolerance_eps = 0.0035; // 0.15% 0.0015; 0.35% 27.10.2018

	//for (i = 0; i < lb; ++i) {
	for (i = lb-1; i >=0 ; i--) {
		// С учетом приоритетов, мы не обрабатываем дважды одно и тоже, и экономим время.

		//if (b[i].g.itypegeom == PRISM) {

		// polygon (b[i].g.itypegeom == POLYGON)
		// Значения прямоугольной призмы xS, xE, yS, yE, zS, zE - хранят окаймляющую
		// полигон прямоугольную призму, что позволит проверять принадлежность точки полигону
		// только для ячеек сетки находящихся внутри данной прямоугольной призмы, что сильно 
		// ускоряет обработку.
		//if ((b[i].g.itypegeom == PRISM)||(b[i].g.itypegeom == POLYGON) || (b[i].g.itypegeom == CYLINDER)) 
		{
			//  (b[i].g.itypegeom == CYLINDER) cylinder

			doublereal x4 = b[i].g.xS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE) || (b[i].g.iPlane == XZ_PLANE))) {
				x4 = b[i].g.xC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == YZ_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.xC;
				}
				else {
					x4 = b[i].g.xC+ b[i].g.Hcyl;
				}
			}
			bool bfound = false;
			/*
			for (j = 0; j <= inx; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (xpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].iL = j;
						bfound = true;
						break;
					}
					else {
						// Абсолютная погрешность.
						if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX) {
							printf("ERROR: if (fabs(100 * (xpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps)\n");
							system("PAUSE");
							block_indexes[i_1].iL = j;
							bfound = true;
							break;
						}
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4)) {
						block_indexes[i_1].iL = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= inx; ++j_63) {
					if ((!bfound) && (fabs(x4 - xpos[j_63])< min_zazor_x /*shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)*/))
					{

#pragma omp critical
						{

							if (!bfound) {

								// Точное совпадение !!!
								block_indexes[i_1].iL = j_63;
								bfound = true;
								//break;
							}
						}
					}
				}
			}
			doublereal tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= inx; ++j_63) {
					if (fabs(x4 - xpos[j_63]) < tzazor) {

#pragma omp critical
						{
							if (fabs(x4 - xpos[j_63]) < tzazor) {

								tzazor = fabs(x4 - xpos[j_63]);
								block_indexes[i_1].iL = j_63;
								bfound = true;
							}
						}
					}
				}
			}

			if (!bfound) {
#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= inx; ++j_63) {
					if ((!bfound) && (x4 < xpos[j_63])) {

#pragma omp critical
						{
							if (!bfound) {
								// Нет точного совпаднения первая встреча.
								block_indexes[i_1].iL = j_63;
								bfound = true;
								//break;
							}
						}
					}
				}
			}
			x4 = b[i].g.xE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE) || (b[i].g.iPlane == XZ_PLANE))) {
				x4 = b[i].g.xC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == YZ_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.xC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.xC;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= inx; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (xpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].iR = j;
						bfound = true;
						break;
					}
					else {
						// Абсолютная погрешность.
						if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX) {
							printf("ERROR: if (fabs(100 * (xpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps)\n");
							system("PAUSE");
							block_indexes[i_1].iR = j;
							bfound = true;
							break;
						}
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4)) {
						block_indexes[i_1].iR = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (integer j_63 = inx; j_63 >= 0; --j_63) {
					if ((!bfound) && (fabs(x4 - xpos[j_63])< min_zazor_x /*shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)*/)) {

#pragma omp critical
						{
							if (!bfound) {
								// Точное совпадение !!!
								block_indexes[i_1].iR = j_63;
								bfound = true;
								//break;
							}
						}
					}
				}
			}
			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= inx; ++j_63) {
					if (fabs(x4 - xpos[j_63]) < tzazor) {

#pragma omp critical
						{

							if (fabs(x4 - xpos[j_63]) < tzazor) {
								tzazor = fabs(x4 - xpos[j_63]);
								block_indexes[i_1].iR = j_63;
								bfound = true;
							}
						}
					}
				}
			}

			if (!bfound) {
				//std::cout << "Error incomming xE\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");

#pragma omp parallel for
				for (integer j_63 = inx; j_63 >= 0; --j_63) {
					if ((!bfound) && (x4 > xpos[j_63])) {

#pragma omp critical
						{

							if (!bfound) {
								// Нет точного совпаднения первая встреча.
								block_indexes[i_1].iR = j_63;
								bfound = true;
								//break;
							}
						}
					}
				}
			}
			x4 = b[i].g.yS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE) || (b[i].g.iPlane == YZ_PLANE))) {
				x4 = b[i].g.yC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.yC;
				}
				else {
					x4 = b[i].g.yC + b[i].g.Hcyl;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= iny; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (ypos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].jL = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4)) {
						block_indexes[i_1].jL = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= iny; ++j_63) {
					if ((!bfound) && (fabs(x4 - ypos[j_63])< min_zazor_y /*shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)*/)) {

#pragma omp critical
						{
							if (!bfound) {

								// Точное совпадение !!!
								block_indexes[i_1].jL = j_63;
								bfound = true;
								//break;
							}
						}
					}
				}
			}

			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.

#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= iny; ++j_63) {
					if (fabs(x4 - ypos[j_63]) < tzazor) {

#pragma omp critical
						{
							if (fabs(x4 - ypos[j_63]) < tzazor) {
								tzazor = fabs(x4 - ypos[j_63]);
								block_indexes[i_1].jL = j_63;
								bfound = true;
							}
						}

					}
				}
			}
			if (!bfound) {
				//std::cout << "Error incomming yS\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");

#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= iny; ++j_63) {
					if ((!bfound) && (x4 < ypos[j_63])) {

#pragma omp critical
						{
							if (!bfound) {
								// Нет точного совпаднения первая встреча.
								block_indexes[i_1].jL = j_63;
								bfound = true;
								//break;
							}

						}
					}
				}
			}
			x4 = b[i].g.yE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE) || (b[i].g.iPlane == YZ_PLANE))) {
				x4 = b[i].g.yC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.yC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.yC ;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= iny; ++j) {
				
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100*(ypos[j] - x4)/fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].jR = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4)) {
						block_indexes[i_1].jR = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (integer j_63 = iny; j_63 >= 0; --j_63) {
					if ((!bfound) && (fabs(x4 - ypos[j_63])< min_zazor_y /*shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)*/)) {
					
#pragma omp critical
						{

							if (!bfound) {

								// Точное совпадение !!!
								block_indexes[i_1].jR = j_63;
								bfound = true;
								//break;
							}
						}
					}
				}
			}
			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= iny; ++j_63) {
					if (fabs(x4 - ypos[j_63]) < tzazor) {

#pragma omp critical
						{
							if (fabs(x4 - ypos[j_63]) < tzazor) {

								tzazor = fabs(x4 - ypos[j_63]);
								block_indexes[i_1].jR = j_63;
								bfound = true;
							}
						}
					}
				}
			}
			if (!bfound) {
				//std::cout << "Error incomming yE\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");

#pragma omp parallel for
				for (integer j_63 = iny; j_63 >= 0; --j_63) {
					if ((!bfound) && (x4 > ypos[j_63])) {

#pragma omp critical
						{
							if (!bfound) {

								// Нет точного совпаднения первая встреча.
								block_indexes[i_1].jR = j_63;
								bfound = true;
								//break;
							}

						}
					}
				}
			}
			x4 = b[i].g.zS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ_PLANE) || (b[i].g.iPlane == YZ_PLANE))) {
				x4 = b[i].g.zC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.zC;
				}
				else {
					x4 = b[i].g.zC + b[i].g.Hcyl;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= inz; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (zpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].kL = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4)) {
						block_indexes[i_1].kL = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {
#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= inz; ++j_63) {
					if ((!bfound) && (fabs(x4 - zpos[j_63])< min_zazor_z /*shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)*/)) {
						
#pragma omp critical
						{

							if (!bfound) {

								// Точное совпадение !!!
								block_indexes[i_1].kL = j_63;
								bfound = true;
								//break;
							}

						}

					}
				}
			}
			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= inz; ++j_63) {
					if (fabs(x4 - zpos[j_63]) < tzazor) {

#pragma omp critical
						{

							if (fabs(x4 - zpos[j_63]) < tzazor) {

								tzazor = fabs(x4 - zpos[j_63]);
								block_indexes[i_1].kL = j_63;
								bfound = true;

							}

						}
					}
				}
			}
			if (!bfound) {

				//std::cout << "Error incomming zS\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");

#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= inz; ++j_63) {
					if ((!bfound) && (x4 < zpos[j_63])) {

#pragma omp critical
						{

							if (!bfound) {

								// Нет точного совпаднения первая встреча.
								block_indexes[i_1].kL = j_63;
								bfound = true;
								//break;
							}
						}
					}
				}
			}
			x4 = b[i].g.zE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ_PLANE) || (b[i].g.iPlane == YZ_PLANE))) {
				x4 = b[i].g.zC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.zC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.zC;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= inz; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (zpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].kR = j;
						bfound = true;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4)) {
						block_indexes[i_1].kR = j;
						bfound = true;
						break;
					}
				}
			}
			*/
			if (!bfound) {

#pragma omp parallel for
				for (integer j_63 = inz; j_63 >= 0; --j_63) {
					if ((!bfound) && (fabs(x4 - zpos[j_63])< min_zazor_z /* shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)*/)) {

#pragma omp critical
						{

							if (!bfound) {

								// Точное совпадение !!!
								block_indexes[i_1].kR = j_63;
								bfound = true;
								//break;
							}

						}
					}
				}
			}
			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
#pragma omp parallel for
				for (integer j_63 = 0; j_63 <= inz; ++j_63) {
					if (fabs(x4 - zpos[j_63]) < tzazor) {

#pragma omp critical
						{

							if (fabs(x4 - zpos[j_63]) < tzazor) {

								tzazor = fabs(x4 - zpos[j_63]);
								block_indexes[i_1].kR = j_63;
								bfound = true;
							}
						}
					}
				}
			}
			if (!bfound) {

				//std::cout << "Error incomming zE\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");
#pragma omp parallel for
				for (integer j_63 = inz; j_63 >= 0; --j_63) {
					if ((!bfound) && (x4 > zpos[j_63])) {

#pragma omp critical
						{
							if (!bfound) {

								// Нет точного совпаднения первая встреча.
								block_indexes[i_1].kR = j_63;
								bfound = true;
								//break;
							}
						}
					}
				}
			}


			if ((block_indexes[i_1].iL >= 0) &&
				(block_indexes[i_1].iR >= 0) &&
				(block_indexes[i_1].iL >= block_indexes[i_1].iR)) {
				printf("enumerate_volume_improved_obobshenie function\n");
				printf("violation of the order block_indexes\n");
				std::cout << b[i].name << std::endl;
				std::cout << "xS=" << b[i].g.xS << " xE=" << b[i].g.xE << std::endl;
				printf("i=%d iL=%lld iR=%lld\n", i, block_indexes[i_1].iL,
					block_indexes[i_1].iR);
				system("pause");
				for (int i92 = 0; i92 <= inx; i92++) {
					std::cout << i92 << " " << xpos[i92] << " \n";
					system("pause");
				}
			}

			if ((block_indexes[i_1].iL >= 0) &&
				(block_indexes[i_1].iR >= 0) &&
				(block_indexes[i_1].iL + 1 == block_indexes[i_1].iR)) {
				// Всего одна клетка на блок.
				//iX_one_CELL_count_statistic++;
			}

			if ((block_indexes[i_1].jL >= 0) &&
				(block_indexes[i_1].jR >= 0) &&
				(block_indexes[i_1].jL >= block_indexes[i_1].jR)) {
				printf("enumerate_volume_improved_obobshenie function\n");
				printf("violation of the order block_indexes\n");
				std::cout << b[i].name << std::endl;
				std::cout << "yS=" << b[i].g.yS << " yE=" << b[i].g.yE << std::endl;
				printf("i=%d jL=%lld jR=%lld\n", i, block_indexes[i_1].jL,
					block_indexes[i_1].jR);
				system("pause");
				for (int i92 = 0; i92 <= iny; i92++) {
					std::cout << i92 << " " << ypos[i92] << " \n";
					system("pause");
				}
			}


			if ((block_indexes[i_1].jL >= 0) &&
				(block_indexes[i_1].jR >= 0) &&
				(block_indexes[i_1].jL + 1 == block_indexes[i_1].jR)) {
				// Всего одна клетка на блок.
				//iY_one_CELL_count_statistic++;
			}

			if ((block_indexes[i_1].kL >= 0) &&
				(block_indexes[i_1].kR >= 0) &&
				(block_indexes[i_1].kL >= block_indexes[i_1].kR)) {
				printf("enumerate_volume_improved_obobshenie function\n");
				printf("violation of the order block_indexes\n");
				std::cout << b[i].name << std::endl;				
				std::cout << "zS=" << b[i].g.zS << " zE=" << b[i].g.zE << std::endl;
				printf("i=%d kL=%lld kR=%lld\n", i, block_indexes[i_1].kL,
					block_indexes[i_1].kR);
				system("pause");
				for (int i92 = 0; i92 <= inz; i92++) {
					std::cout << i92 << " " << zpos[i92] << " \n";
					system("pause");
				}
			}


			if ((block_indexes[i_1].kL >= 0) &&
				(block_indexes[i_1].kR >= 0) &&
				(block_indexes[i_1].kL + 1 == block_indexes[i_1].kR)) {
				// Всего одна клетка на блок.
				//iZ_one_CELL_count_statistic++;
			}

			if (0 && (i == 116)) {
				// debug
				printf("i=%d iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
					block_indexes[i_1].iL,
					block_indexes[i_1].iR, block_indexes[i_1].jL,
					block_indexes[i_1].jR, block_indexes[i_1].kL,
					block_indexes[i_1].kR);
				printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[i].g.xS, b[i].g.xE,
					b[i].g.yS, b[i].g.yE, b[i].g.zS, b[i].g.zE);
				system("pause");
			}

			if ((block_indexes[i_1].iL < 0) ||
				(block_indexes[i_1].iR < 0) ||
				(block_indexes[i_1].jL < 0) ||
				(block_indexes[i_1].jR < 0) ||
				(block_indexes[i_1].kL < 0) ||
				(block_indexes[i_1].kR < 0)) {
				printf("enumerate_volume_improved_obobshenie function\n");
				printf("i=%d iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
					block_indexes[i_1].iL,
					block_indexes[i_1].iR, block_indexes[i_1].jL,
					block_indexes[i_1].jR, block_indexes[i_1].kL,
					block_indexes[i_1].kR);
				printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[i].g.xS, b[i].g.xE,
					b[i].g.yS, b[i].g.yE, b[i].g.zS, b[i].g.zE);
				printf("cabinet: xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[0].g.xS, b[0].g.xE,
					b[0].g.yS, b[0].g.yE, b[0].g.zS, b[0].g.zE);
				printf("ERROR: may be your geometry out of cabinet...\n");
				system("pause");
			}


			//i_1++;
			i_1--;
		}

	}

	// Тотальная проверка block_indexes

	/*for (int i94 = 0; i94 < lb; i94++) {

		std::cout << b[i94].name << " ";
		std::cout << block_indexes[i94].iL << " ";
		std::cout << block_indexes[i94].iR << " ";
		std::cout << block_indexes[i94].jL << " ";
		std::cout << block_indexes[i94].jR << " ";
		std::cout << block_indexes[i94].kL << " ";
		std::cout << block_indexes[i94].kR << " ";
		std::cout << b[i94].g.xS<<" " << b[i94].g.xE << " " << b[i94].g.yS << " " << b[i94].g.yE << " " << b[i94].g.zS << " " << b[i94].g.zE << " "	<< std::endl;

	}
	system("pause");*/

	// Количество проходов существенно сократилось и в итоге это приводит к существенному
	// увеличению быстродействия.
	int m7= lb-1, m8;//  0 m7= lb - 1,
//#pragma omp parallel for
	//for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
		//integer iP = i1 + j1 * inx + k1 * inx*iny;
#pragma omp parallel for
	for (integer iP=0; iP < isize_xyz; ++iP) {
		evt[iP] = -1;
		bvisit[iP] = false;
	}

// integer iP_k1, iP_j1, iP;

   // for (m8 = 0; m8 < lb; m8++) {
	for (m8 = lb-1; m8 >= 0; --m8) {
		// Сначала идет блок с наивысшим приоритетом, мы не рассматриваем ячейки дважды.
		if (b[m8].g.itypegeom == PRISM) {
#pragma omp parallel for
			for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; ++k1) 
			{	
				integer iP_k1 = k1 * inx * iny;
				if ((k1 < 0) || (k1 > inz)) {
					// ERROR
					printf("ERROR PRISM\n");
					//printf("inx=%lld iny=%lld inz=%lld \n", inx, iny, inz);
					std::cout << "inx=" << inx << " iny=" << iny << " inz=" << inz << std::endl;
					printf(" k1=%lld \n",  k1);
					printf("iP_k1=%lld m8=%d", iP_k1, m8);
					printf("kL=%lld kR=%lld\n",  block_indexes[m7].kL, block_indexes[m7].kR);
					system("PAUSE");
				}
				for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; ++j1) 
				{
					integer iP_j1 = j1 * inx + iP_k1;
					if ((j1 < 0) || (j1 > iny)) {
						// ERROR
						printf("ERROR PRISM\n");
						//printf("inx=%lld iny=%lld inz=%lld \n", inx, iny, inz);
						std::cout << "inx=" << inx << " iny=" << iny << " inz=" << inz << std::endl;
						printf(" j1=%lld k1=%lld\n", j1, k1);
						printf("iP_j1=%lld m8=%d", iP_j1, m8);
						printf("jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
						system("PAUSE");
					}
					for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; ++i1)
					{
						integer iP = i1 + iP_j1;

						if ((i1 < 0) || (i1 > inx)) {
							// ERROR
							printf("ERROR PRISM\n");
							//printf("inx=%lld iny=%lld inz=%lld \n", inx, iny, inz);
							std::cout << "inx=" << inx << " iny=" << iny << " inz=" << inz << std::endl;
							printf("i1=%lld j1=%lld k1=%lld \n", i1, j1, k1);
							printf("iP=%lld m8=%d", iP, m8);
							printf("iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
							system("PAUSE");
						}


						if (bvisit[iP] == false)
						{

							bvisit[iP] = true;

							switch (iflag) {
							case TEMPERATURE:
								if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
									evt[iP] = -1;
								}
								else {
									evt[iP] = m8;
								}
								break;
							case HYDRODINAMIC:
								if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
									evt[iP] = -1;
								}
								else {
									evt[iP] = m8;
								}
								break;
							}
						}
					}
				}
			}
			m7--;
			
		}
		else if (b[m8].g.itypegeom == CYLINDER) {

			// как был сформирован призматический объект для цилиндра ? 
			// Надо также сократить число проверяемых точек.
			// Cylinder
			//for (integer i1 = 0; i1 < inx; ++i1) for (integer j1 = 0; j1 < iny; ++j1) for (integer k1 = 0; k1 < inz; ++k1) {
//#pragma omp parallel for
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; ++i1) 
				for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; ++j1)
					for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; ++k1) {
				integer  iP = i1 + j1 * inx + k1 * inx*iny;
				
				if ((i1 < 0) || (i1 > inx) || (j1 < 0) || (j1 > iny) || (k1 < 0) || (k1 > inz)) {
					// ERROR
					printf("ERROR CYLINDER\n");
					printf("iplane=%lld", b[m8].g.iPlane);
					printf("xC=%e yC=%e zC=%e Hcyl=%e\n", b[m8].g.xC, b[m8].g.yC, b[m8].g.zC, b[m8].g.Hcyl);
					printf("Rin=%e Rout=%e\n", b[m8].g.R_in_cyl, b[m8].g.R_out_cyl);
					//printf("inx=%lld iny=%lld inz=%lld \n", inx, iny, inz);
					std::cout << "inx=" << inx << " iny=" << iny << " inz=" << inz << std::endl;
					printf("i1=%lld j1=%lld k1=%lld \n", i1, j1, k1);
					printf("iP=%lld m8=%d", iP, m8);
					printf("iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
					system("PAUSE");
				}


				if (bvisit[iP] == false)
				{

					TOCHKA p;
					p.x = 0.5*(xpos[i1] + xpos[i1 + 1]);
					p.y = 0.5*(ypos[j1] + ypos[j1 + 1]);
					p.z = 0.5*(zpos[k1] + zpos[k1 + 1]);

					switch (b[m8].g.iPlane) {
					case XY_PLANE:
						if (fabs(b[m8].g.R_in_cyl) < 1.0e-36) {
							if ((p.z > b[m8].g.zC) && (p.z < b[m8].g.zC + b[m8].g.Hcyl)) {
								if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.yC - p.y)*(b[m8].g.yC - p.y)) < b[m8].g.R_out_cyl) {

									bvisit[iP] = true;

									switch (iflag) {
									case TEMPERATURE:
										if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
											evt[iP] = -1;
										}
										else {
											evt[iP] = m8;
										}
										break;
									case HYDRODINAMIC:
										if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
											evt[iP] = -1;
										}
										else {
											evt[iP] = m8;
										}
										break;
									}
								}
							}
						}
						else {
							if ((p.z > b[m8].g.zC) && (p.z < b[m8].g.zC + b[m8].g.Hcyl)) {
								if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.yC - p.y)*(b[m8].g.yC - p.y)) < b[m8].g.R_out_cyl) {
									if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.yC - p.y)*(b[m8].g.yC - p.y)) > b[m8].g.R_in_cyl) {

										bvisit[iP] = true;

										switch (iflag) {
										case TEMPERATURE:
											if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										case HYDRODINAMIC:
											if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										}
									}
								}
							}
						}
						break;
					case XZ_PLANE:
						if (fabs(b[m8].g.R_in_cyl) < 1.0e-36) {
							if ((p.y > b[m8].g.yC) && (p.y < b[m8].g.yC + b[m8].g.Hcyl)) {
								if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {

									bvisit[iP] = true;

									switch (iflag) {
									case TEMPERATURE:
										if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
											evt[iP] = -1;
										}
										else {
											evt[iP] = m8;
										}
										break;
									case HYDRODINAMIC:
										if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
											evt[iP] = -1;
										}
										else {
											evt[iP] = m8;
										}
										break;
									}
								}
							}
						}
						else {
							if ((p.y > b[m8].g.yC) && (p.y < b[m8].g.yC + b[m8].g.Hcyl)) {
								if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {
									if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) > b[m8].g.R_in_cyl) {

										bvisit[iP] = true;

										switch (iflag) {
										case TEMPERATURE:
											if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										case HYDRODINAMIC:
											if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										}
									}
								}
							}
						}
						break;
					case YZ_PLANE:
						if (fabs(b[m8].g.R_in_cyl) < 1.0e-36) {
							if ((p.x > b[m8].g.xC) && (p.x < b[m8].g.xC + b[m8].g.Hcyl)) {
								if (sqrt((b[m8].g.yC - p.y)*(b[m8].g.yC - p.y) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {

									bvisit[iP] = true;

									switch (iflag) {
									case TEMPERATURE:
										if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
											evt[iP] = -1;
										}
										else {
											evt[iP] = m8;
										}
										break;
									case HYDRODINAMIC:
										if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
											evt[iP] = -1;
										}
										else {
											evt[iP] = m8;
										}
										break;
									}
								}
							}
						}
						else {
							if ((p.x > b[m8].g.xC) && (p.x < b[m8].g.xC + b[m8].g.Hcyl)) {
								if (sqrt((b[m8].g.yC - p.y)*(b[m8].g.yC - p.y) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {
									if (sqrt((b[m8].g.yC - p.y)*(b[m8].g.yC - p.y) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) > b[m8].g.R_in_cyl) {
										bvisit[iP] = true;

										switch (iflag) {
										case TEMPERATURE:
											if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										case HYDRODINAMIC:
											if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										}
									}
								}
							}
						}
						break;
					}
				}
			}
			m7--;
		}
		else if (b[m8].g.itypegeom == POLYGON) {
			
			// polygon
			// Мы сокращаем число проверяемых точек 
			// рассматривая только точки внутри окаймляющей прямоугольной призмы.

			

            #pragma omp parallel for
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; ++i1) 
				for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; ++j1) 
					for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; ++k1) {

				integer  iP = i1 + j1 * inx + k1 * inx*iny;

				if ((i1 < 0) || (i1 > inx) || (j1 < 0) || (j1 > iny) || (k1 < 0) || (k1 > inz)) {
					// ERROR
					printf("ERROR POLYGON\n");
					//printf("inx=%lld iny=%lld inz=%lld \n", inx, iny, inz);
					std::cout << "inx=" << inx << " iny=" << iny << " inz=" << inz << std::endl;
					printf("i1=%lld j1=%lld k1=%lld \n", i1, j1, k1);
					printf("iP=%lld m8=%d", iP, m8);
					printf("iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
					system("PAUSE");
				}

				if (bvisit[iP] == false)
				{

					//for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
					TOCHKA p;
					p.x = 0.5*(xpos[i1] + xpos[i1 + 1]);
					p.y = 0.5*(ypos[j1] + ypos[j1 + 1]);
					p.z = 0.5*(zpos[k1] + zpos[k1 + 1]);

					int k74 = -1;
					if (in_polygon(p, b[m8].g.nsizei, b[m8].g.xi, b[m8].g.yi, b[m8].g.zi, b[m8].g.hi, b[m8].g.iPlane_obj2, k74, m8)) {
						//printf("i1=%d j1=%d k1=%d inx*iny*inz=%d\n",i1,j1,k1, inx*iny*inz);
						//printf("iL=%d iR=%d jL=%d jR=%d kL=%d kR=%d\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);

						bvisit[iP] = true;

						switch (iflag) {
						case TEMPERATURE:
							if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
								evt[iP] = -1;
							}
							else {
								evt[iP] = m8;
							}
							break;
						case HYDRODINAMIC:
							if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
								evt[iP] = -1;
							}
							else {
								evt[iP] = m8;
							}
							break;
						}
					}

				}
			}
			m7--;
		}
		else if (b[m8].g.itypegeom == CAD_STL) {

			// CAD_STL
			// Мы сокращаем число проверяемых точек 
			// рассматривая только точки внутри окаймляющей прямоугольной призмы.



#pragma omp parallel for
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; ++i1)
				for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; ++j1)
					for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; ++k1) {

				integer  iP = i1 + j1 * inx + k1 * inx * iny;

				if ((i1 < 0) || (i1 > inx) || (j1 < 0) || (j1 > iny) || (k1 < 0) || (k1 > inz)) {
					// ERROR
					printf("ERROR CAD_STL\n");
					//printf("inx=%lld iny=%lld inz=%lld \n", inx, iny, inz);
					std::cout << "inx=" << inx << " iny=" << iny << " inz=" << inz << std::endl;
					printf("i1=%lld j1=%lld k1=%lld \n", i1, j1, k1);
					printf("iP=%lld m8=%d", iP, m8);
					printf("iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
					system("PAUSE");
				}

				if (bvisit[iP] == false)
				{

					//for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
					TOCHKA p;
					p.x = 0.5 * (xpos[i1] + xpos[i1 + 1]);
					p.y = 0.5 * (ypos[j1] + ypos[j1 + 1]);
					p.z = 0.5 * (zpos[k1] + zpos[k1 + 1]);

					int k74 = -1;

					if (b[m8].g.in_CAD_STL_check(p, k74, m8))
					{
						//printf("i1=%d j1=%d k1=%d inx*iny*inz=%d\n",i1,j1,k1, inx*iny*inz);
						//printf("iL=%d iR=%d jL=%d jR=%d kL=%d kR=%d\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);

						bvisit[iP] = true;

						switch (iflag) {
						case TEMPERATURE:
							if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
								evt[iP] = -1;
							}
							else {
								evt[iP] = m8;
							}
							break;
						case HYDRODINAMIC:
							if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
								evt[iP] = -1;
							}
							else {
								evt[iP] = m8;
							}
							break;
						}
					}

				}
			}
			m7--;
		}
		
	}

	//if (bvisit != nullptr) {
	// Оператор delete может быть вызван повторно и даже к null указателю.
		delete[] bvisit;
		bvisit = nullptr;
	//}
	
		if (print_log_message) {
			printf("enumerate_volume_improved 80 procent.\n");
		}
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	int l = 1;
	

	const int iSIZE = inz * iny*inx;
	int* lind = new int[iSIZE];

#pragma omp parallel for
	for (int k_63 = 0; k_63 < iSIZE; ++k_63)
	{
		lind[k_63] = -1;
	}

	for (k = 0; k < inz; ++k)
	{
		integer iP_k = k * inx * iny;
		for (j = 0; j < iny; ++j)
		{
			integer iP_j = j * inx + iP_k;
			for (i = 0; i < inx; ++i)
			{
				integer  iP = i + iP_j;
				if (evt[iP] > -1) {
					lind[iP] = l-1;
					l++;
				}
			}
		}
	}


	// integer iP_k, iP_j;
#pragma omp parallel for
	for (int k_63 = 0; k_63 < inz; ++k_63)
	{
		integer iP_k = k_63 * inx * iny;
		for (int j_63 = 0; j_63 < iny; ++j_63)
		{
			integer iP_j = j_63 * inx + iP_k;
			for (int i_63 = 0; i_63 < inx; ++i_63)
			{
				integer  iP = i_63 + iP_j;
				if (evt[iP] > -1) {
					int ib_63 = evt[iP]; // номер блока был сохранён ранее.
								 // Это очень нужно для записи репорта.
					whot_is_block[lind[iP]] = ib_63; // номер блока которому принадлежит точка (p.x,p.y,p.z).
					tck_int_list[lind[iP]].i = i_63;
					tck_int_list[lind[iP]].j = j_63;
					tck_int_list[lind[iP]].k = k_63;
					evt[iP] = lind[iP]+1;
					//l++;
				}
				else
				{   // не принадлежит расчётной области
					evt[iP] = 0;
				}
			}
		}
	}	
	maxelm = l - 1;

	delete[] lind;

	//if (block_indexes != nullptr) {
	// Оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] block_indexes;
		block_indexes = nullptr;
	//}
		if (print_log_message) {
			printf("enumerate_volume_improved_obobshenie end.\n");
		}

	/*
	bool inDomain = false;
	// нумерация контрольных объёмов начинается с единицы
	// если контрольный объём не принадлежит расчётной области то ставится 0.
	for (i = 0; i<inx; ++i) for (j = 0; j<iny; ++j) for (k = 0; k<inz; ++k) {
	p.x = 0.5*(xpos[i] + xpos[i + 1]);
	p.y = 0.5*(ypos[j] + ypos[j + 1]);
	p.z = 0.5*(zpos[k] + zpos[k + 1]);
	switch (iflag) {
	case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
	case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
	}
	if (inDomain) {
	// принадлежит расчётной области
	evt[i + j*inx + k*inx*iny] = l;
	// Это очень нужно для записи репорта.
	whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
	l++;
	}
	else
	{   // не принадлежит расчётной области
	evt[i + j*inx + k*inx*iny] = 0;
	}
	}

	// maxelm - число контрольных объёмов принадлежащих расчётной области
	maxelm = l - 1;
	*/

#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif

} // enumerate_volume_improved_obobshenie


// глобальная нумерация контрольных объёмов
// для задач теплопроводности
void enumerate_volume(int* &evt, int &maxelm, integer iflag,
	doublereal* &xpos, doublereal* &ypos, doublereal* &zpos, int* &whot_is_block,
	int inx, int iny, int inz, BLOCK* &b, int lb, 
	int lu, UNION* &my_union, integer &iunion_id_p1, TOCKA_SHORT_INT* &tck_int_list,
	WALL* &w, int &lw, SOURCE* &s, int &ls) {

	

	if (lu==0) {
		// 29.12.2017
		// Работает также в случае полигонов и цилиндров.
		enumerate_volume_improved_obobshenie(evt, maxelm, iflag, xpos, ypos, zpos, whot_is_block, inx, iny, inz, b, lb, tck_int_list,w,lw,s,ls);
	}
	else {

		if (lu > 0) {

			// 25.04.2018
			// Работает и для АСЕБЛЕСОВ.

			// Присутствуют также и цилиндры.

			tck_int_list=new TOCKA_SHORT_INT[inx*iny*inz];
			if (tck_int_list == nullptr) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment for evt tck_int_list...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}

			evt = nullptr;
			evt = new int[inx*iny*inz];
			if (evt == nullptr) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment for evt constr struct...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
			whot_is_block = nullptr;
			whot_is_block = new int[inx*iny*inz];
			if (whot_is_block == nullptr) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment for whot_is_block constr struct...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
			TOCHKA p;
			// нумерация в evt начиная с единицы.
			// если не принадлежит расчётной области то стоит 0.
			int l = 1;
			int ib;
			//integer i, j, k;
			bool inDomain = false;
			// нумерация контрольных объёмов начинается с единицы
			// если контрольный объём не принадлежит расчётной области то ставится 0.
			
			// integer iP_k, iP_j, iP;
			for (int k = 0; k < inz; ++k)
			{
				integer iP_k = k * inx * iny;
				for (int j = 0; j < iny; ++j)
				{
					integer iP_j = j * inx + iP_k;
					for (int i = 0; i < inx; ++i)
					{
						integer iP = i + iP_j;

						p.x = 0.5 * (xpos[i] + xpos[i + 1]);
						p.y = 0.5 * (ypos[j] + ypos[j + 1]);
						p.z = 0.5 * (zpos[k] + zpos[k + 1]);
						switch (iflag) {
						case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
						case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
						}
						if (iunion_id_p1 == 0) {
							for (int iu = 0; iu < lu; iu++) {
								if (my_union[iu].active) {
									if (my_union[iu].iunion_parent == -1) {
										if ((p.x > my_union[iu].xS) && (p.x < my_union[iu].xE) &&
											(p.y > my_union[iu].yS) && (p.y < my_union[iu].yE) &&
											(p.z > my_union[iu].zS) && (p.z < my_union[iu].zE)) {
											// Если рассматривать кабинет с асемблесами, то область внутри
											// Асемблесов есть Hollow блок.
											inDomain = false;
										}
									}
								}
							}
						}
						else {

							for (int iu = 0; iu < lu; iu++) {

								if (my_union[iu].active) {

									if (my_union[iu].iunion_parent == iunion_id_p1 - 1) {

										if ((p.x > my_union[iu].xS) && (p.x < my_union[iu].xE) &&
											(p.y > my_union[iu].yS) && (p.y < my_union[iu].yE) &&
											(p.z > my_union[iu].zS) && (p.z < my_union[iu].zE)) {
											// Если рассматривать кабинет с асемблесами, то область внутри
											// Асемблесов есть Hollow блок.
											//inDomain = false;
											// Оставляем так всё верно.
											if (inDomain) {
												//	printf("in our assembles\n");
												//	system("PAUSE");

												// Если рассматривать юнион  my_union[iunion_id_p1-1].name с асемблесами, то область внутри
												// Асемблесов есть Hollow блок.
												inDomain = false;
											}
										}
									}
								}

							}
						}


						if (inDomain) {
							// принадлежит расчётной области
							evt[iP] = l;
							// Это очень нужно для записи репорта.
							whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
							tck_int_list[l - 1].i = i;
							tck_int_list[l - 1].j = j;
							tck_int_list[l - 1].k = k;
							l++;
						}
						else
						{   // не принадлежит расчётной области
							evt[iP] = 0;
						}
					}
				}
			}

			// maxelm - число контрольных объёмов принадлежащих расчётной области
			maxelm = l - 1;
		}

		
		else {
			bool b_only_Prism_object = true;
			for (integer m7 = 0; m7 < lb; m7++) {
				if (b[m7].g.itypegeom != PRISM) {
					b_only_Prism_object = false;
					break;
				}
			}

			if (b_only_Prism_object) {
				// Вся модель построена только на прямоугольных призмах.
				enumerate_volume_improved(evt, maxelm, iflag, xpos, ypos, zpos, whot_is_block, inx, iny, inz, b, lb, w, lw, s, ls);
			}
			else {
				// Присутствуют также и цилиндры.

				evt = nullptr;
				evt = new int[inx*iny*inz];
				if (evt == nullptr) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment for evt constr struct...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
				whot_is_block = nullptr;
				whot_is_block = new int[inx*iny*inz];
				if (whot_is_block == nullptr) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem: not enough memory on your equipment for whot_is_block constr struct...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
				TOCHKA p;
				// нумерация в evt начиная с единицы.
				// если не принадлежит расчётной области то стоит 0.
				int l = 1;
				int ib;
				//integer i, j, k;
				bool inDomain = false;
				// нумерация контрольных объёмов начинается с единицы
				// если контрольный объём не принадлежит расчётной области то ставится 0.
				
				// integer iP_k, iP_j, iP;

				for (int k = 0; k < inz; ++k)
				{
					integer iP_k = k * inx * iny;
					for (int j = 0; j < iny; ++j)
					{
						integer iP_j = j * inx + iP_k;
						for (int i = 0; i < inx; ++i)
						{
							integer iP = i + iP_j;
							p.x = 0.5 * (xpos[i] + xpos[i + 1]);
							p.y = 0.5 * (ypos[j] + ypos[j + 1]);
							p.z = 0.5 * (zpos[k] + zpos[k + 1]);
							switch (iflag) {
							case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
							case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
							}

							if (inDomain) {
								// принадлежит расчётной области
								evt[iP] = l;
								// Это очень нужно для записи репорта.
								whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
								tck_int_list[l - 1].i = i;
								tck_int_list[l - 1].j = j;
								tck_int_list[l - 1].k = k;
								l++;
							}
							else
							{   // не принадлежит расчётной области
								evt[iP] = 0;
							}
						}
					}
				}

				// maxelm - число контрольных объёмов принадлежащих расчётной области
				maxelm = l - 1;
			}
		}

		
	}
} // enumerate_volume

  // Инициализация evt_f. Он нужен для заливки, через которую определяется функция цвета.
  // А функция цвета нужна обязательным образом для корректировки массового баланса при нескольких
  // FLUID областях связанных лишь уравнением теплопередачи, иначе не будет сходимости.
void init_evt_f_alice_improved(int* &evt, integer iflag,
	doublereal* &xpos, doublereal* &ypos, doublereal* &zpos,
	int inx, int iny, int inz, 
	BLOCK* b, int lb, TOCKA_SHORT_INT* &tck_int_list,
	WALL* &w, int &lw, SOURCE* &s, int &ls) {

	tck_int_list = new TOCKA_SHORT_INT[inx*iny*inz];
	if (tck_int_list == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for evt tck_int_list...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	evt = nullptr;
	evt = new int[inx*iny*inz];
	if (evt == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for evt constr struct...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	//TOCHKA p;

	if (print_log_message) {
		printf("enumerate_volume_improved start.\n");
	}
	Block_indexes* block_indexes = new Block_indexes[lb];
	// оператор new не требует проверки на null.
	//if (block_indexes == nullptr) {
		//printf("error in allocation memory for block_indexes in init_evt_f_alice_improved.\n");
		//system("pause");
		//exit(1);
	//}

	// 07.08.2021 Изменение допусков т.к. на большой сетке
	// структурированой для гидродинамики в 8млн узлов не распознались входная и выходная границы.
	// Для того чтобы входная и выходная границы распознавались введены допуски min_zazor_* вместо
	// допусков shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls).
	doublereal min_zazor_x = 1.0e36;
	doublereal min_zazor_y = 1.0e36;
	doublereal min_zazor_z = 1.0e36;

	for (int i = 0; i < inx; ++i) {
		if (fabs(xpos[i + 1] - xpos[i]) < min_zazor_x) {

			min_zazor_x = fabs(xpos[i + 1] - xpos[i]);
		}
	}
	min_zazor_x *= 0.25;

	for (int i = 0; i < iny; ++i) {
		if (fabs(ypos[i + 1] - ypos[i]) < min_zazor_y) {

			min_zazor_y = fabs(ypos[i + 1] - ypos[i]);
		}
	}
	min_zazor_y *= 0.25;

	for (int i = 0; i < inz; ++i) {
		if (fabs(zpos[i + 1] - zpos[i]) < min_zazor_z) {

			min_zazor_z = fabs(zpos[i + 1] - zpos[i]);
		}
	}
	min_zazor_z *= 0.25;


	// 08.04.2018
	for (integer i = 0; i < lb; ++i) {
		// инициализация, на случай если блоки не будут распознаны.
		block_indexes[i].iL = -1;
		block_indexes[i].iR = -2;
		block_indexes[i].jL = -1;
		block_indexes[i].jR = -2;
		block_indexes[i].kL = -1;
		block_indexes[i].kR = -2;
	}

	//integer i=0, j=0, k=0;

	for (integer i = 0; i < lb; ++i) {
		integer i_1 = i;

		doublereal x4 = b[i].g.xS;
		/*for (integer j = 0; j <= inx; ++j) {
			if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].iL = j;
				break;
			}
		}*/
		bool bfound = false;

		if (!bfound) {
			for (integer j = 0; j <= inx; ++j) {
				if ((!bfound) && (fabs(x4 - xpos[j]) < min_zazor_x /*shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)*/))
				{
					// Точное совпадение !!!
					block_indexes[i_1].iL = j;
					bfound = true;
					break;
				}
			}
		}
		doublereal tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (integer  j = 0; j <= inx; ++j) {
				if (fabs(x4 - xpos[j]) < tzazor) {
					tzazor = fabs(x4 - xpos[j]);
					block_indexes[i_1].iL = j;
					bfound = true;
				}
			}
		}

		if (!bfound) {
			for (integer  j = 0; j <= inx; ++j) {
				if ((!bfound) && (x4 < xpos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].iL = j;
					bfound = true;
					break;
				}
			}
		}

		x4 = b[i].g.xE;
		/*for (integer j = 0; j <= inx; ++j) {
			if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].iR = j;
				break;
			}
		}*/
		bfound = false;

		if (!bfound) {
			for (integer j = inx; j >= 0; j--) {
				if ((!bfound) && (fabs(x4 - xpos[j]) < min_zazor_x /*shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].iR = j;
					bfound = true;
					break;
				}
			}
		}
		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (integer j = 0; j <= inx; ++j) {
				if (fabs(x4 - xpos[j]) < tzazor) {
					tzazor = fabs(x4 - xpos[j]);
					block_indexes[i_1].iR = j;
					bfound = true;
				}
			}
		}

		if (!bfound) {
			//std::cout << "Error incomming xE\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (integer j = inx; j >= 0; j--) {
				if ((!bfound) && (x4 > xpos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].iR = j;
					bfound = true;
					break;
				}
			}
		}

		x4 = b[i].g.yS;
		bfound = false;

		/*for (integer j = 0; j <= iny; ++j) {
			if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].jL = j;
				break;
			}
		}*/


		if (!bfound) {
			for (integer j = 0; j <= iny; ++j) {
				if ((!bfound) && (fabs(x4 - ypos[j]) < min_zazor_y /*shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].jL = j;
					bfound = true;
					break;
				}
			}
		}

		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (integer j = 0; j <= iny; ++j) {
				if (fabs(x4 - ypos[j]) < tzazor) {
					tzazor = fabs(x4 - ypos[j]);
					block_indexes[i_1].jL = j;
					bfound = true;
				}
			}
		}
		if (!bfound) {
			//std::cout << "Error incomming yS\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (integer j = 0; j <= iny; ++j) {
				if ((!bfound) && (x4 < ypos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].jL = j;
					bfound = true;
					break;
				}
			}
		}


		x4 = b[i].g.yE;

		bfound = false;

		/*for (integer j = 0; j <= iny; ++j) {
			if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].jR = j;
				break;
			}
		}*/
		if (!bfound) {
			for (integer  j = iny; j >= 0; j--) {
				if ((!bfound) && (fabs(x4 - ypos[j]) < min_zazor_y /*shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].jR = j;
					bfound = true;
					break;
				}
			}
		}
		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (integer j = 0; j <= iny; ++j) {
				if (fabs(x4 - ypos[j]) < tzazor) {
					tzazor = fabs(x4 - ypos[j]);
					block_indexes[i_1].jR = j;
					bfound = true;
				}
			}
		}
		if (!bfound) {
			//std::cout << "Error incomming yE\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (integer  j = iny; j >= 0; j--) {
				if ((!bfound) && (x4 > ypos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].jR = j;
					bfound = true;
					break;
				}
			}
		}

		x4 = b[i].g.zS;

		bfound = false;

		/*for (integer j = 0; j <= inz; ++j) {
			if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].kL = j;
				break;
			}
		}*/
		if (!bfound) {
			for (integer j = 0; j <= inz; ++j) {
				if ((!bfound) && (fabs(x4 - zpos[j]) < min_zazor_z /*shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].kL = j;
					bfound = true;
					break;
				}
			}
		}
		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (integer j = 0; j <= inz; ++j) {
				if (fabs(x4 - zpos[j]) < tzazor) {
					tzazor = fabs(x4 - zpos[j]);
					block_indexes[i_1].kL = j;
					bfound = true;
				}
			}
		}
		if (!bfound) {

			//std::cout << "Error incomming zS\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (integer  j = 0; j <= inz; ++j) {
				if ((!bfound) && (x4 < zpos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].kL = j;
					bfound = true;
					break;
				}
			}
		}

		x4 = b[i].g.zE;

		bfound = false;

		/*for (integer j = 0; j <= inz; ++j) {
			if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].kR = j;
				break;
			}
		}*/
		if (!bfound) {
			for (integer j = inz; j >= 0; j--) {
				if ((!bfound) && (fabs(x4 - zpos[j]) < min_zazor_z /* shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)*/)) {
					// Точное совпадение !!!
					block_indexes[i_1].kR = j;
					bfound = true;
					break;
				}
			}
		}
		tzazor = 1.0e36;
		if (!bfound) {
			// Ищем наиболее близко расположенную границу.
			for (integer j = 0; j <= inz; ++j) {
				if (fabs(x4 - zpos[j]) < tzazor) {
					tzazor = fabs(x4 - zpos[j]);
					block_indexes[i_1].kR = j;
					bfound = true;
				}
			}
		}
		if (!bfound) {

			//std::cout << "Error incomming zE\n";
			//std::cout << "position=" << x4 << "\n";
			//std::cout << b[i].name << std::endl;
			//system("pause");

			for (integer j = inz; j >= 0; j--) {
				if ((!bfound) && (x4 > zpos[j])) {
					// Нет точного совпаднения первая встреча.
					block_indexes[i_1].kR = j;
					bfound = true;
					break;
				}
			}
		}


		if ((block_indexes[i].iL >= 0) &&
			(block_indexes[i].iR >= 0) &&
			(block_indexes[i].iL >= block_indexes[i].iR)) {
			printf("init_evt_f_alice_improved  function\n");
			printf("violation of the order block_indexes\n");
			printf("i=%lld iL=%lld iR=%lld\n", i, block_indexes[i].iL,
				block_indexes[i].iR);
			system("pause");
			for (int i92 = 0; i92 <= inx; i92++) {
				std::cout << i92 << " " << xpos[i92] << " \n";
				system("pause");
			}
		}

		if ((block_indexes[i].iL >= 0) &&
			(block_indexes[i].iR >= 0) &&
			(block_indexes[i].iL + 1 == block_indexes[i].iR)) {
			// Всего одна клетка на блок.
			//iX_one_CELL_count_statistic++;
		}

		if ((block_indexes[i].jL >= 0) &&
			(block_indexes[i].jR >= 0) &&
			(block_indexes[i].jL >= block_indexes[i].jR)) {
			printf("init_evt_f_alice_improved  function\n");
			printf("violation of the order block_indexes\n");
			printf("i=%lld jL=%lld jR=%lld\n", i, block_indexes[i].jL,
				block_indexes[i].jR);
			system("pause");
			for (int i92 = 0; i92 <= iny; i92++) {
				std::cout << i92 << " " << ypos[i92] << " \n";
				system("pause");
			}
		}


		if ((block_indexes[i].jL >= 0) &&
			(block_indexes[i].jR >= 0) &&
			(block_indexes[i].jL + 1 == block_indexes[i].jR)) {
			// Всего одна клетка на блок.
			//iY_one_CELL_count_statistic++;
		}

		if ((block_indexes[i].kL >= 0) &&
			(block_indexes[i].kR >= 0) &&
			(block_indexes[i].kL >= block_indexes[i].kR)) {
			printf("init_evt_f_alice_improved  function\n");
			printf("violation of the order block_indexes\n");
			printf("i=%lld kL=%lld kR=%lld\n", i, block_indexes[i].kL,
				block_indexes[i].kR);
			system("pause");
			for (int i92 = 0; i92 <= inz; i92++) {
				std::cout << i92 << " " << zpos[i92] << " \n";
				system("pause");
			}
		}


		if ((block_indexes[i].kL >= 0) &&
			(block_indexes[i].kR >= 0) &&
			(block_indexes[i].kL + 1 == block_indexes[i].kR)) {
			// Всего одна клетка на блок.
			//iZ_one_CELL_count_statistic++;
		}

		if ((block_indexes[i].iL < 0) ||
			(block_indexes[i].iR < 0) ||
			(block_indexes[i].jL < 0) ||
			(block_indexes[i].jR < 0) ||
			(block_indexes[i].kL < 0) ||
			(block_indexes[i].kR < 0)) {
			printf("init_evt_f_alice_improved function\n");
			printf("i=%lld iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
				block_indexes[i].iL,
				block_indexes[i].iR, block_indexes[i].jL,
				block_indexes[i].jR, block_indexes[i].kL,
				block_indexes[i].kR);
			printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[i].g.xS, b[i].g.xE,
				b[i].g.yS, b[i].g.yE, b[i].g.zS, b[i].g.zE);
			printf("cabinet: xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[0].g.xS, b[0].g.xE,
				b[0].g.yS, b[0].g.yE, b[0].g.zS, b[0].g.zE);
			printf("ERROR: may be your geometry out of cabinet...\n");
			system("pause");
		}


	}

	// Количество проходов существенно сократилось и в итоге это приводит к существенному
	// увеличению быстродействия.
	//integer m7;
	const integer size_xyz = inx * iny * inz;
	for (integer i = 0; i < size_xyz; ++i) {
		evt[i] = -1;
	}

	// integer iP_k, iP_j, iP;

	for (int m7 = 0; m7 < lb; ++m7) {
		for (integer k = block_indexes[m7].kL; k < block_indexes[m7].kR; ++k)
		{
			integer iP_k = k * inx * iny;
			for (integer j = block_indexes[m7].jL; j < block_indexes[m7].jR; ++j)
			{
				integer iP_j = j * inx + iP_k;
				for (integer i = block_indexes[m7].iL; i < block_indexes[m7].iR; ++i)
				{
					integer iP = i + iP_j;
					switch (iflag) {
					case TEMPERATURE:
						if (b[m7].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
							evt[iP] = -1;
						}
						else {
							evt[iP] = m7;
						}
						break;
					case HYDRODINAMIC:
						if ((b[m7].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m7].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
							evt[iP] = -1;
						}
						else {
							evt[iP] = m7;
						}
						break;
					}

				}
			}
		}
	}
	if (print_log_message) {
		printf("enumerate_volume_improved 80 procent.\n");
	}
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	int l = 1;
	for (int k = 0; k < inz; ++k)
	{
		integer iP_k = k * inx * iny;
		for (int j = 0; j < iny; ++j)
		{
			integer iP_j = j * inx + iP_k;
			for (int i = 0; i < inx; ++i)
			{
				integer iP = i + iP_j;
				if (evt[iP] > -1) {
					//integer ib = evt[iP]; // номер блока был сохранён ранее.
					// Это очень нужно для записи репорта.
					//whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
					
					evt[iP] = l;

					tck_int_list[l - 1].i = i;
					tck_int_list[l - 1].j = j;
					tck_int_list[l - 1].k = k;
					l++;
				}
				else {
					// не принадлежит расчётной области
					evt[iP] = 0;
				}
			}
		}
	}


	delete[] block_indexes;

	if (print_log_message) {
		printf("enumerate_volume_improved end.\n");
	}
	/*
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	integer l = 1, ib;
	integer i, j, k;
	bool inDomain = false;
	// нумерация контрольных объёмов начинается с единицы
	// если контрольный объём не принадлежит расчётной области то ставится 0.
	for (i = 0; i<inx; ++i) for (j = 0; j<iny; ++j) for (k = 0; k<inz; ++k) {
	p.x = 0.5*(xpos[i] + xpos[i + 1]);
	p.y = 0.5*(ypos[j] + ypos[j + 1]);
	p.z = 0.5*(zpos[k] + zpos[k + 1]);
	switch (iflag) {
	case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
	case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
	}
	if (inDomain) {
	// принадлежит расчётной области
	evt[i + j*inx + k*inx*iny] = l;
	l++;
	}
	else
	{   // не принадлежит расчётной области
	evt[i + j*inx + k*inx*iny] = 0;
	}
	}
	*/

} // init_evt_f_alice_improved.


  // Инициализация evt_f. Он нужен для заливки, через которую определяется функция цвета.
  // А функция цвета нужна обязательным образом для корректировки массового баланса при нескольких
  // FLUID областях связанных лишь уравнением теплопередачи, иначе не будет сходимости.
void init_evt_f_alice_improved_obobshenie(int* &evt, integer iflag, 
	doublereal* xpos, doublereal* ypos, doublereal* zpos,
	int inx, int iny, int inz,
	BLOCK* b, int lb, TOCKA_SHORT_INT* &tck_int_list,
	WALL* &w, int &lw, SOURCE* &s, int &ls) {




	tck_int_list = new TOCKA_SHORT_INT[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (tck_int_list == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for evt tck_int_list...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	bool* bvisit = nullptr;
	bvisit = new bool[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (bvisit == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	evt = nullptr;
	evt = new int[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (evt == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for evt constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	//TOCHKA p;

	if (print_log_message) {
		printf("enumerate_volume_improved start.\n");
	}
	Block_indexes* block_indexes = new Block_indexes[lb];
	// оператор new не требует проверки на null.
	//if (block_indexes == nullptr) {
		//printf("error in allocation memory for block_indexes in init_evt_f_alice_improved.\n");
		//system("pause");
		//exit(1);
	//}


	// 07.08.2021 Изменение допусков т.к. на большой сетке
	// структурированой для гидродинамики в 8млн узлов не распознались входная и выходная границы.
	// Для того чтобы входная и выходная границы распознавались введены допуски min_zazor_* вместо
	// допусков shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls).
	doublereal min_zazor_x = 1.0e36;
	doublereal min_zazor_y = 1.0e36;
	doublereal min_zazor_z = 1.0e36;

	for (int i = 0; i < inx; ++i) {
		if (fabs(xpos[i + 1] - xpos[i]) < min_zazor_x) {

			min_zazor_x = fabs(xpos[i + 1] - xpos[i]);
		}
	}
	min_zazor_x *= 0.25;

	for (int i = 0; i < iny; ++i) {
		if (fabs(ypos[i + 1] - ypos[i]) < min_zazor_y) {

			min_zazor_y = fabs(ypos[i + 1] - ypos[i]);
		}
	}
	min_zazor_y *= 0.25;

	for (int i = 0; i < inz; ++i) {
		if (fabs(zpos[i + 1] - zpos[i]) < min_zazor_z) {

			min_zazor_z = fabs(zpos[i + 1] - zpos[i]);
		}
	}
	min_zazor_z *= 0.25;




	// 08.04.2018
	for (int i = 0; i < lb; ++i) {
		// инициализация, на случай если блоки не будут распознаны.
		block_indexes[i].iL = -1;
		block_indexes[i].iR = -2;
		block_indexes[i].jL = -1;
		block_indexes[i].jR = -2;
		block_indexes[i].kL = -1;
		block_indexes[i].kR = -2;
	}


	int i=0, j=0, k=0, i_1 = lb-1;

	// Погрешность бывает абсолютная и относительная.
	// Вещественные числа в ЭВМ представляются с конечной точностью.
	// Лучше использовать относительную погрешность в 0.15%.
	//const doublereal otnositelnaq_tolerance_eps = 0.0015; // 0.15%

	for (i = lb-1; i >= 0; --i) {
		//if (b[i].g.itypegeom == PRISM) {

		// polygon (b[i].g.itypegeom == POLYGON)
		// Значения прямоугольной призмы xS, xE, yS, yE, zS, zE - хранят окаймляющую
		// полигон прямоугольную призму, что позволит проверять принадлежность точки полигону
		// только для ячеек сетки находящихся внутри данной прямоугольной призмы, что сильно 
		// ускоряет обработку.
		//if ((b[i].g.itypegeom == PRISM) || (b[i].g.itypegeom == CYLINDER) || (b[i].g.itypegeom == POLYGON))
		{

			doublereal x4 = b[i].g.xS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE) || (b[i].g.iPlane == XZ_PLANE))) {
				x4 = b[i].g.xC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == YZ_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.xC;
				}
				else {
					x4 = b[i].g.xC + b[i].g.Hcyl;
				}
			}
			/*
			for (j = 0; j <= inx; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (xpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].iL = j;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4,b,lb)) {
						block_indexes[i_1].iL = j;
						break;
					}
				}
			}
			*/
			bool bfound = false;

			if (!bfound) {
				for (j = 0; j <= inx; ++j) {
					// Абсолютная погрешность.
					if (fabs(xpos[j] - x4) < min_zazor_x /*shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)*/) {
						block_indexes[i_1].iL = j;
						bfound = true;
						break;
					}
				}
			}

			doublereal tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
				for (j = 0; j <= inx; ++j) {
					if (fabs(x4 - xpos[j]) < tzazor) {
						tzazor = fabs(x4 - xpos[j]);
						block_indexes[i_1].iL = j;
						bfound = true;
					}
				}
			}

			if (!bfound) {
				for (j = 0; j <= inx; ++j) {
					if ((!bfound) && (x4 < xpos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].iL = j;
						bfound = true;
						break;
					}
				}
			}


			x4 = b[i].g.xE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE) || (b[i].g.iPlane == XZ_PLANE))) {
				x4 = b[i].g.xC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == YZ_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.xC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.xC;
				}
			}
			/*
			for (j = 0; j <= inx; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (xpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].iR = j;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4,b,lb, w, lw, s, ls)) {
						block_indexes[i_1].iR = j;
						break;
					}
				}
			}
			*/
			bfound = false;

			if (!bfound) {
				for (j = 0; j <= inx; ++j) {
					// Абсолютная погрешность.
					if (fabs(xpos[j] - x4) < min_zazor_x /*shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)*/) {
						block_indexes[i_1].iR = j;
						bfound = true;
						break;
					}
				}
			}
			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
				for (j = 0; j <= inx; ++j) {
					if (fabs(x4 - xpos[j]) < tzazor) {
						tzazor = fabs(x4 - xpos[j]);
						block_indexes[i_1].iR = j;
						bfound = true;
					}
				}
			}

			if (!bfound) {
				//std::cout << "Error incomming xE\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");

				for (j = inx; j >= 0; j--) {
					if ((!bfound) && (x4 > xpos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].iR = j;
						bfound = true;
						break;
					}
				}
			}

			x4 = b[i].g.yS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE) || (b[i].g.iPlane == YZ_PLANE))) {
				x4 = b[i].g.yC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.yC;
				}
				else {
					x4 = b[i].g.yC + b[i].g.Hcyl;
				}
			}
			/*
			for (j = 0; j <= iny; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (ypos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].jL = j;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4,b,lb, w, lw, s, ls)) {
						block_indexes[i_1].jL = j;
						break;
					}
				}
			}
			*/
			bfound = false;

			if (!bfound) {
				for (j = 0; j <= iny; ++j) {
					// Абсолютная погрешность.
					if (fabs(ypos[j] - x4) < min_zazor_y /*shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)*/) {
						block_indexes[i_1].jL = j;
						bfound = true;
						break;
					}
				}
			}

			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
				for (j = 0; j <= iny; ++j) {
					if (fabs(x4 - ypos[j]) < tzazor) {
						tzazor = fabs(x4 - ypos[j]);
						block_indexes[i_1].jL = j;
						bfound = true;
					}
				}
			}
			if (!bfound) {
				//std::cout << "Error incomming yS\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");

				for (j = 0; j <= iny; ++j) {
					if ((!bfound) && (x4 < ypos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].jL = j;
						bfound = true;
						break;
					}
				}
			}

			x4 = b[i].g.yE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE) || (b[i].g.iPlane == YZ_PLANE))) {
				x4 = b[i].g.yC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.yC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.yC;
				}
			}
			/*
			for (j = 0; j <= iny; ++j) {

				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (ypos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].jR = j;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4,b,lb, w, lw, s, ls)) {
						block_indexes[i_1].jR = j;
						break;
					}
				}
			}
			*/
			bfound = false;

			if (!bfound) {
				for (j = 0; j <= iny; ++j) {
					// Абсолютная погрешность.
					if (fabs(ypos[j] - x4) < min_zazor_y /*shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)*/) {
						block_indexes[i_1].jR = j;
						bfound = true;
						break;
					}
				}
			}

			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
				for (j = 0; j <= iny; ++j) {
					if (fabs(x4 - ypos[j]) < tzazor) {
						tzazor = fabs(x4 - ypos[j]);
						block_indexes[i_1].jR = j;
						bfound = true;
					}
				}
			}
			if (!bfound) {
				//std::cout << "Error incomming yE\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");

				for (j = iny; j >= 0; j--) {
					if ((!bfound) && (x4 > ypos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].jR = j;
						bfound = true;
						break;
					}
				}
			}

			x4 = b[i].g.zS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ_PLANE) || (b[i].g.iPlane == YZ_PLANE))) {
				x4 = b[i].g.zC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.zC;
				}
				else {
					x4 = b[i].g.zC + b[i].g.Hcyl;
				}
			}
			/*
			for (j = 0; j <= inz; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (zpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].kL = j;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4,b,lb, w, lw, s, ls)) {
						block_indexes[i_1].kL = j;
						break;
					}
				}
			}
			*/
			bfound = false;

			if (!bfound) {
				for (j = 0; j <= inz; ++j) {
					// Абсолютная погрешность.
					if (fabs(zpos[j] - x4) < min_zazor_z /*shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)*/) {
						block_indexes[i_1].kL = j;
						bfound = true;
						break;
					}
				}
			}
			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
				for (j = 0; j <= inz; ++j) {
					if (fabs(x4 - zpos[j]) < tzazor) {
						tzazor = fabs(x4 - zpos[j]);
						block_indexes[i_1].kL = j;
						bfound = true;
					}
				}
			}
			if (!bfound) {

				//std::cout << "Error incomming zS\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");

				for (j = 0; j <= inz; ++j) {
					if ((!bfound) && (x4 < zpos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].kL = j;
						bfound = true;
						break;
					}
				}
			}


			x4 = b[i].g.zE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ_PLANE) || (b[i].g.iPlane == YZ_PLANE))) {
				x4 = b[i].g.zC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY_PLANE))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.zC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.zC;
				}
			}
			/*
			for (j = 0; j <= inz; ++j) {
				if (fabs(x4) > 0.0) {
					// Относительная погрешность менее 0.15%.
					if (fabs(100 * (zpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps) {
						block_indexes[i_1].kR = j;
						break;
					}
				}
				else {
					// Абсолютная погрешность.
					if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4,b,lb, w, lw, s, ls)) {
						block_indexes[i_1].kR = j;
						break;
					}
				}
			}
			*/
			bfound = false;

			if (!bfound) {
				for (j = 0; j <= inz; ++j) {
					// Абсолютная погрешность.
					if (fabs(zpos[j] - x4) < min_zazor_z /*shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)*/) {
						block_indexes[i_1].kR = j;
						bfound = true;
						break;
					}
				}
			}
			tzazor = 1.0e36;
			if (!bfound) {
				// Ищем наиболее близко расположенную границу.
				for (j = 0; j <= inz; ++j) {
					if (fabs(x4 - zpos[j]) < tzazor) {
						tzazor = fabs(x4 - zpos[j]);
						block_indexes[i_1].kR = j;
						bfound = true;
					}
				}
			}
			if (!bfound) {

				//std::cout << "Error incomming zE\n";
				//std::cout << "position=" << x4 << "\n";
				//std::cout << b[i].name << std::endl;
				//system("pause");

				for (j = inz; j >= 0; j--) {
					if ((!bfound) && (x4 > zpos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].kR = j;
						bfound = true;
						break;
					}
				}
			}



			if ((block_indexes[i_1].iL >= 0) &&
				(block_indexes[i_1].iR >= 0) &&
				(block_indexes[i_1].iL >= block_indexes[i_1].iR)) {
				printf("init_evt_f_alice_improved_obobshenie  function\n");
				printf("violation of the order block_indexes\n");
				printf("i=%d iL=%lld iR=%lld\n", i, block_indexes[i_1].iL,
					block_indexes[i_1].iR);
				system("pause");
				for (int i92 = 0; i92 <= inx; i92++) {
					std::cout << i92 << " " << xpos[i92] << " \n";
					system("pause");
				}
			}

			if ((block_indexes[i_1].iL >= 0) &&
				(block_indexes[i_1].iR >= 0) &&
				(block_indexes[i_1].iL + 1 == block_indexes[i_1].iR)) {
				// Всего одна клетка на блок.
				//iX_one_CELL_count_statistic++;
			}

			if ((block_indexes[i_1].jL >= 0) &&
				(block_indexes[i_1].jR >= 0) &&
				(block_indexes[i_1].jL >= block_indexes[i_1].jR)) {
				printf("init_evt_f_alice_improved_obobshenie  function\n");
				printf("violation of the order block_indexes\n");
				printf("i=%d jL=%lld jR=%lld\n", i, block_indexes[i_1].jL,
					block_indexes[i_1].jR);
				system("pause");
				for (int i92 = 0; i92 <= iny; i92++) {
					std::cout << i92 << " " << ypos[i92] << " \n";
					system("pause");
				}
			}


			if ((block_indexes[i_1].jL >= 0) &&
				(block_indexes[i_1].jR >= 0) &&
				(block_indexes[i_1].jL + 1 == block_indexes[i_1].jR)) {
				// Всего одна клетка на блок.
				//iY_one_CELL_count_statistic++;
			}

			if ((block_indexes[i_1].kL >= 0) &&
				(block_indexes[i_1].kR >= 0) &&
				(block_indexes[i_1].kL >= block_indexes[i_1].kR)) {
				printf("init_evt_f_alice_improved_obobshenie  function\n");
				printf("violation of the order block_indexes\n");
				printf("i=%d kL=%lld kR=%lld\n", i, block_indexes[i_1].kL,
					block_indexes[i_1].kR);
				system("pause");
				for (int i92 = 0; i92 <= inz; i92++) {
					std::cout << i92 << " " << zpos[i92] << " \n";
					system("pause");
				}
			}


			if ((block_indexes[i_1].kL >= 0) &&
				(block_indexes[i_1].kR >= 0) &&
				(block_indexes[i_1].kL + 1 == block_indexes[i_1].kR)) {
				// Всего одна клетка на блок.
				//iZ_one_CELL_count_statistic++;
			}

			if ((block_indexes[i_1].iL < 0) ||
				(block_indexes[i_1].iR < 0) ||
				(block_indexes[i_1].jL < 0) ||
				(block_indexes[i_1].jR < 0) ||
				(block_indexes[i_1].kL < 0) ||
				(block_indexes[i_1].kR < 0)) {
				if (b[i].g.itypegeom == PRISM) {
					printf("print for PRISM object...\n");
					printf("init_evt_f_alice_improved_obobshenie function\n");
					printf("i=%d iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
						block_indexes[i_1].iL,
						block_indexes[i_1].iR, block_indexes[i_1].jL,
						block_indexes[i_1].jR, block_indexes[i_1].kL,
						block_indexes[i_1].kR);
					printf("xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[i].g.xS, b[i].g.xE,
						b[i].g.yS, b[i].g.yE, b[i].g.zS, b[i].g.zE);
					printf("cabinet: xS=%e xE=%e yS=%e yE=%e zS=%e zE=%e\n", b[0].g.xS, b[0].g.xE,
						b[0].g.yS, b[0].g.yE, b[0].g.zS, b[0].g.zE);
					printf("ERROR: may be your geometry out of cabinet...\n");
					system("pause");
				}
			}

			i_1--;
		}

	}

	// Количество проходов существенно сократилось и в итоге это приводит к существенному
	// увеличению быстродействия.

	// Обязательная проверка !!!
	// И ногда предыдущий метод не срабатывает и это в случае
	// отсутствия исправления приводит к сбою.
	// Здесь приведена коррекция она медленней но работает в 100% случаев.
	// 28.07.2019
	for (integer i_a = lb - 1; i_a >= 0; i_a--) {
		if ((block_indexes[i_a].iL <= -1) || (block_indexes[i_a].iR <= -1) ||
			(block_indexes[i_a].jL <= -1) || (block_indexes[i_a].jR <= -1) ||
			(block_indexes[i_a].kL <= -1) || (block_indexes[i_a].kR <= -1)) {
			// Проблема признана, теперь работаем с каждым проявлением индивидуально.
			if (block_indexes[i_a].iL == -1) {
				doublereal x4 = b[i_a].g.xS;
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY_PLANE) || (b[i_a].g.iPlane == XZ_PLANE))) {
					x4 = b[i_a].g.xC - b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == YZ_PLANE))) {
					if (b[i_a].g.Hcyl > 0.0) {
						x4 = b[i_a].g.xC;
					}
					else {
						x4 = b[i_a].g.xC + b[i_a].g.Hcyl;
					}
				}
				integer j_amin = -1;
				doublereal t_min = 1.0e30;
				for (integer j_a = 0; j_a <= inx; j_a++) {
					// Абсолютная погрешность.
					if (fabs(xpos[j_a] - x4) < t_min) {
						t_min = fabs(xpos[j_a] - x4);
						j_amin = j_a;							
					}					
				}
				if (j_amin > -1) {
					block_indexes[i_a].iL = j_amin;
				}
			}
			if (block_indexes[i_a].iR == -1) {
				doublereal x4 = b[i_a].g.xE;
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY_PLANE) || (b[i_a].g.iPlane == XZ_PLANE))) {
					x4 = b[i_a].g.xC + b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == YZ_PLANE))) {
					if (b[i_a].g.Hcyl > 0.0) {
						x4 = b[i_a].g.xC + b[i_a].g.Hcyl;
					}
					else {
						x4 = b[i_a].g.xC;
					}
				}
				integer j_amin = -1;
				doublereal t_min = 1.0e30;
				for (integer j_a = 0; j_a <= inx; j_a++) {					
					// Абсолютная погрешность.
					if (fabs(xpos[j_a] - x4) < t_min) {
						t_min = fabs(xpos[j_a] - x4);
						j_amin = j_a;
					}
				}
				if (j_amin > -1) {
					block_indexes[i_a].iR = j_amin;
				}

			}
			if (block_indexes[i_a].jL == -1) {
				doublereal x4 = b[i_a].g.yS;
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY_PLANE) || (b[i_a].g.iPlane == YZ_PLANE))) {
					x4 = b[i_a].g.yC - b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XZ_PLANE))) {
					if (b[i_a].g.Hcyl > 0.0) {
						x4 = b[i_a].g.yC;
					}
					else {
						x4 = b[i_a].g.yC + b[i_a].g.Hcyl;
					}
				}
				integer j_amin = -1;
				doublereal t_min = 1.0e30;
				for (integer j_a = 0; j_a <= iny; j_a++) {
					// Абсолютная погрешность.
					if (fabs(ypos[j_a] - x4) < t_min) {
						t_min = fabs(ypos[j_a] - x4);
						j_amin = j_a;
					}
				}
				if (j_amin > -1) {
					block_indexes[i_a].jL = j_amin;
				}
			}
			if (block_indexes[i_a].jR == -1) {
				doublereal x4 = b[i_a].g.yE;
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY_PLANE) || (b[i_a].g.iPlane == YZ_PLANE))) {
					x4 = b[i_a].g.yC + b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XZ_PLANE))) {
					if (b[i_a].g.Hcyl > 0.0) {
						x4 = b[i_a].g.yC + b[i_a].g.Hcyl;
					}
					else {
						x4 = b[i_a].g.yC;
					}
				}
				integer j_amin = -1;
				doublereal t_min = 1.0e30;
				for (integer j_a = 0; j_a <= iny; j_a++) {
					// Абсолютная погрешность.
					if (fabs(ypos[j_a] - x4) < t_min) {
						t_min = fabs(ypos[j_a] - x4);
						j_amin = j_a;
					}
				}
				if (j_amin > -1) {
					block_indexes[i_a].jR = j_amin;
				}

			}
			if (block_indexes[i_a].kL == -1) {
				doublereal x4 = b[i_a].g.zS;
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == YZ_PLANE) || (b[i_a].g.iPlane == XZ_PLANE))) {
					x4 = b[i_a].g.zC - b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY_PLANE))) {
					if (b[i_a].g.Hcyl > 0.0) {
						x4 = b[i_a].g.zC;
					}
					else {
						x4 = b[i_a].g.zC + b[i_a].g.Hcyl;
					}
				}
				integer j_amin = -1;
				doublereal t_min = 1.0e30;
				for (integer j_a = 0; j_a <= inz; j_a++) {
					// Абсолютная погрешность.
					if (fabs(zpos[j_a] - x4) < t_min) {
						t_min = fabs(zpos[j_a] - x4);
						j_amin = j_a;
					}
				}
				if (j_amin > -1) {
					block_indexes[i_a].kL = j_amin;
				}
			}
			if (block_indexes[i_a].kR == -1) {
				doublereal x4 = b[i_a].g.zE;
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == YZ_PLANE) || (b[i_a].g.iPlane == XZ_PLANE))) {
					x4 = b[i_a].g.zC + b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY_PLANE))) {
					if (b[i_a].g.Hcyl > 0.0) {
						x4 = b[i_a].g.zC + b[i_a].g.Hcyl;
					}
					else {
						x4 = b[i_a].g.zC;
					}
				}
				integer j_amin = -1;
				doublereal t_min = 1.0e30;
				for (integer j_a = 0; j_a <= inz; j_a++) {
					// Абсолютная погрешность.
					if (fabs(zpos[j_a] - x4) < t_min) {
						t_min = fabs(zpos[j_a] - x4);
						j_amin = j_a;
					}
				}
				if (j_amin > -1) {
					block_indexes[i_a].kR = j_amin;
				}

			}
		}
	}
	
	const integer size_xyz=inx*iny*inz;
#pragma omp parallel for
	for (integer iP = 0; iP < size_xyz; ++iP) {
		evt[iP] = -1;
		bvisit[iP] = false;
	}
	 integer m7 = lb-1;
	// integer iP_k1, iP_j1, iP, m8;

	 if ((b_on_adaptive_local_refinement_mesh)&&(hash_for_droblenie_xyz!=nullptr)&&
		 (itype_ALICE_Mesh == TYPE_ALICE_MESH::ONE_PASS_COARSE_ALICE_MESH)) {


		 for (integer k1 = 0; k1 < inz; k1++)
		 {
			 integer iP_k1 = k1 * inx * iny;
			 for (integer j1 = 0; j1 < iny; j1++)
			 {
				 integer iP_j1 = j1 * inx + iP_k1;
				 for (integer i1 = 0; i1 < inx; i1++)
				 {
					 integer iP = i1 + iP_j1;

					 evt[iP] = hash_for_droblenie_xyz[i1][j1][k1];
				 }
			 }
		 }


		// Освобождение оперативной памяти из под хеш-таблицы.
		for (int i_54 = 0; i_54 < inx; i_54++) {
			 for (int i_55 = 0; i_55 < iny; i_55++) {
				 delete[] hash_for_droblenie_xyz[i_54][i_55];
				 hash_for_droblenie_xyz[i_54][i_55] = nullptr;
			 }
		 }

		 for (int i_54 = 0; i_54 < inx; i_54++) {
			 delete[] hash_for_droblenie_xyz[i_54];
			 hash_for_droblenie_xyz[i_54] = nullptr;
		 }

		 delete[] hash_for_droblenie_xyz;
		 hash_for_droblenie_xyz = nullptr;

	 }
	 else {

		 for (int m8 = lb - 1; m8 >= 0; --m8) {
			 if (b[m8].g.itypegeom == PRISM) {
#pragma omp parallel for

				 for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; ++k1)
				 {
					 integer iP_k1 = k1 * inx * iny;
					 for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; ++j1)
					 {
						 integer iP_j1 = j1 * inx + iP_k1;
						 for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; ++i1)
						 {
							 integer iP = i1 + iP_j1;

							 if (bvisit[iP] == false)
							 {

								 bvisit[iP] = true;

								 switch (iflag) {
								 case TEMPERATURE:
									 if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
										 evt[iP] = -1;
									 }
									 else {
										 evt[iP] = m8;
									 }
									 break;
								 case HYDRODINAMIC:
									 if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
										 evt[iP] = -1;
									 }
									 else {
										 evt[iP] = m8;
									 }
									 break;
								 }
							 }
						 }

					 }
				 }
				 m7--;
			 }
			 else if (b[m8].g.itypegeom == CYLINDER) {

				 // как был сформирован призматический объект для цилиндра ? 
				 // Надо также сократить число проверяемых точек.
				 // Cylinder
				 //for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {

				 for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; ++k1)
				 {
					 integer iP_k1 = k1 * inx * iny;
					 for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; ++j1)
					 {
						 integer  iP_j1 = j1 * inx + iP_k1;
						 for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; ++i1)
						 {
							 integer  iP = i1 + iP_j1;

							 if (bvisit[iP] == false)
							 {

								 TOCHKA p;
								 p.x = 0.5 * (xpos[i1] + xpos[i1 + 1]);
								 p.y = 0.5 * (ypos[j1] + ypos[j1 + 1]);
								 p.z = 0.5 * (zpos[k1] + zpos[k1 + 1]);

								 switch (b[m8].g.iPlane) {
								 case XY_PLANE:
									 if (fabs(b[m8].g.R_in_cyl) < 1.0e-36) {
										 if ((p.z > b[m8].g.zC) && (p.z < b[m8].g.zC + b[m8].g.Hcyl)) {
											 if (sqrt((b[m8].g.xC - p.x) * (b[m8].g.xC - p.x) + (b[m8].g.yC - p.y) * (b[m8].g.yC - p.y)) < b[m8].g.R_out_cyl) {

												 bvisit[iP] = true;

												 switch (iflag) {
												 case TEMPERATURE:
													 if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
														 evt[iP] = -1;
													 }
													 else {
														 evt[iP] = m8;
													 }
													 break;
												 case HYDRODINAMIC:
													 if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
														 evt[iP] = -1;
													 }
													 else {
														 evt[iP] = m8;
													 }
													 break;
												 }

											 }
										 }
									 }
									 else {
										 if ((p.z > b[m8].g.zC) && (p.z < b[m8].g.zC + b[m8].g.Hcyl)) {
											 if (sqrt((b[m8].g.xC - p.x) * (b[m8].g.xC - p.x) + (b[m8].g.yC - p.y) * (b[m8].g.yC - p.y)) < b[m8].g.R_out_cyl) {
												 if (sqrt((b[m8].g.xC - p.x) * (b[m8].g.xC - p.x) + (b[m8].g.yC - p.y) * (b[m8].g.yC - p.y)) > b[m8].g.R_in_cyl) {


													 bvisit[iP] = true;

													 switch (iflag) {
													 case TEMPERATURE:
														 if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
															 evt[iP] = -1;
														 }
														 else {
															 evt[iP] = m8;
														 }
														 break;
													 case HYDRODINAMIC:
														 if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
															 evt[iP] = -1;
														 }
														 else {
															 evt[iP] = m8;
														 }
														 break;
													 }
												 }
											 }
										 }
									 }
									 break;
								 case XZ_PLANE:
									 if (fabs(b[m8].g.R_in_cyl) < 1.0e-36) {
										 if ((p.y > b[m8].g.yC) && (p.y < b[m8].g.yC + b[m8].g.Hcyl)) {
											 if (sqrt((b[m8].g.xC - p.x) * (b[m8].g.xC - p.x) + (b[m8].g.zC - p.z) * (b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {

												 bvisit[iP] = true;

												 switch (iflag) {
												 case TEMPERATURE:
													 if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
														 evt[iP] = -1;
													 }
													 else {
														 evt[iP] = m8;
													 }
													 break;
												 case HYDRODINAMIC:
													 if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
														 evt[iP] = -1;
													 }
													 else {
														 evt[iP] = m8;
													 }
													 break;
												 }
											 }
										 }
									 }
									 else {
										 if ((p.y > b[m8].g.yC) && (p.y < b[m8].g.yC + b[m8].g.Hcyl)) {
											 if (sqrt((b[m8].g.xC - p.x) * (b[m8].g.xC - p.x) + (b[m8].g.zC - p.z) * (b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {
												 if (sqrt((b[m8].g.xC - p.x) * (b[m8].g.xC - p.x) + (b[m8].g.zC - p.z) * (b[m8].g.zC - p.z)) > b[m8].g.R_in_cyl) {


													 bvisit[iP] = true;

													 switch (iflag) {
													 case TEMPERATURE:
														 if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
															 evt[iP] = -1;
														 }
														 else {
															 evt[iP] = m8;
														 }
														 break;
													 case HYDRODINAMIC:
														 if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
															 evt[iP] = -1;
														 }
														 else {
															 evt[iP] = m8;
														 }
														 break;
													 }
												 }
											 }
										 }
									 }
									 break;
								 case YZ_PLANE:
									 if (fabs(b[m8].g.R_in_cyl) < 1.0e-36) {
										 if ((p.x > b[m8].g.xC) && (p.x < b[m8].g.xC + b[m8].g.Hcyl)) {
											 if (sqrt((b[m8].g.yC - p.y) * (b[m8].g.yC - p.y) + (b[m8].g.zC - p.z) * (b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {



												 bvisit[iP] = true;

												 switch (iflag) {
												 case TEMPERATURE:
													 if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
														 evt[iP] = -1;
													 }
													 else {
														 evt[iP] = m8;
													 }
													 break;
												 case HYDRODINAMIC:
													 if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
														 evt[iP] = -1;
													 }
													 else {
														 evt[iP] = m8;
													 }
													 break;
												 }

											 }
										 }
									 }
									 else {
										 if ((p.x > b[m8].g.xC) && (p.x < b[m8].g.xC + b[m8].g.Hcyl)) {
											 if (sqrt((b[m8].g.yC - p.y) * (b[m8].g.yC - p.y) + (b[m8].g.zC - p.z) * (b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {
												 if (sqrt((b[m8].g.yC - p.y) * (b[m8].g.yC - p.y) + (b[m8].g.zC - p.z) * (b[m8].g.zC - p.z)) > b[m8].g.R_in_cyl) {



													 bvisit[iP] = true;

													 switch (iflag) {
													 case TEMPERATURE:
														 if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
															 evt[iP] = -1;
														 }
														 else {
															 evt[iP] = m8;
														 }
														 break;
													 case HYDRODINAMIC:
														 if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
															 evt[iP] = -1;
														 }
														 else {
															 evt[iP] = m8;
														 }
														 break;
													 }
												 }
											 }
										 }
									 }
									 break;
								 }
							 }
						 }
					 }
				 }
				 m7--;
			 }
			 else if (b[m8].g.itypegeom == POLYGON) {

				 // polygon
				 // Мы сокращаем число проверяемых точек 
				 // рассматривая только точки внутри окаймляющей прямоугольной призмы.



#pragma omp parallel for
				 for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; ++k1)
				 {
					 integer  iP_k1 = k1 * inx * iny;
					 for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; ++j1)
					 {
						 integer  iP_j1 = j1 * inx + iP_k1;
						 for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; ++i1)
						 {

							 integer  iP = i1 + iP_j1;

							 if (bvisit[iP] == false)
							 {

								 //for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
								 TOCHKA p;
								 p.x = 0.5 * (xpos[i1] + xpos[i1 + 1]);
								 p.y = 0.5 * (ypos[j1] + ypos[j1 + 1]);
								 p.z = 0.5 * (zpos[k1] + zpos[k1 + 1]);

								 int k74 = -1;
								 if (in_polygon(p, b[m8].g.nsizei, b[m8].g.xi, b[m8].g.yi, b[m8].g.zi, b[m8].g.hi, b[m8].g.iPlane_obj2, k74, m8)) {
									 //printf("i1=%d j1=%d k1=%d inx*iny*inz=%d\n",i1,j1,k1, inx*iny*inz);
									 //printf("iL=%d iR=%d jL=%d jR=%d kL=%d kR=%d\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);

									 bvisit[iP] = true;

									 switch (iflag) {
									 case TEMPERATURE:
										 if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
											 evt[iP] = -1;
										 }
										 else {
											 evt[iP] = m8;
										 }
										 break;
									 case HYDRODINAMIC:
										 if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
											 evt[iP] = -1;
										 }
										 else {
											 evt[iP] = m8;
										 }
										 break;
									 }
								 }
							 }
						 }
					 }

				 }
				 m7--;
			 }
			 else if (b[m8].g.itypegeom == CAD_STL) {

				 // CAD_STL
				 // Мы сокращаем число проверяемых точек 
				 // рассматривая только точки внутри окаймляющей прямоугольной призмы.



#pragma omp parallel for
				 for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; ++k1)
				 {
					 integer  iP_k1 = k1 * inx * iny;
					 for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; ++j1)
					 {
						 integer  iP_j1 = j1 * inx + iP_k1;
						 for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; ++i1)
						 {

							 integer  iP = i1 + iP_j1;

							 if (bvisit[iP] == false)
							 {

								 //for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
								 TOCHKA p;
								 p.x = 0.5 * (xpos[i1] + xpos[i1 + 1]);
								 p.y = 0.5 * (ypos[j1] + ypos[j1 + 1]);
								 p.z = 0.5 * (zpos[k1] + zpos[k1 + 1]);

								 int k74 = -1;

								 if (b[m8].g.in_CAD_STL_check(p, k74, m8))
								 {
									 //printf("i1=%d j1=%d k1=%d inx*iny*inz=%d\n",i1,j1,k1, inx*iny*inz);
									 //printf("iL=%d iR=%d jL=%d jR=%d kL=%d kR=%d\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);

									 bvisit[iP] = true;

									 switch (iflag) {
									 case TEMPERATURE:
										 if (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
											 evt[iP] = -1;
										 }
										 else {
											 evt[iP] = m8;
										 }
										 break;
									 case HYDRODINAMIC:
										 if ((b[m8].itype == PHYSICS_TYPE_IN_BODY::SOLID) || (b[m8].itype == PHYSICS_TYPE_IN_BODY::HOLLOW)) {
											 evt[iP] = -1;
										 }
										 else {
											 evt[iP] = m8;
										 }
										 break;
									 }
								 }
							 }
						 }
					 }

				 }
				 m7--;
			 }
		 }

	 }

	//if (bvisit != nullptr) {
	// Оператор delete может быть вызван повторно и даже к null указателю.
		delete[] bvisit;
		bvisit = nullptr;
	//}

		if (print_log_message) {
			printf("enumerate_volume_improved 80 procent.\n");
		}
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	//  integer iP_k, iP_j;
	 int l = 1;
	for (k = 0; k < inz; ++k)
	{
		integer  iP_k = k * inx * iny;
		for (j = 0; j < iny; ++j)
		{
			integer  iP_j = j * inx + iP_k;
			for (i = 0; i < inx; ++i)
			{
				integer  iP = i + iP_j;
				if (evt[iP] > -1) {
					//integer ib = evt[iP]; // номер блока был сохранён ранее.
					// Это очень нужно для записи репорта.
					//whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
					
					evt[iP] = l;
					tck_int_list[l - 1].i = i;
					tck_int_list[l - 1].j = j;
					tck_int_list[l - 1].k = k;
					++l;
				}
				else {
					// не принадлежит расчётной области
					evt[iP] = 0;
				}
			}
		}
	}


	delete[] block_indexes;

	if (print_log_message) {
		printf("enumerate_volume_improved end.\n");
	}
	/*
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	integer l = 1, ib;
	integer i, j, k;
	bool inDomain = false;
	// нумерация контрольных объёмов начинается с единицы
	// если контрольный объём не принадлежит расчётной области то ставится 0.
	for (i = 0; i<inx; ++i) for (j = 0; j<iny; ++j) for (k = 0; k<inz; ++k) {
	p.x = 0.5*(xpos[i] + xpos[i + 1]);
	p.y = 0.5*(ypos[j] + ypos[j + 1]);
	p.z = 0.5*(zpos[k] + zpos[k + 1]);
	switch (iflag) {
	case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
	case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
	}
	if (inDomain) {
	// принадлежит расчётной области
	evt[i + j*inx + k*inx*iny] = l;
	l++;
	}
	else
	{   // не принадлежит расчётной области
	evt[i + j*inx + k*inx*iny] = 0;
	}
	}
	*/

} // init_evt_f_alice_improved_obobshenie.

// Инициализация evt_f. Он нужен для заливки, через которую определяется функция цвета.
// А функция цвета нужна обязательным образом для корректировки массового баланса при нескольких
// FLUID областях связанных лишь уравнением теплопередачи, иначе не будет сходимости.
void init_evt_f_alice(int* &evt,  integer iflag, doublereal* xpos, doublereal* ypos, doublereal* zpos, 
	int inx, int iny, int inz,
	BLOCK* b, int lb, TOCKA_SHORT_INT* &tck_int_list,
	WALL* &w, int &lw, SOURCE* &s, int &ls) {


	bool b_only_Prism_object = true;
	for (integer m7 = 0; m7 < lb; m7++) {
		if (b[m7].g.itypegeom != PRISM) {
			b_only_Prism_object = false;
			break;
		}
	}

	//if (1) {
	if (!CAD_GEOMETRY_OCTREE_MESHGEN) {
		// 2.01.2018
		// Обобщение, работающее не только на прямоугольных призмах.
		// Работает на цилиндрах, прямоугольных призмах и полигонах.
		init_evt_f_alice_improved_obobshenie(evt, iflag, xpos, ypos, zpos, inx, iny, inz, b, lb, tck_int_list,w,lw,s,ls);
	}
	else {
		if (b_only_Prism_object) {
			// Модель полностью состоит из прямоугольных призм.

			init_evt_f_alice_improved(evt, iflag, xpos, ypos, zpos, inx, iny, inz, b, lb, tck_int_list, w,lw,s,ls);
		}
		else {
			// В модели присутствуют не только прямоугольные призмы.

			tck_int_list = new TOCKA_SHORT_INT[inx*iny*inz];
			if (tck_int_list == nullptr) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment for evt tck_int_list...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}

			evt = nullptr;
			evt = new int[inx*iny*inz];
			if (evt == nullptr) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment for evt constr struct...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}

			TOCHKA p;
			// нумерация в evt начиная с единицы.
			// если не принадлежит расчётной области то стоит 0.
			int l = 1;
			int ib;
			// integer i, j, k;
			// integer iP_k, iP_j, iP;
			bool inDomain = false;
			// нумерация контрольных объёмов начинается с единицы
			// если контрольный объём не принадлежит расчётной области то ставится 0.
			for (int  k = 0; k < inz; ++k)
			{
				integer  iP_k=k*inx*iny;
				for (int  j = 0; j < iny; ++j) 
				{
					 integer  iP_j=j*inx +iP_k;
					 for (int  i = 0; i < inx; ++i)
					{
						integer  iP=i +iP_j;
						p.x = 0.5*(xpos[i] + xpos[i + 1]);
						p.y = 0.5*(ypos[j] + ypos[j + 1]);
						p.z = 0.5*(zpos[k] + zpos[k + 1]);
						switch (iflag) {
							case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
							case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
						}
						if (inDomain) {
							// принадлежит расчётной области
							evt[iP] = l;
							tck_int_list[l - 1].i = i;
							tck_int_list[l - 1].j = j;
							tck_int_list[l - 1].k = k;
							l++;
						}
						else
						{   // не принадлежит расчётной области
							evt[iP] = 0;
						}
					}
				}
			}
		}
	}

} // init_evt_f_alice.



// находит соседей для каждого внутреннего контрольного
// объёма или 0 если соседа нет.
void constr_neighbour(int* evt, int* ent, int** &neighbour, integer maxelm,
				  int inx, int iny, int inz, TOCKA_SHORT_INT* &tck_int_list) {
    integer i;
	neighbour=nullptr;
	neighbour = new int*[12];
	if (neighbour==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for neighbour constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (i=0; i<12; ++i) {
		neighbour[i]=nullptr;
	}
	for (i=0; i<12; ++i) {
		neighbour[i]=new int[maxelm];
		if (neighbour[i]==nullptr) {
	       // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem: not enough memory on your equipment for neighbour[%lld] constr struct...\n", i);
#else
			printf("Problem: not enough memory on your equipment for neighbour[%d] constr struct...\n", i);
#endif
		   
		   printf("Please any key to exit...\n");
		   //system("PAUSE");
		   system("pause");
		   exit(1);
	    }
	}

	integer ic,l=0;
	integer iE, iN, iT, iW, iS, iB;
    integer iEE, iNN, iTT, iWW, iSS, iBB;
	// проход по всем контрольным объёмам принадлежащим
	// расчётной области.
   // for (i=0; i<(inx); ++i) for (j=0; j<(iny); ++j) for (k=0; k<(inz); ++k) {
	//	ic=i+j*inx+k*inx*iny;
		//if (evt[ic]>0) {

	// Константа для ускорения проверок. 
	const integer iend_size = inx*iny*inz;
	const integer iend_size1 = inx*iny;

	for (integer iscan = 0; iscan < maxelm; iscan++) {
		{
			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;
			ic = i + j * inx + k * inx*iny;
			


				if (l < maxelm) {

					//iE = (i + 1) + j*inx + k*inx*iny + 1;
					iE = ic + 2;
					if (((i + 1) < inx) && (iE > 0) && (iE <= iend_size) && (evt[iE - 1] > 0)) neighbour[E_SIDE][l] = evt[iE - 1];
					else {
						neighbour[E_SIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iW = (i - 1) + j*inx + k*inx*iny + 1;
					iW = ic;
					if (((i - 1) >= 0) && (iW > 0) && (iW <= iend_size) && (evt[iW - 1] > 0)) neighbour[W_SIDE][l] = evt[iW - 1];
					else {
						neighbour[W_SIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iN = i + (j + 1)*inx + k*inx*iny + 1;
					iN = ic + 1 + inx;
					if (((j + 1) < iny) && (iN > 0) && (iN <= iend_size) && (evt[iN - 1] > 0)) neighbour[N_SIDE][l] = evt[iN - 1];
					else {
						neighbour[N_SIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iS = i + (j - 1)*inx + k*inx*iny + 1;
					iS = ic - inx + 1;
					if (((j - 1) >= 0) && (iS > 0) && (iS <= iend_size) && (evt[iS - 1] > 0)) neighbour[S_SIDE][l] = evt[iS - 1];
					else {
						neighbour[S_SIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iT = i + j*inx + (k + 1)*inx*iny + 1;
					iT = ic + iend_size1 + 1;
					if (((k + 1) < inz) && (iT > 0) && (iT <= iend_size) && (evt[iT - 1] > 0)) neighbour[T_SIDE][l] = evt[iT - 1];
					else {
						neighbour[T_SIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iB = i + j*inx + (k - 1)*inx*iny + 1;
					iB = ic - iend_size1 + 1;
					if (((k - 1) >= 0) && (iB > 0) && (iB <= iend_size) && (evt[iB - 1] > 0)) neighbour[B_SIDE][l] = evt[iB - 1];
					else {
						neighbour[B_SIDE][l] = 0; // точка не принадлежит расчётной области
					}

					//iEE = (i + 2) + j*inx + k*inx*iny + 1;
					iEE = ic + 3;
					if (((i + 2) < inx) && (iEE > 0) && (iEE <= iend_size) && (evt[iEE - 1] > 0)) neighbour[EE_SIDE][l] = evt[iEE - 1];
					else neighbour[EE_SIDE][l] = 0; // точка не принадлежит расчётной области
					//iWW = (i - 2) + j*inx + k*inx*iny + 1;
					iWW = ic - 1;
					if (((i - 2) >= 0) && (iWW > 0) && (iWW <= iend_size) && (evt[iWW - 1] > 0)) neighbour[WW_SIDE][l] = evt[iWW - 1];
					else neighbour[WW_SIDE][l] = 0; // точка не принадлежит расчётной области
					//iNN = i + (j + 2)*inx + k*inx*iny + 1;
					iNN = ic + 2 * inx + 1;
					if (((j + 2) < iny) && (iNN > 0) && (iNN <= iend_size) && (evt[iNN - 1] > 0)) neighbour[NN_SIDE][l] = evt[iNN - 1];
					else neighbour[NN_SIDE][l] = 0; // точка не принадлежит расчётной области
					//iSS = i + (j - 2)*inx + k*inx*iny + 1;
					iSS = ic - 2 * inx + 1;
					if (((j - 2) >= 0) && (iSS > 0) && (iSS <= iend_size) && (evt[iSS - 1] > 0)) neighbour[SS_SIDE][l] = evt[iSS - 1];
					else neighbour[SS_SIDE][l] = 0; // точка не принадлежит расчётной области
					//iTT = i + j*inx + (k + 2)*inx*iny + 1;
					iTT = ic + 2 * iend_size1 + 1;
					if (((k + 2) < inz) && (iTT > 0) && (iTT <= iend_size) && (evt[iTT - 1] > 0)) neighbour[TT_SIDE][l] = evt[iTT - 1];
					else neighbour[TT_SIDE][l] = 0; // точка не принадлежит расчётной области
					//iBB = i + j*inx + (k - 2)*inx*iny + 1;
					iBB = ic - 2 * iend_size1 + 1;
					if (((k - 2) >= 0) && (iBB > 0) && (iBB <= iend_size) && (evt[iBB - 1] > 0)) neighbour[BB_SIDE][l] = evt[iBB - 1];
					else neighbour[BB_SIDE][l] = 0; // точка не принадлежит расчётной области

					l++;
				}
				else {
					printf("model error: l>maxelm\n");
					system("pause");
					exit(1);
				}
		}
	}

} // constr_neighbour

// находит соседей для каждого внутреннего контрольного
// объёма или 0 если соседа нет.
void constr_neighbour_flow(int **evt_f2, integer iDom, 
					   int** &neighbour, integer maxelm,
				       int inx, int iny, int inz) {
    //integer i,j,k;	
	neighbour=nullptr;
	neighbour = new int*[12];
	if (neighbour==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for neighbour constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer i=0; i<12; ++i) neighbour[i]=nullptr;
	for (integer i=0; i<12; ++i) {
		neighbour[i]=new int[maxelm];
		if (neighbour[i]==nullptr) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem: not enough memory on your equipment for neighbour[%lld] constr struct...\n", i);
#else
			printf("Problem: not enough memory on your equipment for neighbour[%d] constr struct...\n", i);
#endif
		    
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}

	// Константа для ускорения проверок.
	const integer iend_size = inx*iny*inz;

	 integer l=0;
	//integer ic;
	// integer iE, iN, iT, iW, iS, iB;
    // integer iEE, iNN, iTT, iWW, iSS, iBB;
	// integer ic_k, ic_k_m1, ic_k_p1, ic_k_m2, ic_k_p2;
	// integer ic_j, ic_j_p1, ic_j_m1, ic_j_p2, ic_j_m2;
	// проход по всем контрольным объёмам принадлежащим
	// расчётной области.
	for (integer  k=0; k<(inz); ++k)
	{
		integer ic_k=k*inx*iny;
		integer ic_k_m1=(k-1)*inx*iny;
        integer ic_k_p1=(k+1)*inx*iny;
		integer ic_k_m2=(k-2)*inx*iny;
        integer ic_k_p2=(k+2)*inx*iny;
		for (integer  j=0; j<(iny); ++j)
		{
			 integer ic_j=j*inx+ic_k;
			 integer ic_j_p1=ic_k_p1+j*inx+1;
			 integer ic_j_m1=ic_k_m1+j*inx+1;
			 integer ic_j_p2=ic_k_p2+j*inx+1;
			 integer ic_j_m2=ic_k_m2+j*inx+1;
			 for (integer  i=0; i<(inx); ++i)  
			 {

				integer ic=i+ic_j;
				if ((evt_f2[ENUMERATECONTVOL][ic]>0) && (evt_f2[MASKDOMAINFLUID][ic] == (iDom + 1))) {
			
					if (l < maxelm) {
							 integer iE, iN, iT, iW, iS, iB;
							 integer iEE, iNN, iTT, iWW, iSS, iBB;

							//iE = (i + 1) + j*inx + k*inx*iny + 1;
							iE = ic + 2;
							if (((i + 1) < inx) && (iE > 0) && (iE <= iend_size) && (evt_f2[0][iE - 1] > 0) && (evt_f2[1][iE - 1] == (iDom + 1))) neighbour[E_SIDE][l] = evt_f2[0][iE - 1];
							else neighbour[E_SIDE][l] = 0; // точка не принадлежит расчётной области
							//iW = (i - 1) + j*inx + k*inx*iny + 1;
							iW = ic;
							if (((i - 1) >= 0) && (iW > 0) && (iW <= iend_size) && (evt_f2[0][iW - 1] > 0) && (evt_f2[1][iW - 1] == (iDom + 1))) neighbour[W_SIDE][l] = evt_f2[0][iW - 1];
							else neighbour[W_SIDE][l] = 0; // точка не принадлежит расчётной области
							//iN = i + (j + 1)*inx + k*inx*iny + 1;
							iN = ic + inx + 1;
							if (((j + 1) < iny) && (iN > 0) && (iN <= iend_size) && (evt_f2[0][iN - 1] > 0) && (evt_f2[1][iN - 1] == (iDom + 1))) neighbour[N_SIDE][l] = evt_f2[0][iN - 1];
							else neighbour[N_SIDE][l] = 0; // точка не принадлежит расчётной области
							//iS = i + (j - 1)*inx + k*inx*iny + 1;
							iS = ic - inx + 1;
							if (((j - 1) >= 0) && (iS > 0) && (iS <= iend_size) && (evt_f2[0][iS - 1] > 0) && (evt_f2[1][iS - 1] == (iDom + 1))) neighbour[S_SIDE][l] = evt_f2[0][iS - 1];
							else neighbour[S_SIDE][l] = 0; // точка не принадлежит расчётной области
							iT = i + ic_j_p1;
							if (((k + 1) < inz) && (iT > 0) && (iT <= iend_size) && (evt_f2[0][iT - 1] > 0) && (evt_f2[1][iT - 1] == (iDom + 1))) neighbour[T_SIDE][l] = evt_f2[0][iT - 1];
							else neighbour[T_SIDE][l] = 0; // точка не принадлежит расчётной области
							iB = i + ic_j_m1;
							if (((k - 1) >= 0) && (iB > 0) && (iB <= iend_size) && (evt_f2[0][iB - 1] > 0) && (evt_f2[1][iB - 1] == (iDom + 1))) neighbour[B_SIDE][l] = evt_f2[0][iB - 1];
							else neighbour[B_SIDE][l] = 0; // точка не принадлежит расчётной области

							//iEE = (i + 2) + j*inx + k*inx*iny + 1;
							iEE = ic + 3;
							if (((i + 2) < inx) && (iEE > 0) && (iEE <= iend_size) && (evt_f2[0][iEE - 1] > 0) && (evt_f2[1][iEE - 1] == (iDom + 1))) neighbour[EE_SIDE][l] = evt_f2[0][iEE - 1];
							else neighbour[EE_SIDE][l] = 0; // точка не принадлежит расчётной области
							//iWW = (i - 2) + j*inx + k*inx*iny + 1;
							iWW = ic - 1;
							if (((i - 2) >= 0) && (iWW > 0) && (iWW <= iend_size) && (evt_f2[0][iWW - 1] > 0) && (evt_f2[1][iWW - 1] == (iDom + 1))) neighbour[WW_SIDE][l] = evt_f2[0][iWW - 1];
							else neighbour[WW_SIDE][l] = 0; // точка не принадлежит расчётной области
							//iNN = i + (j + 2)*inx + k*inx*iny + 1;
							iNN = ic + 2 * inx + 1;
							if (((j + 2) < iny) && (iNN > 0) && (iNN <= iend_size) && (evt_f2[0][iNN - 1] > 0) && (evt_f2[1][iNN - 1] == (iDom + 1))) neighbour[NN_SIDE][l] = evt_f2[0][iNN - 1];
							else neighbour[NN_SIDE][l] = 0; // точка не принадлежит расчётной области
							//iSS = i + (j - 2)*inx + k*inx*iny + 1;
							iSS = ic - 2 * inx + 1;
							if (((j - 2) >= 0) && (iSS > 0) && (iSS <= iend_size) && (evt_f2[0][iSS - 1] > 0) && (evt_f2[1][iSS - 1] == (iDom + 1))) neighbour[SS_SIDE][l] = evt_f2[0][iSS - 1];
							else neighbour[SS_SIDE][l] = 0; // точка не принадлежит расчётной области
							iTT = i + ic_j_p2;
							if (((k + 2) < inz) && (iTT > 0) && (iTT <= iend_size) && (evt_f2[0][iTT - 1] > 0) && (evt_f2[1][iTT - 1] == (iDom + 1))) neighbour[TT_SIDE][l] = evt_f2[0][iTT - 1];
							else neighbour[TT_SIDE][l] = 0; // точка не принадлежит расчётной области
							iBB = i + ic_j_m2;
							if (((k - 2) >= 0) && (iBB > 0) && (iBB <= iend_size) && (evt_f2[0][iBB - 1] > 0) && (evt_f2[1][iBB - 1] == (iDom + 1))) neighbour[BB_SIDE][l] = evt_f2[0][iBB - 1];
							else neighbour[BB_SIDE][l] = 0; // точка не принадлежит расчётной области

							l++;
						}
						else {
							printf("model is incorrect: constr struct neighbour index >= maxelm\n");
							system("pause");
							exit(1);
						}
					}
				}
			}
		}

} // constr_neighbour_flow

// Всем соседям присваивается идентификатор центрального узла 
void my_fill_Domain(int* &evt_f, int i, int j, int k, 
					int inx, int iny, int inz) {
        integer iP;
		
		iP=i+j*inx+k*inx*iny;
        //iE=(i+1)+j*inx+k*inx*iny+1;
		//iW=(i-1)+j*inx+k*inx*iny+1;
        //iN=i+(j+1)*inx+k*inx*iny+1;
		//iS=i+(j-1)*inx+k*inx*iny+1;
		//iT=i+j*inx+(k+1)*inx*iny+1;
		//iB=i+j*inx+(k-1)*inx*iny+1;

		//if (((i)<inx) && (iP>=0) && (iP < isize) && evt_f[iP]>0) {
		if (evt_f[iP]>0) {
			
			integer  iE, iW, iN, iS, iT, iB;
			integer isize = inx*iny*inz;
			iE = iP + 2;
			iW = iP;
			iN = iP + 1 + inx;
			iS = iP + 1 - inx;
			integer isize_loc = inx*iny;
			iT = iP + 1 + isize_loc;
			iB = iP + 1 - isize_loc;
            // Узел принадлежит зоне FLUID.
			// Всем существующим соседям присвоить этот номер. 
            if (((i+1)<inx) && (iE>0) && (iE<= isize) && (evt_f[iE-1]>0) && (evt_f[iE - 1]!= evt_f[iP])) evt_f[iE-1]=evt_f[iP];
            if (((i-1)>=0) && (iW>0) && (iW<= isize) && (evt_f[iW-1]>0) && (evt_f[iW - 1] != evt_f[iP])) evt_f[iW-1]=evt_f[iP];
            if (((j+1)<iny) && (iN>0) && (iN<= isize) && (evt_f[iN-1]>0) && (evt_f[iN - 1] != evt_f[iP])) evt_f[iN-1]=evt_f[iP];
            if (((j-1)>=0) && (iS>0) && (iS<= isize) && (evt_f[iS-1]>0) && (evt_f[iS - 1] != evt_f[iP])) evt_f[iS-1]=evt_f[iP];
            if (((k+1)<inz) && (iT>0) && (iT<= isize) && (evt_f[iT-1]>0) && (evt_f[iT - 1] != evt_f[iP])) evt_f[iT-1]=evt_f[iP];
            if (((k-1)>=0) && (iB>0) && (iB<= isize) && (evt_f[iB-1]>0) && (evt_f[iB - 1] != evt_f[iP])) evt_f[iB-1]=evt_f[iP];
		}
} // my_fill_Domain


// Всем соседям присваивается идентификатор центрального узла 
void my_fill_Domain_recursive(int* &evt_f, int i, int j, int k, 
					int inx, int iny, int inz, STEK* &stek) {
        
		integer iP, iE, iW, iN, iS, iT, iB;


		/* Рекурсивная процедура заливки распознаёт области произвольной сложности,
		а чтобы не было переполнения стека рекурсивных вызовов рекурсия полностью заменена
		на итерацию в цикле с помощью стека.
		*/

		
		
		
		// оператор new не требует проверки на null.
		//if (stek==nullptr) {
			//printf("no allocate memory for my_fill_Domain_recursive: stek==nullptr\n");
			//system("PAUSE");
			//system("pause");
			//exit(1);
		//}
		integer itop=0; // маркер вершины стека.

		stek[itop].i=i;
		stek[itop].j=j;
		stek[itop].k=k;
		itop++;

		// Большое количество рекурсивных вызовов не поддерживается компилятором,
		// надо заменить рекурсию на итерацию.
		

		while (itop>0) {

			   // Вынимаем вершину стека.
               i=stek[itop-1].i;
		       j=stek[itop-1].j;
		       k=stek[itop-1].k;
			   itop--;

			   const integer size_xyz=inx*iny*inz;
               const integer iP_k=k*inx*iny;
			   const integer iP_j=iP_k+1+j*inx;
			   const integer iP_k_p1=(k+1)*inx*iny;
			   const integer iP_k_m1=(k-1)*inx*iny;

	           iP=i+j*inx+iP_k;
               iE=(i+1)+iP_j;
		       iW=(i-1)+iP_j;
               iN=i+(j+1)*inx+iP_k+1;
		       iS=i+(j-1)*inx+iP_k+1;
		       iT=i+j*inx+iP_k_p1+1;
		       iB=i+j*inx+iP_k_m1+1;



				if (evt_f[iP]>0) {

					// Узел принадлежит зоне FLUID.
					// Всем существующим соседям присвоить этот номер. 
           

					 // Узел принадлежит зоне FLUID.
					// Всем существующим соседям присвоить этот номер. 
					if (((i+1)<inx) && (iE>0) && (iE<=size_xyz) && (evt_f[iE-1]>0)&&(evt_f[iE-1]!=evt_f[iP])) {
						evt_f[iE-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i+1, j, k, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i+1;
						stek[itop].j=j;
						stek[itop].k=k;
						itop++;
					}
					if (((i-1)>=0) && (iW>0) && (iW<=size_xyz) && (evt_f[iW-1]>0)&&(evt_f[iW-1]!=evt_f[iP])) {
						evt_f[iW-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i-1, j, k, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i-1;
						stek[itop].j=j;
						stek[itop].k=k;
						itop++;			
					}
					if (((j+1)<iny) && (iN>0) && (iN<=size_xyz) && (evt_f[iN-1]>0)&&(evt_f[iN-1]!=evt_f[iP])) {
						evt_f[iN-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i, j+1, k, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i;
	        			stek[itop].j=j+1;
						stek[itop].k=k;
						itop++;
					}
					if (((j-1)>=0) && (iS>0) && (iS<=size_xyz) && (evt_f[iS-1]>0)&&(evt_f[iS-1]!=evt_f[iP])) {
						evt_f[iS-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i, j-1, k, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i;
						stek[itop].j=j-1;
						stek[itop].k=k;
						itop++;			
					}
					if (((k+1)<inz) && (iT>0) && (iT<=size_xyz) && (evt_f[iT-1]>0)&&(evt_f[iT-1]!=evt_f[iP])) {
						evt_f[iT-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i, j, k+1, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i;
						stek[itop].j=j;
						stek[itop].k=k+1;
						itop++;
					}
					if (((k-1)>=0) && (iB>0) && (iB<=size_xyz) && (evt_f[iB-1]>0)&&(evt_f[iB-1]!=evt_f[iP])) {
						evt_f[iB-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i, j, k-1, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i;
						stek[itop].j=j;
	         			stek[itop].k=k-1;
						itop++;
					}
				}	

		}

		
} // my_fill_Domain_recursive


// создание связей гидродинамики с теплопроводностью.
// Первая часть отделена 22 сентября 2016 в связи с появлением АЛИС сетки.
// Эта часть универсальна и подходит и для АЛИС сетки тоже.
void constr_ptr_temp_part1(int &flow_interior, 
	int * &evt_f, int** &evt_f2, integer* &domain_id, 
	int inx, int iny, int inz, integer &icount_part
	) {

#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_max_threads();
	//omp_set_num_threads(8); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
	unsigned int nthreads = number_cores();
	omp_set_num_threads(nthreads); // установка числа потоков
#else 
	int i_my_num_core_parallelesation = 1;
#endif

	int i = 0, j = 0, k = 0;

	// 9 мая 2013 логика проходов может не сработать на очень сложной геометрии и нужен 
	// скорее рекурсивный алгоритм заливки.
	// Рекурсия (без пользовательского стека) испытывает проблемы с STECK OVERFLOW,
	// т.е. глубина рекурсивных вызовов ограничена и её нехватает на большеразмерных задачах.



	// 26 июля 2015 рекурсивный алгоритм прохода связывающий область 
	// любой сложности, если он не встретит ограничений по памяти.
	// Не рекурсивный алгоритм минимален по памяти, но не справляется с очень сложными 
	// конфигурациями расчётной области. 
	// 22 сентября 2016. Теперь справляется т.к. там сделан цикл в 1000 проходов которых
	// точно должно хватить для любой геометрии и досрочный прерыватель если заливка уже
	// успешно выполнена. Прерыватель очень сильно сокращает время построения заливки.

	bool brecursiv_algorithm = true; // 01.03.2019.
	bool bfirst_visit = true;
	integer iPgold = 0; // модификация 26,03,2019
	const integer size_xyz=inx*iny*inz;


	if (brecursiv_algorithm) {
		// Рекурсивный алгоритм.
		// Рекурсия имитируется СТЕКОМ. Там теперь свой стековый алгоритм внутри.
		// Проблемы на геометрии большой блок вум с оребрением и воздушным охлаждением 26,03,2019.
		// integer iP_k,iP_k_p1,iP_k_m1;
        // integer iP_j, iP;
       
		STEK *stek = nullptr;
		stek = new STEK[inx*iny*inz];


		 for ( int k = 0; k<inz; ++k)
		{
			integer  iP_k=k*inx*iny;
            integer  iP_k_p1=(k+1)*inx*iny;
            integer  iP_k_m1=(k-1)*inx*iny;

			for (int  j = 0; j<iny; ++j)
			{
				 integer  iP_j=j*inx +iP_k;
				 for (int  i = 0; i<inx; ++i)
				{

					integer  iP= i + iP_j;
			//if (bfirst_visit)
			if (bfirst_visit||(evt_f[iP]!= evt_f[iPgold])) {
				  integer iE, iW, iN, iS, iT, iB;
				iE = iP +2; //(i + 1) + j*inx + k*inx*iny + 1;
				iW = iP; //(i - 1) + j*inx + k*inx*iny + 1;
				iN = iP+inx+1; //i + (j + 1)*inx + k*inx*iny + 1;
				iS = iP-inx+1; //i + (j - 1)*inx + k*inx*iny + 1;
				iT = i+j*inx+iP_k_p1+1; //i + j*inx + (k + 1)*inx*iny + 1;
				iB = i+j*inx+iP_k_m1+1;//i + j*inx + (k - 1)*inx*iny + 1;
				if ((iP >= 0) && (iP<size_xyz) && (evt_f[iP]>0)) {
					if (((i + 1)<inx) && (iE>0) && (iE <= size_xyz) && (evt_f[iE - 1]>0)) {
						if (((i - 1) >= 0) && (iW>0) && (iW <= size_xyz) && (evt_f[iW - 1]>0)) {
							if (((j + 1)<iny) && (iN>0) && (iN <= size_xyz) && (evt_f[iN - 1]>0)) {
								if (((j - 1) >= 0) && (iS>0) && (iS <= size_xyz) && (evt_f[iS - 1]>0)) {
									if (((k + 1)<inz) && (iT>0) && (iT <= size_xyz) && (evt_f[iT - 1]>0)) {
										if (((k - 1) >= 0) && (iB>0) && (iB <= size_xyz) && (evt_f[iB - 1]>0)) {
#if doubleintprecision == 1
											// printf("in comming %lld\n",evt_f[iP]); // debug.
#else
											 //printf("in comming %d\n",evt_f[iP]); // debug.
#endif
											
											//system("PAUSE");
											
											// Это только первый вызов, вся работа делается внутри алгоритма my_fill_Domain_recursive.
											bfirst_visit = false;// досрочный выход из цикла for.
											// Там внутри рекурсия имитирована стеком. На модели из 576 блоков работает в разы быстрее.
											my_fill_Domain_recursive(evt_f, i, j, k, inx, iny, inz, stek);
											iPgold = iP; // модификация 26,03,2019
										}
									}
								}
							}
						}
					}
				}
			}
				}
			}
		}
	
		 //if (stek != nullptr) {
			 // Оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		 delete[] stek;
		 stek = nullptr;
		 //}
	}
	else {

		bfirst_visit = false;
		bool *evt_etalon = nullptr;
		evt_etalon = new bool[size_xyz+2];
		if (evt_etalon == nullptr) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for evt_etalon constr struct...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		bool *evt_etalon1 = nullptr;
		evt_etalon1 = new bool[size_xyz+2];
		if (evt_etalon1 == nullptr) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for evt_etalon1 constr struct...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
			evt_etalon[i_1] = false; // initialisation.
			evt_etalon1[i_1] = false; // initialisation.
		}

		evt_etalon[0] = true;
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
			if (evt_f[i_1] > 0) {
				evt_etalon[evt_f[i_1]] = true;
			}
		}
		

		if (print_log_message) {
			printf("\n                                          ");
			printf("progress bar:\n	");
			printf("\\");
		}
#ifdef _OPENMP
		omp_set_num_threads(1); // СБОЙ параллельной обработки 29_07_2017
#endif
		// гипотеза - образуются гонки разных нитей и они не могут залить одним цветом расчётную область.
		// 5 августа 2016. (рекурсивный алгоритм вылетает из-за переполнения стека помоему).
		// И переполнение стека есть серьёзная проблема.
		for (integer ilk = 0; ilk < 1000; ilk++) {
			bool bisequal = true;

			//if (ilk % 100 == 0) printf("progress bar %e procent\n", ilk / 10.0);
			if (print_log_message) {
				printf("\b");
				printf("|");
			}
			// Первый проход: Заливка.
			// 6*4 == 24 обхода
			// Меняем всевозможные направления обхода расчётной области,
			// чтобы справится со случаями "заглядывания за угол".
			// Также нужен прерыватель т.к. 1000 проходов это очень долго для tgf20.

#pragma omp parallel for
			for (i = 0; i < inx; i+=2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i+=2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("/");
			}
			if (size_xyz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}

#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("--");
			}
			if (size_xyz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b\b");
				printf("\\");
			}
			if (size_xyz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("|");
			}
			if (size_xyz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("/");
			}
			if (size_xyz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}

#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("--");
			}
			if (size_xyz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b\b");
				printf("\\");
			}
			if (size_xyz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("|");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("/");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("--");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b\b");
				printf("\\");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("|");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("/");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("--");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b\b");
				printf("\\");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("|");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("/");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("--");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b\b");
				printf("\\");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("|");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("/");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b");
				printf("--");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				printf("\b\b");
				printf("\\");
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			if (print_log_message) {
				//printf("\b");
				//printf("|");
			}

			// Код досрочного прерывателя, если всё готово.
			// Иначе процесс сканирования может занять слишком много времени.
			// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
			for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
				evt_etalon1[evt_f[i_1]]=true;
			}

			bisequal = true;
			for (integer i_1 = 0; i_1 < size_xyz; ++i_1) {
				if (evt_etalon[i_1] != evt_etalon1[i_1]) {
					bisequal = false;
				}
			}
			if (bisequal) {
				break;
			}
			else {
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx * iny * inz; ++i_1) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}

			//printf("ilk==%d\n",ilk);
		}
#ifdef _OPENMP
		int inum_core = number_cores();
		omp_set_num_threads(inum_core);// i_my_num_core_parallelesation
#endif
		if (print_log_message) {
			printf("\n");
		}
		//system("PAUSE");

		
		delete[] evt_etalon;
		evt_etalon = nullptr;

		if (evt_etalon1 != nullptr) {
			delete[] evt_etalon1;
		}
		evt_etalon1 = nullptr;
	}

	if (0&&bfirst_visit) {
		printf("ERROR!!!: NOT FILLING\n");
		system("PAUSE");
	}

	if (print_log_message) {
		//printf("part 11\n");
#if doubleintprecision == 1
		printf("part %lld my_fill_Domain constr_ptr_temp_part1 \n", icount_part++); //22
#else
		printf("part %d my_fill_Domain constr_ptr_temp_part1 \n", icount_part++); //22
#endif
	}

	//integer iP=0;


	/*
	// debug GOOD
	// Проверка правильности распознавания по слоям.
	for (k=0; k<inz; ++k) {
	    for (i=0; i<inx; ++i) {
	        for (j=0; j<iny; ++j) {
	             iP=i+j*inx+k*inx*iny;
				 #if doubleintprecision == 1
						printf("%lld ",evt_f[iP]);
				 #else
						printf("%d ",evt_f[iP]);
				 #endif
	             
	        }
	        printf("\n");
	    }
	    system("PAUSE");
	    printf("\n\n");
	}
	*/
#if doubleintprecision == 1
	//printf("evt_f[86]=%lld, evt_f[124]=%lld\n",evt_f[86],evt_f[124]); // debug
#else
	//printf("evt_f[86]=%d, evt_f[124]=%d\n",evt_f[86],evt_f[124]); // debug
#endif
	
	//system("PAUSE");

	// Нам понадобтся знать информацию о различных подобластях
	// для процедуры mass_balance.
	// Т.к. процедуру балансировки масс надо делать для каждой 
	// отдельной связанной гидродинамической
	// подобласти иначе не будет сходимости.
	// Чтобы распознать цвет связаной подобласти он будет хранится в evt_f2[2][iP] позиции. 
	int* evt_f_shadow = nullptr;
	evt_f_shadow = new int[size_xyz];
	// оператор new не требует проверки на null
	//if (evt_f_shadow == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for evt_f_shadow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// Копирование.
	//for (i = 0; i < inx; ++i) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
		//iP = i + j*inx + k*inx*iny;
#pragma omp parallel for
	for (integer iP_loc=0; iP_loc < size_xyz; iP_loc++) {
		evt_f_shadow[iP_loc] = evt_f[iP_loc];
	}

	// Второй проход:
	// Если возникает ошибка, то значение max_domain нужно увеличить.	
	const integer max_domain = 2048;// 256; // максимальное количество зон FLUID
	//integer *domain_id = nullptr;
	//domain_id = nullptr;
	domain_id = new integer[max_domain];
	// оператор new не требует проверки на null.
	//if (domain_id == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for domain_id constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	for (i = 0; i<max_domain; ++i) domain_id[i] = 0; // инициализация

	integer *domain_id_shadow = nullptr;
	domain_id_shadow = new integer[max_domain];
	// оператор new не требует проверки на null.
	//if (domain_id_shadow == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for domain_id_shadow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	for (integer  i = 0; i<max_domain; ++i) domain_id_shadow[i] = 0; // инициализация


	integer id;
	bool bfind = false;
	bool bfind_shadow = false;
	//integer l;
	 int ic = 0; // счётчик по domain_id
	int ic_shadow = 0; // счётчик по domain_id_shadow
	flow_interior = 0;
	// Если жидкие ячейки вообще были.
	if (!bfirst_visit) {
		// integer iP_k, iP_j;

		STEK *stek = nullptr;
		stek = new STEK[inx*iny*inz];

		for (int  k = 0; k < inz; ++k)
		{
			integer  iP_k=k * inx * iny;
			for (int  j = 0; j < iny; ++j)
			{
				integer  iP_j=iP_k+j * inx;
				for (int  i = 0; i < inx; ++i) 
				{

					integer  iP = i + iP_j;
					if (evt_f[iP] > 0) {
						id = evt_f[iP]; // идентификатор связанной FLUID зоны.
						bfind = false;
						// Внимание!!! возможно это медленный участок кода!.
						for (integer l = 0; ((l < ic) && (l < max_domain)); l++) if (domain_id[l] == id) bfind = true;
						if (!bfind) {
							// Пробуем выполнить дозаливку. Добавка 29.04.2019.
							

							my_fill_Domain_recursive(evt_f, i, j, k, inx, iny, inz, stek);

							if (print_log_message) {

								std::cout << "patch 29.04.2019. fluid zone id number = " << ic << ". start zone control volume number = " << id << " " << std::endl;
							}
							if (ic >= max_domain - 1) {
								std::cout <<  "error! nado uvelichit max_domain count..." << std::endl;
								std::cout <<  "icount domain =="<< ic << std::endl; 
								//system("PAUSE");
								system("PAUSE");
								exit(1);
							}
							// Запоминаем идентификатор связанной FLUID зоны.
							domain_id[ic++] = id;
							// Количество изолированных FLUID зон.
							flow_interior = ic;
						}
						bfind_shadow = false;
						// Внимание!!! возможно это медленный участок кода!.
						for (integer l = 0; ((l < ic_shadow) && (l < max_domain)); l++) if (domain_id_shadow[l] == id) bfind_shadow = true;
						if (!bfind_shadow) {
						if (ic_shadow >= max_domain - 1) {
							std::cout <<  "error! nado uvelichit max_domain count..." << std::endl;
							std::cout <<  "icount domain ==" << ic << std::endl;
							//system("PAUSE");
							system("PAUSE");
							exit(1);
						}
						domain_id_shadow[ic_shadow++] = id;
						//flow_interior = ic;
						}
					}
				}
			}
		}
	
		//if (stek != nullptr) {
			// Оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] stek;
		stek = nullptr;
		//}
	}
	if (flow_interior > 1) {
		std::cout <<  "WARNING: flow_interior count = " << flow_interior << std::endl;
		std::cout <<  "Your model contains multiple fluid zones." << std::endl;
		//system("PAUSE");
	}
	if (print_log_message) {
		//printf("part 12\n");
		std::cout << "part " << icount_part++ << "  my fill Domain recursive constr_ptr_temp_part1 " << std::endl; //22
	}

	if (1) {
		// В результате domain_id содержит единственную FLUID зону и её идентификатор 2.

		if (print_log_message) {
			std::cout << "5_08_2016: several fluid domain combine into one fluid domain." << std::endl;
		}
		// Смысл в том что даже если у нас несколько гидрдинамических подобластей то мы делаем одну единственную.
		//ic = 0;
		/*
		// Мутный код, логика ПЛОХАЯ. Исправлено 22 сентября 2016.
		for (i = 0; i < inx; ++i) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
		    iP = i + j*inx + k*inx*iny;
		    if (evt_f[iP] > 0) {
		       evt_f[iP] = 2;
		       id = evt_f[iP];
		       bfind = false;
		       // Внимание!!! возможно это медленный участок кода!.
		       for (l = 0; l < max_domain; l++) if (domain_id[l] == id) bfind = true;
		       if (!bfind) {
		          if (ic >= max_domain - 1) {
		             printf("error! nado uvelichit max_domain count...\n");
		             system("PAUSE");
		             exit(1);
		          }
		          domain_id[ic++] = id;
		          flow_interior = ic;
		       }
		    }
		}
		*/
		const integer size_xyz=inx * iny*inz;
		// В результате domain_id содержит единственную FLUID зону и её идентификатор 2.
		// Если жидкие ячейки вообще были.
		if (!bfirst_visit) {
			//for (i = 0; i < inx; ++i) for (j = 0; j < iny; ++j) for (k = 0; k < inz; ++k) {
				//iP = i + j * inx + k * inx * iny;
#pragma omp parallel for
            for (integer iP_loc=0; iP_loc< size_xyz; iP_loc++)
			{
				if (evt_f[iP_loc] > 0) {
					evt_f[iP_loc] = 2;
				}
			}
		}

		for (integer l = 0; l < max_domain; l++) domain_id[l] = 0;
		ic = 0;
		domain_id[ic++] = 2;
		flow_interior = ic;

		if (print_log_message) {

			std::cout << "part " << icount_part++ << "  evt_f init  constr_ptr_temp_part1 ";
			std::cout << " 05_08_2016; 22_09_2016." << std::endl;
		}
	}



	//if (domain_id != nullptr) {
	// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] domain_id;
		domain_id = nullptr;
	//}

	//if (domain_id_shadow != nullptr) {
	// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] domain_id_shadow;
		domain_id_shadow = nullptr;
	//}

	// Освобождение памяти.
	//if (evt_f_shadow != nullptr) {
	// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] evt_f_shadow;
	//}
	evt_f_shadow = nullptr;

#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif
} // constr_ptr_temp_part1

//  откладываем вычисление evt_f2 как можно дальше т.к. оно не проходит 
// по памяти на больших по размеру  моделях.
// Здесь содержится код вычисляющий только evt_f2.
void constr_ptr_temp_part2(int &flow_interior,
	int * &evt_f, int** &evt_f2, integer* &domain_id,
	int inx, int iny, int inz, 
	TOCKA_SHORT_INT* &tck_int_list, int &maxelm, integer &icount_part) {

	//integer i = 0, j = 0, k = 0;	
	
	// Более быстрое копирование.
	const integer isize = inx * iny*inz;	

	// Второй проход:
	// Если возникает ошибка, то значение max_domain нужно увеличить.
	int max_domain = 2048; // максимальное количество зон FLUID
	//integer *domain_id = nullptr;

	
	domain_id = new integer[max_domain];
	// оператор new не требует после себя проверки на null
	for (int i = 0; i<max_domain; ++i) domain_id[i] = 0; // инициализация

	STEK *stek = nullptr;
	stek = new STEK[inx*iny*inz];

	integer id;
	bool bfind;
	integer l;
	integer ic = 0; // счётчик по domain_id
	flow_interior = 0;
	//for (i = 0; i<inx; ++i) for (j = 0; j<iny; ++j) for (k = 0; k<inz; ++k) {
		//integer iP = i + j*inx + k*inx*iny;
		//if (evt_f[iP]>0) {
	for (integer iscan = 0; iscan<maxelm; iscan++) 
		{
			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;
			integer iP = i + j * inx + k * inx*iny;
			{

			id = evt_f[iP]; // идентификатор связанной FLUID зоны.
			bfind = false;
			// Внимание!!! возможно это медленный участок кода!.
			for (l = 0; ((l < ic) && (l<max_domain)); l++) if (domain_id[l] == id) bfind = true;
			if (!bfind) {
				// Пробуем выполнить дозаливку. Добавка 29.04.2019.
				// 2
				

				my_fill_Domain_recursive(evt_f, i, j, k, inx, iny, inz, stek);

				if (print_log_message) {

					std::cout << "patch 29.04.2019. fluid zone id number = " << ic << ". start zone control volume number = " << id << " " << std::endl;
				}
				if (ic >= max_domain - 1) {
					std::cout << "error! nado uvelichit max_domain count..." << std::endl;
					std::cout << "icount domain ==" << ic << std::endl;
					system("PAUSE");
					exit(1);
				}
				// Запоминаем идентификатор связанной FLUID зоны.
				domain_id[ic++] = id;
				// Количество изолированных FLUID зон.
				flow_interior = (int)(ic);
			}
			
		}
	}

	//if (stek != nullptr) {
		// Оператор delete может быть вызван повторно в том числе и к нулевому указателю.
	delete[] stek;
	stek = nullptr;
	//}

	if (flow_interior > 1) {
		std::cout << "WARNING: flow_interior count = " << flow_interior << std::endl;
		std::cout << "Your model contains multiple fluid zones." << std::endl;
		//system("PAUSE");
	}

	if (print_log_message) {
		//printf("part 12.2\n");
		std::cout << "part " << icount_part++ << " constr domain_id constr_ptr_temp_part2 " << std::endl;
	}
	if (1) {
		// В результате domain_id содержит единственную FLUID зону и её идентификатор 2.

		if (print_log_message) {
			std::cout << "5_08_2016: several fluid domain combine into one fluid domain." << std::endl;
			// Смысл в том что даже если у нас несколько гидрдинамических подобластей то мы делаем одну единственную.
		}
		
		// В результате domain_id содержит единственную FLUID зону и её идентификатор 2.
		// Если жидкие ячейки вообще были.
		
			
#pragma omp parallel for
		for (integer iP_loc = 0; iP_loc < isize; iP_loc++)
		{
			if (evt_f[iP_loc] > 0) {
				evt_f[iP_loc] = 2;
			}
		}
		

		for (integer l = 0; l < max_domain; l++) domain_id[l] = 0;
		ic = 0;
		domain_id[ic++] = 2;
		flow_interior = (int)(ic);

		if (print_log_message) {
			std::cout << "part " << icount_part++ << "  evt_f init  constr_ptr_temp_part2 ";
			std::cout << " 05_08_2016; 22_09_2016." << std::endl;
		}
	}


#if doubleintprecision == 1
	// printf("flow interior=%lld\n",flow_interior); //debug.
#else
	// printf("flow interior=%d\n",flow_interior); //debug.
#endif
	
	// system("PAUSE");

	

	// Выделение оперативной памяти
	evt_f2 = nullptr;
	evt_f2 = new int*[2];
	if (evt_f2 == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout <<  "Problem: not enough memory on your equipment for evt_f2 constr struct..." << std::endl;
		std::cout <<  "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (int i = 0; i<2; ++i) evt_f2[i] = nullptr;
	for (int i = 0; i<2; ++i) {
		evt_f2[i] = new int[isize];
		if (evt_f2[i] == nullptr) {
			// недостаточно памяти на данном оборудовании.
			std::cout <<  "Problem: not enough memory on your equipment for evt_f2["<< i <<"] constr struct..." << std::endl;
			std::cout <<  "Please any key to exit..." << std::endl;
			// system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// Проход по всем контрольным объёмам
	//integer ipos = -1;

	for (integer iP = 0; iP < isize; ++iP) {
		// инициализация.
		evt_f2[MASKDOMAINFLUID][iP] = -1; // не принадлежит расчётной области.
	}
	
	for (integer iscan = 0; iscan<maxelm; iscan++) 
	{
			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;
			integer iP = i + j * inx + k * inx*iny;
			// iP находитсятолько внутри гидродинамической области.

			//12.01.2020 Здесь была какая то непонятная логика. Видимо остаток от 
			// предыдущих устаревших вариантов кода.
			evt_f2[MASKDOMAINFLUID][iP] = 1; // Только одна Fluid зона. номер 1		
	}

	if (print_log_message) {
		//printf("part 13\n");
		std::cout << "part " << icount_part++ << "  constr evt_f2 constr_ptr_temp_part2 " << std::endl;
	}

	/*
	// debug GOOD
	// Проверка правильности распознавания по слоям.
	for (k=0; k<inz; ++k) {
	for (i=0; i<inx; ++i) {
	for (j=0; j<iny; ++j) {
	iP=i+j*inx+k*inx*iny;
	#if doubleintprecision == 1
		printf("%lld ",evt_f2[MASKDOMAINFLUID][iP]);
	#else
		printf("%d ",evt_f2[MASKDOMAINFLUID][iP]);
	#endif
	
	}
	printf("\n");
	}
	system("PAUSE");
	printf("\n\n");
	}
	//*/

	// Освобождение памяти.

	// Здесь нельзя освобождать озу из под domain_id.
	// domain_id используется дальше по коду.
	// domain_id передавалась из вне внутрь данной процедуры.
	//delete[] domain_id;
	//domain_id = nullptr;

	
	

} // constr_ptr_temp_part2

// Выделение памяти под FLOW.
// Дата создания 22 сентября 2016.
// Данная функция лишь несёт нагрузку выделения оперативной памяти, а так из неё всё выгружено при
// АЛИС сетке.
void constr_ptr_temp_allocation_memory_alice(int &flow_interior, FLOW* &f, integer &icount_part)
{
	if (f == nullptr) {
		if (flow_interior != 1) {
			std::cout <<  "error!!! flow_interior!=1" << std::endl;
			system("PAUSE");
			exit(1);
		}
		f = new FLOW[flow_interior];

		// проход по всем жидким зонам.
		for (integer l = 0; l < flow_interior; l++) {
			// Число ненулевых КО.
			// Заглушка которую возможно можно будет оживить в дальнейшем.
			// Самое простое это не использовать информацию whot_is_block для гидродинамических подобластей.
			f[l].whot_is_block = nullptr;
			// Количество элементов в модели будет подсчитано в другом месте в специальной функции.
			//f[l].maxelm = domain_counter[l];
#if doubleintprecision == 1
		//printf("FLUID%lld maxelm=%lld\n", l, f[l].maxelm);
#else
		//printf("FLUID%d maxelm=%d\n", l, f[l].maxelm);
#endif

		// Выделение оперативной памяти.
			f[l].ptr = nullptr;
			// f[l].ptr[fluid_elm_id]=temper_elm_id; // Смысл ptr.
			// Все закомментированное заполняется в другом месте !!!.
			/*
			f[l].ptr = new integer[f[l].maxelm];
			#if doubleintprecision == 1
				printf("delete flow %lld ptr\n", l);
			#else
				printf("delete flow %d ptr\n", l);
			#endif

			if (f[l].ptr == nullptr) { // -2N
				// недостаточно памяти на данном оборудовании.
				#if doubleintprecision == 1
					printf("Problem: not enough memory on your equipment for ptr flow%lld constr struct...\n", l);
				#else
					printf("Problem: not enough memory on your equipment for ptr flow%d constr struct...\n", l);
				#endif

				printf("Please any key to exit...\n");
				//system("PAUSE");
				system("pause");
				exit(1);
			}
			for (j = 0; j < maxelm_t; ++j) {
				if (ptr[MASKDOMAINFLUID][j] == l) {
					f[l].ptr[ptr[ENUMERATECONTVOL][j]] = j;
				}
			}
			*/
		}
	}
	else {
		std::cout <<  "error!!! realloc f" << std::endl;
		system("PAUSE");
		exit(1);
	}
	if (print_log_message) {
		//printf("part 14\n");
		//printf("part 15\n");
		//printf("part 16\n");
		std::cout << "part " << icount_part++ << " constr_ptr_temp_allocation_memory_alice" << std::endl;
		std::cout << "part " << icount_part++ << " constr_ptr_temp_allocation_memory_alice" << std::endl;
		std::cout << "part " << icount_part++ << " constr_ptr_temp_allocation_memory_alice" << std::endl;
	}
} // constr_ptr_temp_allocation_memory_alice


// создание связей гидродинамики с теплопроводностью.
void constr_ptr_temp(int &flow_interior, FLOW* &f, integer maxelm_t,
					 int** &ptr, int *evt_t, int * &evt_f, int** &evt_f2,
					 integer* &domain_id, int inx, int iny, int inz, 
	                 bool breconstruct, integer& icount_part) {
	
    

	
	//integer iP;

	// Если возникает ошибка, то значение max_domain нужно увеличить.
	integer max_domain = 256; // максимальное количество зон FLUID

	int id;
	//bool bfind;
	//bool bfind_shadow;
	int l;
	//integer ic = 0; // счётчик по domain_id
	//integer ic_shadow = 0; // счётчик по domain_id_shadow
	//flow_interior = 0; // Заполняется в constr_ptr_temp_part1.
	
	

	 int i1; // счётчик
	 	
	 int ipos = -1;
	 //integer ipos_shadow = -1;

	// 22 сентября 2016.
	// До части 13 можно использовать и при АЛИС сетке.

	// Выделение памяти:

	 //printf("delete temperature ptr\n");
	 if (ptr != nullptr) {
		 for (unsigned char i = 0; i<2; ++i) {
			 if (ptr[i] != nullptr) {
				 delete[] ptr[i];
				 ptr[i] = nullptr;
			 }
		 }
		 delete[] ptr;
		 ptr = nullptr;
	 }	 

	ptr=nullptr;
    ptr = new int*[2];
	if (ptr==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout <<  "Problem: not enough memory on your equipment for ptr constr struct..." << std::endl;
		std::cout <<  "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	//for (unsigned char  i=0; i<2; ++i) ptr[i]=nullptr;

	for (unsigned char  i=0; i<2; ++i) {
		ptr[i]=new int[maxelm_t];
		// оператор new не требует после себя проверки на null
		//if (ptr[i]==nullptr) {
	        // недостаточно памяти на данном оборудовании.
//#if doubleintprecision == 1
	//		printf("Problem: not enough memory on your equipment for ptr[%d] constr struct...\n", i);
//#else
	//		printf("Problem: not enough memory on your equipment for ptr[%d] constr struct...\n", i);
//#endif
		    
	//	    printf("Please any key to exit...\n");
			//system("PAUSE");
		//	system("pause");
		  //  exit(1);
	    //}   
	}

	for (unsigned char i = 0; i < 2; ++i) {
		for (integer i87 = 0; i87 < maxelm_t; i87++) {
			ptr[i][i87] = -1; // инициализация твердым телом.
		}
	}

	//integer i = 0, j = 0, k = 0;


    //#define ENUMERATECONTVOL 0 // слой нумерации КО
    //#define MASKDOMAINFLUID 1 // слой маска определяет зону FLUID

	if ((ptr != nullptr) && (ptr[0] != nullptr) && (ptr[1] != nullptr)) {

		// Третий проход. Маркировка отдельных зон FLUID:
		l = 0;
		ipos = -1;
		//  integer iP_k, iP_j;
		for (int  k = 0; k < inz; ++k)
		{
			integer  iP_k = k * inx * iny;
			for (int  j = 0; j < iny; ++j)
			{
				integer  iP_j = j * inx + iP_k;
				for (int  i = 0; i < inx; ++i)
				{

					integer iP = i + iP_j;
					if (evt_t[iP] > 0) {
						// Если КО принадлежит модели теплопроводности.
						if (evt_f[iP] > 0) {
							// FLUID Domain
							id = evt_f[iP];
							// Внимание !!! осторожно возможно медленный участок кода.
							for (i1 = 0; i1 < max_domain; i1++) if (domain_id[i1] == id) {
								ipos = i1; // нумерация начинается с нуля.
								break;
							}
							if (l < maxelm_t) {
								// там инициализировано -1, поэтому нумерация начинается с 0.
								ptr[MASKDOMAINFLUID][l] = ipos;// нумерация начинается с нуля.
							}
							else {
								std::cout << "perepolnenie ptr[MASKDOMAINFLUID][l] l>=" << maxelm_t << std::endl;
								system("pause");
								exit(1);
							}
						}
						else {
							// не fluid но солид (и не холлоу)
							if (ptr != nullptr) {
								if (ptr[MASKDOMAINFLUID] != nullptr) {
									if (l < maxelm_t) {
										ptr[MASKDOMAINFLUID][l] = -1;
									}
									else {
										std::cout << "perepolnenie 2 ptr[MASKDOMAINFLUID][l] l>=" << maxelm_t << std::endl;
										system("pause");
										exit(1);
									}
								}
								else {
									std::cout << "error: null ptr[MASKDOMAINFLUID] in constr_ptr_temp." << std::endl;
									system("pause");
									exit(1);
								}
							}
							else {
								//ptr == nullptr;
								std::cout << "error: null ptr in constr_ptr_temp." << std::endl;
								system("pause");
								exit(1);
							}
						}

						l++;
					}
				}
			}
		}

		if (print_log_message) {
			//printf("part 14\n");
			std::cout << "part " << icount_part++ << "  construct ptr[MASKDOMAINFLUID][] constr_ptr_temp " << std::endl;
		}

		// Счётчик внутренних КО в каждой из domain FLUID.
		int *domain_counter = nullptr;
		domain_counter = new int[max_domain];
		// результат применения оператора new не надо проверять на null.
		//if (domain_counter == nullptr) {
			// недостаточно памяти на данном оборудовании.
			//printf("Problem: not enough memory on your equipment for domain_counter constr struct...\n");
			//printf("Please any key to exit...\n");
			//system("PAUSE");
			//system("pause");
			//exit(1);
		//}
		//else 
		{

			for (integer  i = 0; i < max_domain; ++i) domain_counter[i] = 0; // инициализация

			if (ptr != nullptr) {
				if ((ptr[MASKDOMAINFLUID] != nullptr)&&(ptr[ENUMERATECONTVOL]!=nullptr)) {
					for (integer  i = 0; i < maxelm_t; ++i) {
						if (ptr[MASKDOMAINFLUID][i] == -1) {
							ptr[ENUMERATECONTVOL][i] = -1;
						}
						else {
							ptr[ENUMERATECONTVOL][i] = domain_counter[ptr[MASKDOMAINFLUID][i]];
							domain_counter[ptr[MASKDOMAINFLUID][i]]++;
						}
					}
				}
				else {
					std::cout <<  "memory for ptr[MASKDOMAINFLUID] no allocate. error." << std::endl;
					std::cout <<  "see function constr_ptr_temp in module constr struct" << std::endl;
					system("pause");
					exit(1);
				}
			}
			else {
				std::cout <<  "memory for ptr no allocate. error." << std::endl;
				std::cout <<  "see function constr_ptr_temp in module constr struct." << std::endl;
				system("pause");
				exit(1);
			}


			if (!breconstruct) {
				if (flow_interior > 0) {
					if (f == nullptr) {
						f = new FLOW[flow_interior];						
					}
					else {
						std::cout <<  "error: f != nullptr" << std::endl;
						std::cout <<  "file constr_struct. function constr_ptr_temp." << std::endl;
						std::cout <<  "There will be no reallocation of RAM. What is" << std::endl;
						std::cout <<  "already highlighted will be used." << std::endl;
						delete[] f;
						f = new FLOW[flow_interior];
						//system("PAUSE");
						//exit(1);
					}
					for (l = 0; l < flow_interior; l++) {
						f[l].ptr = nullptr;
						f[l].nvtx = nullptr;
						f[l].pa = nullptr;
						f[l].neighbors_for_the_internal_node = nullptr;
						f[l].border_neighbor = nullptr;
						f[l].whot_is_block = nullptr;
						f[l].potent = nullptr;
						f[l].prop = nullptr;
						f[l].prop_b = nullptr;
						f[l].alpha = nullptr;
						f[l].slau = nullptr;
						f[l].slau_bon = nullptr;
						f[l].diag_coef = nullptr;
						f[l].iN = nullptr;
						f[l].id = nullptr;
						f[l].rdistWall = nullptr;
						f[l].SInvariantStrainRateTensor = nullptr;
						f[l].mf = nullptr;
						f[l].icolor_different_fluid_domain = nullptr;
						f[l].ifrontregulationgl = nullptr;
						f[l].ibackregulationgl = nullptr;
					}
				}
				else {
					f = nullptr;
				}
			}
			if (f != nullptr) {
				// проход по всем жидким зонам.
				for (l = 0; l < flow_interior; ++l) {
					// Число ненулевых КО.
					// Заглушка которую возможно можно будет оживить в дальнейшем.
					// Самое простое это не использовать информацию whot_is_block для гидродинамических подобластей.
					f[l].whot_is_block = nullptr;
					f[l].maxelm = domain_counter[l];
					if (print_log_message) {
						std::cout << "FLUID" << l << " maxelm=" << f[l].maxelm << std::endl;
					}

					// Выделение оперативной памяти.
					if (f[l].ptr != nullptr) {
						delete[] f[l].ptr;
					}
					f[l].ptr = nullptr;
					f[l].ptr = new int[f[l].maxelm];

					if (print_log_message) {
						std::cout << "delete flow " << l << " ptr" << std::endl;
					}
					
					if (f[l].ptr == nullptr) { // -2N
						// недостаточно памяти на данном оборудовании.

						std::cout << "Problem: not enough memory on your equipment for ptr flow" << l << " constr struct..." << std::endl;
						std::cout << "Please any key to exit..." << std::endl;
						//system("PAUSE");
						system("pause");
						exit(1);
					}
					for (int j = 0; j < maxelm_t; ++j) {
						if (ptr[MASKDOMAINFLUID][j] == l) {
							f[l].ptr[ptr[ENUMERATECONTVOL][j]] = j;
						}
					}
				}
			}

			if (print_log_message) {
				//printf("part 15\n");
				std::cout << "part " << icount_part++ << " construct ptr[ENUMERATECONTVOL][] constr_ptr_temp " << std::endl;
			}

			// Счётчик внутренних КО на дополнительной структуре.
			for (integer  i = 0; i < max_domain; ++i) domain_counter[i] = 0; // инициализация

			//  integer iP_k, iP_j;

			for (int  k = 0; k < inz; ++k)
			{
				integer  iP_k = k * inx * iny;
				for (int  j = 0; j < iny; ++j)
				{
					integer  iP_j = j * inx + iP_k;
					for (int  i = 0; i < inx; ++i)
					{
						integer  iP = i + iP_j;
						if (evt_f2[MASKDOMAINFLUID][iP] <= 0) {
							evt_f2[ENUMERATECONTVOL][iP] = 0; // не принадлежит расчётной области
						}
						else {
							// Уникальная нумерация узлов гидродинамической подобласти.
							// Нумерация начинается с единицы.
							evt_f2[ENUMERATECONTVOL][iP] = domain_counter[evt_f2[MASKDOMAINFLUID][iP] - 1] + 1;
							domain_counter[evt_f2[MASKDOMAINFLUID][iP] - 1] = domain_counter[evt_f2[MASKDOMAINFLUID][iP] - 1] + 1;
						}
					}
				}
			}


			/*
			// debug GOOD
			// Проверка правильности распознавания по слоям.
			for (k=0; k<inz; ++k) {
				for (i=0; i<inx; ++i) {
					for (j=0; j<iny; ++j) {
						 iP=i+j*inx+k*inx*iny;
						 #if doubleintprecision == 1
								printf("%lld ",evt_f2[ENUMERATECONTVOL][iP]);
						 #else
								printf("%d ",evt_f2[ENUMERATECONTVOL][iP]);
						 #endif
						 
					}
					printf("\n");
				}
				system("PAUSE");
				printf("\n\n");
			}
			//*/


			//if (domain_counter != nullptr) {
			// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
				delete[] domain_counter;
				domain_counter = nullptr;
			//}


		}

		//if (domain_id != nullptr) {
		// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
			delete[] domain_id;
			domain_id = nullptr;
		//}
		
			if (print_log_message) {
				//printf("part 16\n");
				std::cout << "part " << icount_part++ << "  construct evt_f2  constr_ptr_temp " << std::endl;
			}
	}
} // constr_ptr_temp

// реализация функции SetLength из Delphi 
// для изменения размеров динамического массива ra. 
// Возвращает модифицированный массив ra необходимого 
//  размера. 
// Во время работы функции утечки памяти не происходит.
void SetLength_point(TOCHKA* &ra, integer isizeold, integer isize) 
{
	// isize - новая длина динамического массива.
	TOCHKA *temp=nullptr;
	temp = new TOCHKA[isize];
	// оператор new не требует после себя проверки на null
	//if (temp==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for SetLength_pointeger constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
    integer i=0;
	for (i=0; i<isize; ++i) {
		temp[i].x=0.0; // инициализация
		temp[i].y=0.0;
		temp[i].z=0.0;
	}
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
	for (i=0; i<isizeold; ++i) temp[i]=ra[i];
	
	//if (ra != nullptr) {
		// Оператор delete может быть применен повторно в том числе и к нулевому указателю.
		delete[]  ra; // уничтожение объекта
		ra = nullptr;
	//}
	ra = new TOCHKA[isize]; // выделение памяти
	for (i=0; i<isize; ++i) ra[i]=temp[i]; // копирование
	//if (temp != nullptr) {
		// Оператор delete может быть применен повторно в том числе и к нулевому указателю.
		delete[] temp; // освобождение памяти
		temp = nullptr;
	//}
	
} // SetLength_point

// Содержался линейный поиск. Внимание медленная версия.
// добавляет уникальную точку к массиву pa
// длины il
void addpoint_old(TOCHKA* &pa,integer &maxnode,TOCHKA pnew, integer* &ent, integer node, bool* &bvisit) {
	//integer i=0;
	bool bfind=false;
	/* // медленный участок кода
	for (i=0; i<maxnode; ++i) { 
		if ((fabs(pa[i].x-pnew.x)<admission) && (fabs(pa[i].y-pnew.y)<admission) && (fabs(pa[i].z-pnew.z)<admission)) {
			bfind=true;
		}
	}*/
	bfind=bvisit[node];

	if (!bfind) {
		// добавление, точка ещё не встречалась
		//SetLength_point(pa, maxnode, maxnode+1);
		pa[maxnode]=pnew;
		ent[node]=maxnode+1; // нумерация начинается с единицы
		bvisit[node]=true;
		maxnode++;
#if doubleintprecision == 1
		//printf("maxnode=%lld\n",maxnode);
#else
		//printf("maxnode=%d\n",maxnode);
#endif
		
	}
	else {
		// добавок 16.04.2017.
		/*
		// медленный участок кода
		integer ifound = -1;
		for (integer i = 0; i<maxnode; ++i) {
			if ((fabs(pa[i].x - pnew.x)<admission) && (fabs(pa[i].y - pnew.y)<admission) && (fabs(pa[i].z - pnew.z)<admission)) {
				bfind = true;
				ifound = i;
				break;
			}
		}
		ent[node] = ifound+1;
		*/
	}
} // addpoint_old


  // Компактная версия 4 июня 2017.
  // добавляет уникальную точку к массиву pa
void addpoint(TOCHKA* &pa, int &maxnode, TOCHKA_FLOAT pnew, int* &ent, integer node, bool* &bvisit) {

	if (!bvisit[node]) {
		// добавление, точка ещё не встречалась

		TOCHKA pnew_tmp;
		pnew_tmp.x = pnew.x;
		pnew_tmp.y = pnew.y;
		pnew_tmp.z = pnew.z;


		pa[maxnode] = pnew_tmp;
		maxnode++;
		ent[node] = maxnode; // нумерация начинается с единицы
		bvisit[node] = true; // индексация хеша.

	}

} // addpoint

// Компактная версия 4 июня 2017.
  // добавляет уникальную точку к массиву pa
void addpoint(TOCHKA* &pa, int &maxnode, TOCHKA pnew, int* &ent, integer node, bool* &bvisit) {

	if (!bvisit[node]) {
		// добавление, точка ещё не встречалась
		
		pa[maxnode] = pnew;
		maxnode++;
		ent[node] = maxnode; // нумерация начинается с единицы
		bvisit[node] = true; // индексация хеша.
		
	}
	
} // addpoint


// создаёт массив узлов принадлежащих расчётной области
void constr_nodes(TOCHKA* &pa, int &maxnode, int* &ent,
				  integer iflag, int* &whot_is_block, int* &evt,
	int inx, int iny, int inz,
	doublereal *xpos, doublereal *ypos, doublereal *zpos, 
	BLOCK* b, int lb, TOCKA_SHORT_INT* &tck_int_list, int &maxelm) {
	// iflag - принимает два значения: TEMPERATURE или HYDRODINAMIC и 
	// указывает с какие уравнения будут решаться и в какой расчётной области.

	// ent - глобальная нумерация узлов.
	// нумерация узлов начинается с единицы, если узел не принадлежит
	// расчётной области то ставится значение 0.
	const integer iSIZE = (inx + 1)*(iny + 1)*(inz + 1);
	ent = nullptr;
	ent = new int[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (ent == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for ent constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// был ли добавлен узел с номером node
	bool *bvisit = nullptr;
	bvisit = new bool[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (bvisit == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// Выделение памяти. Сразу много для наискорейшего выполнения:
	pa = nullptr;
	pa = new TOCHKA[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (pa == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for pa constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}


	//TOCHKA p;
#pragma omp parallel for
	for (integer i = 0; i<iSIZE; ++i) {
		ent[i] = 0;
		bvisit[i] = false;
	}

	bool inDomain = false;
	integer ib, node;
	maxnode = 0;// maxnode - текущая длина массива pa
	// Цикл по всем контрольным объёмам
	//for (k = 0; k<inz; ++k) for (j = 0; j<iny; ++j)  for (i = 0; i<inx; ++i) {
	for (integer iscan = 0; iscan<maxelm; iscan++) {

			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;

			integer k_1 = k*inx*iny;
			integer k_2 = k*(inx + 1)*(iny + 1);
			integer k_3= (k + 1)*(inx + 1)*(iny + 1);
			integer j_2 = k_2 + j*(inx + 1);
			integer j_3 = (j + 1)*(inx + 1) + k_2;
			integer j_4 = j*(inx + 1) + k_3;
			integer j_5 = (j + 1)*(inx + 1) + k_3;
			integer j_1 = k_1 + j*inx;
			integer i_1 = i + j_1;

		//p.x = 0.5*(xpos[i] + xpos[i + 1]);
		//p.y = 0.5*(ypos[j] + ypos[j + 1]);
		//p.z = 0.5*(zpos[k] + zpos[k + 1]);
		integer iP = 0;
		switch (iflag) {
		case TEMPERATURE:
			// закомментирован медленный участок кода.
			//inDomain = in_model_temp(p, ib, b, lb);
		    iP = evt[i_1];
			if (iP > 0) {
				ib = whot_is_block[iP-1];
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
					inDomain = false;
				}
				else inDomain = true;
			}
			else inDomain = false;			
			break;
		case HYDRODINAMIC: 
			// закомментирован медленный участок кода.
			//inDomain = in_model_flow(p, ib, b, lb);
			iP = evt[i_1];
			if (iP > 0) {
				ib = whot_is_block[iP - 1];
				if (b[ib].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
					inDomain = false;
				}
				else inDomain = true;
			}
			else inDomain = false;
			break;
		}

		if (inDomain) {
			TOCHKA pnew;
			// контрольный объём принадлежит модели
			pnew.x = xpos[i];
			pnew.y = ypos[j];
			pnew.z = zpos[k];
			node = i + j_2;
			// На вызов функции тратится немного времени. Сделаем операции явно без вызова функции.
			//addpoint(pa, maxnode, pnew, ent, node, bvisit);
			if (!bvisit[node]) {
				// добавление, точка ещё не встречалась

				pa[maxnode++] = pnew;
				ent[node] = maxnode; // нумерация начинается с единицы
				bvisit[node] = true; // индексация хеша.

			}

			pnew.x = xpos[i + 1];
			pnew.y = ypos[j];
			pnew.z = zpos[k];
			node = (i + 1) + j_2;
			//addpoint(pa, maxnode, pnew, ent, node, bvisit);
			if (!bvisit[node]) {
				// добавление, точка ещё не встречалась

				pa[maxnode++] = pnew;
				ent[node] = maxnode; // нумерация начинается с единицы
				bvisit[node] = true; // индексация хеша.

			}

			pnew.x = xpos[i];
			pnew.y = ypos[j + 1];
			pnew.z = zpos[k];
			node = i + j_3;
			//addpoint(pa, maxnode, pnew, ent, node, bvisit);
			if (!bvisit[node]) {
				// добавление, точка ещё не встречалась

				pa[maxnode++] = pnew;
				ent[node] = maxnode; // нумерация начинается с единицы
				bvisit[node] = true; // индексация хеша.

			}

			pnew.x = xpos[i + 1];
			pnew.y = ypos[j + 1];
			pnew.z = zpos[k];
			node = (i + 1) + j_3;
			//addpoint(pa, maxnode, pnew, ent, node, bvisit);
			if (!bvisit[node]) {
				// добавление, точка ещё не встречалась

				pa[maxnode++] = pnew;
				ent[node] = maxnode; // нумерация начинается с единицы
				bvisit[node] = true; // индексация хеша.

			}

			pnew.x = xpos[i];
			pnew.y = ypos[j];
			pnew.z = zpos[k + 1];
			node = i + j_4;
			//addpoint(pa, maxnode, pnew, ent, node, bvisit);
			if (!bvisit[node]) {
				// добавление, точка ещё не встречалась

				pa[maxnode++] = pnew;
				ent[node] = maxnode; // нумерация начинается с единицы
				bvisit[node] = true; // индексация хеша.

			}

			pnew.x = xpos[i + 1];
			pnew.y = ypos[j];
			pnew.z = zpos[k + 1];
			node = (i + 1) + j_4;
			//addpoint(pa, maxnode, pnew, ent, node, bvisit);
			if (!bvisit[node]) {
				// добавление, точка ещё не встречалась

				pa[maxnode++] = pnew;
				ent[node] = maxnode; // нумерация начинается с единицы
				bvisit[node] = true; // индексация хеша.

			}

			pnew.x = xpos[i];
			pnew.y = ypos[j + 1];
			pnew.z = zpos[k + 1];
			node = i + j_5;
			//addpoint(pa, maxnode, pnew, ent, node, bvisit);
			if (!bvisit[node]) {
				// добавление, точка ещё не встречалась

				pa[maxnode++] = pnew;
				ent[node] = maxnode; // нумерация начинается с единицы
				bvisit[node] = true; // индексация хеша.

			}

			pnew.x = xpos[i + 1];
			pnew.y = ypos[j + 1];
			pnew.z = zpos[k + 1];
			node = (i + 1) + j_5;
			//addpoint(pa, maxnode, pnew, ent, node, bvisit);
			if (!bvisit[node]) {
				// добавление, точка ещё не встречалась

				pa[maxnode++] = pnew;
				ent[node] = maxnode; // нумерация начинается с единицы
				bvisit[node] = true; // индексация хеша.

			}

		}
	}

	// Освобождение памяти 
	// уменьшение размера массива
	SetLength_point(pa, (inx + 1)*(iny + 1)*(inz + 1), maxnode);
	
	delete[] bvisit;
	bvisit = nullptr;
	
} // constr_nodes


  // создаёт массив узлов принадлежащих расчётной области
void constr_nodes_not_optimaze(TOCHKA* &pa, integer &maxnode, integer* &ent,
	integer iflag, int* &whot_is_block, int* &evt,
	integer inx, integer iny, integer inz,
	doublereal *xpos, doublereal *ypos, doublereal *zpos,
	BLOCK* b, integer lb, TOCKA_SHORT_INT* &tck_int_list, integer &maxelm) {
	// iflag - принимает два значения: TEMPERATURE или HYDRODINAMIC и 
	// указывает с какие уравнения будут решаться и в какой расчётной области.

	// ent - глобальная нумерация узлов.
	// нумерация узлов начинается с единицы, если узел не принадлежит
	// расчётной области то ставится значение 0.
	const integer iSIZE = (inx + 1)*(iny + 1)*(inz + 1);

	bool *in_dom = new bool[iSIZE];

	ent = nullptr;
	ent = new integer[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (ent == nullptr) {
	// недостаточно памяти на данном оборудовании.
	//printf("Problem: not enough memory on your equipment for ent constr struct...\n");
	//printf("Please any key to exit...\n");
	//system("PAUSE");
	//system("pause");
	//exit(1);
	//}
	// был ли добавлен узел с номером node
	bool *bvisit = nullptr;
	bvisit = new bool[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (bvisit == nullptr) {
	// недостаточно памяти на данном оборудовании.
	//printf("Problem: not enough memory on your equipment for bvisit constr struct...\n");
	//printf("Please any key to exit...\n");
	//system("PAUSE");
	//system("pause");
	//exit(1);
	//}
	// Выделение памяти. Сразу много для наискорейшего выполнения:
	pa = nullptr;
	pa = new TOCHKA[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (pa == nullptr) {
	// недостаточно памяти на данном оборудовании.
	//printf("Problem: not enough memory on your equipment for pa constr struct...\n");
	//printf("Please any key to exit...\n");
	//system("PAUSE");
	//system("pause");
	//exit(1);
	//}


	//TOCHKA p;
#pragma omp parallel for
	for (integer i = 0; i<iSIZE; ++i) {
		ent[i] = 0;
		bvisit[i] = false;
		in_dom[i] = false;
	}

	bool inDomain = false;
	integer ib, node;
	maxnode = 0;// maxnode - текущая длина массива pa
				// Цикл по всем контрольным объёмам
				//for (k = 0; k<inz; ++k) for (j = 0; j<iny; ++j)  for (i = 0; i<inx; ++i) {


	for (integer iscan = 0; iscan < maxelm; iscan++) {

		int i = tck_int_list[iscan].i;
		int j = tck_int_list[iscan].j;
		int k = tck_int_list[iscan].k;

		integer k_1 = k*inx*iny;
		integer j_1 = k_1 + j*inx;
		integer i_1 = i + j_1;

		integer iP = 0;
		switch (iflag) {
		case TEMPERATURE:
			// закомментирован медленный участок кода.
			//inDomain = in_model_temp(p, ib, b, lb);
			iP = evt[i_1];
			if (iP > 0) {
				ib = whot_is_block[iP - 1];
				if (b[ib].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
					inDomain = false;
				}
				else inDomain = true;
			}
			else inDomain = false;
			break;
		case HYDRODINAMIC:
			// закомментирован медленный участок кода.
			//inDomain = in_model_flow(p, ib, b, lb);
			iP = evt[i_1];
			if (iP > 0) {
				ib = whot_is_block[iP - 1];
				if (b[ib].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
					inDomain = false;
				}
				else inDomain = true;
			}
			else inDomain = false;
			break;
		}

		if (inDomain) {
			in_dom[i_1] = true;
		}
	}

	
		for (integer iscan = 0; iscan < maxelm; iscan++) {

			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;			

			integer k_1 = k*inx*iny;
			integer j_1 = k_1 + j*inx;
			integer i_1 = i + j_1;

			if (in_dom[i_1]) {
				TOCHKA pnew;
				
					// контрольный объём принадлежит модели
					pnew.x = xpos[i];
					pnew.y = ypos[j];
					pnew.z = zpos[k];
					node = i + j*(inx + 1)+ k*(inx + 1)*(iny + 1);
					//addpoint(pa, maxnode, pnew, ent, node, bvisit);
					// Раскрываем вызов функции. Есть гипотеза что на вызов функции 
					// тратися время.
					if (!bvisit[node]) {
						// добавление, точка ещё не встречалась

						pa[maxnode++] = pnew;
						ent[node] = maxnode; // нумерация начинается с единицы
						bvisit[node] = true; // индексация хеша.

					}

			}
		}
	

		for (integer iscan = 0; iscan < maxelm; iscan++) {

			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;

			integer k_1 = k*inx*iny;
			integer j_1 = k_1 + j*inx;
			integer i_1 = i + j_1;

			if (in_dom[i_1]) {
				TOCHKA pnew;
				pnew.x = xpos[i + 1];
				pnew.y = ypos[j];
				pnew.z = zpos[k];
				node = (i + 1) + j*(inx + 1) + k*(inx + 1)*(iny + 1);
				//addpoint(pa, maxnode, pnew, ent, node, bvisit);
				if (!bvisit[node]) {
					// добавление, точка ещё не встречалась

					pa[maxnode++] = pnew;
					ent[node] = maxnode; // нумерация начинается с единицы
					bvisit[node] = true; // индексация хеша.

				}
			}
		}

		for (integer iscan = 0; iscan < maxelm; iscan++) {

			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;

			integer k_1 = k*inx*iny;
			integer j_1 = k_1 + j*inx;
			integer i_1 = i + j_1;

			if (in_dom[i_1]) {
				TOCHKA pnew;
				pnew.x = xpos[i];
				pnew.y = ypos[j + 1];
				pnew.z = zpos[k];
				node = i + (j + 1)*(inx + 1) + k*(inx + 1)*(iny + 1);
				//addpoint(pa, maxnode, pnew, ent, node, bvisit);
				if (!bvisit[node]) {
					// добавление, точка ещё не встречалась

					pa[maxnode++] = pnew;
					ent[node] = maxnode; // нумерация начинается с единицы
					bvisit[node] = true; // индексация хеша.

				}
			}
		}

		for (integer iscan = 0; iscan < maxelm; iscan++) {

			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;

			integer k_1 = k*inx*iny;
			integer j_1 = k_1 + j*inx;
			integer i_1 = i + j_1;

			if (in_dom[i_1]) {
				TOCHKA pnew;
				pnew.x = xpos[i + 1];
				pnew.y = ypos[j + 1];
				pnew.z = zpos[k];
				node = (i + 1) + (j + 1)*(inx + 1) + k*(inx + 1)*(iny + 1);
				//addpoint(pa, maxnode, pnew, ent, node, bvisit);
				if (!bvisit[node]) {
					// добавление, точка ещё не встречалась

					pa[maxnode++] = pnew;
					ent[node] = maxnode; // нумерация начинается с единицы
					bvisit[node] = true; // индексация хеша.

				}
			}
		}

		for (integer iscan = 0; iscan < maxelm; iscan++) {

			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;

			integer k_1 = k*inx*iny;
			integer j_1 = k_1 + j*inx;
			integer i_1 = i + j_1;

			if (in_dom[i_1]) {
				TOCHKA pnew;
				pnew.x = xpos[i];
				pnew.y = ypos[j];
				pnew.z = zpos[k + 1];
				node = i + j*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
				//addpoint(pa, maxnode, pnew, ent, node, bvisit);
				if (!bvisit[node]) {
					// добавление, точка ещё не встречалась

					pa[maxnode++] = pnew;
					ent[node] = maxnode; // нумерация начинается с единицы
					bvisit[node] = true; // индексация хеша.

				}
			}
		}

		for (integer iscan = 0; iscan < maxelm; iscan++) {

			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;

			integer k_1 = k*inx*iny;
			integer j_1 = k_1 + j*inx;
			integer i_1 = i + j_1;

			if (in_dom[i_1]) {
				TOCHKA pnew;
				pnew.x = xpos[i + 1];
				pnew.y = ypos[j];
				pnew.z = zpos[k + 1];
				node = (i + 1) + j*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
				//addpoint(pa, maxnode, pnew, ent, node, bvisit);
				if (!bvisit[node]) {
					// добавление, точка ещё не встречалась

					pa[maxnode++] = pnew;
					ent[node] = maxnode; // нумерация начинается с единицы
					bvisit[node] = true; // индексация хеша.

				}
			}
		}

		for (integer iscan = 0; iscan < maxelm; iscan++) {

			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;

			integer k_1 = k*inx*iny;
			integer j_1 = k_1 + j*inx;
			integer i_1 = i + j_1;

			if (in_dom[i_1]) {
				TOCHKA pnew;
				pnew.x = xpos[i];
				pnew.y = ypos[j + 1];
				pnew.z = zpos[k + 1];
				node = i + (j + 1)*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
				//addpoint(pa, maxnode, pnew, ent, node, bvisit);
				if (!bvisit[node]) {
					// добавление, точка ещё не встречалась

					pa[maxnode++] = pnew;
					ent[node] = maxnode; // нумерация начинается с единицы
					bvisit[node] = true; // индексация хеша.

				}
			}
		}

		for (integer iscan = 0; iscan < maxelm; iscan++) {

			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;

			integer k_1 = k*inx*iny;
			integer j_1 = k_1 + j*inx;
			integer i_1 = i + j_1;

			if (in_dom[i_1]) {
				TOCHKA pnew;
				pnew.x = xpos[i + 1];
				pnew.y = ypos[j + 1];
				pnew.z = zpos[k + 1];
				node = (i + 1) + (j + 1)*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
				//addpoint(pa, maxnode, pnew, ent, node, bvisit);
				if (!bvisit[node]) {
					// добавление, точка ещё не встречалась

					pa[maxnode++] = pnew;
					ent[node] = maxnode; // нумерация начинается с единицы
					bvisit[node] = true; // индексация хеша.

				}
			}
		}

	delete[] in_dom;

	// Освобождение памяти 
	// уменьшение размера массива
	SetLength_point(pa, (inx + 1)*(iny + 1)*(inz + 1), maxnode);
	if (bvisit != nullptr) {
		delete[] bvisit;
		bvisit = nullptr;
	}
} // constr_nodes_not_optimaze


// создаёт массив узлов принадлежащих расчётной области
void constr_nodes_flow(TOCHKA* &pa, int &maxnode, int* &ent,
	int **evt_f2, integer iDom,
	int inx, int iny, int inz,
	doublereal *xpos, doublereal *ypos, doublereal *zpos) {

	const integer iSIZE = (inx + 1)*(iny + 1)*(inz + 1);

	// ent - глобальная нумерация узлов.
	// нумерация узлов начинается с единицы, если узел не принадлежит
	// расчётной области то ставится значение 0.
	//ent = nullptr;
	ent = new int[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (ent == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for ent flow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// был ли добавлен узел с номером node
	bool *bvisit = nullptr;
	bvisit = new bool[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (bvisit == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for bvisit flow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// Выделение памяти. Сразу много для наискорейшего выполнения:
	//pa = nullptr;
	pa = new TOCHKA[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (pa == nullptr) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for pa flow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}


	// Инициализация:
	integer i=0, j=0, k=0;
	for (i = 0; i<iSIZE; ++i) {
		ent[i] = 0;
		bvisit[i] = false;
	}

	bool inDomain = false;
	integer node, ivol;
	TOCHKA pnew;
	maxnode = 0;// maxnode - текущая длина массива pa
	// Цикл по всем контрольным объёмам
	for (k = 0; k<inz; ++k) for (j = 0; j<iny; ++j) for (i = 0; i<inx; ++i) {

		ivol = i + j*inx + k*inx*iny; // номер контрольного объёма
		if (evt_f2[1][ivol] == (iDom + 1)) inDomain = true;

		if (inDomain) {
			// контрольный объём принадлежит модели
			pnew.x = xpos[i];
			pnew.y = ypos[j];
			pnew.z = zpos[k];
			node = i + j*(inx + 1) + k*(inx + 1)*(iny + 1);
			addpoint(pa, maxnode, pnew, ent, node, bvisit);
			pnew.x = xpos[i + 1];
			pnew.y = ypos[j];
			pnew.z = zpos[k];
			node = (i + 1) + j*(inx + 1) + k*(inx + 1)*(iny + 1);
			addpoint(pa, maxnode, pnew, ent, node, bvisit);
			pnew.x = xpos[i];
			pnew.y = ypos[j + 1];
			pnew.z = zpos[k];
			node = i + (j + 1)*(inx + 1) + k*(inx + 1)*(iny + 1);
			addpoint(pa, maxnode, pnew, ent, node, bvisit);
			pnew.x = xpos[i + 1];
			pnew.y = ypos[j + 1];
			pnew.z = zpos[k];
			node = (i + 1) + (j + 1)*(inx + 1) + k*(inx + 1)*(iny + 1);
			addpoint(pa, maxnode, pnew, ent, node, bvisit);
			pnew.x = xpos[i];
			pnew.y = ypos[j];
			pnew.z = zpos[k + 1];
			node = i + j*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			addpoint(pa, maxnode, pnew, ent, node, bvisit);
			pnew.x = xpos[i + 1];
			pnew.y = ypos[j];
			pnew.z = zpos[k + 1];
			node = (i + 1) + j*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			addpoint(pa, maxnode, pnew, ent, node, bvisit);
			pnew.x = xpos[i];
			pnew.y = ypos[j + 1];
			pnew.z = zpos[k + 1];
			node = i + (j + 1)*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			addpoint(pa, maxnode, pnew, ent, node, bvisit);
			pnew.x = xpos[i + 1];
			pnew.y = ypos[j + 1];
			pnew.z = zpos[k + 1];
			node = (i + 1) + (j + 1)*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			addpoint(pa, maxnode, pnew, ent, node, bvisit);
		}
	}

	// Освобождение памяти 
	// уменьшение размера массива
	SetLength_point(pa, iSIZE, maxnode);
	//if (bvisit != nullptr) {
	// оператор delete может применяться повторно в том числе и к null указателю.
		delete[] bvisit;
		bvisit = nullptr;
	//}
} // constr_nodes_flow



// для каждого контрольного объёма принадлежащему
// расчётной области определяет номера его вершин.
void constr_nvtx(int* evt, int* ent, int** &nvtx, int &maxelm,
	int inx, int iny, int inz, TOCKA_SHORT_INT* &tck_int_list) {

	// нумерация ent начинается с единицы.

	// проити по всем контрольным объёмам принадлежащим
	// расчётной области
	// maxelm - число контрольных объёмов принадлежащих расчётной области
	//maxelm = 0;
	integer l = 0;
	
	// подсчёт количества контрольных объёмов 
	// принадлежащих расчётной области.
	/*
	for (i = 0; i<(inx); ++i) for (j = 0; j<(iny); ++j) for (k = 0; k<(inz); ++k) {
		ic = i + j*inx + k*inx*iny;
		if (evt[ic] > 0) {
			maxelm++;
		}
	}
	*/
	nvtx = nullptr;
	nvtx = new int*[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
	if (nvtx == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for nvtx constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer i = 0; i<NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) nvtx[i] = nullptr;
	for (integer i = 0; i<NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) {
		nvtx[i] = new int[maxelm];
		if (nvtx[i] == nullptr) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem: not enough memory on your equipment for nvtx[%lld] constr struct...\n", i);
#else
			printf("Problem: not enough memory on your equipment for nvtx[%d] constr struct...\n", i);
#endif
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// Обход по всем контрольным объёмам.
	//for (i = 0; i<(inx); ++i) for (j = 0; j<(iny); ++j) for (k = 0; k<(inz); ++k) {
		//ic = i + j*inx + k*inx*iny;
		//if (evt[ic]>0) {
	// integer i,j,k;
   //  integer iP_k, iP_j, iP_k_p1, iP_j_p1;
	for (integer iscan=0; iscan<maxelm; iscan++) {
		{
			 int i = tck_int_list[iscan].i;
			 int j = tck_int_list[iscan].j;
			 int  k = tck_int_list[iscan].k;
			// integer ic = i + j * inx + k * inx*iny;

			 // контрольный объём принадлежит расчётной области.
			 integer iP_k=k*(inx + 1)*(iny + 1);
			 integer iP_j=j*(inx + 1);
             integer iP_k_p1=(k+1)*(inx + 1)*(iny + 1);
             integer iP_j_p1=(j+1)*(inx + 1);
			 //integer i1, i2, i3, i4, i5, i6, i7, i8;

			integer i1 = i + iP_j + iP_k;
			integer i2 = (i + 1) + iP_j + iP_k;
			integer i3 = i + iP_j_p1 + iP_k;
			integer i4 = (i + 1) + iP_j_p1 + iP_k;
			integer i5 = i + iP_j + iP_k_p1;
			integer i6 = (i + 1) + iP_j + iP_k_p1;
			integer i7 = i + iP_j_p1 + iP_k_p1;
			integer i8 = (i + 1) + iP_j_p1 + iP_k_p1;
			// Добавка 15 апреля 2017.
			// ent - нумерация начинается с единицы.
			if ((ent[i1] > 0) && (ent[i2] > 0) && (ent[i3] > 0) && (ent[i4] > 0) && (ent[i5] > 0) && (ent[i6] > 0) && (ent[i7] > 0) && (ent[i8] > 0)) {
				nvtx[0][l] = ent[i1];
				nvtx[1][l] = ent[i2];
				nvtx[2][l] = ent[i3];
				nvtx[3][l] = ent[i4];
				nvtx[4][l] = ent[i5];
				nvtx[5][l] = ent[i6];
				nvtx[6][l] = ent[i7];
				nvtx[7][l] = ent[i8];

				l++;
				/*
				// Проверка сделанная выше показывает что это условие никогда не выполняется. 11.01.2020
				if ((ent[i1] == 0) || (ent[i2] == 0) || (ent[i3] == 0) || (ent[i4] == 0) || (ent[i5] == 0) || (ent[i6] == 0) || (ent[i7] == 0) || (ent[i8] == 0)) {
#if doubleintprecision == 1
					printf("error: ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#else
					printf("error: ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#endif
					//system("PAUSE");
					system("pause");
				}
				*/
			}
			else {
#if doubleintprecision == 1
				printf("error TEMP ent <= 0: ent[%lld]=%d ent[%lld]=%d ent[%lld]=%d ent[%lld]=%d ent[%lld]=%d ent[%lld]=%d ent[%lld]=%d ent[%lld]=%d \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#else
				printf("error TEMP ent <= 0: ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#endif
				//system("PAUSE");
				system("pause");
			}
			
			
		}
	}
} // constr_nvtx

// Функция с быстрым временем выполнения.
// 9.08.2017.
void walk_in_octree_icolor_different_fluid_domain(octree* &oc,
	int inx, int iny, int inz,
	integer* &icolor_different_fluid_domain,
	int** evt_f2, integer maxelm_flow, 
	TOCKA_SHORT_INT*& tck_int_list) {

	//integer idiagnostic_error_coloc_counter1 = 0, idiagnostic_error_coloc_counter2 = 0;
	top_ALICE_STACK = 0;
	if (oc->link0 != nullptr) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != nullptr) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != nullptr) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != nullptr) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != nullptr) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != nullptr) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != nullptr) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != nullptr) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link7);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link7->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link7->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link7->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link7->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link7->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link7->minz;
		top_ALICE_STACK++;
	}
	while (top_ALICE_STACK > 0) {
		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != nullptr) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist  ) {



				octree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;

				// разбиение на 8.
				//int minx = my_ALICE_STACK[top_ALICE_STACK - 1].minx;
				//int maxx = my_ALICE_STACK[top_ALICE_STACK - 1].maxx;
				//int miny = my_ALICE_STACK[top_ALICE_STACK - 1].miny;
				//int maxy = my_ALICE_STACK[top_ALICE_STACK - 1].maxy;
				//int minz = my_ALICE_STACK[top_ALICE_STACK - 1].minz;
				//int maxz = my_ALICE_STACK[top_ALICE_STACK - 1].maxz;

				// Дробление  вызывается.
				my_ALICE_STACK[top_ALICE_STACK - 1].link = nullptr;
				top_ALICE_STACK--;

				// Обработка узла.
				if (octree1->inum_FD > 0) {
					integer ielm = octree1->inum_FD - 1; // номер nvtx.
					
					//integer ic = 0;
					//ic = minx + miny*inx + minz*inx*iny;
					//12.01.2020
					//ic = tck_int_list[ielm].i + tck_int_list[ielm].j * inx + tck_int_list[ielm].k * inx * iny;
					integer id_found = ielm;

					icolor_different_fluid_domain[id_found] = 1; // Гидродинамическая ячейка. Присваиваем единицу.
						
				}

			}
			else {
				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != nullptr) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != nullptr) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != nullptr) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != nullptr) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != nullptr) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != nullptr) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != nullptr) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != nullptr) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link7);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link7->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link7->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link7->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link7->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link7->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link7->minz;
					top_ALICE_STACK++;
				}
			}
		}
		//}
		//system("PAUSE");
	}

} // walk_in_octree_icolor_different_fluid_domain

// 22 сентября 2016 вычисление функции цвета для АЛИС сетки.
void constr_icolor_different_fluid_domain_alice(integer maxelm_flow, integer*& icolor_different_fluid_domain,
	int inx, int iny, int inz, int** evt_f2, doublereal*& x_pos, doublereal*& y_pos, doublereal*& z_pos, TOCHKA* pa,
	int** nvtx, octree*& oc,
	TOCKA_SHORT_INT*& tck_int_list, BLOCK* b, integer lb)
{
	

	//printf("maxelm_flow=%lld maxelm=%lld\n", maxelm_flow, maxelm);// они равны.

	// Выделение памяти.
	icolor_different_fluid_domain = nullptr;
	icolor_different_fluid_domain = new integer[maxelm_flow];
	if (icolor_different_fluid_domain == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for icolor_different_fluid_domain constr_nvtx_flow in constr_struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		exit(1);
	}

#pragma omp parallel for 
	for (integer i = 0; i < maxelm_flow; ++i) {
		// инициализация не существующим значением.
		// Нумерация гидродинамических областей начинается с единицы.
		icolor_different_fluid_domain[i] = -1;// заглушка.
	}
	

	if (print_log_message) {
		std::cout << "maxelm_flow=" << maxelm_flow << " inx=" << inx << ", iny=" << iny << " inz=" << inz << std::endl;
	}

	/*
	// УСТАРЕЛО. Полностью заменено на функцию walk_in_octree_icolor_different_fluid_domain.
	// Чрезвычайно медленный кусок кода. Устарело 13.01.2020. Оставить для истории закомментированным.
	//for (integer i = 0; i<(inx); ++i) for (integer j = 0; j<(iny); ++j) for (integer k = 0; k<(inz); ++k) {
	for (integer iscan = 0; iscan < maxelm_flow; iscan++)
	{
		i = tck_int_list[iscan].i;
		j = tck_int_list[iscan].j;
		k = tck_int_list[iscan].k;
		integer ic = i + j * inx + k * inx * iny;//iP
		
		//integer ic = 0;
		//ic = i + j*inx + k*inx*iny;
		// Координаты центра первичной ячейки.
		doublereal x_c = 0.5*(x_pos[i]+x_pos[i+1]);
		doublereal y_c = 0.5*(y_pos[j] + y_pos[j + 1]);
		doublereal z_c = 0.5*(z_pos[k] + z_pos[k + 1]);

		/*
		TOCHKA p_23;
		p_23.x = x_c; p_23.y = y_c; p_23.z = z_c;
		integer ib_23;
		// Проверка показала что они все гидродинамические.
		if (!in_model_flow(p_23, ib_23, b, lb)) {
			printf("it's not a liquid cell\n");
			system("pause");
		}
		*/
	    /*
		// Поиск нужной ic64.
		bool bfound = false;
		integer id_found = -1;
		//begin
		//printf("begin "); //debug
		for (integer i_1 = 0; i_1 < maxelm_flow; ++i_1) {
			if (pa[nvtx[0][i_1] - 1].x>=pa[nvtx[1][i_1] - 1].x) {
				printf("fatall error! ");
				printf("maxelm_flow=%d\n", maxelm_flow);
				printf("pa[nvtx[0][i_1] - 1].x=%e pa[nvtx[1][i_1] - 1].x=%e\n", pa[nvtx[0][i_1] - 1].x, pa[nvtx[1][i_1] - 1].x);
				printf("%d %d %d %d %d %d %d %d\n", nvtx[0][i_1], nvtx[1][i_1], nvtx[2][i_1], nvtx[3][i_1], nvtx[4][i_1], nvtx[5][i_1], nvtx[6][i_1], nvtx[7][i_1]);
				system("pause");
			}
			if (pa[nvtx[0][i_1] - 1].y >= pa[nvtx[3][i_1] - 1].y) {
				printf("fatall error! ");
				printf("maxelm_flow=%d\n", maxelm_flow);
				printf("pa[nvtx[0][i_1] - 1].y=%e pa[nvtx[3][i_1] - 1].y=%e\n", pa[nvtx[0][i_1] - 1].y, pa[nvtx[3][i_1] - 1].y);
				printf("%d %d %d %d %d %d %d %d\n", nvtx[0][i_1], nvtx[1][i_1], nvtx[2][i_1], nvtx[3][i_1], nvtx[4][i_1], nvtx[5][i_1], nvtx[6][i_1], nvtx[7][i_1]);
				system("pause");
			}
			if (pa[nvtx[0][i_1] - 1].z >= pa[nvtx[4][i_1] - 1].z) {
				printf("fatall error! ");
				printf("maxelm_flow=%d\n", maxelm_flow);
				printf("pa[nvtx[0][i_1] - 1].z=%e pa[nvtx[4][i_1] - 1].z=%e\n", pa[nvtx[0][i_1] - 1].z, pa[nvtx[4][i_1] - 1].z);
				printf("%d %d %d %d %d %d %d %d\n", nvtx[0][i_1], nvtx[1][i_1], nvtx[2][i_1], nvtx[3][i_1], nvtx[4][i_1], nvtx[5][i_1], nvtx[6][i_1], nvtx[7][i_1]);
				system("pause");
			}

			if ((x_c>=pa[nvtx[0][i_1] - 1].x) && (x_c<=pa[nvtx[1][i_1] - 1].x) && 
			   (y_c>=pa[nvtx[0][i_1] - 1].y) && (y_c<=pa[nvtx[3][i_1] - 1].y) &&
			   (z_c>=pa[nvtx[0][i_1] - 1].z) && (z_c<=pa[nvtx[4][i_1] - 1].z))
			{
				bfound = true;
				id_found = i_1;
				//printf("id_found=%lld %lld",i_1, evt_f2[MASKDOMAINFLUIDCOLOR][ic]);// debug
				// В случае нахождения досрочный выход из цикла for.
				break;
			}
		}
		//printf("end\n");
		//system("pause");
		if (bfound) {
			// Находит одно и тоже значение id_found несколько раз.
			icolor_different_fluid_domain[id_found] = evt_f2[MASKDOMAINFLUIDCOLOR][ic]; // Присвоение цвета который хранится в третьем поле evt_f2.
			if (evt_f2[MASKDOMAINFLUIDCOLOR][ic] <= 0) {
				printf("evt_f2[MASKDOMAINFLUIDCOLOR][%lld]=%lld\n",ic, evt_f2[MASKDOMAINFLUIDCOLOR][ic]);
				system("pause");
			}
		}
	}
	*/
	// 9.08.2017 он заменён на более быстрый работающий за один проход.
	walk_in_octree_icolor_different_fluid_domain(oc, inx, iny, inz, icolor_different_fluid_domain, evt_f2, maxelm_flow, tck_int_list);

	// инвариантная проверка - заполнены ли все жидкие ячейки.
	integer ic_err1 = 0, ic_err2 = 0;

#pragma omp parallel for reduction(+ : ic_err1, ic_err2) 
	for (integer i = 0; i < maxelm_flow; ++i) {
		// инициализация не существующим значением.
		// Нумерация гидродинамических областей начинается с единицы.
		if (icolor_different_fluid_domain[i] <= -1) {
			ic_err1++;
		}
		if (icolor_different_fluid_domain[i] == 0) {
			ic_err2++;
		}
	}
	if ((ic_err1 > 0)||(ic_err2>0)) {
		std::cout <<  "file constr_struct. function constr_icolor_different_fluid_domain" << std::endl;
		std::cout <<  "ERROR!!! icolor_different_fluid_domain[i] == -1 value is empty" << std::endl;
		std::cout <<  "number errors <=-1="<< ic_err1 <<" 0="<< ic_err2 <<" sum="<< ic_err1+ ic_err2 << std::endl;
		//system("PAUSE");
	}

} //  constr_icolor_different_fluid_domain_alice

// Вычисляет maxelm;
void calculate_max_elm(int &maxelm, int** evt_f2, 
	int inx, int iny, int inz, integer iDom) {

	maxelm = 0;
	int im = 0;
	// integer i, j, k;
	// integer ic, ic_k, ic_j;
	// подсчёт количества контрольных объёмов 
	// принадлежащих расчётной области.
#pragma omp parallel for reduction(+ : im)
	for (int k = 0; k < inz; ++k) 
	{
		integer  ic_k=k * inx * iny;
		for (int  j = 0; j < iny; ++j)
		{
			integer  ic_j=j * inx +ic_k;
			for (int i = 0; i < inx; ++i)
	        {
				integer ic = i +  ic_j;
				if ((evt_f2[ENUMERATECONTVOL][ic] > 0) && (evt_f2[MASKDOMAINFLUID][ic] == (iDom + 1))) {
					im++;
				}
			}
		}
	}

	maxelm = im;
}//calculate_max_elm


// для каждого контрольного объёма принадлежащему
// расчётной области определяет номера его вершин.
void constr_nvtx_flow(int** &evt_f2, integer* &icolor_different_fluid_domain,
	integer iDom, int* &ent,
	int** &nvtx, int maxelm,
	int inx, int iny, int inz) {

	// проити по всем контрольным объёмам принадлежащим
	// расчётной области
	// maxelm - число контрольных объёмов принадлежащих расчётной области
	//maxelm = 0;
	//integer i, j, k, ic;
	integer l = 0;
	
	
	// Выделение памяти.
	icolor_different_fluid_domain = nullptr;
	icolor_different_fluid_domain = new integer[maxelm];
	if (icolor_different_fluid_domain == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment ";
		std::cout << "for icolor_different_fluid_domain constr_nvtx_flow";
		std::cout << "in constr_struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		exit(1);
	}

#pragma omp parallel for
	for (integer i = 0; i < maxelm; ++i) {
		// Для структурированной сетки.
		icolor_different_fluid_domain[i] = -1;// заглушка инициализация не существующим значением.
	}
	// integer ic_k, ic_j;

	integer ic64 = 0;

#pragma omp parallel for reduction(+ : ic64)	
	 for (int  k = 0; k<inz; ++k)
	 {
		integer  ic_k=k*inx*iny;
		for (int  j = 0; j<iny; ++j)
		{
            integer ic_j=j*inx +ic_k;
			for (int  i = 0; i<inx; ++i)
			{

				integer  ic = i +  ic_j;
				if ((evt_f2[ENUMERATECONTVOL][ic] > 0) && (evt_f2[MASKDOMAINFLUID][ic] == (iDom + 1))) {
		
					ic64++;
				}
			}
		}
	}

#pragma omp parallel for
	 for (integer i = 0; i < ic64; ++i) {
		 // Для структурированной сетки.
		 // Это точно активная жидкая ячейка. У нас только одна слитая воедино гидродинамическая подобласть, поэтому просто присвоим 1
		 icolor_different_fluid_domain[i] = 1;
	 }

	nvtx = nullptr;
	nvtx = new int*[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
	if (nvtx == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for nvtx flow constr struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer  i = 0; i<NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) nvtx[i] = nullptr;
	for (integer   i = 0; i<NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) {
		nvtx[i] = new int[maxelm];
		if (nvtx[i] == nullptr) {
			// недостаточно памяти на данном оборудовании.
			std::cout <<  "Problem: not enough memory on your equipment for nvtx[" << i << "] constr struct..." << std::endl; 
			
			std::cout <<  "Please any key to exit..." << std::endl;
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// integer ic_k1, ic_k_p1, ic_j1, ic_j_P1;
	// Обход по всем контрольным объёмам.
	for (int   k = 0; k<(inz); ++k)
	 {
		integer   ic_k=k*inx*iny;
		integer  ic_k1=k*(inx+1)*(iny+1);
		integer  ic_k_p1=(k+1)*(inx+1)*(iny+1);
		for (int  j = 0; j<(iny); ++j)
		{
            integer  ic_j=j*inx +ic_k;
			integer  ic_j1=j*(inx+1);
			integer  ic_j_P1=(j+1)*(inx+1);
			for (int   i = 0; i<(inx); ++i)
			{

				integer  ic = i +  ic_j;
				if ((evt_f2[ENUMERATECONTVOL][ic]>0) && (evt_f2[MASKDOMAINFLUID][ic] == (iDom + 1))) {
					// контрольный объём принадлежит расчётной области.
					//integer i1, i2, i3, i4, i5, i6, i7, i8;
					integer  i1 = i + ic_j1 + ic_k1;
					integer  i2 = (i + 1) + ic_j1 + ic_k1;
					integer  i3 = i + ic_j_P1 + ic_k1;
					integer  i4 = (i + 1) + ic_j_P1 + ic_k1;
					integer  i5 = i + ic_j1 + ic_k_p1;
					integer  i6 = (i + 1) + ic_j1 + ic_k_p1;
					integer  i7 = i + ic_j_P1 + ic_k_p1;
					integer  i8 = (i + 1) + ic_j_P1 + ic_k_p1;

			

					// Добавка 15 апреля 2017.
					// ent - нумерация начинается с единицы.
					if ((ent[i1] > 0) && (ent[i2] > 0) && (ent[i3] > 0) && 
						(ent[i4] > 0) && (ent[i5] > 0) && (ent[i6] > 0) && 
						(ent[i7] > 0) && (ent[i8] > 0)) {


						nvtx[0][l] = ent[i1];
						nvtx[1][l] = ent[i2];
						nvtx[2][l] = ent[i3];
						nvtx[3][l] = ent[i4];
						nvtx[4][l] = ent[i5];
						nvtx[5][l] = ent[i6];
						nvtx[6][l] = ent[i7];
						nvtx[7][l] = ent[i8];
						l++;

						/*
						// Этого условия не может быть. Проверка сделана выше.
						if ((ent[i1] == 0) || (ent[i2] == 0) || (ent[i3] == 0) || (ent[i4] == 0) || (ent[i5] == 0) || (ent[i6] == 0) || (ent[i7] == 0) || (ent[i8] == 0)) {
#if doubleintprecision == 1
							printf("error: ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#else
							printf("error: ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#endif
							//system("PAUSE");
							system("pause");
						}
						*/
					}
					else {
						std::cout <<  "error FLOW ent <= 0: ent[" <<i1<<"]="<<ent[i1]<<" ent["<<i2<<"]="<<ent[i2]<<" ent["<<13<<"]=";
						std::cout << ent[i3] <<"ent["<<i4<<"]="<<ent[i4]<<" ent["<<i5<<"]="<<ent[i5]<<" ent["<<i6<<"]="<<ent[i6]<<" ent["<<i7<<"]="<<ent[i7]<<" ent["<<i8<<"]="<<ent[i8]<<" "<<std::endl;
						//system("PAUSE");
						system("pause");
					}

				}
			}
		}
	}
} // constr_nvtx_flow

// Заносит свойства материалов в структуру
void constr_prop(int* &evt, int* &whot_is_block, int* &ent, float** &prop, integer maxelm, integer iflag, BLOCK* b,
	integer lb, int inx, int iny, int inz, doublereal* &Sc, POWER_TIME_DEPEND* &ipower_time_depend,
				 doublereal* &xpos, doublereal* &ypos, doublereal* &zpos, TPROP* matlist, TOCKA_SHORT_INT* &tck_int_list,
	bool* &bActiveShearModule) {

	if (bActiveShearModule != nullptr) {
		delete[] bActiveShearModule;
		bActiveShearModule = nullptr;
	}
	bActiveShearModule = new bool[maxelm];

	
	prop=nullptr;
	prop=new float*[SIZE_PROPERTIES_ARRAY];
	if (prop==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for prop constr struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer i=0; i<SIZE_PROPERTIES_ARRAY; ++i) prop[i]=nullptr;
	for (integer i=0; i<SIZE_PROPERTIES_ARRAY; ++i) {
		if (steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
			if (i <= 5) {
				prop[i] = new float[maxelm];
				if (prop[i] == nullptr) {
					// недостаточно памяти на данном оборудовании.
					std::cout << "Problem: not enough memory on your equipment for prop[" << i << "] constr struct..." << std::endl;

					std::cout << "Please any key to exit..." << std::endl;
					//system("PAUSE");
					system("pause");
					exit(1);
				}
			}
		}
		else {
			prop[i] = new float[maxelm];
			if (prop[i] == nullptr) {
				// недостаточно памяти на данном оборудовании.
				std::cout << "Problem: not enough memory on your equipment for prop[" << i << "] constr struct..." << std::endl;

				std::cout << "Please any key to exit..." << std::endl;
				//system("PAUSE");
				system("pause");
				exit(1);
			}
		}
		
	}
	Sc=nullptr;
	Sc=new doublereal[maxelm];
	if (Sc==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for Sc constr struct..."<<std::endl;
		std::cout << "Please any key to exit..."<<std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	ipower_time_depend = nullptr;
	ipower_time_depend = new POWER_TIME_DEPEND[maxelm];
	if (ipower_time_depend == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for ipower_time_depend constr struct..."<<std::endl;
		std::cout << "Please any key to exit..."<<std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}

	
   
	integer l=0; 
	
	// Проход по всем контрольным объёмам в порядке обхода:
	//for (integer i=0; i<inx; ++i) for (integer j=0; j<iny; ++j) for (integer k=0; k<inz; ++k) {

		//integer ic = i + j*inx + k*inx*iny;
		//if (evt[ic] > 0) {
	for (integer iscan = 0; iscan<maxelm; iscan++) {
		{
			int i = tck_int_list[iscan].i;
			int j = tck_int_list[iscan].j;
			int k = tck_int_list[iscan].k;
			integer ic = i + j * inx + k * inx*iny;

			TOCHKA p;

				p.x = 0.5*(xpos[i] + xpos[i + 1]);
				p.y = 0.5*(ypos[j] + ypos[j + 1]);
				p.z = 0.5*(zpos[k] + zpos[k + 1]);
				/*
				switch (iflag) {
				  case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
				  case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
				}
				*/

				bool inDomain = false;
				integer ib;

				integer iP = 0;
				switch (iflag) {
				case TEMPERATURE:
					// закомментирован медленный участок кода.
					//inDomain = in_model_temp(p, ib, b, lb);
					iP = evt[ic];
					if (iP > 0) {
						ib = whot_is_block[iP - 1];
						if (b[ib].itype == PHYSICS_TYPE_IN_BODY::HOLLOW) {
							inDomain = false;
						}
						else inDomain = true;
					}
					else inDomain = false;
					break;
				case HYDRODINAMIC:
					// закомментирован медленный участок кода.
					//inDomain = in_model_flow(p, ib, b, lb);
					iP = evt[ic];
					if (iP > 0) {
						ib = whot_is_block[iP - 1];
						if (b[ib].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
							inDomain = false;
						}
						else inDomain = true;
					}
					else inDomain = false;
					break;
				}

				

				if (inDomain) {
					if (l < maxelm)
					{

						integer imat_id_loc = b[ib].imatid;
						switch (iflag) {
						case TEMPERATURE:  prop[RHO][l] = matlist[imat_id_loc].rho;
							//prop[HEAT_CAPACITY][l] = matlist[b[ib].imatid].cp;
							//prop[LAM][l] = matlist[b[ib].imatid].lam;
							
							if (matlist[imat_id_loc].n_cp == 1) {
								prop[HEAT_CAPACITY][l] = matlist[imat_id_loc].arr_cp[0];
							}
							else {
								prop[HEAT_CAPACITY][l] = get_cp(matlist[imat_id_loc].n_cp, matlist[imat_id_loc].temp_cp, matlist[imat_id_loc].arr_cp, 25.0);
							}
							prop[LAM][l] = get_lam(matlist[imat_id_loc].n_lam, matlist[imat_id_loc].temp_lam, matlist[imat_id_loc].arr_lam, 25.0);

							prop[MULT_LAM_X][l] = matlist[imat_id_loc].orthotropy_multiplyer_x;
							prop[MULT_LAM_Y][l] = matlist[imat_id_loc].orthotropy_multiplyer_y;
							prop[MULT_LAM_Z][l] = matlist[imat_id_loc].orthotropy_multiplyer_z;
							
							if ((steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY)&&
								(((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
									(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
									(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
									(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE)))) {

								if (matlist[imat_id_loc].n_beta_t_solid == 1) {
									prop[BETA_T_MECHANICAL][l] = matlist[imat_id_loc].arr_beta_t_solid[0];
								}
								else {
									prop[BETA_T_MECHANICAL][l] = get_beta_t_solid(matlist[imat_id_loc].n_beta_t_solid, matlist[imat_id_loc].temp_beta_t_solid, matlist[imat_id_loc].arr_beta_t_solid, 25.0);
								}

								if (matlist[imat_id_loc].n_Poisson_ratio == 1) {
									prop[POISSON_RATIO][l] = matlist[imat_id_loc].arr_Poisson_ratio[0];
								}
								else {
									prop[POISSON_RATIO][l] = get_Poisson_ratio(matlist[imat_id_loc].n_Poisson_ratio, matlist[imat_id_loc].temp_Poisson_ratio, matlist[imat_id_loc].arr_Poisson_ratio, 25.0);
								}

								if (matlist[imat_id_loc].n_YoungModule == 1) {
									prop[YOUNG_MODULE][l] = matlist[imat_id_loc].arr_Young_Module[0];
								}
								else {
									prop[YOUNG_MODULE][l] = get_Young_Module(matlist[imat_id_loc].n_YoungModule, matlist[imat_id_loc].temp_Young_Module, matlist[imat_id_loc].arr_Young_Module, 25.0);
								}
								prop[MULT_BETA_T_MECHANICAL_X][l] = matlist[imat_id_loc].orthotropy_multiplyer_x_beta_t_solid;
								prop[MULT_BETA_T_MECHANICAL_Y][l] = matlist[imat_id_loc].orthotropy_multiplyer_y_beta_t_solid;
								prop[MULT_BETA_T_MECHANICAL_Z][l] = matlist[imat_id_loc].orthotropy_multiplyer_z_beta_t_solid;
								prop[MULT_YOUNG_MODULE_X][l] = matlist[imat_id_loc].orthotropy_multiplyer_x_Young_Module;
								prop[MULT_YOUNG_MODULE_Y][l] = matlist[imat_id_loc].orthotropy_multiplyer_y_Young_Module;
								prop[MULT_YOUNG_MODULE_Z][l] = matlist[imat_id_loc].orthotropy_multiplyer_z_Young_Module;
								
								// множитель коэффициента Пуассона.
								prop[MULT_POISSON_RATIO_YZ][l] = matlist[imat_id_loc].orthotropy_multiplyer_Poisson_ratio_yz;
								prop[MULT_POISSON_RATIO_XZ][l] = matlist[imat_id_loc].orthotropy_multiplyer_Poisson_ratio_xz;
								prop[MULT_POISSON_RATIO_XY][l] = matlist[imat_id_loc].orthotropy_multiplyer_Poisson_ratio_xy;
								prop[MULT_POISSON_RATIO_ZY][l] = matlist[imat_id_loc].orthotropy_multiplyer_Poisson_ratio_zy;
								prop[MULT_POISSON_RATIO_ZX][l] = matlist[imat_id_loc].orthotropy_multiplyer_Poisson_ratio_zx;
								prop[MULT_POISSON_RATIO_YX][l] = matlist[imat_id_loc].orthotropy_multiplyer_Poisson_ratio_yx;

								bActiveShearModule[l] = matlist[imat_id_loc].bActive_ShearModule;
								prop[SHEAR_MODULE_YZ][l] = matlist[imat_id_loc].ShearModule_yz;
								prop[SHEAR_MODULE_XZ][l] = matlist[imat_id_loc].ShearModule_xz;
								prop[SHEAR_MODULE_XY][l] = matlist[imat_id_loc].ShearModule_xy;
							}
							//Sc[evt[ic] - 1] = b[ib].Sc;
							if (b[ib].n_Sc == 1) {
								Sc[evt[ic] - 1] = b[ib].arr_Sc[0];
							}
							else {
								Sc[evt[ic] - 1] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, 20.0);
							}
							ipower_time_depend[evt[ic] - 1] = b[ib].ipower_time_depend;

							//Sc[l]= get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, 20.0);
							//ipower_time_depend[l] = b[ib].ipower_time_depend;
							break;
						case HYDRODINAMIC: prop[RHO][l] = matlist[imat_id_loc].rho;
							prop[MU_DYNAMIC_VISCOSITY][l] = matlist[imat_id_loc].mu;
							prop[BETA_T][l] = matlist[imat_id_loc].beta_t;
							break;
						}

						l++;
					}
					else {
						std::cout << "model is incorrect: prop indx >=maxelm"<<std::endl;
						std::cout << "see constr struct file."<<std::endl;
						system("pause");
						exit(1);
					}
				}
		}
	}
} // constr_prop

// Заносит свойства материалов в структуру
// для жидкой зоны с номером iDom.
void constr_prop_flow(int* &evt, int* &whot_is_block, 
	int** &evt_f2, integer iDom, float** &prop,
					  integer maxelm, BLOCK* &b, integer lb, 
					  int inx, int iny, int inz, 
				      doublereal* &xpos, doublereal* &ypos, doublereal* &zpos,
					  TPROP* &matlist) {
	integer i;
	prop=nullptr;
	prop=new float*[3];
	if (prop==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for prop flow constr struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (i=0; i<3; ++i) prop[i]=nullptr;
	for (i=0; i<3; ++i) {
		prop[i]=new float[maxelm];
		if (prop[i]==nullptr) {
	        // недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for prop["<< i <<"] constr struct..." << std::endl;
		    
		    std::cout << "Please any key to exit..." << std::endl;
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}

	//integer j,k,ic;
	//TOCHKA p;
    bool inDomain=false;
	 integer ib,l=0; 
	// integer ic_k, ic_j;

	// Проход по всем контрольным объёмам в порядке обхода:
	for (int k=0; k<inz; ++k)
	{
		integer  ic_k=k*inx*iny;
		for (int  j=0; j<iny; ++j) 
		{
			integer   ic_j=j*inx+ic_k;
			for (int  i=0; i<inx; ++i) 
			{

				integer  ic=i+ic_j; // номер КО
        
				if ((evt_f2[ENUMERATECONTVOL][ic]>0) && (evt_f2[MASKDOMAINFLUID][ic]==(iDom+1))) {

						//p.x = 0.5*(xpos[i] + xpos[i + 1]);
						//p.y = 0.5*(ypos[j] + ypos[j + 1]);
						//p.z = 0.5*(zpos[k] + zpos[k + 1]);

						//inDomain = in_model_flow(p, ib, b, lb); // определяем номер блока жидкой зоны

						integer iP = evt[ic];
						if (iP > 0) {
							ib = whot_is_block[iP - 1];
							if (b[ib].itype != PHYSICS_TYPE_IN_BODY::FLUID) {
								inDomain = false;
							}
							else inDomain = true;
						}
						else inDomain = false;


						if (inDomain) {
							if (l < maxelm) {
								prop[RHO][l] = matlist[b[ib].imatid].rho;
								prop[MU_DYNAMIC_VISCOSITY][l] = matlist[b[ib].imatid].mu;
								prop[BETA_T][l] = matlist[b[ib].imatid].beta_t;

								l++;
							}
							else {
								std::cout << "model is incorrect: prop 2 indx>=maxelm" << std::endl;
								std::cout << "see file constr struct in function constr_prop_flow" << std::endl;
								system("pause");
								exit(1);
							}
					}
				}
			}
		}
	}
} // constr_prop_flow

// Алгоритм уникальной нумерации граней.
// В результате каждая граничная грань,
// принадлежащая расчётной области имеет
// уникальный номер.
// Внимание !!! Грани внутренних источников тепла также пронумерованы.
void enumerate_gran_temp(int** &gran, int maxelm, int** &nvtx,
						 int &maxbound, int** &neighbour, TOCHKA* &pa,
						 SOURCE* &s, integer ls) {
	// Выделение оперативной памяти.
	if (gran != nullptr) {
		for (integer i = 0; i < 6; ++i) {
			if (gran[i] != nullptr) {
				delete[] gran[i];
				gran[i] = nullptr;
			}
		}
		delete[] gran;
	}
	//gran = nullptr;

	gran=nullptr;
	gran=new int*[6];
	if (gran==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for gran constr struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	
	for (int i=0; i<6; ++i) gran[i]=nullptr;
	for (int i=0; i<6; ++i) {
		gran[i]=new int[maxelm];
		if (gran[i]==nullptr) {
	        // недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for gran[" << i <<"] constr struct..." << std::endl;
		    
		    std::cout << "Please any key to exit..." << std::endl;
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}
	// инициализация:
	for (integer i = 0; i < 6; ++i) {

#pragma omp parallel for
		for (integer j = 0; j < maxelm; ++j) {
			gran[i][j] = -1;
		}
	}

	bool *bvisit=nullptr;
	bvisit=new bool[maxelm];
	// Результат работы оператора new не требует проверки на null.
	//if (bvisit==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
#pragma omp parallel for 
	for (integer j = 0; j < maxelm; ++j) {
		bvisit[j] = false; // признак посещения узла.
	}

	maxbound=0; // общее число граней
	
	bool bfind=false;
	
	for (integer i=0; i<maxelm; ++i) {
		// цикл по всем соседям
		// G текущая грань.
		for (integer G=0; G<6; G++) {
			if (neighbour[G][i]==0) {
				// сосед не существует
				gran[G][i]=maxbound;
				maxbound++;
				// записываем идентификатор грани и
				// увеличиваем счётчик граней.
			}
			else {
				// сосед существует, грань внутренняя
				if (bvisit[neighbour[G][i]-1]) {

					integer GR; // G Reverse
					// узел уже был посещён
					switch (G) {
						case E_SIDE: GR=W_SIDE; break;
						case W_SIDE: GR=E_SIDE; break;
						case N_SIDE: GR=S_SIDE; break;
						case S_SIDE: GR=N_SIDE; break;
						case T_SIDE: GR=B_SIDE; break;
						case B_SIDE: GR=T_SIDE; break;
					}
					gran[G][i]=gran[GR][neighbour[G][i]-1];
				}
				else {

					doublereal x_c = 0.0, y_c = 0.0, z_c = 0.0;
					integer iplane;

					// узнать координаты центра грани и ориентацию в пространстве
					switch (G) {
						case E_SIDE: x_c=pa[nvtx[1][i]-1].x;
							     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
								 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
								 iplane=YZ_PLANE;
							     break;
						case W_SIDE: x_c=pa[nvtx[0][i]-1].x;
							     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
								 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
								 iplane=YZ_PLANE;
							     break;
						case N_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
							     y_c=pa[nvtx[2][i]-1].y;
								 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
								 iplane=XZ_PLANE;
							     break;
						case S_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
							     y_c=pa[nvtx[0][i]-1].y;
								 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
								 iplane=XZ_PLANE;
							     break;
						case T_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                                 y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
								 z_c=pa[nvtx[4][i]-1].z;
                                 iplane=XY_PLANE;
							     break;
						case B_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                                 y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
								 z_c=pa[nvtx[0][i]-1].z;
                                 iplane=XY_PLANE;
							     break;
					} // end case
					bfind=false;
					for (integer j=0; j<ls; ++j) {
						if (s[j].iPlane==iplane) {
							switch (iplane) {
							case XY_PLANE: s[j].g.zE = s[j].g.zS;
									if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(z_c-s[j].g.zE)<admission)) bfind=true;
									      break;
							case XZ_PLANE:  s[j].g.yE = s[j].g.yS;
									if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (fabs(y_c-s[j].g.yE)<admission)) bfind=true;
									      break;
							case YZ_PLANE: s[j].g.xE = s[j].g.xS;
									if ((z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(x_c-s[j].g.xE)<admission)) bfind=true;
									      break;
							}	
						}
						if (bfind) break; // досрочный выход из цикла for.
					}
					if (bfind) {
						// нужно присвоить грани соответствующий номер maxbound.
                        gran[G][i]=maxbound;
						maxbound++;
					}
				}
			}
		}
		bvisit[i]=true; // узел посещён
	}

	//if (bvisit != nullptr) {
		delete[] bvisit;
		bvisit = nullptr;
	//}
	// После этого обхода имеем количество граней maxbound.
	// Для каждой грани знаем её позицию в сетке.

} // enumerate_gran_temp

// Алгоритм уникальной нумерации граней.
// В результате каждая граничная грань,
// принадлежащая расчётной области имеет
// уникальный номер.
void enumerate_gran_flow(int** &gran, int maxelm, int** &nvtx,
						 int &maxbound, int** &neighbour, TOCHKA* &pa) {
	// Выделение оперативной памяти.
	gran=nullptr;
	gran=new int*[6];
	if (gran==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for gran flow constr struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	integer i=0, j=0;
	for (i=0; i<6; ++i) gran[i]=nullptr;
	for (i=0; i<6; ++i) {
		gran[i]=new int[maxelm];
		if (gran[i]==nullptr) {
	        // недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for gran["<< i << "] constr struct..." << std::endl;
		   
		    std::cout << "Please any key to exit..." << std::endl;
			//system("PAUSE");
			system("pause");
			exit(1);
	    }
	}
	// инициализация:
	for (i = 0; i < 6; ++i) {
#pragma omp for private(j)
		for (j = 0; j < maxelm; ++j) {
			gran[i][j] = -1;
		}
	}

	maxbound=0; // общее число граней
	integer G; // текущая грань

	for (i=0; i<maxelm; ++i) {
		// цикл по всем соседям
		for (G=0; G<6; ++G) {
			if (neighbour[G][i]==0) {
				// сосед не существует
				gran[G][i]=maxbound++;
				// записываем идентификатор грани и
				// увеличиваем счётчик граней.
			}
		}	
	}

   
	// После этого обхода имеем количество граней maxbound.
	// Для каждой грани знаем её позицию в сетке.

} // enumerate_gran_flow

// Вычисление соседей для каждого внутреннего узла. 
// Причём множество соседей выбирается и среди внутренних КО и 
// среди граничных КО.
// Универсальность: подходит и для TEMPER и для FLOW
void constr_neighbors_for_the_internal_node(
	int*** &neighbors_for_the_internal_node,
	int maxelm, int** &gran,
						int** &neighbour) 
{
	

	int inumcor = number_cores();

#ifdef _OPENMP
	omp_set_num_threads(inumcor); // установка числа потоков
#else
	inumcor = 1; // один поток.
#endif

	// Многократное выделение оперативной памяти для neighbors_for_the_internal_node занимает очень много времени.
	//printf("incomming maxelm=%lld\n",maxelm);
	//system("pause");

	//integer i=0;
	// Выделение оперативной памяти.
	neighbors_for_the_internal_node=nullptr;
	neighbors_for_the_internal_node = new int**[12];
	if (neighbors_for_the_internal_node==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for neighbors_for_the_internal_node constr struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer i=0; i<12; ++i) neighbors_for_the_internal_node[i]=nullptr;
	// Возможно выделение оперативной памяти параллельно не самая лучшая идея из-за гонки за ресурс.
	// Распределением ресурса должна заниматься операционная система.
//#pragma omp parallel for
	for (integer i=0; i<12; ++i) {

		neighbors_for_the_internal_node[i] = new int* [4];

		if (neighbors_for_the_internal_node[i] == nullptr) {
			// недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for neighbors_for_the_internal_node[" << i << "] constr struct..." << std::endl;

			std::cout << "Please any key to exit..." << std::endl;
			//system("PAUSE");
			system("pause");
			exit(1);
		}

		if (b_on_adaptive_local_refinement_mesh) {
			// Адаптивная локальная измельченная сетка.

			for (integer j = 0; j < 4; ++j) {
				neighbors_for_the_internal_node[i][j] = new int[maxelm];
			}
		}
		else {
			// Структурированная сетка.

			neighbors_for_the_internal_node[i][0] = new int[maxelm];

			neighbors_for_the_internal_node[i][1] = nullptr;
			neighbors_for_the_internal_node[i][2] = nullptr;
			neighbors_for_the_internal_node[i][3] = nullptr;

		}

		
		
	}

	
	


		//integer G, GG; // текущая грань
		// проход по всем граням принадлежащим текущему КО в определённом порядке.
		for (integer G=0; G<6; G++) {

			integer GG;

			// Вычисление дальнего соседа
			switch (G) {
			case E_SIDE: GG = EE_SIDE; break; // E
			case N_SIDE: GG = NN_SIDE; break; // N
			case T_SIDE: GG = TT_SIDE; break; // T
			case W_SIDE: GG = WW_SIDE; break; // W
			case S_SIDE: GG = SS_SIDE; break; // S
			case B_SIDE: GG = BB_SIDE; break; // B
			}

			// Обход в цикле по всем внутренним КО принадлежащим расчётной области
#pragma omp parallel for
			for (integer i = 0; i<maxelm; ++i) {

			if (gran[G][i]>-1) {
				// граничная грань:
				neighbors_for_the_internal_node[G][0][i]=maxelm+gran[G][i]; // граничные КО нумеруются в последнюю очередь,
				if (neighbors_for_the_internal_node[G][1] != nullptr) {
					neighbors_for_the_internal_node[G][1][i] = -1;
				}
				if (neighbors_for_the_internal_node[G][2] != nullptr) {
					neighbors_for_the_internal_node[G][2][i] = -1;
				}
				if (neighbors_for_the_internal_node[G][3] != nullptr) {
					neighbors_for_the_internal_node[G][3][i] = -1;
				}
				// после maxelm внутренних КО.
								
				neighbors_for_the_internal_node[GG][0][i] = -1; // соседа нет.
				if (neighbors_for_the_internal_node[GG][1] != nullptr) {
					neighbors_for_the_internal_node[GG][1][i] = -1;
				}
				if (neighbors_for_the_internal_node[GG][2] != nullptr) {
					neighbors_for_the_internal_node[GG][2][i] = -1;
				}
				if (neighbors_for_the_internal_node[GG][3] != nullptr) {
					neighbors_for_the_internal_node[GG][3][i] = -1;
				}
				

			}
			else {
				// соседом является внутренний КО.
				neighbors_for_the_internal_node[G][0][i]=neighbour[G][i]-1;
				if (neighbors_for_the_internal_node[G][1] != nullptr) {
					neighbors_for_the_internal_node[G][1][i] = -1;
				}
				if (neighbors_for_the_internal_node[G][2] != nullptr) {
					neighbors_for_the_internal_node[G][2][i] = -1;
				}
				if (neighbors_for_the_internal_node[G][3] != nullptr) {
					neighbors_for_the_internal_node[G][3][i] = -1;
				}
				

				if (gran[G][neighbour[G][i]-1]>-1) {
					// граничный КО
					neighbors_for_the_internal_node[GG][0][i] = maxelm+gran[G][neighbour[G][i]-1];
					if (neighbors_for_the_internal_node[GG][1] != nullptr) {
						neighbors_for_the_internal_node[GG][1][i] = -1;
					}
					if (neighbors_for_the_internal_node[GG][2] != nullptr) {
						neighbors_for_the_internal_node[GG][2][i] = -1;
					}
					if (neighbors_for_the_internal_node[GG][3] != nullptr) {
						neighbors_for_the_internal_node[GG][3][i] = -1;
					}
				}
				else {
					neighbors_for_the_internal_node[GG][0][i] = neighbour[GG][i] - 1; // внутренний КО.
					if (neighbors_for_the_internal_node[GG][1] != nullptr) {
						neighbors_for_the_internal_node[GG][1][i] = -1;
					}
					if (neighbors_for_the_internal_node[GG][2] != nullptr) {
						neighbors_for_the_internal_node[GG][2][i] = -1;
					}
					if (neighbors_for_the_internal_node[GG][3] != nullptr) {
						neighbors_for_the_internal_node[GG][3][i] = -1;
					}
				}

			}
		}
	}

	// neighbors_for_the_internal_node: maxelm -> maxp == maxelm + maxbound;

	// Признак того что узел граничный: номер КО >= maxelm.
	// Массив boundary (определяющий граничные узлы) становится ненужен.

	//omp_set_num_threads(1);

} // constr_neighbors_for_the_internal_node

// соседи в плоскости границы для граничных
// узлов. Это требуется для симметризации 
// портрета матрицы СЛАУ.
void constr_boundary_neighbour(integer G, integer i,
						   BOUND* &border_neighbor, 
						   int** &gran, int** &neighbour,
						   int*** &neighbors_for_the_internal_node,
						   integer maxelm) {
	int j=0; // счётчик цикла for

    for (j=0; j<6; ++j) border_neighbor[gran[G][i]].iW[j]=-1; // инициализация 

	// 13 10 2016 Выключил данный код т.к. он считается ненужным.
	/*
	// Эти значения требуются для симметризации портрета матрицы СЛАУ:
	switch (G) {
	case ESIDE: if ((neighbour[NSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[NSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[NSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[NSIDE][i].iNODE1];
		if ((neighbour[SSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[SSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[SSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[SSIDE][i].iNODE1];
		if ((neighbour[TSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[TSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[TSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[TSIDE][i].iNODE1];
		if ((neighbour[BSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[BSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[BSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[BSIDE][i].iNODE1];
			  break;
	case WSIDE: if ((neighbour[NSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[NSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[NSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[NSIDE][i].iNODE1];
		if ((neighbour[SSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[SSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[SSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[SSIDE][i].iNODE1];
		if ((neighbour[TSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[TSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[TSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[TSIDE][i].iNODE1];
		if ((neighbour[BSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[BSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[BSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[BSIDE][i].iNODE1];
			  break;
	case NSIDE: if ((neighbour[ESIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[ESIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[ESIDE] = maxelm + gran[G][neighbors_for_the_internal_node[ESIDE][i].iNODE1];
		if ((neighbour[WSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[WSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[WSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[WSIDE][i].iNODE1];
		if ((neighbour[TSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[TSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[TSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[TSIDE][i].iNODE1];
		if ((neighbour[BSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[BSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[BSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[BSIDE][i].iNODE1];
			  break;
	case SSIDE: if ((neighbour[ESIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[ESIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[ESIDE] = maxelm + gran[G][neighbors_for_the_internal_node[ESIDE][i].iNODE1];
		if ((neighbour[WSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[WSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[WSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[WSIDE][i].iNODE1];
		if ((neighbour[TSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[TSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[TSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[TSIDE][i].iNODE1];
		if ((neighbour[BSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[BSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[BSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[BSIDE][i].iNODE1];
			  break;
	case TSIDE: if ((neighbour[ESIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[ESIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[ESIDE] = maxelm + gran[G][neighbors_for_the_internal_node[ESIDE][i].iNODE1];
		if ((neighbour[WSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[WSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[WSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[WSIDE][i].iNODE1];
		if ((neighbour[NSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[NSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[NSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[NSIDE][i].iNODE1];
		if ((neighbour[SSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[SSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[SSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[SSIDE][i].iNODE1];
		      break;
	case BSIDE: if ((neighbour[ESIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[ESIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[ESIDE] = maxelm + gran[G][neighbors_for_the_internal_node[ESIDE][i].iNODE1];
		if ((neighbour[WSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[WSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[WSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[WSIDE][i].iNODE1];
		if ((neighbour[NSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[NSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[NSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[NSIDE][i].iNODE1];
		if ((neighbour[SSIDE][i]>0) && (gran[G][neighbors_for_the_internal_node[SSIDE][i].iNODE1]>-1)) border_neighbor[gran[G][i]].iW[SSIDE] = maxelm + gran[G][neighbors_for_the_internal_node[SSIDE][i].iNODE1];
		      break;
	} // end сбор информации для симметризации СЛАУ.

	*/

} // constr_boundary_neighbour

// Заполнение информации о граничных узлах:
// Заполняются border_neighbor и binternalsource.
void constr_border_neighbor_temp(BOUND* &border_neighbor, int* &whot_is_block,
	bool* &binternalsource, int maxelm, int maxbound, 
	int** &gran, int** &neighbour,
	int*** &neighbors_for_the_internal_node,
	int** &nvtx, TOCHKA* &pa,
	BLOCK* &b, int lb, int lw, WALL* &w, SOURCE* &s, int ls,
	TOCHKA pavg) {

	// Алгоритм:
    // Выделение оперативной памяти:
    // gran[0..5][0..maxbound-1]
	border_neighbor=nullptr;
	border_neighbor= new BOUND[maxbound];
	// оператор new не требует после себя проверки на null
	//if (border_neighbor==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for border_neighbor constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	binternalsource=nullptr;
    binternalsource = new bool[maxbound];
	// оператор new не требует после себя проверки на null
	//if (binternalsource==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for binternalsource constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	//bool bfluidsolid=false;

    // Проход по всем внутренним КО принадлежащим расчётной области.
	
#pragma omp parallel for 
	for (integer i = 0; i < maxbound; ++i) {
		binternalsource[i] = false; // инициализация.
	}
	    	

	bool *bvisit=nullptr;
	bvisit=new bool[maxbound];
	// оператор new не требует после себя проверки на null
	//if (bvisit==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
#pragma omp parallel for 
	for (integer i = 0; i < maxbound; ++i) {
		bvisit[i] = false;
	}

	//integer il=0; // debug
    //for (i=0; i<maxelm; ++i) for (G=0; G<6; G++) if (gran[G][i]>-1) il++;
    //for (i=0; i<maxelm; ++i) for (G=0; G<6; G++) if (neighbour[G][i]==0) il++;
#if doubleintprecision == 1
	//printf("il=%lld\n",il); system("PAUSE");
#else
	//printf("il=%d\n",il); system("PAUSE");
#endif
   
	// Препроцессинг.
	integer*** wlist = new integer**[4];
	integer** iwsize = new integer*[4];
	for (integer j = 0; j < 4; ++j) {
		if (j == 0) {
			wlist[j] = nullptr;
			iwsize[j] = nullptr;
		}
		else {			
			wlist[j] = new integer*[4];
			iwsize[j] = new integer[4];
			for (integer j_1 = 0; j_1 < 4; ++j_1) {
				wlist[j][j_1] = new integer[lw];
				iwsize[j][j_1] = 0;
			}
		}
		
	}


	for (integer j = 0; j < lw; ++j) {
		switch (w[j].iPlane) {
		case XY_PLANE:
			if (w[j].g.xS <= pavg.x) {
				if (w[j].g.yS <= pavg.y) {
					wlist[XY_PLANE][0][iwsize[XY_PLANE][0]] = j;
					iwsize[XY_PLANE][0]++;
				}
				if (w[j].g.yE >= pavg.y) {
					wlist[XY_PLANE][1][iwsize[XY_PLANE][1]] = j;
					iwsize[XY_PLANE][1]++;
				}
			}
			if (w[j].g.xE >= pavg.x)  {
				if (w[j].g.yS <= pavg.y) {
					wlist[XY_PLANE][2][iwsize[XY_PLANE][2]] = j;
					iwsize[XY_PLANE][2]++;
				}
				if (w[j].g.yE >= pavg.y) {
					wlist[XY_PLANE][3][iwsize[XY_PLANE][3]] = j;
					iwsize[XY_PLANE][3]++;
				}
			}

			
			break;
		case YZ_PLANE:
			if (w[j].g.yS <= pavg.y) {
				if (w[j].g.zS <= pavg.z) {
					wlist[YZ_PLANE][0][iwsize[YZ_PLANE][0]] = j;
					iwsize[YZ_PLANE][0]++;
				}
				if (w[j].g.zE >= pavg.z)  {
					wlist[YZ_PLANE][1][iwsize[YZ_PLANE][1]] = j;
					iwsize[YZ_PLANE][1]++;
				}
			}
			if (w[j].g.yE >= pavg.y) {
				if (w[j].g.zS <= pavg.z) {
					wlist[YZ_PLANE][2][iwsize[YZ_PLANE][2]] = j;
					iwsize[YZ_PLANE][2]++;
				}
				if (w[j].g.zE >= pavg.z) {
					wlist[YZ_PLANE][3][iwsize[YZ_PLANE][3]] = j;
					iwsize[YZ_PLANE][3]++;
				}
			}

			
			break;
		case XZ_PLANE:
			if (w[j].g.xS <= pavg.x) {
				if (w[j].g.zS <= pavg.z) {
					wlist[XZ_PLANE][0][iwsize[XZ_PLANE][0]] = j;
					iwsize[XZ_PLANE][0]++;
				}
				if (w[j].g.zE >= pavg.z) {
					wlist[XZ_PLANE][1][iwsize[XZ_PLANE][1]] = j;
					iwsize[XZ_PLANE][1]++;
				}
			}
			if (w[j].g.xE >= pavg.x) {
				if (w[j].g.zS <= pavg.z) {
					wlist[XZ_PLANE][2][iwsize[XZ_PLANE][2]] = j;
					iwsize[XZ_PLANE][2]++;
				}
				if (w[j].g.zE >= pavg.z)  {
					wlist[XZ_PLANE][3][iwsize[XZ_PLANE][3]] = j;
					iwsize[XZ_PLANE][3]++;
				}
			}
			
			break;
		}
	}


#pragma omp parallel for
	for (int i=0; i<maxelm; ++i) {
		// проход по всем граням внутреннего КО в особом порядке.
		for (integer G=0; G<6; G++) {
			if (gran[G][i]>-1) {

				doublereal x_c = 0.0, y_c = 0.0, z_c = 0.0; // координаты центра грани
				integer iplane; // плоскость в которой лежит грань.

				// граничная грань
				doublereal dS = 0.0;
				TOCHKA control_volume_size_dopusk;
                // узнать координаты центра грани и ориентацию в пространстве
				switch (G) {
					case E_SIDE: x_c=pa[nvtx[1][i]-1].x;
						     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 control_volume_size_dopusk.x = 0.25*fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x);
							 dS = fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=YZ_PLANE;
						     break;
					case W_SIDE: x_c=pa[nvtx[0][i]-1].x;
						     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 control_volume_size_dopusk.x = 0.25*fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x);
							 dS = fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=YZ_PLANE;
						     break;
					case N_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						     y_c=pa[nvtx[2][i]-1].y;
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 control_volume_size_dopusk.y = 0.25*fabs(pa[nvtx[0][i] - 1].y- pa[nvtx[2][i] - 1].y);
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=XZ_PLANE;
						     break;
					case S_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						     y_c=pa[nvtx[0][i]-1].y;
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 control_volume_size_dopusk.y = 0.25*fabs(pa[nvtx[0][i] - 1].y - pa[nvtx[2][i] - 1].y);
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=XZ_PLANE;
						     break;
					case T_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                             y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=pa[nvtx[4][i]-1].z;
							 control_volume_size_dopusk.z = 0.25*fabs(pa[nvtx[4][i] - 1].z- pa[nvtx[0][i] - 1].z);
                             iplane=XY_PLANE;
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y);
						     break;
				    case B_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                             y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
					    	 z_c=pa[nvtx[0][i]-1].z;
							 control_volume_size_dopusk.z = 0.25*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y);
                             iplane=XY_PLANE;
						     break;
				} // end case

				if (neighbour[G][i]==0) { // грань лежит на границе расчётной области.

					// соседи в плоскости границы для граничных
                    // узлов. Это требуется для симметризации 
                    // портрета матрицы СЛАУ.
                    constr_boundary_neighbour(G,i,border_neighbor, gran, neighbour, neighbors_for_the_internal_node, maxelm);

					

					// грань лежит на границе расчётной области.
					switch (G) {
						case E_SIDE: border_neighbor[gran[G][i]].Norm=W_SIDE;
							     border_neighbor[gran[G][i]].iII=neighbors_for_the_internal_node[W_SIDE][0][i];
							     break;
						case W_SIDE: border_neighbor[gran[G][i]].Norm=E_SIDE;
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[E_SIDE][0][i];
							     break;
						case N_SIDE: border_neighbor[gran[G][i]].Norm=S_SIDE;
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[S_SIDE][0][i];
							     break;
						case S_SIDE: border_neighbor[gran[G][i]].Norm = N_SIDE;
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[N_SIDE][0][i];
							     break;
						case T_SIDE: border_neighbor[gran[G][i]].Norm=B_SIDE; 
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[B_SIDE][0][i];
							     break;
						case B_SIDE: border_neighbor[gran[G][i]].Norm=T_SIDE; 
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[T_SIDE][0][i];
							     break;
					} // end определение внутренней нормали
                    border_neighbor[gran[G][i]].iB=maxelm+gran[G][i];
					border_neighbor[gran[G][i]].iI=i;
					border_neighbor[gran[G][i]].iI1 = -1;
					border_neighbor[gran[G][i]].iI2 = -1;
					bvisit[gran[G][i]]=true; // грань была посещена
					border_neighbor[gran[G][i]].dS = dS;
					// координаты центра грани.
					border_neighbor[gran[G][i]].p_c.x = x_c;
					border_neighbor[gran[G][i]].p_c.y = y_c;
					border_neighbor[gran[G][i]].p_c.z = z_c;
					// Вычисляем emissivity:
					integer ibfound = -1;
					ibfound = whot_is_block[i];
					if (ibfound >= 0) {
						// блок найден.
						// Определяем внутреннюю нормаль:
						switch (border_neighbor[gran[G][i]].Norm) {
						case W_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissE;
							break;
						case E_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissW;
							break;
						case S_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissN;
							break;
						case N_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissS;
							break;
						case B_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissT;
							break;
						case T_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissB;
							break;
						}
					}
					else {
						std::cout << "error: emissivity calculation block unfound. in constr_border_neighbor_temp." << std::endl;
						system("PAUSE");
						exit(1);
					}

					bool bfind=false;
					integer jpos=-1;
					for (integer j=0; j<ls; ++j) {
						if (s[j].iPlane==iplane) {
							// Важно не только попадание в фокус объекта но и нахождение с объектом на одном уровне:
							switch (iplane) {
								case XY_PLANE: if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(z_c-s[j].g.zE)<admission_bon_con)) { bfind=true; jpos=j; } break;
								case YZ_PLANE: if ((z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(x_c-s[j].g.xE)<admission_bon_con)) { bfind=true; jpos=j; } break;
								case XZ_PLANE: if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (fabs(y_c-s[j].g.yE)<admission_bon_con)) { bfind=true; jpos=j; } break;
							}
						}
					}
					if (bfind) {
						// внешний источник тепла
						//printf("source out found...\n"); // debug
						//system("PAUSE");
						border_neighbor[gran[G][i]].MCB=jpos;
					}
					else {

						
						if (1) {
							// Важно не только попадание в фокус объекта но и нахождение с объектом на одном уровне:
							switch (iplane) {
							case XY_PLANE:
								if (x_c <= pavg.x) {
									if (y_c <= pavg.y) {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][0]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{
												integer j_ind = wlist[iplane][0][j];
												if ((x_c > w[j_ind].g.xS) &&
													(x_c < w[j_ind].g.xE) &&
													(y_c > w[j_ind].g.yS) &&
													(y_c < w[j_ind].g.yE) &&
													(fabs(z_c - w[j_ind].g.zE) < control_volume_size_dopusk.z))
												{
													bfind = true; jpos = j_ind;
												}
												
											}
										}
									}
									else {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][1]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{
												integer j_ind = wlist[iplane][1][j];
												if ((x_c > w[j_ind].g.xS) &&
													(x_c < w[j_ind].g.xE) &&
													(y_c > w[j_ind].g.yS) &&
													(y_c < w[j_ind].g.yE) &&
													(fabs(z_c - w[j_ind].g.zE) < control_volume_size_dopusk.z))
												{
													bfind = true; jpos = j_ind;
												}
											}
										}
									}
								}
								else {
									if (y_c <= pavg.y) {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][2]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{
												integer j_ind = wlist[iplane][2][j];
												if ((x_c > w[j_ind].g.xS) &&
													(x_c < w[j_ind].g.xE) &&
													(y_c > w[j_ind].g.yS) &&
													(y_c < w[j_ind].g.yE) &&
													(fabs(z_c - w[j_ind].g.zE) < control_volume_size_dopusk.z))
												{
													bfind = true; jpos = j_ind;
												}
											}
										}
									}
									else {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][3]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{
												integer j_ind = wlist[iplane][3][j];
												if ((x_c > w[j_ind].g.xS) &&
													(x_c < w[j_ind].g.xE) &&
													(y_c > w[j_ind].g.yS) &&
													(y_c < w[j_ind].g.yE) &&
													(fabs(z_c - w[j_ind].g.zE) < control_volume_size_dopusk.z))
												{
													bfind = true; jpos = j_ind;
												}
											}
										}
									}
								}


								break;
							case YZ_PLANE:
								if (y_c <= pavg.y) {
									if (z_c <= pavg.z) {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][0]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{

												integer j_ind = wlist[iplane][0][j];
												if ((z_c > w[j_ind].g.zS) &&
													(z_c < w[j_ind].g.zE) &&
													(y_c > w[j_ind].g.yS) &&
													(y_c < w[j_ind].g.yE) &&
													(fabs(x_c - w[j_ind].g.xE) < control_volume_size_dopusk.x))
												{
													bfind = true; jpos = j_ind;
												}

											}
										}
									}
									else {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][1]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{

												integer j_ind = wlist[iplane][1][j];
												if ((z_c > w[j_ind].g.zS) &&
													(z_c < w[j_ind].g.zE) &&
													(y_c > w[j_ind].g.yS) &&
													(y_c < w[j_ind].g.yE) &&
													(fabs(x_c - w[j_ind].g.xE) < control_volume_size_dopusk.x))
												{
													bfind = true; jpos = j_ind;
												}

											}
										}
									}
								}
								else {
									if (z_c <= pavg.z) {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][2]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{

												integer j_ind = wlist[iplane][2][j];
												if ((z_c > w[j_ind].g.zS) &&
													(z_c < w[j_ind].g.zE) &&
													(y_c > w[j_ind].g.yS) &&
													(y_c < w[j_ind].g.yE) &&
													(fabs(x_c - w[j_ind].g.xE) < control_volume_size_dopusk.x))
												{
													bfind = true; jpos = j_ind;
												}

											}
										}
									}
									else {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][3]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{

												integer j_ind = wlist[iplane][3][j];
												if ((z_c > w[j_ind].g.zS) &&
													(z_c < w[j_ind].g.zE) &&
													(y_c > w[j_ind].g.yS) &&
													(y_c < w[j_ind].g.yE) &&
													(fabs(x_c - w[j_ind].g.xE) < control_volume_size_dopusk.x))
												{
													bfind = true; jpos = j_ind;
												}

											}
										}
									}
								}


								break;
							case XZ_PLANE:
								if (x_c <= pavg.x) {
									if (z_c <= pavg.z) {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][0]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{
												integer j_ind = wlist[iplane][0][j];
												if ((x_c > w[j_ind].g.xS) &&
													(x_c < w[j_ind].g.xE) &&
													(z_c > w[j_ind].g.zS) &&
													(z_c < w[j_ind].g.zE) &&
													(fabs(y_c - w[j_ind].g.yE) < control_volume_size_dopusk.y))
												{
													bfind = true; jpos = j_ind;
												}
											}
										}
									}
									else {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][1]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{
												integer j_ind = wlist[iplane][1][j];
												if ((x_c > w[j_ind].g.xS) &&
													(x_c < w[j_ind].g.xE) &&
													(z_c > w[j_ind].g.zS) &&
													(z_c < w[j_ind].g.zE) &&
													(fabs(y_c - w[j_ind].g.yE) < control_volume_size_dopusk.y))
												{
													bfind = true; jpos = j_ind;
												}
											}
										}
									}
								}
								else {
									if (z_c <= pavg.z) {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][2]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{
												integer j_ind = wlist[iplane][2][j];
												if ((x_c > w[j_ind].g.xS) &&
													(x_c < w[j_ind].g.xE) &&
													(z_c > w[j_ind].g.zS) &&
													(z_c < w[j_ind].g.zE) &&
													(fabs(y_c - w[j_ind].g.yE) < control_volume_size_dopusk.y))
												{
													bfind = true; jpos = j_ind;
												}
											}
										}
									}
									else {
										for (integer j = 0; (!bfind) && (j < iwsize[iplane][3]); ++j) {
											//if (w[wlist[iplane][j]].iPlane==iplane)
											{
												integer j_ind = wlist[iplane][3][j];
												if ((x_c > w[j_ind].g.xS) &&
													(x_c < w[j_ind].g.xE) &&
													(z_c > w[j_ind].g.zS) &&
													(z_c < w[j_ind].g.zE) &&
													(fabs(y_c - w[j_ind].g.yE) < control_volume_size_dopusk.y))
												{
													bfind = true; jpos = j_ind;
												}
											}
										}
									}
								}


								break;
							}
						}
						else {
							for (integer j = 0; j < lw; ++j) {
								if (w[j].iPlane == iplane) {
									// Важно не только попадание в фокус объекта но и нахождение с объектом на одном уровне:
									switch (iplane) {
									case XY_PLANE: if ((x_c > w[j].g.xS) && (x_c < w[j].g.xE) && (y_c > w[j].g.yS) && (y_c < w[j].g.yE) && (fabs(z_c - w[j].g.zE) < admission_bon_con)) { bfind = true; jpos = j; } break;
									case YZ_PLANE: if ((z_c > w[j].g.zS) && (z_c < w[j].g.zE) && (y_c > w[j].g.yS) && (y_c < w[j].g.yE) && (fabs(x_c - w[j].g.xE) < admission_bon_con)) { bfind = true; jpos = j; } break;
									case XZ_PLANE: if ((x_c > w[j].g.xS) && (x_c < w[j].g.xE) && (z_c > w[j].g.zS) && (z_c < w[j].g.zE) && (fabs(y_c - w[j].g.yE) < admission_bon_con)) { bfind = true; jpos = j; } break;
									}
								}
							}
						}
					   
					   if (bfind) {
						   // найдена стека
						   //printf("ambient wall found...\n"); // debug
						   //system("PAUSE");
                           border_neighbor[gran[G][i]].MCB=ls+jpos;
					   }
					   else {
						   // граничное условие по умолчанию
						   border_neighbor[gran[G][i]].MCB=ls+lw;
					   }

					}

				}
				else 
				{
					// внутренняя грань - это внутренний источник тепла.
					if (!bvisit[gran[G][i]]) {
						// если внутренняя грань ещё не была посещена.

						// Нужно проверить где находится твёрдое тело:
						// Нормаль будет направлена в сторону твёрдого тела.
						// координаты центра контрольного объёма

						TOCHKA p_c;
						p_c.x=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						p_c.y=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
						p_c.z=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);

						int ib; // номер блока
						bool bi_fluid=in_model_flow(p_c,ib,b,lb); // КО i принадлежит жидкой зоне.

						if (bi_fluid) {
							// КО i принадлежит жидкой зоне

							// Записываем в соседи тот КО который принадлежит твёрдому телу
							// с противоположной от жидкого КО стороны рассматриваемой грани.
                            border_neighbor[gran[G][i]].Norm=G; // внутренняя нормаль
							border_neighbor[gran[G][i]].iB=maxelm+gran[G][i];
							border_neighbor[gran[G][i]].iI=neighbour[G][i]-1;
							border_neighbor[gran[G][i]].iI1 = neighbour[G][i] - 1;
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[G][0][neighbour[G][i] - 1];
							border_neighbor[gran[G][i]].dS = dS;
							// координаты центра грани.
							border_neighbor[gran[G][i]].p_c.x = x_c;
							border_neighbor[gran[G][i]].p_c.y = y_c;
							border_neighbor[gran[G][i]].p_c.z = z_c;
							//bfluidsolid=true;							

							// Вычисляем emissivity:
							integer ibfound = -1;
							ibfound = whot_is_block[neighbour[G][i] - 1];
							if (ibfound >= 0) {
								// блок найден.
								// Определяем внутреннюю нормаль:
								switch (border_neighbor[gran[G][i]].Norm) {
								case W_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissE;
									break;
								case E_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissW;
									break;
								case S_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissN;
									break;
								case N_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissS;
									break;
								case B_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissT;
									break;
								case T_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissB;
									break;
								}
							}
							else {
								std::cout << "error: emissivity calculation block unfound. in constr_border_neighbor_temp." << std::endl;
								system("PAUSE");
								exit(1);
							}


						}
						else {
							// В будущем здесь также можно предусмотреть ситуацию
							// когда плоский источник окружают два твёрдых тела.
							// Пока это не сделано. Там нужно будет выбрать тогда
							// из двух твёрдых тел одно с наибольшей теплопроводностью.

							// КО i принадлежит SOLID

							switch (G) {
						       case E_SIDE: border_neighbor[gran[G][i]].Norm=W_SIDE;
								   border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[W_SIDE][0][i];
							            break;
						       case W_SIDE: border_neighbor[gran[G][i]].Norm=E_SIDE;
								   border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[E_SIDE][0][i];
							            break;
						       case N_SIDE: border_neighbor[gran[G][i]].Norm=S_SIDE;
								   border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[S_SIDE][0][i];
							            break;
						       case S_SIDE: border_neighbor[gran[G][i]].Norm=N_SIDE;
								   border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[N_SIDE][0][i];
							            break;
						       case T_SIDE: border_neighbor[gran[G][i]].Norm=B_SIDE; 
								   border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[B_SIDE][0][i];
							            break;
						       case B_SIDE: border_neighbor[gran[G][i]].Norm=T_SIDE; 
								   border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[T_SIDE][0][i];
							            break;
					         } // end определение внутренней нормали
                            border_neighbor[gran[G][i]].iB=maxelm+gran[G][i];
					        border_neighbor[gran[G][i]].iI=i;
							border_neighbor[gran[G][i]].iI1 = i;
							border_neighbor[gran[G][i]].dS = dS;
							// координаты центра грани.
							border_neighbor[gran[G][i]].p_c.x = x_c;
							border_neighbor[gran[G][i]].p_c.y = y_c;
							border_neighbor[gran[G][i]].p_c.z = z_c;
							//bfluidsolid=true;
							// здесь подразумевается что для всех внутренних граней
							// имеем ситуацию когда источник посередине между SOLID и FLUID.

							// Вычисляем emissivity:
							integer ibfound = -1;
							ibfound = whot_is_block[i];
							if (ibfound >= 0) {
								// блок найден.
								// Определяем внутреннюю нормаль:
								switch (border_neighbor[gran[G][i]].Norm) {
								case W_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissE;
									break;
								case E_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissW;
									break;
								case S_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissN;
									break;
								case N_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissS;
									break;
								case B_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissT;
									break;
								case T_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissB;
									break;
								}
							}
							else {
								std::cout << "error: emissivity calculation block unfound. in constr_border_neighbor_temp." << std::endl;
								system("PAUSE");
								exit(1);
							}

						}

						// соседи в плоскости границы для граничных
                        // узлов. Это требуется для симметризации 
                        // портрета матрицы СЛАУ.
                        constr_boundary_neighbour(G,i,border_neighbor, gran, neighbour, neighbors_for_the_internal_node, maxelm);

						integer jpos=-1;

						// Нужно определить MCB:
                        bool bfind=false;
					    for (integer j=0; j<ls; ++j) {
						   if (s[j].iPlane==iplane) {
							  switch (iplane) {
									case XY_PLANE: if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(z_c-s[j].g.zE)<admission_bon_con)) { bfind=true; jpos=j; } break;
								    case YZ_PLANE: if ((z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(x_c-s[j].g.xE)<admission_bon_con)) { bfind=true; jpos=j; } break;
								    case XZ_PLANE: if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (fabs(y_c-s[j].g.yE)<admission_bon_con)) { bfind=true; jpos=j; } break;
							  }
						   }
					    }
						
					    if (bfind) {
						    // внутрений источник тепла
							//printf("internal source found...\n"); // debug
							//system("PAUSE");
						    border_neighbor[gran[G][i]].MCB=jpos;
							binternalsource[gran[G][i]]=true; // внутренний источник тепла на границе жидкости и твёрдого тела
						}
						else {
							// Ошибка внутренний источник тепла не обнаружен!
							std::cout << "Error !!! internal source not found..." << std::endl;
							//system("PAUSE");
							system("pause");
							exit(0);
							// закомментировал 11.01.2020
							// как недостижимый код.
							//border_neighbor[gran[G][i]].MCB=ls+lw;
						}

					}
					else
					{

						TOCHKA p_c;
						p_c.x = 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x);
						p_c.y = 0.5*(pa[nvtx[1][i] - 1].y + pa[nvtx[3][i] - 1].y);
						p_c.z = 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z);

						int ib; // номер блока
						bool bi_fluid = in_model_flow(p_c, ib, b, lb);

						if (bi_fluid) {
							// КО i принадлежит жидкой зоне

							// Записываем в соседи тот КО который принадлежит твёрдому телу
							// с противоположной от жидкого КО стороны рассматриваемой грани.
							
							border_neighbor[gran[G][i]].iI2 = neighbour[G][i] - 1;
							

						}
						else {
							border_neighbor[gran[G][i]].iI2 = i;
						}
					}

					bvisit[gran[G][i]]=true; // грань была посещена
				}
			}
		}
	}

	// Освобождение оперативной памяти
	delete[] bvisit; 
	bvisit= nullptr;
	for (integer j = 0; j < 4; ++j) {
		for (integer j_1 = 0; j_1 < 4; ++j_1) {
			if (j != 0) {
				delete[] wlist[j][j_1];
				wlist[j][j_1] = nullptr;
			}
		}
		if (j != 0) {
			delete[] wlist[j];
			wlist[j] = nullptr;
			delete[] iwsize[j];
			iwsize[j] = nullptr;
		}
	}
	delete[] wlist;
	delete[] iwsize;
	wlist= nullptr;
	iwsize= nullptr;

	// border_neighbor: maxbound -> maxp == maxelm + maxbound;

} // constr_border_neighbor_temp 

// Заполнение информации о граничных узлах:
void constr_border_neighbor_flow(BOUND* &border_neighbor,  int* &whot_is_block,
								 int* &ptr, int maxelm, int maxbound,
	                     int** &gran, int** &neighbour, 
						 int*** &neighbors_for_the_internal_node,
						 int** &nvtx, TOCHKA* &pa,
						int lw, WALL* &w, int ls, BLOCK* &b) 
{

	// Алгоритм:
    // Выделение оперативной памяти:
    // gran[0..5][0..maxbound-1]
	border_neighbor=nullptr;
	border_neighbor= new BOUND[maxbound];
	if (border_neighbor==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for border_neighbor constr struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}

    // Проход по всем внутренним КО принадлежащим расчётной области.
	    

#pragma omp parallel for
	for (int i=0; i<maxelm; ++i) {
		// проход по всем граням внутреннего КО в особом порядке.
		// G - текущая грань.
		for (integer G=0; G<6; ++G) {
			if (gran[G][i]>-1) {
				// граничная грань
				
				doublereal dS = 0.0;

				doublereal x_c = 0.0, y_c = 0.0, z_c = 0.0; // координаты центра грани
				integer iplane; // плоскость в которой лежит грань.

                // узнать координаты центра грани и ориентацию в пространстве
				switch (G) {
					case E_SIDE: x_c=pa[nvtx[1][i]-1].x;
						     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=YZ_PLANE;
						     break;
					case W_SIDE: x_c=pa[nvtx[0][i]-1].x;
						     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=YZ_PLANE;
						     break;
					case N_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						     y_c=pa[nvtx[2][i]-1].y;
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=XZ_PLANE;
						     break;
					case S_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						     y_c=pa[nvtx[0][i]-1].y;
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=XZ_PLANE;
						     break;
					case T_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                             y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=pa[nvtx[4][i]-1].z;
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y);
                             iplane=XY_PLANE;
						     break;
				    case B_SIDE: x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                             y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
					    	 z_c=pa[nvtx[0][i]-1].z;
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y);
                             iplane=XY_PLANE;
						     break;
				} // end case

				if (neighbour[G][i]==0) {

                    // соседи в плоскости границы для граничных
                    // узлов. Это требуется для симметризации 
                    // портрета матрицы СЛАУ.
                    constr_boundary_neighbour(G,i,border_neighbor, gran, neighbour, neighbors_for_the_internal_node, maxelm);


					// грань лежит на границе расчётной области.
					switch (G) {
						case E_SIDE: border_neighbor[gran[G][i]].Norm=W_SIDE;
							     border_neighbor[gran[G][i]].iII=neighbors_for_the_internal_node[W_SIDE][0][i];
							     break;
						case W_SIDE: border_neighbor[gran[G][i]].Norm=E_SIDE;
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[E_SIDE][0][i];
							     break;
						case N_SIDE: border_neighbor[gran[G][i]].Norm=S_SIDE;
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[S_SIDE][0][i];
							     break;
						case S_SIDE: border_neighbor[gran[G][i]].Norm=N_SIDE;
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[N_SIDE][0][i];
							     break;
						case T_SIDE: border_neighbor[gran[G][i]].Norm=B_SIDE; 
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[B_SIDE][0][i];
							     break;
						case B_SIDE: border_neighbor[gran[G][i]].Norm=T_SIDE; 
							border_neighbor[gran[G][i]].iII = neighbors_for_the_internal_node[T_SIDE][0][i];
							     break;
					} // end определение внутренней нормали
                    border_neighbor[gran[G][i]].iB=maxelm+gran[G][i];
					border_neighbor[gran[G][i]].iI=i;
					border_neighbor[gran[G][i]].dS = dS; // площадь грани.
					// координаты центра грани.
					border_neighbor[gran[G][i]].p_c.x = x_c;
					border_neighbor[gran[G][i]].p_c.y = y_c;
					border_neighbor[gran[G][i]].p_c.z = z_c;
					// Вычисляем emissivity:
					integer ibfound = -1;
					ibfound = whot_is_block[ptr[i]];
					if (ibfound >= 0) {
						// блок найден.
						// Определяем внутреннюю нормаль:
						switch (border_neighbor[gran[G][i]].Norm) {
						case W_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissE;
							break;
						case E_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissW;
							break;
						case S_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissN;
							break;
						case N_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissS;
							break;
						case B_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissT;
							break;
						case T_SIDE: border_neighbor[gran[G][i]].emissivity = b[ibfound].radiation.emissB;
							break;
						}
					}
					else {
						std::cout << "error: emissivity calculation block unfound. in constr_border_neighbor_temp." << std::endl;
						system("PAUSE");
						exit(1);
					}
					
					integer jpos = 0; // счётчик цикла.
					bool bfind=false;
					for (integer j=0; j<lw; ++j) {
						if (w[j].iPlane==iplane) {
							switch (iplane) {
								case XY_PLANE: if ((x_c>w[j].g.xS) && (x_c<w[j].g.xE) && (y_c>w[j].g.yS) && (y_c<w[j].g.yE) && (fabs(z_c-w[j].g.zE)<admission_bon_con)) 
										  {
											  //printf("z_c=%f, zE=%f\n",z_c, w[j].g.zE);
											  //system("PAUSE"); // debug
											  bfind=true; jpos=j;
										  }
										  break;
								case YZ_PLANE: if ((z_c>w[j].g.zS) && (z_c<w[j].g.zE) && (y_c>w[j].g.yS) && (y_c<w[j].g.yE) && (fabs(x_c-w[j].g.xE)<admission_bon_con)) { bfind=true; jpos=j; } break;
								case XZ_PLANE: if ((x_c>w[j].g.xS) && (x_c<w[j].g.xE) && (z_c>w[j].g.zS) && (z_c<w[j].g.zE) && (fabs(y_c-w[j].g.yE)<admission_bon_con)) { bfind=true; jpos=j; } break;
							}
						}
					}
					if (bfind) {
						// внешний источник тепла
						border_neighbor[gran[G][i]].MCB=ls+jpos;
#if doubleintprecision == 1
						//printf("MCB=%lld\n",border_neighbor[gran[G][i]].MCB);
#else
						//printf("MCB=%d\n",border_neighbor[gran[G][i]].MCB);
#endif
						
						// system("PAUSE"); // debug
					}
					else {
                         // граничное условие по умолчанию
					     border_neighbor[gran[G][i]].MCB=ls+lw;
						 //printf("default\n");
						 //system("PAUSE"); // debug
					}
				} 
			}
		}
	}

	
	// border_neighbor: maxbound -> maxp == maxelm + maxbound;

} // constr_border_neighbor_flow

// Свойства материала на границе твердотельной области.
void constr_prop_bound(float** &prop, float** &prop_b, int maxelm, int maxbound,
					   int** &gran, int** &neighbour, int** &nvtx, TOCHKA* &pa,
					   BLOCK* &b, int lb) {
	// Выделение оперативной памяти:
	prop_b=nullptr;
	prop_b=new float*[6];
	if (prop_b==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for prop_b constr struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	integer i=0; // счётчик
	for (i=0; i<6; ++i) prop_b[i]=nullptr;
	for (i=0; i<6; ++i) {
		prop_b[i]=new float[maxbound];
		if (prop_b[i]==nullptr) {
	        // недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for prop_b["<< i <<"] constr struct..."<< std::endl;
		    std::cout << "Please any key to exit..."<< std::endl;
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}
	//integer G; // текущая грань
	bool *bvisit=nullptr;
	bvisit=new bool[maxbound];
	// Результат работы оператора new не требует проверки на null.
	//if (bvisit==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem: not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	for (i=0; i<maxbound; ++i) bvisit[i]=false; // признак посещаемости узлов
	TOCHKA p_c;
	int ib;
	bool bi_fluid;


    // в цикле по всем внутренним КО принадлежащим расчётной области:
	for (i=0; i<maxelm; ++i) {
		// проход по всем граням внутреннего КО в особом порядке.

		for (int G = 0; G < 6; ++G) {
			if (gran[G][i] > -1) {
				// граничная грань
				if (neighbour[G][i] == 0) {
					// грань лежит на границе расчётной области

					// В граничный узел сносятся свойства прилегающего КО:
					prop_b[RHO][gran[G][i]] = prop[RHO][i];
					prop_b[HEAT_CAPACITY][gran[G][i]] = prop[HEAT_CAPACITY][i];
					prop_b[LAM][gran[G][i]] = prop[LAM][i];
					prop_b[MULT_LAM_X][gran[G][i]] = prop[MULT_LAM_X][i];
					prop_b[MULT_LAM_Y][gran[G][i]] = prop[MULT_LAM_Y][i];
					prop_b[MULT_LAM_Z][gran[G][i]] = prop[MULT_LAM_Z][i];
					bvisit[gran[G][i]] = true;
				}
				else
				{
					// внутренняя грань
					// Это случай расположения плоского бесконечно тонкого источника тепла.
					// Есть два пути. Первый присвоить грани свойства твёрдого тела. Второй 
					// построить среднегармоническую аппроксимацию с воздухом.
					// Выберем первый путь.

					if (!bvisit[gran[G][i]]) {

						// Нужно проверить где находится твёрдое тело:
						// Нормаль будет направлена в сторону твёрдого тела.
						p_c.x = 0.5 * (pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x);
						p_c.y = 0.5 * (pa[nvtx[1][i] - 1].y + pa[nvtx[3][i] - 1].y);
						p_c.z = 0.5 * (pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z);

						bi_fluid = in_model_flow(p_c, ib, b, lb);

						if (bi_fluid) {
							prop_b[RHO][gran[G][i]] = prop[RHO][neighbour[G][i] - 1];
							prop_b[HEAT_CAPACITY][gran[G][i]] = prop[HEAT_CAPACITY][neighbour[G][i] - 1];
							prop_b[LAM][gran[G][i]] = prop[LAM][neighbour[G][i] - 1];
							prop_b[MULT_LAM_X][gran[G][i]] = prop[MULT_LAM_X][neighbour[G][i] - 1];
							prop_b[MULT_LAM_Y][gran[G][i]] = prop[MULT_LAM_Y][neighbour[G][i] - 1];
							prop_b[MULT_LAM_Z][gran[G][i]] = prop[MULT_LAM_Z][neighbour[G][i] - 1];
						}
						else {
							prop_b[RHO][gran[G][i]] = prop[RHO][i];
							prop_b[HEAT_CAPACITY][gran[G][i]] = prop[HEAT_CAPACITY][i];
							prop_b[LAM][gran[G][i]] = prop[LAM][i];
							prop_b[MULT_LAM_X][gran[G][i]] = prop[MULT_LAM_X][i];
							prop_b[MULT_LAM_Y][gran[G][i]] = prop[MULT_LAM_Y][i];
							prop_b[MULT_LAM_Z][gran[G][i]] = prop[MULT_LAM_Z][i];
						}

					}

					bvisit[gran[G][i]] = true;
				}
			}
	}


		
	}

	delete[] bvisit;
	bvisit=nullptr;

} // constr_prop_bound 

// Свойства материала на границе жидкой области.
void constr_prop_bound_flow(float** &prop, float** &prop_b, 
							integer maxelm, integer maxbound,
					   int** &gran, int** &neighbour)
{
	// Выделение оперативной памяти:
	prop_b=nullptr;
	prop_b=new float*[3];
	if (prop_b==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for prop_b constr struct..." << std::endl;
		std::cout << "Please any key to exit..." << std::endl;
	//	system("PAUSE");
		system("pause");
		exit(1);
	}
	integer i=0; // счётчик
	for (i=0; i<3; ++i) prop_b[i]=nullptr;
	for (i=0; i<3; ++i) {
		prop_b[i]=new float[maxbound];
		if (prop_b[i]==nullptr) {
	        // недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for prop_b["<< i<<"] constr struct..." << std::endl;		   
		    std::cout << "Please any key to exit..." << std::endl;
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}
	integer G; // текущая грань
	


    // в цикле по всем внутренним КО принадлежащим расчётной области:
	for (i=0; i<maxelm; ++i) {
		// проход по всем граням внутреннего КО в особом порядке.
		for (G=0; G<6; G++) {
			if (gran[G][i]>-1) {
				// граничная грань
				if (neighbour[G][i]==0) {
					// грань лежит на границе расчётной области

					// В граничный узел сносятся свойства прилегающего КО:
					prop_b[RHO][gran[G][i]]=prop[RHO][i];
					prop_b[MU_DYNAMIC_VISCOSITY][gran[G][i]]=prop[MU_DYNAMIC_VISCOSITY][i];
					prop_b[BETA_T][gran[G][i]]=prop[BETA_T][i];

				}
				
			}
		}
	}

} // constr_prop_bound_flow

// Нужно ли фиксировать давление в одной точке.
// По умолчанию используется последняя граничная точка.
// Данная функция универсальна и работает также и на АЛИС сетке.
// 24 сентября 2016.
void determination_of_activity_flow(BOUND* &border_neighbor, integer maxbound, 
									integer ls, integer lw, WALL* &w, 
									bool &bactiv, bool &bpressure_fix, bool &bLR1free) 
{
	integer i=0; // счётчик цикла for
	bool bpfix=true; // по умолчанию давление нужно фиксировать в одной точке

	// 10.02.2017
	//bpfix = false; так делать нельзя т.к. решение развалится.
	
	for (i=0; i<maxbound; ++i) {
		if ((border_neighbor[i].MCB>=ls) && (border_neighbor[i].MCB<(ls+lw))) {
			if (w[border_neighbor[i].MCB - ls].bpressure || 
				w[border_neighbor[i].MCB - ls].bopening) {
				// давление фиксировать ненужно
				bpfix=false;
				break; // досрочный выход из цикла for
			}
		}
	}

	bLR1free=bpfix; // применяем плавающий полилинейный солвер для поправки давления, т.к. он показывает более быструю сходимость.
	//bLR1free=false; // ->
	
	//->//bpfix=false; // если поправку не фиксировать в одной точке то алгоритм SIMPLE покажет более быструю сходимость.
    //->//bactiv=!bpfix; // Но это только в том случае если на всей границе расчётной области задана нулевая скорость.
	bactiv=true;
	/*
	bool bnonzeromag=false;

	// Если на границе области скорость ненулевая то нужно присвоить активность расчётной области:
	// Это как течение в каверне или тест Валь Девиса.
    for (i=0; i<maxbound; ++i) {
		if ((border_neighbor[i].MCB>=ls) && (border_neighbor[i].MCB<(ls+lw))) {
            if (!w[border_neighbor[i].MCB-ls].bpressure) {
				if ((fabs(w[border_neighbor[i].MCB-ls].Vx)>admission) || (fabs(w[border_neighbor[i].MCB-ls].Vy)>admission) || (fabs(w[border_neighbor[i].MCB-ls].Vz)>admission)) bnonzeromag=true;
			}
		}
	}
	if (!bactiv) bactiv=bnonzeromag;
	*/

    bpressure_fix=bpfix;

} // determination_of_activity_flow

// выделение оперативной памяти для задачи теплопроводности.
void allocation_memory_temp(doublereal* &potent, doublereal** &total_deformation, 
	equation3D* &sl, equation3D_bon* &slb, BOUND* &border_neighbor, 
	integer maxelm, integer maxbound, integer ls, integer lw, WALL* &w,
	doublereal operatingtemperature) {

    // выделение памяти под искомые полевые величины.
	potent=nullptr;
    potent=new doublereal[maxelm+maxbound];
	if (potent==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for potent constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}

	total_deformation = nullptr;
	if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
		(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
	{
		total_deformation = new doublereal*[SIZE_DEFORMATION_ARRAY];
		if (total_deformation == nullptr) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for total_deformation constr struct...\n");
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
		for (integer i_1 = 0; i_1 < SIZE_DEFORMATION_ARRAY; ++i_1) {
			total_deformation[i_1] = nullptr;
		}
		for (integer i_1 = 0; i_1 < SIZE_DEFORMATION_ARRAY; ++i_1) {
			total_deformation[i_1] = new doublereal[maxelm + maxbound];
			if (total_deformation[i_1] == nullptr) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment for total_deformation[%lld] constr struct...\n", i_1);
				printf("Please any key to exit...\n");
				//system("PAUSE");
				system("pause");
				exit(1);
			}
		}
	}

    doublereal Tamb=0.0;
    bool bTamb=false;
       // Давление должно быть инициализировано с учётом граничных условий.
       if (lw>0) {
          Tamb=1.0e+20;
	      for (integer i=0; i<maxbound; ++i) {
		      if (((border_neighbor[i].MCB<(ls+lw))&&
				  (border_neighbor[i].MCB>=ls)) && 
				  (w[border_neighbor[i].MCB-ls].ifamily== WALL_BOUNDARY_CONDITION::DIRICHLET_FAMILY)) {
			      // на границе задана температура
			      Tamb=fmin(Tamb,w[border_neighbor[i].MCB-ls].Tamb);
                  bTamb=true;
		      }
	      }
	   }

	   // Внимание инициализация значением operatingtemperature !!!.
       if (!bTamb) Tamb= operatingtemperature; // 0.0;

	// обнуление
	integer j=0;
#pragma omp parallel for
	for (j = 0; j < (maxelm + maxbound); ++j) {
		potent[j] = Tamb; // ноль градусов цельсия начальная температура.
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL_AND_TEMPERATURE) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL_AND_TEMPERATURE))
		{
			if (total_deformation != nullptr) {
				for (integer i_1 = 0; i_1 < SIZE_DEFORMATION_ARRAY; ++i_1) {
					total_deformation[i_1][j] = 0.0; // нулевая полная деформация
				}
			}
		}
	}

	sl=nullptr;
	if ((steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY)&&
		(steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T)) {

		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL))
		{
			// Отключаем выделение оперативной памяти под СЛАУ при решении чистой механики.
			// Т.к. для механики своя СЛАУ.
		}
		else {
			sl = new equation3D[maxelm]; // коэффициенты матрицы СЛАУ для внутренних КО.
			if (sl == nullptr) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment for slau temperature constr struct...\n");
				printf("Please any key to exit...\n");
				//system("PAUSE");
				system("pause");
				exit(1);
			}
		}
	}

	slb=nullptr;
	if ((steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) &&
		(steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T)) {
		if ((steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::STEADY_STATIC_STRUCTURAL) ||
			(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::UNSTEADY_STATIC_STRUCTURAL))
		{
			// Отключаем выделение оперативной памяти под СЛАУ при решении чистой механики.
			// Т.к. для механики своя СЛАУ.
		}
		else {
			slb = new equation3D_bon[maxbound]; // коэффициенты матрицы СЛАУ для граничных КО
			if (slb == nullptr) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem: not enough memory on your equipment for slau boundary temperature constr struct...\n");
				printf("Please any key to exit...\n");
				//system("PAUSE");
				system("pause");
				exit(1);
			}
		}
	}
} // allocation_memory_temp


// Возвращает true если внешние габариты расчётной области
// совпадают с внешними габаритами гидродинамической подобласти,
// считанной из файла load.txt.
bool is_EXTERNAL_FLOW(TOCHKA* pa1, int** nvtx1, integer maxelm1,
	doublereal* &x47, doublereal* &y47, doublereal* &z47, int** nvtx2, integer maxelm2,
	doublereal eps_minx, doublereal eps_miny, doublereal eps_minz,
	doublereal eps_maxx, doublereal eps_maxy, doublereal eps_maxz) {

	doublereal xmin1 = 1.0e37;
	doublereal ymin1 = 1.0e37;
	doublereal zmin1 = 1.0e37;
	doublereal xmax1 = 1.0e37;
	doublereal ymax1 = 1.0e37;
	doublereal zmax1 = 1.0e37;
	for (integer i = 0; i < maxelm1; ++i) {
		for (integer j = 0; j < 8; ++j) {
			if (pa1[nvtx1[j][i]-1].x < xmin1) xmin1 = pa1[nvtx1[j][i]-1].x;
			if (pa1[nvtx1[j][i]-1].y < ymin1) ymin1 = pa1[nvtx1[j][i]-1].y;
			if (pa1[nvtx1[j][i]-1].z < zmin1) zmin1 = pa1[nvtx1[j][i]-1].z;

			if (pa1[nvtx1[j][i]-1].x > xmax1) xmax1 = pa1[nvtx1[j][i]-1].x;
			if (pa1[nvtx1[j][i]-1].y > ymax1) ymax1 = pa1[nvtx1[j][i]-1].y;
			if (pa1[nvtx1[j][i]-1].z > zmax1) zmax1 = pa1[nvtx1[j][i]-1].z;
		}
	}

	doublereal xmin2 = 1.0e37;
	doublereal ymin2 = 1.0e37;
	doublereal zmin2 = 1.0e37;
	doublereal xmax2 = 1.0e37;
	doublereal ymax2 = 1.0e37;
	doublereal zmax2 = 1.0e37;
	for (integer i = 0; i < maxelm2; ++i) {
		for (integer j = 0; j < 8; ++j) {
			if (x47[nvtx2[j][i]] < xmin2) xmin2 = x47[nvtx2[j][i]];
			if (y47[nvtx2[j][i]] < ymin2) ymin2 = y47[nvtx2[j][i]];
			if (z47[nvtx2[j][i]] < zmin2) zmin2 = z47[nvtx2[j][i]];

			if (x47[nvtx2[j][i]] > xmax2) xmax2 = x47[nvtx2[j][i]];
			if (y47[nvtx2[j][i]] > ymax2) ymax2 = y47[nvtx2[j][i]];
			if (z47[nvtx2[j][i]] > zmax2) zmax2 = z47[nvtx2[j][i]];
		}
	}

	return ((fabs(xmin1-xmin2)<eps_minx)&&
		(fabs(ymin1 - ymin2)<eps_miny)&&
		(fabs(zmin1 - zmin2)<eps_minz)&&
		(fabs(xmax1 - xmax2)<eps_maxx)&&
		(fabs(ymax1 - ymax2)<eps_maxy)&&
		(fabs(zmax1 - zmax2)<eps_maxz));
}//is_EXTERNAL_FLOW


  // возвращает скорректированный массовый поток.
  // скорректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
  // скорректированной скорости или просто скорости на основе простейшей интерполяции.
  // Это нужно для отдельного решения уравнения конвекции-диффузии где задана пользовательская скорость,
  // никакой поправки Рхи-Чоу просто интерполяция.
  // 26.03.2017 15.09.2018 Декларация для использования уже здесь. Реализация в файле pamendment3.c.
void return_calc_correct_mass_flux_only_interpolation(integer iP, doublereal** potent, TOCHKA* pa, float** prop, float** prop_b,
	int** nvtx, int*** neighbors_for_the_internal_node, integer maxelm,
	doublereal* &mfcurrentretune, BOUND* &border_neighbor, int &ls, int &lw, 
	integer* ilevel_alice, int* ptr);

// выделение оперативной памяти для задачи гидродинамики.
void allocation_memory_flow(doublereal** &potent, 
							equation3D** &sl, equation3D_bon** &slb,
	BOUND* &border_neighbor, integer maxelm, integer maxbound,
	doublereal* &alpha, integer ls, integer lw, WALL* &w,
	doublereal dgx, doublereal dgy, doublereal dgz, int** &nvtx,
	TOCHKA* &pa, float** prop, int*** &neighbors_for_the_internal_node,
	doublereal eps_minx, doublereal eps_miny, doublereal eps_minz,
	doublereal eps_maxx, doublereal eps_maxy, doublereal eps_maxz,
	VISCOSITY_MODEL &iflowregime, float** &prop_b) {

	int inumcor = number_cores();

#ifdef _OPENMP
	omp_set_num_threads(inumcor); // установка числа потоков
#else 
	inumcor = 1;// один поток.
#endif

	// выделение памяти под искомые полевые величины.
	potent=nullptr;
    potent=new doublereal*[SIZE_FLOW_POTENT_ARRAY];
	if (potent==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for potent constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer i=0; i<SIZE_FLOW_POTENT_ARRAY; ++i) potent[i]=nullptr;
	for (integer i=0; i<SIZE_FLOW_POTENT_ARRAY; ++i) {
		potent[i]=new doublereal[maxelm+maxbound];
		if (potent[i]==nullptr) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem: not enough memory on your equipment for potent[%lld] constr struct...\n", i);
#else
			printf("Problem: not enough memory on your equipment for potent[%d] constr struct...\n", i);
#endif
		    
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}

	// обнуление
	for (integer i = 0; i < SIZE_FLOW_POTENT_ARRAY; ++i) {
		if (i == TURBULENT_KINETIK_ENERGY) {
#pragma omp parallel for
			for (integer j = 0; j < (maxelm + maxbound); ++j) {
				potent[TURBULENT_KINETIK_ENERGY][j] = 4.0e-12;
			}
		}
		else if (i== TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA) {
#pragma omp parallel for
			for (integer j = 0; j < (maxelm + maxbound); ++j) {
				potent[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA][j] = 200.0;
			}
		}
		else {
#pragma omp parallel for
			for (integer j = 0; j < (maxelm + maxbound); ++j) {
				potent[i][j] = 0.0;
			}
		}
	}

	// Паскалевская составляющая давления при естественной конвекции 
	// не вычисляется из уравнений Навье-Стокса а просто добавляется при 
	// инициализации, хотя она влияет на скорости давления в Навье-Стоксе.
	
	doublereal minX1 = 1.0e20;
	doublereal minY1 = 1.0e20;
	doublereal minZ1 = 1.0e20;
	minX1 = 0.0;
	minY1 = 0.0;
	minZ1 = 0.0;

	for(integer iP = 0; iP < maxelm; ++iP) {
		TOCHKA p;
		center_cord3D(iP, nvtx, pa, p,100);
		// вычисление размеров текущего контрольного объёма:
		//doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
		//volume3D(iP, nvtx, pa, dx, dy, dz);
		//if (p.x < minX1) minX1 = p.x-0.5*dx;
		//if (p.y < minY1) minY1 = p.y-0.5*dy;
		//if (p.z < minZ1) minZ1 = p.z-0.5*dz;
		minX1 += p.x;
		minY1 += p.y;
		minZ1 += p.z;
	}
	minX1 /= maxelm;
	minY1 /= maxelm;
	minZ1 /= maxelm;
#pragma omp parallel for
	for (integer iP = 0; iP < maxelm; ++iP) {
		TOCHKA p;
		center_cord3D(iP, nvtx, pa, p,100);
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
		volume3D(iP, nvtx, pa, dx, dy, dz);
		potent[PRESS][iP] += prop[RHO][iP] * (
			dgx*(p.x - minX1) + 
			dgy*(p.y - minY1) +
			dgz*(p.z - minZ1));
		integer iE1=-1, iN1=-1, iT1=-1, iW1=-1, iS1=-1, iB1=-1; // номера соседних контрольных объёмов
		integer iE2=-1, iN2=-1, iT2=-1, iW2=-1, iS2=-1, iB2=-1;
		integer iE3=-1, iN3=-1, iT3=-1, iW3=-1, iS3=-1, iB3=-1;
		integer iE4=-1, iN4=-1, iT4=-1, iW4=-1, iS4=-1, iB4=-1;

		iE1 = neighbors_for_the_internal_node[E_SIDE][0][iP]; 
		iN1 = neighbors_for_the_internal_node[N_SIDE][0][iP];
		iT1 = neighbors_for_the_internal_node[T_SIDE][0][iP]; 
		iW1 = neighbors_for_the_internal_node[W_SIDE][0][iP]; 
		iS1 = neighbors_for_the_internal_node[S_SIDE][0][iP]; 
		iB1 = neighbors_for_the_internal_node[B_SIDE][0][iP]; 

		if (b_on_adaptive_local_refinement_mesh) {
			iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP]; iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP];
			iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP];
			iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP]; iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP]; iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
			iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP];
			iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP];
			iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP]; iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP]; iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];
		}
		/*
		// Так делать ни в коем случае нельзя т.к. возникает нефизичность скорости на границе расчётной области.
		if (0) {
			if (iE >= maxelm) {
				potent[PRESS][iE] += prop[RHO][iP] * (dgx*(p.x + 0.5*dx - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iW >= maxelm) {
				potent[PRESS][iW] += prop[RHO][iP] * (dgx*(p.x - 0.5*dx - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iN >= maxelm) {
				potent[PRESS][iN] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y + 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iS >= maxelm) {
				potent[PRESS][iS] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iT >= maxelm) {
				potent[PRESS][iT] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z + 0.5*dz - minZ1));
			}
			if (iB >= maxelm) {
				potent[PRESS][iB] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z - 0.5*dz - minZ1));
			}
		}
		else */ 
		{
			if (iE1 >= maxelm) {
				potent[PRESS][iE1] += prop[RHO][iP] * (dgx*(p.x + 0.5*dx - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iW1 >= maxelm) {
				potent[PRESS][iW1] += prop[RHO][iP] * (dgx*(p.x - 0.5*dx - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iN1 >= maxelm) {
				potent[PRESS][iN1] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y + 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iS1 >= maxelm) {
				potent[PRESS][iS1] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iT1 >= maxelm) {
				potent[PRESS][iT1] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z + 0.5*dz - minZ1));
			}
			if (iB1 >= maxelm) {
				potent[PRESS][iB1] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z - 0.5*dz - minZ1));
			}

			if (iE2 >= maxelm) {
				potent[PRESS][iE2] += prop[RHO][iP] * (dgx*(p.x + 0.5*dx - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iW2 >= maxelm) {
				potent[PRESS][iW2] += prop[RHO][iP] * (dgx*(p.x - 0.5*dx - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iN2 >= maxelm) {
				potent[PRESS][iN2] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y + 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iS2 >= maxelm) {
				potent[PRESS][iS2] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iT2 >= maxelm) {
				potent[PRESS][iT2] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z + 0.5*dz - minZ1));
			}
			if (iB2 >= maxelm) {
				potent[PRESS][iB2] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z - 0.5*dz - minZ1));
			}

			if (iE3 >= maxelm) {
				potent[PRESS][iE3] += prop[RHO][iP] * (dgx*(p.x + 0.5*dx - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iW3 >= maxelm) {
				potent[PRESS][iW3] += prop[RHO][iP] * (dgx*(p.x - 0.5*dx - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iN3 >= maxelm) {
				potent[PRESS][iN3] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y + 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iS3 >= maxelm) {
				potent[PRESS][iS3] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iT3 >= maxelm) {
				potent[PRESS][iT3] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z + 0.5*dz - minZ1));
			}
			if (iB3 >= maxelm) {
				potent[PRESS][iB3] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z - 0.5*dz - minZ1));
			}

			if (iE4 >= maxelm) {
				potent[PRESS][iE4] += prop[RHO][iP] * (dgx*(p.x + 0.5*dx  - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iW4 >= maxelm) {
				potent[PRESS][iW4] += prop[RHO][iP] * (dgx*(p.x - 0.5*dx - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
			}
			if (iN4 >= maxelm) {
				potent[PRESS][iN4] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y + 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iS4 >= maxelm) {
				potent[PRESS][iS4] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - 0.5*dy - minY1) + dgz*(p.z - minZ1));
			}
			if (iT4 >= maxelm) {
				potent[PRESS][iT4] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z + 0.5*dz - minZ1));
			}
			if (iB4 >= maxelm) {
				potent[PRESS][iB4] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z - 0.5*dz - minZ1));
			}
		}
	}
	
	// 10.02.2017
#pragma omp parallel for
	for (integer iP = 0; iP < maxelm+maxbound; ++iP) {
		//potent[PRESS][iP] *= -1.0;
		// ANSYS Icepak не учитывает гидростатический перепад давления.
		// Для проверки этой гипотезы здесь произведено зануление.
		potent[PRESS][iP] = 0.0;
	}


	/*
	// Давление должно быть инициализировано с учётом граничных условий.
	for (i=0; i<maxbound; ++i) {
		if (((border_neighbor[i].MCB<(ls+lw))&&(border_neighbor[i].MCB>=ls)) && (w[border_neighbor[i].MCB-ls].bpressure)) {
			// на границе задано давление
			potent[PRESS][border_neighbor[i].iB]=w[border_neighbor[i].MCB-ls].P;
		}
	}
	*/


	// Загрузка распределения начальной скорости.
	
	FILE* fp_inicialization_data=NULL;	
#ifdef MINGW_COMPILLER
	int err_inicialization_data = 0;
	fp_inicialization_data=fopen64("load.txt", "r");
	if (fp_inicialization_data == NULL) {
		err_inicialization_data = 1; // Файла несуществует.
	}
#else
	errno_t err_inicialization_data = 0;
    err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif

	if (maxelm > 0) {
		// Только если есть жидкие ячейки.

		if (err_inicialization_data != 0) {
			// файл load.txt отсутствует

			// 04.04.2019
		    // Если мы решаем гидродинамическую задачу то считывать файл скоростей ненадо.

			// открытие неудачно или файл отсутствует.



			// Инициализация компонент скорости во внутренности расчётной области.
			// 26 марта 2017.
			for (integer i = 0; i < maxelm; ++i) {
				potent[VELOCITY_X_COMPONENT][i] = starting_speed_Vx;
				potent[VELOCITY_Y_COMPONENT][i] = starting_speed_Vy;
				potent[VELOCITY_Z_COMPONENT][i] = starting_speed_Vz;
				// скорректированнное поле скорости должно удовлетворять уравнению неразрывности.
				potent[VXCOR][i] = starting_speed_Vx;
				potent[VYCOR][i] = starting_speed_Vy;
				potent[VZCOR][i] = starting_speed_Vz;
			}

		}
		else {
			// файл load.txt присутствует
			if ((fabs(starting_speed_Vx) > 1.0e-20) ||
				(fabs(starting_speed_Vy) > 1.0e-20) ||
				(fabs(starting_speed_Vz) > 1.0e-20)) 
			{
				// Аналитическое задание скоростей приоритетней наличия файла load.txt.

				// Инициализация компонент скорости во внутренности расчётной области.
			    // 26 марта 2017.
				for (integer i = 0; i < maxelm; ++i) {
					potent[VELOCITY_X_COMPONENT][i] = starting_speed_Vx;
					potent[VELOCITY_Y_COMPONENT][i] = starting_speed_Vy;
					potent[VELOCITY_Z_COMPONENT][i] = starting_speed_Vz;
					// скорректированнное поле скорости должно удовлетворять уравнению неразрывности.
					potent[VXCOR][i] = starting_speed_Vx;
					potent[VYCOR][i] = starting_speed_Vy;
					potent[VZCOR][i] = starting_speed_Vz;
				}
			}
			else {

			// 23 июля 2017.
			// Гидродинамические распределения найдены отдельным cfd расчётом и 
			// сохранены в текстовый файл load.txt.
			// Здесь же решается чисто тепловая задача с полем скорости считанным из
			// предварительно подготовленного файла load.txt. 
			// При этом т.к. в тепловой модели сетка изменилась то осуществляется переинтерполяция с 
			// сетки на сетку гидродинамических характеристик с помощью метода наименьших квадратов К.Ф.Гаусса. 

			// удачное открытие 
			printf("load start...\n");

			// Структура файла:
			// maxnode
			// maxelm
			// x
			// y
			// z
			// Vx
			// Vy
			// Vz
			// MUT турбулентная динамическая вязкость.
			// nvtx

			integer maxnode47 = 0;
			integer maxelm47 = 0;
			integer din47 = 0;
#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
			fscanf(fp_inicialization_data, "%lld", &din47);
			maxnode47 = din47;
			fscanf(fp_inicialization_data, "%lld", &din47);
			maxelm47 = din47;
			fscanf(fp_inicialization_data, "%lld", &din47);
			// Режим течения: ламинарный или конкретная модель турбулентности.
			switch (din47) {
			case 0: iflowregime = VISCOSITY_MODEL::LAMINAR;
				break;
			case 1: iflowregime = VISCOSITY_MODEL::ZEROEQMOD;
				break;
			case 2: iflowregime = VISCOSITY_MODEL::SMAGORINSKY;
				break;
			case 3: iflowregime = VISCOSITY_MODEL::RNG_LES;
				break;
			case 4: iflowregime = VISCOSITY_MODEL::RANS_SPALART_ALLMARES;
				break;
			case 5: iflowregime = VISCOSITY_MODEL::RANS_MENTER_SST;
				break;
			case 6: iflowregime = VISCOSITY_MODEL::RANS_STANDART_K_EPS;
				break;
			case 7: iflowregime = VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST;
				break;
			default: iflowregime = VISCOSITY_MODEL::LAMINAR;
				break;
			}
#else
			fscanf(fp_inicialization_data, "%d", &din47);
			maxnode47 = din47;
			fscanf(fp_inicialization_data, "%d", &din47);
			maxelm47 = din47;
			fscanf(fp_inicialization_data, "%d", &din47);
			// Режим течения: ламинарный или конкретная модель турбулентности.
			switch (din47) {
			case 0: iflowregime = VISCOSITY_MODEL::LAMINAR;
				break;
			case 1: iflowregime = VISCOSITY_MODEL::ZEROEQMOD;
				break;
			case 2: iflowregime = VISCOSITY_MODEL::SMAGORINSKY;
				break;
			case 3: iflowregime = VISCOSITY_MODEL::RNG_LES;
				break;
			case 4: iflowregime = VISCOSITY_MODEL::RANS_SPALART_ALLMARES;
				break;
			case 5: iflowregime = VISCOSITY_MODEL::RANS_MENTER_SST;
				break;
			case 6: iflowregime = VISCOSITY_MODEL::RANS_STANDART_K_EPS;
				break;
			case 7: iflowregime = VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST;
				break;
			default: iflowregime = VISCOSITY_MODEL::LAMINAR;
				break;
			}
#endif	
#else
#if doubleintprecision == 1
			fscanf_s(fp_inicialization_data, "%lld", &din47);
			maxnode47 = din47;
			fscanf_s(fp_inicialization_data, "%lld", &din47);
			maxelm47 = din47;
			fscanf_s(fp_inicialization_data, "%lld", &din47);
			// Режим течения: ламинарный или конкретная модель турбулентности.
			switch (din47) {
			case 0: iflowregime = VISCOSITY_MODEL::LAMINAR;
				break;
			case 1: iflowregime = VISCOSITY_MODEL::ZEROEQMOD;
				break;
			case 2: iflowregime = VISCOSITY_MODEL::SMAGORINSKY;
				break;
			case 3: iflowregime = VISCOSITY_MODEL::RNG_LES;
				break;
			case 4: iflowregime = VISCOSITY_MODEL::RANS_SPALART_ALLMARES;
				break;
			case 5: iflowregime = VISCOSITY_MODEL::RANS_MENTER_SST;
				break;
			case 6: iflowregime = VISCOSITY_MODEL::RANS_STANDART_K_EPS;
				break;
			case 7: iflowregime = VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST;
				break;
			default : iflowregime = VISCOSITY_MODEL::LAMINAR;
				break;
			}
			  
#else
			fscanf_s(fp_inicialization_data, "%d", &din47);
			maxnode47 = din47;
			fscanf_s(fp_inicialization_data, "%d", &din47);
			maxelm47 = din47;
			fscanf_s(fp_inicialization_data, "%d", &din47);
			// Режим течения: ламинарный или конкретная модель турбулентности.
			switch (din47) {
			case 0: iflowregime = VISCOSITY_MODEL::LAMINAR;
				break;
			case 1: iflowregime = VISCOSITY_MODEL::ZEROEQMOD;
				break;
			case 2: iflowregime = VISCOSITY_MODEL::SMAGORINSKY;
				break;
			case 3: iflowregime = VISCOSITY_MODEL::RNG_LES;
				break;
			case 4: iflowregime = VISCOSITY_MODEL::RANS_SPALART_ALLMARES;
				break;
			case 5: iflowregime = VISCOSITY_MODEL::RANS_MENTER_SST;
				break;
			case 6: iflowregime = VISCOSITY_MODEL::RANS_STANDART_K_EPS;
				break;
			case 7: iflowregime = VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST;
				break;
			default: iflowregime = VISCOSITY_MODEL::LAMINAR;
				break;
			}
#endif	
#endif
			


			float fin47 = 0.0;
			doublereal *x47 = nullptr,
				       *y47 = nullptr,
				       *z47 = nullptr,
				       *Vx47 = nullptr,
				       *Vy47 = nullptr,
				       *Vz47 = nullptr,
				       *Mut47 = nullptr;
			int** nvtx47 = nullptr;
			// Везде нумерация с нуля.
			x47 = new doublereal[maxnode47];
			y47 = new doublereal[maxnode47];
			z47 = new doublereal[maxnode47];
			Vx47 = new doublereal[maxnode47];
			Vy47 = new doublereal[maxnode47];
			Vz47 = new doublereal[maxnode47];
			Mut47= new doublereal[maxnode47];
			// Считывание данных.
#ifdef MINGW_COMPILLER
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf(fp_inicialization_data, "%f", &fin47);
				x47[i_47] = fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf(fp_inicialization_data, "%f", &fin47);
				y47[i_47] = fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf(fp_inicialization_data, "%f", &fin47);
				z47[i_47] = fin47;
			}

			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf(fp_inicialization_data, "%f", &fin47);
				Vx47[i_47] = my_multiplyer_velocity_load*fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf(fp_inicialization_data, "%f", &fin47);
				Vy47[i_47] = my_multiplyer_velocity_load*fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf(fp_inicialization_data, "%f", &fin47);
				Vz47[i_47] = my_multiplyer_velocity_load*fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf(fp_inicialization_data, "%f", &fin47);
				Mut47[i_47] = fmax(0.0, my_multiplyer_velocity_load*fin47);
			}
#else
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf_s(fp_inicialization_data, "%f", &fin47);
				x47[i_47] = fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf_s(fp_inicialization_data, "%f", &fin47);
				y47[i_47] = fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf_s(fp_inicialization_data, "%f", &fin47);
				z47[i_47] = fin47;
			}

			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf_s(fp_inicialization_data, "%f", &fin47);
				Vx47[i_47] = my_multiplyer_velocity_load*fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf_s(fp_inicialization_data, "%f", &fin47);
				Vy47[i_47] = my_multiplyer_velocity_load*fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf_s(fp_inicialization_data, "%f", &fin47);
				Vz47[i_47] = my_multiplyer_velocity_load*fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf_s(fp_inicialization_data, "%f", &fin47);
				Mut47[i_47] = fmax(0.0, my_multiplyer_velocity_load*fin47);
			}
#endif
			doublereal Speed_min = 1.0e37, Speed_max=-1.0e37;
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				doublereal mag_speed = sqrt((Vx47[i_47])*(Vx47[i_47])+
					                        (Vy47[i_47])*(Vy47[i_47])+
					                        (Vz47[i_47])*(Vz47[i_47]));

				if (mag_speed < Speed_min) {
					Speed_min = mag_speed;
				}
				if (mag_speed > Speed_max) {
					Speed_max = mag_speed;
				}
			}
			printf("LOAD: Minimum speed = %e, Maximum speed = %e \n", Speed_min, Speed_max);

			nvtx47 = new int*[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
			for (integer i_47 = 0; i_47 < NUMBER_OF_VERTEX_FINITE_ELEMENT(); i_47++) {
				nvtx47[i_47] = new int[maxelm47];
			}

			for (integer j_47 = 0; j_47 < maxelm47; j_47++) {
				for (integer i_47 = 0; i_47 < NUMBER_OF_VERTEX_FINITE_ELEMENT(); i_47++) {

					int din48 = -1;

#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
					fscanf(fp_inicialization_data, "%d", &din48);
#else
					fscanf(fp_inicialization_data, "%d", &din48);
#endif
#else
#if doubleintprecision == 1
					fscanf_s(fp_inicialization_data, "%d", &din48);
#else
					fscanf_s(fp_inicialization_data, "%d", &din48);
#endif
#endif


					nvtx47[i_47][j_47] = din48;
					if ((i_47 > 0) && (nvtx47[i_47-1][j_47]==din48)) {
						printf("FILE load.txt is incorrect. \n");
						printf("nvtx is error. Two identical values in a row.\n");
						system("PAUSE");
						exit(1);
					}
				}
			}
			fclose(fp_inicialization_data);

			printf("load done...\n");
			printf("interpolation start...\n");

			// Инициализация компонент скорости во внутренности расчётной области.
			// 26 марта 2017.

			// Метод линейного порядка.
			doublereal min_x = 1e37;
			doublereal min_y = 1e37;
			doublereal min_z = 1e37;
			doublereal max_x = -1e37;
			doublereal max_y = -1e37;
			doublereal max_z = -1e37;

			for (integer i_48 = 0; i_48 < maxelm47; i_48++) {
				for (integer i_47 = 0; i_47 < NUMBER_OF_VERTEX_FINITE_ELEMENT(); i_47++) {
					if (x47[nvtx47[i_47][i_48]] < min_x) {
						min_x = x47[nvtx47[i_47][i_48]];
					}
					if (y47[nvtx47[i_47][i_48]] < min_y) {
						min_y = y47[nvtx47[i_47][i_48]];
					}
					if (z47[nvtx47[i_47][i_48]] < min_z) {
						min_z = z47[nvtx47[i_47][i_48]];
					}
					if (x47[nvtx47[i_47][i_48]] > max_x) {
						max_x = x47[nvtx47[i_47][i_48]];
					}
					if (y47[nvtx47[i_47][i_48]] > max_y) {
						max_y = y47[nvtx47[i_47][i_48]];
					}
					if (z47[nvtx47[i_47][i_48]] > max_z) {
						max_z = z47[nvtx47[i_47][i_48]];
					}
				}
			}

			doublereal min_x1 = min_x;
			doublereal min_y1 = min_y;
			doublereal min_z1 = min_z;
			doublereal max_x1 = max_x;
			doublereal max_y1 = max_y;
			doublereal max_z1 = max_z;

			// Ищем сеточный узел максимально близкий к среднему значению.
			doublereal avgx_separator = 0.5*(min_x+ max_x);
			doublereal avgx_SH = avgx_separator;// avgx_SH - запоминаемая позиция сепаратора avgx_separator.
			{
				doublereal delta = 1.0e37;
				for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
					if (fabs(avgx_separator - x47[i_47]) < delta) {
						delta = fabs(avgx_separator - x47[i_47]);
						avgx_SH = x47[i_47];
					}
				}
			}
			avgx_separator = avgx_SH;

			// Ищем сеточный узел максимально близкий к среднему значению.
			doublereal avgy_separator = 0.5*(min_y + max_y);
			doublereal avgy_SH = avgy_separator;
			{
				doublereal delta = 1.0e37;
				for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
					if (fabs(avgy_separator - y47[i_47]) < delta) {
						delta = fabs(avgy_separator - y47[i_47]);
						avgy_SH = y47[i_47];
					}
				}
			}
			avgy_separator = avgy_SH;

			// Ищем сеточный узел максимально близкий к среднему значению.
			doublereal avgz_separator = 0.5*(min_z + max_z);
			doublereal avgz_SH = avgz_separator;
			{
				doublereal delta = 1.0e37;
				for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
					if (fabs(avgz_separator - z47[i_47]) < delta) {
						delta = fabs(avgz_separator - z47[i_47]);
						avgz_SH = z47[i_47];
					}
				}
			}
			avgz_separator = avgz_SH;

			min_x = 1.05*fabs(max_x - min_x);
			//if (min_x < 1.0e-36) {
				//min_x = 1.05*fabs(max_x);
			//}
			min_y = 1.05*fabs(max_y - min_y);
			//if (min_y < 1.0e-36) {
				//min_y = 1.05*fabs(max_y);
			//}
			min_z = 1.05*fabs(max_z - min_z);
			//if (min_z < 1.0e-36) {
				//min_z = 1.05*fabs(max_z);
			//}

			min_x = min_y= min_z= 0.0;

			TOCHKA** pointerlist = new TOCHKA*[maxelm];
			doublereal** rthdsd_Gauss = new doublereal*[maxelm];
			for (integer i_47 = 0; i_47 < maxelm; i_47++) {
				pointerlist[i_47] = new TOCHKA[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
				rthdsd_Gauss[i_47] = new doublereal[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
			}


			
			integer i_11[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
			integer*** oct_load1 = nullptr;
			doublereal avgx[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			doublereal avgy[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			doublereal avgz[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			integer** i_22 = nullptr;
			i_22 = new integer*[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
			for (integer i_1 = 0; i_1 < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i_1) {
				i_22[i_1] = new integer[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
				for (integer i_2 = 0; i_2 < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i_2) {
					i_22[i_1][i_2] = 0;
				}
			}

			// Возвращает true если внешние габариты расчётной области
            // совпадают с внешними габаритами гидродинамической подобласти,
            // считанной из файла load.txt.
			bool bext_flow=is_EXTERNAL_FLOW(pa, nvtx, maxelm,
				x47, y47, z47, nvtx47, maxelm47,
				eps_minx, eps_miny, eps_minz,
				eps_maxx, eps_maxy, eps_maxz );

			if (bext_flow) printf("EXTERNAL FLOW ACELERATOR is true..\n");

			// Были ли вычислены компоненты скорости на АЛИС ? 22,03,2019
			bool b_load_has_been_saved_on_ALICE = !bext_flow;
			if (b_load_has_been_saved_on_ALICE) {
				if ((maxelm > 6000000) || (maxelm47 > 100000)) {
					if (!b_on_adaptive_local_refinement_mesh) {
						// Сетка структурированная, число тепловых узлов
						// более 8М, число гидродинамических узлов более 100К
						// Т.е. в случае очень больших сеток мы пытаемся ускорится
						// с помощью octree дерева.
						b_load_has_been_saved_on_ALICE = false;
					}
				}			
			}

			// Был активен запрещающий оператор до 03.09.2020
			//b_load_has_been_saved_on_ALICE = true;// true ускорение не применяется.

			if (!b_load_has_been_saved_on_ALICE) {

				printf("maxelm CFD=%lld \n", maxelm47);	
				integer** oct_load = nullptr;
				oct_load = new integer*[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
				if (oct_load == nullptr) {
					printf("ERROR ALLOCATION MEMORY: array oct_load is nullptr.");
					system("PAUSE");
				}
				for (integer i = 0; i < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) {
					//oct_load[i] = nullptr;
					oct_load[i] = new integer[maxelm47];
					if (oct_load[i] == nullptr) {
						printf("ERROR ALLOCATION MEMORY: array oct_load[%lld] is nullptr.", i);
						system("PAUSE");
					}
					for (integer i_1 = 0; i_1 < maxelm47; ++i_1) {
						oct_load[i][i_1] = -1;
					}
				}
				
				for (integer i_47 = 0; i_47 < maxelm47; i_47++) {
					// Вычисляем координаты центра ячейки i47.
					doublereal xc47 = 0.0, yc47 = 0.0, zc47 = 0.0;
					for (integer i_9 = 0; i_9 < NUMBER_OF_VERTEX_FINITE_ELEMENT(); i_9++) {
						xc47 += 0.125*(x47[nvtx47[i_9][i_47]]);
						yc47 += 0.125*(y47[nvtx47[i_9][i_47]]);
						zc47 += 0.125*(z47[nvtx47[i_9][i_47]]);
					}


					// Сортируем:
					if (xc47 < avgx_separator) {
						if (yc47 < avgy_separator) {
							if (zc47 < avgz_separator) {
								oct_load[0][i_11[0]] = i_47;
								i_11[0]++;
							}
							else {
								oct_load[1][i_11[1]] = i_47;
								i_11[1]++;
							}
						}
						else {
							if (zc47 < avgz_separator) {
								oct_load[2][i_11[2]] = i_47;
								i_11[2]++;
							}
							else {
								oct_load[3][i_11[3]] = i_47;
								i_11[3]++;
							}
						}
					}
					else {
						if (yc47 < avgy_separator) {
							if (zc47 < avgz_separator) {
								oct_load[4][i_11[4]] = i_47;
								i_11[4]++;
							}
							else {
								oct_load[5][i_11[5]] = i_47;
								i_11[5]++;
							}
						}
						else {
							if (zc47 < avgz_separator) {
								oct_load[6][i_11[6]] = i_47;
								i_11[6]++;
							}
							else {
								oct_load[7][i_11[7]] = i_47;
								i_11[7]++;
							}
						}
					}
				}
				// Для задач внешнего обтекания одного уровня octree дерева мало.
				// Чтобы ещё ускорить моделирование нужен ещё один уровень octree дерева.
				
				oct_load1 = new integer**[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
				// Если new не может выделить оперативную память он направляет исключение.
				//if (oct_load1 == nullptr) {				
					//printf("ERROR ALLOCATION MEMORY: array oct_load1 is nullptr.");
					//system("PAUSE");
				//}
				for (integer i = 0; i < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) {
					//oct_load1[i] = nullptr;
					oct_load1[i] = new integer*[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
					//if (oct_load1 == nullptr) {
						//printf("ERROR ALLOCATION MEMORY: array oct_load1 is nullptr.");
						//system("PAUSE");
					//}
					for (integer j = 0; j < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++j) {
						oct_load1[i][j] = new integer[i_11[i]];
						if (oct_load1[i][j] == nullptr) {
							printf("ERROR ALLOCATION MEMORY: array oct_load1[%lld][%lld] is nullptr.", i, j);
							system("PAUSE");
						}

						for (integer i_1 = 0; i_1 < i_11[i]; ++i_1) {
							oct_load1[i][j][i_1] = -1;
						}
					}
				}

				printf("balance octree level 1: %e %e %e %e %e %e %e %e\n", (double)(i_11[0]) / (double)(maxelm47 + 1), (double)(i_11[1]) / (double)(maxelm47 + 1), (double)(i_11[2]) / (double)(maxelm47 + 1), (double)(i_11[3]) / (double)(maxelm47 + 1), (double)(i_11[4]) / (double)(maxelm47 + 1), (double)(i_11[5]) / (double)(maxelm47 + 1), (double)(i_11[6]) / (double)(maxelm47 + 1), (double)(i_11[7]) / (double)(maxelm47 + 1));


				
				for (integer i = 0; i < NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) {

					// Разбиение максимально близкое 
					// к способу построения octree дерева
					switch (i) {
					case 0: avgx[i] = min_x1+0.25*fabs(max_x1-min_x1);
						avgy[i] = min_y1 + 0.25*fabs(max_y1-min_y1);
						avgz[i] = min_z1 + 0.25*fabs(max_z1-min_z1);
						break;
					case 1: avgx[i] = min_x1 + 0.25* fabs(max_x1 - min_x1);
						avgy[i] = min_y1 + 0.25* fabs(max_y1 - min_y1);
						avgz[i] = min_z1 + 0.75* fabs(max_z1 - min_z1);
						break;
					case 2: avgx[i] = min_x1 + 0.25* fabs(max_x1 - min_x1);
						avgy[i] = min_y1 + 0.75* fabs(max_y1 - min_y1);
						avgz[i] = min_z1 + 0.25* fabs(max_z1 - min_z1);
						break;
					case 3: avgx[i] = min_x1 + 0.25* fabs(max_x1 - min_x1);;
						avgy[i] = min_y1 + 0.75* fabs(max_y1 - min_y1);
						avgz[i] = min_z1 + 0.75* fabs(max_z1 - min_z1);
						break;
					case 4: avgx[i] = min_x1 + 0.75* fabs(max_x1 - min_x1);;
						avgy[i] = min_y1 + 0.25* fabs(max_y1 - min_y1);
						avgz[i] = min_z1 + 0.25* fabs(max_z1 - min_z1);
						break;
					case 5: avgx[i] = min_x1 + 0.75* fabs(max_x1 - min_x1);;
						avgy[i] = min_y1 + 0.25* fabs(max_y1 - min_y1);
						avgz[i] = min_z1 + 0.75* fabs(max_z1 - min_z1);
						break;
					case 6: avgx[i] = min_x1 + 0.75* fabs(max_x1 - min_x1);;
						avgy[i] = min_y1 + 0.75* fabs(max_y1 - min_y1);
						avgz[i] = min_z1 + 0.25* fabs(max_z1 - min_z1);
						break;
					case 7: avgx[i] = min_x1 + 0.75* fabs(max_x1 - min_x1);;
						avgy[i] = min_y1 + 0.75* fabs(max_y1 - min_y1);
						avgz[i] = min_z1 + 0.75* fabs(max_z1 - min_z1);
						break;
					}

					/*
					// Новое разбиение.
					for (integer i_1 = 0; i_1 < i_11[i]; ++i_1) {
						for (integer i_9 = 0; i_9 < 8; i_9++) {
							avgx[i] += (0.125*x47[nvtx47[i_9][oct_load[i][i_1]]]) / i_11[i];
							avgy[i] += (0.125*y47[nvtx47[i_9][oct_load[i][i_1]]]) / i_11[i];
							avgz[i] += (0.125*z47[nvtx47[i_9][oct_load[i][i_1]]]) / i_11[i];
						}
					}
					*/

					// 23.03.2019
					// Ищем сеточный узел максимально близкий к среднему значению.
					doublereal avgx_separator1 = avgx[i];
					doublereal avgx_SH1 = avgx_separator1;
					{
						doublereal delta = 1.0e37;
						for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
							if (fabs(avgx_separator1 - x47[i_47]) < delta) {
								delta = fabs(avgx_separator1 - x47[i_47]);
								avgx_SH1 = x47[i_47];
							}
						}
					}
					avgx_separator1 = avgx_SH1;
					avgx[i] = avgx_separator1; // Мы точно попали на сеточную линию.

					// Ищем сеточный узел максимально близкий к среднему значению.
					doublereal avgy_separator1 = avgy[i];
					doublereal avgy_SH1 = avgy_separator1;
					{
						doublereal delta = 1.0e37;
						for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
							if (fabs(avgy_separator1 - y47[i_47]) < delta) {
								delta = fabs(avgy_separator1 - y47[i_47]);
								avgy_SH1 = y47[i_47];
							}
						}
					}
					avgy_separator1 = avgy_SH1;
					avgy[i] = avgy_separator1; // Мы точно попали на сеточную линию.

											   // Ищем сеточный узел максимально близкий к среднему значению.
					doublereal avgz_separator1 = avgz[i];
					doublereal avgz_SH1 = avgz_separator1;
					{
						doublereal delta = 1.0e37;
						for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
							if (fabs(avgz_separator1 - z47[i_47]) < delta) {
								delta = fabs(avgz_separator1 - z47[i_47]);
								avgz_SH1 = z47[i_47];
							}
						}
					}
					avgz_separator1 = avgz_SH1;
					avgz[i] = avgz_separator1; // Мы точно попали на сеточную линию.

				}

				for (integer i = 0; i < 8; ++i) {
					delete[] oct_load[i];
					oct_load[i] = nullptr;
				}
				delete[] oct_load;
				oct_load = nullptr;

				for (integer i_47 = 0; i_47 < maxelm47; i_47++) {
					// Вычисляем координаты центра ячейки i47.
					doublereal xc47 = 0.0, yc47 = 0.0, zc47 = 0.0;
					for (integer i_9 = 0; i_9 < 8; i_9++) {
						xc47 += 0.125*(x47[nvtx47[i_9][i_47]]);
						yc47 += 0.125*(y47[nvtx47[i_9][i_47]]);
						zc47 += 0.125*(z47[nvtx47[i_9][i_47]]);
					}
					// Сортируем:
					if (xc47 < avgx_separator) {
						if (yc47 < avgy_separator) {
							if (zc47 < avgz_separator) {
								//oct_load[0][i_11[0]] = i_47;
								//i_11[0]++;
								// Распределяем по octree дереву:
								//integer identifikator = 0;
								CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(0, xc47, yc47, zc47, i_47, oct_load1, i_22, avgx, avgy, avgz);
							}
							else {
								//oct_load[1][i_11[1]] = i_47;
								//i_11[1]++;
								// Распределяем по octree дереву:
								//integer identifikator = 1;
								CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(1, xc47, yc47, zc47, i_47, oct_load1, i_22, avgx, avgy, avgz);
							}
						}
						else {
							if (zc47 < avgz_separator) {
								//oct_load[2][i_11[2]] = i_47;
								//i_11[2]++;
								// Распределяем по octree дереву:
								//integer identifikator = 2;
								CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(2, xc47, yc47, zc47, i_47, oct_load1, i_22, avgx, avgy, avgz);
							}
							else {
								//oct_load[3][i_11[3]] = i_47;
								//i_11[3]++;
								// Распределяем по octree дереву:
								//integer identifikator = 3;
								CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(3, xc47, yc47, zc47, i_47, oct_load1, i_22, avgx, avgy, avgz);
							}
						}
					}
					else {
						if (yc47 < avgy_separator) {
							if (zc47 < avgz_separator) {
								//oct_load[4][i_11[4]] = i_47;
								//i_11[4]++;
								// Распределяем по octree дереву:
								//integer identifikator = 4;
								CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(4, xc47, yc47, zc47, i_47, oct_load1, i_22, avgx, avgy, avgz);
							}
							else {
								//oct_load[5][i_11[5]] = i_47;
								//i_11[5]++;
								// Распределяем по octree дереву:
								//integer identifikator = 5;
								CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(5, xc47, yc47, zc47, i_47, oct_load1, i_22, avgx, avgy, avgz);
							}
						}
						else {
							if (zc47 < avgz_separator) {
								//oct_load[6][i_11[6]] = i_47;
								//i_11[6]++;
								// Распределяем по octree дереву:
								//integer identifikator = 6;
								CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(6, xc47, yc47, zc47, i_47, oct_load1, i_22, avgx, avgy, avgz);
							}
							else {
								//oct_load[7][i_11[7]] = i_47;
								//i_11[7]++;
								// Распределяем по octree дереву:
								//integer identifikator = 7;
								CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD(7, xc47, yc47, zc47, i_47, oct_load1, i_22, avgx, avgy, avgz);
							}
						}
					}
				}
			}

			printf("balance 8 octree level 2 is construct: Ok\n");

			if (b_load_has_been_saved_on_ALICE) {
				//printf("b_load_has_been_saved_on_ALICE\n");
				//system("pause");

				bool* bmask = new bool[maxelm];
				for (integer i = 0; i < maxelm; ++i) {
					bmask[i] = false;
				}


				for (integer i = 0; i < maxelm; ++i) {
					doublereal xc47, yc47, zc47;

					TOCHKA p;
					center_cord3D(i, nvtx, pa, p, 100);
					xc47 = p.x;
					yc47 = p.y;
					zc47 = p.z;


					// Найдем номер конечного изопараметрического элемента 
			        // в котором содержится точка p.
					bool bfound = false;
					integer ifound = -1;

					// Только если точка лежит внутри жидкой зоны мы запускаем поиск её местоположения.
					// Эффективность этого предиктора сильно снижена если жидких зон несколько.
					if ((xc47 >= min_x1) && (xc47 <= max_x1) &&
						(yc47 >= min_y1) && (yc47 <= max_y1) &&
						(zc47 >= min_z1) && (zc47 <= max_z1))
					{



						// 22.03.2019 Для АЛИС сетки мы не можем разбить на три уровня octree дерева так как 
						// границы ячеек сетки не вдоль координатных линий а уступами, ступеньками сильно неравномерные
						// отклоняющиеся от линии разбиения тов меньшую то в большую сторону.
#pragma omp parallel
						{

							bool bfound_loc = false;
							integer ifound_loc = -1;

#pragma omp for 
							for (integer i_47 = 0; i_47 < maxelm47; i_47++) {
								if (bfound_loc) continue;

								if (x47[nvtx47[0][i_47]] >= x47[nvtx47[1][i_47]]) {
									printf("ERROR logic X !!! if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]])\n");
									//printf("0=%e 1=%e %lld %lld\n", x47[nvtx47[0][i_47]], x47[nvtx47[1][i_47]], nvtx47[0][i_47], nvtx47[1][i_47]);
									std::cout << "0=" << x47[nvtx47[0][i_47]] << " 1=" << x47[nvtx47[1][i_47]] << " " << nvtx47[0][i_47] << " " << nvtx47[1][i_47] << std::endl;

									printf("in function allocation memory flow!!! constr_struct.cpp\n");
									system("PAUSE");
								}
								if (y47[nvtx47[0][i_47]] >= y47[nvtx47[3][i_47]]) {
									printf("ERROR logic Y if ((yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) !!!\n");
									printf("0=%e 3=%e\n", y47[nvtx47[0][i_47]], y47[nvtx47[3][i_47]]);
									printf("in function allocation memory flow!!! constr_struct.cpp\n");
									system("PAUSE");
								}
								if (z47[nvtx47[0][i_47]] >= z47[nvtx47[4][i_47]]) {
									printf("ERROR logic Z if (zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))!!!\n");
									printf("0=%e 4=%e\n", z47[nvtx47[0][i_47]], z47[nvtx47[4][i_47]]);
									printf("in function allocation memory flow!!! constr_struct.cpp\n");
									system("PAUSE");
								}
								if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) &&
									(yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
									(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
								{
									ifound_loc = i_47;
									bfound_loc = true;
									//break;
								}
							}

#pragma omp critical
							{
								if (bfound_loc) {
									bfound = true;
									ifound = ifound_loc;
								}
							}
						}
					}
				
					// МНК
					// Три раза, для каждой компоненты скорости.
					doublereal starting_speed_Vx_47 = 0.0;
					doublereal starting_speed_Vy_47 = 0.0;
					doublereal starting_speed_Vz_47 = 0.0;
					doublereal starting_Mut_47 = 0.0;

					if (bfound) {

						bmask[i] = true;

						//TOCHKA p;
						//center_cord3D(i, nvtx, pa, p, 100);
						p.x = p.x + min_x;
						p.y = p.y + min_y;
						p.z = p.z + min_z;

						// VX

						for (integer j = 0; j <= 7; ++j) {
							TOCHKA p1;
							p1.x = x47[nvtx47[j][ifound]];
							p1.y = y47[nvtx47[j][ifound]];
							p1.z = z47[nvtx47[j][ifound]];
							p1.x = p1.x + min_x;
							p1.y = p1.y + min_y;
							p1.z = p1.z + min_z;

							pointerlist[i][j] = p1;
							if (fabs(p1.x) < -1.0e-36) {
								printf("problem x=%e\n", p1.x);
								system("PAUSE");
							}
							if (fabs(p1.y) < -1.0e-36) {
								printf("problem y=%e\n", p1.y);
								system("PAUSE");
							}
							if (fabs(p1.z) < -1.0e-36) {
								printf("problem z=%e\n", p1.z);
								system("PAUSE");
							}
							rthdsd_Gauss[i][j] = Vx47[nvtx47[j][ifound]];
						}

						doublereal** Xmatr = new doublereal * [4];
						for (integer j = 0; j <= 3; ++j) {
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




						for (integer j = 0; j < 8; ++j) {

							Xmatr[0][0] += 1.0;
							Xmatr[0][1] += pointerlist[i][j].x;
							Xmatr[0][2] += pointerlist[i][j].y;
							Xmatr[0][3] += pointerlist[i][j].z;

							Xmatr[1][0] += pointerlist[i][j].x;
							Xmatr[1][1] += pointerlist[i][j].x * pointerlist[i][j].x;
							Xmatr[1][2] += pointerlist[i][j].x * pointerlist[i][j].y;
							Xmatr[1][3] += pointerlist[i][j].x * pointerlist[i][j].z;

							Xmatr[2][0] += pointerlist[i][j].y;
							Xmatr[2][1] += pointerlist[i][j].y * pointerlist[i][j].x;
							Xmatr[2][2] += pointerlist[i][j].y * pointerlist[i][j].y;
							Xmatr[2][3] += pointerlist[i][j].y * pointerlist[i][j].z;

							Xmatr[3][0] += pointerlist[i][j].z;
							Xmatr[3][1] += pointerlist[i][j].z * pointerlist[i][j].x;
							Xmatr[3][2] += pointerlist[i][j].z * pointerlist[i][j].y;
							Xmatr[3][3] += pointerlist[i][j].z * pointerlist[i][j].z;

							bmatr[0] += rthdsd_Gauss[i][j];
							bmatr[1] += pointerlist[i][j].x * rthdsd_Gauss[i][j];
							bmatr[2] += pointerlist[i][j].y * rthdsd_Gauss[i][j];
							bmatr[3] += pointerlist[i][j].z * rthdsd_Gauss[i][j];
						}


						for (integer j1 = 0; j1 <= 100; j1++) {
							koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
							koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
							koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
							koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
						}

						starting_speed_Vx_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));

						// VY

						for (integer j = 0; j <= 7; ++j) {
							rthdsd_Gauss[i][j] = Vy47[nvtx47[j][ifound]];
						}

						for (integer j1 = 0; j1 <= 3; j1++) {
							for (integer j2 = 0; j2 <= 3; j2++) {
								Xmatr[j1][j2] = 0.0;
							}
							bmatr[j1] = 0.0;
							koefmatr[j1] = 0.0;
						}




						for (integer j = 0; j < 8; ++j) {

							Xmatr[0][0] += 1.0;
							Xmatr[0][1] += pointerlist[i][j].x;
							Xmatr[0][2] += pointerlist[i][j].y;
							Xmatr[0][3] += pointerlist[i][j].z;

							Xmatr[1][0] += pointerlist[i][j].x;
							Xmatr[1][1] += pointerlist[i][j].x * pointerlist[i][j].x;
							Xmatr[1][2] += pointerlist[i][j].x * pointerlist[i][j].y;
							Xmatr[1][3] += pointerlist[i][j].x * pointerlist[i][j].z;

							Xmatr[2][0] += pointerlist[i][j].y;
							Xmatr[2][1] += pointerlist[i][j].y * pointerlist[i][j].x;
							Xmatr[2][2] += pointerlist[i][j].y * pointerlist[i][j].y;
							Xmatr[2][3] += pointerlist[i][j].y * pointerlist[i][j].z;

							Xmatr[3][0] += pointerlist[i][j].z;
							Xmatr[3][1] += pointerlist[i][j].z * pointerlist[i][j].x;
							Xmatr[3][2] += pointerlist[i][j].z * pointerlist[i][j].y;
							Xmatr[3][3] += pointerlist[i][j].z * pointerlist[i][j].z;

							bmatr[0] += rthdsd_Gauss[i][j];
							bmatr[1] += pointerlist[i][j].x * rthdsd_Gauss[i][j];
							bmatr[2] += pointerlist[i][j].y * rthdsd_Gauss[i][j];
							bmatr[3] += pointerlist[i][j].z * rthdsd_Gauss[i][j];
						}


						for (integer j1 = 0; j1 <= 100; j1++) {
							koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
							koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
							koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
							koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
						}

						starting_speed_Vy_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));


						// VZ

						for (integer j = 0; j <= 7; ++j) {
							rthdsd_Gauss[i][j] = Vz47[nvtx47[j][ifound]];
						}

						for (integer j1 = 0; j1 <= 3; j1++) {
							for (integer j2 = 0; j2 <= 3; j2++) {
								Xmatr[j1][j2] = 0.0;
							}
							bmatr[j1] = 0.0;
							koefmatr[j1] = 0.0;
						}




						for (integer j = 0; j < 8; ++j) {

							Xmatr[0][0] += 1.0;
							Xmatr[0][1] += pointerlist[i][j].x;
							Xmatr[0][2] += pointerlist[i][j].y;
							Xmatr[0][3] += pointerlist[i][j].z;

							Xmatr[1][0] += pointerlist[i][j].x;
							Xmatr[1][1] += pointerlist[i][j].x * pointerlist[i][j].x;
							Xmatr[1][2] += pointerlist[i][j].x * pointerlist[i][j].y;
							Xmatr[1][3] += pointerlist[i][j].x * pointerlist[i][j].z;

							Xmatr[2][0] += pointerlist[i][j].y;
							Xmatr[2][1] += pointerlist[i][j].y * pointerlist[i][j].x;
							Xmatr[2][2] += pointerlist[i][j].y * pointerlist[i][j].y;
							Xmatr[2][3] += pointerlist[i][j].y * pointerlist[i][j].z;

							Xmatr[3][0] += pointerlist[i][j].z;
							Xmatr[3][1] += pointerlist[i][j].z * pointerlist[i][j].x;
							Xmatr[3][2] += pointerlist[i][j].z * pointerlist[i][j].y;
							Xmatr[3][3] += pointerlist[i][j].z * pointerlist[i][j].z;

							bmatr[0] += rthdsd_Gauss[i][j];
							bmatr[1] += pointerlist[i][j].x * rthdsd_Gauss[i][j];
							bmatr[2] += pointerlist[i][j].y * rthdsd_Gauss[i][j];
							bmatr[3] += pointerlist[i][j].z * rthdsd_Gauss[i][j];
						}


						for (integer j1 = 0; j1 <= 100; j1++) {
							koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
							koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
							koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
							koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
						}

						starting_speed_Vz_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));


						// Mut

						for (integer j = 0; j <= 7; ++j) {
							rthdsd_Gauss[i][j] = Mut47[nvtx47[j][ifound]];
						}

						for (integer j1 = 0; j1 <= 3; j1++) {
							for (integer j2 = 0; j2 <= 3; j2++) {
								Xmatr[j1][j2] = 0.0;
							}
							bmatr[j1] = 0.0;
							koefmatr[j1] = 0.0;
						}




						for (integer j = 0; j < 8; ++j) {

							Xmatr[0][0] += 1.0;
							Xmatr[0][1] += pointerlist[i][j].x;
							Xmatr[0][2] += pointerlist[i][j].y;
							Xmatr[0][3] += pointerlist[i][j].z;

							Xmatr[1][0] += pointerlist[i][j].x;
							Xmatr[1][1] += pointerlist[i][j].x * pointerlist[i][j].x;
							Xmatr[1][2] += pointerlist[i][j].x * pointerlist[i][j].y;
							Xmatr[1][3] += pointerlist[i][j].x * pointerlist[i][j].z;

							Xmatr[2][0] += pointerlist[i][j].y;
							Xmatr[2][1] += pointerlist[i][j].y * pointerlist[i][j].x;
							Xmatr[2][2] += pointerlist[i][j].y * pointerlist[i][j].y;
							Xmatr[2][3] += pointerlist[i][j].y * pointerlist[i][j].z;

							Xmatr[3][0] += pointerlist[i][j].z;
							Xmatr[3][1] += pointerlist[i][j].z * pointerlist[i][j].x;
							Xmatr[3][2] += pointerlist[i][j].z * pointerlist[i][j].y;
							Xmatr[3][3] += pointerlist[i][j].z * pointerlist[i][j].z;

							bmatr[0] += rthdsd_Gauss[i][j];
							bmatr[1] += pointerlist[i][j].x * rthdsd_Gauss[i][j];
							bmatr[2] += pointerlist[i][j].y * rthdsd_Gauss[i][j];
							bmatr[3] += pointerlist[i][j].z * rthdsd_Gauss[i][j];
						}


						for (integer j1 = 0; j1 <= 100; j1++) {
							koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
							koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
							koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
							koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
						}

						starting_Mut_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));



						for (integer j = 0; j <= 3; ++j) {
							delete[] Xmatr[j];
						}
						delete[] Xmatr;
						delete[] bmatr;
						delete[] koefmatr;
					}

					if (starting_speed_Vx_47 != starting_speed_Vx_47) {
						printf("NAN pri interpolation Vx %lld in allocation_memory_flow in constr_struct.cpp\n", i);
					}
					if (starting_speed_Vy_47 != starting_speed_Vy_47) {
						printf("NAN pri interpolation Vy %lld in allocation_memory_flow in constr_struct.cpp\n", i);
					}
					if (starting_speed_Vz_47 != starting_speed_Vz_47) {
						printf("NAN pri interpolation Vz %lld in allocation_memory_flow in constr_struct.cpp\n", i);
					}
					if (starting_Mut_47 != starting_Mut_47) {
						printf("NAN pri interpolation Mut %lld in allocation_memory_flow in constr_struct.cpp\n", i);
					}

					potent[VELOCITY_X_COMPONENT][i] = starting_speed_Vx_47;
					potent[VELOCITY_Y_COMPONENT][i] = starting_speed_Vy_47;
					potent[VELOCITY_Z_COMPONENT][i] = starting_speed_Vz_47;
					potent[MUT][i] = fmax(0.0,starting_Mut_47);
					// скорректированнное поле скорости должно удовлетворять уравнению неразрывности.
					potent[VXCOR][i] = starting_speed_Vx_47;
					potent[VYCOR][i] = starting_speed_Vy_47;
					potent[VZCOR][i] = starting_speed_Vz_47;

					if ((i!=0)&&(i % 100 == 0)) {
						// Хронология выполнения. Сколько осталось до конца.
						printf("\b\b\b%2d%%", (int)(100 * (maxelm - i) / maxelm));
					}
				
				}

				// По видимому сохранённый периметр гидродинамической подобласти немного уже чем реально начерченная
				// гидродинамическая подобласть. Образуется что то типо щели в приграничной области на стыке сохранеённой гидродинамической
				// области и реально начерченной. В этой щели нулевая скорость. Щель с нулевой скоростью приводит к неверным результатам.
				// Здесь мы смыкаем данную щель заполняя ненулевой скоростью по непрерывности с основной считанной гидродинамичсеской
				// подобластью. 16.12.2021

				bool  bideptch_15_12_2021 = true;

				while (bideptch_15_12_2021) {

					bideptch_15_12_2021=false;

					for (integer i = 0; i < maxelm; ++i) {
						for (int G = 0; G < 6; ++G) {
							if ((bmask[i]) && ((neighbors_for_the_internal_node[G][0][i] >= maxelm) ||
								((neighbors_for_the_internal_node[G][0][i] < maxelm) && (!bmask[neighbors_for_the_internal_node[G][0][i]])))) {

								potent[VELOCITY_X_COMPONENT][neighbors_for_the_internal_node[G][0][i]] = potent[VELOCITY_X_COMPONENT][i];
								potent[VELOCITY_Y_COMPONENT][neighbors_for_the_internal_node[G][0][i]] = potent[VELOCITY_Y_COMPONENT][i];
								potent[VELOCITY_Z_COMPONENT][neighbors_for_the_internal_node[G][0][i]] = potent[VELOCITY_Z_COMPONENT][i];
								potent[MUT][neighbors_for_the_internal_node[G][0][i]] = potent[MUT][i];
								// скорректированнное поле скорости должно удовлетворять уравнению неразрывности.
								potent[VXCOR][neighbors_for_the_internal_node[G][0][i]] = potent[VXCOR][i];
								potent[VYCOR][neighbors_for_the_internal_node[G][0][i]] = potent[VYCOR][i];
								potent[VZCOR][neighbors_for_the_internal_node[G][0][i]] = potent[VZCOR][i];
								if (neighbors_for_the_internal_node[G][0][i] < maxelm) {
									bmask[neighbors_for_the_internal_node[G][0][i]] = true;
									bideptch_15_12_2021 = true;
								}
							}

						}
					}
				}

				delete[] bmask;
			}
			else {
				//printf("not b_load_has_been_saved_on_ALICE\n");
				//system("pause");

				bool* bmask = new bool[maxelm];
				for (integer i = 0; i < maxelm; ++i) {
					bmask[i] = false;
				}

				for (integer i = 0; i < maxelm; ++i) {
					doublereal xc47, yc47, zc47;

					TOCHKA p;
					center_cord3D(i, nvtx, pa, p, 100);
					xc47 = p.x;
					yc47 = p.y;
					zc47 = p.z;


					// Найдем номер конечного изопараметрического элемента 
					// в котором содержится точка p.
					
					bool bfound = false;
					integer ifound = -1;

					// Только если точка лежит внутри жидкой зоны мы запускаем поиск её местоположения.
					// Эффективность этого предиктора сильно снижена если жидких зон несколько.
					if ((xc47 >= min_x1) && (xc47 <= max_x1) &&
						(yc47 >= min_y1) && (yc47 <= max_y1) &&
						(zc47 >= min_z1) && (zc47 <= max_z1))
					{



						// Редуцируем (уменьшим) сложность поиска в 8 раз.
						// За счёт одного уровня octree дерева.
						if (xc47 < avgx_separator) {
							if (yc47 < avgy_separator) {
								if (zc47 < avgz_separator) {

									CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(0, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);

								}
								else {

									CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(1, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);

								}
							}
							else {
								if (zc47 < avgz_separator) {

									CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(2, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);

								}
								else {

									CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(3, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);

								}
							}
						}
						else {
							if (yc47 < avgy_separator) {
								if (zc47 < avgz_separator) {

									CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(4, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);

								}
								else {

									CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(5, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);

								}
							}
							else {
								if (zc47 < avgz_separator) {

									CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(6, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);

								}
								else {

									CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(7, xc47, yc47, zc47, oct_load1, i_22, avgx, avgy, avgz, x47, y47, z47, nvtx47, bfound, ifound);

								}
							}
						}


					}
				
				
					// МНК
						// Три раза, для каждой компоненты скорости.
					doublereal starting_speed_Vx_47 = 0.0;
					doublereal starting_speed_Vy_47 = 0.0;
					doublereal starting_speed_Vz_47 = 0.0;
					doublereal starting_Mut_47 = 0.0;

					if (bfound) {

						bmask[i] = true;

						//TOCHKA p;
						//center_cord3D(i, nvtx, pa, p, 100);
						p.x = p.x + min_x;
						p.y = p.y + min_y;
						p.z = p.z + min_z;

						// VX

						for (integer j = 0; j <= 7; ++j) {
							TOCHKA p1;
							p1.x = x47[nvtx47[j][ifound]];
							p1.y = y47[nvtx47[j][ifound]];
							p1.z = z47[nvtx47[j][ifound]];
							p1.x = p1.x + min_x;
							p1.y = p1.y + min_y;
							p1.z = p1.z + min_z;

							pointerlist[i][j] = p1;
							if (fabs(p1.x) < -1.0e-36) {
								printf("problem x=%e\n", p1.x);
								system("PAUSE");
							}
							if (fabs(p1.y) < -1.0e-36) {
								printf("problem y=%e\n", p1.y);
								system("PAUSE");
							}
							if (fabs(p1.z) < -1.0e-36) {
								printf("problem z=%e\n", p1.z);
								system("PAUSE");
							}
							rthdsd_Gauss[i][j] = Vx47[nvtx47[j][ifound]];
						}

						doublereal** Xmatr = new doublereal * [4];
						for (integer j = 0; j <= 3; ++j) {
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




						for (integer j = 0; j < 8; ++j) {

							Xmatr[0][0] += 1.0;
							Xmatr[0][1] += pointerlist[i][j].x;
							Xmatr[0][2] += pointerlist[i][j].y;
							Xmatr[0][3] += pointerlist[i][j].z;

							Xmatr[1][0] += pointerlist[i][j].x;
							Xmatr[1][1] += pointerlist[i][j].x * pointerlist[i][j].x;
							Xmatr[1][2] += pointerlist[i][j].x * pointerlist[i][j].y;
							Xmatr[1][3] += pointerlist[i][j].x * pointerlist[i][j].z;

							Xmatr[2][0] += pointerlist[i][j].y;
							Xmatr[2][1] += pointerlist[i][j].y * pointerlist[i][j].x;
							Xmatr[2][2] += pointerlist[i][j].y * pointerlist[i][j].y;
							Xmatr[2][3] += pointerlist[i][j].y * pointerlist[i][j].z;

							Xmatr[3][0] += pointerlist[i][j].z;
							Xmatr[3][1] += pointerlist[i][j].z * pointerlist[i][j].x;
							Xmatr[3][2] += pointerlist[i][j].z * pointerlist[i][j].y;
							Xmatr[3][3] += pointerlist[i][j].z * pointerlist[i][j].z;

							bmatr[0] += rthdsd_Gauss[i][j];
							bmatr[1] += pointerlist[i][j].x * rthdsd_Gauss[i][j];
							bmatr[2] += pointerlist[i][j].y * rthdsd_Gauss[i][j];
							bmatr[3] += pointerlist[i][j].z * rthdsd_Gauss[i][j];
						}


						for (integer j1 = 0; j1 <= 16; j1++) {
#pragma omp parallel sections 
{
#pragma omp section
	{
		koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
	}
#pragma omp section
	{
		koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
	}
#pragma omp section
	{
		koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
	}
#pragma omp section
	{
		koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
	}							
}
						}

						starting_speed_Vx_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));

						// VY

						for (integer j = 0; j <= 7; ++j) {
							rthdsd_Gauss[i][j] = Vy47[nvtx47[j][ifound]];
						}

						for (integer j1 = 0; j1 <= 3; j1++) {
							for (integer j2 = 0; j2 <= 3; j2++) {
								Xmatr[j1][j2] = 0.0;
							}
							bmatr[j1] = 0.0;
							koefmatr[j1] = 0.0;
						}




						for (integer j = 0; j < 8; ++j) {

							Xmatr[0][0] += 1.0;
							Xmatr[0][1] += pointerlist[i][j].x;
							Xmatr[0][2] += pointerlist[i][j].y;
							Xmatr[0][3] += pointerlist[i][j].z;

							Xmatr[1][0] += pointerlist[i][j].x;
							Xmatr[1][1] += pointerlist[i][j].x * pointerlist[i][j].x;
							Xmatr[1][2] += pointerlist[i][j].x * pointerlist[i][j].y;
							Xmatr[1][3] += pointerlist[i][j].x * pointerlist[i][j].z;

							Xmatr[2][0] += pointerlist[i][j].y;
							Xmatr[2][1] += pointerlist[i][j].y * pointerlist[i][j].x;
							Xmatr[2][2] += pointerlist[i][j].y * pointerlist[i][j].y;
							Xmatr[2][3] += pointerlist[i][j].y * pointerlist[i][j].z;

							Xmatr[3][0] += pointerlist[i][j].z;
							Xmatr[3][1] += pointerlist[i][j].z * pointerlist[i][j].x;
							Xmatr[3][2] += pointerlist[i][j].z * pointerlist[i][j].y;
							Xmatr[3][3] += pointerlist[i][j].z * pointerlist[i][j].z;

							bmatr[0] += rthdsd_Gauss[i][j];
							bmatr[1] += pointerlist[i][j].x * rthdsd_Gauss[i][j];
							bmatr[2] += pointerlist[i][j].y * rthdsd_Gauss[i][j];
							bmatr[3] += pointerlist[i][j].z * rthdsd_Gauss[i][j];
						}


						for (integer j1 = 0; j1 <= 100; j1++) {
							koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
							koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
							koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
							koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
						}

						starting_speed_Vy_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));


						// VZ

						for (integer j = 0; j <= 7; ++j) {
							rthdsd_Gauss[i][j] = Vz47[nvtx47[j][ifound]];
						}

						for (integer j1 = 0; j1 <= 3; j1++) {
							for (integer j2 = 0; j2 <= 3; j2++) {
								Xmatr[j1][j2] = 0.0;
							}
							bmatr[j1] = 0.0;
							koefmatr[j1] = 0.0;
						}




						for (integer j = 0; j < 8; ++j) {

							Xmatr[0][0] += 1.0;
							Xmatr[0][1] += pointerlist[i][j].x;
							Xmatr[0][2] += pointerlist[i][j].y;
							Xmatr[0][3] += pointerlist[i][j].z;

							Xmatr[1][0] += pointerlist[i][j].x;
							Xmatr[1][1] += pointerlist[i][j].x * pointerlist[i][j].x;
							Xmatr[1][2] += pointerlist[i][j].x * pointerlist[i][j].y;
							Xmatr[1][3] += pointerlist[i][j].x * pointerlist[i][j].z;

							Xmatr[2][0] += pointerlist[i][j].y;
							Xmatr[2][1] += pointerlist[i][j].y * pointerlist[i][j].x;
							Xmatr[2][2] += pointerlist[i][j].y * pointerlist[i][j].y;
							Xmatr[2][3] += pointerlist[i][j].y * pointerlist[i][j].z;

							Xmatr[3][0] += pointerlist[i][j].z;
							Xmatr[3][1] += pointerlist[i][j].z * pointerlist[i][j].x;
							Xmatr[3][2] += pointerlist[i][j].z * pointerlist[i][j].y;
							Xmatr[3][3] += pointerlist[i][j].z * pointerlist[i][j].z;

							bmatr[0] += rthdsd_Gauss[i][j];
							bmatr[1] += pointerlist[i][j].x * rthdsd_Gauss[i][j];
							bmatr[2] += pointerlist[i][j].y * rthdsd_Gauss[i][j];
							bmatr[3] += pointerlist[i][j].z * rthdsd_Gauss[i][j];
						}


						for (integer j1 = 0; j1 <= 100; j1++) {
							koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
							koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
							koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
							koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
						}

						starting_speed_Vz_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));


						// Mut

						for (integer j = 0; j <= 7; ++j) {
							rthdsd_Gauss[i][j] = Mut47[nvtx47[j][ifound]];
						}

						for (integer j1 = 0; j1 <= 3; j1++) {
							for (integer j2 = 0; j2 <= 3; j2++) {
								Xmatr[j1][j2] = 0.0;
							}
							bmatr[j1] = 0.0;
							koefmatr[j1] = 0.0;
						}




						for (integer j = 0; j < 8; ++j) {

							Xmatr[0][0] += 1.0;
							Xmatr[0][1] += pointerlist[i][j].x;
							Xmatr[0][2] += pointerlist[i][j].y;
							Xmatr[0][3] += pointerlist[i][j].z;

							Xmatr[1][0] += pointerlist[i][j].x;
							Xmatr[1][1] += pointerlist[i][j].x * pointerlist[i][j].x;
							Xmatr[1][2] += pointerlist[i][j].x * pointerlist[i][j].y;
							Xmatr[1][3] += pointerlist[i][j].x * pointerlist[i][j].z;

							Xmatr[2][0] += pointerlist[i][j].y;
							Xmatr[2][1] += pointerlist[i][j].y * pointerlist[i][j].x;
							Xmatr[2][2] += pointerlist[i][j].y * pointerlist[i][j].y;
							Xmatr[2][3] += pointerlist[i][j].y * pointerlist[i][j].z;

							Xmatr[3][0] += pointerlist[i][j].z;
							Xmatr[3][1] += pointerlist[i][j].z * pointerlist[i][j].x;
							Xmatr[3][2] += pointerlist[i][j].z * pointerlist[i][j].y;
							Xmatr[3][3] += pointerlist[i][j].z * pointerlist[i][j].z;

							bmatr[0] += rthdsd_Gauss[i][j];
							bmatr[1] += pointerlist[i][j].x * rthdsd_Gauss[i][j];
							bmatr[2] += pointerlist[i][j].y * rthdsd_Gauss[i][j];
							bmatr[3] += pointerlist[i][j].z * rthdsd_Gauss[i][j];
						}


						for (integer j1 = 0; j1 <= 100; j1++) {
							koefmatr[0] = (bmatr[0] - Xmatr[0][1] * koefmatr[1] - Xmatr[0][2] * koefmatr[2] - Xmatr[0][3] * koefmatr[3]) / Xmatr[0][0];
							koefmatr[1] = (bmatr[1] - Xmatr[1][0] * koefmatr[0] - Xmatr[1][2] * koefmatr[2] - Xmatr[1][3] * koefmatr[3]) / Xmatr[1][1];
							koefmatr[2] = (bmatr[2] - Xmatr[2][0] * koefmatr[0] - Xmatr[2][1] * koefmatr[1] - Xmatr[2][3] * koefmatr[3]) / Xmatr[2][2];
							koefmatr[3] = (bmatr[3] - Xmatr[3][0] * koefmatr[0] - Xmatr[3][1] * koefmatr[1] - Xmatr[3][2] * koefmatr[2]) / Xmatr[3][3];
						}

						starting_Mut_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));



						for (integer j = 0; j <= 3; ++j) {
							delete[] Xmatr[j];
						}
						delete[] Xmatr;
						delete[] bmatr;
						delete[] koefmatr;
					}

					if (starting_speed_Vx_47 != starting_speed_Vx_47) {
						printf("NAN pri interpolation Vx %lld in allocation_memory_flow in constr_struct.cpp\n", i);
					}
					if (starting_speed_Vy_47 != starting_speed_Vy_47) {
						printf("NAN pri interpolation Vy %lld in allocation_memory_flow in constr_struct.cpp\n", i);
					}
					if (starting_speed_Vz_47 != starting_speed_Vz_47) {
						printf("NAN pri interpolation Vz %lld in allocation_memory_flow in constr_struct.cpp\n", i);
					}
					if (starting_Mut_47 != starting_Mut_47) {
						printf("NAN pri interpolation Mut %lld in allocation_memory_flow in constr_struct.cpp\n", i);
					}

					potent[VELOCITY_X_COMPONENT][i] = starting_speed_Vx_47;
					potent[VELOCITY_Y_COMPONENT][i] = starting_speed_Vy_47;
					potent[VELOCITY_Z_COMPONENT][i] = starting_speed_Vz_47;
					potent[MUT][i] = fmax(0.0, starting_Mut_47);
					// скорректированнное поле скорости должно удовлетворять уравнению неразрывности.
					potent[VXCOR][i] = starting_speed_Vx_47;
					potent[VYCOR][i] = starting_speed_Vy_47;
					potent[VZCOR][i] = starting_speed_Vz_47;

					if ((i!=0)&&(i % 100 == 0)) {
						// Хронология выполнения. Сколько осталось до конца.
						printf("\b\b\b%2d%%", (int)(100 * (maxelm - i) / maxelm));
					}
				
				}


				// По видимому сохранённый периметр гидродинамической подобласти немного уже чем реально начерченная
				// гидродинамическая подобласть. Образуется что то типо щели в приграничной области на стыке сохранеённой гидродинамической
				// области и реально начерченной. В этой щели нулевая скорость. Щель с нулевой скоростью приводит к неверным результатам.
				// Здесь мы смыкаем данную щель заполняя ненулевой скоростью по непрерывности с основной считанной гидродинамичсеской
				// подобластью. 16.12.2021

				bool bideptch_15_12_2021 = true;

				while (bideptch_15_12_2021) {

					bideptch_15_12_2021=false;

					for (integer i = 0; i < maxelm; ++i) {
						for (int G = 0; G < 6; ++G) {
							if ((bmask[i]) && ((neighbors_for_the_internal_node[G][0][i]>=maxelm)||
								((neighbors_for_the_internal_node[G][0][i]<maxelm)&&(!bmask[neighbors_for_the_internal_node[G][0][i]])))) {

								potent[VELOCITY_X_COMPONENT][neighbors_for_the_internal_node[G][0][i]] = potent[VELOCITY_X_COMPONENT][i];
								potent[VELOCITY_Y_COMPONENT][neighbors_for_the_internal_node[G][0][i]] = potent[VELOCITY_Y_COMPONENT][i];
								potent[VELOCITY_Z_COMPONENT][neighbors_for_the_internal_node[G][0][i]] = potent[VELOCITY_Z_COMPONENT][i];
								potent[MUT][neighbors_for_the_internal_node[G][0][i]] = potent[MUT][i];
								// скорректированнное поле скорости должно удовлетворять уравнению неразрывности.
								potent[VXCOR][neighbors_for_the_internal_node[G][0][i]] = potent[VXCOR][i];
								potent[VYCOR][neighbors_for_the_internal_node[G][0][i]] = potent[VYCOR][i];
								potent[VZCOR][neighbors_for_the_internal_node[G][0][i]] = potent[VZCOR][i];
								if (neighbors_for_the_internal_node[G][0][i] < maxelm) {
									bmask[neighbors_for_the_internal_node[G][0][i]] = true;
									bideptch_15_12_2021 = true;
								}
							}

						}
					}
				}

				delete[] bmask;

			}


				
			


			doublereal Speed_min1 = 1.0e60, Speed_max1 = -1.0e37;
			for (integer i_47 = 0; i_47 < maxelm; i_47++) {
				doublereal mag_speed = sqrt((potent[VELOCITY_X_COMPONENT][i_47])*(potent[VELOCITY_X_COMPONENT][i_47]) + (potent[VELOCITY_Y_COMPONENT][i_47])*(potent[VELOCITY_Y_COMPONENT][i_47]) + (potent[VELOCITY_Z_COMPONENT][i_47])*(potent[VELOCITY_Z_COMPONENT][i_47]));
				if (mag_speed < Speed_min1) {
					Speed_min1 = mag_speed;
				}
				if (mag_speed > Speed_max1) {
					Speed_max1 = mag_speed;
				}
			}
			printf("Interpolation: Minimum speed = %e, Maximum speed = %e \n", Speed_min1, Speed_max1);
			//system("PAUSE");

			for (integer i_47 = 0; i_47 < maxelm; i_47++) {
				// Сохраняем модуль скорости неизменным. 29.03.2019
				potent[VELOCITY_X_COMPONENT][i_47] *= Speed_max / Speed_max1;
				potent[VELOCITY_Y_COMPONENT][i_47] *= Speed_max / Speed_max1;
				potent[VELOCITY_Z_COMPONENT][i_47] *= Speed_max / Speed_max1;
			}

			for (integer i_1 = 0; i_1 < 8; ++i_1) {
				delete[] i_22[i_1];
			}
			delete[] i_22;

			if (oct_load1 != nullptr) {
				for (integer i = 0; i < 8; ++i) {
					for (integer j = 0; j < 8; ++j) {
						delete[] oct_load1[i][j];
						oct_load1[i][j] = nullptr;
					}
				}

				for (integer i = 0; i < 8; ++i) {
					delete[] oct_load1[i];
					oct_load1[i] = nullptr;
				}
				delete[] oct_load1;
				oct_load1 = nullptr;
			}

			for (integer i = 0; i < maxelm; ++i) {
				delete[] pointerlist[i];
				delete[] rthdsd_Gauss[i];
			}
			delete[] pointerlist;
			delete[] rthdsd_Gauss;

			// Освобождение оперативной памяти.
			delete[] x47;
			delete[] y47;
			delete[] z47;
			delete[] Vx47;
			delete[] Vy47;
			delete[] Vz47;
			delete[] Mut47;
			for (integer i_47 = 0; i_47 < 8; i_47++) {
				delete[] nvtx47[i_47];
			}
			delete[] nvtx47;

			printf("done.\n");
			//system("PAUSE");
		}
	    }
	}
	// Лучше начинать с нулевого поля скорости.
	// на границе должны быть выполнены граничные условия 8 мая 2013г (revised 2 апреля 2019г).
	// Скорость тоже должна быть инициализирована с учётом граничных условий Дирихле.
#pragma omp parallel for
	for (integer i = 0; i < maxbound; ++i) {		
		
		if ((border_neighbor[i].MCB < (ls + lw)) && (border_neighbor[i].MCB >= ls)) {

			// Добавил УСЛОВИЕ СИММЕТРИИ 02.04.2019.

			if ((w[border_neighbor[i].MCB - ls].bpressure) || 
				(w[border_neighbor[i].MCB - ls].bopening)) {
				// Выходная граница потока на которой стоит условие Неймана.
				// Присваиваем скорость из ближайшей внутренней точки расчётной области.
				potent[VELOCITY_X_COMPONENT][border_neighbor[i].iB] = potent[VELOCITY_X_COMPONENT][border_neighbor[i].iI];
				potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iB] = potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iI];
				potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iB] = potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iI];
				potent[MUT][border_neighbor[i].iB] = potent[MUT][border_neighbor[i].iI];
				// скорректированнное поле скорости должно удовлетворять уравнению неразрывности.
				potent[VXCOR][border_neighbor[i].iB] = potent[VELOCITY_X_COMPONENT][border_neighbor[i].iI];
				potent[VYCOR][border_neighbor[i].iB] = potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iI];
				potent[VZCOR][border_neighbor[i].iB] = potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iI];
				//printf("Vx=%e, Vy=%e, Vz=%e\n",w[border_neighbor[i].MCB-ls].Vx, w[border_neighbor[i].MCB-ls].Vy, w[border_neighbor[i].MCB-ls].Vz);
				//system("PAUSE"); // debug
			}
			else if (w[border_neighbor[i].MCB - ls].bsymmetry) {
				// Тангенсальные компоненты скорости условие Неймана,
				// нормальная компонента скорости равна нулю.

				potent[MUT][border_neighbor[i].iB] = potent[MUT][border_neighbor[i].iI];

				switch (w[border_neighbor[i].MCB - ls].iPlane) {
				case XY_PLANE:
					potent[VELOCITY_X_COMPONENT][border_neighbor[i].iB] = potent[VELOCITY_X_COMPONENT][border_neighbor[i].iI];
					potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iB] = potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iI];
					potent[VXCOR][border_neighbor[i].iB] = potent[VELOCITY_X_COMPONENT][border_neighbor[i].iI];
					potent[VYCOR][border_neighbor[i].iB] = potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iI];
					potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iB] = 0.0;
					potent[VZCOR][border_neighbor[i].iB] = 0.0;
					break;
				case XZ_PLANE:
					potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iB] = 0.0;
					potent[VYCOR][border_neighbor[i].iB] = 0.0;
					potent[VELOCITY_X_COMPONENT][border_neighbor[i].iB] = potent[VELOCITY_X_COMPONENT][border_neighbor[i].iI];
					potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iB] = potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iI];
					potent[VXCOR][border_neighbor[i].iB] = potent[VELOCITY_X_COMPONENT][border_neighbor[i].iI];
					potent[VZCOR][border_neighbor[i].iB] = potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iI];
					break;
				case YZ_PLANE: 
					potent[VELOCITY_X_COMPONENT][border_neighbor[i].iB] = 0.0;
					potent[VXCOR][border_neighbor[i].iB] = 0.0;
					potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iB] = potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iI];
					potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iB] = potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iI];
					potent[VYCOR][border_neighbor[i].iB] = potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iI];
					potent[VZCOR][border_neighbor[i].iB] = potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iI];
					break;
				default:
					printf("ERROR: Undefined iPlane wall...\n");
					system("PAUSE");
					exit(1);
					break;
				}
			}
			else {
				// на границе задана скорость
				potent[VELOCITY_X_COMPONENT][border_neighbor[i].iB] = w[border_neighbor[i].MCB - ls].Vx;
				potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iB] = w[border_neighbor[i].MCB - ls].Vy;
				potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iB] = w[border_neighbor[i].MCB - ls].Vz;
				potent[MUT][border_neighbor[i].iB] = 2.0*prop_b[MU_DYNAMIC_VISCOSITY][i] / prop_b[RHO][i];
				// скорректированное поле скорости должно удовлетворять уравнению неразрывности.
				potent[VXCOR][border_neighbor[i].iB] = w[border_neighbor[i].MCB - ls].Vx;
				potent[VYCOR][border_neighbor[i].iB] = w[border_neighbor[i].MCB - ls].Vy;
				potent[VZCOR][border_neighbor[i].iB] = w[border_neighbor[i].MCB - ls].Vz;
				//printf("Vx=%e, Vy=%e, Vz=%e\n",w[border_neighbor[i].MCB-ls].Vx, w[border_neighbor[i].MCB-ls].Vy, w[border_neighbor[i].MCB-ls].Vz);
				//system("PAUSE"); // debug
			}
		}
		else {
			// Абсолютно твёрдая неподвижная стенка на которой выполнено условие прилипания,
			// скорость равна нулю.
			// 02.04.2019
			potent[VELOCITY_X_COMPONENT][border_neighbor[i].iB] = 0.0;
			potent[VELOCITY_Y_COMPONENT][border_neighbor[i].iB] = 0.0;
			potent[VELOCITY_Z_COMPONENT][border_neighbor[i].iB] = 0.0;
			potent[MUT][border_neighbor[i].iB] = 0.0;
			// скорректированное поле скорости должно удовлетворять уравнению неразрывности.
			potent[VXCOR][border_neighbor[i].iB] = 0.0;
			potent[VYCOR][border_neighbor[i].iB] = 0.0;
			potent[VZCOR][border_neighbor[i].iB] = 0.0;
		}
		
	}
	//*/
	 
	alpha=nullptr;
	alpha = new doublereal[12];
	if (alpha==nullptr) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for alpha constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	
	// Патанкар Fluent Гаврилов Андрей
    alpha[VELOCITY_X_COMPONENT]=0.7; // 0.5 0.7 0.8
    alpha[VELOCITY_Y_COMPONENT]=0.7;
    alpha[VELOCITY_Z_COMPONENT]=0.7;
    alpha[PRESS]=0.3; //0.8 0.3 (0.2, но вообще говоря он использует SIMPLEC по умолчанию).
	alpha[NUSHA_SL] = 0.8;// 1.0;// 0.7; // Модифицированная кинематическая турбулентная вязкость.
	alpha[TURBULENT_KINETIK_ENERGY_SL] = 0.8;// 1.0;// 0.7; // Кинетическая энергия турбулентных пульсаций.
	alpha[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL] = 0.8;// 1.0;// 0.7; // удельная скорость её диссипации
	alpha[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL] = 0.8;
	alpha[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL] = 0.8;
	alpha[GAMMA_LANGTRY_MENTER_SL] = 0.8;
	alpha[RE_THETA_LANGTRY_MENTER_SL] = 0.8;

   // коэффициенты матрицы СЛАУ для внутренних КО.
   sl=nullptr;
   if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
	   sl = new equation3D*[12];
	   if (sl == nullptr) {
		   // недостаточно памяти на данном оборудовании.
		   printf("Problem: not enough memory on your equipment for slau flow constr struct...\n");
		   printf("Please any key to exit...\n");
		   //system("PAUSE");
		   system("pause");
		   exit(1);
	   }


	   for (integer i = 0; i < 12; ++i) {
		   sl[i] = nullptr;
	   }

#pragma omp parallel for
	   for (integer i = 0; i < 12; ++i) {
		   switch (i) {
		   case VELOCITY_X_COMPONENT: sl[VELOCITY_X_COMPONENT] = new equation3D[maxelm];
			   if (sl[VELOCITY_X_COMPONENT] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[VX] constr struct...\n");
				   printf("Please any key to exit...\n");
				   // system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case VELOCITY_Y_COMPONENT: sl[VELOCITY_Y_COMPONENT] = new equation3D[maxelm];
			   if (sl[VELOCITY_Y_COMPONENT] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[VY] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case VELOCITY_Z_COMPONENT: sl[VELOCITY_Z_COMPONENT] = new equation3D[maxelm];
			   if (sl[VELOCITY_Z_COMPONENT] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[VZ] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case PRESS: sl[PRESS] = nullptr; break;
		   case PAM: sl[PAM] = new equation3D[maxelm];
			   if (sl[PAM] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[PAM] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case NUSHA_SL:
			   // Модифицированная кинематическая турбулентная вязкость.
			   sl[NUSHA_SL] = new equation3D[maxelm];
			   if (sl[NUSHA_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[NUSHA_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case TURBULENT_KINETIK_ENERGY_SL: // k SST
			   // Кинетическая энергия турбулентных пульсаций
			   sl[TURBULENT_KINETIK_ENERGY_SL] = new equation3D[maxelm];
			   if (sl[TURBULENT_KINETIK_ENERGY_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[TURBULENT_KINETIK_ENERGY_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL: // omega SST
			   // Удельная скорость её диссипации.
			   sl[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL] = new equation3D[maxelm];
			   if (sl[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL: // Двухслойная модель на основе Standart KE model
			   // turbulence kinetik energy
			   sl[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL] = new equation3D[maxelm];
			   if (sl[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL: // Двухслойная модель на основе Standart KE model
			   // turbulent dissipation rate (epsilon)
			   sl[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL] = new equation3D[maxelm];
			   if (sl[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case GAMMA_LANGTRY_MENTER_SL:  // Модель ламинарно турбулентного перехода.
			   sl[GAMMA_LANGTRY_MENTER_SL] = new equation3D[maxelm];
			   if (sl[GAMMA_LANGTRY_MENTER_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[GAMMA_LANGTRY_MENTER_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case RE_THETA_LANGTRY_MENTER_SL: // Модель ламинарно турбулентного перехода.
			   sl[RE_THETA_LANGTRY_MENTER_SL] = new equation3D[maxelm];
			   if (sl[RE_THETA_LANGTRY_MENTER_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau[RE_THETA_LANGTRY_MENTER_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   }
	   }
   }

   // коэффициенты матрицы СЛАУ для граничных КО
   slb=nullptr;
   if (steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::NETWORK_T_UNSTEADY) {
	   slb = new equation3D_bon*[12];
	   if (slb == nullptr) {
		   // недостаточно памяти на данном оборудовании.
		   printf("Problem: not enough memory on your equipment for slau boundary constr struct...\n");
		   printf("Please any key to exit...\n");
		   //system("PAUSE");
		   system("pause");
		   exit(1);
	   }

	   for (integer i = 0; i < 12; ++i) {
		   slb[i] = nullptr;
	   }

#pragma omp parallel for
	   for (int i = 0; i < 12; ++i) {
		   switch (i) {
		   case VELOCITY_X_COMPONENT: slb[VELOCITY_X_COMPONENT] = new equation3D_bon[maxbound];
			   if (slb[VELOCITY_X_COMPONENT] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[VX] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case VELOCITY_Y_COMPONENT: slb[VELOCITY_Y_COMPONENT] = new equation3D_bon[maxbound];
			   if (slb[VELOCITY_Y_COMPONENT] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[VY] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case VELOCITY_Z_COMPONENT: slb[VELOCITY_Z_COMPONENT] = new equation3D_bon[maxbound];
			   if (slb[VELOCITY_Z_COMPONENT] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[VZ] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case PRESS: slb[PRESS] = nullptr; break;
		   case PAM: slb[PAM] = new equation3D_bon[maxbound];
			   if (slb[PAM] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[PAM] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case NUSHA_SL: slb[NUSHA_SL] = new equation3D_bon[maxbound];
			   if (slb[NUSHA_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[NUSHA_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case TURBULENT_KINETIK_ENERGY_SL: // k SST
			   slb[TURBULENT_KINETIK_ENERGY_SL] = new equation3D_bon[maxbound];
			   if (slb[TURBULENT_KINETIK_ENERGY_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[TURBULENT_KINETIK_ENERGY_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL: // omega SST
			   slb[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL] = new equation3D_bon[maxbound];
			   if (slb[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[TURBULENT_SPECIFIC_DISSIPATION_RATE_OMEGA_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL: // Двухслойная модель на основе Standart KE model
				// turbulence kinetik energy
			   slb[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL] = new equation3D_bon[maxbound];
			   if (slb[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[TURBULENT_KINETIK_ENERGY_STD_K_EPS_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL: // Двухслойная модель на основе Standart KE model
				// turbulent dissipation rate (epsilon)
			   slb[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL] = new equation3D_bon[maxbound];
			   if (slb[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[TURBULENT_DISSIPATION_RATE_EPSILON_STD_K_EPS_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case GAMMA_LANGTRY_MENTER_SL : // Модель Ламинарно турбулентного перехода.
			   slb[GAMMA_LANGTRY_MENTER_SL] = new equation3D_bon[maxbound];
			   if (slb[GAMMA_LANGTRY_MENTER_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[GAMMA_LANGTRY_MENTER_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   case RE_THETA_LANGTRY_MENTER_SL :  // Модель Ламинарно турбулентного перехода.
			   slb[RE_THETA_LANGTRY_MENTER_SL] = new equation3D_bon[maxbound];
			   if (slb[RE_THETA_LANGTRY_MENTER_SL] == nullptr) {
				   // недостаточно памяти на данном оборудовании.
				   printf("Problem: not enough memory on your equipment for slau_bon[RE_THETA_LANGTRY_MENTER_SL] constr struct...\n");
				   printf("Please any key to exit...\n");
				   //system("PAUSE");
				   system("pause");
				   exit(1);
			   }
			   break;
		   }
	   }
   }
				
  
  // omp_set_num_threads(1);

} // allocation_memory_flow


  // выделение оперативной памяти для задачи гидродинамики.
void allocation_memory_flow_2(
	doublereal** &diag_coef,  integer maxelm, integer maxbound,
	 doublereal* &SInvariantStrainRateTensor,
	doublereal** &mf) {

	// диагональные коэффициенты компонент скорости 
	// необходимо хранить для вычисления поправки Рхи-Чоу.
	diag_coef = nullptr;
	diag_coef = new doublereal*[3];
	if (diag_coef == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for diag_coef constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}

	for (integer i = 0; i<3; ++i) {
		diag_coef[i] = nullptr;
		diag_coef[i] = new doublereal[maxelm + maxbound];
		if (diag_coef[i] == nullptr) {
			// недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for diag_coef["<< i <<"] constr struct..." << std::endl;
			
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// инициализация
	for (integer i = 0; i<3; ++i) for (integer j = 0; j<(maxelm + maxbound); ++j) diag_coef[i][j] = 1.0;

	// Выделение памяти под S инвариант тензора скоростей-деформаций.
	SInvariantStrainRateTensor = nullptr;
	SInvariantStrainRateTensor = new doublereal[maxelm + maxbound];
	if (SInvariantStrainRateTensor == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for SInvariantStrainRateTensor constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	// инициализация S инварианта strain-rate тензора:
	for (integer i = 0; i<(maxelm + maxbound); ++i) {
		SInvariantStrainRateTensor[i] = 0.0; // инициализация
	}

	mf = nullptr;
	mf = new doublereal*[maxelm];
	if (mf == nullptr) {
		// недостаточно памяти на данном оборудовании.
		std::cout << "Problem: not enough memory on your equipment for mf constr struct..."<< std::endl;
		std::cout << "Please any key to exit..."<< std::endl;
		//system("PAUSE");
		system("pause");
		exit(1);
	}

	for (integer i = 0; i<maxelm; ++i) {
		mf[i] = nullptr;
		mf[i] = new doublereal[6];
		if (mf[i] == nullptr) {
			// недостаточно памяти на данном оборудовании.
			std::cout << "Problem: not enough memory on your equipment for mf["<< i <<"] constr struct..." << std::endl;
			
			std::cout << "Please any key to exit..."<< std::endl;
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// инициализация:
	for (integer i = 0; i<maxelm; ++i) {
		for (integer j = 0; j<6; ++j) {
			mf[i][j] = 0.0;
		}
	}
}

// создаёт связи между контрольными объёмами для графической 
// визуализации.
// универсальна: подходит и для температуры и для течения
// 29.12.2021 - переписал функцию в параллельном режиме с использованием OpenMP. Теперь данная функция работает параллельно.
void constr_nvtxcell(int* &evt, BOUND* &border_neighbor, int maxbound, 
	int maxelm, bool bextendedprint, int** &nvtxcell, integer &ncell_gl,
	int inx, int iny, int inz, TOCKA_SHORT_INT* &tck_int_list) {

	// Если bextendedprinteger = true то мы имеем дело с расширенной печатью включая граничные значения узлов.

	// для графической визуализации
	// nvtxcell - связи для контрольных объёмов.
	// ncell - количество связей для контрольных объёмов.
	ncell_gl = 0;
	// ncell > maxelm.
	//integer i1, i2, i3, i4, i5, i6, i7, i8;
	//for (i = 0; i<(inx - 1); ++i) for (j = 0; j<(iny - 1); ++j) for (k = 0; k<(inz - 1); ++k) {
	// Сокращает число просмотров на больших моделях.

	int incore = number_cores();

#ifdef _OPENMP
	omp_set_num_threads(incore); // установка числа потоков
#else
	incore = 1; // один поток.
#endif

	integer ncell = 0;
	// integer i,j,k;
	// integer iP_k, iP_k_p1, iP_j, iP_j_p1;

/*#pragma omp parallel for reduction(+:ncell)
	for (integer iscan = 0; iscan < maxelm; ++iscan) {

		int i = tck_int_list[iscan].i;
		int j = tck_int_list[iscan].j;
		int k = tck_int_list[iscan].k;

		if ((i < (inx - 1)) && (j < (iny - 1)) && (k < (inz - 1)))
		{

			// Подсчёт количества связей между внутренними КО в цикле
			 integer iP_k=k * inx*iny;
			 integer iP_k_p1=(k+1) * inx*iny;
			 integer iP_j=j * inx;
			 integer iP_j_p1=(j+1) * inx;

			integer i1 = i + iP_j + iP_k;
			integer i2 = (i + 1) + iP_j + iP_k;
			integer i3 = (i + 1) + iP_j_p1 + iP_k;
			integer i4 = i + iP_j_p1 + iP_k;
			integer i5 = i + iP_j + iP_k_p1;
			integer i6 = (i + 1) + iP_j + iP_k_p1;
			integer i7 = (i + 1) + iP_j_p1 + iP_k_p1;
			integer i8 = i + iP_j_p1 + iP_k_p1;
		
		if ((evt[i1] > 0) && (evt[i2] > 0) && (evt[i3] > 0) && (evt[i4] > 0) && (evt[i5] > 0) && (evt[i6] > 0) && (evt[i7] > 0) && (evt[i8] > 0)) {
			ncell++;

			if (bextendedprint) {

				// сканируем каждую внутреннюю нормаль в поисках граничных 
				// плоских бесконечно-тонких объёмов на которых выставлено граничное условие.
				integer inorm = E_SIDE;
				bool bfound1 = false, bfound2 = false, bfound3 = false, bfound4 = false;

				// поиск.
				for (int isel = 0; isel < maxbound; ++isel) {
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i1] - 1)) {
						bfound1 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i4] - 1)) {
						bfound2 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i5] - 1)) {
						bfound3 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i8] - 1)) {
						bfound4 = true;
					}
				}

				if (bfound1&&bfound2&&bfound3&&bfound4) {
					ncell++;
				}
				bfound1 = false;
				bfound2 = false;
				bfound3 = false;
				bfound4 = false;

				inorm = W_SIDE;
				// поиск.
				for (int isel = 0; isel < maxbound; ++isel) {
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i2] - 1)) {
						bfound1 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i3] - 1)) {
						bfound2 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i6] - 1)) {
						bfound3 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i7] - 1)) {
						bfound4 = true;
					}
				}

				if (bfound1&&bfound2&&bfound3&&bfound4) {
					ncell++;
				}

				bfound1 = false;
				bfound2 = false;
				bfound3 = false;
				bfound4 = false;

				inorm = N_SIDE;
				// поиск.
				for (int isel = 0; isel < maxbound; ++isel) {
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i1] - 1)) {
						bfound1 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i2] - 1)) {
						bfound2 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i5] - 1)) {
						bfound3 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i6] - 1)) {
						bfound4 = true;
					}
				}

				if (bfound1&&bfound2&&bfound3&&bfound4) {
					ncell++;
				}

				bfound1 = false;
				bfound2 = false;
				bfound3 = false;
				bfound4 = false;
				inorm = S_SIDE;
				// поиск.
				for (int isel = 0; isel < maxbound; ++isel) {
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i3] - 1)) {
						bfound1 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i4] - 1)) {
						bfound2 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i7] - 1)) {
						bfound3 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i8] - 1)) {
						bfound4 = true;
					}
				}

				if (bfound1&&bfound2&&bfound3&&bfound4) {
					ncell++;
				}

				bfound1 = false;
				bfound2 = false;
				bfound3 = false;
				bfound4 = false;
				inorm = T_SIDE;
				// поиск.
				for (int isel = 0; isel < maxbound; ++isel) {
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i1] - 1)) {
						bfound1 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i2] - 1)) {
						bfound2 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i3] - 1)) {
						bfound3 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i4] - 1)) {
						bfound4 = true;
					}
				}

				if (bfound1&&bfound2&&bfound3&&bfound4) {
					ncell++;
				}

				bfound1 = false;
				bfound2 = false;
				bfound3 = false;
				bfound4 = false;
				inorm = B_SIDE;
				// поиск.
				for (int isel = 0; isel < maxbound; ++isel) {
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i5] - 1)) {
						bfound1 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i6] - 1)) {
						bfound2 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i7] - 1)) {
						bfound3 = true;
					}
					if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i8] - 1)) {
						bfound4 = true;
					}
				}

				if (bfound1&&bfound2&&bfound3&&bfound4) {
					ncell++;
				}

			} // bextendedprint

		}
	  }
	}
	*/

	integer*** idx = new integer**[incore];
	for (integer tid = 0; tid < incore; ++tid) {
		idx[tid] = new integer*[7];
		for (integer i_63 = 0; i_63 < 7; ++i_63) {
			idx[tid][i_63] = new integer[maxelm];
			for (integer iscan = 0; iscan < maxelm; ++iscan) {
				idx[tid][i_63][iscan] = -1; // инициализация.
			}
		}
	}
	integer* tid_count = new integer[incore];
	integer* tid_countSum = new integer[incore];
	for (integer tid = 0; tid < incore; ++tid) {
		tid_count[tid] = 0;
		tid_countSum[tid] = 0;
	}

	ncell = 0;
#pragma omp parallel for reduction(+:ncell)
	for (integer iscan = 0; iscan < maxelm; ++iscan) {


		int tid = omp_get_thread_num();

		int i = tck_int_list[iscan].i;
		int j = tck_int_list[iscan].j;
		int k = tck_int_list[iscan].k;

		if ((i < (inx - 1)) && (j < (iny - 1)) && (k < (inz - 1)))
		{

			// Подсчёт количества связей между внутренними КО в цикле
			integer iP_k = k * inx*iny;
			integer iP_k_p1 = (k + 1) * inx*iny;
			integer iP_j = j * inx;
			integer iP_j_p1 = (j + 1) * inx;

			integer i1 = i + iP_j + iP_k;
			integer i2 = (i + 1) + iP_j + iP_k;
			integer i3 = (i + 1) + iP_j_p1 + iP_k;
			integer i4 = i + iP_j_p1 + iP_k;
			integer i5 = i + iP_j + iP_k_p1;
			integer i6 = (i + 1) + iP_j + iP_k_p1;
			integer i7 = (i + 1) + iP_j_p1 + iP_k_p1;
			integer i8 = i + iP_j_p1 + iP_k_p1;

			if ((evt[i1] > 0) && (evt[i2] > 0) && (evt[i3] > 0) && (evt[i4] > 0) && (evt[i5] > 0) && (evt[i6] > 0) && (evt[i7] > 0) && (evt[i8] > 0)) {

				idx[tid][6][iscan] = tid_count[tid];
				++tid_count[tid];
				ncell++;

				if (bextendedprint) {

					// сканируем каждую внутреннюю нормаль в поисках граничных 
					// плоских бесконечно-тонких объёмов на которых выставлено граничное условие.
					integer inorm = E_SIDE;
					bool bfound1 = false, bfound2 = false, bfound3 = false, bfound4 = false;

					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i1] - 1)) {
							bfound1 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i4] - 1)) {
							bfound2 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i5] - 1)) {
							bfound3 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i8] - 1)) {
							bfound4 = true;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						idx[tid][E_SIDE][iscan] = tid_count[tid];
						++tid_count[tid];
						ncell++;
					}
					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;

					inorm = W_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i2] - 1)) {
							bfound1 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i3] - 1)) {
							bfound2 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i6] - 1)) {
							bfound3 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i7] - 1)) {
							bfound4 = true;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						idx[tid][W_SIDE][iscan] = tid_count[tid];
						++tid_count[tid];
						ncell++;
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;

					inorm = N_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i1] - 1)) {
							bfound1 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i2] - 1)) {
							bfound2 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i5] - 1)) {
							bfound3 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i6] - 1)) {
							bfound4 = true;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						idx[tid][N_SIDE][iscan] = tid_count[tid];
						++tid_count[tid];
						ncell++;
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;
					inorm = S_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i3] - 1)) {
							bfound1 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i4] - 1)) {
							bfound2 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i7] - 1)) {
							bfound3 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i8] - 1)) {
							bfound4 = true;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						idx[tid][S_SIDE][iscan] = tid_count[tid];
						++tid_count[tid];
						ncell++;
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;
					inorm = T_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i1] - 1)) {
							bfound1 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i2] - 1)) {
							bfound2 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i3] - 1)) {
							bfound3 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i4] - 1)) {
							bfound4 = true;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						idx[tid][T_SIDE][iscan] = tid_count[tid];
						++tid_count[tid];
						ncell++;
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;
					inorm = B_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i5] - 1)) {
							bfound1 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i6] - 1)) {
							bfound2 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i7] - 1)) {
							bfound3 = true;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i8] - 1)) {
							bfound4 = true;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						idx[tid][B_SIDE][iscan] = tid_count[tid];
						++tid_count[tid];
						ncell++;
					}

				} // bextendedprint

			}
		}
	}

	for (integer tid = 1; tid < incore; ++tid) {
		tid_countSum[tid] = tid_countSum[tid-1]+ tid_count[tid-1];
	}

	// Выделение ОП 
	ncell_gl = ncell;
	nvtxcell = nullptr;
	nvtxcell = new int*[NUMBER_OF_VERTEX_FINITE_ELEMENT()];
	if (nvtxcell == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for nvtxcell constr struct...\n");
		printf("Please any key to exit...\n");
		// system("PAUSE");
		system("pause");
		exit(1);
	}
	//integer l = 0;
	for (integer i = 0; i<NUMBER_OF_VERTEX_FINITE_ELEMENT(); ++i) {
		nvtxcell[i] = nullptr;
		nvtxcell[i] = new int[ncell];
		if (nvtxcell[i] == nullptr) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem: not enough memory on your equipment for nvtxcell[%lld] constr struct...\n", i);
#else
			printf("Problem: not enough memory on your equipment for nvtxcell[%d] constr struct...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	//for (i = 0; i<(inx - 1); ++i) for (j = 0; j<(iny - 1); ++j) for (k = 0; k<(inz - 1); ++k) {
	// Сокращает число просмотров на больших моделях.

#pragma omp parallel for 
	for (integer iscan = 0; iscan < maxelm; ++iscan) {

		int tid = omp_get_thread_num();

		 int i = tck_int_list[iscan].i;
		 int j = tck_int_list[iscan].j;
		 int k = tck_int_list[iscan].k;

		if ((i < (inx - 1)) && (j < (iny - 1)) && (k < (inz - 1)))
		{

			integer iP_k=k * inx*iny;
			integer iP_k_p1=(k+1) * inx*iny;
			integer iP_j=j * inx;
			integer iP_j_p1=(j+1) * inx;

			// Установка восьмиточечных связей для связывания КО.
			integer i1 = i + iP_j + iP_k;
			integer i2 = (i + 1) + iP_j + iP_k;
			integer i3 = (i + 1) + iP_j_p1 + iP_k;
			integer i4 = i + iP_j_p1 + iP_k;
			integer i5 = i + iP_j + iP_k_p1;
			integer i6 = (i + 1) + iP_j + iP_k_p1;
			integer i7 = (i + 1) + iP_j_p1 + iP_k_p1;
			integer i8 = i + iP_j_p1 + iP_k_p1;

			
			if ((evt[i1] > 0) && (evt[i2] > 0) && (evt[i3] > 0) && (evt[i4] > 0) && (evt[i5] > 0) && (evt[i6] > 0) && (evt[i7] > 0) && (evt[i8] > 0)) {

				integer l = idx[tid][6][iscan] + tid_countSum[tid];

				nvtxcell[0][l] = evt[i1]; nvtxcell[1][l] = evt[i2]; nvtxcell[2][l] = evt[i3]; nvtxcell[3][l] = evt[i4];
				nvtxcell[4][l] = evt[i5]; nvtxcell[5][l] = evt[i6]; nvtxcell[6][l] = evt[i7]; nvtxcell[7][l] = evt[i8];
				

				if (bextendedprint) {

					// сканируем каждую внутреннюю нормаль в поисках граничных 
					// плоских бесконечно-тонких объёмов на которых вытавлено граничное условие.
					integer inorm = E_SIDE;
					bool bfound1 = false, bfound2 = false, bfound3 = false, bfound4 = false;
					int im1 = -1, im2 = -1, im3 = -1, im4 = -1;


					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i1] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
							//im1=border_neighbor[isel].iB+1;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i4] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
							//im2=border_neighbor[isel].iB+1;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i5] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
							//im3=border_neighbor[isel].iB+1;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i8] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
							//im4=border_neighbor[isel].iB+1;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {

						l = idx[tid][E_SIDE][iscan] + tid_countSum[tid];

						nvtxcell[0][l] = im1; nvtxcell[1][l] = evt[i2]; nvtxcell[2][l] = evt[i3]; nvtxcell[3][l] = im2;
						nvtxcell[4][l] = im3; nvtxcell[5][l] = evt[i6]; nvtxcell[6][l] = evt[i7]; nvtxcell[7][l] = im4;
						

					}
					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;

					inorm = W_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i2] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
							//im1=border_neighbor[isel].iB+1;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i3] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
							//im2=border_neighbor[isel].iB+1;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i6] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i7] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {

						l = idx[tid][W_SIDE][iscan] + tid_countSum[tid];

						nvtxcell[0][l] = evt[i1]; nvtxcell[1][l] = im1; nvtxcell[2][l] = im2; nvtxcell[3][l] = evt[i4];
						nvtxcell[4][l] = evt[i5]; nvtxcell[5][l] = im3; nvtxcell[6][l] = im4; nvtxcell[7][l] = evt[i8];
						
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;

					inorm = N_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i1] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i2] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i5] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i6] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {

						l = idx[tid][N_SIDE][iscan] + tid_countSum[tid];

						nvtxcell[0][l] = im1; nvtxcell[1][l] = im2; nvtxcell[2][l] = evt[i3]; nvtxcell[3][l] = evt[i4];
						nvtxcell[4][l] = im3; nvtxcell[5][l] = im4; nvtxcell[6][l] = evt[i7]; nvtxcell[7][l] = evt[i8];
						
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;
					inorm = S_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i3] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i4] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i7] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i8] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {

						l = idx[tid][S_SIDE][iscan] + tid_countSum[tid];

						nvtxcell[0][l] = evt[i1]; nvtxcell[1][l] = evt[i2]; nvtxcell[2][l] = im1; nvtxcell[3][l] = im2;
						nvtxcell[4][l] = evt[i5]; nvtxcell[5][l] = evt[i6]; nvtxcell[6][l] = im3; nvtxcell[7][l] = im4;
						
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;
					inorm = T_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i1] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i2] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i3] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i4] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {

						l = idx[tid][T_SIDE][iscan] + tid_countSum[tid];

						nvtxcell[0][l] = im1; nvtxcell[1][l] = im2; nvtxcell[2][l] = im3; nvtxcell[3][l] = im4;
						nvtxcell[4][l] = evt[i5]; nvtxcell[5][l] = evt[i6]; nvtxcell[6][l] = evt[i7]; nvtxcell[7][l] = evt[i8];
						
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;
					inorm = B_SIDE;
					// поиск.
					for (int isel = 0; isel < maxbound; ++isel) {
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i5] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i6] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i7] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((border_neighbor[isel].Norm == inorm) && (border_neighbor[isel].iI == evt[i8] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {

						l = idx[tid][B_SIDE][iscan] + tid_countSum[tid];

						nvtxcell[0][l] = evt[i1]; nvtxcell[1][l] = evt[i2]; nvtxcell[2][l] = evt[i3]; nvtxcell[3][l] = evt[i4];
						nvtxcell[4][l] = im1; nvtxcell[5][l] = im2; nvtxcell[6][l] = im3; nvtxcell[7][l] = im4;
						
					}

				} // bextendedprint

			}
		}
	}


	
	for (integer tid = 0; tid < incore; ++tid) {
		for (integer i_63 = 0; i_63 < 7; ++i_63) {
			delete[] idx[tid][i_63];			
		}
	}
	for (integer tid = 0; tid < incore; ++tid) {
		delete[] idx[tid];
	}
	delete[] idx;
	delete[] tid_count;
	delete[] tid_countSum;
	

	//omp_set_num_threads(1);
} // constr_nvtxcell

// Вычисление расстояния до ближайшей стенки
void calcdistwall(FLOW &f, integer ls, integer lw, WALL* &w) {
	// Вычисление расстояния до ближайшей стенки требуется во многих моделях турбулентности.
    f.rdistWall=nullptr;
	f.rdistWall=new doublereal[f.maxelm+f.maxbound];
	if (f.rdistWall==nullptr) {
	      // недостаточно памяти на данном оборудовании.
		  printf("Problem: not enough memory on your equipment for rdistWall constr struct...\n");
		  printf("Please any key to exit...\n");
		  //system("PAUSE");
		  system("pause");
		  exit(1);
	}
#pragma omp parallel for
	for (integer i = 0; i < (f.maxelm + f.maxbound); ++i) {
		f.rdistWall[i] = 1e30; // присваиваем очень большое вещественное число.
	}

	bool *bwall=nullptr;
	bwall=new bool[f.maxbound]; 
	if (bwall==nullptr) {
	      // недостаточно памяти на данном оборудовании.
		  printf("Problem: not enough memory on your equipment for bwall constr struct...\n");
		  printf("Please any key to exit...\n");
		  //system("PAUSE");
		  system("pause");
		  exit(1);
	}
#pragma omp parallel for
	for (integer i = 0; i < f.maxbound; ++i) {
		bwall[i] = false; // граница не принадлежит твёрдой стенке.
	}
	doublereal **rposition=nullptr;
	rposition=new doublereal*[3];
	if (rposition==nullptr) {
	      // недостаточно памяти на данном оборудовании.
		   std::cout << "Problem: not enough memory on your equipment for rposition constr struct..." << std::endl;
		   std::cout << "Please any key to exit..." << std::endl;
		  //system("PAUSE");
		  system("pause");
		  exit(1);
	}
	for (integer i=0; i<3; ++i) {
		rposition[i]=nullptr;
		rposition[i]=new doublereal[f.maxbound];
		if (rposition[i]==nullptr) {
	      // недостаточно памяти на данном оборудовании.
          std::cout << "Problem: not enough memory on your equipment for rposition[" << i <<"] constr struct..." << std::endl;
		  
		  std::cout << "Please any key to exit..." << std::endl;
		  //system("PAUSE");
		  system("pause");
		  exit(1);
	    }
	}
	const integer X=0;
	const integer Y=1;
	const integer Z=2;
	for (integer i=0; i<3; ++i) {
#pragma omp parallel for
		for (integer j=0; j<(f.maxbound); ++j) {
			rposition[i][j]=0.0; // инициализация.
		}
	}

	for (integer i=0; i<f.maxelm; ++i) {
		// цикл по всем внутренним контрольным объёмам
		for (integer G=0; G<6; G++) {

			
			if (f.neighbors_for_the_internal_node[G][0][i] >= f.maxelm) {
				// это граничный узел
				integer inumber=f.neighbors_for_the_internal_node[G][0][i] - f.maxelm; // номер граничного узла
				//TOCHKA p; // координаты центра КО.
				//center_cord3D(i, f.nvtx, f.pa, p,100); // вычисление координат центра КО.
				// вычисление размеров текущего контрольного объёма:
	            //doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контрольного объёма
	           // volume3D(i, f.nvtx, f.pa, dx, dy, dz);

				if ((f.border_neighbor[inumber].MCB==(ls+lw))||
					(f.border_neighbor[inumber].MCB<ls)) 
				{
					bwall[inumber]=true;					
				}

				if ((f.border_neighbor[inumber].MCB<(ls + lw)) && 
					(f.border_neighbor[inumber].MCB >= ls) && 
					(!w[f.border_neighbor[inumber].MCB - ls].bsymmetry) && 
					(!w[f.border_neighbor[inumber].MCB - ls].bopening) &&
					(!w[f.border_neighbor[inumber].MCB - ls].bpressure)) 
				{
					doublereal Magnitude=fabs(w[f.border_neighbor[inumber].MCB-ls].Vx)+fabs(w[f.border_neighbor[inumber].MCB-ls].Vy)+fabs(w[f.border_neighbor[inumber].MCB-ls].Vz);
					if (Magnitude<admission) {
						// Стенка на которой задана нулевая скорость
						// Следовательно это твёрдая стенка
						bwall[inumber]=true;
					} // Magnitude
				} // WALL

				/*
				switch (G) {
					case ESIDE: rposition[X][inumber]=p.x+0.5*dx;
						     rposition[Y][inumber]=p.y;
							 rposition[Z][inumber]=p.z;
						     break;
					case WSIDE: rposition[X][inumber]=p.x-0.5*dx;
						     rposition[Y][inumber]=p.y;
							 rposition[Z][inumber]=p.z;
						     break;
					case NSIDE: rposition[X][inumber]=p.x;
						     rposition[Y][inumber]=p.y+0.5*dy;
							 rposition[Z][inumber]=p.z;
						     break;
					case SSIDE: rposition[X][inumber]=p.x;
						     rposition[Y][inumber]=p.y-0.5*dy;
							 rposition[Z][inumber]=p.z;
						     break;
					case TSIDE: rposition[X][inumber]=p.x;
						     rposition[Y][inumber]=p.y;
							 rposition[Z][inumber]=p.z+0.5*dz;
						     break;
					case BSIDE: rposition[X][inumber]=p.x;
						     rposition[Y][inumber]=p.y;
							 rposition[Z][inumber]=p.z-0.5*dz;
						     break;
					}
					*/
					rposition[X][inumber]= f.border_neighbor[inumber].p_c.x;
					rposition[Y][inumber]= f.border_neighbor[inumber].p_c.y;
					rposition[Z][inumber]= f.border_neighbor[inumber].p_c.z;
					
			}
		
			if (b_on_adaptive_local_refinement_mesh) {

				if (f.neighbors_for_the_internal_node[G][1][i] >= f.maxelm) {
					// это граничный узел
					integer inumber = f.neighbors_for_the_internal_node[G][1][i] - f.maxelm; // номер граничного узла
					//TOCHKA p; // координаты центра КО.
					//center_cord3D(i, f.nvtx, f.pa, p, 100); // вычисление координат центра КО.
															// вычисление размеров текущего контрольного объёма:
					//doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
					//volume3D(i, f.nvtx, f.pa, dx, dy, dz);

					if ((f.border_neighbor[inumber].MCB == (ls + lw)) ||
						(f.border_neighbor[inumber].MCB < ls))
					{
						bwall[inumber] = true;
					}

					if ((f.border_neighbor[inumber].MCB < (ls + lw)) &&
						(f.border_neighbor[inumber].MCB >= ls) &&
						(!w[f.border_neighbor[inumber].MCB - ls].bsymmetry) &&
						(!w[f.border_neighbor[inumber].MCB - ls].bopening) &&
						(!w[f.border_neighbor[inumber].MCB - ls].bpressure))
					{
						doublereal Magnitude = fabs(w[f.border_neighbor[inumber].MCB - ls].Vx) + fabs(w[f.border_neighbor[inumber].MCB - ls].Vy) + fabs(w[f.border_neighbor[inumber].MCB - ls].Vz);
						if (Magnitude < admission) {
							// Стенка на которой задана нулевая скорость
							// Следовательно это твёрдая стенка
							bwall[inumber] = true;
						} // Magnitude
					} // WALL

					rposition[X][inumber] = f.border_neighbor[inumber].p_c.x;
					rposition[Y][inumber] = f.border_neighbor[inumber].p_c.y;
					rposition[Z][inumber] = f.border_neighbor[inumber].p_c.z;
					//printf("%e %e %e\n", rposition[X][inumber], rposition[Y][inumber], rposition[Z][inumber]);
					//system("pause");

				}

				if (f.neighbors_for_the_internal_node[G][2][i] >= f.maxelm) {
					// это граничный узел
					integer inumber = f.neighbors_for_the_internal_node[G][2][i] - f.maxelm; // номер граничного узла
																		//TOCHKA p; // координаты центра КО.
																		//center_cord3D(i, f.nvtx, f.pa, p, 100); // вычисление координат центра КО.
																		// вычисление размеров текущего контрольного объёма:
																		//doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
																		//volume3D(i, f.nvtx, f.pa, dx, dy, dz);

					if ((f.border_neighbor[inumber].MCB == (ls + lw)) ||
						(f.border_neighbor[inumber].MCB < ls)) {
						bwall[inumber] = true;
					}

					if ((f.border_neighbor[inumber].MCB < (ls + lw)) &&
						(f.border_neighbor[inumber].MCB >= ls) &&
						(!w[f.border_neighbor[inumber].MCB - ls].bsymmetry) &&
						(!w[f.border_neighbor[inumber].MCB - ls].bopening) &&
						(!w[f.border_neighbor[inumber].MCB - ls].bpressure)) {
						doublereal Magnitude = fabs(w[f.border_neighbor[inumber].MCB - ls].Vx) + fabs(w[f.border_neighbor[inumber].MCB - ls].Vy) + fabs(w[f.border_neighbor[inumber].MCB - ls].Vz);
						if (Magnitude < admission) {
							// Стенка на которой задана нулевая скорость
							// Следовательно это твёрдая стенка
							bwall[inumber] = true;
						} // Magnitude
					} // WALL

					rposition[X][inumber] = f.border_neighbor[inumber].p_c.x;
					rposition[Y][inumber] = f.border_neighbor[inumber].p_c.y;
					rposition[Z][inumber] = f.border_neighbor[inumber].p_c.z;

					//printf("%e %e %e\n", rposition[X][inumber], rposition[Y][inumber], rposition[Z][inumber]);
					//system("pause");

				}

				if (f.neighbors_for_the_internal_node[G][3][i] >= f.maxelm) {
					// это граничный узел
					integer inumber = f.neighbors_for_the_internal_node[G][3][i] - f.maxelm; // номер граничного узла
																		//TOCHKA p; // координаты центра КО.
																		//center_cord3D(i, f.nvtx, f.pa, p, 100); // вычисление координат центра КО.
																		// вычисление размеров текущего контрольного объёма:
																		//doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
																		//volume3D(i, f.nvtx, f.pa, dx, dy, dz);

					if ((f.border_neighbor[inumber].MCB == (ls + lw)) ||
						(f.border_neighbor[inumber].MCB < ls))
					{
						bwall[inumber] = true;
					}

					if ((f.border_neighbor[inumber].MCB < (ls + lw)) &&
						(f.border_neighbor[inumber].MCB >= ls) &&
						(!w[f.border_neighbor[inumber].MCB - ls].bsymmetry) &&
						(!w[f.border_neighbor[inumber].MCB - ls].bopening) &&
						(!w[f.border_neighbor[inumber].MCB - ls].bpressure))
					{
						doublereal Magnitude = fabs(w[f.border_neighbor[inumber].MCB - ls].Vx) + fabs(w[f.border_neighbor[inumber].MCB - ls].Vy) + fabs(w[f.border_neighbor[inumber].MCB - ls].Vz);
						if (Magnitude < admission) {
							// Стенка на которой задана нулевая скорость
							// Следовательно это твёрдая стенка
							bwall[inumber] = true;
						} // Magnitude
					} // WALL

					rposition[X][inumber] = f.border_neighbor[inumber].p_c.x;
					rposition[Y][inumber] = f.border_neighbor[inumber].p_c.y;
					rposition[Z][inumber] = f.border_neighbor[inumber].p_c.z;

					//printf("%e %e %e\n", rposition[X][inumber], rposition[Y][inumber], rposition[Z][inumber]);
					//system("pause");

				}
			}
		}
	}

	// Непосредственное вычисление расстояния:
	// Только для внутренних КО.
#pragma omp parallel for
	for (integer i=0; i<f.maxelm; ++i) {
		TOCHKA p; // координаты центра КО.
		center_cord3D(i, f.nvtx, f.pa, p,100); // вычисление координат центра КО.
		for (integer j=0; j<f.maxbound; ++j) {
			if (bwall[j]) {
				doublereal dist=sqrt((p.x-rposition[X][j])*
					(p.x-rposition[X][j])+(p.y-rposition[Y][j])*
					(p.y-rposition[Y][j])+(p.z-rposition[Z][j])*
					(p.z-rposition[Z][j]));
				if (dist<f.rdistWall[i]) f.rdistWall[i]=dist;
			}
		}
	}

	// Непосредственное вычисление расстояния для граничных КО:
	// принадлежащих входной границе, выходной границе или границе с условиями симметрии.
#pragma omp parallel for
	for (integer i=f.maxelm; i<f.maxelm+f.maxbound; ++i) {
		integer inumber=i-f.maxelm; // номер граничного узла
		if (bwall[inumber]) {
			// расстояние равно нулю
			f.rdistWall[i]=0.0;
		}
		else {
            for (integer j=0; j<f.maxbound; ++j) {
		 	    if (bwall[j]) {
                    doublereal dist=sqrt((rposition[X][inumber]-rposition[X][j])*
						(rposition[X][inumber]-rposition[X][j])+
						(rposition[Y][inumber]-rposition[Y][j])*
						(rposition[Y][inumber]-rposition[Y][j])+
						(rposition[Z][inumber]-rposition[Z][j])*
						(rposition[Z][inumber]-rposition[Z][j]));
				    if (dist<f.rdistWall[i]) f.rdistWall[i]=dist;
			    }
			}
		}
	}

	// Вычисление толщины пограничного слоя в формуле Эскудиера
	//f.rdistWallmax=-1e20; // очень маленькое значение
	//for (integer i=0; i<(f.maxelm+f.maxbound); ++i) {
		//if (f.rdistWall[i]>f.rdistWallmax) f.rdistWallmax=f.rdistWall[i];
	//}

	f.rdistWallmax = get_max_array_elm(f.rdistWall, f.maxelm + f.maxbound);

	// Освобождение оперативной памяти.
	delete[] bwall;
	bwall=nullptr;
	for (integer i=0; i<3; ++i) {
		delete[] rposition[i];
		rposition[i]=nullptr;
	}
	delete[] rposition;
	rposition=nullptr;
} // calcdistwall

void constr_fluid_equation(FLOW* &f,  integer flow_interior,
	                       integer ls, integer lw, WALL* w) {

	// Если bCFXcalcdistancetotheWALL==true то будет использован
	// способ вычисления расстояния до ближайшей стенки используемый в ANSYS CFX.
	// А если bCFXcalcdistancetotheWALL==false то будет организован простой двойной цикл
    // вычисления кратчайшего расстояния до ближайшей стенки.
	bool bCFXcalcdistancetotheWALL=true;//true;

	if (b_on_adaptive_local_refinement_mesh) {
		// Обычный двойной цикл благо ячеек мало.
		bCFXcalcdistancetotheWALL = false;
	}

	if (flow_interior==0) {
		printf("f[0].maxelm=%d\n", f[0].maxelm);
		if (f[0].maxelm > 0) {
			if (eqin.imaxflD > 0) {
				// flow_interior==0 никак не согласуется с наличием жидких ячеек в сетке. 6.04.2019
				printf("Error: liquid zones should not be...\n");
				printf("Will be carried out of the program.\n");
				printf("please, press any key to exit...\n");
				//system("PAUSE");
				system("pause");
				exit(0); // выход из приложения
			}
		}
	}
	else {
		if (flow_interior!=eqin.imaxflD) {
			if ((eqin.itemper==1)&&(eqin.imaxflD==0)) {
				for (integer i=0; i<flow_interior; ++i) {
				    f[i].rdistWall=nullptr; // инициализация
					f[i].iflowregime= VISCOSITY_MODEL::LAMINAR; // VISCOSITY_MODEL::LAMINAR
				}
				printf("pure heat calculation.\n");
			}
			else {

			     printf("number of recognized program of liquid\n");
			     printf("isolated areas do not correspond to the\n");
			     printf("number of isolated areas of a user-defined liquid.\n");
			     printf("Will be carried out of the program.\n");
                 printf("please, press any key to exit...\n");
#if doubleintprecision == 1
				 printf("debug! flow_interior=%lld eqin.imaxflD=%lld\n", flow_interior, eqin.imaxflD);
#else
				 printf("debug! flow_interior=%d eqin.imaxflD=%d\n", flow_interior, eqin.imaxflD);
#endif
			     
			     //system("PAUSE");
				 system("pause");
				 exit(0); // выход из приложения
			}
		}
		else {

			integer i=0; // счётчик цикла for

			// всё в порядке количество изолированных жидких зон распознанных программой
			// равно количеству изолированных жидких зон заданных пользователем в интерфейсе.

			// нужно расположить жидкие зоны заданные пользователем в том порядке в котором они были 
			// распознаны программой. Это можно сделать с помощью матрицы расстояний от опорной точки заданной пользователем
			// до каждой программной жидкой зоны (минимальное расстояние).
			if (flow_interior==1) {
				// всё впорядке т.к. имеем только одну жидкую зону.
				
				if (eqin.fluidinfo[0].iflowregime== FLOW_REGIME::LAMINAR) {
					// Ламинарный режим
					f[0].iflowregime= VISCOSITY_MODEL::LAMINAR; // LAMINAR
				}
				else {
					if (eqin.fluidinfo[0].iturbmodel== TURBULENT_MODEL::ZEROEQMOD) {
                       f[0].iflowregime= VISCOSITY_MODEL::ZEROEQMOD; // ZEROEQMOD
					}
					if (eqin.fluidinfo[0].iturbmodel== TURBULENT_MODEL::SMAGORINSKY) {
                       f[0].iflowregime= VISCOSITY_MODEL::SMAGORINSKY; // SMAGORINSKY
					}
					if (eqin.fluidinfo[0].iturbmodel== TURBULENT_MODEL::RNG_LES) {
                       f[0].iflowregime= VISCOSITY_MODEL::RNG_LES; // RNG_LES
					}
					if (eqin.fluidinfo[0].iturbmodel == TURBULENT_MODEL::RANS_SPALART_ALLMARES) {
					   f[0].iflowregime = VISCOSITY_MODEL::RANS_SPALART_ALLMARES; // RANS_SPALART_ALLMARES
					}
					if (eqin.fluidinfo[0].iturbmodel == TURBULENT_MODEL::RANS_MENTER_SST) {
						f[0].iflowregime = VISCOSITY_MODEL::RANS_MENTER_SST; // K-Omega SST Menter
					}
					if (eqin.fluidinfo[0].iturbmodel == TURBULENT_MODEL::RANS_STANDART_K_EPS) {
						// Двухслойная модель на основе стандартной K-epsilon модели [2001]
						f[0].iflowregime = VISCOSITY_MODEL::RANS_STANDART_K_EPS;
					}
					if (eqin.fluidinfo[0].iturbmodel == TURBULENT_MODEL::RANS_LANGTRY_MENTOR_SST) {
						f[0].iflowregime = VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST; // Модель Ментора Лантгрии.
					}
					
					// заполнение параметров модели Смагоринского.
                       /*
					   // debug.
					   // свойства модели задаются явно прямо здесь. Отладочный вариант.
					   f[0].smaginfo.Cs=0.151; 
					   f[0].smaginfo.bfdelta=true; // учёт неравномерности сетки
					   f[0].smaginfo.bSmagorinsky_Lilly=true; // модель Смагоринского-Лиллу.
					   f[0].smaginfo.bsurface_roughness=false; // не учитывать шероховатость стенки.
					   f[0].smaginfo.ipowerroughness=2; // наиболее рекомендуемый вариант.
					   f[0].smaginfo.roughness=10e-6; // максимальная шероховатость поверхности 10 мкм.
					   */
					   f[0].smaginfo.Cs=eqin.fluidinfo[0].Cs; // постоянная Смагоринского.
					   // Использовать ли динамическую модель Германо для определения
					   // квадрата постоянной Смагоринского.
					   f[0].smaginfo.bDynamic_Stress=eqin.fluidinfo[0].bDynamic_Stress;
					   // использовать ли минимаксное ограничение на постоянную Смагоринского.
					   f[0].smaginfo.bLimiters_Cs=eqin.fluidinfo[0].bLimiters_Cs;
					   f[0].smaginfo.minCs=eqin.fluidinfo[0].rminCs; // минимальное значение постоянной Смагоринского
					   f[0].smaginfo.maxCs=eqin.fluidinfo[0].rmaxCs; // максимальное значение постоянной Смагоринского
					   // Тип фильтра который используется в модели Германо.
					   f[0].smaginfo.itypeFiltrGermano=eqin.fluidinfo[0].itypeFiltrGermano;
					   // учёт неравномерности сетки
					   if (eqin.fluidinfo[0].bfdelta==1) {
						   f[0].smaginfo.bfdelta=true; // включён. 
					   }
					   else {
						   f[0].smaginfo.bfdelta=false; // выключен.
					   }
					   // модель Смагоринского-Лиллу.
					   if (eqin.fluidinfo[0].bSmagorinskyLilly==1) {
						   f[0].smaginfo.bSmagorinsky_Lilly=true; // включён.
					   }
					   else {
						   f[0].smaginfo.bSmagorinsky_Lilly=false; // выключен.
					   }
					   // учёт шероховатости стенки.
					   if (eqin.fluidinfo[0].bsurface_roughness==1) {
						   f[0].smaginfo.bsurface_roughness=true; // учитывать шероховатость стенки.
					   }
					   else {
						   f[0].smaginfo.bsurface_roughness=false; // не учитывать шероховатость стенки.
					   }

					   f[0].smaginfo.bRichardsonCorrect=eqin.fluidinfo[0].bSwirlAmendment; // поправка на течения с кривизной линий тока.
					   f[0].smaginfo.bSelectiveSmagorinsky=eqin.fluidinfo[0].bSelectiveSmagorinsky; // Избирательность в модели Смагоринского.
					   
					   f[0].smaginfo.ipowerroughness=eqin.fluidinfo[0].ipowerroughness; // наиболее рекомендуемый вариант ipowerroughness==2.
					   f[0].smaginfo.roughness=eqin.fluidinfo[0].roughness; // максимальная шероховатость поверхности 10 мкм.

					   f[0].smaginfo.rRichardsonMultiplyer=eqin.fluidinfo[0].rRi_mult; // множитель корректирующий турбулентное число Ричардсона.
					   f[0].smaginfo.SSangle=eqin.fluidinfo[0].rSelectiveAngle; // пороговое значение угла в модели Selective Smagorinsky.
					   f[0].smaginfo.itypeFILTRSelectiveSmagorinsky=eqin.fluidinfo[0].itypeSelectiveSmagorinsky_filtr; // тип фильтра в модели Selective Smagorinsky.
				}

				// Вычисление вспомогательных величин для моделей турбулентности:
			    bool bdist=false; // вычислять ли расстояние до стенки.
			    switch (f[0].iflowregime) {
			      case VISCOSITY_MODEL::LAMINAR: bdist=false; break;  // Ламинарное течение.
			      case VISCOSITY_MODEL::ZEROEQMOD: bdist=true; break; // Zero Equation Turbulence model (RANS)
				  case VISCOSITY_MODEL::SMAGORINSKY: bdist=true; break; // Smagorinsky model (LES) (содержит как опцию модель Германо 1991 года).
				  case VISCOSITY_MODEL::RNG_LES: bdist=false; break; // RNG_LES Renormalization Group Theory (LES)
				  case VISCOSITY_MODEL::RANS_SPALART_ALLMARES: bdist = true;  break; // Spallart Allmares (RANS).
				  case VISCOSITY_MODEL::RANS_MENTER_SST: bdist = true; break; // K-Omega SST (RANS)
				  case VISCOSITY_MODEL::RANS_STANDART_K_EPS: bdist = true; break; // Стандартная KE модель (RANS) на основе двухслойной модели.
				  case VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST: bdist = true; break; // Модель ламинарно турбулентного перехода Ментора Лантгрии.
				}
				f[0].rdistWall=nullptr; // инициализация
			    if (bdist) {
				   // Вычисление расстояния до ближайшей стенки
				   printf("please wait! program calculate distance to the wall...\n");
				   if (!bCFXcalcdistancetotheWALL) {
					   // простой переборный алгоритм.
					   calcdistwall(f[0], ls, lw, w);
				   }
				   else {

					   if (number_cores() > 1) {
						   // Алгоритм на основе ICCG солвера не работает в параллельном режиме. 28.12.2021

						   // простой переборный алгоритм.
						   calcdistwall(f[0], ls, lw, w);
					   }
					   else {
						   // способ используемый в ANSYS CFX.
						   calcdistwallCFX(f[0], ls, lw, w);
					   }
				   }
			    }
			}
			else {
				// имеем две жидкие зоны или больше.
			    doublereal** distzone=nullptr;
				// distzone[i][j] - кратчайшее расстояние от опорной точки i-ого жидкого отсека пользователя до j-ой распознанной программой жидкой зоны.
				// выделение памяти оперативной памяти
				distzone=new doublereal*[flow_interior];
				if (distzone==nullptr) {
	                // недостаточно памяти на данном оборудовании.
		            printf("Problem: not enough memory on your equipment for distzone constr struct...\n");
		            printf("Please any key to exit...\n");
					//system("PAUSE");
					system("pause");
		            exit(1);
	            }
				for (i=0; i<flow_interior; ++i) {
					distzone[i]=nullptr;
					distzone[i]=new doublereal[flow_interior];
					if (distzone[i]==nullptr) {
	                     // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
						printf("Problem: not enough memory on your equipment for distzone[%lld] constr struct...\n", i);
#else
						printf("Problem: not enough memory on your equipment for distzone[%d] constr struct...\n", i);
#endif
		                
		                 printf("Please any key to exit...\n");
						 //system("PAUSE");
						 system("pause");
		                 exit(1);
	                }
				}
				// инициализация:
				for (i=0; i<flow_interior; ++i) for (integer j=0; j<flow_interior; ++j) distzone[i][j]=1e30; // присваиваем максимально большое число.
				// Вычисление расстояний:
				for (i=0; i<flow_interior; ++i) {
					for (integer j=0; j<flow_interior; ++j) {
						// только до внутренних контрольных объёмов
						TOCHKA p; // координаты центра КО.
						for (integer iP=0; iP<f[j].maxelm; ++iP) {
							center_cord3D(iP, f[j].nvtx, f[j].pa, p,100); // вычисление координат центра КО.
							distzone[i][j]=fmin(distzone[i][j],sqrt((eqin.fluidinfo[i].xc-p.x)*(eqin.fluidinfo[i].xc-p.x)+(eqin.fluidinfo[i].yc-p.y)*(eqin.fluidinfo[i].yc-p.y)+(eqin.fluidinfo[i].zc-p.z)*(eqin.fluidinfo[i].zc-p.z)));
						}
					}
				}
				integer *key=nullptr;
				key=new integer[flow_interior]; // ключ для сортировки
				if (key==nullptr) {
	                // недостаточно памяти на данном оборудовании.
		            printf("Problem: not enough memory on your equipment for key constr struct...\n");
		            printf("Please any key to exit...\n");
					//system("PAUSE");
					system("pause");
		            exit(1);
	            }
				// инициализация.
				for (i=0; i<flow_interior; ++i) key[i]=flow_interior; // присваиваем максимальное значение.
				for (i=0; i<flow_interior; ++i) {
					// поиск минимума в i-ой строке distzone
					doublereal rmin=1e30;
					integer imin=flow_interior;
					for (integer j=0; j<flow_interior; ++j) {
						if (distzone[i][j]<rmin) {
							rmin=distzone[i][j];
							imin=j;
						}
					}
					// минимальное расстояние в i-ой строке находится в ячейке imin
					key[i]=imin;
				}

				bool berror=false;
				for (i=0; i<flow_interior; ++i) if (key[i]>=flow_interior) {
					// ошибка распознавания жидких зон заданных пользоваателем.
					printf("Error detection of liquid zones defined by the user!\n");
					printf("Will be carried out of the program.\n");
					printf("Please, press any key to exit...\n");
					//system("PAUSE");
					system("pause");
					exit(0);
				}

				// сортировка массива eqin.fluidinfo по ключу key.
                // пузырьковая сортировка BubbleSort
				integer numberOfPairs=flow_interior;
				bool swappedElements=true;
				while (swappedElements) {
                    numberOfPairs--;
					swappedElements=false;
					for (i=0; i<numberOfPairs; ++i) {
						if (key[i]>key[i+1]) {
							// Swap i <-> i+1
							TFLOWINFO termflow;
							termflow=eqin.fluidinfo[i];
							eqin.fluidinfo[i]=eqin.fluidinfo[i+1];
							eqin.fluidinfo[i+1]=termflow;
							// swap key i <-> i+1
							integer ikeyterm;
							ikeyterm=key[i];
							key[i]=key[i+1];
							key[i+1]=ikeyterm;
							swappedElements=true;
						}
					}
				}

				// Проверка корректности, т.к. в массиве key[] не должно быть повторяющихся значений ключей.
				berror=false; // ошибки не обнаружено.
				for (i=0; i<(flow_interior-1); ++i) {
					// следует помнить что ключи key[] теперь отсортированы в порядке возрастания 
					if (key[i]==key[i+1]) berror=true;
				}

				if (berror) {
					// ошибка! две или более непересекающиеся по жидкости  жидкие зоны
					// ссылаются на одну и ту же обособленную жидкую зону распознанную программой. 
                    printf("Error! two or more non-overlapping zones \n");
					printf("of fluid liquid refer to the same\n");
					printf("isolated area of the liquid the recognized program.\n");
					printf("Will be carried out of the program.\n");
					printf("Please, press any key to exit...\n");
					//system("PAUSE");
					system("pause");
					exit(0);
				}

				// Освобождение оперативной памяти
				for (i=0; i<flow_interior; ++i) {
					delete[] distzone[i];
					distzone[i]=nullptr;
				}
				delete[] distzone;
				distzone=nullptr;
				delete[] key;
				key=nullptr;

				for (i=0; i<flow_interior; ++i) {
					if (eqin.fluidinfo[i].iflowregime== FLOW_REGIME::LAMINAR) {
					    // Ламинарный режим
					    f[i].iflowregime= VISCOSITY_MODEL::LAMINAR; // LAMINAR
				    }
				    else {
					    if (eqin.fluidinfo[i].iturbmodel== TURBULENT_MODEL::ZEROEQMOD) {
                            f[i].iflowregime= VISCOSITY_MODEL::ZEROEQMOD; // ZEROEQMOD
					    }
						if (eqin.fluidinfo[i].iturbmodel== TURBULENT_MODEL::SMAGORINSKY) {
                            f[i].iflowregime= VISCOSITY_MODEL::SMAGORINSKY; // SMAGORINSKY	
					    }
                        if (eqin.fluidinfo[i].iturbmodel== TURBULENT_MODEL::RNG_LES) {
                            f[i].iflowregime= VISCOSITY_MODEL::RNG_LES; // RNG_LES	
					    }						
						if (eqin.fluidinfo[i].iturbmodel == TURBULENT_MODEL::RANS_SPALART_ALLMARES) {
							f[i].iflowregime = VISCOSITY_MODEL::RANS_SPALART_ALLMARES; // RANS_SPALART_ALLMARES
						}
						if (eqin.fluidinfo[i].iturbmodel == TURBULENT_MODEL::RANS_MENTER_SST) {
							f[i].iflowregime = VISCOSITY_MODEL::RANS_MENTER_SST; // K-Omega SST
						}						
						if (eqin.fluidinfo[i].iturbmodel == TURBULENT_MODEL::RANS_STANDART_K_EPS) {
							// Двухслойная модель на основе стандартной k-epsilon модели [2001].
							f[i].iflowregime = VISCOSITY_MODEL::RANS_STANDART_K_EPS;
						}
						if (eqin.fluidinfo[i].iturbmodel == TURBULENT_MODEL::RANS_LANGTRY_MENTOR_SST) {
							f[i].iflowregime = VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST; // модель Ментора Лантгрии.
						}

						f[i].smaginfo.Cs=eqin.fluidinfo[i].Cs; // постоянная Смагоринского.
						// Использовать ли динамическую модель Германо для определения
					    // квадрата постоянной Смагоринского.
					    f[i].smaginfo.bDynamic_Stress=eqin.fluidinfo[i].bDynamic_Stress;
						// использовать ли минимаксное ограничение на постоянную Смагоринского.
					    f[i].smaginfo.bLimiters_Cs=eqin.fluidinfo[i].bLimiters_Cs;
					    f[i].smaginfo.minCs=eqin.fluidinfo[i].rminCs; // минимальное значение постоянной Смагоринского
					    f[i].smaginfo.maxCs=eqin.fluidinfo[i].rmaxCs; // максимальное значение постоянной Смагоринского
					    // Тип фильтра который используется в модели Германо.
					    f[i].smaginfo.itypeFiltrGermano=eqin.fluidinfo[i].itypeFiltrGermano;
					    // учёт неравномерности сетки
					    if (eqin.fluidinfo[i].bfdelta==1) {
						    f[i].smaginfo.bfdelta=true; // включён. 
					    }
					    else {
						    f[i].smaginfo.bfdelta=false; // выключен.
					    }
					    // модель Смагоринского-Лиллу.
					    if (eqin.fluidinfo[i].bSmagorinskyLilly==1) {
						    f[i].smaginfo.bSmagorinsky_Lilly=true; // включён.
					    }
					    else {
						    f[i].smaginfo.bSmagorinsky_Lilly=false; // выключен.
					    }
					    // учёт шероховатости стенки.
					    if (eqin.fluidinfo[i].bsurface_roughness==1) {
						    f[i].smaginfo.bsurface_roughness=true; // учитывать шероховатость стенки.
					    }
					    else {
						    f[i].smaginfo.bsurface_roughness=false; // не учитывать шероховатость стенки.
					    }
					   
						f[i].smaginfo.bRichardsonCorrect=eqin.fluidinfo[i].bSwirlAmendment; // поправка на течения с кривизной линий тока.
					    f[i].smaginfo.bSelectiveSmagorinsky=eqin.fluidinfo[i].bSelectiveSmagorinsky; // Избирательность в модели Смагоринского.

						f[i].smaginfo.ipowerroughness=eqin.fluidinfo[i].ipowerroughness; // наиболее рекомендуемый вариант ipowerroughness==2.
					    f[i].smaginfo.roughness=eqin.fluidinfo[i].roughness; // максимальная шероховатость поверхности 10 мкм.

						f[i].smaginfo.rRichardsonMultiplyer=eqin.fluidinfo[i].rRi_mult; // множитель корректирующий турбулентное число Ричардсона.
					    f[i].smaginfo.SSangle=eqin.fluidinfo[i].rSelectiveAngle; // пороговое значение угла в модели Selective Smagorinsky.
					    f[i].smaginfo.itypeFILTRSelectiveSmagorinsky=eqin.fluidinfo[i].itypeSelectiveSmagorinsky_filtr; // тип фильтра в модели Selective Smagorinsky.

				    }

				    // Вычисление вспомогательных величин для моделей турбулентности:
			        bool bdist=false; // вычислять ли расстояние до стенки.
			        switch (f[i].iflowregime) {
			           case VISCOSITY_MODEL::LAMINAR: bdist=false; break;
			           case VISCOSITY_MODEL::ZEROEQMOD: bdist=true; break; // Zero Equation Model (RANS)
					   case VISCOSITY_MODEL::SMAGORINSKY: bdist=true; break; // Smagorinsky model (LES) (содержит как опцию модель Германо 1991 года).
					   case VISCOSITY_MODEL::RNG_LES: bdist=false; break; // RNG_LES Renormalization Group Theory (LES)
					   case VISCOSITY_MODEL::RANS_SPALART_ALLMARES: bdist = true; break; // RANS Спалларт Аллмарес (Spallart Allmares)
					   case VISCOSITY_MODEL::RANS_MENTER_SST: bdist = true; break; // RANS k-Omega SST Menter [1993].
					   case VISCOSITY_MODEL::RANS_STANDART_K_EPS: bdist = true; break; // RANS двухслойная модель на основе стандартной KE модели
					   case VISCOSITY_MODEL::RANS_LANGTRY_MENTOR_SST: bdist = true; break; // Модель ламинарно турбулентного перехода Ментора Лантгрии.
					}
					f[i].rdistWall=nullptr; // инициализация
			        if (bdist) {
				       // Вычисление расстояния до ближайшей стенки
				       printf("please wait! program calculate distance to the wall...\n");
					   if (!bCFXcalcdistancetotheWALL) {
					        // простой переборный алгоритм.
					        calcdistwall(f[i], ls, lw, w);
				       }
				       else {
					        // способ используемый в ANSYS CFX.
					        calcdistwallCFX(f[i], ls, lw, w);
				       }
			        }
				}
				
			}
		}
	}
} // constr_fluid_equation

// Создаёт необходимые связи для модели излучения внутри блока.
void constr_link_on_surface_for_radiation_model(integer maxelm, int* &whot_is_block, int*** neighbors_for_the_internal_node,
	                                            int** nvtx, TOCHKA* pa,  BLOCK* &b, integer lb)
{
	if (print_log_message) {
		std::cout << "number threads =" << number_cores() << std::endl;
	}

#ifdef _OPENMP 
	omp_set_num_threads(number_cores()); // установка числа потоков
#endif

	integer* lb_size = new integer[lb];
#pragma omp parallel for
	for (integer i = 0; i < lb; ++i)
	{
		lb_size[i] = 0;
	}

	
	for (integer iP = 0; iP < maxelm; ++iP) {
		integer ib = whot_is_block[iP];
		if (b[ib].radiation.binternalRadiation) {
			//isize++;
			lb_size[ib]++;
		}
	}

	integer isize = 0;

#pragma omp parallel for
	for (integer i = 0; i < lb; ++i)
	{
		if (lb_size[i] > isize) {
#pragma omp critical
			{
				if (lb_size[i] > isize) {
					isize = lb_size[i];
				}
			}
		}

	}


	isize += 2;

	delete[] lb_size;

	MY_PAIR* list = new MY_PAIR[isize];

	MY_PAIR** list_omp = new MY_PAIR*[number_cores()];

	integer* list_candidate = new integer[isize];


	integer** list_candidateomp = new integer*[number_cores()];

	for (integer i = 0; i < number_cores(); ++i) {
		list_candidateomp[i] = new integer[isize];
		list_omp[i]= new MY_PAIR[isize];

#pragma omp parallel for
		for (integer iP = 0; iP < isize; ++iP) {
			list_candidateomp[i][iP] = -1;
		}
	}

	for (integer i = 0; i < lb; ++i) {
		if (b[i].radiation.binternalRadiation) {
			integer itek = i;

			integer icandidate_limit = 0;
			integer* icandidate_limitomp = new integer[number_cores()];


			for (integer i = 0; i < number_cores(); ++i) {
				icandidate_limitomp[i] = 0;
			}

#pragma omp parallel for
			for (integer iP = 0; iP < maxelm; ++iP) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
						int tid = omp_get_thread_num();
						list_candidateomp[tid][icandidate_limitomp[tid]] = iP;
						icandidate_limitomp[tid]++;
				}
			}

			for (integer i = 0; i < number_cores(); ++i) {
#pragma omp parallel for
				for (integer j = 0; j < icandidate_limitomp[i]; ++j) 
				{
					list_candidate[icandidate_limit + j] = list_candidateomp[i][j];
				}
				icandidate_limit += icandidate_limitomp[i];
			}

			

			delete[] icandidate_limitomp;

#pragma omp parallel sections
			{

#pragma omp section
				{

					integer ic = 0;
					int tid = omp_get_thread_num();

					//for (integer iP = 0; iP < maxelm; ++iP)
					for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
					{
						integer iP = list_candidate[iP_1];
						//integer ib = whot_is_block[iP];
						//if (ib == itek) 
						{
							//integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
							//iE = neighbors_for_the_internal_node[E_SIDE][0][iP];
							//iN = neighbors_for_the_internal_node[N_SIDE][0][iP];
							//iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
							integer iW = neighbors_for_the_internal_node[W_SIDE][0][iP];
							//iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; 
							//iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

							// вычисление размеров текущего контрольного объёма:
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
							volume3D_q(iP, nvtx, pa, dx, dy, dz);

							if (iW >= maxelm) {
								// граничный узел.
								list_omp[tid][ic].node1 = iP;
								list_omp[tid][ic].node21 = iW;
								list_omp[tid][ic].node22 = -1;
								list_omp[tid][ic].node23 = -1;
								list_omp[tid][ic].node24 = -1;
								//list_omp[tid][ic].dS = dy*dz; // Площадь ячейки.
								list_omp[tid][ic].dS1 = dy * dz;
								list_omp[tid][ic].dS2 = 0.0;
								list_omp[tid][ic].dS3 = 0.0;
								list_omp[tid][ic].dS4 = 0.0;
								ic++;
							}
							else
							{
								// Внутренняя грань, но соседи принадлежат разным блокам.
								if (itek != whot_is_block[iW]) {
									list_omp[tid][ic].node1 = iP;
									list_omp[tid][ic].node21 = iW;
									list_omp[tid][ic].node22 = -1;
									list_omp[tid][ic].node23 = -1;
									list_omp[tid][ic].node24 = -1;
									//list_omp[tid][ic].dS = dy*dz; // Площадь ячейки.
									list_omp[tid][ic].dS1 = dy * dz;
									list_omp[tid][ic].dS2 = 0.0;
									list_omp[tid][ic].dS3 = 0.0;
									list_omp[tid][ic].dS4 = 0.0;
									ic++;
								}
							}

						}
					}

#pragma omp critical
					{
						b[itek].radiation.nodelistW = new MY_PAIR[ic];
					}
					b[itek].radiation.nodelistWsize = ic;
					for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistW[j] = list_omp[tid][j];

				}

#pragma omp section
				{
					integer ic = 0;
					int tid = omp_get_thread_num();


					//for (integer iP = 0; iP < maxelm; ++iP)
					for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
					{
						integer iP = list_candidate[iP_1];
						//integer ib = whot_is_block[iP];
						//if (ib == itek)
						{
							//integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
							integer iE = neighbors_for_the_internal_node[E_SIDE][0][iP];
							//iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; 
							//iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
							//iW = neighbors_for_the_internal_node[W_SIDE][0][iP];
							//iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; 
							//iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

							// вычисление размеров текущего контрольного объёма:
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
							volume3D_q(iP, nvtx, pa, dx, dy, dz);

							if (iE >= maxelm) {
								list_omp[tid][ic].node1 = iP;
								list_omp[tid][ic].node21 = iE;
								list_omp[tid][ic].node22 = -1;
								list_omp[tid][ic].node23 = -1;
								list_omp[tid][ic].node24 = -1;
								//list_omp[tid][ic].dS = dy*dz; // Площадь ячейки.
								list_omp[tid][ic].dS1 = dy * dz;
								list_omp[tid][ic].dS2 = 0.0;
								list_omp[tid][ic].dS3 = 0.0;
								list_omp[tid][ic].dS4 = 0.0;
								ic++;
							}
							else
							{
								if (itek != whot_is_block[iE]) {
									list_omp[tid][ic].node1 = iP;
									list_omp[tid][ic].node21 = iE;
									list_omp[tid][ic].node22 = -1;
									list_omp[tid][ic].node23 = -1;
									list_omp[tid][ic].node24 = -1;
									//list_omp[tid][ic].dS = dy*dz; // Площадь ячейки.
									list_omp[tid][ic].dS1 = dy * dz;
									list_omp[tid][ic].dS2 = 0.0;
									list_omp[tid][ic].dS3 = 0.0;
									list_omp[tid][ic].dS4 = 0.0;
									ic++;
								}
							}

						}
					}

#pragma omp critical
					{
						b[itek].radiation.nodelistE = new MY_PAIR[ic];
					}
					b[itek].radiation.nodelistEsize = ic;
					for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistE[j] = list_omp[tid][j];
				}

#pragma omp section
				{

					integer ic = 0;
					int tid = omp_get_thread_num();

					//for (integer iP = 0; iP < maxelm; ++iP) 
					for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
					{
						integer iP = list_candidate[iP_1];
						//integer ib = whot_is_block[iP];
						//if (ib == itek)
						{
							//integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
							//iE = neighbors_for_the_internal_node[E_SIDE][0][iP];
							//iN = neighbors_for_the_internal_node[N_SIDE][0][iP];
							//iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
							//iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; 
							integer iS = neighbors_for_the_internal_node[S_SIDE][0][iP];
							//iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

							// вычисление размеров текущего контрольного объёма:
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
							volume3D_q(iP, nvtx, pa, dx, dy, dz);

							if (iS >= maxelm) {
								list_omp[tid][ic].node1 = iP;
								list_omp[tid][ic].node21 = iS;
								list_omp[tid][ic].node22 = -1;
								list_omp[tid][ic].node23 = -1;
								list_omp[tid][ic].node24 = -1;
								//list_omp[tid][ic].dS = dx*dz; // Площадь ячейки.
								list_omp[tid][ic].dS1 = dx * dz;
								list_omp[tid][ic].dS2 = 0.0;
								list_omp[tid][ic].dS3 = 0.0;
								list_omp[tid][ic].dS4 = 0.0;
								ic++;
							}
							else
							{
								if (itek != whot_is_block[iS]) {
									list_omp[tid][ic].node1 = iP;
									list_omp[tid][ic].node21 = iS;
									list_omp[tid][ic].node22 = -1;
									list_omp[tid][ic].node23 = -1;
									list_omp[tid][ic].node24 = -1;
									//list_omp[tid][ic].dS = dx*dz; // Площадь ячейки.
									list_omp[tid][ic].dS1 = dx * dz;
									list_omp[tid][ic].dS2 = 0.0;
									list_omp[tid][ic].dS3 = 0.0;
									list_omp[tid][ic].dS4 = 0.0;
									ic++;
								}
							}

						}
					}

#pragma omp critical
					{
						b[itek].radiation.nodelistS = new MY_PAIR[ic];
					}
					b[itek].radiation.nodelistSsize = ic;
					for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistS[j] = list_omp[tid][j];

				}

#pragma omp section
				{
					integer ic = 0;
					int tid = omp_get_thread_num();

					//for (integer iP = 0; iP < maxelm; ++iP)
					for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
					{
						integer iP = list_candidate[iP_1];
						//integer ib = whot_is_block[iP];
						//if (ib == itek)
						{
							// integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
							//iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; 
							integer iN = neighbors_for_the_internal_node[N_SIDE][0][iP];
							//iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
							//iW = neighbors_for_the_internal_node[W_SIDE][0][iP];
							//iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; 
							//iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

							// вычисление размеров текущего контрольного объёма:
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
							volume3D_q(iP, nvtx, pa, dx, dy, dz);

							if (iN >= maxelm) {
								list_omp[tid][ic].node1 = iP;
								list_omp[tid][ic].node21 = iN;
								list_omp[tid][ic].node22 = -1;
								list_omp[tid][ic].node23 = -1;
								list_omp[tid][ic].node24 = -1;
								//list_omp[tid][ic].dS = dx*dz; // Площадь ячейки.
								list_omp[tid][ic].dS1 = dx * dz;
								list_omp[tid][ic].dS2 = 0.0;
								list_omp[tid][ic].dS3 = 0.0;
								list_omp[tid][ic].dS4 = 0.0;
								ic++;
							}
							else
							{
								if (itek != whot_is_block[iN]) {
									list_omp[tid][ic].node1 = iP;
									list_omp[tid][ic].node21 = iN;
									list_omp[tid][ic].node22 = -1;
									list_omp[tid][ic].node23 = -1;
									list_omp[tid][ic].node24 = -1;
									//list_omp[tid][ic].dS = dx*dz; // Площадь ячейки.
									list_omp[tid][ic].dS1 = dx * dz;
									list_omp[tid][ic].dS2 = 0.0;
									list_omp[tid][ic].dS3 = 0.0;
									list_omp[tid][ic].dS4 = 0.0;
									ic++;
								}
							}

						}
					}

#pragma omp critical
					{
						b[itek].radiation.nodelistN = new MY_PAIR[ic];
					}
					b[itek].radiation.nodelistNsize = ic;
					for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistN[j] = list_omp[tid][j];

				}

#pragma omp section
				{
					integer ic = 0;
					int tid = omp_get_thread_num();

					//for (integer iP = 0; iP < maxelm; ++iP)
					for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
					{
						integer iP = list_candidate[iP_1];
						//integer ib = whot_is_block[iP];
						//if (ib == itek)
						{
							//integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
							//iE = neighbors_for_the_internal_node[E_SIDE][0][iP];
							//iN = neighbors_for_the_internal_node[N_SIDE][0][iP]; 
							//iT = neighbors_for_the_internal_node[T_SIDE][0][iP]; 
							//iW = neighbors_for_the_internal_node[W_SIDE][0][iP]; 
							//iS = neighbors_for_the_internal_node[S_SIDE][0][iP]; 
							integer iB = neighbors_for_the_internal_node[B_SIDE][0][iP];

							// вычисление размеров текущего контрольного объёма:
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
							volume3D_q(iP, nvtx, pa, dx, dy, dz);

							if (iB >= maxelm) {
								list_omp[tid][ic].node1 = iP;
								list_omp[tid][ic].node21 = iB;
								list_omp[tid][ic].node22 = -1;
								list_omp[tid][ic].node23 = -1;
								list_omp[tid][ic].node24 = -1;
								//list_omp[tid][ic].dS = dx*dy; // Площадь ячейки.
								list_omp[tid][ic].dS1 = dx * dy;
								list_omp[tid][ic].dS2 = 0.0;
								list_omp[tid][ic].dS3 = 0.0;
								list_omp[tid][ic].dS4 = 0.0;
								ic++;
							}
							else
							{
								if (itek != whot_is_block[iB]) {
									list_omp[tid][ic].node1 = iP;
									list_omp[tid][ic].node21 = iB;
									list_omp[tid][ic].node22 = -1;
									list_omp[tid][ic].node23 = -1;
									list_omp[tid][ic].node24 = -1;
									//list_omp[tid][ic].dS = dx*dy; // Площадь ячейки.
									list_omp[tid][ic].dS1 = dx * dy;
									list_omp[tid][ic].dS2 = 0.0;
									list_omp[tid][ic].dS3 = 0.0;
									list_omp[tid][ic].dS4 = 0.0;
									ic++;
								}
							}

						}
					}

#pragma omp critical
					{
						b[itek].radiation.nodelistB = new MY_PAIR[ic];
					}
					b[itek].radiation.nodelistBsize = ic;
					for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistB[j] = list_omp[tid][j];
				}

#pragma omp section
				{
					integer ic = 0;
					int tid = omp_get_thread_num();


					//for (integer iP = 0; iP < maxelm; ++iP) 
					for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
					{
						integer iP = list_candidate[iP_1];
						//integer ib = whot_is_block[iP];
						//if (ib == itek)
						{
							//integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
							//iE = neighbors_for_the_internal_node[E_SIDE][0][iP]; 
							//iN = neighbors_for_the_internal_node[N_SIDE][0][iP];
							integer iT = neighbors_for_the_internal_node[T_SIDE][0][iP];
							//iW = neighbors_for_the_internal_node[W_SIDE][0][iP];
							//iS = neighbors_for_the_internal_node[S_SIDE][0][iP];
							//iB = neighbors_for_the_internal_node[B_SIDE][0][iP];


							// вычисление размеров текущего контрольного объёма:
							doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
							volume3D_q(iP, nvtx, pa, dx, dy, dz);

							if (iT >= maxelm) {
								list_omp[tid][ic].node1 = iP;
								list_omp[tid][ic].node21 = iT;
								list_omp[tid][ic].node22 = -1;
								list_omp[tid][ic].node23 = -1;
								list_omp[tid][ic].node24 = -1;
								//list_omp[tid][ic].dS = dx*dy; // Площадь ячейки.
								list_omp[tid][ic].dS1 = dx * dy;
								list_omp[tid][ic].dS2 = 0.0;
								list_omp[tid][ic].dS3 = 0.0;
								list_omp[tid][ic].dS4 = 0.0;
								ic++;
							}
							else
							{
								if (itek != whot_is_block[iT]) {
									list_omp[tid][ic].node1 = iP;
									list_omp[tid][ic].node21 = iT;
									list_omp[tid][ic].node22 = -1;
									list_omp[tid][ic].node23 = -1;
									list_omp[tid][ic].node24 = -1;
									//list_omp[tid][ic].dS = dx*dy; // Площадь ячейки.
									list_omp[tid][ic].dS1 = dx * dy;
									list_omp[tid][ic].dS2 = 0.0;
									list_omp[tid][ic].dS3 = 0.0;
									list_omp[tid][ic].dS4 = 0.0;
									ic++;
								}
							}

						}
					}

#pragma omp critical
					{
						b[itek].radiation.nodelistT = new MY_PAIR[ic];
					}
					b[itek].radiation.nodelistTsize = ic;
					for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistT[j] = list_omp[tid][j];
				}
			}

		}
	}

	delete[] list;
	list = nullptr;

	delete[] list_candidate;
	list_candidate = nullptr;

	for (integer i = 0; i < number_cores(); ++i) {
		delete[] list_candidateomp[i];
		delete[] list_omp[i];
	}
	delete[] list_candidateomp;
	delete[] list_omp;

} // constr_link_on_surface_for_radiation_model


// Начало написания 14 сентября 2016. 15_06.
// Продолжение написания 20 сентября 2016. 14_45.
// Правка 12 июня 2019.
// Также необходимо менять код в модуле my_solver_v0_03.cpp,
// и в модуле my_elmatr_quad.cpp
// Создаёт необходимые связи для модели излучения внутри блока.
void constr_link_on_surface_for_radiation_model_alice(integer maxelm, BOUND* &border_neighbor, int* &whot_is_block, int*** neighbors_for_the_internal_node,
	int** nvtx, TOCHKA* pa, BLOCK* &b, integer lb)
{

	integer isize = 0;
	for (integer iP = 0; iP < maxelm; ++iP) {
		integer ib = whot_is_block[iP];
		if (b[ib].radiation.binternalRadiation) {
			isize++;
		}
	}
	isize += 2;

	MY_PAIR* list = new MY_PAIR[isize];
	integer* list_candidate = new integer[maxelm];

	for (integer i = 0; i < lb; ++i) {
		if (b[i].radiation.binternalRadiation) {
			integer itek = i;

			integer icandidate_limit = 0;

			for (integer iP = 0; iP < maxelm; ++iP) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
					list_candidate[icandidate_limit] = iP;
					icandidate_limit++;
				}
			}

			
			integer ic = 0;

			//for (integer iP = 0; iP < maxelm; ++iP)
			for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
			{
				integer iP = list_candidate[iP_1];
				//integer ib = whot_is_block[iP];
				//if (ib == itek)
				{
					//integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					//integer iE2, iN2, iT2, iW2, iS2, iB2;
					//integer iE3, iN3, iT3, iW3, iS3, iB3;
					//integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					//iE1 = neighbors_for_the_internal_node[E_SIDE][0][iP]; 
					//iN1 = neighbors_for_the_internal_node[N_SIDE][0][iP];
					//iT1 = neighbors_for_the_internal_node[T_SIDE][0][iP]; 
					integer iW1 = neighbors_for_the_internal_node[W_SIDE][0][iP]; 
					//iS1 = neighbors_for_the_internal_node[S_SIDE][0][iP];
					//iB1 = neighbors_for_the_internal_node[B_SIDE][0][iP];

					//iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP];
					//iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP];
					//iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
					integer iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP]; 
					//iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP];
					//iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

					//iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP];
					//iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP];
					//iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
					integer iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; 
					//iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP];
					//iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

					//iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP];
					//iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP];
					//iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
					integer iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP];
					//iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP];
					//iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					bool binc = false;
					list[ic].node21 = -1;
					list[ic].node22 = -1;
					list[ic].node23 = -1;
					list[ic].node24 = -1;
					list[ic].dS1 = 0.0;
					list[ic].dS2 = 0.0;
					list[ic].dS3 = 0.0;
					list[ic].dS4 = 0.0;

					if (iW1 == -1) {
						//printf("iW1 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iW1 != -1) {

						if (iW1 >= maxelm) {
							// граничный узел.
							list[ic].node1 = iP;
							list[ic].node21 = iW1;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							//list[ic].dS1 = dy*dz; // Площадь ячейки.
							list[ic].dS1 = border_neighbor[iW1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							// Внутренняя грань, но соседи принадлежат разным блокам.
							if (itek != whot_is_block[iW1]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iW1, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node21 = iW1;
								//list[ic].dS = dy*dz; // Площадь ячейки.
								list[ic].dS1 = dy_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iW2 == -1) {
						//printf("iW2 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iW2 != -1) {
						if (iW2 >= maxelm) {
							// граничный узел.
							list[ic].node1 = iP;
							list[ic].node22 = iW2;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							//list[ic].dS2 = dy*dz; // Площадь ячейки.
							list[ic].dS2 = border_neighbor[iW2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							// Внутренняя грань, но соседи принадлежат разным блокам.
							if (itek != whot_is_block[iW2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iW2, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node22 = iW2;
								//list[ic].dS = dy*dz; // Площадь ячейки.
								list[ic].dS2 = dy_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iW3 == -1) {
						//printf("iW3 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iW3 != -1) {
						if (iW3 >= maxelm) {
							// граничный узел.
							list[ic].node1 = iP;
							list[ic].node23 = iW3;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							//list[ic].dS3 = dy*dz; // Площадь ячейки.
							list[ic].dS3 = border_neighbor[iW3 - maxelm].dS; // Площадь ячейки.

							binc = true;
						}
						else
						{
							// Внутренняя грань, но соседи принадлежат разным блокам.
							if (itek != whot_is_block[iW3]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iW3, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node23 = iW3;
								//list[ic].dS = dy*dz; // Площадь ячейки.
								list[ic].dS3 = dy_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iW4 == -1) {
						//printf("iW4 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iW4 != -1) {
						if (iW4 >= maxelm) {
							// граничный узел.
							list[ic].node1 = iP;
							list[ic].node24 = iW4;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							//list[ic].dS4 = dy*dz; // Площадь ячейки.
							list[ic].dS4 = border_neighbor[iW4 - maxelm].dS; // Площадь ячейки.

							binc = true;
						}
						else
						{
							// Внутренняя грань, но соседи принадлежат разным блокам.
							if (itek != whot_is_block[iW4]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iW4, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node24 = iW4;
								//list[ic].dS = dy*dz; // Площадь ячейки.
								list[ic].dS4 = dy_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (binc) {
						ic++;
					}
				}
			}

			b[itek].radiation.nodelistW = new MY_PAIR[ic];
			b[itek].radiation.nodelistWsize = ic;
			for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistW[j] = list[j];

			ic = 0;

			//for (integer iP = 0; iP < maxelm; ++iP) 
			for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
			{
				integer iP = list_candidate[iP_1];
				//integer ib = whot_is_block[iP];
				//if (ib == itek)
				{
				

					//integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					//integer iE2, iN2, iT2, iW2, iS2, iB2;
					//integer iE3, iN3, iT3, iW3, iS3, iB3;
					//integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					integer iE1 = neighbors_for_the_internal_node[E_SIDE][0][iP];
					//iN1 = neighbors_for_the_internal_node[N_SIDE][0][iP]; 
					//iT1 = neighbors_for_the_internal_node[T_SIDE][0][iP];
					//iW1 = neighbors_for_the_internal_node[W_SIDE][0][iP];
					//iS1 = neighbors_for_the_internal_node[S_SIDE][0][iP]; 
					//iB1 = neighbors_for_the_internal_node[B_SIDE][0][iP];

					integer iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; 
					//iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP];
					//iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
					//iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP];
					//iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP];
					//iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

					integer iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP];
					//iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP]; 
					//iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
					//iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP];
					//iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; 
					//iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

					integer iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP];
					//iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP]; 
					//iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
					//iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; 
					//iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP]; 
					//iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];


					bool binc = false;
					list[ic].node21 = -1;
					list[ic].node22 = -1;
					list[ic].node23 = -1;
					list[ic].node24 = -1;
					list[ic].dS1 = 0.0;
					list[ic].dS2 = 0.0;
					list[ic].dS3 = 0.0;
					list[ic].dS4 = 0.0;

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iE1 == -1) {
						//printf("iE1 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iE1 != -1) {
						if (iE1 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node21 = iE1;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							//list[ic].dS1 = dy*dz; // Площадь ячейки.
							list[ic].dS1 = border_neighbor[iE1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iE1]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iE1, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node21 = iE1;
								//list[ic].dS = dy*dz; // Площадь ячейки.
								list[ic].dS1 = dy_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iE2 == -1) {
						//printf("iE2 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iE2 != -1) {
						if (iE2 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node22 = iE2;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							//list[ic].dS2 = dy*dz; // Площадь ячейки.
							list[ic].dS2 = border_neighbor[iE2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iE2]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iE2, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node22 = iE2;
								//list[ic].dS = dy*dz; // Площадь ячейки.
								list[ic].dS2 = dy_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iE3 == -1) {
						//printf("iE3 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iE3 != -1) {
						if (iE3 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node23 = iE3;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							//list[ic].dS3 = dy*dz; // Площадь ячейки.
							list[ic].dS3 = border_neighbor[iE3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iE3]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iE3, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node23 = iE3;
								//list[ic].dS = dy*dz; // Площадь ячейки.
								list[ic].dS3 = dy_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iE4 == -1) {
						//printf("iE4 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iE4 != -1) {
						if (iE4 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node24 = iE4;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							//list[ic].dS4 = dy*dz; // Площадь ячейки.
							list[ic].dS4 = border_neighbor[iE4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iE4]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iE4, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node24 = iE4;
								//list[ic].dS = dy*dz; // Площадь ячейки.
								list[ic].dS4 = dy_loc * dz_loc; // Площадь ячейки
								binc = true;
							}
						}
					}

					if (binc) {
						ic++;
					}

				}
			}

			b[itek].radiation.nodelistE = new MY_PAIR[ic];
			b[itek].radiation.nodelistEsize = ic;
			for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistE[j] = list[j];

			ic = 0;

			//for (integer iP = 0; iP < maxelm; ++iP)
			for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
			{
				integer iP = list_candidate[iP_1];
				//integer ib = whot_is_block[iP];
				//if (ib == itek)
				{
					
					//integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					//integer iE2, iN2, iT2, iW2, iS2, iB2;
					//integer iE3, iN3, iT3, iW3, iS3, iB3;
					//integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					//iE1 = neighbors_for_the_internal_node[E_SIDE][0][iP];
					//iN1 = neighbors_for_the_internal_node[N_SIDE][0][iP];
					//iT1 = neighbors_for_the_internal_node[T_SIDE][0][iP];
					//iW1 = neighbors_for_the_internal_node[W_SIDE][0][iP]; 
					integer iS1 = neighbors_for_the_internal_node[S_SIDE][0][iP];
					//iB1 = neighbors_for_the_internal_node[B_SIDE][0][iP];

					//iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP];
					//iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP];
					//iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
					//iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP];
					integer iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP];
					//iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

					//iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP];
					//iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP];
					//iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
					//iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP];
					integer iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP];
					//iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

					//iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP];
					//iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP];
					//iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
					//iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP];
					integer iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP];
					//iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];


					bool binc = false;
					list[ic].node21 = -1;
					list[ic].node22 = -1;
					list[ic].node23 = -1;
					list[ic].node24 = -1;
					list[ic].dS1 = 0.0;
					list[ic].dS2 = 0.0;
					list[ic].dS3 = 0.0;
					list[ic].dS4 = 0.0;

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iS1 == -1) {
						//printf("iS1 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iS1 != -1) {
						if (iS1 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node21 = iS1;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							//list[ic].dS1 = dx*dz; // Площадь ячейки.
							list[ic].dS1 = border_neighbor[iS1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iS1]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iS1, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node21 = iS1;
								//list[ic].dS = dx*dz; // Площадь ячейки.
								list[ic].dS1 = dx_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iS2 == -1) {
						//printf("iS2 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iS2 != -1) {
						if (iS2 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node22 = iS2;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							//list[ic].dS2 = dx*dz; // Площадь ячейки.
							list[ic].dS2 = border_neighbor[iS2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iS2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iS2, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node22 = iS2;
								//list[ic].dS = dx*dz; // Площадь ячейки.
								list[ic].dS2 = dx_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iS3 == -1) {
						//printf("iS3 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iS3 != -1) {
						if (iS3 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node23 = iS3;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							//list[ic].dS3 = dx*dz; // Площадь ячейки.
							list[ic].dS3 = border_neighbor[iS3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iS3]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iS3, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node23 = iS3;
								//list[ic].dS = dx*dz; // Площадь ячейки.
								list[ic].dS3 = dx_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iS4 == -1) {
						//printf("iS4 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iS4 != -1) {
						if (iS4 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node24 = iS4;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							//list[ic].dS4 = dx*dz; // Площадь ячейки.
							list[ic].dS4 = border_neighbor[iS4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iS4]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iS4, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node24 = iS4;
								//list[ic].dS = dx*dz; // Площадь ячейки.
								list[ic].dS4 = dx_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (binc) {
						ic++;
					}

				}
			}

			b[itek].radiation.nodelistS = new MY_PAIR[ic];
			b[itek].radiation.nodelistSsize = ic;
			for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistS[j] = list[j];

			
			ic = 0;

			//for (integer iP = 0; iP < maxelm; ++iP) 
			for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
			{
				integer iP = list_candidate[iP_1];
				//integer ib = whot_is_block[iP];
				//if (ib == itek)
				{
					
					//integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					//integer iE2, iN2, iT2, iW2, iS2, iB2;
					//integer iE3, iN3, iT3, iW3, iS3, iB3;
					//integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					//iE1 = neighbors_for_the_internal_node[E_SIDE][0][iP];
					integer iN1 = neighbors_for_the_internal_node[N_SIDE][0][iP]; 
					//iT1 = neighbors_for_the_internal_node[T_SIDE][0][iP];
					//iW1 = neighbors_for_the_internal_node[W_SIDE][0][iP];
					//iS1 = neighbors_for_the_internal_node[S_SIDE][0][iP]; 
					//iB1 = neighbors_for_the_internal_node[B_SIDE][0][iP];

					//iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; 
					integer iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP]; 
					//iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
					//iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP];
					//iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; 
					//iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

					//iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP];
					integer iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP];
					//iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
					//iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP];
					//iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; 
					//iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

					//iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP]; 
					integer iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP];
					//iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
					//iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP];
					//iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP];
					//iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];


					bool binc = false;
					list[ic].node21 = -1;
					list[ic].node22 = -1;
					list[ic].node23 = -1;
					list[ic].node24 = -1;
					list[ic].dS1 = 0.0;
					list[ic].dS2 = 0.0;
					list[ic].dS3 = 0.0;
					list[ic].dS4 = 0.0;


					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);


					if (iN1 == -1) {
						//printf("iN1 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iN1 != -1) {
						if (iN1 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node21 = iN1;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							//list[ic].dS1 = dx*dz; // Площадь ячейки.
							list[ic].dS1 = border_neighbor[iN1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iN1]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iN1, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node21 = iN1;
								//list[ic].dS = dx*dz; // Площадь ячейки.
								list[ic].dS1 = dx_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iN2 == -1) {
						//printf("iN2 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iN2 != -1) {
						if (iN2 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node22 = iN2;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							//list[ic].dS2 = dx*dz; // Площадь ячейки.
							list[ic].dS2 = border_neighbor[iN2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iN2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iN2, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node22 = iN2;
								//list[ic].dS = dx*dz; // Площадь ячейки.
								list[ic].dS2 = dx_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iN3 == -1) {
						//printf("iN3 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iN3 != -1) {
						if (iN3 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node23 = iN3;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							//list[ic].dS3 = dx*dz; // Площадь ячейки.
							list[ic].dS3 = border_neighbor[iN3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iN3]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iN3, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node23 = iN3;
								//list[ic].dS = dx*dz; // Площадь ячейки.
								list[ic].dS3 = dx_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iN4 == -1) {
						//printf("iN4 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iN4 != -1) {
						if (iN4 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node24 = iN4;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							//list[ic].dS4 = dx*dz; // Площадь ячейки.
							list[ic].dS4 = border_neighbor[iN4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iN4]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iN4, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node24 = iN4;
								//list[ic].dS = dx*dz; // Площадь ячейки.
								list[ic].dS4 = dx_loc * dz_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (binc) {
						ic++;
					}

				}
			}

			b[itek].radiation.nodelistN = new MY_PAIR[ic];
			b[itek].radiation.nodelistNsize = ic;
			for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistN[j] = list[j];

			
			ic = 0;

			//for (integer iP = 0; iP < maxelm; ++iP) 
			for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
			{
				integer iP = list_candidate[iP_1];
				//integer ib = whot_is_block[iP];
				//if (ib == itek)
				{
				
					//integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					//integer iE2, iN2, iT2, iW2, iS2, iB2;
					//integer iE3, iN3, iT3, iW3, iS3, iB3;
					//integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					//iE1 = neighbors_for_the_internal_node[E_SIDE][0][iP];
					//iN1 = neighbors_for_the_internal_node[N_SIDE][0][iP];
					//iT1 = neighbors_for_the_internal_node[T_SIDE][0][iP];
					//iW1 = neighbors_for_the_internal_node[W_SIDE][0][iP];
					//iS1 = neighbors_for_the_internal_node[S_SIDE][0][iP];
					integer iB1 = neighbors_for_the_internal_node[B_SIDE][0][iP];

					//iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP];
					//iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP];
					//iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
					//iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP];
					//iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP];
					integer iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

					//iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP];
					//iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP];
					//iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
					//iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; 
					//iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP];
					integer iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

					//iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP];
					//iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP];
					//iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
					//iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP];
					//iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP];
					integer iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];


					bool binc = false;
					list[ic].node21 = -1;
					list[ic].node22 = -1;
					list[ic].node23 = -1;
					list[ic].node24 = -1;
					list[ic].dS1 = 0.0;
					list[ic].dS2 = 0.0;
					list[ic].dS3 = 0.0;
					list[ic].dS4 = 0.0;


					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iB1 == -1) {
						//printf("iB1 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iB1 != -1) {
						if (iB1 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node21 = iB1;
							//list[ic].dS1 = dx*dy; // Площадь ячейки.
							list[ic].dS1 = border_neighbor[iB1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iB1]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iB1, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node21 = iB1;
								list[ic].dS1 = dx_loc * dy_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iB2 == -1) {
						//printf("iB2 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iB2 != -1) {
						if (iB2 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node22 = iB2;
							//list[ic].dS2 = dx*dy; // Площадь ячейки.
							list[ic].dS2 = border_neighbor[iB2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iB2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iB2, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node22 = iB2;
								list[ic].dS2 = dx_loc * dy_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iB3 == -1) {
						//printf("iB3 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iB3 != -1) {
						if (iB3 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node23 = iB3;
							//list[ic].dS3 = dx*dy; // Площадь ячейки.
							list[ic].dS3 = border_neighbor[iB3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iB3]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iB3, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node23 = iB3;
								list[ic].dS3 = dx_loc * dy_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iB4 == -1) {
						//printf("iB4 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iB4 != -1) {
						if (iB4 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node24 = iB4;
							//list[ic].dS4 = dx*dy; // Площадь ячейки.
							list[ic].dS4 = border_neighbor[iB4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iB4]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iB4, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node24 = iB4;
								list[ic].dS4 = dx_loc * dy_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (binc) {
						ic++;
					}

				}
			}

			b[itek].radiation.nodelistB = new MY_PAIR[ic];
			b[itek].radiation.nodelistBsize = ic;
			for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistB[j] = list[j];

			ic = 0;

			//for (integer iP = 0; iP < maxelm; ++iP)
			for (integer iP_1 = 0; iP_1 < icandidate_limit; iP_1++)
			{
				integer iP = list_candidate[iP_1];
				//integer ib = whot_is_block[iP];
				//if (ib == itek)
				{
				
					//integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					//integer iE2, iN2, iT2, iW2, iS2, iB2;
					//integer iE3, iN3, iT3, iW3, iS3, iB3;
					//integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					//iE1 = neighbors_for_the_internal_node[E_SIDE][0][iP];
					//iN1 = neighbors_for_the_internal_node[N_SIDE][0][iP];
					integer iT1 = neighbors_for_the_internal_node[T_SIDE][0][iP];
					//iW1 = neighbors_for_the_internal_node[W_SIDE][0][iP];
					//iS1 = neighbors_for_the_internal_node[S_SIDE][0][iP];
					//iB1 = neighbors_for_the_internal_node[B_SIDE][0][iP];

					//iE2 = neighbors_for_the_internal_node[E_SIDE][1][iP]; 
					//iN2 = neighbors_for_the_internal_node[N_SIDE][1][iP];
					integer iT2 = neighbors_for_the_internal_node[T_SIDE][1][iP];
					//iW2 = neighbors_for_the_internal_node[W_SIDE][1][iP];
					//iS2 = neighbors_for_the_internal_node[S_SIDE][1][iP]; 
					//iB2 = neighbors_for_the_internal_node[B_SIDE][1][iP];

					//iE3 = neighbors_for_the_internal_node[E_SIDE][2][iP];
					//iN3 = neighbors_for_the_internal_node[N_SIDE][2][iP];
					integer iT3 = neighbors_for_the_internal_node[T_SIDE][2][iP];
					//iW3 = neighbors_for_the_internal_node[W_SIDE][2][iP]; 
					//iS3 = neighbors_for_the_internal_node[S_SIDE][2][iP]; 
					//iB3 = neighbors_for_the_internal_node[B_SIDE][2][iP];

					//iE4 = neighbors_for_the_internal_node[E_SIDE][3][iP];
					//iN4 = neighbors_for_the_internal_node[N_SIDE][3][iP];
					integer iT4 = neighbors_for_the_internal_node[T_SIDE][3][iP];
					//iW4 = neighbors_for_the_internal_node[W_SIDE][3][iP]; 
					//iS4 = neighbors_for_the_internal_node[S_SIDE][3][iP];
					//iB4 = neighbors_for_the_internal_node[B_SIDE][3][iP];


					bool binc = false;
					list[ic].node21 = -1;
					list[ic].node22 = -1;
					list[ic].node23 = -1;
					list[ic].node24 = -1;
					list[ic].dS1 = 0.0;
					list[ic].dS2 = 0.0;
					list[ic].dS3 = 0.0;
					list[ic].dS4 = 0.0;

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iT1 == -1) {
						//printf("iT1 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iT1 != -1) {
						if (iT1 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node21 = iT1;
							//list[ic].dS1 = dx*dy; // Площадь ячейки.
							list[ic].dS1 = border_neighbor[iT1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iT1]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iT1, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node21 = iT1;
								list[ic].dS1 = dx_loc * dy_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iT2 == -1) {
						//printf("iT2 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iT2 != -1) {
						if (iT2 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node22 = iT2;
							//list[ic].dS2 = dx*dy; // Площадь ячейки.
							list[ic].dS2 = border_neighbor[iT2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iT2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iT2, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node22 = iT2;
								list[ic].dS2 = dx_loc * dy_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iT3 == -1) {
						//printf("iT3 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iT3 != -1) {
						if (iT3 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node23 = iT3;
							//list[ic].dS3 = dx*dy; // Площадь ячейки.
							list[ic].dS3 = border_neighbor[iT3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iT3]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iT3, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node23 = iT3;
								list[ic].dS3 = dx_loc * dy_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (iT4 == -1) {
						//printf("iT4 == -1\n");
						//system("PAUSE");
						//system("PAUSE");
					}

					if (iT4 != -1) {
						if (iT4 >= maxelm) {
							list[ic].node1 = iP;
							list[ic].node24 = iT4;
							//list[ic].dS4 = dx*dy; // Площадь ячейки.
							list[ic].dS4 = border_neighbor[iT4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iT4]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контрольного объёма
								volume3D(iT4, nvtx, pa, dx_loc, dy_loc, dz_loc);

								list[ic].node1 = iP;
								list[ic].node24 = iT4;
								list[ic].dS4 = dx_loc * dy_loc; // Площадь ячейки.
								binc = true;
							}
						}
					}

					if (binc) {
						ic++;
					}

				}
			}

			b[itek].radiation.nodelistT = new MY_PAIR[ic];
			b[itek].radiation.nodelistTsize = ic;
			for (integer j = 0; j < ic; ++j) b[itek].radiation.nodelistT[j] = list[j];

		}
	}

	delete[] list;
	list = nullptr;
	delete[] list_candidate;
	list_candidate = nullptr;

} // constr_link_on_surface_for_radiation_model_alice


  // освобождение памяти уровня 1
  // для уравнений гидродинамики.
void free_level1_flow(FLOW* &fglobal, int &flow_interior) {
	
	for (integer iflow = 0; iflow<flow_interior; iflow++) {
		
		if (fglobal != nullptr) {

			if (print_log_message) {
#if doubleintprecision == 1
				printf("delete flow %lld icolor_diferent_fluid_domain\n", iflow);
#else
				printf("delete flow %d icolor_diferent_fluid_domain\n", iflow);
#endif
			}

			if (fglobal[iflow].whot_is_block != nullptr) {
				delete[] fglobal[iflow].whot_is_block;
				fglobal[iflow].whot_is_block = nullptr;
			}

			delete[] fglobal[iflow].center_coord;
			delete[] fglobal[iflow].volume;

			if (fglobal[iflow].icolor_different_fluid_domain != nullptr) {
				delete[] fglobal[iflow].icolor_different_fluid_domain;
				fglobal[iflow].icolor_different_fluid_domain = nullptr;
			}

			if (print_log_message) {
#if doubleintprecision == 1
				printf("delete flow %lld neighbors in boundary nodes\n", iflow);
				// border_neighbor
#else
				printf("delete flow %d neighbors in boundary nodes\n", iflow);
#endif
			}

			if (fglobal[iflow].border_neighbor != nullptr) {
				delete[] fglobal[iflow].border_neighbor;
				fglobal[iflow].border_neighbor = nullptr;
			}

			if (print_log_message) {
#if doubleintprecision == 1
				printf("delete flow %lld rdistWall\n", iflow);
#else
				printf("delete flow %d rdistWall\n", iflow);
#endif
			}

			if (fglobal[iflow].rdistWall != nullptr) { // -1N
				delete[] fglobal[iflow].rdistWall;
				fglobal[iflow].rdistWall = nullptr;
			}

			if (print_log_message) {
#if doubleintprecision == 1
				printf("delete flow %lld neighbors in internal nodes\n", iflow);
				// neighbors_for_the_internal_node
#else
				printf("delete flow %d neighbors in internal nodes\n", iflow);
#endif
			}

			if (fglobal[iflow].neighbors_for_the_internal_node != nullptr) {
				for (integer i = 0; i < 12; ++i) {
					if (fglobal[iflow].neighbors_for_the_internal_node[i] != nullptr) {
						delete[] fglobal[iflow].neighbors_for_the_internal_node[i]; // -12N
					}
				}
			}
			if (fglobal[iflow].neighbors_for_the_internal_node != nullptr) {
				delete[] fglobal[iflow].neighbors_for_the_internal_node;
				fglobal[iflow].neighbors_for_the_internal_node = nullptr;
			}

			if (print_log_message) {
#if doubleintprecision == 1
				printf("delete flow %lld nvtx\n", iflow);
#else
				printf("delete flow %d nvtx\n", iflow);
#endif
			}

			for (integer i = 0; i < 8; ++i) { // -8N
#if doubleintprecision == 1
											//printf("%lld\n",i); system("PAUSE"); // debug
#else
											//printf("%d\n",i); system("PAUSE"); // debug
#endif

				if (fglobal[iflow].nvtx != nullptr) {
					if (fglobal[iflow].nvtx[i] != nullptr) {
						delete[] fglobal[iflow].nvtx[i];
					}
				}

			}
			if (fglobal[iflow].nvtx != nullptr) {
				delete[] fglobal[iflow].nvtx;
				fglobal[iflow].nvtx = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld properties in boundary nodes\n", iflow);
				// prop_b
#else
				printf("delete flow %d properties in boundary nodes\n", iflow);
#endif
			}

			if (fglobal[iflow].prop_b != nullptr) {
				for (integer i = 0; i < 3; ++i) {
					if (fglobal[iflow].prop_b[i] != nullptr) {
						delete[] fglobal[iflow].prop_b[i]; // -3N
					}
				}
			}
			if (fglobal[iflow].prop_b != nullptr) {
				delete[] fglobal[iflow].prop_b;
				fglobal[iflow].prop_b = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld properties in internal nodes\n", iflow);
				// prop
#else
				printf("delete flow %d properties in internal nodes\n", iflow);
#endif
		}
			
			if (fglobal[iflow].prop != nullptr) {
				for (integer i = 0; i<3; ++i) {
					if (fglobal[iflow].prop[i] != nullptr) {
						delete[] fglobal[iflow].prop[i]; // -3N
					}
				}
			}
			if (fglobal[iflow].prop != nullptr) {
				delete[] fglobal[iflow].prop;
				fglobal[iflow].prop = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld mass flux\n", iflow);
				// mf
#else
				printf("delete flow %d mass flux\n", iflow);
#endif
			}
			
			if (fglobal[iflow].mf != nullptr) {
				for (integer i = 0; i<fglobal[iflow].maxelm; ++i) {
					if (fglobal[iflow].mf[i] != nullptr) {
						delete[] fglobal[iflow].mf[i]; // -6N
						fglobal[iflow].mf[i] = nullptr;
					}
				}
			}
			if (fglobal[iflow].mf != nullptr) {
				delete[] fglobal[iflow].mf;
				fglobal[iflow].mf = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld SInvariantStrainRateTensor\n", iflow);
#else
				printf("delete flow %d SInvariantStrainRateTensor\n", iflow);
#endif
			}
			
			if (fglobal[iflow].SInvariantStrainRateTensor != nullptr) { // -1N
				delete[] fglobal[iflow].SInvariantStrainRateTensor;
				fglobal[iflow].SInvariantStrainRateTensor = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld diagonal coefficients\n", iflow);
				// diag_coef
#else
				printf("delete flow %d diagonal coefficients\n", iflow);
#endif
			}
			
			if (fglobal[iflow].diag_coef != nullptr) {
				for (integer i = 0; i<3; ++i) {
					if (fglobal[iflow].diag_coef[i] != nullptr) {
						delete[] fglobal[iflow].diag_coef[i]; // -3N
						fglobal[iflow].diag_coef[i] = nullptr;
					}
				}
			}
			if (fglobal[iflow].diag_coef != nullptr) {
				delete[] fglobal[iflow].diag_coef;
				fglobal[iflow].diag_coef = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld pa\n", iflow);
#else
				printf("delete flow %d pa\n", iflow);
#endif
			}
			
			if (fglobal[iflow].pa != nullptr) { // -3N
				delete[] fglobal[iflow].pa;
				fglobal[iflow].pa = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld ptr\n", iflow);
#else
				printf("delete flow %d ptr\n", iflow);
#endif
			}
			
			if (fglobal[iflow].ptr != nullptr) { // -N
				delete[] fglobal[iflow].ptr;
				fglobal[iflow].ptr = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld iN\n", iflow);
#else
				printf("delete flow %d iN\n", iflow);
#endif
			}
			
			if (fglobal[iflow].iN != nullptr) {
				for (integer i = 0; i<3; ++i) {
					if (fglobal[iflow].iN[i] != nullptr) {
						delete fglobal[iflow].iN[i]; // -max3(iWE, iSN, iBT)
						fglobal[iflow].iN[i] = nullptr;
					}
				}
			}
			if (fglobal[iflow].iN != nullptr) {
				delete fglobal[iflow].iN;
				fglobal[iflow].iN = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld id\n", iflow);
#else
				printf("delete flow %d id\n", iflow);
#endif
			}
			
			if (fglobal[iflow].id != nullptr) {
				for (integer i = 0; i<3; ++i) {
					if (fglobal[iflow].id[i] != nullptr) {
						integer n1 = 0;
						switch (i) {
						case 0: n1 = fglobal[iflow].iWE;  break;
						case 1: n1 = fglobal[iflow].iSN;  break;
						case 2: n1 = fglobal[iflow].iBT;  break;
						}
						for (integer j = 0; j<n1; ++j) {
							if (fglobal[iflow].id[i][j] != nullptr) {
								delete fglobal[iflow].id[i][j];
								fglobal[iflow].id[i][j] = nullptr;
							}
						}
					}
				}
			}
			if (fglobal[iflow].id != nullptr) {
				for (integer i = 0; i<3; ++i) {
					if (fglobal[iflow].id[i] != nullptr) {
						delete fglobal[iflow].id[i]; // 
						fglobal[iflow].id[i] = nullptr;
					}
				}
			}
			if (fglobal[iflow].id != nullptr) {
				delete fglobal[iflow].id;
				fglobal[iflow].id = nullptr;
			}

		}
	}
} // free_level1_flow


  // освобождение памяти уровня 2
  // для уравнений гидродинамики.
void free_level2_flow(FLOW* &fglobal, int &flow_interior) {
	
	for (integer iflow = 0; iflow<flow_interior; iflow++) {
		integer j = 0;
		if (fglobal != nullptr) {

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow%lld slau\n", iflow);
#else
				printf("delete flow%d slau\n", iflow);
#endif
			}
			
			if (fglobal[iflow].slau != nullptr) {
				for (j = 0; j<12; ++j) {
					if (fglobal[iflow].slau[j] != nullptr) {
						delete[] fglobal[iflow].slau[j];
						fglobal[iflow].slau[j] = nullptr;
					}
				}
				delete[] fglobal[iflow].slau;
				fglobal[iflow].slau = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow%lld slau_bon\n", iflow);
#else
				printf("delete flow%d slau_bon\n", iflow);
#endif
			}
			
			if (fglobal[iflow].slau_bon != nullptr) {
				for (j = 0; j<12; ++j) {
					if (fglobal[iflow].slau_bon[j] != nullptr) {
						delete[] fglobal[iflow].slau_bon[j];
						fglobal[iflow].slau_bon[j] = nullptr;
					}
				}
				delete[] fglobal[iflow].slau_bon;
				fglobal[iflow].slau_bon = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete fluid %lld potent\n", iflow);
#else
				printf("delete fluid %d potent\n", iflow);
#endif
			}
		
			if (fglobal[iflow].potent != nullptr) { // -26N
				for (integer i = 0; i<SIZE_FLOW_POTENT_ARRAY; ++i) {
					if (fglobal[iflow].potent[i] != nullptr) {
						delete[] fglobal[iflow].potent[i];
						fglobal[iflow].potent[i] = nullptr;
					}
				}
				delete[] fglobal[iflow].potent;
				fglobal[iflow].potent = nullptr;
			}

			if (print_log_message) {

#if doubleintprecision == 1
				printf("delete flow %lld alpha\n", iflow);
#else
				printf("delete flow %d alpha\n", iflow);
#endif
			}
			
			if (fglobal[iflow].alpha != nullptr) { // -5
				delete[] fglobal[iflow].alpha;
				fglobal[iflow].alpha = nullptr;
			}

		}

	}
	//flow_interior=0;
} // free_level2_flow


bool bfirst_additional_count_CAD = true;

/* Заполняет структуру данных для задачи теплопроводности TEMPER
* и структуры данных для задачи гидродинамики i*FLOW. 
* Возможны варианты: 
* 1. Отсеков гидродинамики вообще нет, т.е.
*    везде только SOLID.
* 2. Отсек гидродинамики есть но он пассивен, т.е. нулевые
*    гидродинамические компоненты есть решение. В этом случае
*    нужна специальная галочка что решать не нужно. Этот 
*    случай возникает когда нет вынужденной конвекции и когда 
*    нет естественной конвекции.
* 3. Есть отсек с вынужденной конвекцией он должен решаться 
*    в первую очередь.  
* 4. Есть отсек с естественной конвекцией. Тогда возникает 
*    взаимосогласованная задача.
*/
void load_TEMPER_and_FLOW(TEMPER& t, FLOW*& f, int& inx, int& iny, int& inz,
	doublereal*& xpos, doublereal*& ypos, doublereal*& zpos, int& flow_interior,
	BLOCK*& b, int lb, int lw, WALL*& w, SOURCE*& s, int ls,
	int lu, UNION*& my_union, doublereal temp_ref,
	TPROP*& matlist, bool bextendedprint,
	doublereal dgx, doublereal dgy, doublereal dgz, bool bALICEflag,
	bool breconstruct, integer& iunion_id_p1) {

	bool b_first_constr_neighbors_for_the_internal_node = true;

	// eqin - информация о решаемом наборе уравнений.

	integer icount_part = 1; // 19 счётчик хода выполнения

	integer ipolygon_count_limit = 0;

	integer iPnodes_count_shadow = 0, inx_shadow = 0, iny_shadow = 0, inz_shadow = 0;

	//инициализация.
	inx_shadow = inx;
	iny_shadow = iny;
	inz_shadow = inz;
	iPnodes_count_shadow = iPnodes_count_shadow_memo;

RESTART_CONSTRUCT: // возвращаемся к перестроению всех структур данных:

	ipolygon_count_limit++;

	if (print_log_message) {

		//system("pause");
		std::cout << "iunion_id_p1 = " << iunion_id_p1 << std::endl;

#if doubleintprecision == 1
		printf("inx=%d, iny=%d, inz=%d. Nodes=%lld\n", inx + 1, iny + 1, inz + 1, static_cast<integer>((inx + 1) * (iny + 1) * (inz + 1)));
#else
		printf("inx=%d, iny=%d, inz=%d. Nodes=%d\n", inx + 1, iny + 1, inz + 1, (inx + 1) * (iny + 1) * (inz + 1));
#endif
	}

	inx_shadow = inx;
	iny_shadow = iny;
	inz_shadow = inz;


	//for (integer i_53 = 0; i_53 <= inz; i_53++) {
#if doubleintprecision == 1
	//printf("i=%lld zpos[%lld]=%e ", i_53, i_53, zpos[i_53]);
#else
	//printf("i=%d zpos[%d]=%e ", i_53, i_53, zpos[i_53]);
#endif

		//if (i_53 < inz) printf("h=%e",zpos[i_53+1]-zpos[i_53]);
		//printf("\n");
	//}
	//system("PAUSE");

	// ТЕПЛОПРОВОДНОСТЬ:
	t.ilevel_alice = nullptr;

	int* evt_t = nullptr; // глобальная нумерация контрольных объёмов
	TOCKA_SHORT_INT* tck_int_list = nullptr;

	if (!bALICEflag) {
		enumerate_volume(evt_t, t.maxelm, TEMPERATURE, xpos, ypos, zpos, t.whot_is_block, inx, iny, inz, b, lb, lu, my_union, iunion_id_p1, tck_int_list, w, lw, s, ls);
		t.ilevel_alice = new integer[t.maxelm];
		// При конформной сетке всюду один и тот-же уровень - первый.
		for (integer i_3 = 0; i_3 < t.maxelm; ++i_3) {
			t.ilevel_alice[i_3] = 1;
		}
	}
	else {
		// evt_t отсутствует.
		// Построено t.maxelm, t.whot_is_block. Замена enumerate_volume.
		// ent_t отсутствует.
		// pa[0..maxnod-1]
		// Построено t.pa, t.maxnod. Замена constr_nodes.
		// nvtx[0..7][0..maxelm-1]
		// Построено t.nvtx. Замена constr_nvtx.
		// Построено t.prop, t.Sc, t.ipower_time_depend. Замена constr_prop.
		// Это закодировано 8 сентября 2016.
		// Осуществляет уникальную нумерацию листьев принадлежащих модели. Замена constr_neighbour.
		// Оставляет пустым neighbour=nullptr.
		// Эта функция за раз вычисляет части 1-6. part 1..6.
		constr_nodes_nvtx_prop_alice(oc_global, inx, iny, inz, t.maxelm, xpos, ypos, zpos, TEMPERATURE, b, lb, t.whot_is_block, t.pa, t.maxnod, t.nvtx,
			t.prop, t.Sc, t.ipower_time_depend, matlist, t.ilevel_alice, t.bActiveShearModule);

		if (b_on_adaptive_local_refinement_mesh) {
			if (print_log_message) {

				// Инвариант корректности АЛИС сетки.
				printf("the invariant correctness 8 november 2016...\n");
				printf("ANES_ALICE_CORRECT in constr_struct_alice.cpp file\n");
			}
			ANES_ALICE_CORRECT(t.maxnod, t.pa, t.maxelm, t.nvtx);
		}
	}
	if (print_log_message) {

		//printf("part 1\n");
		std::cout << "part " << icount_part++ << " enumerate_volume" << std::endl; //22
	}

	int* ent_t = nullptr; // глобальная нумерация узлов.
	if (!bALICEflag) {
		// pa[0..maxnod-1]
		constr_nodes(t.pa, t.maxnod, ent_t, TEMPERATURE, t.whot_is_block, evt_t, inx, iny, inz, xpos, ypos, zpos, b, lb, tck_int_list, t.maxelm);
	}

	if (print_log_message) {

		//printf("part 2\n");    
		std::cout << "part " << icount_part++ << " constr_nodes" << std::endl; //22
	}

	int** neighbour = nullptr;
	if (!bALICEflag) {
		/*
		#ifdef _OPENMP
				int inumcor = number_cores();
				omp_set_num_threads(inumcor);
				//omp_set_num_threads(3); // установка числа потоков

		#pragma omp parallel
				{

		#pragma omp parallel sections
					{
		#pragma omp section
						{
							// находит соседей только среди внутренних КО
							// для каждого внутреннего контрольного
							// объёма или 0 если соседа нет.
							// neighbour[0..11][0..maxelm-1]
							constr_neighbour(evt_t, ent_t, neighbour, t.maxelm, inx, iny, inz, tck_int_list);
						}

		#pragma omp section
						{
							// для каждого контрольного объёма принадлежащему
							// расчётной области определяет номера его вершин.
							// nvtx[0..7][0..maxelm-1]
							// Нумерация начинается с единицы.
							// maxelm передается уже посчитанным заранее, его не надо вычислять заново.
							constr_nvtx(evt_t, ent_t, t.nvtx, t.maxelm, inx, iny, inz, tck_int_list);
						}
		#pragma omp section
						{
							// Чистая гидродинамика будет экспортироваться если поле температур равно константе.
							// Поле температур равно константе если нет источников тепла, и имеется только одна температура Tamb.
							// Также необходимо чтобы в этом случае было ненулевое вынужденное течение. При соблюдении всех этих условий
							// экспортировать нужно с последним параметром равным 2. Аналогично если активны и температурные эффекты и
							// гидродинамика то нужно экспортировать с последним параметром равным 3. В целях уменьшения расхода оперативной
							// памяти нужно экспортировать части 1 и 3 прямо здесь и сейчас, т.к. позже память под необходимые для экспорта
							// структуры данных будет уничтожена. Но на данный момент ничего неизвестно о втором и третьем случае, т.е. экспортом
							// с последним параметром 2 и 3. Поэтому принято решение во время экспорта последней (средняя часть итогового файла)
							// части заменить вторую строку первого файла, отвечающую за перечень экспортируемых переменных.
							//printf("part 5\n");

							// Заносит свойства материалов в структуру prop для внутренних КО.
							// prop[0..2][0..maxelm-1]
							constr_prop(evt_t, t.whot_is_block, ent_t, t.prop, t.maxelm, TEMPERATURE, b, lb, inx, iny, inz, t.Sc, t.ipower_time_depend, xpos, ypos, zpos, matlist, tck_int_list, t.bActiveShearModule);
						}
					}
				}

				//omp_set_num_threads(1); // установка числа потоков
		#else
		*/
		// один поток
		// находит соседей только среди внутренних КО 
		// для каждого внутреннего контрольного
		// объёма или 0 если соседа нет.
		// neighbour[0..11][0..maxelm-1]
		constr_neighbour(evt_t, ent_t, neighbour, t.maxelm, inx, iny, inz, tck_int_list);
		// для каждого контрольного объёма принадлежащему
		// расчётной области определяет номера его вершин.
		// nvtx[0..7][0..maxelm-1]
		// Нумерация начинается с единицы.
		// maxelm передается уже посчитанным заранее, его не надо вычислять заново.
		constr_nvtx(evt_t, ent_t, t.nvtx, t.maxelm, inx, iny, inz, tck_int_list);

		// Чистая гидродинамика будет экспортироваться если поле температур равно константе.
		// Поле температур равно константе если нет источников тепла, и имеется только одна температура Tamb.
		// Также необходимо чтобы в этом случае было ненулевое вынужденное течение. При соблюдении всех этих условий
		// экспортировать нужно с последним параметром равным 2. Аналогично если активны и температурные эффекты и 
		// гидродинамика то нужно экспортировать с последним параметром равным 3. В целях уменьшения расхода оперативной 
		// памяти нужно экспортировать части 1 и 3 прямо здесь и сейчас, т.к. позже память под необходимые для экспорта
		// структуры данных будет уничтожена. Но на данный момент ничего неизвестно о втором и третьем случае, т.е. экспортом 
		// с последним параметром 2 и 3. Поэтому принято решение во время экспорта последней (средняя часть итогового файла) 
		// части заменить вторую строку первого файла, отвечающую за перечень экспортируемых переменных.
		//printf("part 5\n");

		// Заносит свойства материалов в структуру prop для внутренних КО.
		// prop[0..2][0..maxelm-1]
		constr_prop(evt_t, t.whot_is_block, ent_t, t.prop, t.maxelm, TEMPERATURE, b, lb, inx, iny, inz, t.Sc, t.ipower_time_depend, xpos, ypos, zpos, matlist, tck_int_list, t.bActiveShearModule);

		//#endif
	}

	if (print_log_message) {

		//printf("part 3\n");
		std::cout << "part " << icount_part++ << " constr_neighbours" << std::endl; //22


		//if (!bALICEflag) {

		//}
		//printf("part 4\n");
		std::cout << "part " << icount_part++ << " constr_nvtx" << std::endl; //22

	}

	//if (!bALICEflag) {

	//}
	//printf("part 6\n");
	if (print_log_message) {

#if doubleintprecision == 1
		printf("part %lld constr_prop\n", icount_part++); //22
#else
		printf("part %d constr_prop\n", icount_part++); //22
#endif
	}

	/*
	*  neighbors_for_the_internal_node, border_neighbor ( i - internal, b - boundary).
	*/

	// Алгоритм уникальной нумерации граней.
	// В результате каждая граничная грань,
	// принадлежащая расчётной области имеет
	// уникальный номер.
	int** gran_t = nullptr;
	if (!bALICEflag) {
		enumerate_gran_temp(gran_t, t.maxelm, t.nvtx, t.maxbound, neighbour, t.pa, s, ls);
	}
	else {
		// calculate_max_bound_temp закодирована 10 сентября 2016. в 10_34.
		// gran_t==nullptr не используется.
		calculate_max_bound_temp(oc_global, t.maxbound, (inx + 1) * (iny + 1) * (inz + 1), b, lb, s, ls);

		if (print_log_message) {

#if doubleintprecision == 1
			printf("maxbound=%d\n", t.maxbound);
#else
			printf("maxbound=%d\n", t.maxbound);
#endif
		}

		// constr_neighbors_for_the_internal_node_alice закодирована 14 сентября 2016. в 12_36.
		// заполнение prop_b. Замена constr_prop_bound. закодирована 14 сентября 2016. в 14_24.
		// nvtxcell остаётся тождественный nullptr. Визуализация в tecplot на адаптивной локально-измельчённой сетке
		// будет производится как описано в документации ANES решателя.
		// Заполнение t.border_neighbor, t.binternalsource. Замена constr_border_neighbor_temp закодирована 20 сентября 2016. в 11_24.
		constr_neighbors_for_the_internal_node_prop_b_alice(oc_global, t.border_neighbor, t.binternalsource, t.neighbors_for_the_internal_node,
			t.prop, t.prop_b, t.maxbound, t.maxelm, b, lb, s, ls, w, lw, t.whot_is_block, t.nvtx, t.pa);
	}
	t.maxp = t.maxelm + t.maxbound;


	if (print_log_message) {

#if doubleintprecision == 1
		printf("maxelm=%d maxbound=%d maxp=%d\n", t.maxelm, t.maxbound, t.maxp);
#else
		printf("maxelm=%d maxbound=%d maxp=%d\n", t.maxelm, t.maxbound, t.maxp);
#endif

		//printf("part 7\n");
#if doubleintprecision == 1
		printf("part %lld enumerate_gran_temp\n", icount_part++);
#else
		printf("part %d enumerate_gran_temp\n", icount_part++);
#endif

	}

	if (!bALICEflag) {
		// Вычисление соседей для каждого внутреннего узла. 
		// Причём множество соседей выбирается и среди  
		// внутренних КО и среди граничных КО.
		if (!b_first_constr_neighbors_for_the_internal_node) {
			if (t.neighbors_for_the_internal_node != nullptr) {

				if (print_log_message) {

					std::cout << "free t.neighbors_for_the_internal_node\n";
				}

				for (int i = 0; i < 12; ++i) {
					if (t.neighbors_for_the_internal_node[i] != nullptr) {
						for (int j = 0; j < 4; ++j) {
							if (t.neighbors_for_the_internal_node[i][j] != nullptr) {
								delete[] t.neighbors_for_the_internal_node[i][j];
								t.neighbors_for_the_internal_node[i][j] = nullptr;
							}
						}
					}
				}

				for (int i = 0; i < 12; ++i) {
					if (t.neighbors_for_the_internal_node[i] != nullptr) {
						delete[] t.neighbors_for_the_internal_node[i]; // -12N
						t.neighbors_for_the_internal_node[i] = nullptr;
					}
				}
			}
			if (t.neighbors_for_the_internal_node != nullptr) {
				delete[] t.neighbors_for_the_internal_node;
				t.neighbors_for_the_internal_node = nullptr;
			}
		}
		b_first_constr_neighbors_for_the_internal_node = false;
		constr_neighbors_for_the_internal_node(t.neighbors_for_the_internal_node, t.maxelm, gran_t, neighbour);
	}

	if (print_log_message) {

		//printf("part 8\n");
#if doubleintprecision == 1
		printf("part %lld constr_neighbors_for_the_internal_node\n", icount_part++); //8
#else
		printf("part %d constr_neighbors_for_the_internal_node\n", icount_part++); //8
#endif
	}

	// Заполнение информации о граничных узлах:
	if (!bALICEflag) {



		TOCHKA pavg;
		//pavg.x = 0.0;
		//pavg.y = 0.0;
		//pavg.z = 0.0;
		/*
		for (integer i_53 = 0; i_53 <= inx; i_53++) {
			pavg.x += xpos[i_53];
		}
		pavg.x /= (inx + 1.0);
		for (integer i_53 = 0; i_53 <= iny; i_53++) {
			pavg.y += ypos[i_53];
		}
		pavg.y /= (iny + 1.0);
		for (integer i_53 = 0; i_53 <= inz; i_53++) {
			pavg.z += zpos[i_53];
		}
		pavg.z /= (inz + 1.0);
		*/
		pavg.x = xpos[(inx + 1) / 2];
		pavg.y = ypos[(iny + 1) / 2];
		pavg.z = zpos[(inz + 1) / 2];
		// Внимание !!! не заполнено поле dS.
		// pavg - середина расчётной области. 
		// Требуется для ускорения вычислений.
		constr_border_neighbor_temp(t.border_neighbor, t.whot_is_block, t.binternalsource, t.maxelm, t.maxbound, gran_t, neighbour, t.neighbors_for_the_internal_node, t.nvtx, t.pa, b, lb, lw, w, s, ls, pavg);


	}


	if (print_log_message) {

		//printf("part 9\n");
#if doubleintprecision == 1
		printf("part %lld constr_border_neighbor_temp\n", icount_part++); //9
#else
		printf("part %d constr_border_neighbor_temp\n", icount_part++); //9
#endif
	}

	// Свойства материала на границе твердотельной области.
	if (!bALICEflag) {


		/*
		#ifdef _OPENMP
				int inumcor = number_cores();
				omp_set_num_threads(inumcor); // установка числа потоков
				//omp_set_num_threads(4);

		#pragma omp parallel
				{
		#pragma omp parallel sections
					{
		#pragma omp section
						{
							constr_prop_bound(t.prop, t.prop_b, t.maxelm, t.maxbound, gran_t, neighbour, t.nvtx, t.pa, b, lb);
						}
		#pragma omp section
						{
							// Нужно не забыть выделить память под potent.
							//t.alpha=0.9; // параметр нижней релаксации для уравнения теплопроводности
							// На прошлой неделе 12-16 сентября 2016, было выяснено что для стабильности сходимости алгоритма
							// должно стоять значение именно 1.0, а нижняя релаксация обеспечивается с помощью константы
							// bHORF == 0.25 как рекомендовано в Theory Guide ANSYS Fluent.
							t.alpha = 1.0; // ВНИМАНИЕ !!! только значение t.alpha = 1.0 стабильно.

							// выделение оперативной памяти для задачи теплопроводности.
							// Инициализация также производится. Производится автоматическое определение минимальной температуры на существующих стенках WALL.
							// Данная функция универсальна  подходит также и для AЛИС сетки.
							if (!breconstruct) {
								allocation_memory_temp(t.potent, t.total_deformation, t.slau, t.slau_bon, t.border_neighbor, t.maxelm, t.maxbound, ls, lw, w, temp_ref);
								//printf("part 11\n");
							}
						}
		#pragma omp section
						{
							//if (!bALICEflag) {
								constr_link_on_surface_for_radiation_model(t.maxelm, t.whot_is_block, t.neighbors_for_the_internal_node, t.nvtx, t.pa, b, lb);
							//}
						}
		#pragma omp section
						{
							//if (!bALICEflag) {
								if (bextendedprint) {
									printf("extended print. please wait...\n");
								}
								constr_nvtxcell(evt_t, t.border_neighbor, t.maxbound, t.maxelm, bextendedprint, t.nvtxcell, t.ncell, inx, iny, inz, tck_int_list);
							//}
						}


					}
				}

				//omp_set_num_threads(1);

		#else
		*/

		constr_prop_bound(t.prop, t.prop_b, t.maxelm, t.maxbound, gran_t, neighbour, t.nvtx, t.pa, b, lb);

		// Нужно не забыть выделить память под potent.
		//t.alpha=0.9; // параметр нижней релаксации для уравнения теплопроводности
		// На прошлой неделе 12-16 сентября 2016, было выяснено что для стабильности сходимости алгоритма
		// должно стоять значение именно 1.0, а нижняя релаксация обеспечивается с помощью константы 
		// bHORF == 0.25 как рекомендовано в Theory Guide ANSYS Fluent. 
		t.alpha = 1.0; // ВНИМАНИЕ !!! только значение t.alpha = 1.0 стабильно. 

		// выделение оперативной памяти для задачи теплопроводности.
		// Инициализация также производится. Производится автоматическое определение минимальной температуры на существующих стенках WALL.
		// Данная функция универсальна  подходит также и для AЛИС сетки.
		if (!breconstruct) {
			allocation_memory_temp(t.potent, t.total_deformation, t.slau, t.slau_bon, t.border_neighbor, t.maxelm, t.maxbound, ls, lw, w, temp_ref);
			//printf("part 11\n");
		}

		//if (!bALICEflag) {
		constr_link_on_surface_for_radiation_model(t.maxelm, t.whot_is_block, t.neighbors_for_the_internal_node, t.nvtx, t.pa, b, lb);

		if (print_log_message) {

			if (bextendedprint) {
				printf("extended print. please wait...\n");
			}
		}
		constr_nvtxcell(evt_t, t.border_neighbor, t.maxbound, t.maxelm, bextendedprint, t.nvtxcell, t.ncell, inx, iny, inz, tck_int_list);
		//}
//#endif
	}

	if (print_log_message) {

		//printf("part 10\n");
#if doubleintprecision == 1
		printf("part %lld constr_prop_bound\n", icount_part++); //10
#else
		printf("part %d constr_prop_bound\n", icount_part++); //10
#endif
	}
	//exit(1);


	// выделение оперативной памяти для задачи теплопроводности.
	// Инициализация также производится. Производится автоматическое определение минимальной температуры на существующих стенках WALL.
	// Данная функция универсальна  подходит также и для AЛИС сетки.
	if (!breconstruct) {
		if (print_log_message) {

			//printf("part 11\n");
#if doubleintprecision == 1
			printf("part %lld allocation_memory_temp\n", icount_part++); //11
#else
			printf("part %d allocation_memory_temp\n", icount_part++); //11
#endif
		}
	}




	// Создаёт необходимые связи для модели излучения внутри блока.
	if ((!breconstruct) && (bALICEflag)) {

		// Нужно не забыть выделить память под potent.
		//t.alpha=0.9; // параметр нижней релаксации для уравнения теплопроводности
		// На прошлой неделе 12-16 сентября 2016, было выяснено что для стабильности сходимости алгоритма
		// должно стоять значение именно 1.0, а нижняя релаксация обеспечивается с помощью константы 
		// bHORF == 0.25 как рекомендовано в Theory Guide ANSYS Fluent. 
		t.alpha = 1.0; // ВНИМАНИЕ !!! только значение t.alpha = 1.0 стабильно. 

		// выделение оперативной памяти для задачи теплопроводности.
		// Инициализация также производится. Производится автоматическое определение минимальной температуры на существующих стенках WALL.
		// Данная функция универсальна  подходит также и для AЛИС сетки.

		allocation_memory_temp(t.potent, t.total_deformation, t.slau, t.slau_bon, t.border_neighbor, t.maxelm, t.maxbound, ls, lw, w, temp_ref);
		//printf("part 11\n");


		// border_neighbor - введено для хранения площади грани граничного узла на АЛИС сетке.
		constr_link_on_surface_for_radiation_model_alice(t.maxelm, t.border_neighbor, t.whot_is_block, t.neighbors_for_the_internal_node, t.nvtx, t.pa, b, lb);
	}
	// создаёт связи между контрольными объёмами для графической 
	// визуализации.

	if (print_log_message) {

		//printf("part 11\n");
#if doubleintprecision == 1
		printf("part %lld constr_link_on_surface_for_radiation_model\n", icount_part++); //12
#else
		printf("part %d constr_link_on_surface_for_radiation_model\n", icount_part++); //12
#endif
	}



	if (!bALICEflag) {
		if (bextendedprint) {
			if (print_log_message) {

				printf("extended print. please wait...\n");
			}
		}
		//constr_nvtxcell(evt_t, t.border_neighbor, t.maxbound, t.maxelm, bextendedprint, t.nvtxcell, t.ncell, inx, iny, inz, tck_int_list);
	}

	//if (tck_int_list != nullptr) {
	delete[] tck_int_list;
	tck_int_list = nullptr;
	//}

#if doubleintprecision == 1
	/*
	for (integer i = 0; i <= inx; ++i) {
		printf("xpos[%lld]=%e\n", i, xpos[i]);
	}
	// индексация начинается с нуля и заканчивается значением inx
	for (integer i = 0; i <= iny; ++i) {
		printf("ypos[%lld]=%e\n", i, ypos[i]);
	}
	// индексация начинается с нуля и заканчивается значением inx
	for (integer i = 0; i <= inz; ++i) {
		printf("zpos[%lld]=%e\n", i, zpos[i]);
	}
	system("PAUSE");
	*/
#else
	/*
	for (integer i = 0; i <= inx; ++i) {
		printf("xpos[%d]=%e\n", i, xpos[i]);
	}
	// индексация начинается с нуля и заканчивается значением inx
	for (integer i = 0; i <= iny; ++i) {
		printf("ypos[%d]=%e\n", i, ypos[i]);
	}
	// индексация начинается с нуля и заканчивается значением inx
	for (integer i = 0; i <= inz; ++i) {
		printf("zpos[%d]=%e\n", i, zpos[i]);
	}
	system("PAUSE");
	*/
#endif

	if (print_log_message) {

		//printf("part 12\n");
#if doubleintprecision == 1
		printf("part %lld constr_nvtxcell\n", icount_part++); //12
#else
		printf("part %d constr_nvtxcell\n", icount_part++); //12
#endif
	}

	// При использовании АЛИС сетки визуализация производится в соответствии с документацией 
	// ANES солвера, поэтому добавление сеточных линий совсем не требуется.
	if ((!bALICEflag)) {

		// 15 августа. В случае нехватки узлов для качественной визуализации на основе 
		// Coarce Mesh сетки произойдёт автоматическое добавление недостающих линий и возврат.
		// Проверка nvtxcell
		integer iPnodes_count = 0;
		bool* found_hash_table = new bool[t.maxelm + 1];
		for (integer i = 1; i <= t.maxelm; ++i) {

			if (bfirst_additional_count_CAD) {

				found_hash_table[i] = false; // инициализация.
				integer ib = t.whot_is_block[i - 1];

				if (b[ib].g.itypegeom == POLYGON) {
					// 30.08.2017. на полигонах мы не используем этот метод.
					found_hash_table[i] = true; // инициализация.
				}
				if (b[ib].g.itypegeom == CAD_STL) {

					// 14.11.2020. на CAD STL мы не используем этот метод.
					// Для улучшения качества Coarce сетки используем. 21.11.2020
					if (!b[ib].g.CAD_is_PRISM()) {
						// CAD объект не порожден прямой прямоугольной призмой параллельной координатным осям.
						found_hash_table[i] = true; // инициализация.
					}
				}
			}
			else {
				// Не добавляем сеточные линии.
				found_hash_table[i] = true; // инициализация.
			}
		}

		bfirst_additional_count_CAD = false;

		for (integer i_1 = 0; i_1 < t.ncell; ++i_1) {
			found_hash_table[t.nvtxcell[0][i_1]] = true;
			found_hash_table[t.nvtxcell[1][i_1]] = true;
			found_hash_table[t.nvtxcell[2][i_1]] = true;
			found_hash_table[t.nvtxcell[3][i_1]] = true;
			found_hash_table[t.nvtxcell[4][i_1]] = true;
			found_hash_table[t.nvtxcell[5][i_1]] = true;
			found_hash_table[t.nvtxcell[6][i_1]] = true;
			found_hash_table[t.nvtxcell[7][i_1]] = true;
		}


		iPnodes_count = 0;
		for (integer i = 1; i <= t.maxelm; ++i) {
			if (!found_hash_table[i]) {
				// узел i не найден в nvtxcell
				iPnodes_count++;
			}
		}

		if (print_log_message) {

#if doubleintprecision == 1
			printf("problem visualisation nodes count iPnodes_count=%lld\n", iPnodes_count);
#else
			printf("problem visualisation nodes count iPnodes_count=%d\n", iPnodes_count);
#endif
		}
		/*
		bool bpolygon = false;
		for (integer i83 = 1; i83 < lb; i83++) {
			if (b[i83].g.itypegeom == POLYGON) {
				bpolygon = true;
			}
			if (b[i83].g.itypegeom == CAD_STL) {
				bpolygon = true;
			}
		}*/

		//system("PAUSE");
		/*
		bool bcondition = true;
		//if (ipolygon_count_limit > 1)
		{
			bcondition = bglobal_restart_06_10_2018;
		}*/
		if (/*(bcondition)&&*/(iPnodes_count > 0) && (((ireconstruction_free_construct_alloc == 1) && (ipolygon_count_limit < iGLOBAL_RESTART_LIMIT)) ||
			((ireconstruction_free_construct_alloc == 0) && (ipolygon_count_limit < 2))/*&&(!bpolygon)*/)) {

			for (integer i = 1; i <= t.maxelm; ++i) {
				if (!found_hash_table[i]) {
					integer iP = i - 1;
					// iP - номер центрального контрольного объёма
					// iP внутренний КО 0..maxelm-1
					integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
					iE = t.neighbors_for_the_internal_node[E_SIDE][0][iP];
					iN = t.neighbors_for_the_internal_node[N_SIDE][0][iP];
					iT = t.neighbors_for_the_internal_node[T_SIDE][0][iP];
					iW = t.neighbors_for_the_internal_node[W_SIDE][0][iP];
					iS = t.neighbors_for_the_internal_node[S_SIDE][0][iP];
					iB = t.neighbors_for_the_internal_node[B_SIDE][0][iP];

					// вычисление размеров текущего контрольного объёма:
					//doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контрольного объёма
					//volume3D(iP, t.nvtx, t.pa, dx, dy, dz);
					TOCHKA pointP;
					center_cord3D(iP, t.nvtx, t.pa, pointP, 100);


					if ((iE >= t.maxelm) && (iW >= t.maxelm)) {
						// добавить в xpos
						addboundary(xpos, inx, pointP.x, YZ_PLANE, b, lb, w, lw, s, ls);
					}
					if ((iN >= t.maxelm) && (iS >= t.maxelm)) {
						// добавить в ypos
						addboundary(ypos, iny, pointP.y, XZ_PLANE, b, lb, w, lw, s, ls);
					}
					//if (!b_one_cell_z)
					{
						// Только если у нас невырожденная задача по оси Oz.
						// код не готов к одной клетке.
						if ((iT >= t.maxelm) && (iB >= t.maxelm)) {
							// добавить в zpos
							addboundary(zpos, inz, pointP.z, XY_PLANE, b, lb, w, lw, s, ls);
						}
					}

				}
			}

			// упорядочивание по возрастанию
			Sort_method<doublereal>(xpos, inx);
			Sort_method<doublereal>(ypos, iny);
			Sort_method<doublereal>(zpos, inz);

			// Нужно освободить оперативную память из под всех структур данных:
			free_level1_temp(t);
			if (!breconstruct) {
				free_level2_temp(t); // освобождение памяти из под матриц.
			}

			delete[] found_hash_table;

			// 22.01.2017
			if (neighbour != nullptr) {
				for (integer i_97 = 0; i_97 < 12; i_97++) {
					if (neighbour[i_97] != nullptr) {
						delete[] neighbour[i_97];
						neighbour[i_97] = nullptr;
					}
				}
				delete[] neighbour;
			}
			neighbour = nullptr;

			if (gran_t != nullptr) {
				for (integer i_97 = 0; i_97 < 6; i_97++) {
					if (gran_t[i_97] != nullptr) {
						delete[] gran_t[i_97];
						gran_t[i_97] = nullptr;
					}
				}
				delete[] gran_t;
			}
			gran_t = nullptr;

			//if (ent_t != nullptr) {
			delete[] ent_t;
			//}
			ent_t = nullptr;

			//if (evt_t != nullptr) {
			delete[] evt_t;
			evt_t = nullptr;
			//}

			iPnodes_count_shadow_memo = iPnodes_count;

			if (print_log_message) {

				std::cout << "progon number is ipolygon_count_limit=" << ipolygon_count_limit << std::endl;
				std::cout << "iPnodes_count = " << iPnodes_count_shadow << "=" << iPnodes_count << ", x= " << inx_shadow << "=" << inx << ", y= " << iny_shadow << "=" << iny << ", z= " << inz_shadow << "=" << inz << "\n";
			}

			if (abs(iPnodes_count_shadow - iPnodes_count) +
				abs(inx_shadow - inx) + abs(iny_shadow - iny) + abs(inz_shadow - inz) == 0) {
				bglobal_restart_06_10_2018_stagnation[ipolygon_count_limit] = true;
				if (print_log_message) {

					printf("Stagnation situation...\n");
					//system("pause");
					//bglobal_restart_06_10_2018 = false;
					printf("bglobal_restart_06_10_2018 = false\n");
				}
			}
			else {
				bglobal_restart_06_10_2018_stagnation[ipolygon_count_limit] = false;
				//bglobal_restart_06_10_2018 = true;
				if (print_log_message) {

					printf("bglobal_restart_06_10_2018 = true;\n");
				}
			}
			//system("pause");
			/*
			if (ipolygon_count_limit == iGLOBAL_RESTART_LIMIT) {
				bool bcheck = true;
				for (integer i_4 = 1; i_4 <= iGLOBAL_RESTART_LIMIT; i_4++) {
					bcheck = bcheck && bglobal_restart_06_10_2018_stagnation[i_4];
				}
				if (bcheck) {
					bglobal_restart_06_10_2018 = false;
					printf("bglobal_restart_06_10_2018 = false\n");
				}
				else {
					bglobal_restart_06_10_2018 = true;
					printf("bglobal_restart_06_10_2018 = true;\n");
					//bglobal_restart_06_10_2018 = false;
				}
			}
			*/
			iPnodes_count_shadow = iPnodes_count;
			// горячий рестарт рестарт.
			goto RESTART_CONSTRUCT;

		}
		delete[] found_hash_table;
	}




	// ГИДРОДИНАМИКА:
	integer i = 0;

	int maxelm_global_flow = 0; // максимальное суммарное количество внутренних КО принадлежащих расчётной области.
	int* evt_f = nullptr; // глобальная нумерация контрольных объёмов
	int* whot_is_block_fl = nullptr; // побочная информация.

	if (print_log_message) {

		//printf("part 13\n");
#if doubleintprecision == 1
		printf("part %lld FLUID DYNAMIC start construct\n", icount_part++); //13
#else
		printf("part %d FLUID DYNAMIC start construct\n", icount_part++); //13
#endif
	}


	TOCKA_SHORT_INT* tck_int_list_flow = nullptr;

	if (!bALICEflag) {
		enumerate_volume(evt_f, maxelm_global_flow, HYDRODINAMIC, xpos, ypos, zpos, whot_is_block_fl, inx, iny, inz, b, lb, lu, my_union, iunion_id_p1, tck_int_list_flow, w, lw, s, ls);
	}
	else {
		calculate_max_elm(oc_global, maxelm_global_flow, HYDRODINAMIC, b, lb, false);
		if (print_log_message) {

			printf("maxelm calculate succseful.\n");
		}
		// Это потребуется для алгоритма заливки.
		init_evt_f_alice(evt_f, HYDRODINAMIC, xpos, ypos, zpos, inx, iny, inz, b, lb, tck_int_list_flow, w, lw, s, ls);
	}

	//if (whot_is_block_fl != nullptr) {
		// информация найдена преждевременно и мы её уничтожим.
	delete[] whot_is_block_fl;
	//}

	if (print_log_message) {

#if doubleintprecision == 1
		printf("maxelm_global_flow=%d\n", maxelm_global_flow);
#else
		printf("maxelm_global_flow=%d\n", maxelm_global_flow);
#endif

		//printf("part 14\n");
#if doubleintprecision == 1
		printf("part %lld enumerate_volume flow\n", icount_part++); //22
#else
		printf("part %d enumerate_volume flow\n", icount_part++); //22
#endif
}
	t.rootWE=nullptr; t.rootSN=nullptr; t.rootBT=nullptr;

	
	// 22 сентября 2016 была выделена первая часть, она является
	// универсальной и подходит для АЛИС сетки тоже.
	// Части 11, 12, 12.5, 13.
	// Заполнение evt_f, domain_id, evt_f2.
	int** evt_f2 = nullptr;
	integer* domain_id = nullptr;
	constr_ptr_temp_part1(flow_interior, evt_f, evt_f2, domain_id, inx, iny, inz, icount_part);

	

	   // Установка связей между
	   // гидродинамикой и теплообменом.
       // Части   14, 15, 16, 17, 18
	   // После 22 сентября 2016 эта функция несёт нагрузку второй части,
	   // первая часть универсальна и подходит и для АЛИС сетки.
	   if (!bALICEflag) {
		   if (domain_id != nullptr) {
			   delete[] domain_id;
			   domain_id = nullptr;
		   }
		   constr_ptr_temp_part2(flow_interior, evt_f, evt_f2, domain_id, inx, iny, inz, tck_int_list_flow, maxelm_global_flow, icount_part);

		   constr_ptr_temp(flow_interior, f, t.maxelm, t.ptr, evt_t, evt_f, evt_f2, domain_id, inx, iny, inz, breconstruct, icount_part);

		   if (evt_f != nullptr) {
			   delete[] evt_f;
		   }
		   evt_f = nullptr;
	   }
	   else {
		   // Построение структуры ptr на АЛИС сетке.
		   // Выделение памяти под FLOW.
		   // Дата создания 22 сентября 2016.
		   // Данная функция лишь несёт нагрузку выделения оперативной памяти, а так из неё всё выгружено при
		   // АЛИС сетке.
		   constr_ptr_temp_allocation_memory_alice(flow_interior, f, icount_part);
	   }
		
	   // Гарантированно освобождаем память для domain_id.
	   // для структурированной сетки память освобождена внутри функции constr_ptr_temp.
	   delete[] domain_id;
	   domain_id = nullptr;
		
		//STUB 17.11.2019
	   // Экспорт в графический визуализатор.
	   if (!bALICEflag) {
		   if (print_log_message) {

			   printf("please wait...\n");
		   }
		   // записывает первую и третью части экспортируемого в tecplot360
		   // файла на диск. Это делается для того чтобы выгрузить в дальнейшем 
		   // 8*t.maxelm*sizeof(doublereal) из оперативной памяти компьютера.
		   // Эти же файлы можно будет использовать и для печати гидродинамики.
		   // Т.к. печать планируется производить только для внутренних КО.
		   if (iunion_id_p1 == 0) {
			   exporttecplotxy360T_3D_part1and3(t, t.maxelm, t.maxbound, bextendedprint, t.ncell, t.nvtx, t.nvtxcell, t.pa, t.border_neighbor, 1, t.ptr);
		   }
		   else {
			   exporttecplotxy360T_3D_part1and3(my_union[iunion_id_p1 - 1].t, my_union[iunion_id_p1 - 1].t.maxelm,
				   my_union[iunion_id_p1 - 1].t.maxbound, bextendedprint, my_union[iunion_id_p1 - 1].t.ncell,
				   my_union[iunion_id_p1 - 1].t.nvtx, my_union[iunion_id_p1 - 1].t.nvtxcell,
				   my_union[iunion_id_p1 - 1].t.pa, my_union[iunion_id_p1 - 1].t.border_neighbor, 1, my_union[iunion_id_p1 - 1].t.ptr);

			   
		   }
	   }
	   if (t.nvtxcell != nullptr) {
		   for (integer i74 = 0; i74 < 8; i74++) { // -8N
			   if (t.nvtxcell[i74] != nullptr) {
				   delete[] t.nvtxcell[i74];
				   t.nvtxcell[i74] = nullptr;
			   }
		   }
		   delete[] t.nvtxcell;
	   }
	   t.nvtxcell = nullptr; // если указатель равен nullptr то это предохраняет от повторного удаления памяти.


		if (neighbour != nullptr) {
			for (i = 0; i < 12; ++i) {
				if (neighbour[i] != nullptr) {
					delete[] neighbour[i];
					neighbour[i] = nullptr;
				}
			}
			delete[] neighbour;
		}
		neighbour = nullptr;

		// освобождение оперативной памяти TEMPER
		if (ent_t != nullptr) {
			delete[] ent_t;
		}
		ent_t = nullptr;
		//if (evt_t != nullptr) {
			//delete evt_t;
		//}
		//evt_t = nullptr;
		if (neighbour != nullptr) {
			for (i = 0; i < 12; ++i) {
				if (neighbour[i] != nullptr) {
					delete[] neighbour[i];
					neighbour[i] = nullptr;
				}
			}
			delete[] neighbour;
		}
		neighbour = nullptr;

		if (gran_t != nullptr) {
			for (i = 0; i < 6; ++i) {
				if (gran_t[i] != nullptr) {
					delete[] gran_t[i];
					gran_t[i] = nullptr;
				}
			}
			delete[] gran_t;
		}
		gran_t = nullptr;

		if (maxelm_global_flow > 0) {
			// есть зоны FLUID

		// освобождение оперативной памяти TEMPER
			if (ent_t != nullptr) {
				delete[] ent_t;
			}
			ent_t = nullptr;

			if (neighbour != nullptr) {
				for (i = 0; i < 12; ++i) {
					if (neighbour[i] != nullptr) {
						delete[] neighbour[i];
						neighbour[i] = nullptr;
					}
				}
				delete[] neighbour;
			}
			neighbour = nullptr;

			if (gran_t != nullptr) {
				for (i = 0; i < 6; ++i) {
					if (gran_t[i] != nullptr) {
						delete[] gran_t[i];
						gran_t[i] = nullptr;
					}
				}
				delete[] gran_t;
			}
			gran_t = nullptr;

			integer j = 0; // счётчики
			int* ent_f = nullptr; // глобальная нумерация узлов.
			int** neighbour_f = nullptr; // соседи для внутренних КО только среди внутренних КО.
			int** gran_f = nullptr; // уникальная нумерация граней



			// Проход по всем жидким зонам:
			for (i = 0; i < flow_interior; ++i) {
				// своё построение для каждой зоны FLUID.


				if (bALICEflag) {
					// Замена частей 19, 20, 21 и 22.
					// Построение f[i].pa, f[i].maxnod. Замена constr_nodes_flow.
					// ВНИМАНИЕ !!! icolor_different_fluid_domain НЕ построено.
					// Построение f[i].nvtx, f[i].maxelm. Замена constr_nvtx_flow.
					// Построение f[i].prop. Замена constr_prop_flow.
					// 21.09.2016 Начало написания. Все сделано кроме icolor_different_fluid_domain.
					// Специально введено inum_FD в octree дерево.
					// 22.09.2016 Сделано заполнение f[0].ptr[fluid_elm_id]=соответствующий temper_elm_id.
					// 25.09.2016 Заполнение t.ptr. 
					constr_nodes_nvtx_prop_flow_alice(oc_global, inx, iny, inz, f[i].maxelm, xpos, ypos, zpos,
						HYDRODINAMIC, b, lb, f[i].pa, f[i].maxnod, f[i].nvtx, f[i].prop, matlist, f[i].ptr,
						f[i].whot_is_block, t.ptr, t.maxelm, t.bActiveShearModule);

					// 22 сентября 2016 вычисление функции цвета для АЛИС сетки.
					// Вычисление icolor_different_fluid_domain.
					constr_ptr_temp_part2(flow_interior, evt_f, evt_f2, domain_id, inx, iny, inz, tck_int_list_flow, maxelm_global_flow, icount_part);
					if (evt_f != nullptr) {
						delete[] evt_f;
					}
					evt_f = nullptr;
					// maxelm_global_flow в точности равно f[i].maxelm.
					// b, lb - передаются только для дебага.
					constr_icolor_different_fluid_domain_alice(f[i].maxelm, f[i].icolor_different_fluid_domain,
						inx, iny, inz, evt_f2, xpos, ypos, zpos, f[i].pa, f[i].nvtx, oc_global, tck_int_list_flow, b, lb);
				}




				// integer *ent_f; - глобальная нумерация узлов.
				if (!bALICEflag) {
					constr_nodes_flow(f[i].pa, f[i].maxnod, ent_f, evt_f2, i, inx, iny, inz, xpos, ypos, zpos);
				}


				//#endif

				if (print_log_message) {

#if doubleintprecision == 1
					printf("part %lld constr_nodes_flow\n", icount_part++); //19
#else
					printf("part %d constr_nodes_flow\n", icount_part++); //19
#endif
				}


				// для каждого контрольного объёма принадлежащему
				// расчётной области определяет номера его вершин.
				// icolor_different_fluid_domain, nvtx, 
				// вычисляет maxelm.
				if (!bALICEflag) {
					// Вычисляет maxelm;
					calculate_max_elm(f[i].maxelm, evt_f2,
						inx, iny, inz, i);
				}


#ifdef _OPENMP
				int inumcor = number_cores();
				omp_set_num_threads(inumcor); // установка числа потоков
#endif

/*
* // Не работает в параллельном режиме !!! отключил 24,07,2021
#ifdef _OPENMP
			int inumcor = number_cores();
			omp_set_num_threads(inumcor); // установка числа потоков
			//omp_set_num_threads(3);
#pragma omp parallel
			{
#pragma omp parallel sections
				{

#pragma omp section
					{
						if (!bALICEflag) {
							constr_nvtx_flow(evt_f2, f[i].icolor_different_fluid_domain,
								i, ent_f,
								f[i].nvtx, f[i].maxelm, inx, iny, inz);
						}
					}

#pragma omp section
					{
						// находит соседей для каждого внутреннего контрольного
									// объёма или 0 если соседа нет.
						if (!bALICEflag) {
							// Это просто пропускаем т.к. информация о neighbour_f
							//при сетках АЛИС не требуется, она хранится в octree дереве.
							constr_neighbour_flow(evt_f2, i, neighbour_f, f[i].maxelm, inx, iny, inz);
						}
					}
#pragma omp section
					{
						if (!bALICEflag) {
							constr_prop_flow(evt_t, t.whot_is_block, evt_f2,
								i, f[i].prop, f[i].maxelm,
								b, lb, inx, iny, inz, xpos, ypos, zpos, matlist);
						}
					}


				}
			}


			//omp_set_num_threads(1);

#else
*/

				if (!bALICEflag) {
					constr_nvtx_flow(evt_f2, f[i].icolor_different_fluid_domain, i, ent_f,
						f[i].nvtx, f[i].maxelm, inx, iny, inz);
				}

				// находит соседей для каждого внутреннего контрольного
				// объёма или 0 если соседа нет.
				if (!bALICEflag) {
					// Это просто пропускаем т.к. информация о neighbour_f
					//при сетках АЛИС не требуется, она хранится в octree дереве.
					constr_neighbour_flow(evt_f2, i, neighbour_f, f[i].maxelm, inx, iny, inz);
				}

				if (!bALICEflag) {
					constr_prop_flow(evt_t, t.whot_is_block, evt_f2, i, f[i].prop, f[i].maxelm, b, lb, inx, iny, inz, xpos, ypos, zpos, matlist);
				}

				//#endif

				if (print_log_message) {

#if doubleintprecision == 1
					printf("part %lld constr_nvtx_flow\n", icount_part++); //20
#else
					printf("part %d constr_nvtx_flow\n", icount_part++); //20
#endif
				}

				// использует neighbour_f
				if (!bALICEflag) {
					enumerate_gran_flow(gran_f, f[i].maxelm, f[i].nvtx, f[i].maxbound, neighbour_f, f[i].pa);
				}

				if (print_log_message) {

#if doubleintprecision == 1
					printf("part %lld constr_neighbour_flow\n", icount_part++); //21
#else
					printf("part %d constr_neighbour_flow\n", icount_part++); //21
#endif
				}

				// Создаёт связи между контрольными объёмами для графической 
				// визуализации. Своей графической визуализации у FLOW нет, она 
				// использует визуализацию предложенную в части TEMPER. Графическая
				// визуализация делается только для внутренних КО не затрагивая 
				// граничные узлы.
				/*
				* Внимание ! оставлен синтаксис старого вызова.
				* constr_nvtxcell(evt_f, f.nvtxcell, f.ncell, inx, iny, inz);
				* printf("part not_number\n");
				*/


				// Заносит свойства материалов внутренних КО в структуру prop.
				if (print_log_message) {

#if doubleintprecision == 1
					printf("%d\n", f[i].maxelm);
#else
					printf("%d\n", f[i].maxelm);
#endif
				}

				if (bwait) {
					//system("PAUSE");
					system("pause");
				}

				if (evt_t != nullptr) {
					delete[] evt_t;
				}
				evt_t = nullptr;
				if (print_log_message) {

#if doubleintprecision == 1
					printf("part %lld constr_properties_flow\n", icount_part++); //22
#else
					printf("part %d constr_properties_flow\n", icount_part++); //22
#endif
				}




				/*
				 *  neighbors_for_the_internal_node, border_neighbor ( i - internal, b - boundary).
				 */

				 // Алгоритм уникальной нумерации граней.
				 // В результате каждая граничная грань,
				 // принадлежащая расчётной области имеет
				 // уникальный номер.
				if (!bALICEflag) {

				}
				else {
					calculate_max_bound_flow(oc_global, f[i].maxbound, (inx + 1) * (iny + 1) * (inz + 1), b, lb, s, ls);
					if (print_log_message) {

#if doubleintprecision == 1
						printf("maxbound=%d\n", f[i].maxbound);
#else
						printf("maxbound=%d\n", f[i].maxbound);
#endif
					}

				}
				f[i].maxp = f[i].maxelm + f[i].maxbound;

				if (print_log_message) {

#if doubleintprecision == 1
					printf("maxelm=%d maxbound=%d maxp=%d\n", f[i].maxelm, f[i].maxbound, f[i].maxp);
					printf("part %lld enumerate_gran_flow\n", icount_part++);
#else
					printf("maxelm=%d maxbound=%d maxp=%d\n", f[i].maxelm, f[i].maxbound, f[i].maxp);
					printf("part %d enumerate_gran_flow\n", icount_part++);
#endif
				}





				if (bALICEflag) {
					// constr_neighbors_for_the_internal_node_prop_b_flow_alice закодирована  22 сентября 2016 на основе
					// constr_neighbors_for_the_internal_node_prop_b_alice которая была закодирована 14 сентября 2016. в 12_36.
					// заполнение prop_b. Замена constr_prop_bound_flow. закодирована 14 сентября 2016. в 14_24.
					// nvtxcell остаётся тождественный nullptr. Визуализация в tecplot на адаптивной локально-измельчённой сетке
					// будет производится как описано в документации ANES решателя.
					// Заполнение f[i].border_neighbor Замена constr_border_neighbor_flow закодирована 20 сентября 2016. в 11_24.
					// Работоспособна 24 сентября 2016.
					constr_neighbors_for_the_internal_node_prop_b_flow_alice(oc_global, f[i].border_neighbor,
						f[i].neighbors_for_the_internal_node, f[i].prop, f[i].prop_b,
						f[i].maxbound, f[i].maxelm,
						b, lb, s, ls, w, lw, f[i].whot_is_block, f[i].nvtx, f[i].pa, matlist);

					/*
					// 26.09.2016 Добавок для АЛИС сетки.
					// Убрано 26.03.2019 Как  бессмысленная проверка.
					integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
					integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
					integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов
					integer iP = 13; // !!!! внимание почему именно 13. Проверка бессмысленна.

					// -1 если не используется и [0..maxelm+maxbound-1] если используется.

					iE2 = f[i].neighbors_for_the_internal_node[ESIDE][iP].iNODE2; iN2 = f[i].neighbors_for_the_internal_node[NSIDE][iP].iNODE2; iT2 = f[i].neighbors_for_the_internal_node[TSIDE][iP].iNODE2;
					iW2 = f[i].neighbors_for_the_internal_node[WSIDE][iP].iNODE2; iS2 = f[i].neighbors_for_the_internal_node[SSIDE][iP].iNODE2; iB2 = f[i].neighbors_for_the_internal_node[BSIDE][iP].iNODE2;
					iE3 = f[i].neighbors_for_the_internal_node[ESIDE][iP].iNODE3; iN3 = f[i].neighbors_for_the_internal_node[NSIDE][iP].iNODE3; iT3 = f[i].neighbors_for_the_internal_node[TSIDE][iP].iNODE3;
					iW3 = f[i].neighbors_for_the_internal_node[WSIDE][iP].iNODE3; iS3 = f[i].neighbors_for_the_internal_node[SSIDE][iP].iNODE3; iB3 = f[i].neighbors_for_the_internal_node[BSIDE][iP].iNODE3;
					iE4 = f[i].neighbors_for_the_internal_node[ESIDE][iP].iNODE4; iN4 = f[i].neighbors_for_the_internal_node[NSIDE][iP].iNODE4; iT4 = f[i].neighbors_for_the_internal_node[TSIDE][iP].iNODE4;
					iW4 = f[i].neighbors_for_the_internal_node[WSIDE][iP].iNODE4; iS4 = f[i].neighbors_for_the_internal_node[SSIDE][iP].iNODE4; iB4 = f[i].neighbors_for_the_internal_node[BSIDE][iP].iNODE4;

					if (iE2 > -1 || iW2 > -1 || iN2 > -1 || iS2 > -1 || iT2 > -1 || iB2 > -1 ||
						iE3 > -1 || iW3 > -1 || iN3 > -1 || iS3 > -1 || iT3 > -1 || iB3 > -1 ||
						iE4 > -1 || iW4 > -1 || iN4 > -1 || iS4 > -1 || iT4 > -1 || iB4 > -1) {
						printf("correct\n");
						//system("PAUSE");
					}
					else {
						printf("incorrect\n");
						system("PAUSE");
					}
					*/

				}

				// Вычисление соседей для каждого внутреннего узла. 
				// Причём множество соседей выбирается и среди внутренних КО и 
				// среди граничных КО.
				if (!bALICEflag) {
					if (f[i].neighbors_for_the_internal_node != nullptr) {

						if (print_log_message) {

							std::cout << "free f[i].neighbors_for_the_internal_node\n";
						}

						for (int i_72 = 0; i_72 < 12; i_72++) {
							if (f[i].neighbors_for_the_internal_node[i_72] != nullptr) {
								for (int j_72 = 0; j_72 < 4; j_72++) {
									if (f[i].neighbors_for_the_internal_node[i_72][j_72] != nullptr) {
										delete[] f[i].neighbors_for_the_internal_node[i_72][j_72];
										f[i].neighbors_for_the_internal_node[i_72][j_72] = nullptr;
									}
								}
							}
						}

						for (int i_72 = 0; i_72 < 12; i_72++) {
							if (f[i].neighbors_for_the_internal_node[i_72] != nullptr) {
								delete[] f[i].neighbors_for_the_internal_node[i_72]; // -12N
								f[i].neighbors_for_the_internal_node[i_72] = nullptr;
							}
						}

						if (f[i].neighbors_for_the_internal_node != nullptr) {
							delete[] f[i].neighbors_for_the_internal_node;
							f[i].neighbors_for_the_internal_node = nullptr;
						}
					}

					constr_neighbors_for_the_internal_node(f[i].neighbors_for_the_internal_node, f[i].maxelm, gran_f, neighbour_f);
				}

				if (print_log_message) {

#if doubleintprecision == 1
					printf("part %lld constr_neighbour\n", icount_part++); //24
#else
					printf("part %d constr_neighbour\n", icount_part++); //24
#endif
				}


#ifdef _OPENMP 
				inumcor = number_cores();
				omp_set_num_threads(inumcor); // установка числа потоков
				//omp_set_num_threads(2); // установка числа потоков
#endif

/*#ifdef _OPENMP
			inumcor = number_cores();
			omp_set_num_threads(inumcor); // установка числа потоков
			//omp_set_num_threads(2); // установка числа потоков

#pragma omp parallel
			{
#pragma omp parallel sections
				{
#pragma omp section
					{
						// Заполнение информации о граничных узлах:
						if (!bALICEflag) {
							constr_border_neighbor_flow(f[i].border_neighbor, t.whot_is_block, f[i].ptr, f[i].maxelm, f[i].maxbound, gran_f, neighbour_f, f[i].neighbors_for_the_internal_node, f[i].nvtx, f[i].pa, lw, w, ls, b);
						}

						// 24 сентября 2016 (данная функция универсальна и подходит также и для АЛИС сетки).
						// Нужно ли фиксировать давление в одной точке.
						// По умолчанию используется последняя граничная точка.
						determination_of_activity_flow(f[i].border_neighbor, f[i].maxbound, ls, lw, w, f[i].bactive, f[i].bPressureFix, f[i].bLR1free);

					}
#pragma omp section
					{
						// Свойства жидкого материала на границе жидкой области.
						if (!bALICEflag) {
							constr_prop_bound_flow(f[i].prop, f[i].prop_b, f[i].maxelm, f[i].maxbound, gran_f, neighbour_f);
						}

					}
				}
			}

#ifdef _OPENMP
			omp_set_num_threads(1);
#endif

#else*/
// Заполнение информации о граничных узлах:
				if (!bALICEflag) {
					constr_border_neighbor_flow(f[i].border_neighbor, t.whot_is_block, f[i].ptr, f[i].maxelm, f[i].maxbound, gran_f, neighbour_f, f[i].neighbors_for_the_internal_node, f[i].nvtx, f[i].pa, lw, w, ls, b);
				}

				// Свойства жидкого материала на границе жидкой области.
				if (!bALICEflag) {
					constr_prop_bound_flow(f[i].prop, f[i].prop_b, f[i].maxelm, f[i].maxbound, gran_f, neighbour_f);
				}

				// 24 сентября 2016 (данная функция универсальна и подходит также и для АЛИС сетки).
				// Нужно ли фиксировать давление в одной точке.
				// По умолчанию используется последняя граничная точка.
				determination_of_activity_flow(f[i].border_neighbor, f[i].maxbound, ls, lw, w, f[i].bactive, f[i].bPressureFix, f[i].bLR1free);


				//#endif

				if (print_log_message) {

#if doubleintprecision == 1
					printf("part %lld constr_neighbour_db_flow\n", icount_part++); //25
#else
					printf("part %d constr_neighbour_db_flow\n", icount_part++); //25
#endif
				}



				if (print_log_message) {

#if doubleintprecision == 1
					printf("part %lld constr_prop_bound_flow\n", icount_part++); //26
#else
					printf("part %d constr_prop_bound_flow\n", icount_part++); //26
#endif



#if doubleintprecision == 1
					printf("part %lld determination_of_activity_flow\n", icount_part++); //27
#else
					printf("part %d determination_of_activity_flow\n", icount_part++); //27
#endif
				}





				// Нужно не забыть выделить память под potent.
				// Также необходимо установить параметры релаксации.

				// выделение оперативной памяти для задачи гидродинамики.
				// инициализирует правильным образом вектор potent.
				// загружает компоненты скорости и турбулентную вязкость из файла load.txt
				// 24 сентября 2016 (теперь этот метод работает и на АЛИС сетке).
				if (!breconstruct) {
					// 31.03.2019 Передаём внутрь также значения допусков
					// для функции is_EXTERNAL_FLOW
					allocation_memory_flow(f[i].potent, f[i].slau, f[i].slau_bon, f[i].border_neighbor,
						f[i].maxelm, f[i].maxbound, f[i].alpha, ls, lw, w, dgx, dgy, dgz, f[i].nvtx,
						f[i].pa, f[i].prop, f[i].neighbors_for_the_internal_node, fabs(xpos[1] - xpos[0]), fabs(ypos[1] - ypos[0]), fabs(zpos[1] - zpos[0]),
						fabs(xpos[inx] - xpos[inx - 1]), fabs(ypos[iny] - ypos[iny - 1]), fabs(zpos[inz] - zpos[inz - 1]),
						f[i].iflowregime, f[i].prop_b);
				}
				allocation_memory_flow_2(
					f[i].diag_coef, f[i].maxelm, f[i].maxbound,
					f[i].SInvariantStrainRateTensor,
					f[i].mf);

				f[i].center_coord = new TOCHKA[f[i].maxelm];
				f[i].volume = new TOCHKA[f[i].maxelm];
				for (integer i60 = 0; i60 < f[i].maxelm; i60++) {

					TOCHKA p_loc;
					center_cord3D(i60, f[i].nvtx, f[i].pa, p_loc, 100);

					f[i].center_coord[i60] = p_loc;

					doublereal dx = 0.0, dy = 0.0, dz = 0.0; // размеры контрольного объёма
					volume3D(i60, f[i].nvtx, f[i].pa, dx, dy, dz);

					f[i].volume[i60].x = fabs(dx);
					f[i].volume[i60].y = fabs(dy);
					f[i].volume[i60].z = fabs(dz);


				}

				{
					// Загрузка распределения начальной скорости.

					FILE* fp_inicialization_data;
#ifdef MINGW_COMPILLER
					int err_inicialization_data = 0;
					fp_inicialization_data = fopen64("load.txt", "r");
					if (fp_inicialization_data == NULL) {
						err_inicialization_data = 1;
				}
#else
					errno_t err_inicialization_data = 0;
					err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif
					if ((f[i].maxelm > 0) && ((err_inicialization_data == 0) || 
						(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_STEADY) ||
						(steady_or_unsteady_global_determinant == PHYSICAL_MODEL_SWITCH::CFD_UNSTEADY) ||
						((starting_speed_Vx * starting_speed_Vx +
							starting_speed_Vy * starting_speed_Vy +
							starting_speed_Vz * starting_speed_Vz > 1.0e-30) &&
							(steady_or_unsteady_global_determinant != PHYSICAL_MODEL_SWITCH::CFD_STEADY)))) {
						// Только если присутствуют жидкие ячейки.


						// открытие успешно и файл присутствует.
						// Либо мы считаем модуль "ПИОНЕР" аналитически заданная пользователем 
						// скорость не равна нулю и мы не считаем гидродинамику на основе
						// SIMPLE алгоритма.



						// Инициализация компонент скорости во внутренности расчётной области.
						// 26 марта 2017.
						for (integer i32 = 0; i32 < f[i].maxelm; i32++) {
							// 15.09.2018	


							// вычисляем скорректированный массовый поток через грани КО.
							// Массовый поток вычисляется по обычным формулам но в данном
							// случае без монотонизирующей поправки Рхи-Чоу. При его вычислении используются
							// простая линейная интерполяция скорости на грань КО.

							//bool bsimplelinearinterpol = true; // выполняется простая линейная интерполяция скорости на грань.


							return_calc_correct_mass_flux_only_interpolation(i32,
								f[i].potent,
								f[i].pa,
								f[i].prop,
								f[i].prop_b,
								f[i].nvtx,
								f[i].neighbors_for_the_internal_node,
								f[i].maxelm,
								f[i].mf[i32],
								f[i].border_neighbor,
								ls, lw,
								t.ilevel_alice, f[i].ptr);

							//if (fabs(f[i].mf[i32][BSIDE]) + fabs(f[i].mf[i32][TSIDE]) > 0.0) {
							//	printf("non zero mf. Ok.\n");
							//	system("PAUSE");
							//}

							//printf("incomming\n");
							//printf("%e %e %e %d\n", starting_speed_Vx, starting_speed_Vy, starting_speed_Vz, steady_or_unsteady_global_determinant);
							//system("PAUSE");

							//printf("%e\n", f[i].mf[i32][TSIDE]);
							//doublereal ts = f[i].mf[i32][TSIDE] + f[i].mf[i32][BSIDE] + f[i].mf[i32][ESIDE] + f[i].mf[i32][WSIDE] + f[i].mf[i32][NSIDE] + f[i].mf[i32][SSIDE];
							//if (ts != ts) {
								//printf("%d %e %e %e %e %e %e\n",i32, f[i].mf[i32][TSIDE], f[i].mf[i32][BSIDE], f[i].mf[i32][ESIDE], f[i].mf[i32][WSIDE], f[i].mf[i32][NSIDE], f[i].mf[i32][SSIDE]);
							//}

						}
						if (err_inicialization_data == 0) {
							// файл точно был успешно открыт до этого.
							fclose(fp_inicialization_data);
						}
						//system("PAUSE");
					}
			}

				if (print_log_message) {

#if doubleintprecision == 1
					printf("part %lld allocation_memory_flow_2.\n", icount_part++); //28
#else
					printf("part %d allocation_memory_flow_2\n", icount_part++); //28
#endif
		}
																		 
																		 
																		 
			//system("PAUSE");
			 if (0) {
		        xyplot( f, flow_interior, t);
		        printf("after allocate memory. OK.\n");
				system("pause");
	            //system("PAUSE"); // debug 
	         }


			// Освобождение оперативной памяти для каждой зоны FLUID.
			 if (ent_f != nullptr) {
				 delete[] ent_f;
			 }
			ent_f=nullptr;
			if (neighbour_f != nullptr) {
				for (j = 0; j < 12; ++j) {
					if (neighbour_f[j] != nullptr) {
						delete[] neighbour_f[j];
						neighbour_f[j] = nullptr;
					}
				}
				delete[] neighbour_f;
			}
			neighbour_f=nullptr;

			if (gran_f != nullptr) {
				for (j = 0; j < 6; ++j) {
					if (gran_f[j] != nullptr) {
						delete[] gran_f[j];
						gran_f[j] = nullptr;
					}
				}
				delete[] gran_f;
			}			
			gran_f=nullptr;

			// инициализация для полилинейного метода:
			// 

			f[i].OpTemp=temp_ref; // Operating Temperature

			if (print_log_message) {

#if doubleintprecision == 1
				printf("end FLUID %lld\n", i);
#else
				printf("end FLUID %d\n", i);
#endif
			}
			
			//system("PAUSE");

		} // конец цикла по всем FLUID зонам.


		// Вычисление расстояния до стенки не предназначено для 
		// АЛИС сетки т.к. там содержится сборка матрицы не учитывающая  
		// возможного наличия четырёх соседей.
			if (print_log_message) {

				printf("share information about the turbulence model\n");
			}
		constr_fluid_equation(f, flow_interior, ls, lw, w);
		if (print_log_message) {

			printf("end fluid\n");
		}
		
		if (evt_f2 != nullptr) {
			for (i = 0; i < 2; ++i) {
				if (evt_f2[i] != nullptr) {
					delete[] evt_f2[i];
					evt_f2[i] = nullptr;
				}
			}
			delete[] evt_f2;
		}
		evt_f2=nullptr;

	} else {
       
	   if (!breconstruct) {
		   if (f != nullptr) {
			   // Освобождение памяти для структуры flow.
			   free_level1_flow(f, flow_interior);
			   free_level2_flow(f, flow_interior); // освобождение памяти из под матриц.
			   delete[] f;
			   f = nullptr;
		   }
		   f = new FLOW[1];
		   f[0].ptr = nullptr;
	   }
	   flow_interior = 0; // нет зон FLUID
	   // Связи гидродинамики с теплообменом нет, т.к. 
	   // нет гидродинамической составляющей.
	   if (!breconstruct) {
		   // Если память занята то мы её освободим.
		   if (t.ptr != nullptr) {
			   for (i = 0; i<2; ++i) {
				   if (t.ptr[i] != nullptr) {
					   delete[] t.ptr[i];
					   t.ptr[i] = nullptr;
				   }
			   }
			   delete[] t.ptr;
			   t.ptr = nullptr;
		   }
		   t.ptr = nullptr;
	   }
	   // обнуляющая инициализация ссылок.
	   f[0].nvtx=nullptr;
	   f[0].maxnod=0;
	   f[0].maxelm=0;
	   f[0].pa=nullptr;
	   f[0].neighbors_for_the_internal_node=nullptr;
	   f[0].maxbound=0;
	   f[0].maxp=0;
	   f[0].border_neighbor=nullptr;
	   f[0].icolor_different_fluid_domain = nullptr;
	   if (f[0].ptr != nullptr) {
		   delete[] f[0].ptr;
	   }
	   f[0].ptr=nullptr;
	   f[0].potent=nullptr;
	   f[0].prop=nullptr;
	   f[0].prop_b=nullptr;
	   f[0].alpha=nullptr;
	   if (!breconstruct) {
		   f[0].slau = nullptr;
		   f[0].slau_bon = nullptr;
	   }
	   f[0].mf=nullptr;
	   f[0].rdistWall=nullptr;
	   f[0].SInvariantStrainRateTensor=nullptr;
	   f[0].diag_coef=nullptr;
	   f[0].bactive=false;
	   f[0].OpTemp=temp_ref;
	   f[0].bLR1free=true;
	   f[0].bPressureFix=true;
	   if (!breconstruct) {
		   f[0].iN = nullptr;
		   f[0].id = nullptr;
	   }
	   f[0].iflowregime= VISCOSITY_MODEL::LAMINAR;

	   // Вся выделенная память должна быть освобождена.
	   if (evt_f2 != nullptr) {
		   for (i = 0; i < 2; ++i) {
			   if (evt_f2[i] != nullptr) {
				   delete[] evt_f2[i];
				   evt_f2[i] = nullptr;
			   }
		   }
		   delete[] evt_f2;
	   }
	   evt_f2 = nullptr;
	}

	if (tck_int_list_flow != nullptr) {
		delete[] tck_int_list_flow;
		tck_int_list_flow = nullptr;
	}

	if (evt_t != nullptr) {
		delete[] evt_t;
	}
	evt_t = nullptr;
	
	if (maxelm_global_flow == 0) {
        // освобождение оперативной памяти TEMPER
		if (ent_t != nullptr) {
			delete[] ent_t;
		}
		ent_t=nullptr;
		if (evt_t != nullptr) {
			delete[] evt_t;
		}
		evt_t=nullptr;
		if (neighbour != nullptr) {
			for (i = 0; i < 6; ++i) {
				if (neighbour[i] != nullptr) {
					delete[] neighbour[i];
					neighbour[i] = nullptr;
				}
			}
			delete[] neighbour;
		}
		neighbour=nullptr;
		if (gran_t != nullptr) {
			for (i = 0; i < 6; ++i) {
				if (gran_t[i] != nullptr) {
					delete[] gran_t[i];
					gran_t[i] = nullptr;
				}
			}
			delete[] gran_t;
		}
		gran_t=nullptr;

		// освобождение оперативной памяти FLOW
		if (evt_f != nullptr) {
			delete[] evt_f;
		}
		evt_f=nullptr;
	}

	if (0) {
		// Старая версия xyplot на АЛИС сетке использовать невозможно и нужно прописать интерполяцию.
		xyplot( f, flow_interior, t);
		printf("load temper and flow. OK.\n");
		system("pause");
	    //system("PAUSE"); // debug avtosave
	}
    	
} // load_TEMPER_and_FLOW


/*
// записывает первую часть файла meshin.txt,
// а именно координаты узлов.
// Для задач теплопроводности.
void meshin_nodes_temp(char* fname, TEMPER t) {
   
    FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, fname, "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		integer i=0;
		// первая секция - координат узлов.
		for (i=0; i<t.maxnod; ++i) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld %e %e %e\n",i+1,t.pa[i].x,t.pa[i].y,t.pa[i].z);
		#else
			fprintf(fp,"%d %e %e %e\n",i+1,t.pa[i].x,t.pa[i].y,t.pa[i].z);
		#endif
			
		}
        fprintf(fp,"\n");
		// вторая секция - контрольные объёмы.
		for (i=0; i<t.maxelm; ++i) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld %lld %lld %lld %lld",i+1,t.nvtx[0][i],t.nvtx[1][i],t.nvtx[2][i],t.nvtx[3][i]);
			fprintf(fp," %lld %lld %lld %lld", t.nvtx[4][i],t.nvtx[5][i],t.nvtx[6][i],t.nvtx[7][i]);
			fprintf(fp," %lld %lld %lld", t.neighbors_for_the_internal_node[ESIDE][i].iNODE1, t.neighbors_for_the_internal_node[NSIDE][i].iNODE1, t.neighbors_for_the_internal_node[TSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", t.neighbors_for_the_internal_node[WSIDE][i].iNODE1, t.neighbors_for_the_internal_node[SSIDE][i].iNODE1, t.neighbors_for_the_internal_node[BSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", t.neighbors_for_the_internal_node[EE][i].iNODE1, t.neighbors_for_the_internal_node[NN][i].iNODE1, t.neighbors_for_the_internal_node[TTSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", t.neighbors_for_the_internal_node[WW][i].iNODE1, t.neighbors_for_the_internal_node[SS][i].iNODE1, t.neighbors_for_the_internal_node[BB][i].iNODE1);
		#else
			fprintf(fp,"%d %d %d %d %d",i+1,t.nvtx[0][i],t.nvtx[1][i],t.nvtx[2][i],t.nvtx[3][i]);
			fprintf(fp," %d %d %d %d", t.nvtx[4][i],t.nvtx[5][i],t.nvtx[6][i],t.nvtx[7][i]);
			fprintf(fp," %d %d %d", t.neighbors_for_the_internal_node[ESIDE][i].iNODE1, t.neighbors_for_the_internal_node[NSIDE][i].iNODE1, t.neighbors_for_the_internal_node[TSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", t.neighbors_for_the_internal_node[WSIDE][i].iNODE1, t.neighbors_for_the_internal_node[SSIDE][i].iNODE1, t.neighbors_for_the_internal_node[BSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", t.neighbors_for_the_internal_node[EE][i].iNODE1, t.neighbors_for_the_internal_node[NN][i].iNODE1, t.neighbors_for_the_internal_node[TTSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", t.neighbors_for_the_internal_node[WW][i].iNODE1, t.neighbors_for_the_internal_node[SS][i].iNODE1, t.neighbors_for_the_internal_node[BB][i].iNODE1);
		#endif
			
			fprintf(fp," %e %e %e", t.prop[RHO][i], t.prop[HEAT_CAPACITY][i], t.prop[LAM][i]);
			fprintf(fp," %e %e %e", t.prop[MULT_LAM_X][i], t.prop[MULT_LAM_Y][i], t.prop[MULT_LAM_Z][i]);
		}
        fprintf(fp,"\n");
        // четвёртая секция - связей контрольных объёмов для графической визуализации.
		for (i=0; i<t.ncell; ++i) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld %lld %lld %lld ",t.nvtxcell[0][i], t.nvtxcell[1][i], t.nvtxcell[2][i], t.nvtxcell[3][i]);
			fprintf(fp,"%lld %lld %lld %lld\n",t.nvtxcell[4][i], t.nvtxcell[5][i], t.nvtxcell[6][i], t.nvtxcell[7][i]);
		#else
			fprintf(fp,"%d %d %d %d ",t.nvtxcell[0][i], t.nvtxcell[1][i], t.nvtxcell[2][i], t.nvtxcell[3][i]);
			fprintf(fp,"%d %d %d %d\n",t.nvtxcell[4][i], t.nvtxcell[5][i], t.nvtxcell[6][i], t.nvtxcell[7][i]);
		#endif
			
		}
		fprintf(fp,"\n");
       fclose(fp); // закрытие файла
	}
} // meshin_nodes_temp
*/

/*
// записывает первую часть файла meshin.txt,
// а именно координаты узлов.
// Для задач гидродинамики.
void meshin_nodes_flow(char* fname, FLOW f) {
   
    FILE *fp;
	errno_t err;
	// создание файла для записи.
	if ((err = fopen_s( &fp, fname, "w")) != 0) {
		printf("Create File Error\n");
	}
	else {
		integer i=0;
        // первая секция - координат узлов.
		for (i=0; i<f.maxnod; ++i) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld %e %e %e\n",i+1,f.pa[i].x,f.pa[i].y,f.pa[i].z);
		#else
			fprintf(fp,"%d %e %e %e\n",i+1,f.pa[i].x,f.pa[i].y,f.pa[i].z);
		#endif
			
		}
		fprintf(fp,"\n");
		// вторая секция - контрольные объёмы.
		for (i=0; i<f.maxelm; ++i) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld %lld %lld %lld %lld",i+1,f.nvtx[0][i],f.nvtx[1][i],f.nvtx[2][i],f.nvtx[3][i]);
			fprintf(fp," %lld %lld %lld %lld",f.nvtx[4][i],f.nvtx[5][i],f.nvtx[6][i],f.nvtx[7][i]);
			fprintf(fp," %lld %lld %lld", f.neighbors_for_the_internal_node[ESIDE][i].iNODE1, f.neighbors_for_the_internal_node[NSIDE][i].iNODE1, f.neighbors_for_the_internal_node[TSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", f.neighbors_for_the_internal_node[WSIDE][i].iNODE1, f.neighbors_for_the_internal_node[SSIDE][i].iNODE1, f.neighbors_for_the_internal_node[BSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", f.neighbors_for_the_internal_node[EE][i].iNODE1, f.neighbors_for_the_internal_node[NN][i].iNODE1, f.neighbors_for_the_internal_node[TTSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", f.neighbors_for_the_internal_node[WW][i].iNODE1, f.neighbors_for_the_internal_node[SS][i].iNODE1, f.neighbors_for_the_internal_node[BB][i].iNODE1);
		#else
			fprintf(fp,"%d %d %d %d %d",i+1,f.nvtx[0][i],f.nvtx[1][i],f.nvtx[2][i],f.nvtx[3][i]);
			fprintf(fp," %d %d %d %d",f.nvtx[4][i],f.nvtx[5][i],f.nvtx[6][i],f.nvtx[7][i]);
			fprintf(fp," %d %d %d", f.neighbors_for_the_internal_node[ESIDE][i].iNODE1, f.neighbors_for_the_internal_node[NSIDE][i].iNODE1, f.neighbors_for_the_internal_node[TSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", f.neighbors_for_the_internal_node[WSIDE][i].iNODE1, f.neighbors_for_the_internal_node[SSIDE][i].iNODE1, f.neighbors_for_the_internal_node[BSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", f.neighbors_for_the_internal_node[EE][i].iNODE1, f.neighbors_for_the_internal_node[NN][i].iNODE1, f.neighbors_for_the_internal_node[TTSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", f.neighbors_for_the_internal_node[WW][i].iNODE1, f.neighbors_for_the_internal_node[SS][i].iNODE1, f.neighbors_for_the_internal_node[BB][i].iNODE1);
		#endif
			
			fprintf(fp," %e %e %e\n", f.prop[RHO][i], f.prop[MU][i], f.prop[BETA_T][i]);
		}
        fprintf(fp,"\n");
		
        fclose(fp); // закрытие файла
	}
} // meshin_nodes_flow
*/


#endif
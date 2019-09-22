// Файл constr_struct.cpp собирает структуры
// для решения задачи сопряжённого теплообмена.
// Данная версия работает на декартовой структурированной сетке с семиточечным шаблоном.
// Функционал данного модуля универсальный, он также предназначен также для вызова
// соответствующих функций из модуля constr_struct_alice.cpp на адаптивной локально измельченной сетке АЛИС.
// 26.09.2016 Написание данного модуля в основном завершено.

#pragma once
#ifndef _CONSTR_STRUCT_CPP_
#define _CONSTR_STRUCT_CPP_ 1

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
void CONSTRUCT_SECONDARY_LEVEL_OCTREE_LOAD1(integer identifikator,
	doublereal xc47, doublereal yc47, doublereal zc47, 
	integer*** &oct_load1, integer** &i_22,
	doublereal avgx[], doublereal avgy[], doublereal avgz[],
	doublereal* &x47, doublereal* &y47, doublereal* &z47,
	integer** &nvtx47, bool &bfound, integer &ifound) {

if (xc47 < avgx[identifikator]) {
	if (yc47 < avgy[identifikator]) {
		if (zc47 < avgz[identifikator]) {
			for (integer i_471 = 0; i_471 < i_22[identifikator][0]; i_471++) {
				integer i_47 = oct_load1[identifikator][0][i_471];
				if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) && (yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
					(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
				{
					ifound = i_47;
					bfound = true;
					break;
				}
			}
		}
		else {
			for (integer i_471 = 0; i_471 < i_22[identifikator][1]; i_471++) {
				integer i_47 = oct_load1[identifikator][1][i_471];
				if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) && (yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
					(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
				{
					ifound = i_47;
					bfound = true;
					break;
				}
			}
		}
	}
	else {
		if (zc47 < avgz[identifikator]) {
			for (integer i_471 = 0; i_471 < i_22[identifikator][2]; i_471++) {
				integer i_47 = oct_load1[identifikator][2][i_471];
				if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) && (yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
					(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
				{
					ifound = i_47;
					bfound = true;
					break; 
				}
			}
		}
		else {
			for (integer i_471 = 0; i_471 < i_22[identifikator][3]; i_471++) {
				integer i_47 = oct_load1[identifikator][3][i_471];
				if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) && (yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
					(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
				{
					ifound = i_47;
					bfound = true;
					break;
				}
			}
		}
	}
}
else {
	if (yc47 < avgy[identifikator]) {
		if (zc47 < avgz[identifikator]) {
			for (integer i_471 = 0; i_471 < i_22[identifikator][4]; i_471++) {
				integer i_47 = oct_load1[identifikator][4][i_471];
				if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) && (yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
					(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
				{
					ifound = i_47; 
					bfound = true; 
					break; 
				}
			}
		}
		else {
			for (integer i_471 = 0; i_471 < i_22[identifikator][5]; i_471++) {
				integer i_47 = oct_load1[identifikator][5][i_471];
				if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) && (yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
					(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
				{
					ifound = i_47;
					bfound = true;
					break;
				}
			}
		}
	}
	else {
		if (zc47 < avgz[identifikator]) {
			for (integer i_471 = 0; i_471 < i_22[identifikator][6]; i_471++) {
				integer i_47 = oct_load1[identifikator][6][i_471];
				if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) && (yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
					(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
				{
					ifound = i_47;
					bfound = true;
					break;
				}
			}
		}
		else {
					
						for (integer i_471 = 0; i_471 < i_22[identifikator][7]; i_471++) {
							
								integer i_47 = oct_load1[identifikator][7][i_471]; 
								if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) && (yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) && 
									(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
								{
									ifound = i_47; 
									bfound = true; 
									break; 
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
const unsigned char MU_LAME = 6;  // Коэффициент Лямэ
const unsigned char LAMBDA_LAME = 7; // Коэффициент Лямэ.
const unsigned char BETA_T_MECHANICAL = 8; // коэффициент линейного теплового расширения для Механики.


const unsigned char MU = 1; // динамическая вязкость
const unsigned char BETA_T = 2; // коэффициент линейного температурного расширения

const unsigned char VX = 0; // горизонтальная скорость (обязательно должен начинаться с нуля)
const unsigned char VY = 1; // вертикальная скорость
const unsigned char VZ = 2;
const unsigned char PRESS = 3; // давление
const unsigned char PAM = 4; // поправка давления
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
const unsigned char TOTALDEFORMATIONVAR = 27;
const unsigned char SPEED = 28; // Модуль скорости теплоносителя.
// Spalart-Allmares (RANS). 19.04.2019
const unsigned char NUSHA = 29; // Модифицированная кинематическая турбулентная вязкость.
// Градиенты модифицированной кинематической турбулентной вязкости.
const unsigned char GRADXNUSHA = 30;
const unsigned char GRADYNUSHA = 31;
const unsigned char GRADZNUSHA = 32;
// Static Structural
// Деформации:
const unsigned char TOTALDEFORMATION = 0; // полная
const unsigned char XDEFORMATION = 1; // х компонента
const unsigned char YDEFORMATION = 2; // y компонента
const unsigned char ZDEFORMATION = 3; // z компонента
// von-Mises


//const unsigned char SPEED = 99; // velocity magnitude (модуль скорости) константа используется при создании анимации.
const unsigned char TEMP = 100; // уравнение теплопроводности


const unsigned char ENUMERATECONTVOL = 0; // слой нумерации КО
const unsigned char MASKDOMAINFLUID = 1; // слой маска определяет зону FLUID
const unsigned char MASKDOMAINFLUIDCOLOR = 2; // В evt_f2 цвет гидродинамической подобласти.

const unsigned char LAMINAR = 0; // ламинарное течение
const unsigned char ZEROEQMOD = 1; // турбулентное течение - алгебраическая модель турбулентности Zero Equation Model
const unsigned char SMAGORINSKY = 2; // турбулентное течение - LES моделирование одна из разновидностей модели Смагоринского.
const unsigned char RNG_LES = 3; // Based on Renormalization Group Theory. (модель соответствует описанию CFD-Wiki).
// Динамическая модель Германо 1991 года. (основывается на модели Смагоринского и реализуется в виде её опции - bDynamic_Stress).

typedef struct TTOCKA_INT {
	integer i=0, j=0, k=0;
} TOCKA_INT;

// Объявление функции код которой будет реализован ниже.
// проверка построеной сетки
// экспорт результата расчёта в программу tecplot360
// части 1 и 3.
void exporttecplotxy360T_3D_part1and3(TEMPER &t, integer maxelm, integer maxbound, bool bextendedprint, integer ncell, integer** nvtx, integer** nvtxcell, TOCHKA* pa, BOUND* sosedb, integer ivarexport, integer** ptr_out);
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


	printf("delete temperature sosedi\n");
	if (t.sosedi != NULL) {
#pragma omp parallel for
		for (int i = 0; i<12; i++) {
			if (t.sosedi[i] != NULL) {
				delete[] t.sosedi[i]; // -12N
				t.sosedi[i] = NULL;
			}
		}
	}
	if (t.sosedi != NULL) {
		delete[] t.sosedi;
		t.sosedi = NULL;
	}

#pragma omp parallel sections
	{

		{
			if (t.whot_is_block != NULL) {
				delete[] t.whot_is_block;
				t.whot_is_block = NULL;
			}
		}

#pragma omp section 
		{
			printf("delete temperature Sc\n");
			if (t.Sc != NULL) {
				delete[] t.Sc; // -N
				t.Sc = NULL;
			}
		}

#pragma omp section 
		{
			printf("delete temperature ipower_time_depend\n");
			if (t.ipower_time_depend != NULL) {
				delete[] t.ipower_time_depend; // -N
				t.ipower_time_depend = NULL;
			}
		}


#pragma omp section 
		{
			printf("delete temperature pa\n");
			if (t.pa != NULL) {
				delete[] t.pa; // -3N
				t.pa = NULL;
			}
		}

#pragma omp section 
		{
			printf("delete temperature sosedb\n");
			if (t.sosedb != NULL) {
				delete[] t.sosedb;
				t.sosedb = NULL;
			}
		}

#pragma omp section 
		{
			printf("delete temperature binternal source\n");
			if (t.binternalsource != NULL) {
				delete[] t.binternalsource;
				t.binternalsource = NULL;
			}
		}

#pragma omp section 
		{
			printf("delete temperature prop\n");
			if (t.prop != NULL) {
				for (int i = 0; i < 9; i++) {
					if (t.prop[i] != NULL) {
						delete[] t.prop[i]; // -3N
					}
				}
			}
			if (t.prop != NULL) {
				delete[] t.prop;
				t.prop = NULL;
			}
		}


#pragma omp section 
		{
			printf("delete temperature ptr\n");
			if (t.ptr != NULL) {
				for (int i = 0; i < 2; i++) {
					if (t.ptr[i] != NULL) {
						delete[] t.ptr[i];
						t.ptr[i] = NULL;
					}
				}
				delete[] t.ptr;
				t.ptr = NULL;
			}
		}

#pragma omp section 
		{
			printf("delete temperature prop_b\n");
			if (t.prop_b != NULL) {
				for (int i = 0; i < 6; i++) {
					if (t.prop_b[i] != NULL) {
						delete[] t.prop_b[i]; // -3N
						t.prop_b[i] = NULL;
					}
				}
			}
			if (t.prop_b != NULL) {
				delete[] t.prop_b;
				t.prop_b = NULL;
			}
		}

#pragma omp section 
		{
			printf("delete temperature nvtx and nvtxcell\n");
			for (int i = 0; i < 8; i++) { // -8N
				if (t.nvtx != NULL) {
					if (t.nvtx[i] != NULL) {
						delete[] t.nvtx[i];
						t.nvtx[i] = NULL;
					}
				}

			}
			if (t.nvtx != NULL) {
				delete[] t.nvtx;
				t.nvtx = NULL;
			}
		}

#pragma omp section 
		{
			for (int i = 0; i < 8; i++) { // -8N

				if (t.nvtxcell != NULL) {
					if (t.nvtxcell[i] != NULL) {
						delete[] t.nvtxcell[i]; // может быть уже удалено
						t.nvtxcell[i] = NULL;
					}
				}
			}
			if (t.nvtxcell != NULL) {
				delete[] t.nvtxcell;
				t.nvtxcell = NULL;
			}
		}

#pragma omp section 
		{
			// 26_09_2016.
			printf("delete ilevel_alice\n");
			if (t.ilevel_alice != NULL) {
				delete[] t.ilevel_alice;
				t.ilevel_alice = NULL;
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
	printf("delete temperature slau\n");
	if (t.slau != NULL) {
		delete[] t.slau;
		t.slau = NULL;
	}
	printf("delete temperature slau_bon\n");
	if (t.slau_bon != NULL) {
		delete[] t.slau_bon;
		t.slau_bon = NULL;
	}

	printf("delete temperature potent\n");
	if (t.potent != NULL) {
		delete[] t.potent;
		t.potent = NULL;
	}

	printf("delete total deformation \n");
	if (t.total_deformation != NULL) {
		for (integer i_1 = 0; i_1 < 4; i_1++) {
			if (t.total_deformation[i_1] != NULL) {
				delete[] t.total_deformation[i_1];
				t.total_deformation[i_1] = NULL;
			}
		}

		delete[] t.total_deformation;
		t.total_deformation = NULL;
	}
} // free_level2_temp

// вычисление центральной точки центрального контрольного объёма
void center_cord3D(integer iP, integer** nvtx, TOCHKA* pa, TOCHKA &p, integer id) {
	// вычисление центральной точки центрального контрольного объёма
	if (iP < 0) {
		p.x = 0.0;
		p.y = 0.0;
		p.z = 0.0;
		printf("fatall error in center_cord3D !!! ");
		switch (id) {
		case WSIDE: printf("W"); break;
		case ESIDE: printf("E"); break;
		case NSIDE: printf("N"); break;
		case SSIDE: printf("S"); break;
		case BSIDE: printf("B"); break;
		case TSIDE: printf("T"); break;
		case WW: printf("WW"); break;
		case EE: printf("EE"); break;
		case NN: printf("NN"); break;
		case SS: printf("SS"); break;
		case BB: printf("BB"); break;
		case TTSIDE: printf("TT"); break;
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
		if (dz < 1.0e-40) {
			printf("ERROR Z : INCORRECT NUMERATION IN NVTX MAY BE...\n");
			//printf("error : center cord 3D slipanie po Z. delta < shorter_length_for_simplificationZ\n");
			printf("error : center cord 3D slipanie po Z. delta < 1.0e-40\n");
#if doubleintprecision == 1
			printf("iP=%lld dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
			printf("nvtx[%lld] : %lld %lld %lld %lld %lld %lld %lld %lld \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#else
			printf("iP=%d dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
			printf("nvtx[%d] : %d %d %d %d %d %d %d %d \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#endif

				//system("PAUSE");
			system("PAUSE");
		}
		//if (dx < shorter_length_for_simplificationX(p.x)) {
		if (dx < 1.0e-40) {
		    printf("ERROR X : INCORRECT NUMERATION IN NVTX MAY BE...\n");
			//printf("error : center cord 3D slipanie po X. delta < shorter_length_for_simplificationX\n");
			printf("error : center cord 3D slipanie po X. delta < 1.0e-40\n");
#if doubleintprecision == 1
			printf("iP=%lld dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
			printf("nvtx[%lld] : %lld %lld %lld %lld %lld %lld %lld %lld \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#else
			printf("iP=%d dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
			printf("nvtx[%d] : %d %d %d %d %d %d %d %d \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#endif
				//system("PAUSE");
			system("PAUSE");
		}
		//if (dy < shorter_length_for_simplificationY(p.y)) {
		if (dy < 1.0e-40) {
		    printf("ERROR Y: INCORRECT NUMERATION IN NVTX MAY BE...\n");
			//printf("error : center cord 3D slipanie po Y. delta < shorter_length_for_simplificationY\n");
			printf("error : center cord 3D slipanie po Y. delta <1.0e-40\n");
#if doubleintprecision == 1
			printf("iP=%lld dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
			printf("nvtx[%lld] : %lld %lld %lld %lld %lld %lld %lld %lld \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#else
			printf("iP=%d dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
			printf("nvtx[%d] : %d %d %d %d %d %d %d %d \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#endif
				//system("PAUSE");
			system("PAUSE");
		}
	}

	//printf("%e %e %e\n",dx,dy,dz); // debug GOOD
	//system("PAUSE");	
} // center_cord3D

// вычисляет размеры контрольного объёма
// 23.03.2019 Делаем более сильную проверку и на отрицательный КО тоже.
void volume3D(integer iP, integer** nvtx, TOCHKA* pa, doublereal &dx, doublereal &dy, doublereal &dz) {
	// вычисление размеров текущего контрольного объёма:
	
	TOCHKA p; // Для локальной привязки.
	// координаты центра контрольного объёма
	p.x = 0.5*(pa[nvtx[1][iP] - 1].x + pa[nvtx[0][iP] - 1].x);
	p.y = 0.5*(pa[nvtx[2][iP] - 1].y + pa[nvtx[0][iP] - 1].y);
	p.z = 0.5*(pa[nvtx[4][iP] - 1].z + pa[nvtx[0][iP] - 1].z);

	//dx = 0.0; dy = 0.0; dz = 0.0;// размеры контрольного объёма
	dx=pa[nvtx[1][iP]-1].x-pa[nvtx[0][iP]-1].x;
	dy=pa[nvtx[2][iP]-1].y-pa[nvtx[0][iP]-1].y;
    dz=pa[nvtx[4][iP]-1].z-pa[nvtx[0][iP]-1].z;
	//if (fabs(pa[nvtx[4][iP] - 1].z - pa[nvtx[0][iP] - 1].z) < 1.0e-40) {
	//if (dz < shorter_length_for_simplificationZ(p.z)) {
	if (dz < 1.0e-40) {
	    printf("ERROR 23.03.2019 NVTX perenumeration is ERROR!!!\n");
		//printf("error : volume 3D slipanie po Z. delta < shorter_length_for_simplificationZ\n");
		printf("error : volume 3D slipanie po Z. delta < 1.0e-40\n");
#if doubleintprecision == 1
		printf("iP=%lld dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
		printf("nvtx[%lld] : %lld %lld %lld %lld %lld %lld %lld %lld \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#else
		printf("iP=%d dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
		printf("nvtx[%d] : %d %d %d %d %d %d %d %d \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#endif
			//system("PAUSE");
		system("PAUSE");
	}
	//if (fabs(pa[nvtx[1][iP] - 1].x - pa[nvtx[0][iP] - 1].x) < 1.0e-40) {
	//if (dx < shorter_length_for_simplificationX(p.x)) {
	if (dx < 1.0e-40) {
	    printf("ERROR 23.03.2019 NVTX perenumeration is ERROR!!!\n");
		//printf("error : volume 3D slipanie po X. delta < shorter_length_for_simplificationX\n");
		printf("error : volume 3D slipanie po X. delta <  1.0e-40\n");
#if doubleintprecision == 1
		printf("iP=%lld dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
		printf("nvtx[%lld] : %lld %lld %lld %lld %lld %lld %lld %lld \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#else
		printf("iP=%d dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
		printf("nvtx[%d] : %d %d %d %d %d %d %d %d \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#endif
		//system("PAUSE");
		system("PAUSE");
	}
	//if (fabs(pa[nvtx[3][iP] - 1].y - pa[nvtx[0][iP] - 1].y) < 1.0e-40) {
//if (dy < shorter_length_for_simplificationY(p.y)) {
if (dy < 1.0e-40) {
	    printf("ERROR 23.03.2019 NVTX perenumeration is ERROR!!!\n");
		//printf("error : volume 3D slipanie po Y. delta < shorter_length_for_simplificationY\n");
		printf("error : volume 3D slipanie po Y. delta < 1.0e-40\n");
#if doubleintprecision == 1
		printf("iP=%lld dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
		printf("nvtx[%lld] : %lld %lld %lld %lld %lld %lld %lld %lld \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#else
		printf("iP=%d dx=%e dy=%e dz=%e\n", iP, dx, dy, dz);
		printf("nvtx[%d] : %d %d %d %d %d %d %d %d \n", iP, nvtx[0][iP] - 1, nvtx[1][iP] - 1, nvtx[2][iP] - 1, nvtx[3][iP] - 1, nvtx[4][iP] - 1, nvtx[5][iP] - 1, nvtx[6][iP] - 1, nvtx[7][iP] - 1);
#endif
		//system("PAUSE");
		system("PAUSE");
	}
	//printf("%e %e %e\n",dx,dy,dz); // debug GOOD
	//system("PAUSE");	
} // volume3D 

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
	doublereal *x11=NULL, *y11=NULL, *z11=NULL;
	integer ix11 = 0, iy11 = 0, iz11 = 0;
	integer *ijk_block_id = NULL;
	doublereal epsTolx = 1.0e40, epsToly = 1.0e40, epsTolz = 1.0e40;
	integer *b_non_prism=NULL;
	integer lb_non_prism=0;
} QuickSearchBlockid;

QuickSearchBlockid QSBid; 

void init_QSBid(integer lb, BLOCK* &b, WALL* &w, integer &lw, SOURCE* &s, integer &ls) {

	// Вычисляем количество непризматических объектов QSBid.lb_non_prism и заносим их в отдельный список QSBid.b_non_prism.
	for (integer i = 0; i < lb; i++) {
		if (b[i].g.itypegeom != PRISM) {
			QSBid.lb_non_prism++;
		}
	}
	if (QSBid.lb_non_prism > 0) {
		QSBid.b_non_prism = new integer[QSBid.lb_non_prism];
	}
	integer j = 0;
	for (integer i = 0; i < lb; i++) {
		if (b[i].g.itypegeom != PRISM) {
			QSBid.b_non_prism[j] = i;
			j++;
		}
	}

	// X
	QSBid.ix11 = 1;
	QSBid.x11 = new doublereal[QSBid.ix11 + 1]; // число границ по оси x

	QSBid.x11[0] = b[0].g.xS; // начало области
	QSBid.x11[QSBid.ix11] = b[0].g.xE; // конец области

	for (integer i = 1; i < lb; i++) {

		if (b[i].g.itypegeom == PRISM) {
			doublereal x_1 = b[i].g.xS;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(QSBid.x11, QSBid.ix11, x_1, YZ,b,lb, w, lw, s, ls);
			}
			x_1 = b[i].g.xE;
			if ((x_1 >= b[0].g.xS) && (x_1 <= b[0].g.xE)) {
				addboundary(QSBid.x11, QSBid.ix11, x_1, YZ, b, lb, w, lw, s, ls);
			}
		}
	}

	Sort_method(QSBid.x11, QSBid.ix11);

	// Y
	QSBid.iy11 = 1;
	QSBid.y11 = new doublereal[QSBid.iy11 + 1]; // число границ по оси y

	QSBid.y11[0] = b[0].g.yS; // начало области
	QSBid.y11[QSBid.iy11] = b[0].g.yE; // конец области

	for (integer i = 1; i < lb; i++) {

		if (b[i].g.itypegeom == PRISM) {
			doublereal y_1 = b[i].g.yS;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(QSBid.y11, QSBid.iy11, y_1,XZ, b, lb, w, lw, s, ls);
			}
			y_1 = b[i].g.yE;
			if ((y_1 >= b[0].g.yS) && (y_1 <= b[0].g.yE)) {
				addboundary(QSBid.y11, QSBid.iy11, y_1, XZ, b, lb, w, lw, s, ls);
			}
		}
	}

	Sort_method(QSBid.y11, QSBid.iy11);

	// Z
	QSBid.iz11 = 1;
	QSBid.z11 = new doublereal[QSBid.iz11 + 1]; // число границ по оси z

	QSBid.z11[0] = b[0].g.zS; // начало области
	QSBid.z11[QSBid.iz11] = b[0].g.zE; // конец области

	for (integer i = 1; i < lb; i++) {

		if (b[i].g.itypegeom == PRISM) {
			doublereal z_1 = b[i].g.zS;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(QSBid.z11, QSBid.iz11, z_1,XY, b, lb, w, lw, s, ls);
			}
			z_1 = b[i].g.zE;
			if ((z_1 >= b[0].g.zS) && (z_1 <= b[0].g.zE)) {
				addboundary(QSBid.z11, QSBid.iz11, z_1,XY, b, lb, w, lw, s, ls);
			}
		}
	}

	Sort_method(QSBid.z11, QSBid.iz11);

	QSBid.ijk_block_id = new integer[(QSBid.ix11 + 1)*(QSBid.iy11 + 1)*(QSBid.iz11 + 1)];
#pragma omp parallel for
	for (integer i_94 = 0; i_94 < (QSBid.ix11 + 1)*(QSBid.iy11 + 1)*(QSBid.iz11 + 1); i_94++) {
		QSBid.ijk_block_id[i_94] = lb;
	}

	// Вычисление допуска.
	QSBid.epsTolx = 1.0e40;
	QSBid.epsToly = 1.0e40;
	QSBid.epsTolz = 1.0e40;
	for (integer i = 0; i < QSBid.ix11; i++) {
		if (fabs(QSBid.x11[i + 1] - QSBid.x11[i]) < QSBid.epsTolx) {
			QSBid.epsTolx = 0.5*fabs(QSBid.x11[i + 1] - QSBid.x11[i]);
		}
	}
	for (integer i = 0; i < QSBid.iy11; i++) {
		if (fabs(QSBid.y11[i + 1] - QSBid.y11[i]) < QSBid.epsToly) {
			QSBid.epsToly = 0.5*fabs(QSBid.y11[i + 1] - QSBid.y11[i]);
		}
	}
	for (integer i = 0; i < QSBid.iz11; i++) {
		if (fabs(QSBid.z11[i + 1] - QSBid.z11[i]) < QSBid.epsTolz) {
			QSBid.epsTolz = 0.5*fabs(QSBid.z11[i + 1] - QSBid.z11[i]);
		}
	}

	/*
	QSBid.epsToolx = -1.0e40;
	QSBid.epsTooly = -1.0e40;
	QSBid.epsToolz = -1.0e40;
	for (integer i = 0; i < inx; i++) {
		if (fabs(xpos[i + 1] - xpos[i]) > QSBid.epsToolx) {
			QSBid.epsToolx = 0.5*fabs(xpos[i + 1] - xpos[i]);
		}
	}
	for (integer i = 0; i < iny; i++) {
		if (fabs(ypos[i + 1] - ypos[i]) > QSBid.epsTooly) {
			QSBid.epsTooly = 0.5*fabs(ypos[i + 1] - ypos[i]);
		}
	}
	for (integer i = 0; i < inz; i++) {
		if (fabs(zpos[i + 1] - zpos[i]) > QSBid.epsToolz) {
			QSBid.epsToolz = 0.5*fabs(zpos[i + 1] - zpos[i]);
		}
	}
	*/

	Block_indexes* block_indexes = new Block_indexes[lb];
	// Оператор new не требует проверки.
		//if (block_indexes == NULL) {
			//printf("error in allocation memory for block_indexes in enumerate_volume_improved.\n");
			//system("pause");
			//exit(1);
		//}

	// initialization 30.07.2019
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
			for (j = 0; j <= QSBid.ix11; j++) {
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
				for (j = 0; j <= QSBid.ix11; j++) {
					if ((!bfound) && (fabs(x4 - QSBid.x11[j])<shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].iL = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = 0; j <= QSBid.ix11; j++) {
					if ((!bfound) && (x4 > QSBid.x11[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].iL = j;
						bfound = true;
						break;
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
			for (j = 0; j <= QSBid.ix11; j++) {
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
				for (j = QSBid.ix11; j >= 0; j--) {
					if ((!bfound) && (fabs(x4 - QSBid.x11[j])<shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].iR = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = QSBid.ix11; j >= 0; j--) {
					if ((!bfound) && (x4 <= QSBid.x11[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].iR = j;
						bfound = true;
						break;
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
			for (j = 0; j <= QSBid.iy11; j++) {
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
				for (j = 0; j <= QSBid.iy11; j++) {
					if ((!bfound) && (fabs(x4 - QSBid.y11[j])<shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].jL = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = 0; j <= QSBid.iy11; j++) {
					if ((!bfound) && (x4 > QSBid.y11[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].jL = j;
						bfound = true;
						break;
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
			for (j = 0; j <= QSBid.iy11; j++) {

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
				for (j = QSBid.iy11; j >= 0; j--) {
					if ((!bfound) && (fabs(x4 - QSBid.y11[j])<shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].jR = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = QSBid.iy11; j >= 0; j--) {
					if ((!bfound) && (x4 <= QSBid.y11[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].jR = j;
						bfound = true;
						break;
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
			for (j = 0; j <= QSBid.iz11; j++) {
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
				for (j = 0; j <= QSBid.iz11; j++) {
					if ((!bfound) && (fabs(x4 - QSBid.z11[j])<shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].kL = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = 0; j <= QSBid.iz11; j++) {
					if ((!bfound) && (x4 > QSBid.z11[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].kL = j;
						bfound = true;
						break;
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
			for (j = 0; j <= QSBid.iz11; j++) {
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
				for (j = QSBid.iz11; j >= 0; j--) {
					if ((!bfound) && (fabs(x4 - QSBid.z11[j])<shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls))) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].kR = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = QSBid.iz11; j >= 0; j--) {
					if ((!bfound) && (x4 <= QSBid.z11[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i].kR = j;
						bfound = true;
						break;
					}
				}
			}

			if ((block_indexes[i].iL>=0)&&
				(block_indexes[i].iR>=0)&&
				(block_indexes[i].iL >= block_indexes[i].iR)) {
				printf("init_QSBid function\n");
				printf("violation of the order\n");
				printf("i=%lld iL=%lld iR=%lld delta=%e\n", i, block_indexes[i].iL,
					block_indexes[i].iR, shorter_length_for_simplificationX(QSBid.x11[block_indexes[i].iR], b, lb, w, lw, s, ls));
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
				printf("i=%lld jL=%lld jR=%lld delta=%e\n", i, block_indexes[i].jL,
					block_indexes[i].jR, shorter_length_for_simplificationY(QSBid.y11[block_indexes[i].jR], b, lb, w, lw, s, ls));
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
		printf("STATISTICS: X=%lld Y=%lld Z=%lld\n", iX_one_CELL_count_statistic, iY_one_CELL_count_statistic, iZ_one_CELL_count_statistic);
		// Работает и при такой ситуации. Проверено на одной из моделей.
		//system("pause");
	}


	// block_indexes подготовлен.

	bool* bvisit = NULL;
	bvisit = new bool[(QSBid.ix11 + 1)*(QSBid.iy11 + 1)*(QSBid.iz11 + 1)];

	// Количество проходов существенно сократилось и в итоге это приводит к существенному
		// увеличению быстродействия.
	integer m7 = lb - 1, m8;

#pragma omp parallel for
	for (integer iP = 0; iP < (QSBid.ix11 + 1)*(QSBid.iy11 + 1)*(QSBid.iz11 + 1); iP++) {
		bvisit[iP] = false;
	}




	for (m8 = lb - 1; m8 >= 0; m8--) {
		m7 = m8;
		if (b[m8].g.itypegeom == PRISM)
		{
			//#pragma omp parallel for
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++)
				for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++)
					for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++)
					{
						integer iP = i1 + j1 * (QSBid.ix11 + 1) + k1 * (QSBid.ix11 + 1)*(QSBid.iy11 + 1);

						if ((i1 < 0) || (i1 > QSBid.ix11) || (j1 < 0) || (j1 > QSBid.iy11) || (k1 < 0) || (k1 > QSBid.iz11)) {
							// ERROR
							printf("ERROR PRISM\n");
							printf("inx=%lld iny=%lld inz=%lld \n", QSBid.ix11, QSBid.iy11, QSBid.iz11);
							printf("i1=%lld j1=%lld k1=%lld \n", i1, j1, k1);
							printf("iP=%lld m8=%lld", iP, m8);
							printf("iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
							system("PAUSE");
						}

						if (bvisit[iP] == false)
						{

							bvisit[iP] = true;

							QSBid.ijk_block_id[iP] = m8;
						}
					}
			//m7--;
		}
	}

	
	// Так делать нельзя т.к. мы заполняли только на призмах и не рассматривали цилиндры и полигоны.
	// Поэтому какие-то позиции из ijk_block_id[iP] отстанутся равными lb т.е. незаполненными.
	// Это нормально. НО!!! здесь мы добавили проверку if (bvisit[iP] == true)  и поэтому гарантируется
	// Что блоки только призматические и только.
	for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++)
	{
		for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++)
		{
			for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++)
			{
				integer iP = i1 + j1 * (QSBid.ix11 + 1) + k1 * (QSBid.ix11 + 1)*(QSBid.iy11 + 1);
				if (bvisit[iP] == true) {
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
	bvisit = NULL;

	delete[] block_indexes;
	block_indexes = NULL;

}

void free_QSBid() {

	if (QSBid.b_non_prism!=NULL) {
	    delete[] QSBid.b_non_prism;
	    QSBid.b_non_prism = NULL;
    }

	if (QSBid.x11 != NULL) {
		delete[] QSBid.x11;
		QSBid.x11 = NULL;
	}
	if (QSBid.y11 != NULL) {
		delete[] QSBid.y11;
		QSBid.y11 = NULL;
	}
	if (QSBid.z11 != NULL) {
		delete[] QSBid.z11;
		QSBid.z11 = NULL;
	}

	if (QSBid.ijk_block_id != NULL) {
		delete[] QSBid.ijk_block_id;
		QSBid.ijk_block_id = NULL;
	}

	QSBid.epsTolx = 1.0e40;
	QSBid.epsToly = 1.0e40;
	QSBid.epsTolz = 1.0e40;

	QSBid.ix11 = 0;
	QSBid.iy11 = 0;
	QSBid.iz11 = 0;
}

// Какому прямоугольному параллелепипеду принадлежит точка p.
integer myisblock_id_PRISM_only(integer lb, BLOCK* &b, TOCHKA p) {
	// return i_vacant + (inx + 1)*j_vacant + (inx + 1)*(iny + 1)*k_vacant;
	integer ikey = hash_key_alice33Q(QSBid.ix11, QSBid.iy11, QSBid.iz11, QSBid.x11, QSBid.y11, QSBid.z11, p, QSBid.epsTolx, QSBid.epsToly, QSBid.epsTolz);
	integer ib= QSBid.ijk_block_id[ikey];
	if ((ib < 0) || (ib >= lb)) {
		printf("error in function myisblock_id_PRISM_only.\n");
		printf("p.x=%e p.y=%e p.z=%e ikey==%lld\n",p.x,p.y,p.z,ikey);
		printf("ikey=%lld size=%lld ib=%lld\n",ikey, (QSBid.ix11 + 1)*(QSBid.iy11 + 1)*(QSBid.iz11 + 1), ib);
		printf("lx=%lld ly=%lld iz=%lld\n", QSBid.ix11, QSBid.iy11, QSBid.iz11);
		
		// Вызываем очень медленный наивный алгоритм.
		ib = lb;
		for (integer i_1 = lb - 1; i_1 >= 0; i_1--) {
			if (b[i_1].g.itypegeom == PRISM) {
				if ((b[i_1].g.xS < p.x) && (b[i_1].g.xE > p.x) && (b[i_1].g.yS < p.y) && (b[i_1].g.yE > p.y) && (b[i_1].g.zS < p.z) && (b[i_1].g.zE > p.z)) {
					ib = i_1;
					// Если блок найден то сканирование сразу прекращается.
					break;
				}
			}
		}
		printf("ib==%lld lb=%lld\n",ib,lb);
		//system("PAUSE");
	}
	return ib;
} // myisblock_id_PRISM_only

// Полнофункциональная медленная версия работающая и для цилиндров и для полигонов,
// а также для прямоугольных параллелепипедов.
// 01.05.2019
integer myisblock_id(integer lb, BLOCK* &b, TOCHKA p) {
	integer ib = 0;
	//for (integer i_1 = lb - 1; i_1 >= 0; i_1--) {
		//if ((b[i_1].g.xS < p.x) && (b[i_1].g.xE > p.x) && (b[i_1].g.yS < p.y) && (b[i_1].g.yE > p.y) && (b[i_1].g.zS < p.z) && (b[i_1].g.zE > p.z)) {
			//ib = i_1;
			// Если блок найден то сканирование сразу прекращается.
			//break;
		//}
	//}

	integer i = 0, k = 0;

	// Параллельный алгоритм существенно медленней однопоточного.
	// Однопоточный вариант полностью пересмотрен и использует двоичный поиск и хеш таблицы.
	// Откажемся от многопоточного варианта т.к. он неэффективен.
	//integer lb4 = lb / 4; // Это очень плохое деление пригодное только в случае,
	// если все блоки равноправны : содержат примерно одинаковое число ячеек.

/*
#ifdef _OPENMP
	integer k1 = 0, k2 = 0, k3 = 0, k4 = 0;
#pragma omp parallel shared(k1,k2,k3,k4) num_threads(4)
	{
		integer tid = omp_get_thread_num();
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

				if (b[i1].g.itypegeom == CYLINDER) {
					// Cylinder
					switch (b[i1].g.iPlane) {
					case XY:
						if (fabs(b[i1].g.R_in_cyl) < 1.0e-40) {
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
						if (fabs(b[i1].g.R_in_cyl) < 1.0e-40) {
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
						if (fabs(b[i1].g.R_in_cyl) < 1.0e-40) {
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

				if (b[i2].g.itypegeom == CYLINDER) {
					// Cylinder
					switch (b[i2].g.iPlane) {
					case XY:
						if (fabs(b[i2].g.R_in_cyl) < 1.0e-40) {
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
						if (fabs(b[i2].g.R_in_cyl) < 1.0e-40) {
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
						if (fabs(b[i2].g.R_in_cyl) < 1.0e-40) {
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

				if (b[i3].g.itypegeom == CYLINDER) {
					// Cylinder
					switch (b[i3].g.iPlane) {
					case XY:
						if (fabs(b[i3].g.R_in_cyl) < 1.0e-40) {
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
						if (fabs(b[i3].g.R_in_cyl) < 1.0e-40) {
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
						if (fabs(b[i3].g.R_in_cyl) < 1.0e-40) {
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

				if (b[i4].g.itypegeom == CYLINDER) {
					// Cylinder
					switch (b[i4].g.iPlane) {
					case XY:
						if (fabs(b[i4].g.R_in_cyl) < 1.0e-40) {
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
						if (fabs(b[i4].g.R_in_cyl) < 1.0e-40) {
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
						if (fabs(b[i4].g.R_in_cyl) < 1.0e-40) {
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
   integer kprism = myisblock_id_PRISM_only(lb,b,p); 
   k = lb; ib = lb;

   bool bfound_out = false;
   //integer *b_non_prism = NULL;
   //integer lb_non_prism = 0;
   for (integer i_74 = QSBid.lb_non_prism - 1; i_74 >= 0; i_74--) {
	   i = QSBid.b_non_prism[i_74]; // Номер непризматического объекта.
	   if (b[i].g.itypegeom == POLYGON) {

		   if ((p.x > b[i].g.xS) && (p.x < b[i].g.xE) && (p.y > b[i].g.yS) && (p.y < b[i].g.yE) && (p.z > b[i].g.zS) && (p.z < b[i].g.zE)) {

			   bool bfound = false;
			   // определяет принадлежность точки полигону.
			   bfound = in_polygon(p, b[i].g.nsizei, b[i].g.xi, b[i].g.yi, b[i].g.zi, b[i].g.hi, b[i].g.iPlane_obj2, k, i);
			   if (bfound) {
				   // Нашли и сразу завершили проверку в случае успеха.
				   bfound_out = true;
				   goto OUTOF_IN_MODEL_TEMP1;
			   }
		   }

	   }
	   if (b[i].g.itypegeom == CYLINDER) {
		   // Cylinder
		   switch (b[i].g.iPlane) {
		   case XY:
			   if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
				   if ((p.z >= b[i].g.zC) && (p.z <= b[i].g.zC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.yC - p.y)*(b[i].g.yC - p.y)) < b[i].g.R_out_cyl) {
						   k = i;
						   bfound_out = true;
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
							   bfound_out = true;
							   // Нашли и сразу завершили проверку в случае успеха.
							   goto OUTOF_IN_MODEL_TEMP1;
						   }
					   }
				   }
			   }
			   break;
		   case XZ:
			   if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
				   if ((p.y >= b[i].g.yC) && (p.y <= b[i].g.yC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.xC - p.x)*(b[i].g.xC - p.x) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
						   k = i;
						   bfound_out = true;
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
							   bfound_out = true;
							   // Нашли и сразу завершили проверку в случае успеха.
							   goto OUTOF_IN_MODEL_TEMP1;
						   }
					   }
				   }
			   }
			   break;
		   case YZ:
			   if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
				   if ((p.x >= b[i].g.xC) && (p.x <= b[i].g.xC + b[i].g.Hcyl)) {
					   if (sqrt((b[i].g.yC - p.y)*(b[i].g.yC - p.y) + (b[i].g.zC - p.z)*(b[i].g.zC - p.z)) < b[i].g.R_out_cyl) {
						   k = i;
						   bfound_out = true;
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
							   bfound_out = true;
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
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
	

OUTOF_IN_MODEL_TEMP1:

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
integer myisblock_id(integer lb, BLOCK* &b, doublereal x11, doublereal y11, doublereal z11) {
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
bool in_model_temp(TOCHKA p, integer &ib, BLOCK* b, integer lb) {
	
	bool ret = true;// по умолчанию принадлежит модели
	ib = myisblock_id(lb, b, p);
	
	if (ib == lb) {
		ret = false;
	}
	else {
		if (b[ib].itype == HOLLOW) ret = false;
	}
	return ret;

} // in_model_temp

// проверяет принадлежит ли контрольный объём
// гидродинамической модели.
// Возвращает параметр ib равный номеру блока
// которому принадлежит контрольный объём.
bool in_model_flow(TOCHKA p, integer &ib, BLOCK* b, integer lb) {
	
	bool ret=true;// по умолчанию принадлежит модели
	ib = myisblock_id(lb, b, p);
	   	
	if (ib == lb) {
		ret = false;
	}
	else {
		if ((b[ib].itype == SOLID) || (b[ib].itype == HOLLOW)) ret = false;
	}
	return ret;

} // in_model_flow

// проверяет принадлежит ли контрольный объём
// тепловой модели.
// Возвращает параметр ib равный номеру блока
// которому принадлежит контрольный объём.
bool in_model_temp_stab(TOCHKA p, integer &ib, BLOCK* b, integer lb) {

	integer i = 0, k = 0;
	bool ret = true;// по умолчанию принадлежит модели
	integer lb4 = lb / 4; // Это очень плохое деление пригодное только в случае,
	// если все блоки равноправны : содержат примерно одинаковое число ячеек.

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
	if (bYES_PARALLELISM_ON) {
		omp_set_num_threads(4);
		integer k1 = 0, k2 = 0, k3 = 0, k4 = 0;
#pragma omp parallel shared(k1,k2,k3,k4) num_threads(4)
		{
			integer tid = omp_get_thread_num();
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

					if (b[i1].g.itypegeom == CYLINDER) {
						// Cylinder
						switch (b[i1].g.iPlane) {
						case XY:
							if (fabs(b[i1].g.R_in_cyl) < 1.0e-40) {
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
							if (fabs(b[i1].g.R_in_cyl) < 1.0e-40) {
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
							if (fabs(b[i1].g.R_in_cyl) < 1.0e-40) {
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

					if (b[i2].g.itypegeom == CYLINDER) {
						// Cylinder
						switch (b[i2].g.iPlane) {
						case XY:
							if (fabs(b[i2].g.R_in_cyl) < 1.0e-40) {
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
							if (fabs(b[i2].g.R_in_cyl) < 1.0e-40) {
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
							if (fabs(b[i2].g.R_in_cyl) < 1.0e-40) {
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

					if (b[i3].g.itypegeom == CYLINDER) {
						// Cylinder
						switch (b[i3].g.iPlane) {
						case XY:
							if (fabs(b[i3].g.R_in_cyl) < 1.0e-40) {
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
							if (fabs(b[i3].g.R_in_cyl) < 1.0e-40) {
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
							if (fabs(b[i3].g.R_in_cyl) < 1.0e-40) {
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

					if (b[i4].g.itypegeom == CYLINDER) {
						// Cylinder
						switch (b[i4].g.iPlane) {
						case XY:
							if (fabs(b[i4].g.R_in_cyl) < 1.0e-40) {
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
							if (fabs(b[i4].g.R_in_cyl) < 1.0e-40) {
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
							if (fabs(b[i4].g.R_in_cyl) < 1.0e-40) {
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

		omp_set_num_threads(1);
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
		if (b[i].g.itypegeom == CYLINDER) {
			// Cylinder
			switch (b[i].g.iPlane) {
			case XY:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
		if (b[i].g.itypegeom == CYLINDER) {
			// Cylinder
			switch (b[i].g.iPlane) {
			case XY:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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

	if (b[k].itype == HOLLOW) ret = false;
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
bool in_model_flow_stab(TOCHKA p, integer &ib, BLOCK* b, integer lb) {
	integer i = 0, k = 0;
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

		if (b[i].g.itypegeom == CYLINDER) {
			// Cylinder
			switch (b[i].g.iPlane) {
			case XY:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
			case XZ:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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
			case YZ:
				if (fabs(b[i].g.R_in_cyl) < 1.0e-40) {
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

	if ((b[k].itype == SOLID) || (b[k].itype == HOLLOW)) ret = false;
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
void enumerate_volume_improved(integer* &evt, integer &maxelm, integer iflag,
	doublereal* xpos, doublereal* ypos, doublereal* zpos, integer* &whot_is_block,
	integer inx, integer iny, integer inz, BLOCK* b, integer lb,
	WALL* &w, integer &lw, SOURCE* &s, integer &ls) {

#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_num_threads();
	omp_set_num_threads(8); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
#else
	int i_my_num_core_parallelesation = 1;
#endif
	

	evt = NULL;
	evt = new integer[inx*iny*inz];
	if (evt == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for evt constr struct...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}
	whot_is_block = NULL;
	whot_is_block = new integer[inx*iny*inz];
	if (whot_is_block == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for whot_is_block constr struct...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}
	

	printf("enumerate_volume_improved start.\n");
	Block_indexes* block_indexes = new Block_indexes[lb];
	// оператор new не требует проверки на null.
	//if (block_indexes==NULL) {
		//printf("error in allocation memory for block_indexes in enumerate_volume_improved.\n");
		//system("pause");
		//exit(1);
	//}
	integer i, j, k;

	// 08.04.2018
	for (i = 0; i < lb; i++) {
		// инициализация, на случай если блоки не будут распознаны.
		block_indexes[i].iL = -1;
		block_indexes[i].iR = -2;
		block_indexes[i].jL = -1;
		block_indexes[i].jR = -2;
		block_indexes[i].kL = -1;
		block_indexes[i].kR = -2;
	}

	for (i = 0; i < lb; i++) {
		doublereal x4 = b[i].g.xS;
		for (j = 0; j <= inx; j++) {
			if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].iL = j;
				break;
			}
		}
		x4 = b[i].g.xE;
		for (j = 0; j <= inx; j++) {
			if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].iR = j;
				break;
			}
		}
		x4 = b[i].g.yS;
		for (j = 0; j <= iny; j++) {
			if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].jL = j;
				break;
			}
		}
		x4 = b[i].g.yE;
		for (j = 0; j <= iny; j++) {
			if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].jR = j;
				break;
			}
		}
		x4 = b[i].g.zS;
		for (j = 0; j <= inz; j++) {
			if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].kL = j;
				break;
			}
		}
		x4 = b[i].g.zE;
		for (j = 0; j <= inz; j++) {
			if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].kR = j;
				break;
			}
		}

		if ((block_indexes[i].iL >= 0) &&
			(block_indexes[i].iR >= 0) &&
			(block_indexes[i].iL >= block_indexes[i].iR)) {
			printf("enumerate_volume_improved function\n");
			printf("violation of the order block_indexes\n");
			printf("i=%lld iL=%lld iR=%lld\n", i, block_indexes[i].iL,
				block_indexes[i].iR);
			system("pause");
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
			printf("i=%lld jL=%lld jR=%lld\n", i, block_indexes[i].jL,
				block_indexes[i].jR);
			system("pause");
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
			printf("i=%lld kL=%lld kR=%lld\n", i, block_indexes[i].kL,
				block_indexes[i].kR);
			system("pause");
		}


		if ((block_indexes[i].kL >= 0) &&
			(block_indexes[i].kR >= 0) &&
			(block_indexes[i].kL + 1 == block_indexes[i].kR)) {
			// Всего одна клетка на блок.
			//iZ_one_CELL_count_statistic++;
		}

		if (0 && (i == 116)) {
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
			printf("enumerate_volume_improved function\n");
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
	integer m7;

	/*
	integer ib_stub = -1;
	// Мы найдем самый большой по размеру Hollow block, иначе будет просто кабинет.
	ib_stub = 0;
	doublereal vol_stub = -1.0;
	for (i = 0; i < lb; i++) {
		if (b[i].itype == HOLLOW) {
			if (fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS) > vol_stub) {
				ib_stub = i;
				vol_stub = fabs(b[i].g.xE - b[i].g.xS)*fabs(b[i].g.yE - b[i].g.yS)*fabs(b[i].g.zE - b[i].g.zS);
			}
		}
	}
	*/

#pragma omp parallel for
	for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
		evt[i1 + j1*inx + k1*inx*iny] = -1;// -1
	}
	for (m7 = 0; m7 < lb; m7++) {

#pragma omp parallel for
		for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++) for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++) for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++) {
			switch (iflag) {
			case TEMPERATURE:
				if (b[m7].itype == HOLLOW) {
					evt[i1 + j1*inx + k1*inx*iny] = -1;
				}
				else {
					evt[i1 + j1*inx + k1*inx*iny] = m7;
				}
					break;
			case HYDRODINAMIC:
				if ((b[m7].itype == SOLID) || (b[m7].itype == HOLLOW)) {
					evt[i1 + j1*inx + k1*inx*iny] = -1;
				}
				else {
					evt[i1 + j1*inx + k1*inx*iny] = m7;
				}
				break;
			}
		
		}
	}
	printf("enumerate_volume_improved 80 procent.\n");
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	integer l = 1, ib;
	for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
		if (evt[i + j*inx + k*inx*iny] > -1) {
			ib = evt[i + j*inx + k*inx*iny]; // номер блока был сохранён ранее.
			// Это очень нужно для записи репорта.
			whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
			evt[i + j*inx + k*inx*iny] = l;
			l++;
		}
		else
		{   // не принадлежит расчётной области
			evt[i + j*inx + k*inx*iny] = 0;
		}
	}
	maxelm = l - 1;
	
	delete[] block_indexes;

	printf("enumerate_volume_improved end.\n");

	/*
	bool inDomain = false;
	// нумерация контрольных объёмов начинается с единицы
	// если контрольный объём не принадлежит расчётной области то ставится 0.
	for (i = 0; i<inx; i++) for (j = 0; j<iny; j++) for (k = 0; k<inz; k++) {
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
void enumerate_volume_improved_obobshenie(integer* &evt, integer &maxelm, integer iflag,
	doublereal* xpos, doublereal* ypos, doublereal* zpos, integer* &whot_is_block,
	integer inx, integer iny, integer inz, BLOCK* b, integer lb, TOCKA_INT* &tck_int_list,
	WALL* &w, integer &lw, SOURCE* &s, integer &ls) {

#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_num_threads();
	omp_set_num_threads(8); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
#else
	int i_my_num_core_parallelesation = 1;
#endif

	tck_int_list = new TOCKA_INT[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (tck_int_list == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for evt tck_int_list...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	bool* bvisit = NULL;
	bvisit = new bool[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (bvisit == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	evt = NULL;
	evt = new integer[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (evt == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for evt constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}
	whot_is_block = NULL;
	whot_is_block = new integer[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (whot_is_block == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for whot_is_block constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}


	printf("enumerate_volume_improved_obobshenie start.\n");
	Block_indexes* block_indexes = new Block_indexes[lb];
	// оператор new не требует проверки на null.
	//if (block_indexes == NULL) {
		//printf("error in allocation memory for block_indexes in enumerate_volume_improved.\n");
		//system("pause");
		//exit(1);
	//}
	integer i=0, j=0, k=0, i_1=lb-1;//lb-1

	// 08.04.2018
	for (i = 0; i < lb; i++) {
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

	//for (i = 0; i < lb; i++) {
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
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == XZ))) {
				x4 = b[i].g.xC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == YZ))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.xC;
				}
				else {
					x4 = b[i].g.xC+ b[i].g.Hcyl;
				}
			}
			bool bfound = false;
			/*
			for (j = 0; j <= inx; j++) {
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
							printf("ERROR : if (fabs(100 * (xpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps)\n");
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
				for (j = 0; j <= inx; j++) {
					if ((!bfound) && (fabs(x4 - xpos[j])<shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)))
					{
						// Точное совпадение !!!
						block_indexes[i_1].iL = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = 0; j <= inx; j++) {
					if ((!bfound) && (x4 > xpos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].iL = j;
						bfound = true;
						break;
					}
				}
			}
			x4 = b[i].g.xE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == XZ))) {
				x4 = b[i].g.xC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == YZ))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.xC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.xC;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= inx; j++) {
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
							printf("ERROR : if (fabs(100 * (xpos[j] - x4) / fabs(x4)) < otnositelnaq_tolerance_eps)\n");
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
				for (j = inx; j >= 0; j--) {
					if ((!bfound) && (fabs(x4 - xpos[j])<shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls))) {
						// Точное совпадение !!!
						block_indexes[i_1].iR = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = inx; j >= 0; j--) {
					if ((!bfound) && (x4 < xpos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].iR = j;
						bfound = true;
						break;
					}
				}
			}
			x4 = b[i].g.yS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == YZ))) {
				x4 = b[i].g.yC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.yC;
				}
				else {
					x4 = b[i].g.yC + b[i].g.Hcyl;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= iny; j++) {
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
				for (j = 0; j <= iny; j++) {
					if ((!bfound) && (fabs(x4 - ypos[j])<shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls))) {
						// Точное совпадение !!!
						block_indexes[i_1].jL = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = 0; j <= iny; j++) {
					if ((!bfound) && (x4 > ypos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].jL = j;
						bfound = true;
						break;
					}
				}
			}
			x4 = b[i].g.yE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == YZ))) {
				x4 = b[i].g.yC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.yC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.yC ;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= iny; j++) {
				
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
				for (j = iny; j >= 0; j--) {
					if ((!bfound) && (fabs(x4 - ypos[j])<shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls))) {
						// Точное совпадение !!!
						block_indexes[i_1].jR = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = iny; j >= 0; j--) {
					if ((!bfound) && (x4 < ypos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].jR = j;
						bfound = true;
						break;
					}
				}
			}
			x4 = b[i].g.zS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ) || (b[i].g.iPlane == YZ))) {
				x4 = b[i].g.zC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.zC;
				}
				else {
					x4 = b[i].g.zC + b[i].g.Hcyl;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= inz; j++) {
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
				for (j = 0; j <= inz; j++) {
					if ((!bfound) && (fabs(x4 - zpos[j])<shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls))) {
						// Точное совпадение !!!
						block_indexes[i_1].kL = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = 0; j <= inz; j++) {
					if ((!bfound) && (x4 > zpos[j])) {
						// Нет точного совпаднения первая встреча.
						block_indexes[i_1].kL = j;
						bfound = true;
						break;
					}
				}
			}
			x4 = b[i].g.zE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ) || (b[i].g.iPlane == YZ))) {
				x4 = b[i].g.zC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.zC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.zC;
				}
			}
			bfound = false;
			/*
			for (j = 0; j <= inz; j++) {
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
				for (j = inz; j >= 0; j--) {
					if ((!bfound) && (fabs(x4 - zpos[j])<shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls))) {
						// Точное совпадение !!!
						block_indexes[i_1].kR = j;
						bfound = true;
						break;
					}
				}
			}
			if (!bfound) {
				for (j = inz; j >= 0; j--) {
					if ((!bfound) && (x4 < zpos[j])) {
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
				printf("enumerate_volume_improved_obobshenie function\n");
				printf("violation of the order block_indexes\n");
				printf("i=%lld iL=%lld iR=%lld\n", i, block_indexes[i_1].iL,
					block_indexes[i_1].iR);
				system("pause");
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
				printf("i=%lld jL=%lld jR=%lld\n", i, block_indexes[i_1].jL,
					block_indexes[i_1].jR);
				system("pause");
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
				printf("i=%lld kL=%lld kR=%lld\n", i, block_indexes[i_1].kL,
					block_indexes[i_1].kR);
				system("pause");
			}


			if ((block_indexes[i_1].kL >= 0) &&
				(block_indexes[i_1].kR >= 0) &&
				(block_indexes[i_1].kL + 1 == block_indexes[i_1].kR)) {
				// Всего одна клетка на блок.
				//iZ_one_CELL_count_statistic++;
			}

			if (0 && (i == 116)) {
				// debug
				printf("i=%lld iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
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
				printf("i=%lld iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
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




	// Количество проходов существенно сократилось и в итоге это приводит к существенному
	// увеличению быстродействия.
	integer m7= lb-1, m8;//  0 m7= lb - 1,
//#pragma omp parallel for
	//for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
		//integer iP = i1 + j1 * inx + k1 * inx*iny;
#pragma omp parallel for
	for (integer iP=0; iP<inx*iny*inz; iP++) {
		evt[iP] = -1;
		bvisit[iP] = false;
	}
   // for (m8 = 0; m8 < lb; m8++) {
	for (m8 = lb-1; m8 >= 0; m8--) {
		// Сначала идет блок с наивысшим приоритетом, мы не рассматриваем ячейки дважды.
		if (b[m8].g.itypegeom == PRISM) {
#pragma omp parallel for
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++) for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++) for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++) {
				integer iP = i1 + j1 * inx + k1 * inx*iny;

				if ((i1 < 0) || (i1 > inx) || (j1 < 0) || (j1 > iny) || (k1 < 0) || (k1 > inz)) {
					// ERROR
					printf("ERROR PRISM\n");
					printf("inx=%lld iny=%lld inz=%lld \n", inx, iny, inz);
					printf("i1=%lld j1=%lld k1=%lld \n", i1, j1, k1);
					printf("iP=%lld m8=%lld", iP, m8);
					printf("iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);
					system("PAUSE");
				}


				if (bvisit[iP] == false)
				{

					bvisit[iP] = true;

					switch (iflag) {
					case TEMPERATURE:
						if (b[m8].itype == HOLLOW) {
							evt[iP] = -1;
						}
						else {
							evt[iP] = m8;
						}
						break;
					case HYDRODINAMIC:
						if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
							evt[iP] = -1;
						}
						else {
							evt[iP] = m8;
						}
						break;
					}
				}
			}
			m7--;
			
		}
		else if (b[m8].g.itypegeom == CYLINDER) {

			// TODO как был сформирован призматический объект для цилиндра ? 
			// Надо также сократить число проверяемых точек.
			// Cylinder
			//for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
//#pragma omp parallel for
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++) for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++) for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++) {
				integer iP = i1 + j1 * inx + k1 * inx*iny;
				
				if ((i1 < 0) || (i1 > inx) || (j1 < 0) || (j1 > iny) || (k1 < 0) || (k1 > inz)) {
					// ERROR
					printf("ERROR CYLINDER\n");
					printf("iplane=%lld", b[m8].g.iPlane);
					printf("xC=%e yC=%e zC=%e Hcyl=%e\n", b[m8].g.xC, b[m8].g.yC, b[m8].g.zC, b[m8].g.Hcyl);
					printf("Rin=%e Rout=%e\n", b[m8].g.R_in_cyl, b[m8].g.R_out_cyl);
					printf("inx=%lld iny=%lld inz=%lld \n", inx, iny, inz);
					printf("i1=%lld j1=%lld k1=%lld \n", i1, j1, k1);
					printf("iP=%lld m8=%lld", iP, m8);
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
					case XY:
						if (fabs(b[m8].g.R_in_cyl) < 1.0e-40) {
							if ((p.z > b[m8].g.zC) && (p.z < b[m8].g.zC + b[m8].g.Hcyl)) {
								if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.yC - p.y)*(b[m8].g.yC - p.y)) < b[m8].g.R_out_cyl) {

									bvisit[iP] = true;

									switch (iflag) {
									case TEMPERATURE:
										if (b[m8].itype == HOLLOW) {
											evt[iP] = -1;
										}
										else {
											evt[iP] = m8;
										}
										break;
									case HYDRODINAMIC:
										if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
											if (b[m8].itype == HOLLOW) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										case HYDRODINAMIC:
											if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
					case XZ:
						if (fabs(b[m8].g.R_in_cyl) < 1.0e-40) {
							if ((p.y > b[m8].g.yC) && (p.y < b[m8].g.yC + b[m8].g.Hcyl)) {
								if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {

									bvisit[iP] = true;

									switch (iflag) {
									case TEMPERATURE:
										if (b[m8].itype == HOLLOW) {
											evt[iP] = -1;
										}
										else {
											evt[iP] = m8;
										}
										break;
									case HYDRODINAMIC:
										if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
											if (b[m8].itype == HOLLOW) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										case HYDRODINAMIC:
											if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
					case YZ:
						if (fabs(b[m8].g.R_in_cyl) < 1.0e-40) {
							if ((p.x > b[m8].g.xC) && (p.x < b[m8].g.xC + b[m8].g.Hcyl)) {
								if (sqrt((b[m8].g.yC - p.y)*(b[m8].g.yC - p.y) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {

									bvisit[iP] = true;

									switch (iflag) {
									case TEMPERATURE:
										if (b[m8].itype == HOLLOW) {
											evt[iP] = -1;
										}
										else {
											evt[iP] = m8;
										}
										break;
									case HYDRODINAMIC:
										if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
											if (b[m8].itype == HOLLOW) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										case HYDRODINAMIC:
											if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++) for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++) for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++) {

				integer iP = i1 + j1 * inx + k1 * inx*iny;

				if ((i1 < 0) || (i1 > inx) || (j1 < 0) || (j1 > iny) || (k1 < 0) || (k1 > inz)) {
					// ERROR
					printf("ERROR POLYGON\n");
					printf("inx=%lld iny=%lld inz=%lld \n", inx, iny, inz);
					printf("i1=%lld j1=%lld k1=%lld \n", i1, j1, k1);
					printf("iP=%lld m8=%lld", iP, m8);
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

					integer k74 = -1;
					if (in_polygon(p, b[m8].g.nsizei, b[m8].g.xi, b[m8].g.yi, b[m8].g.zi, b[m8].g.hi, b[m8].g.iPlane_obj2, k74, m8)) {
						//printf("i1=%d j1=%d k1=%d inx*iny*inz=%d\n",i1,j1,k1, inx*iny*inz);
						//printf("iL=%d iR=%d jL=%d jR=%d kL=%d kR=%d\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);

						bvisit[iP] = true;

						switch (iflag) {
						case TEMPERATURE:
							if (b[m8].itype == HOLLOW) {
								evt[iP] = -1;
							}
							else {
								evt[iP] = m8;
							}
							break;
						case HYDRODINAMIC:
							if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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

	//if (bvisit != NULL) {
	// Оператор delete может быть вызван повторно и даже к null указателю.
		delete[] bvisit;
		bvisit = NULL;
	//}
	

	printf("enumerate_volume_improved 80 procent.\n");
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	integer l = 1;
	integer ib;
	for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
		integer iP = i + j * inx + k * inx * iny;
		if (evt[iP] > -1) {
			ib = evt[iP]; // номер блока был сохранён ранее.
						 // Это очень нужно для записи репорта.
			whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
			tck_int_list[l - 1].i = i;
			tck_int_list[l - 1].j = j;
			tck_int_list[l - 1].k = k;
			evt[iP] = l;
			l++;
		}
		else
		{   // не принадлежит расчётной области
			evt[iP] = 0;
		}
	}	
	maxelm = l - 1;

	//if (block_indexes != NULL) {
	// Оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] block_indexes;
		block_indexes = NULL;
	//}

	printf("enumerate_volume_improved_obobshenie end.\n");

	/*
	bool inDomain = false;
	// нумерация контрольных объёмов начинается с единицы
	// если контрольный объём не принадлежит расчётной области то ставится 0.
	for (i = 0; i<inx; i++) for (j = 0; j<iny; j++) for (k = 0; k<inz; k++) {
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
void enumerate_volume(integer* &evt, integer &maxelm, integer iflag,
	doublereal* xpos, doublereal* ypos, doublereal* zpos, integer* &whot_is_block,
	integer inx, integer iny, integer inz, BLOCK* b, integer lb, 
	integer lu, UNION* &my_union, integer &iunion_id_p1, TOCKA_INT* &tck_int_list,
	WALL* &w, integer &lw, SOURCE* &s, integer &ls) {

	

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

			tck_int_list=new TOCKA_INT[inx*iny*inz];
			if (tck_int_list == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for evt tck_int_list...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}

			evt = NULL;
			evt = new integer[inx*iny*inz];
			if (evt == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for evt constr struct...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
			whot_is_block = NULL;
			whot_is_block = new integer[inx*iny*inz];
			if (whot_is_block == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for whot_is_block constr struct...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}
			TOCHKA p;
			// нумерация в evt начиная с единицы.
			// если не принадлежит расчётной области то стоит 0.
			integer l = 1, ib;
			integer i, j, k;
			bool inDomain = false;
			// нумерация контрольных объёмов начинается с единицы
			// если контрольный объём не принадлежит расчётной области то ставится 0.
			for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				p.x = 0.5*(xpos[i] + xpos[i + 1]);
				p.y = 0.5*(ypos[j] + ypos[j + 1]);
				p.z = 0.5*(zpos[k] + zpos[k + 1]);
				switch (iflag) {
				case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
				case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
				}
				if (iunion_id_p1 == 0) {
					for (integer iu = 0; iu < lu; iu++) {
						if ((p.x > my_union[iu].xS) && (p.x < my_union[iu].xE) &&
							(p.y > my_union[iu].yS) && (p.y < my_union[iu].yE) &&
							(p.z > my_union[iu].zS) && (p.z < my_union[iu].zE)) {
							// Если рассматривать кабинет с асемблесами, то область внутри
							// Асемблесов есть Hollow блок.
							inDomain = false;
						}
					}
				}
				else {
					if ((p.x > my_union[iunion_id_p1-1].xS) && (p.x < my_union[iunion_id_p1-1].xE) &&
						(p.y > my_union[iunion_id_p1-1].yS) && (p.y < my_union[iunion_id_p1-1].yE) &&
						(p.z > my_union[iunion_id_p1-1].zS) && (p.z < my_union[iunion_id_p1-1].zE)) {
						// Если рассматривать кабинет с асемблесами, то область внутри
						// Асемблесов есть Hollow блок.
						//inDomain = false;
						// Оставляем так всё верно.
						if (inDomain) {
						//	printf("in our assembles\n");
						//	system("PAUSE");
						}
					}
				}


				if (inDomain) {
					// принадлежит расчётной области
					evt[i + j*inx + k*inx*iny] = l;
					// Это очень нужно для записи репорта.
					whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
					tck_int_list[l - 1].i = i;
					tck_int_list[l - 1].j = j;
					tck_int_list[l - 1].k = k;
					l++;
				}
				else
				{   // не принадлежит расчётной области
					evt[i + j*inx + k*inx*iny] = 0;
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

				evt = NULL;
				evt = new integer[inx*iny*inz];
				if (evt == NULL) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem : not enough memory on your equipment for evt constr struct...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
				whot_is_block = NULL;
				whot_is_block = new integer[inx*iny*inz];
				if (whot_is_block == NULL) {
					// недостаточно памяти на данном оборудовании.
					printf("Problem : not enough memory on your equipment for whot_is_block constr struct...\n");
					printf("Please any key to exit...\n");
					exit(1);
				}
				TOCHKA p;
				// нумерация в evt начиная с единицы.
				// если не принадлежит расчётной области то стоит 0.
				integer l = 1, ib;
				integer i, j, k;
				bool inDomain = false;
				// нумерация контрольных объёмов начинается с единицы
				// если контрольный объём не принадлежит расчётной области то ставится 0.
				for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
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
						tck_int_list[l - 1].i = i;
						tck_int_list[l - 1].j = j;
						tck_int_list[l - 1].k = k;
						l++;
					}
					else
					{   // не принадлежит расчётной области
						evt[i + j*inx + k*inx*iny] = 0;
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
void init_evt_f_alice_improved(integer* &evt, integer iflag, doublereal* xpos, doublereal* ypos, doublereal* zpos,
	integer inx, integer iny, integer inz, BLOCK* b, integer lb, TOCKA_INT* &tck_int_list,
	WALL* &w, integer &lw, SOURCE* &s, integer &ls) {

	tck_int_list = new TOCKA_INT[inx*iny*inz];
	if (tck_int_list == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for evt tck_int_list...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	evt = NULL;
	evt = new integer[inx*iny*inz];
	if (evt == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for evt constr struct...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	//TOCHKA p;

	printf("enumerate_volume_improved start.\n");
	Block_indexes* block_indexes = new Block_indexes[lb];
	// оператор new не требует проверки на null.
	//if (block_indexes == NULL) {
		//printf("error in allocation memory for block_indexes in init_evt_f_alice_improved.\n");
		//system("pause");
		//exit(1);
	//}

	// 08.04.2018
	for (integer i = 0; i < lb; i++) {
		// инициализация, на случай если блоки не будут распознаны.
		block_indexes[i].iL = -1;
		block_indexes[i].iR = -2;
		block_indexes[i].jL = -1;
		block_indexes[i].jR = -2;
		block_indexes[i].kL = -1;
		block_indexes[i].kR = -2;
	}

	integer i=0, j=0, k=0;

	for (i = 0; i < lb; i++) {
		doublereal x4 = b[i].g.xS;
		for (j = 0; j <= inx; j++) {
			if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].iL = j;
				break;
			}
		}
		x4 = b[i].g.xE;
		for (j = 0; j <= inx; j++) {
			if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].iR = j;
				break;
			}
		}
		x4 = b[i].g.yS;
		for (j = 0; j <= iny; j++) {
			if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].jL = j;
				break;
			}
		}
		x4 = b[i].g.yE;
		for (j = 0; j <= iny; j++) {
			if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].jR = j;
				break;
			}
		}
		x4 = b[i].g.zS;
		for (j = 0; j <= inz; j++) {
			if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].kL = j;
				break;
			}
		}
		x4 = b[i].g.zE;
		for (j = 0; j <= inz; j++) {
			if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
				block_indexes[i].kR = j;
				break;
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
	integer m7;
	for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
		evt[i + j*inx + k*inx*iny] = -1;
	}
	for (m7 = 0; m7 < lb; m7++) {
		for (i = block_indexes[m7].iL; i < block_indexes[m7].iR; i++) for (j = block_indexes[m7].jL; j < block_indexes[m7].jR; j++) for (k = block_indexes[m7].kL; k < block_indexes[m7].kR; k++) {
			switch (iflag) {
			case TEMPERATURE:
				if (b[m7].itype == HOLLOW) {
					evt[i + j*inx + k*inx*iny] = -1;
				}
				else {
					evt[i + j*inx + k*inx*iny] = m7;
				}
				break;
			case HYDRODINAMIC:
				if ((b[m7].itype == SOLID) || (b[m7].itype == HOLLOW)) {
					evt[i + j*inx + k*inx*iny] = -1;
				}
				else {
					evt[i + j*inx + k*inx*iny] = m7;
				}
				break;
			}

		}
	}
	printf("enumerate_volume_improved 80 procent.\n");
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	integer l = 1, ib;
	for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
		if (evt[i + j*inx + k*inx*iny] > -1) {
			ib = evt[i + j*inx + k*inx*iny]; // номер блока был сохранён ранее.
											 // Это очень нужно для записи репорта.
											 //whot_is_block[l - 1] = ib; // номер блока которому принадлежит точка (p.x,p.y,p.z).
			evt[i + j*inx + k*inx*iny] = l;

			tck_int_list[l - 1].i = i;
			tck_int_list[l - 1].j = j;
			tck_int_list[l - 1].k = k;
			l++;
		}
		else {
			// не принадлежит расчётной области
			evt[i + j*inx + k*inx*iny] = 0;
		}
	}


	delete[] block_indexes;

	printf("enumerate_volume_improved end.\n");
	/*
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	integer l = 1, ib;
	integer i, j, k;
	bool inDomain = false;
	// нумерация контрольных объёмов начинается с единицы
	// если контрольный объём не принадлежит расчётной области то ставится 0.
	for (i = 0; i<inx; i++) for (j = 0; j<iny; j++) for (k = 0; k<inz; k++) {
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
void init_evt_f_alice_improved_obobshenie(integer* &evt, integer iflag, 
	doublereal* xpos, doublereal* ypos, doublereal* zpos,
	integer inx, integer iny, integer inz, BLOCK* b, integer lb, TOCKA_INT* &tck_int_list,
	WALL* &w, integer &lw, SOURCE* &s, integer &ls) {




	tck_int_list = new TOCKA_INT[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (tck_int_list == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for evt tck_int_list...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	bool* bvisit = NULL;
	bvisit = new bool[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (bvisit == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	evt = NULL;
	evt = new integer[inx*iny*inz];
	// оператор new не требует проверки на null.
	//if (evt == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for evt constr struct...\n");
		//printf("Please any key to exit...\n");
		//exit(1);
	//}

	//TOCHKA p;

	printf("enumerate_volume_improved start.\n");
	Block_indexes* block_indexes = new Block_indexes[lb];
	// оператор new не требует проверки на null.
	//if (block_indexes == NULL) {
		//printf("error in allocation memory for block_indexes in init_evt_f_alice_improved.\n");
		//system("pause");
		//exit(1);
	//}


	// 08.04.2018
	for (integer i = 0; i < lb; i++) {
		// инициализация, на случай если блоки не будут распознаны.
		block_indexes[i].iL = -1;
		block_indexes[i].iR = -2;
		block_indexes[i].jL = -1;
		block_indexes[i].jR = -2;
		block_indexes[i].kL = -1;
		block_indexes[i].kR = -2;
	}


	integer i=0, j=0, k=0, i_1 = lb-1;

	// Погрешность бывает абсолютная и относительная.
	// Вещественные числа в ЭВМ представляются с конечной точностью.
	// Лучше использовать относительную погрешность в 0.15%.
	//const doublereal otnositelnaq_tolerance_eps = 0.0015; // 0.15%

	for (i = lb-1; i >= 0; i--) {
		//if (b[i].g.itypegeom == PRISM) {

		// polygon (b[i].g.itypegeom == POLYGON)
		// Значения прямоугольной призмы xS, xE, yS, yE, zS, zE - хранят окаймляющую
		// полигон прямоугольную призму, что позволит проверять принадлежность точки полигону
		// только для ячеек сетки находящихся внутри данной прямоугольной призмы, что сильно 
		// ускоряет обработку.
		//if ((b[i].g.itypegeom == PRISM) || (b[i].g.itypegeom == CYLINDER) || (b[i].g.itypegeom == POLYGON))
		{

			doublereal x4 = b[i].g.xS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == XZ))) {
				x4 = b[i].g.xC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == YZ))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.xC;
				}
				else {
					x4 = b[i].g.xC + b[i].g.Hcyl;
				}
			}
			/*
			for (j = 0; j <= inx; j++) {
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
			for (j = 0; j <= inx; j++) {
				// Абсолютная погрешность.
				if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i_1].iL = j;
					break;
				}
			}
			x4 = b[i].g.xE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == XZ))) {
				x4 = b[i].g.xC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == YZ))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.xC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.xC;
				}
			}
			/*
			for (j = 0; j <= inx; j++) {
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
			for (j = 0; j <= inx; j++) {
				// Абсолютная погрешность.
				if (fabs(xpos[j] - x4) < shorter_length_for_simplificationX(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i_1].iR = j;
					break;
				}
			}
			x4 = b[i].g.yS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == YZ))) {
				x4 = b[i].g.yC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.yC;
				}
				else {
					x4 = b[i].g.yC + b[i].g.Hcyl;
				}
			}
			/*
			for (j = 0; j <= iny; j++) {
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
			for (j = 0; j <= iny; j++) {
				// Абсолютная погрешность.
				if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i_1].jL = j;
					break;
				}
			}
			x4 = b[i].g.yE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY) || (b[i].g.iPlane == YZ))) {
				x4 = b[i].g.yC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.yC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.yC;
				}
			}
			/*
			for (j = 0; j <= iny; j++) {

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
			for (j = 0; j <= iny; j++) {
				// Абсолютная погрешность.
				if (fabs(ypos[j] - x4) < shorter_length_for_simplificationY(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i_1].jR = j;
					break;
				}
			}
			x4 = b[i].g.zS;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ) || (b[i].g.iPlane == YZ))) {
				x4 = b[i].g.zC - b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.zC;
				}
				else {
					x4 = b[i].g.zC + b[i].g.Hcyl;
				}
			}
			/*
			for (j = 0; j <= inz; j++) {
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
			for (j = 0; j <= inz; j++) {
				// Абсолютная погрешность.
				if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i_1].kL = j;
					break;
				}
			}
			x4 = b[i].g.zE;
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XZ) || (b[i].g.iPlane == YZ))) {
				x4 = b[i].g.zC + b[i].g.R_out_cyl;
			}
			if ((b[i].g.itypegeom == CYLINDER) && ((b[i].g.iPlane == XY))) {
				if (b[i].g.Hcyl > 0.0) {
					x4 = b[i].g.zC + b[i].g.Hcyl;
				}
				else {
					x4 = b[i].g.zC;
				}
			}
			/*
			for (j = 0; j <= inz; j++) {
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
			for (j = 0; j <= inz; j++) {
				// Абсолютная погрешность.
				if (fabs(zpos[j] - x4) < shorter_length_for_simplificationZ(x4, b, lb, w, lw, s, ls)) {
					block_indexes[i_1].kR = j;
					break;
				}
			}

			if ((block_indexes[i_1].iL >= 0) &&
				(block_indexes[i_1].iR >= 0) &&
				(block_indexes[i_1].iL >= block_indexes[i_1].iR)) {
				printf("init_evt_f_alice_improved_obobshenie  function\n");
				printf("violation of the order block_indexes\n");
				printf("i=%lld iL=%lld iR=%lld\n", i, block_indexes[i_1].iL,
					block_indexes[i_1].iR);
				system("pause");
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
				printf("i=%lld jL=%lld jR=%lld\n", i, block_indexes[i_1].jL,
					block_indexes[i_1].jR);
				system("pause");
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
				printf("i=%lld kL=%lld kR=%lld\n", i, block_indexes[i_1].kL,
					block_indexes[i_1].kR);
				system("pause");
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
					printf("i=%lld iL=%lld iR=%lld jL=%lld jR=%lld kL=%lld kR=%lld\n", i,
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
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY) || (b[i_a].g.iPlane == XZ))) {
					x4 = b[i_a].g.xC - b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == YZ))) {
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
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY) || (b[i_a].g.iPlane == XZ))) {
					x4 = b[i_a].g.xC + b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == YZ))) {
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
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY) || (b[i_a].g.iPlane == YZ))) {
					x4 = b[i_a].g.yC - b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XZ))) {
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
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY) || (b[i_a].g.iPlane == YZ))) {
					x4 = b[i_a].g.yC + b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XZ))) {
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
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == YZ) || (b[i_a].g.iPlane == XZ))) {
					x4 = b[i_a].g.zC - b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY))) {
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
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == YZ) || (b[i_a].g.iPlane == XZ))) {
					x4 = b[i_a].g.zC + b[i_a].g.R_out_cyl;
				}
				if ((b[i_a].g.itypegeom == CYLINDER) && ((b[i_a].g.iPlane == XY))) {
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
	
#pragma omp parallel for
	for (integer iP = 0; iP<inx*iny*inz; iP++) {
		evt[iP] = -1;
		bvisit[iP] = false;
	}
	integer m7 = lb-1, m8;

	for (m8 = lb-1; m8 >= 0; m8--) {
		if (b[m8].g.itypegeom == PRISM) {
#pragma omp parallel for
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++) for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++) for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++) {
				integer iP = i1 + j1 * inx + k1 * inx*iny;

				if (bvisit[iP] == false)
				{

					bvisit[iP] = true;

					switch (iflag) {
					case TEMPERATURE:
						if (b[m8].itype == HOLLOW) {
							evt[iP] = -1;
						}
						else {
							evt[iP] = m8;
						}
						break;
					case HYDRODINAMIC:
						if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
							evt[iP] = -1;
						}
						else {
							evt[iP] = m8;
						}
						break;
					}

				}
			}
			m7--;
		}
		else if (b[m8].g.itypegeom == CYLINDER) {

			// TODO как был сформирован призматический объект для цилиндра ? 
			// Надо также сократить число проверяемых точек.
			// Cylinder
			//for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {

			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++) for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++) for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++) {
					integer iP = i1 + j1 * inx + k1 * inx*iny;

					if (bvisit[iP] == false)
					{

						TOCHKA p;
						p.x = 0.5*(xpos[i1] + xpos[i1 + 1]);
						p.y = 0.5*(ypos[j1] + ypos[j1 + 1]);
						p.z = 0.5*(zpos[k1] + zpos[k1 + 1]);

						switch (b[m8].g.iPlane) {
						case XY:
							if (fabs(b[m8].g.R_in_cyl) < 1.0e-40) {
								if ((p.z > b[m8].g.zC) && (p.z < b[m8].g.zC + b[m8].g.Hcyl)) {
									if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.yC - p.y)*(b[m8].g.yC - p.y)) < b[m8].g.R_out_cyl) {
									
										bvisit[iP] = true;

										switch (iflag) {
										case TEMPERATURE:
											if (b[m8].itype == HOLLOW) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										case HYDRODINAMIC:
											if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
												if (b[m8].itype == HOLLOW) {
													evt[iP] = -1;
												}
												else {
													evt[iP] = m8;
												}
												break;
											case HYDRODINAMIC:
												if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
						case XZ:
							if (fabs(b[m8].g.R_in_cyl) < 1.0e-40) {
								if ((p.y > b[m8].g.yC) && (p.y < b[m8].g.yC + b[m8].g.Hcyl)) {
									if (sqrt((b[m8].g.xC - p.x)*(b[m8].g.xC - p.x) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {
										
										bvisit[iP] = true;

										switch (iflag) {
										case TEMPERATURE:
											if (b[m8].itype == HOLLOW) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										case HYDRODINAMIC:
											if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
												if (b[m8].itype == HOLLOW) {
													evt[iP] = -1;
												}
												else {
													evt[iP] = m8;
												}
												break;
											case HYDRODINAMIC:
												if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
						case YZ:
							if (fabs(b[m8].g.R_in_cyl) < 1.0e-40) {
								if ((p.x > b[m8].g.xC) && (p.x < b[m8].g.xC + b[m8].g.Hcyl)) {
									if (sqrt((b[m8].g.yC - p.y)*(b[m8].g.yC - p.y) + (b[m8].g.zC - p.z)*(b[m8].g.zC - p.z)) < b[m8].g.R_out_cyl) {
										
										
										
										bvisit[iP] = true;

										switch (iflag) {
										case TEMPERATURE:
											if (b[m8].itype == HOLLOW) {
												evt[iP] = -1;
											}
											else {
												evt[iP] = m8;
											}
											break;
										case HYDRODINAMIC:
											if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
												if (b[m8].itype == HOLLOW) {
													evt[iP] = -1;
												}
												else {
													evt[iP] = m8;
												}
												break;
											case HYDRODINAMIC:
												if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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
			for (integer i1 = block_indexes[m7].iL; i1 < block_indexes[m7].iR; i1++) for (integer j1 = block_indexes[m7].jL; j1 < block_indexes[m7].jR; j1++) for (integer k1 = block_indexes[m7].kL; k1 < block_indexes[m7].kR; k1++) {

				integer iP = i1 + j1 * inx + k1 * inx*iny;

				if (bvisit[iP] == false)
				{

					//for (integer i1 = 0; i1 < inx; i1++) for (integer j1 = 0; j1 < iny; j1++) for (integer k1 = 0; k1 < inz; k1++) {
					TOCHKA p;
					p.x = 0.5*(xpos[i1] + xpos[i1 + 1]);
					p.y = 0.5*(ypos[j1] + ypos[j1 + 1]);
					p.z = 0.5*(zpos[k1] + zpos[k1 + 1]);

					integer k74 = -1;
					if (in_polygon(p, b[m8].g.nsizei, b[m8].g.xi, b[m8].g.yi, b[m8].g.zi, b[m8].g.hi, b[m8].g.iPlane_obj2, k74, m8)) {
						//printf("i1=%d j1=%d k1=%d inx*iny*inz=%d\n",i1,j1,k1, inx*iny*inz);
						//printf("iL=%d iR=%d jL=%d jR=%d kL=%d kR=%d\n", block_indexes[m7].iL, block_indexes[m7].iR, block_indexes[m7].jL, block_indexes[m7].jR, block_indexes[m7].kL, block_indexes[m7].kR);

						bvisit[iP] = true;

						switch (iflag) {
						case TEMPERATURE:
							if (b[m8].itype == HOLLOW) {
								evt[iP] = -1;
							}
							else {
								evt[iP] = m8;
							}
							break;
						case HYDRODINAMIC:
							if ((b[m8].itype == SOLID) || (b[m8].itype == HOLLOW)) {
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

	//if (bvisit != NULL) {
	// Оператор delete может быть вызван повторно и даже к null указателю.
		delete[] bvisit;
		bvisit = NULL;
	//}

	printf("enumerate_volume_improved 80 procent.\n");
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	integer l = 1, ib;
	for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
		integer iP = i + j * inx + k * inx*iny;
		if (evt[iP] > -1) {
			ib = evt[iP]; // номер блока был сохранён ранее.
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


	delete[] block_indexes;

	printf("enumerate_volume_improved end.\n");
	/*
	// нумерация в evt начиная с единицы.
	// если не принадлежит расчётной области то стоит 0.
	integer l = 1, ib;
	integer i, j, k;
	bool inDomain = false;
	// нумерация контрольных объёмов начинается с единицы
	// если контрольный объём не принадлежит расчётной области то ставится 0.
	for (i = 0; i<inx; i++) for (j = 0; j<iny; j++) for (k = 0; k<inz; k++) {
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
void init_evt_f_alice(integer* &evt,  integer iflag, doublereal* xpos, doublereal* ypos, doublereal* zpos, 
	integer inx, integer iny, integer inz, BLOCK* b, integer lb, TOCKA_INT* &tck_int_list,
	WALL* &w, integer &lw, SOURCE* &s, integer &ls) {


	bool b_only_Prism_object = true;
	for (integer m7 = 0; m7 < lb; m7++) {
		if (b[m7].g.itypegeom != PRISM) {
			b_only_Prism_object = false;
			break;
		}
	}


	if (1) {
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

			tck_int_list = new TOCKA_INT[inx*iny*inz];
			if (tck_int_list == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for evt tck_int_list...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}

			evt = NULL;
			evt = new integer[inx*iny*inz];
			if (evt == NULL) {
				// недостаточно памяти на данном оборудовании.
				printf("Problem : not enough memory on your equipment for evt constr struct...\n");
				printf("Please any key to exit...\n");
				exit(1);
			}

			TOCHKA p;
			// нумерация в evt начиная с единицы.
			// если не принадлежит расчётной области то стоит 0.
			integer l = 1, ib;
			integer i, j, k;
			bool inDomain = false;
			// нумерация контрольных объёмов начинается с единицы
			// если контрольный объём не принадлежит расчётной области то ставится 0.
			for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
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
					tck_int_list[l - 1].i = i;
					tck_int_list[l - 1].j = j;
					tck_int_list[l - 1].k = k;
					l++;
				}
				else
				{   // не принадлежит расчётной области
					evt[i + j*inx + k*inx*iny] = 0;
				}
			}
		}
	}

} // init_evt_f_alice.



// находит соседей для каждого внутреннего контрольного
// объёма или 0 если соседа нет.
void constr_sosed(integer* evt, integer* ent, integer** &sosed, integer maxelm,
				  integer inx, integer iny, integer inz, TOCKA_INT* &tck_int_list) {
    integer i;
	sosed=NULL;
	sosed = new integer*[12];
	if (sosed==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for sosed constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (i=0; i<12; i++) {
		sosed[i]=NULL;
	}
	for (i=0; i<12; i++) {
		sosed[i]=new integer[maxelm];
		if (sosed[i]==NULL) {
	       // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for sosed[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for sosed[%d] constr struct...\n", i);
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
   // for (i=0; i<(inx); i++) for (j=0; j<(iny); j++) for (k=0; k<(inz); k++) {
	//	ic=i+j*inx+k*inx*iny;
		//if (evt[ic]>0) {

	// Константа для ускорения проверок. 
	const integer iend_size = inx*iny*inz;
	const integer iend_size1 = inx*iny;

	for (integer iscan = 0; iscan < maxelm; iscan++) {
		{
			integer i = tck_int_list[iscan].i;
			integer j = tck_int_list[iscan].j;
			integer k = tck_int_list[iscan].k;
			ic = i + j * inx + k * inx*iny;
			


				if (l < maxelm) {

					//iE = (i + 1) + j*inx + k*inx*iny + 1;
					iE = ic + 2;
					if (((i + 1) < inx) && (iE > 0) && (iE <= iend_size) && (evt[iE - 1] > 0)) sosed[ESIDE][l] = evt[iE - 1];
					else {
						sosed[ESIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iW = (i - 1) + j*inx + k*inx*iny + 1;
					iW = ic;
					if (((i - 1) >= 0) && (iW > 0) && (iW <= iend_size) && (evt[iW - 1] > 0)) sosed[WSIDE][l] = evt[iW - 1];
					else {
						sosed[WSIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iN = i + (j + 1)*inx + k*inx*iny + 1;
					iN = ic + 1 + inx;
					if (((j + 1) < iny) && (iN > 0) && (iN <= iend_size) && (evt[iN - 1] > 0)) sosed[NSIDE][l] = evt[iN - 1];
					else {
						sosed[NSIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iS = i + (j - 1)*inx + k*inx*iny + 1;
					iS = ic - inx + 1;
					if (((j - 1) >= 0) && (iS > 0) && (iS <= iend_size) && (evt[iS - 1] > 0)) sosed[SSIDE][l] = evt[iS - 1];
					else {
						sosed[SSIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iT = i + j*inx + (k + 1)*inx*iny + 1;
					iT = ic + iend_size1 + 1;
					if (((k + 1) < inz) && (iT > 0) && (iT <= iend_size) && (evt[iT - 1] > 0)) sosed[TSIDE][l] = evt[iT - 1];
					else {
						sosed[TSIDE][l] = 0; // точка не принадлежит расчётной области
					}
					//iB = i + j*inx + (k - 1)*inx*iny + 1;
					iB = ic - iend_size1 + 1;
					if (((k - 1) >= 0) && (iB > 0) && (iB <= iend_size) && (evt[iB - 1] > 0)) sosed[BSIDE][l] = evt[iB - 1];
					else {
						sosed[BSIDE][l] = 0; // точка не принадлежит расчётной области
					}

					//iEE = (i + 2) + j*inx + k*inx*iny + 1;
					iEE = ic + 3;
					if (((i + 2) < inx) && (iEE > 0) && (iEE <= iend_size) && (evt[iEE - 1] > 0)) sosed[EE][l] = evt[iEE - 1];
					else sosed[EE][l] = 0; // точка не принадлежит расчётной области
					//iWW = (i - 2) + j*inx + k*inx*iny + 1;
					iWW = ic - 1;
					if (((i - 2) >= 0) && (iWW > 0) && (iWW <= iend_size) && (evt[iWW - 1] > 0)) sosed[WW][l] = evt[iWW - 1];
					else sosed[WW][l] = 0; // точка не принадлежит расчётной области
					//iNN = i + (j + 2)*inx + k*inx*iny + 1;
					iNN = ic + 2 * inx + 1;
					if (((j + 2) < iny) && (iNN > 0) && (iNN <= iend_size) && (evt[iNN - 1] > 0)) sosed[NN][l] = evt[iNN - 1];
					else sosed[NN][l] = 0; // точка не принадлежит расчётной области
					//iSS = i + (j - 2)*inx + k*inx*iny + 1;
					iSS = ic - 2 * inx + 1;
					if (((j - 2) >= 0) && (iSS > 0) && (iSS <= iend_size) && (evt[iSS - 1] > 0)) sosed[SS][l] = evt[iSS - 1];
					else sosed[SS][l] = 0; // точка не принадлежит расчётной области
					//iTT = i + j*inx + (k + 2)*inx*iny + 1;
					iTT = ic + 2 * iend_size1 + 1;
					if (((k + 2) < inz) && (iTT > 0) && (iTT <= iend_size) && (evt[iTT - 1] > 0)) sosed[TTSIDE][l] = evt[iTT - 1];
					else sosed[TTSIDE][l] = 0; // точка не принадлежит расчётной области
					//iBB = i + j*inx + (k - 2)*inx*iny + 1;
					iBB = ic - 2 * iend_size1 + 1;
					if (((k - 2) >= 0) && (iBB > 0) && (iBB <= iend_size) && (evt[iBB - 1] > 0)) sosed[BB][l] = evt[iBB - 1];
					else sosed[BB][l] = 0; // точка не принадлежит расчётной области

					l++;
				}
				else {
					printf("model error: l>maxelm\n");
					system("pause");
					exit(1);
				}
		}
	}

} // constr_sosed

// находит соседей для каждого внутреннего контрольного
// объёма или 0 если соседа нет.
void constr_sosed_flow(integer **evt_f2, integer iDom, 
					   integer** &sosed, integer maxelm,
				       integer inx, integer iny, integer inz) {
    integer i,j,k;	
	sosed=NULL;
	sosed = new integer*[12];
	if (sosed==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for sosed constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (i=0; i<12; i++) sosed[i]=NULL;
	for (i=0; i<12; i++) {
		sosed[i]=new integer[maxelm];
		if (sosed[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for sosed[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for sosed[%d] constr struct...\n", i);
#endif
		    
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}

	// Константа для ускорения проверок.
	const integer iend_size = inx*iny*inz;

	integer ic,l=0;
	integer iE, iN, iT, iW, iS, iB;
    integer iEE, iNN, iTT, iWW, iSS, iBB;
	// проход по всем контрольным объёмам принадлежащим
	// расчётной области.
    for (i=0; i<(inx); i++) for (j=0; j<(iny); j++) for (k=0; k<(inz); k++) {
		ic=i+j*inx+k*inx*iny;
			if ((evt_f2[ENUMERATECONTVOL][ic]>0) && (evt_f2[MASKDOMAINFLUID][ic] == (iDom + 1))) {


				
					if (l < maxelm) {
							//iE = (i + 1) + j*inx + k*inx*iny + 1;
							iE = ic + 2;
							if (((i + 1) < inx) && (iE > 0) && (iE <= iend_size) && (evt_f2[0][iE - 1] > 0) && (evt_f2[1][iE - 1] == (iDom + 1))) sosed[ESIDE][l] = evt_f2[0][iE - 1];
							else sosed[ESIDE][l] = 0; // точка не принадлежит расчётной области
							//iW = (i - 1) + j*inx + k*inx*iny + 1;
							iW = ic;
							if (((i - 1) >= 0) && (iW > 0) && (iW <= iend_size) && (evt_f2[0][iW - 1] > 0) && (evt_f2[1][iW - 1] == (iDom + 1))) sosed[WSIDE][l] = evt_f2[0][iW - 1];
							else sosed[WSIDE][l] = 0; // точка не принадлежит расчётной области
							//iN = i + (j + 1)*inx + k*inx*iny + 1;
							iN = ic + inx + 1;
							if (((j + 1) < iny) && (iN > 0) && (iN <= iend_size) && (evt_f2[0][iN - 1] > 0) && (evt_f2[1][iN - 1] == (iDom + 1))) sosed[NSIDE][l] = evt_f2[0][iN - 1];
							else sosed[NSIDE][l] = 0; // точка не принадлежит расчётной области
							//iS = i + (j - 1)*inx + k*inx*iny + 1;
							iS = ic - inx + 1;
							if (((j - 1) >= 0) && (iS > 0) && (iS <= iend_size) && (evt_f2[0][iS - 1] > 0) && (evt_f2[1][iS - 1] == (iDom + 1))) sosed[SSIDE][l] = evt_f2[0][iS - 1];
							else sosed[SSIDE][l] = 0; // точка не принадлежит расчётной области
							iT = i + j*inx + (k + 1)*inx*iny + 1;
							if (((k + 1) < inz) && (iT > 0) && (iT <= iend_size) && (evt_f2[0][iT - 1] > 0) && (evt_f2[1][iT - 1] == (iDom + 1))) sosed[TSIDE][l] = evt_f2[0][iT - 1];
							else sosed[TSIDE][l] = 0; // точка не принадлежит расчётной области
							iB = i + j*inx + (k - 1)*inx*iny + 1;
							if (((k - 1) >= 0) && (iB > 0) && (iB <= iend_size) && (evt_f2[0][iB - 1] > 0) && (evt_f2[1][iB - 1] == (iDom + 1))) sosed[BSIDE][l] = evt_f2[0][iB - 1];
							else sosed[BSIDE][l] = 0; // точка не принадлежит расчётной области

							//iEE = (i + 2) + j*inx + k*inx*iny + 1;
							iEE = ic + 3;
							if (((i + 2) < inx) && (iEE > 0) && (iEE <= iend_size) && (evt_f2[0][iEE - 1] > 0) && (evt_f2[1][iEE - 1] == (iDom + 1))) sosed[EE][l] = evt_f2[0][iEE - 1];
							else sosed[EE][l] = 0; // точка не принадлежит расчётной области
							//iWW = (i - 2) + j*inx + k*inx*iny + 1;
							iWW = ic - 1;
							if (((i - 2) >= 0) && (iWW > 0) && (iWW <= iend_size) && (evt_f2[0][iWW - 1] > 0) && (evt_f2[1][iWW - 1] == (iDom + 1))) sosed[WW][l] = evt_f2[0][iWW - 1];
							else sosed[WW][l] = 0; // точка не принадлежит расчётной области
							//iNN = i + (j + 2)*inx + k*inx*iny + 1;
							iNN = ic + 2 * inx + 1;
							if (((j + 2) < iny) && (iNN > 0) && (iNN <= iend_size) && (evt_f2[0][iNN - 1] > 0) && (evt_f2[1][iNN - 1] == (iDom + 1))) sosed[NN][l] = evt_f2[0][iNN - 1];
							else sosed[NN][l] = 0; // точка не принадлежит расчётной области
							//iSS = i + (j - 2)*inx + k*inx*iny + 1;
							iSS = ic - 2 * inx + 1;
							if (((j - 2) >= 0) && (iSS > 0) && (iSS <= iend_size) && (evt_f2[0][iSS - 1] > 0) && (evt_f2[1][iSS - 1] == (iDom + 1))) sosed[SS][l] = evt_f2[0][iSS - 1];
							else sosed[SS][l] = 0; // точка не принадлежит расчётной области
							iTT = i + j*inx + (k + 2)*inx*iny + 1;
							if (((k + 2) < inz) && (iTT > 0) && (iTT <= iend_size) && (evt_f2[0][iTT - 1] > 0) && (evt_f2[1][iTT - 1] == (iDom + 1))) sosed[TTSIDE][l] = evt_f2[0][iTT - 1];
							else sosed[TTSIDE][l] = 0; // точка не принадлежит расчётной области
							iBB = i + j*inx + (k - 2)*inx*iny + 1;
							if (((k - 2) >= 0) && (iBB > 0) && (iBB <= iend_size) && (evt_f2[0][iBB - 1] > 0) && (evt_f2[1][iBB - 1] == (iDom + 1))) sosed[BB][l] = evt_f2[0][iBB - 1];
							else sosed[BB][l] = 0; // точка не принадлежит расчётной области

							l++;
						}
						else {
							printf("model is incorrect : constr struct sosed index >= maxelm\n");
							system("pause");
							exit(1);
						}
		}
	}

} // constr_sosed_flow

// Всем соседям присваивается идентификатор центрального узла 
void my_fill_Domain(integer* &evt_f, integer i, integer j, integer k, 
					integer inx, integer iny, integer inz) {
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
void my_fill_Domain_recursive(integer* &evt_f, integer i, integer j, integer k, 
					integer inx, integer iny, integer inz) {
        
		integer iP, iE, iW, iN, iS, iT, iB;


		/* Рекурсивная процедура заливки распознаёт области произвольной сложности,
		а чтобы не было переполнения стека рекурсивных вызовов рекурсия полностью заменена
		на итерацию в цикле с помощью стека.
		*/

		// тот самый стек.
		typedef struct TSTEK {
	        integer i;
	        integer j;
	        integer k;
        } STEK;
		
		STEK *stek=NULL;
		stek=new STEK[inx*iny*inz];
		// оператор new не требует проверки на null.
		//if (stek==NULL) {
			//printf("no allocate memory for my_fill_Domain_recursive : stek==NULL\n");
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
			if (itop>0) {
				// Вынимаем вершину стека.
               i=stek[itop-1].i;
		       j=stek[itop-1].j;
		       k=stek[itop-1].k;
			   itop--;

	           iP=i+j*inx+k*inx*iny;
               iE=(i+1)+j*inx+k*inx*iny+1;
		       iW=(i-1)+j*inx+k*inx*iny+1;
               iN=i+(j+1)*inx+k*inx*iny+1;
		       iS=i+(j-1)*inx+k*inx*iny+1;
		       iT=i+j*inx+(k+1)*inx*iny+1;
		       iB=i+j*inx+(k-1)*inx*iny+1;



				if (evt_f[iP]>0) {

					// Узел принадлежит зоне FLUID.
					// Всем существующим соседям присвоить этот номер. 
           

					 // Узел принадлежит зоне FLUID.
					// Всем существующим соседям присвоить этот номер. 
					if (((i+1)<inx) && (iE>0) && (iE<=inx*iny*inz) && (evt_f[iE-1]>0)&&(evt_f[iE-1]!=evt_f[iP])) {
						evt_f[iE-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i+1, j, k, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i+1;
						stek[itop].j=j;
						stek[itop].k=k;
						itop++;
					}
					if (((i-1)>=0) && (iW>0) && (iW<=inx*iny*inz) && (evt_f[iW-1]>0)&&(evt_f[iW-1]!=evt_f[iP])) {
						evt_f[iW-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i-1, j, k, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i-1;
						stek[itop].j=j;
						stek[itop].k=k;
						itop++;			
					}
					if (((j+1)<iny) && (iN>0) && (iN<=inx*iny*inz) && (evt_f[iN-1]>0)&&(evt_f[iN-1]!=evt_f[iP])) {
						evt_f[iN-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i, j+1, k, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i;
	        			stek[itop].j=j+1;
						stek[itop].k=k;
						itop++;
					}
					if (((j-1)>=0) && (iS>0) && (iS<=inx*iny*inz) && (evt_f[iS-1]>0)&&(evt_f[iS-1]!=evt_f[iP])) {
						evt_f[iS-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i, j-1, k, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i;
						stek[itop].j=j-1;
						stek[itop].k=k;
						itop++;			
					}
					if (((k+1)<inz) && (iT>0) && (iT<=inx*iny*inz) && (evt_f[iT-1]>0)&&(evt_f[iT-1]!=evt_f[iP])) {
						evt_f[iT-1]=evt_f[iP];
						//my_fill_Domain_recursive(evt_f, i, j, k+1, inx, iny, inz);
						// PUSH вставить в стек.
						stek[itop].i=i;
						stek[itop].j=j;
						stek[itop].k=k+1;
						itop++;
					}
					if (((k-1)>=0) && (iB>0) && (iB<=inx*iny*inz) && (evt_f[iB-1]>0)&&(evt_f[iB-1]!=evt_f[iP])) {
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

		}

		//if (stek != NULL) {
		// Оператор delete может быть вызван повторно в том числе и к нулевому указателю.
			delete[] stek;
			stek = NULL;
		//}
} // my_fill_Domain_recursive


// создание связей гидродинамики с теплопроводностью.
// Первая часть отделена 22 сентября 2016 в связи с появлением АЛИС сетки.
// Эта часть универсальна и подходит и для АЛИС сетки тоже.
void constr_ptr_temp_part1(integer &flow_interior, 
	integer * &evt_f, integer** &evt_f2, integer* &domain_id, 
	integer inx, integer iny, integer inz, integer &icount_part
	) {

#ifdef _OPENMP
	int i_my_num_core_parallelesation = omp_get_num_threads();
	//omp_set_num_threads(8); // оптимально 8 потоков, 10 потоков уже проигрыш по времени.
	unsigned int nthreads = number_cores();
	omp_set_num_threads(nthreads); // установка числа потоков
#else 
	int i_my_num_core_parallelesation = 1;
#endif

	integer i = 0, j = 0, k = 0;

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

	if (brecursiv_algorithm) {
		// Рекурсивный алгоритм.
		// Рекурсия имитируется СТЕКОМ. Там теперь свой стековый алгоритм внутри.
		// Проблемы на геометрии 1ТТ113 26,03,2019.

		for (i = 0; i<inx; i++) for (j = 0; j<iny; j++) for (k = 0; k<inz; k++) {
			integer iP= i + j*inx + k*inx*iny;
			//if (bfirst_visit)
			if (bfirst_visit||(evt_f[iP]!= evt_f[iPgold])) {
				integer iE, iW, iN, iS, iT, iB;
				iE = (i + 1) + j*inx + k*inx*iny + 1;
				iW = (i - 1) + j*inx + k*inx*iny + 1;
				iN = i + (j + 1)*inx + k*inx*iny + 1;
				iS = i + (j - 1)*inx + k*inx*iny + 1;
				iT = i + j*inx + (k + 1)*inx*iny + 1;
				iB = i + j*inx + (k - 1)*inx*iny + 1;
				if ((iP >= 0) && (iP<inx*iny*inz) && (evt_f[iP]>0)) {
					if (((i + 1)<inx) && (iE>0) && (iE <= inx*iny*inz) && (evt_f[iE - 1]>0)) {
						if (((i - 1) >= 0) && (iW>0) && (iW <= inx*iny*inz) && (evt_f[iW - 1]>0)) {
							if (((j + 1)<iny) && (iN>0) && (iN <= inx*iny*inz) && (evt_f[iN - 1]>0)) {
								if (((j - 1) >= 0) && (iS>0) && (iS <= inx*iny*inz) && (evt_f[iS - 1]>0)) {
									if (((k + 1)<inz) && (iT>0) && (iT <= inx*iny*inz) && (evt_f[iT - 1]>0)) {
										if (((k - 1) >= 0) && (iB>0) && (iB <= inx*iny*inz) && (evt_f[iB - 1]>0)) {
#if doubleintprecision == 1
											// printf("in comming %lld\n",evt_f[iP]); // debug.
#else
											 //printf("in comming %d\n",evt_f[iP]); // debug.
#endif
											
											//system("PAUSE");
											
											// Это только первый вызов, вся работа делается внутри алгоритма my_fill_Domain_recursive.
											bfirst_visit = false;// досрочный выход из цикла for.
											// Там внутри рекурсия имитирована стеком. На модели из 576 блоков работает в разы быстрее.
											my_fill_Domain_recursive(evt_f, i, j, k, inx, iny, inz);
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
	else {

		bfirst_visit = false;
		bool *evt_etalon = NULL;
		evt_etalon = new bool[inx*iny*inz+2];
		if (evt_etalon == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for evt_etalon constr struct...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

		bool *evt_etalon1 = NULL;
		evt_etalon1 = new bool[inx*iny*inz+2];
		if (evt_etalon1 == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for evt_etalon1 constr struct...\n");
			printf("Please any key to exit...\n");
			exit(1);
		}

#pragma omp parallel for
		for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
			evt_etalon[i_1] = false; // initialisation.
			evt_etalon1[i_1] = false; // initialisation.
		}

		evt_etalon[0] = true;
#pragma omp parallel for
		for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
			if (evt_f[i_1] > 0) {
				evt_etalon[evt_f[i_1]] = true;
			}
		}
		

		printf("\n                                          ");
		printf("progress bar :\n	");
		printf("\\");
#ifdef _OPENMP
		omp_set_num_threads(1); // СБОЙ параллельной обработки 29_07_2017
#endif
		// гипотеза - образуются гонки разных нитей и они не могут залить одним цветом расчётную область.
		// 5 августа 2016. (рекурсивный алгоритм вылетает из-за переполнения стека помоему).
		// И переполнение стека есть серьёзная проблема.
		for (integer ilk = 0; ilk < 1000; ilk++) {
			bool bisequal = true;

			//if (ilk % 100 == 0) printf("progress bar %e procent\n", ilk / 10.0);
			printf("\b");
			printf("|");
			// Первый проход: Заливка.
			// 6*4 == 24 обхода
			// Меняем всевозможные направления обхода расчётной области,
			// чтобы справится со случаями "заглядывания за угол".
			// Также нужен прерыватель т.к. 1000 проходов это очень долго для tgf20.

#pragma omp parallel for
			for (i = 0; i < inx; i+=2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i+=2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("/");
			if (inx*iny*inz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}

#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("--");
			if (inx*iny*inz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b\b");
			printf("\\");
			if (inx*iny*inz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("|");
			if (inx*iny*inz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("/");
			if (inx*iny*inz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}

#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("--");
			if (inx*iny*inz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b\b");
			printf("\\");
			if (inx*iny*inz > 9000000) {
				// Очень большая модель, следовательно надо проверять досрочное прерывание.
				// Код досрочного прерывателя, если всё готово.
				// Иначе процесс сканирования может занять слишком много времени.
				// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
				evt_etalon1[0] = true;
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_f[i_1] > 0) {
						evt_etalon1[evt_f[i_1]] = true;
					}
				}

				bisequal = true;
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					if (evt_etalon[i_1] != evt_etalon1[i_1]) {
						bisequal = false;
						break;
					}
				}
				if (bisequal) {
					break;
				}
#pragma omp parallel for
				for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
					evt_etalon[i_1] = evt_etalon1[i_1];
					evt_etalon1[i_1] = false;// initialization.
				}
			}
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("|");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("/");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("--");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b\b");
			printf("\\");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("|");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("/");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("--");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b\b");
			printf("\\");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("|");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("/");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("--");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b\b");
			printf("\\");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("|");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("/");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b");
			printf("--");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			printf("\b\b");
			printf("\\");
#pragma omp parallel for
			for (i = 0; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
#pragma omp parallel for
			for (i = 1; i < inx; i += 2) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				my_fill_Domain(evt_f, i, j, k, inx, iny, inz);
			}
			//printf("\b");
			//printf("|");

			// Код досрочного прерывателя, если всё готово.
			// Иначе процесс сканирования может занять слишком много времени.
			// Недостаток в том что он требует 2*inx*iny*inz целых ячеек памяти.
			for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
				evt_etalon1[evt_f[i_1]]=true;
			}

			bisequal = true;
			for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
				if (evt_etalon[i_1] != evt_etalon1[i_1]) {
					bisequal = false;
				}
			}
			if (bisequal) {
				break;
			}
#pragma omp parallel for
			for (integer i_1 = 0; i_1 < inx*iny*inz; i_1++) {
				evt_etalon[i_1] = evt_etalon1[i_1];
				evt_etalon1[i_1] = false;// initialization.
			}

			//printf("ilk==%d\n",ilk);
		}
#ifdef _OPENMP
		omp_set_num_threads(8);// i_my_num_core_parallelesation
#endif
		printf("\n");
		//system("PAUSE");

		if (evt_etalon != NULL) {
			delete[] evt_etalon;
		}
		evt_etalon = NULL;

		if (evt_etalon1 != NULL) {
			delete[] evt_etalon1;
		}
		evt_etalon1 = NULL;
	}

	if (0&&bfirst_visit) {
		printf("ERROR!!!: NOT FILLING\n");
		system("PAUSE");
	}

	//printf("part 11\n");
#if doubleintprecision == 1
	printf("part %lld my_fill_Domain constr_ptr_temp_part1 \n", icount_part++); //22
#else
	printf("part %d my_fill_Domain constr_ptr_temp_part1 \n", icount_part++); //22
#endif

	integer iP=0;


	/*
	// debug GOOD
	// Проверка правильности распознавания по слоям.
	for (k=0; k<inz; k++) {
	    for (i=0; i<inx; i++) {
	        for (j=0; j<iny; j++) {
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
	integer* evt_f_shadow = NULL;
	evt_f_shadow = new integer[inx*iny*inz];
	// оператор new не требует проверки на null
	//if (evt_f_shadow == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for evt_f_shadow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// Копирование.
	//for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
		//iP = i + j*inx + k*inx*iny;
#pragma omp parallel for
	for (integer iP_loc=0; iP_loc < inx * iny*inz; iP_loc++) {
		evt_f_shadow[iP_loc] = evt_f[iP_loc];
	}

	// Второй проход:
	// Если возникает ошибка, то значение max_domain нужно увеличить.	
	const integer max_domain = 2048;// 256; // максимальное количество зон FLUID
	//integer *domain_id = NULL;
	domain_id = NULL;
	domain_id = new integer[max_domain];
	// оператор new не требует проверки на null.
	//if (domain_id == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for domain_id constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	for (i = 0; i<max_domain; i++) domain_id[i] = 0; // инициализация

	integer *domain_id_shadow = NULL;
	domain_id_shadow = new integer[max_domain];
	// оператор new не требует проверки на null.
	//if (domain_id_shadow == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for domain_id_shadow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	for (i = 0; i<max_domain; i++) domain_id_shadow[i] = 0; // инициализация


	integer id;
	bool bfind = false;
	bool bfind_shadow = false;
	//integer l;
	integer ic = 0; // счётчик по domain_id
	integer ic_shadow = 0; // счётчик по domain_id_shadow
	flow_interior = 0;
	// Если жидкие ячейки вообще были.
	if (!bfirst_visit) {
		for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
			iP = i + j * inx + k * inx * iny;
			if (evt_f[iP] > 0) {
				id = evt_f[iP]; // идентификатор связанной FLUID зоны.
				bfind = false;
				// Внимание!!! возможно это медленный участок кода!.
				for (integer l = 0; ((l < ic) && (l < max_domain)); l++) if (domain_id[l] == id) bfind = true;
				if (!bfind) {
					// Пробуем выполнить дозаливку. Добавка 29.04.2019.
					my_fill_Domain_recursive(evt_f, i, j, k, inx, iny, inz);

					printf("patch 29.04.2019. fluid zone id number = %lld. start zone control volume number = %lld \n", ic, id);
					if (ic >= max_domain - 1) {
						printf("error! nado uvelichit max_domain count...\n");
						printf("icount domain ==%lld\n", ic);
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
						printf("error! nado uvelichit max_domain count...\n");
						printf("icount domain ==%lld\n", ic);
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
	if (flow_interior > 1) {
		printf("WARNING : flow_interior count = %lld\n", flow_interior);
		printf("Your model contains multiple fluid zones.\n");
		//system("PAUSE");
	}
	//printf("part 12\n");
#if doubleintprecision == 1
	printf("part %lld my fill Domain recursive constr_ptr_temp_part1 \n", icount_part++); //22
#else
	printf("part %d  my fill Domain recursive constr_ptr_temp_part1 \n", icount_part++); //22
#endif

	if (1) {
		// В результате domain_id содержит единственную FLUID зону и её идентификатор 2.


		printf("5_08_2016: several fluid domain combine into one fluid domain.\n");
		// Смысл в том что даже если у нас несколько гидрдинамических подобластей то мы делаем одну единственную.
		//ic = 0;
		/*
		// Мутный код, логика ПЛОХАЯ. Исправлено 22 сентября 2016.
		for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
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
		// В результате domain_id содержит единственную FLUID зону и её идентификатор 2.
		// Если жидкие ячейки вообще были.
		if (!bfirst_visit) {
			//for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				//iP = i + j * inx + k * inx * iny;
#pragma omp parallel for
            for (integer iP_loc=0; iP_loc< inx * iny*inz; iP_loc++)
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

#if doubleintprecision == 1
		printf("part %lld evt_f init  constr_ptr_temp_part1 ", icount_part++); 
#else
		printf("part %d  evt_f init  constr_ptr_temp_part1 ", icount_part++); 
#endif
		printf(" 05_08_2016; 22_09_2016.\n");
	}



	//if (domain_id != NULL) {
	// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] domain_id;
		domain_id = NULL;
	//}

	//if (domain_id_shadow != NULL) {
	// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] domain_id_shadow;
		domain_id_shadow = NULL;
	//}

	// Освобождение памяти.
	//if (evt_f_shadow != NULL) {
	// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] evt_f_shadow;
	//}
	evt_f_shadow = NULL;

#ifdef _OPENMP
	omp_set_num_threads(i_my_num_core_parallelesation);
#endif
} // constr_ptr_temp_part1

//  откладываем вычисление evt_f2 как можно дальше т.к. оно не проходит 
// по памяти на большеразмерных  моделях.
// Здесь содержится код вычисляющий только evt_f2.
void constr_ptr_temp_part2(integer &flow_interior,
	integer * &evt_f, integer** &evt_f2, integer* &domain_id,
	integer inx, integer iny, integer inz, 
	TOCKA_INT* &tck_int_list, integer &maxelm, integer &icount_part) {

	integer i = 0, j = 0, k = 0;

	// Нам понадобтся знать информацию о различных подобластях
	// для процедуры mass_balance.
	// Т.к. процедуру балансировки масс надо делать для каждой 
	// отдельной связанной гидродинамической
	// подобласти иначе не будет сходимости.
	// Чтобы распознать цвет связаной подобласти он будет хранится в evt_f2[2][iP] позиции. 
	integer* evt_f_shadow = NULL;
	evt_f_shadow = new integer[inx*iny*inz];
	// Результат работы оператора new не требует проверки на null.
	//if (evt_f_shadow == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for evt_f_shadow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// Копирование.
	//for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
		//integer iP = i + j*inx + k*inx*iny;
		//evt_f_shadow[iP] = evt_f[iP];
	//}
	// Более быстрое копирование.
	const integer isize = inx * iny*inz;
	for (integer iP = 0; iP < isize; iP++) {
		evt_f_shadow[iP] = evt_f[iP];
	}

	// Второй проход:
	// Если возникает ошибка, то значение max_domain нужно увеличить.
	integer max_domain = 256; // максимальное количество зон FLUID
						  //integer *domain_id = NULL;
	domain_id = NULL;
	domain_id = new integer[max_domain];
	// оператор new не требует после себя проверки на null
	//if (domain_id == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for domain_id constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	for (i = 0; i<max_domain; i++) domain_id[i] = 0; // инициализация

	integer *domain_id_shadow = NULL;
	domain_id_shadow = new integer[max_domain];
	// оператор new не требует после себя проверки на null
	//if (domain_id_shadow == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for domain_id_shadow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	for (i = 0; i<max_domain; i++) domain_id_shadow[i] = 0; // инициализация


	integer id;
	bool bfind;
	bool bfind_shadow;
	integer l;
	integer ic = 0; // счётчик по domain_id
	integer ic_shadow = 0; // счётчик по domain_id_shadow
	flow_interior = 0;
	//for (i = 0; i<inx; i++) for (j = 0; j<iny; j++) for (k = 0; k<inz; k++) {
		//integer iP = i + j*inx + k*inx*iny;
		//if (evt_f[iP]>0) {
	for (integer iscan = 0; iscan<maxelm; iscan++) 
		{
			integer i = tck_int_list[iscan].i;
			integer j = tck_int_list[iscan].j;
			integer k = tck_int_list[iscan].k;
			integer iP = i + j * inx + k * inx*iny;
			{

			id = evt_f[iP]; // идентификатор связанной FLUID зоны.
			bfind = false;
			// Внимание!!! возможно это медленный участок кода!.
			for (l = 0; l<max_domain; l++) if (domain_id[l] == id) bfind = true;
			if (!bfind) {
				if (ic >= max_domain - 1) {
					printf("error! nado uvelichit max_domain count...\n");
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
			for (l = 0; l<max_domain; l++) if (domain_id_shadow[l] == id) bfind_shadow = true;
			if (!bfind_shadow) {
				if (ic_shadow >= max_domain - 1) {
					printf("error! nado uvelichit max_domain count...\n");
					//system("PAUSE");
					system("PAUSE");
					exit(1);
				}
				domain_id_shadow[ic_shadow++] = id;
				//flow_interior = ic;
			}
		}
	}
	//printf("part 12.2\n");
#if doubleintprecision == 1
	printf("part %lld constr domain_id constr_ptr_temp_part2 \n", icount_part++);
#else
	printf("part %d constr domain_id constr_ptr_temp_part2 \n", icount_part++);
#endif

#if doubleintprecision == 1
	// printf("flow interior=%lld\n",flow_interior); //debug.
#else
	// printf("flow interior=%d\n",flow_interior); //debug.
#endif
	
	// system("PAUSE");

	integer i1; // счётчик

			// Выделение оперативной памяти
	evt_f2 = NULL;
	evt_f2 = new integer*[3];
	if (evt_f2 == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for evt_f2 constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (i = 0; i<3; i++) evt_f2[i] = NULL;
	for (i = 0; i<3; i++) {
		evt_f2[i] = new integer[inx*iny*inz];
		if (evt_f2[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for evt_f2[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for evt_f2[%d] constr struct...\n", i);
#endif
				printf("Please any key to exit...\n");
			// system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// Проход по всем контрольным объёмам
	integer ipos = -1;
	integer ipos_shadow = -1;

	for (integer iP = 0; iP < isize; iP++) {
		evt_f2[MASKDOMAINFLUID][iP] = 0; // не принадлежит расчётной области.
	}
	//for (i = 0; i<inx; i++) for (j = 0; j<iny; j++) for (k = 0; k<inz; k++) {
		//integer iP = i + j*inx + k*inx*iny;

	for (integer iscan = 0; iscan<maxelm; iscan++) 
		{
			integer i = tck_int_list[iscan].i;
			integer j = tck_int_list[iscan].j;
			integer k = tck_int_list[iscan].k;
			integer iP = i + j * inx + k * inx*iny;
			
		//if (evt_f[iP] == 0) {
			//evt_f2[MASKDOMAINFLUID][iP] = 0; // не принадлежит расчётной области.
		//}
		//else
			{
			// FLUID Domain
			// 22.09.2016 Он может быть только равен 2.
			integer id = evt_f[iP];
			// Внимание !!! осторожно возможно медленный участок кода.
			for (i1 = 0; i1<max_domain; i1++) if (domain_id[i1] == id) {
				// 22.09.2016 Он обязатенльно найдётся и его номер обязательно будет ipos==i1==0.
				ipos = i1;
				break;
			}
			// на единицу больше чем уникальный идентификатор FLUID
			// Здесь всегда присваивается лишь единица.
			evt_f2[MASKDOMAINFLUID][iP] = ipos + 1; // Нумерация FLUID Domain начинается с 1.

													// FLUID Domain
			integer id_shadow = evt_f_shadow[iP];
			// Внимание !!! осторожно возможно медленный участок кода.
			for (i1 = 0; i1<max_domain; i1++) if (domain_id_shadow[i1] == id_shadow) {
				// Присваивается номер связанной FLUID зоны начиная с нуля.
				ipos_shadow = i1;
				break;
			}
			// на единицу больше чем уникальный идентификатор FLUID
			evt_f2[MASKDOMAINFLUIDCOLOR][iP] = ipos_shadow + 1; // Нумерация FLUID Domain начинается с 1.
		}

	}

	//printf("part 13\n");
#if doubleintprecision == 1
	printf("part %lld constr evt_f2 constr_ptr_temp_part2 \n", icount_part++);
#else
	printf("part %d  constr evt_f2 constr_ptr_temp_part2 \n", icount_part++);
#endif

	/*
	// debug GOOD
	// Проверка правильности распознавания по слоям.
	for (k=0; k<inz; k++) {
	for (i=0; i<inx; i++) {
	for (j=0; j<iny; j++) {
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

	//if (domain_id_shadow != NULL) {
	// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] domain_id_shadow;
		domain_id_shadow = NULL;
	//}

		// Здесь нельзя освобождать озу из под domain_id.
		// domain_id используется дальше по коду.
		//delete[] domain_id;
		//domain_id = NULL;

	// Освобождение памяти.
	//if (evt_f_shadow != NULL) {
	// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
		delete[] evt_f_shadow;
	//}
	evt_f_shadow = NULL;

} // constr_ptr_temp_part2

// Выделение памяти под FLOW.
// Дата создания 22 сентября 2016.
// Данная функция лишь несёт нагрузку выделения оперативной памяти, а так из неё всё выгружено при
// АЛИС сетке.
void constr_ptr_temp_allocation_memory_alice(integer &flow_interior, FLOW* &f, integer &icount_part)
{
	f = new FLOW[flow_interior];
	// проход по всем жидким зонам.
	for (integer l = 0; l < flow_interior; l++) {
		// Число ненулевых КО.
		// Заглушка которую возможно можно будет оживить в дальнейшем.
		// Самое простое это не использовать информацию whot_is_block для гидродинамических подобластей.
		f[l].whot_is_block = NULL;
		// Количество элементов в модели будет подсчитано в другом месте в специальной функции.
		//f[l].maxelm = domain_counter[l];
#if doubleintprecision == 1
		//printf("FLUID%lld maxelm=%lld\n", l, f[l].maxelm);
#else
		//printf("FLUID%d maxelm=%d\n", l, f[l].maxelm);
#endif
		
		// Выделение оперативной памяти.
		f[l].ptr = NULL;
		// f[l].ptr[fluid_elm_id]=temper_elm_id; // Смысл ptr.
		// Все закоментированное заполняется в другом месте !!!.
		/*
		f[l].ptr = new integer[f[l].maxelm];
		#if doubleintprecision == 1
			printf("delete flow %lld ptr\n", l);
		#else
			printf("delete flow %d ptr\n", l);
		#endif
		
		if (f[l].ptr == NULL) { // -2N
			// недостаточно памяти на данном оборудовании.
			#if doubleintprecision == 1
				printf("Problem : not enough memory on your equipment for ptr flow%lld constr struct...\n", l);
			#else
				printf("Problem : not enough memory on your equipment for ptr flow%d constr struct...\n", l);
			#endif
			
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
		for (j = 0; j < maxelm_t; j++) {
			if (ptr[MASKDOMAINFLUID][j] == l) {
				f[l].ptr[ptr[ENUMERATECONTVOL][j]] = j;
			}
		}
		*/
	}
	//printf("part 14\n");
	//printf("part 15\n");
	//printf("part 16\n");
#if doubleintprecision == 1
	printf("part %lld constr_ptr_temp_allocation_memory_alice\n", icount_part++);
#else
	printf("part %d constr_ptr_temp_allocation_memory_alice\n", icount_part++);
#endif
#if doubleintprecision == 1
	printf("part %lld constr_ptr_temp_allocation_memory_alice\n", icount_part++);
#else
	printf("part %d constr_ptr_temp_allocation_memory_alice\n", icount_part++);
#endif
#if doubleintprecision == 1
	printf("part %lld constr_ptr_temp_allocation_memory_alice\n", icount_part++);
#else
	printf("part %d constr_ptr_temp_allocation_memory_alice\n", icount_part++);
#endif
} // constr_ptr_temp_allocation_memory_alice


// создание связей гидродинамики с теплопроводностью.
void constr_ptr_temp(integer &flow_interior, FLOW* &f, integer maxelm_t,
					 integer** &ptr, integer *evt_t, integer * &evt_f, integer** &evt_f2,
					 integer* &domain_id, integer inx, integer iny, integer inz, 
	                 bool breconstruct, integer& icount_part) {
	
    

	
	integer iP;

	// Если возникает ошибка, то значение max_domain нужно увеличить.
	integer max_domain = 256; // максимальное количество зон FLUID

	integer id;
	//bool bfind;
	//bool bfind_shadow;
	integer l;
	integer ic = 0; // счётчик по domain_id
	integer ic_shadow = 0; // счётчик по domain_id_shadow
	//flow_interior = 0; // Заполняется в constr_ptr_temp_part1.
	
	

	 integer i1; // счётчик
	 	
	 integer ipos = -1;
	 integer ipos_shadow = -1;

	// 22 сентября 2016.
	// До части 13 можно использовать и при АЛИС сетке.

	// Выделение памяти :

	 //printf("delete temperature ptr\n");
	 if (ptr != NULL) {
		 for (unsigned char i = 0; i<2; i++) {
			 if (ptr[i] != NULL) {
				 delete[] ptr[i];
				 ptr[i] = NULL;
			 }
		 }
		 delete[] ptr;
		 ptr = NULL;
	 }	 

	ptr=NULL;
    ptr = new integer*[2];
	if (ptr==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for ptr constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (unsigned char  i=0; i<2; i++) ptr[i]=NULL;

	for (unsigned char  i=0; i<2; i++) {
		ptr[i]=new integer[maxelm_t];
		// оператор new не требует после себя проверки на null
		//if (ptr[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
//#if doubleintprecision == 1
	//		printf("Problem : not enough memory on your equipment for ptr[%d] constr struct...\n", i);
//#else
	//		printf("Problem : not enough memory on your equipment for ptr[%d] constr struct...\n", i);
//#endif
		    
	//	    printf("Please any key to exit...\n");
			//system("PAUSE");
		//	system("pause");
		  //  exit(1);
	    //}   
	}

	for (unsigned char i = 0; i < 2; i++) {
		for (integer i87 = 0; i87 < maxelm_t; i87++) {
			ptr[i][i87] = -1; // инициализация твердым телом.
		}
	}

	integer i = 0, j = 0, k = 0;


    //#define ENUMERATECONTVOL 0 // слой нумерации КО
    //#define MASKDOMAINFLUID 1 // слой маска определяет зону FLUID
    //#define MASKDOMAINFLUIDCOLOR 2 // В evt_f2 цвет гидродинамической подобласти.

	if ((ptr != NULL) && (ptr[0] != NULL) && (ptr[1] != NULL)) {

		// Третий проход. Маркировка отдельных зон FLUID:
		l = 0;
		ipos = -1;
		for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
			iP = i + j*inx + k*inx*iny;
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
						printf("perepolnenie ptr[MASKDOMAINFLUID][l] l>=maxelm_t\n");
						system("pause");
						exit(1);
					}
				}
				else {
					// не fluid но солид (и не холлоу)
					if (ptr != NULL) {
						if (ptr[MASKDOMAINFLUID] != NULL) {
							if (l < maxelm_t) {
								ptr[MASKDOMAINFLUID][l] = -1;
							}
							else {
								printf("perepolnenie 2 ptr[MASKDOMAINFLUID][l] l>=maxelm_t\n");
								system("pause");
								exit(1);
							}
						}
						else {
							printf("error : null ptr[MASKDOMAINFLUID] in constr_ptr_temp.\n");
							system("pause");
							exit(1);
						}
					}
					else {
						//ptr == NULL;
						printf("error : null ptr in constr_ptr_temp.\n");
						system("pause");
						exit(1);
					}
				}

				l++;
			}
		}

		//printf("part 14\n");
#if doubleintprecision == 1
		printf("part %lld construct ptr[MASKDOMAINFLUID][] constr_ptr_temp \n", icount_part++);
#else
		printf("part %d  construct ptr[MASKDOMAINFLUID][] constr_ptr_temp \n", icount_part++);
#endif

		// Счётчик внутренних КО в каждой из domain FLUID.
		integer *domain_counter = NULL;
		domain_counter = new integer[max_domain];
		// результат применения оператора new не надо проверять на null.
		//if (domain_counter == NULL) {
			// недостаточно памяти на данном оборудовании.
			//printf("Problem : not enough memory on your equipment for domain_counter constr struct...\n");
			//printf("Please any key to exit...\n");
			//system("PAUSE");
			//system("pause");
			//exit(1);
		//}
		//else 
		{

			for (i = 0; i < max_domain; i++) domain_counter[i] = 0; // инициализация

			if (ptr != NULL) {
				if ((ptr[MASKDOMAINFLUID] != NULL)&&(ptr[ENUMERATECONTVOL]!=NULL)) {
					for (i = 0; i < maxelm_t; i++) {
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
					printf("memory for ptr[MASKDOMAINFLUID] no allocate. error.\n");
					printf("see function constr_ptr_temp in module constr struct\n");
					system("pause");
					exit(1);
				}
			}
			else {
				printf("memory for ptr no allocate. error.\n");
				printf("see function constr_ptr_temp in module constr struct\n");
				system("pause");
				exit(1);
			}


			if (!breconstruct) {
				if (flow_interior > 0) {
					f = new FLOW[flow_interior];
					for (l = 0; l < flow_interior; l++) {
						f[l].ptr = NULL;
						f[l].nvtx = NULL;
						f[l].pa = NULL;
						f[l].sosedi = NULL;
						f[l].sosedb = NULL;
						f[l].whot_is_block = NULL;
						f[l].potent = NULL;
						f[l].prop = NULL;
						f[l].prop_b = NULL;
						f[l].alpha = NULL;
						f[l].slau = NULL;
						f[l].slau_bon = NULL;
						f[l].diag_coef = NULL;
						f[l].iN = NULL;
						f[l].id = NULL;
						f[l].rdistWall = NULL;
						f[l].SInvariantStrainRateTensor = NULL;
						f[l].mf = NULL;
						f[l].icolor_different_fluid_domain = NULL;
						f[l].ifrontregulationgl = NULL;
						f[l].ibackregulationgl = NULL;
					}
				}
				else {
					f = NULL;
				}
			}
			if (f != NULL) {
				// проход по всем жидким зонам.
				for (l = 0; l < flow_interior; l++) {
					// Число ненулевых КО.
					// Заглушка которую возможно можно будет оживить в дальнейшем.
					// Самое простое это не использовать информацию whot_is_block для гидродинамических подобластей.
					f[l].whot_is_block = NULL;
					f[l].maxelm = domain_counter[l];
#if doubleintprecision == 1
					printf("FLUID%lld maxelm=%lld\n", l, f[l].maxelm);
#else
					printf("FLUID%d maxelm=%d\n", l, f[l].maxelm);
#endif
					
					// Выделение оперативной памяти.
					if (f[l].ptr != NULL) {
						delete[] f[l].ptr;
					}
					f[l].ptr = NULL;
					f[l].ptr = new integer[f[l].maxelm];
#if doubleintprecision == 1
					printf("delete flow %lld ptr\n", l);
#else
					printf("delete flow %d ptr\n", l);
#endif
					
					if (f[l].ptr == NULL) { // -2N
						// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
						printf("Problem : not enough memory on your equipment for ptr flow%lld constr struct...\n", l);
#else
						printf("Problem : not enough memory on your equipment for ptr flow%d constr struct...\n", l);
#endif
						printf("Please any key to exit...\n");
						//system("PAUSE");
						system("pause");
						exit(1);
					}
					for (j = 0; j < maxelm_t; j++) {
						if (ptr[MASKDOMAINFLUID][j] == l) {
							f[l].ptr[ptr[ENUMERATECONTVOL][j]] = j;
						}
					}
				}
			}

			//printf("part 15\n");
#if doubleintprecision == 1
			printf("part %lld construct ptr[ENUMERATECONTVOL][] constr_ptr_temp \n", icount_part++);
#else
			printf("part %d construct ptr[ENUMERATECONTVOL][] constr_ptr_temp \n", icount_part++);
#endif

			// Счётчик внутренних КО на дополнительной структуре.
			for (i = 0; i < max_domain; i++) domain_counter[i] = 0; // инициализация

			for (i = 0; i < inx; i++) for (j = 0; j < iny; j++) for (k = 0; k < inz; k++) {
				iP = i + j*inx + k*inx*iny;
				if (evt_f2[MASKDOMAINFLUID][iP] == 0) {
					evt_f2[ENUMERATECONTVOL][iP] = 0; // не принадлежит расчётной области
				}
				else {
					// Уникальная нумерация узлов гидродинамической подобласти.
					// Нумерация начинается с единицы.
					evt_f2[ENUMERATECONTVOL][iP] = domain_counter[evt_f2[MASKDOMAINFLUID][iP] - 1] + 1;
					domain_counter[evt_f2[MASKDOMAINFLUID][iP] - 1] = domain_counter[evt_f2[MASKDOMAINFLUID][iP] - 1] + 1;
				}
			}


			/*
			// debug GOOD
			// Проверка правильности распознавания по слоям.
			for (k=0; k<inz; k++) {
				for (i=0; i<inx; i++) {
					for (j=0; j<iny; j++) {
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


			//if (domain_counter != NULL) {
			// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
				delete[] domain_counter;
				domain_counter = NULL;
			//}


		}

		//if (domain_id != NULL) {
		// оператор delete может быть вызван повторно в том числе и к нулевому указателю.
			delete[] domain_id;
			domain_id = NULL;
		//}
		

		//printf("part 16\n");
#if doubleintprecision == 1
			printf("part %lld construct evt_f2 constr_ptr_temp \n", icount_part++);
#else
			printf("part %d  construct evt_f2  constr_ptr_temp \n", icount_part++);
#endif

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
	TOCHKA *temp=NULL;
	temp = new TOCHKA[isize];
	// оператор new не требует после себя проверки на null
	//if (temp==NULL) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for SetLength_pointeger constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
    integer i=0;
	for (i=0; i<isize; i++) {
		temp[i].x=0.0; // инициализация
		temp[i].y=0.0;
		temp[i].z=0.0;
	}
    /*
	integer isizeold;
	if (ra != NULL) {
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
	for (i=0; i<isizeold; i++) temp[i]=ra[i];
	
	//if (ra != NULL) {
		// Оператор delete может быть применен повторно в том числе и к нулевому указателю.
		delete[]  ra; // уничтожение объекта
		ra = NULL;
	//}
	ra = new TOCHKA[isize]; // выделение памяти
	for (i=0; i<isize; i++) ra[i]=temp[i]; // копирование
	//if (temp != NULL) {
		// Оператор delete может быть применен повторно в том числе и к нулевому указателю.
		delete[] temp; // освобождение памяти
		temp = NULL;
	//}
	
} // SetLength_point

// Содержался линейный поиск. Внимание медленная версия.
// добавляет уникальную точку к массиву pa
// длины il
void addpoint_old(TOCHKA* &pa,integer &maxnode,TOCHKA pnew, integer* &ent, integer node, bool* &bvisit) {
	//integer i=0;
	bool bfind=false;
	/* // медленный участок кода
	for (i=0; i<maxnode; i++) { 
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
		for (integer i = 0; i<maxnode; i++) {
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
void addpoint(TOCHKA* &pa, integer &maxnode, TOCHKA pnew, integer* &ent, integer node, bool* &bvisit) {

	if (!bvisit[node]) {
		// добавление, точка ещё не встречалась
		
		pa[maxnode++] = pnew;
		ent[node] = maxnode; // нумерация начинается с единицы
		bvisit[node] = true; // индексация хеша.
		
	}
	
} // addpoint


// создаёт массив узлов принадлежащих расчётной области
void constr_nodes(TOCHKA* &pa, integer &maxnode, integer* &ent, integer iflag, integer* &whot_is_block, integer* &evt,
	integer inx, integer iny, integer inz,
	doublereal *xpos, doublereal *ypos, doublereal *zpos, BLOCK* b, integer lb, TOCKA_INT* &tck_int_list, integer &maxelm) {
	// iflag - принимает два значения : TEMPERATURE или HYDRODINAMIC и 
	// указывает с какие уравнения будут решаться и в какой расчётной области.

	// ent - глобальная нумерация узлов.
	// нумерация узлов начинается с единицы, если узел не принадлежит
	// расчётной области то ставится значение 0.
	const integer iSIZE = (inx + 1)*(iny + 1)*(inz + 1);
	ent = NULL;
	ent = new integer[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (ent == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for ent constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// был ли добавлен узел с номером node
	bool *bvisit = NULL;
	bvisit = new bool[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (bvisit == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// Выделение памяти. Сразу много для наискорейшего выполнения:
	pa = NULL;
	pa = new TOCHKA[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (pa == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for pa constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}


	TOCHKA p, pnew;
#pragma omp parallel for
	for (integer i = 0; i<iSIZE; i++) {
		ent[i] = 0;
		bvisit[i] = false;
	}

	bool inDomain = false;
	integer ib, node;
	maxnode = 0;// maxnode - текущая длина массива pa
	// Цикл по всем контрольным объёмам
	//for (i = 0; i<inx; i++) for (j = 0; j<iny; j++) for (k = 0; k<inz; k++) {
	for (integer iscan = 0; iscan<maxelm; iscan++) {

			integer i = tck_int_list[iscan].i;
			integer j = tck_int_list[iscan].j;
			integer k = tck_int_list[iscan].k;

		p.x = 0.5*(xpos[i] + xpos[i + 1]);
		p.y = 0.5*(ypos[j] + ypos[j + 1]);
		p.z = 0.5*(zpos[k] + zpos[k + 1]);
		integer iP = 0;
		switch (iflag) {
		case TEMPERATURE:
			// закоментирован медленный участок кода.
			//inDomain = in_model_temp(p, ib, b, lb);
		    iP = evt[i + j*inx + k*inx*iny];
			if (iP > 0) {
				ib = whot_is_block[iP-1];
				if (b[ib].itype == HOLLOW) {
					inDomain = false;
				}
				else inDomain = true;
			}
			else inDomain = false;			
			break;
		case HYDRODINAMIC: 
			// закоментирован медленный участок кода.
			//inDomain = in_model_flow(p, ib, b, lb);
			iP = evt[i + j*inx + k*inx*iny];
			if (iP > 0) {
				ib = whot_is_block[iP - 1];
				if (b[ib].itype != FLUID) {
					inDomain = false;
				}
				else inDomain = true;
			}
			else inDomain = false;
			break;
		}

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
	SetLength_point(pa, (inx + 1)*(iny + 1)*(inz + 1), maxnode);
	if (bvisit != NULL) {
		delete[] bvisit;
		bvisit = NULL;
	}
} // constr_nodes





// создаёт массив узлов принадлежащих расчётной области
void constr_nodes_flow(TOCHKA* &pa, integer &maxnode, integer* &ent,
	integer **evt_f2, integer iDom,
	integer inx, integer iny, integer inz,
	doublereal *xpos, doublereal *ypos, doublereal *zpos) {

	const integer iSIZE = (inx + 1)*(iny + 1)*(inz + 1);

	// ent - глобальная нумерация узлов.
	// нумерация узлов начинается с единицы, если узел не принадлежит
	// расчётной области то ставится значение 0.
	//ent = NULL;
	ent = new integer[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (ent == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for ent flow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// был ли добавлен узел с номером node
	bool *bvisit = NULL;
	bvisit = new bool[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (bvisit == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit flow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	// Выделение памяти. Сразу много для наискорейшего выполнения:
	//pa = NULL;
	pa = new TOCHKA[iSIZE];
	// Оператор new не требует после себя проверки на null.
	//if (pa == NULL) {
		// недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for pa flow constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}


	// Инициализация:
	integer i=0, j=0, k=0;
	for (i = 0; i<iSIZE; i++) {
		ent[i] = 0;
		bvisit[i] = false;
	}

	bool inDomain = false;
	integer node, ivol;
	TOCHKA pnew;
	maxnode = 0;// maxnode - текущая длина массива pa
	// Цикл по всем контрольным объёмам
	for (i = 0; i<inx; i++) for (j = 0; j<iny; j++) for (k = 0; k<inz; k++) {

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
	//if (bvisit != NULL) {
	// оператор delete может применяться повторно в том числе и к null указателю.
		delete[] bvisit;
		bvisit = NULL;
	//}
} // constr_nodes_flow



// для каждого контрольного объёма принадлежащему
// расчётной области определяет номера его вершин.
void constr_nvtx(integer* evt, integer* ent, integer** &nvtx, integer &maxelm,
	integer inx, integer iny, integer inz, TOCKA_INT* &tck_int_list) {

	// нумерация ent начинается с единицы.

	// проити по всем контрольным объёмам принадлежащим
	// расчётной области
	// maxelm - число контрольных объёмов принадлежащих расчётной области
	//maxelm = 0;
	integer ic, l = 0;
	integer i1, i2, i3, i4, i5, i6, i7, i8;
	// подсчёт количества контрольных объёмов 
	// принадлежащих расчётной области.
	/*
	for (i = 0; i<(inx); i++) for (j = 0; j<(iny); j++) for (k = 0; k<(inz); k++) {
		ic = i + j*inx + k*inx*iny;
		if (evt[ic] > 0) {
			maxelm++;
		}
	}
	*/
	nvtx = NULL;
	nvtx = new integer*[8];
	if (nvtx == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for nvtx constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer i = 0; i<8; i++) nvtx[i] = NULL;
	for (integer i = 0; i<8; i++) {
		nvtx[i] = new integer[maxelm];
		if (nvtx[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for nvtx[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for nvtx[%d] constr struct...\n", i);
#endif
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// Обход по всем контрольным объёмам.
	//for (i = 0; i<(inx); i++) for (j = 0; j<(iny); j++) for (k = 0; k<(inz); k++) {
		//ic = i + j*inx + k*inx*iny;
		//if (evt[ic]>0) {
	for (integer iscan=0; iscan<maxelm; iscan++) {
		{
			integer i = tck_int_list[iscan].i;
			integer j = tck_int_list[iscan].j;
			integer k = tck_int_list[iscan].k;
			ic = i + j * inx + k * inx*iny;

			// контрольный объём принадлежит расчётной области.
			i1 = i + j*(inx + 1) + k*(inx + 1)*(iny + 1);
			i2 = (i + 1) + j*(inx + 1) + k*(inx + 1)*(iny + 1);
			i3 = i + (j + 1)*(inx + 1) + k*(inx + 1)*(iny + 1);
			i4 = (i + 1) + (j + 1)*(inx + 1) + k*(inx + 1)*(iny + 1);
			i5 = i + j*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			i6 = (i + 1) + j*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			i7 = i + (j + 1)*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			i8 = (i + 1) + (j + 1)*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
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

				if ((ent[i1] == 0) || (ent[i2] == 0) || (ent[i3] == 0) || (ent[i4] == 0) || (ent[i5] == 0) || (ent[i6] == 0) || (ent[i7] == 0) || (ent[i8] == 0)) {
#if doubleintprecision == 1
					printf("error : ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#else
					printf("error : ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#endif
					//system("PAUSE");
					system("pause");
				}
			}
			else {
#if doubleintprecision == 1
				printf("error TEMP ent <= 0: ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
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
void walk_in_octree_icolor_different_fluid_domain(octTree* &oc, integer inx, integer iny, integer inz, integer* &icolor_different_fluid_domain, integer** evt_f2) {


	top_ALICE_STACK = 0;
	if (oc->link0 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link0);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link0->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link0->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link0->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link0->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link0->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link0->minz;
		top_ALICE_STACK++;
	}
	if (oc->link1 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link1);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link1->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link1->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link1->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link1->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link1->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link1->minz;
		top_ALICE_STACK++;
	}
	if (oc->link2 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link2);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link2->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link2->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link2->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link2->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link2->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link2->minz;
		top_ALICE_STACK++;
	}
	if (oc->link3 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link3);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link3->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link3->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link3->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link3->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link3->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link3->minz;
		top_ALICE_STACK++;
	}
	if (oc->link4 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link4);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link4->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link4->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link4->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link4->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link4->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link4->minz;
		top_ALICE_STACK++;
	}
	if (oc->link5 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link5);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link5->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link5->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link5->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link5->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link5->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link5->minz;
		top_ALICE_STACK++;
	}
	if (oc->link6 != NULL) {
		my_ALICE_STACK[top_ALICE_STACK].link = (oc->link6);
		my_ALICE_STACK[top_ALICE_STACK].maxx = oc->link6->maxx;
		my_ALICE_STACK[top_ALICE_STACK].minx = oc->link6->minx;
		my_ALICE_STACK[top_ALICE_STACK].maxy = oc->link6->maxy;
		my_ALICE_STACK[top_ALICE_STACK].miny = oc->link6->miny;
		my_ALICE_STACK[top_ALICE_STACK].maxz = oc->link6->maxz;
		my_ALICE_STACK[top_ALICE_STACK].minz = oc->link6->minz;
		top_ALICE_STACK++;
	}
	if (oc->link7 != NULL) {
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
		if (my_ALICE_STACK[top_ALICE_STACK - 1].link != NULL) {
			if (my_ALICE_STACK[top_ALICE_STACK - 1].link->dlist == true) {



				octTree* octree1 = my_ALICE_STACK[top_ALICE_STACK - 1].link;

				// разбиение на 8.
				integer minx = my_ALICE_STACK[top_ALICE_STACK - 1].minx;
				integer maxx = my_ALICE_STACK[top_ALICE_STACK - 1].maxx;
				integer miny = my_ALICE_STACK[top_ALICE_STACK - 1].miny;
				integer maxy = my_ALICE_STACK[top_ALICE_STACK - 1].maxy;
				integer minz = my_ALICE_STACK[top_ALICE_STACK - 1].minz;
				integer maxz = my_ALICE_STACK[top_ALICE_STACK - 1].maxz;

				// Дробление  вызывается.
				my_ALICE_STACK[top_ALICE_STACK - 1].link = NULL;
				top_ALICE_STACK--;

				// Обработка узла.
				if (octree1->inum_FD > 0) {
					integer ielm = octree1->inum_FD - 1; // номер nvtx.
					integer ic = 0;
					ic = minx + miny*inx + minz*inx*iny;
					integer id_found = ielm;

					icolor_different_fluid_domain[id_found] = evt_f2[MASKDOMAINFLUIDCOLOR][ic]; // Присвоение цвета который хранится в третьем поле evt_f2.
				}

			}
			else {
				// продолжаем добираться до листьев.
				STACK_ALICE buf1 = my_ALICE_STACK[top_ALICE_STACK - 1];
				STACK_ALICE* buf = &buf1;
				top_ALICE_STACK--;
				if (buf->link->link0 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link0);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link0->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link0->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link0->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link0->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link0->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link0->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link1 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link1);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link1->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link1->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link1->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link1->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link1->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link1->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link2 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link2);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link2->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link2->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link2->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link2->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link2->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link2->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link3 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link3);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link3->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link3->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link3->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link3->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link3->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link3->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link4 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link4);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link4->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link4->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link4->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link4->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link4->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link4->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link5 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link5);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link5->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link5->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link5->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link5->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link5->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link5->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link6 != NULL) {
					my_ALICE_STACK[top_ALICE_STACK].link = (buf->link->link6);
					my_ALICE_STACK[top_ALICE_STACK].maxx = buf->link->link6->maxx;
					my_ALICE_STACK[top_ALICE_STACK].minx = buf->link->link6->minx;
					my_ALICE_STACK[top_ALICE_STACK].maxy = buf->link->link6->maxy;
					my_ALICE_STACK[top_ALICE_STACK].miny = buf->link->link6->miny;
					my_ALICE_STACK[top_ALICE_STACK].maxz = buf->link->link6->maxz;
					my_ALICE_STACK[top_ALICE_STACK].minz = buf->link->link6->minz;
					top_ALICE_STACK++;
				}
				if (buf->link->link7 != NULL) {
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
// Данная функция использует поле  evt_f2[MASKDOMAINFLUIDCOLOR] и только.
void constr_icolor_different_fluid_domain_alice(integer maxelm_flow, integer* &icolor_different_fluid_domain,
	integer inx, integer iny, integer inz, integer** evt_f2,  doublereal* &x_pos, doublereal* &y_pos, doublereal* &z_pos, TOCHKA* pa,
	integer** nvtx, octTree* &oc)
{
	integer i = 0, j = 0, k = 0; // Счётчики цикла for.

	// Выделение памяти.
	icolor_different_fluid_domain = NULL;
	icolor_different_fluid_domain = new integer[maxelm_flow];
	if (icolor_different_fluid_domain == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for icolor_different_fluid_domain constr_nvtx_flow in constr_struct...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

#if doubleintprecision == 1
	printf("maxelm_flow=%lld inx=%lld, iny=%lld inz=%lld\n", maxelm_flow, inx, iny, inz);
#else
	printf("maxelm_flow=%d inx=%d, iny=%d inz=%d\n", maxelm_flow, inx, iny, inz);
#endif
	
	/*
	// Чрезвычайно медленный кусок кода.
	for (i = 0; i<(inx); i++) for (j = 0; j<(iny); j++) for (k = 0; k<(inz); k++) {
		integer ic = 0;
		ic = i + j*inx + k*inx*iny;
		// Координаты центра первичной ячейки.
		doublereal x_c = 0.5*(x_pos[i]+x_pos[i+1]);
		doublereal y_c = 0.5*(y_pos[j] + y_pos[j + 1]);
		doublereal z_c = 0.5*(z_pos[k] + z_pos[k + 1]);

		
		// Поиск нужной ic64.
		bool bfound = false;
		integer id_found = -1;
		
		for (integer i_1 = 0; i_1 < maxelm_flow; i_1++) {
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

			if ((x_c>pa[nvtx[0][i_1] - 1].x) && (x_c<pa[nvtx[1][i_1] - 1].x) && (y_c>pa[nvtx[0][i_1] - 1].y) && (y_c<pa[nvtx[3][i_1] - 1].y) && (z_c>pa[nvtx[0][i_1] - 1].z) && (z_c<pa[nvtx[4][i_1] - 1].z)) {
				bfound = true;
				id_found = i_1;
				// В случае нахождения досрочный выход из цикла for.
				break;
			}
		}
		if (bfound) {
			icolor_different_fluid_domain[id_found] = evt_f2[MASKDOMAINFLUIDCOLOR][ic]; // Присвоение цвета который хранится в третьем поле evt_f2.
		}
	}
	*/
	// 9.08.2017 он заменён на более быстрый работающий за один проход.
	walk_in_octree_icolor_different_fluid_domain(oc,  inx, iny, inz, icolor_different_fluid_domain,  evt_f2);
} //  constr_icolor_different_fluid_domain_alice





// для каждого контрольного объёма принадлежащему
// расчётной области определяет номера его вершин.
void constr_nvtx_flow(integer** evt_f2, integer* &icolor_different_fluid_domain, integer iDom, integer* ent,
	integer** &nvtx, integer &maxelm,
	integer inx, integer iny, integer inz) {
	// проити по всем контрольным объёмам принадлежащим
	// расчётной области
	// maxelm - число контрольных объёмов принадлежащих расчётной области
	maxelm = 0;
	integer i, j, k;
	integer ic, l = 0;
	integer i1, i2, i3, i4, i5, i6, i7, i8;
	// подсчёт количества контрольных объёмов 
	// принадлежащих расчётной области.
	for (i = 0; i<(inx); i++) for (j = 0; j<(iny); j++) for (k = 0; k<(inz); k++) {
		ic = i + j*inx + k*inx*iny;
		if ((evt_f2[ENUMERATECONTVOL][ic] > 0) && (evt_f2[MASKDOMAINFLUID][ic] == (iDom + 1))) {   
			maxelm++;
		}
	}
	// Выделение памяти.
	icolor_different_fluid_domain = NULL;
	icolor_different_fluid_domain = new integer[maxelm];
	if (icolor_different_fluid_domain == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for icolor_different_fluid_domain constr_nvtx_flow in constr_struct...\n");
		printf("Please any key to exit...\n");
		exit(1);
	}

	integer ic64 = 0;
	for (i = 0; i<(inx); i++) for (j = 0; j<(iny); j++) for (k = 0; k<(inz); k++) {
		ic = i + j*inx + k*inx*iny;
		if ((evt_f2[ENUMERATECONTVOL][ic] > 0) && (evt_f2[MASKDOMAINFLUID][ic] == (iDom + 1))) {

			icolor_different_fluid_domain[ic64] = evt_f2[MASKDOMAINFLUIDCOLOR][ic]; // Присвоение цвета который хранится в третьем поле evt_f2.
			ic64++;
		}
	}
	nvtx = NULL;
	nvtx = new integer*[8];
	if (nvtx == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for nvtx flow constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (i = 0; i<8; i++) nvtx[i] = NULL;
	for (i = 0; i<8; i++) {
		nvtx[i] = new integer[maxelm];
		if (nvtx[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for nvtx[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for nvtx[%d] constr struct...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// Обход по всем контрольным объёмам.
	for (i = 0; i<(inx); i++) for (j = 0; j<(iny); j++) for (k = 0; k<(inz); k++) {
		ic = i + j*inx + k*inx*iny;
		if ((evt_f2[ENUMERATECONTVOL][ic]>0) && (evt_f2[MASKDOMAINFLUID][ic] == (iDom + 1))) {
			// контрольный объём принадлежит расчётной области.
			i1 = i + j*(inx + 1) + k*(inx + 1)*(iny + 1);
			i2 = (i + 1) + j*(inx + 1) + k*(inx + 1)*(iny + 1);
			i3 = i + (j + 1)*(inx + 1) + k*(inx + 1)*(iny + 1);
			i4 = (i + 1) + (j + 1)*(inx + 1) + k*(inx + 1)*(iny + 1);
			i5 = i + j*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			i6 = (i + 1) + j*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			i7 = i + (j + 1)*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);
			i8 = (i + 1) + (j + 1)*(inx + 1) + (k + 1)*(inx + 1)*(iny + 1);

			

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

					if ((ent[i1] == 0) || (ent[i2] == 0) || (ent[i3] == 0) || (ent[i4] == 0) || (ent[i5] == 0) || (ent[i6] == 0) || (ent[i7] == 0) || (ent[i8] == 0)) {
#if doubleintprecision == 1
						printf("error : ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#else
						printf("error : ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#endif
						//system("PAUSE");
						system("pause");
					}
			}
			else {
#if doubleintprecision == 1
				printf("error FLOW ent <= 0: ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld ent[%lld]=%lld \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#else
				printf("error FLOW ent <= 0: ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d ent[%d]=%d \n", i1, ent[i1], i2, ent[i2], i3, ent[i3], i4, ent[i4], i5, ent[i5], i6, ent[i6], i7, ent[i7], i8, ent[i8]);
#endif
				//system("PAUSE");
				system("pause");
			}

		}
	}
} // constr_nvtx_flow

// Заносит свойства материалов в структуру
void constr_prop(integer* evt, integer* &whot_is_block, integer* ent, doublereal** &prop, integer maxelm, integer iflag, BLOCK* b,
	integer lb, integer inx, integer iny, integer inz, doublereal* &Sc, integer* &ipower_time_depend,
				 doublereal *xpos, doublereal *ypos, doublereal *zpos, TPROP* matlist, TOCKA_INT* &tck_int_list) {
	integer i;
	prop=NULL;
	prop=new doublereal*[9];
	if (prop==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for prop constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (i=0; i<9; i++) prop[i]=NULL;
	for (i=0; i<9; i++) {
		prop[i]=new doublereal[maxelm];
		if (prop[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for prop[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for prop[%d] constr struct...\n", i);
#endif
		    
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
	    }
	}
	Sc=NULL;
	Sc=new doublereal[maxelm];
	if (Sc==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for Sc constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	ipower_time_depend = NULL;
	ipower_time_depend = new integer[maxelm];
	if (ipower_time_depend == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for ipower_time_depend constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}

	TOCHKA p;
    bool inDomain=false;
	integer ib,l=0; 
	
	// Проход по всем контрольным объёмам в порядке обхода:
	//for (integer i=0; i<inx; i++) for (integer j=0; j<iny; j++) for (integer k=0; k<inz; k++) {

		//integer ic = i + j*inx + k*inx*iny;
		//if (evt[ic] > 0) {
	for (integer iscan = 0; iscan<maxelm; iscan++) {
		{
			integer i = tck_int_list[iscan].i;
			integer j = tck_int_list[iscan].j;
			integer k = tck_int_list[iscan].k;
			integer ic = i + j * inx + k * inx*iny;

				p.x = 0.5*(xpos[i] + xpos[i + 1]);
				p.y = 0.5*(ypos[j] + ypos[j + 1]);
				p.z = 0.5*(zpos[k] + zpos[k + 1]);
				/*
				switch (iflag) {
				  case TEMPERATURE: inDomain = in_model_temp(p, ib, b, lb); break;
				  case HYDRODINAMIC: inDomain = in_model_flow(p, ib, b, lb); break;
				}
				*/
				integer iP = 0;
				switch (iflag) {
				case TEMPERATURE:
					// закоментирован медленный участок кода.
					//inDomain = in_model_temp(p, ib, b, lb);
					iP = evt[i + j*inx + k*inx*iny];
					if (iP > 0) {
						ib = whot_is_block[iP - 1];
						if (b[ib].itype == HOLLOW) {
							inDomain = false;
						}
						else inDomain = true;
					}
					else inDomain = false;
					break;
				case HYDRODINAMIC:
					// закоментирован медленный участок кода.
					//inDomain = in_model_flow(p, ib, b, lb);
					iP = evt[i + j*inx + k*inx*iny];
					if (iP > 0) {
						ib = whot_is_block[iP - 1];
						if (b[ib].itype != FLUID) {
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
						switch (iflag) {
						case TEMPERATURE:  prop[RHO][l] = matlist[b[ib].imatid].rho;
							//prop[HEAT_CAPACITY][l] = matlist[b[ib].imatid].cp;
							//prop[LAM][l] = matlist[b[ib].imatid].lam;
							prop[HEAT_CAPACITY][l] = get_cp(matlist[b[ib].imatid].n_cp, matlist[b[ib].imatid].temp_cp, matlist[b[ib].imatid].arr_cp, 25.0);
							prop[LAM][l] = get_lam(matlist[b[ib].imatid].n_lam, matlist[b[ib].imatid].temp_lam, matlist[b[ib].imatid].arr_lam, 25.0);

							prop[MULT_LAM_X][l] = matlist[b[ib].imatid].orthotropy_multiplyer_x;
							prop[MULT_LAM_Y][l] = matlist[b[ib].imatid].orthotropy_multiplyer_y;
							prop[MULT_LAM_Z][l] = matlist[b[ib].imatid].orthotropy_multiplyer_z;
							// Коэффициенты Лямэ.
							prop[MU_LAME][l] = matlist[b[ib].imatid].mu_Lame;
							prop[LAMBDA_LAME][l] = matlist[b[ib].imatid].lambda_Lame;
							prop[BETA_T_MECHANICAL][l] = matlist[b[ib].imatid].beta_t_solid;

							//Sc[evt[i + j*inx + k*inx*iny] - 1] = b[ib].Sc;
							Sc[evt[i + j*inx + k*inx*iny] - 1] = get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, 20.0);
							ipower_time_depend[evt[i + j*inx + k*inx*iny] - 1] = b[ib].ipower_time_depend;

							//Sc[l]= get_power(b[ib].n_Sc, b[ib].temp_Sc, b[ib].arr_Sc, 20.0);
							//ipower_time_depend[l] = b[ib].ipower_time_depend;
							break;
						case HYDRODINAMIC: prop[RHO][l] = matlist[b[ib].imatid].rho;
							prop[MU][l] = matlist[b[ib].imatid].mu;
							prop[BETA_T][l] = matlist[b[ib].imatid].beta_t;
							break;
						}

						l++;
					}
					else {
						printf("model is incorrect : prop indx >=maxelm\n");
						printf("see constr struct file.\n");
						system("pause");
						exit(1);
					}
				}
		}
	}
} // constr_prop

// Заносит свойства материалов в структуру
// для жидкой зоны с номером iDom.
void constr_prop_flow(integer* evt, integer* &whot_is_block, 
	integer** evt_f2, integer iDom, doublereal** &prop,
					  integer maxelm, BLOCK* b, integer lb, 
					  integer inx, integer iny, integer inz, 
				      doublereal *xpos, doublereal *ypos, doublereal *zpos, TPROP* matlist) {
	integer i;
	prop=NULL;
	prop=new doublereal*[3];
	if (prop==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for prop flow constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (i=0; i<3; i++) prop[i]=NULL;
	for (i=0; i<3; i++) {
		prop[i]=new doublereal[maxelm];
		if (prop[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for prop[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for prop[%d] constr struct...\n", i);
#endif
		    
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}

	integer j,k,ic;
	//TOCHKA p;
    bool inDomain=false;
	integer ib,l=0; 
	
	// Проход по всем контрольным объёмам в порядке обхода:
	for (i=0; i<inx; i++) for (j=0; j<iny; j++) for (k=0; k<inz; k++) {
		ic=i+j*inx+k*inx*iny; // номер КО
        
		if ((evt_f2[ENUMERATECONTVOL][ic]>0) && (evt_f2[MASKDOMAINFLUID][ic]==(iDom+1))) {

					//p.x = 0.5*(xpos[i] + xpos[i + 1]);
					//p.y = 0.5*(ypos[j] + ypos[j + 1]);
					//p.z = 0.5*(zpos[k] + zpos[k + 1]);

					//inDomain = in_model_flow(p, ib, b, lb); // определяем номер блока жидкой зоны

					integer iP = evt[i + j*inx + k*inx*iny];
					if (iP > 0) {
						ib = whot_is_block[iP - 1];
						if (b[ib].itype != FLUID) {
							inDomain = false;
						}
						else inDomain = true;
					}
					else inDomain = false;


					if (inDomain) {
						if (l < maxelm) {
							prop[RHO][l] = matlist[b[ib].imatid].rho;
							prop[MU][l] = matlist[b[ib].imatid].mu;
							prop[BETA_T][l] = matlist[b[ib].imatid].beta_t;

							l++;
						}
						else {
							printf("model is incorrect: prop 2 indx>=maxelm\n");
							printf("see file constr struct in function constr_prop_flow\n");
							system("pause");
							exit(1);
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
void enumerate_gran_temp(integer** &gran, integer maxelm, integer** nvtx,
						 integer &maxbound, integer **sosed, TOCHKA* pa,
						 SOURCE* s, integer ls) {
	// Выделение оперативной памяти.
	if (gran != NULL) {
		for (integer i = 0; i < 6; i++) {
			if (gran[i] != NULL) {
				delete[] gran[i];
				gran[i] = NULL;
			}
		}
		delete[] gran;
	}
	gran = NULL;

	gran=NULL;
	gran=new integer*[6];
	if (gran==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for gran constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	integer i=0, j=0;
	for (i=0; i<6; i++) gran[i]=NULL;
	for (i=0; i<6; i++) {
		gran[i]=new integer[maxelm];
		if (gran[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for gran[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for gran[%d] constr struct...\n", i);
#endif
		    
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}
	// инициализация:
    for (i=0; i<6; i++) for (j=0; j<maxelm; j++) gran[i][j]=-1;

	bool *bvisit=NULL;
	bvisit=new bool[maxelm];
	// Результат работы оператора new не требует проверки на null.
	//if (bvisit==NULL) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
    for (j=0; j<maxelm; j++) bvisit[j]=false; // признак посещения узла.

	maxbound=0; // общее число граней
	integer G; // текущая грань
	doublereal x_c=0.0, y_c=0.0, z_c=0.0;
	integer iplane;
	bool bfind=false;
	integer GR; // G Reverse
	for (i=0; i<maxelm; i++) {
		// цикл по всем соседям
		for (G=0; G<6; G++) {
			if (sosed[G][i]==0) {
				// сосед не существует
				gran[G][i]=maxbound++;
				// записываем идентификатор грани и
				// увеличиваем счётчик граней.
			}
			else {
				// сосед существует, грань внутренняя
				if (bvisit[sosed[G][i]-1]) {
					// узел уже был посещён
					switch (G) {
						case ESIDE : GR=WSIDE; break;
						case WSIDE: GR=ESIDE; break;
						case NSIDE : GR=SSIDE; break;
						case SSIDE : GR=NSIDE; break;
						case TSIDE : GR=BSIDE; break;
						case BSIDE : GR=TSIDE; break;
					}
					gran[G][i]=gran[GR][sosed[G][i]-1];
				}
				else {
					// узнать координаты центра грани и ориентацию в пространстве
					switch (G) {
						case ESIDE : x_c=pa[nvtx[1][i]-1].x;
							     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
								 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
								 iplane=YZ;
							     break;
						case WSIDE: x_c=pa[nvtx[0][i]-1].x;
							     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
								 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
								 iplane=YZ;
							     break;
						case NSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
							     y_c=pa[nvtx[2][i]-1].y;
								 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
								 iplane=XZ;
							     break;
						case SSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
							     y_c=pa[nvtx[0][i]-1].y;
								 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
								 iplane=XZ;
							     break;
						case TSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                                 y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
								 z_c=pa[nvtx[4][i]-1].z;
                                 iplane=XY;
							     break;
						case BSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                                 y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
								 z_c=pa[nvtx[0][i]-1].z;
                                 iplane=XY;
							     break;
					} // end case
					bfind=false;
					for (j=0; j<ls; j++) {
						if (s[j].iPlane==iplane) {
							switch (iplane) {
							case XY: s[j].g.zE = s[j].g.zS;
									      if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(z_c-s[j].g.zE)<admission)) bfind=true;
									      break;
							case XZ:  s[j].g.yE = s[j].g.yS;
									if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (fabs(y_c-s[j].g.yE)<admission)) bfind=true;
									      break;
							case YZ: s[j].g.xE = s[j].g.xS;
									if ((z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(x_c-s[j].g.xE)<admission)) bfind=true;
									      break;
							}	
						}
						if (bfind) break; // досрочный выход из цикла for.
					}
					if (bfind) {
						// нужно присвоить грани соответствующий номер maxbound.
                        gran[G][i]=maxbound++;
					}
				}
			}
		}
		bvisit[i]=true; // узел посещён
	}

	//if (bvisit != NULL) {
		delete[] bvisit;
		bvisit = NULL;
	//}
	// После этого обхода имеем количество граней maxbound.
	// Для каждой грани знаем её позицию в сетке.

} // enumerate_gran_temp

// Алгоритм уникальной нумерации граней.
// В результате каждая граничная грань,
// принадлежащая расчётной области имеет
// уникальный номер.
void enumerate_gran_flow(integer** &gran, integer maxelm, integer** nvtx,
						 integer &maxbound, integer **sosed, TOCHKA* pa) {
	// Выделение оперативной памяти.
	gran=NULL;
	gran=new integer*[6];
	if (gran==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for gran flow constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	integer i=0, j=0;
	for (i=0; i<6; i++) gran[i]=NULL;
	for (i=0; i<6; i++) {
		gran[i]=new integer[maxelm];
		if (gran[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for gran[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for gran[%d] constr struct...\n", i);
#endif
		   
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
	    }
	}
	// инициализация:
    for (i=0; i<6; i++) for (j=0; j<maxelm; j++) gran[i][j]=-1;

	maxbound=0; // общее число граней
	integer G; // текущая грань

	for (i=0; i<maxelm; i++) {
		// цикл по всем соседям
		for (G=0; G<6; G++) {
			if (sosed[G][i]==0) {
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
void constr_sosedi(ALICE_PARTITION** &sosedi, integer maxelm, integer** gran,
						integer** sosed) {
	//integer i=0;
	// Выделение оперативной памяти.
	sosedi=NULL;
	sosedi = new ALICE_PARTITION*[12];
	if (sosedi==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for sosedi constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer i=0; i<12; i++) sosedi[i]=NULL;
#pragma omp parallel for
	for (integer i=0; i<12; i++) {
		sosedi[i] = new ALICE_PARTITION[maxelm];
		if (sosedi[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for sosedi[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for sosedi[%d] constr struct...\n", i);
#endif
		   
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}

	
	// Обход в цикле по всем внутренним КО принадлежащим расчётной области
#pragma omp parallel for
	for (integer i=0; i<maxelm; i++) {

		integer G, GG; // текущая грань
		// проход по всем граням принадлежащим текущему КО в определённом порядке.
		for (G=0; G<6; G++) {
			if (gran[G][i]>-1) {
				// граничная грань:
				sosedi[G][i].iNODE1=maxelm+gran[G][i]; // граничные КО нумеруются в последнюю очередь,
				sosedi[G][i].bdroblenie4 = false;
				sosedi[G][i].iNODE2 = -1;
				sosedi[G][i].iNODE3 = -1;
				sosedi[G][i].iNODE4 = -1;
				// после maxelm внутренних КО.

				// Вычисление дальнего соседа
				switch (G) {
					case ESIDE : GG=EE; break; // E
					case NSIDE : GG=NN; break; // N
					case TSIDE : GG=TTSIDE; break; // T
					case WSIDE: GG=WW; break; // W
					case SSIDE : GG=SS; break; // S
					case BSIDE : GG=BB; break; // B
				}
				sosedi[GG][i].iNODE1=-1; // соседа нет.
				sosedi[GG][i].bdroblenie4 = false;
				sosedi[GG][i].iNODE2 = -1;
				sosedi[GG][i].iNODE3 = -1;
				sosedi[GG][i].iNODE4 = -1;

			}
			else {
				// соседом является внутренний КО.
				sosedi[G][i].iNODE1=sosed[G][i]-1;
				sosedi[G][i].bdroblenie4 = false;
				sosedi[G][i].iNODE2 = -1;
				sosedi[G][i].iNODE3 = -1;
				sosedi[G][i].iNODE4 = -1;

				// Вычисление дальнего соседа
				switch (G) {
					case ESIDE : GG=EE; break; // E
					case NSIDE : GG=NN; break; // N
					case TSIDE : GG=TTSIDE; break; // T
					case WSIDE: GG=WW; break; // W
					case SSIDE : GG=SS; break; // S
					case BSIDE : GG=BB; break; // B
				}

				if (gran[G][sosed[G][i]-1]>-1) {
					// граничный КО
					sosedi[GG][i].iNODE1=maxelm+gran[G][sosed[G][i]-1];
					sosedi[GG][i].bdroblenie4 = false;
					sosedi[GG][i].iNODE2 = -1;
					sosedi[GG][i].iNODE3 = -1;
					sosedi[GG][i].iNODE4 = -1;
				}
				else {
					sosedi[GG][i].iNODE1 = sosed[GG][i] - 1; // внутренний КО.
					sosedi[GG][i].bdroblenie4 = false;
					sosedi[GG][i].iNODE2 = -1;
					sosedi[GG][i].iNODE3 = -1;
					sosedi[GG][i].iNODE4 = -1;
				}

			}
		}
	}

	// sosedi: maxelm -> maxp == maxelm + maxbound;

	// Признак того что узел граничный: номер КО >= maxelm.
	// Массив boundary (определяющий граничные узлы) становится ненужен.

} // constr_sosedi

// соседи в плоскости границы для граничных
// узлов. Это требуется для симметризации 
// портрета матрицы СЛАУ.
void constr_boundary_sosed(integer G, integer i, BOUND* &sosedb, integer** gran, integer** sosed, ALICE_PARTITION **sosedi, integer maxelm) {
	integer j=0; // счётчик цикла for

    for (j=0; j<6; j++) sosedb[gran[G][i]].iW[j]=-1; // инициализация 

	// 13 10 2016 Выключил данный код т.к. он считается ненужным.
	/*
	// Эти значения требуются для симметризации портрета матрицы СЛАУ:
	switch (G) {
	case ESIDE: if ((sosed[NSIDE][i]>0) && (gran[G][sosedi[NSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[NSIDE] = maxelm + gran[G][sosedi[NSIDE][i].iNODE1];
		if ((sosed[SSIDE][i]>0) && (gran[G][sosedi[SSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[SSIDE] = maxelm + gran[G][sosedi[SSIDE][i].iNODE1];
		if ((sosed[TSIDE][i]>0) && (gran[G][sosedi[TSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[TSIDE] = maxelm + gran[G][sosedi[TSIDE][i].iNODE1];
		if ((sosed[BSIDE][i]>0) && (gran[G][sosedi[BSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[BSIDE] = maxelm + gran[G][sosedi[BSIDE][i].iNODE1];
			  break;
	case WSIDE: if ((sosed[NSIDE][i]>0) && (gran[G][sosedi[NSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[NSIDE] = maxelm + gran[G][sosedi[NSIDE][i].iNODE1];
		if ((sosed[SSIDE][i]>0) && (gran[G][sosedi[SSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[SSIDE] = maxelm + gran[G][sosedi[SSIDE][i].iNODE1];
		if ((sosed[TSIDE][i]>0) && (gran[G][sosedi[TSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[TSIDE] = maxelm + gran[G][sosedi[TSIDE][i].iNODE1];
		if ((sosed[BSIDE][i]>0) && (gran[G][sosedi[BSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[BSIDE] = maxelm + gran[G][sosedi[BSIDE][i].iNODE1];
			  break;
	case NSIDE: if ((sosed[ESIDE][i]>0) && (gran[G][sosedi[ESIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[ESIDE] = maxelm + gran[G][sosedi[ESIDE][i].iNODE1];
		if ((sosed[WSIDE][i]>0) && (gran[G][sosedi[WSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[WSIDE] = maxelm + gran[G][sosedi[WSIDE][i].iNODE1];
		if ((sosed[TSIDE][i]>0) && (gran[G][sosedi[TSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[TSIDE] = maxelm + gran[G][sosedi[TSIDE][i].iNODE1];
		if ((sosed[BSIDE][i]>0) && (gran[G][sosedi[BSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[BSIDE] = maxelm + gran[G][sosedi[BSIDE][i].iNODE1];
			  break;
	case SSIDE: if ((sosed[ESIDE][i]>0) && (gran[G][sosedi[ESIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[ESIDE] = maxelm + gran[G][sosedi[ESIDE][i].iNODE1];
		if ((sosed[WSIDE][i]>0) && (gran[G][sosedi[WSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[WSIDE] = maxelm + gran[G][sosedi[WSIDE][i].iNODE1];
		if ((sosed[TSIDE][i]>0) && (gran[G][sosedi[TSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[TSIDE] = maxelm + gran[G][sosedi[TSIDE][i].iNODE1];
		if ((sosed[BSIDE][i]>0) && (gran[G][sosedi[BSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[BSIDE] = maxelm + gran[G][sosedi[BSIDE][i].iNODE1];
			  break;
	case TSIDE: if ((sosed[ESIDE][i]>0) && (gran[G][sosedi[ESIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[ESIDE] = maxelm + gran[G][sosedi[ESIDE][i].iNODE1];
		if ((sosed[WSIDE][i]>0) && (gran[G][sosedi[WSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[WSIDE] = maxelm + gran[G][sosedi[WSIDE][i].iNODE1];
		if ((sosed[NSIDE][i]>0) && (gran[G][sosedi[NSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[NSIDE] = maxelm + gran[G][sosedi[NSIDE][i].iNODE1];
		if ((sosed[SSIDE][i]>0) && (gran[G][sosedi[SSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[SSIDE] = maxelm + gran[G][sosedi[SSIDE][i].iNODE1];
		      break;
	case BSIDE : if ((sosed[ESIDE][i]>0) && (gran[G][sosedi[ESIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[ESIDE] = maxelm + gran[G][sosedi[ESIDE][i].iNODE1];
		if ((sosed[WSIDE][i]>0) && (gran[G][sosedi[WSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[WSIDE] = maxelm + gran[G][sosedi[WSIDE][i].iNODE1];
		if ((sosed[NSIDE][i]>0) && (gran[G][sosedi[NSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[NSIDE] = maxelm + gran[G][sosedi[NSIDE][i].iNODE1];
		if ((sosed[SSIDE][i]>0) && (gran[G][sosedi[SSIDE][i].iNODE1]>-1)) sosedb[gran[G][i]].iW[SSIDE] = maxelm + gran[G][sosedi[SSIDE][i].iNODE1];
		      break;
	} // end сбор информации для симметризации СЛАУ.

	*/

} // constr_boundary_sosed

// Заполнение информации о граничных узлах:
// Заполняются sosedb и binternalsource.
void constr_sosedb_temp(BOUND* &sosedb, integer* &whot_is_block, bool* &binternalsource, integer maxelm, integer maxbound, 
	integer** gran, integer** sosed, ALICE_PARTITION **sosedi, integer **nvtx, TOCHKA* pa,
						BLOCK* b, integer lb, integer lw, WALL* w, SOURCE* s, integer ls) {

	// Алгоритм:
    // Выделение оперативной памяти:
    // gran[0..5][0..maxbound-1]
	sosedb=NULL;
	sosedb= new BOUND[maxbound];
	// оператор new не требует после себя проверки на null
	//if (sosedb==NULL) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for sosedb constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	binternalsource=NULL;
    binternalsource = new bool[maxbound];
	// оператор new не требует после себя проверки на null
	//if (binternalsource==NULL) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for binternalsource constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	//bool bfluidsolid=false;

    // Проход по всем внутренним КО принадлежащим расчётной области.
	integer i=0,j,jpos; // счётчик цикла.
	for (i=0; i<maxbound; i++) binternalsource[i]=false; // инициализация.
	bool bfind; 
	integer G=-1; // текущая грань.
    doublereal x_c=0.0, y_c=0.0, z_c=0.0; // координаты центра грани
	integer iplane; // плоскость в которой лежит грань.
	TOCHKA p_c;
	integer ib; // номер блока
	bool bi_fluid; // КО i принадлежит жидкой зоне.

	bool *bvisit=NULL;
	bvisit=new bool[maxbound];
	// оператор new не требует после себя проверки на null
	//if (bvisit==NULL) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	for (i=0; i<maxbound; i++) bvisit[i]=false;

	//integer il=0; // debug
    //for (i=0; i<maxelm; i++) for (G=0; G<6; G++) if (gran[G][i]>-1) il++;
    //for (i=0; i<maxelm; i++) for (G=0; G<6; G++) if (sosed[G][i]==0) il++;
#if doubleintprecision == 1
	//printf("il=%lld\n",il); system("PAUSE");
#else
	//printf("il=%d\n",il); system("PAUSE");
#endif
   

	for (i=0; i<maxelm; i++) {
		// проход по всем граням внутреннего КО в особом порядке.
		for (G=0; G<6; G++) {
			if (gran[G][i]>-1) {
				// граничная грань
				doublereal dS = 0.0;
                // узнать координаты центра грани и ориентацию в пространстве
				switch (G) {
					case ESIDE : x_c=pa[nvtx[1][i]-1].x;
						     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=YZ;
						     break;
					case WSIDE: x_c=pa[nvtx[0][i]-1].x;
						     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=YZ;
						     break;
					case NSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						     y_c=pa[nvtx[2][i]-1].y;
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=XZ;
						     break;
					case SSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						     y_c=pa[nvtx[0][i]-1].y;
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=XZ;
						     break;
					case TSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                             y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=pa[nvtx[4][i]-1].z;
                             iplane=XY;
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y);
						     break;
				    case BSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                             y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
					    	 z_c=pa[nvtx[0][i]-1].z;
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y);
                             iplane=XY;
						     break;
				} // end case

				if (sosed[G][i]==0) { // грань лежит на границе расчётной области.

					// соседи в плоскости границы для граничных
                    // узлов. Это требуется для симметризации 
                    // портрета матрицы СЛАУ.
                    constr_boundary_sosed(G,i,sosedb, gran, sosed, sosedi, maxelm);

					

					// грань лежит на границе расчётной области.
					switch (G) {
						case ESIDE : sosedb[gran[G][i]].Norm=WSIDE;
							     sosedb[gran[G][i]].iII=sosedi[WSIDE][i].iNODE1;
							     break;
						case WSIDE: sosedb[gran[G][i]].Norm=ESIDE;
							sosedb[gran[G][i]].iII = sosedi[ESIDE][i].iNODE1;
							     break;
						case NSIDE : sosedb[gran[G][i]].Norm=SSIDE;
							sosedb[gran[G][i]].iII = sosedi[SSIDE][i].iNODE1;
							     break;
						case SSIDE: sosedb[gran[G][i]].Norm = NSIDE;
							sosedb[gran[G][i]].iII = sosedi[NSIDE][i].iNODE1;
							     break;
						case TSIDE : sosedb[gran[G][i]].Norm=BSIDE; 
							sosedb[gran[G][i]].iII = sosedi[BSIDE][i].iNODE1;
							     break;
						case BSIDE : sosedb[gran[G][i]].Norm=TSIDE; 
							sosedb[gran[G][i]].iII = sosedi[TSIDE][i].iNODE1;
							     break;
					} // end определеиие внутренней нормали
                    sosedb[gran[G][i]].iB=maxelm+gran[G][i];
					sosedb[gran[G][i]].iI=i;
					sosedb[gran[G][i]].iI1 = -1;
					sosedb[gran[G][i]].iI2 = -1;
					bvisit[gran[G][i]]=true; // грань была посещена
					sosedb[gran[G][i]].dS = dS;
					// координаты центра грани.
					sosedb[gran[G][i]].p_c.x = x_c;
					sosedb[gran[G][i]].p_c.y = y_c;
					sosedb[gran[G][i]].p_c.z = z_c;
					// Вычисляем emissivity:
					integer ibfound = -1;
					ibfound = whot_is_block[i];
					if (ibfound >= 0) {
						// блок найден.
						// Определяем внутреннюю нормаль:
						switch (sosedb[gran[G][i]].Norm) {
						case WSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissE;
							break;
						case ESIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissW;
							break;
						case SSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissN;
							break;
						case NSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissS;
							break;
						case BSIDE : sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissT;
							break;
						case TSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissB;
							break;
						}
					}
					else {
						printf("error: emissivity calculation block unfound. in constr_sosedb_temp.\n");
						system("PAUSE");
						exit(1);
					}

					bfind=false;
					for (j=0; j<ls; j++) {
						if (s[j].iPlane==iplane) {
							// Важно не только попадание в фокус объекта но и нахождение с объектом на одном уровне:
							switch (iplane) {
								case XY : if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(z_c-s[j].g.zE)<admission)) { bfind=true; jpos=j; } break;
								case YZ : if ((z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(x_c-s[j].g.xE)<admission)) { bfind=true; jpos=j; } break;
								case XZ : if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (fabs(y_c-s[j].g.yE)<admission)) { bfind=true; jpos=j; } break;
							}
						}
					}
					if (bfind) {
						// внешний источник тепла
						//printf("source out found...\n"); // debug
						//system("PAUSE");
						sosedb[gran[G][i]].MCB=jpos;
					}
					else {
                       for (j=0; j<lw; j++) {
						  if (w[j].iPlane==iplane) {
							  // Важно не только попадание в фокус объекта но и нахождение с объектом на одном уровне:
						      switch (iplane) {
								case XY : if ((x_c>w[j].g.xS) && (x_c<w[j].g.xE) && (y_c>w[j].g.yS) && (y_c<w[j].g.yE) && (fabs(z_c-w[j].g.zE)<admission)) { bfind=true; jpos=j; } break;
								case YZ : if ((z_c>w[j].g.zS) && (z_c<w[j].g.zE) && (y_c>w[j].g.yS) && (y_c<w[j].g.yE) && (fabs(x_c-w[j].g.xE)<admission)) { bfind=true; jpos=j; } break;
								case XZ : if ((x_c>w[j].g.xS) && (x_c<w[j].g.xE) && (z_c>w[j].g.zS) && (z_c<w[j].g.zE) && (fabs(y_c-w[j].g.yE)<admission)) { bfind=true; jpos=j; } break;
							  }
						  }
					   }
					   if (bfind) {
						   // найдена стека
						   //printf("ambient wall found...\n"); // debug
						   //system("PAUSE");
                           sosedb[gran[G][i]].MCB=ls+jpos;
					   }
					   else {
						   // граничное условие по умолчанию
						   sosedb[gran[G][i]].MCB=ls+lw;
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
						
						p_c.x=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						p_c.y=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
						p_c.z=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);

						bi_fluid=in_model_flow(p_c,ib,b,lb);

						if (bi_fluid) {
							// КО i принадлежит жидкой зоне

							// Записываем в соседи тот КО который принадлежит твёрдому телу
							// с противоположной от жидкого КО стороны рассматриваемой грани.
                            sosedb[gran[G][i]].Norm=G; // внутренняя нормаль
							sosedb[gran[G][i]].iB=maxelm+gran[G][i];
							sosedb[gran[G][i]].iI=sosed[G][i]-1;
							sosedb[gran[G][i]].iI1 = sosed[G][i] - 1;
							sosedb[gran[G][i]].iII = sosedi[G][sosed[G][i] - 1].iNODE1;
							sosedb[gran[G][i]].dS = dS;
							// координаты центра грани.
							sosedb[gran[G][i]].p_c.x = x_c;
							sosedb[gran[G][i]].p_c.y = y_c;
							sosedb[gran[G][i]].p_c.z = z_c;
							//bfluidsolid=true;							

							// Вычисляем emissivity:
							integer ibfound = -1;
							ibfound = whot_is_block[sosed[G][i] - 1];
							if (ibfound >= 0) {
								// блок найден.
								// Определяем внутреннюю нормаль:
								switch (sosedb[gran[G][i]].Norm) {
								case WSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissE;
									break;
								case ESIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissW;
									break;
								case SSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissN;
									break;
								case NSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissS;
									break;
								case BSIDE : sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissT;
									break;
								case TSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissB;
									break;
								}
							}
							else {
								printf("error: emissivity calculation block unfound. in constr_sosedb_temp.\n");
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
						       case ESIDE : sosedb[gran[G][i]].Norm=WSIDE;
								   sosedb[gran[G][i]].iII = sosedi[WSIDE][i].iNODE1;
							            break;
						       case WSIDE: sosedb[gran[G][i]].Norm=ESIDE;
								   sosedb[gran[G][i]].iII = sosedi[ESIDE][i].iNODE1;
							            break;
						       case NSIDE : sosedb[gran[G][i]].Norm=SSIDE;
								   sosedb[gran[G][i]].iII = sosedi[SSIDE][i].iNODE1;
							            break;
						       case SSIDE : sosedb[gran[G][i]].Norm=NSIDE;
								   sosedb[gran[G][i]].iII = sosedi[NSIDE][i].iNODE1;
							            break;
						       case TSIDE : sosedb[gran[G][i]].Norm=BSIDE; 
								   sosedb[gran[G][i]].iII = sosedi[BSIDE][i].iNODE1;
							            break;
						       case BSIDE : sosedb[gran[G][i]].Norm=TSIDE; 
								   sosedb[gran[G][i]].iII = sosedi[TSIDE][i].iNODE1;
							            break;
					         } // end определеиие внутренней нормали
                            sosedb[gran[G][i]].iB=maxelm+gran[G][i];
					        sosedb[gran[G][i]].iI=i;
							sosedb[gran[G][i]].iI1 = i;
							sosedb[gran[G][i]].dS = dS;
							// координаты центра грани.
							sosedb[gran[G][i]].p_c.x = x_c;
							sosedb[gran[G][i]].p_c.y = y_c;
							sosedb[gran[G][i]].p_c.z = z_c;
							//bfluidsolid=true;
							// здесь подразумевается что для всех внутренних граней
							// имеем ситуацию когда источник посередине между SOLID и FLUID.

							// Вычисляем emissivity:
							integer ibfound = -1;
							ibfound = whot_is_block[i];
							if (ibfound >= 0) {
								// блок найден.
								// Определяем внутреннюю нормаль:
								switch (sosedb[gran[G][i]].Norm) {
								case WSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissE;
									break;
								case ESIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissW;
									break;
								case SSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissN;
									break;
								case NSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissS;
									break;
								case BSIDE : sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissT;
									break;
								case TSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissB;
									break;
								}
							}
							else {
								printf("error: emissivity calculation block unfound. in constr_sosedb_temp.\n");
								system("PAUSE");
								exit(1);
							}

						}

						// соседи в плоскости границы для граничных
                        // узлов. Это требуется для симметризации 
                        // портрета матрицы СЛАУ.
                        constr_boundary_sosed(G,i,sosedb, gran, sosed, sosedi, maxelm);

						// Нужно определить MCB:
                        bfind=false;
					    for (j=0; j<ls; j++) {
						   if (s[j].iPlane==iplane) {
							  switch (iplane) {
									case XY : if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(z_c-s[j].g.zE)<admission)) { bfind=true; jpos=j; } break;
								    case YZ : if ((z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (y_c>s[j].g.yS) && (y_c<s[j].g.yE) && (fabs(x_c-s[j].g.xE)<admission)) { bfind=true; jpos=j; } break;
								    case XZ : if ((x_c>s[j].g.xS) && (x_c<s[j].g.xE) && (z_c>s[j].g.zS) && (z_c<s[j].g.zE) && (fabs(y_c-s[j].g.yE)<admission)) { bfind=true; jpos=j; } break;
							  }
						   }
					    }
						
					    if (bfind) {
						    // внутрений источник тепла
							//printf("internal source found...\n"); // debug
							//system("PAUSE");
						    sosedb[gran[G][i]].MCB=jpos;
							binternalsource[gran[G][i]]=true; // внутренний источник тепла на границе жидкости и твёрдого тела
						}
						else {
							// Ошибка внутренний источник тепла не обнаружен!
							printf("Error !!! internal source not found...\n");
							//system("PAUSE");
							system("pause");
							exit(0);
							sosedb[gran[G][i]].MCB=ls+lw;
						}

					}
					else
					{
						p_c.x = 0.5*(pa[nvtx[0][i] - 1].x + pa[nvtx[1][i] - 1].x);
						p_c.y = 0.5*(pa[nvtx[1][i] - 1].y + pa[nvtx[3][i] - 1].y);
						p_c.z = 0.5*(pa[nvtx[0][i] - 1].z + pa[nvtx[4][i] - 1].z);

						bi_fluid = in_model_flow(p_c, ib, b, lb);

						if (bi_fluid) {
							// КО i принадлежит жидкой зоне

							// Записываем в соседи тот КО который принадлежит твёрдому телу
							// с противоположной от жидкого КО стороны рассматриваемой грани.
							
							sosedb[gran[G][i]].iI2 = sosed[G][i] - 1;
							

						}
						else {
							sosedb[gran[G][i]].iI2 = i;
						}
					}

					bvisit[gran[G][i]]=true; // грань была посещена
				}
			}
		}
	}

	// Освобождение оперативной памяти
	delete[] bvisit; 
	bvisit=NULL;

	// sosedb : maxbound -> maxp == maxelm + maxbound;

} // constr_sosedb_temp 

// Заполнение информации о граничных узлах:
void constr_sosedb_flow(BOUND* &sosedb,  integer* &whot_is_block, integer* ptr, integer maxelm, integer maxbound,
	                     integer** gran, integer** sosed, ALICE_PARTITION **sosedi, integer **nvtx, TOCHKA* pa,
						integer lw, WALL* w, integer ls, BLOCK* b) {

	// Алгоритм:
    // Выделение оперативной памяти:
    // gran[0..5][0..maxbound-1]
	sosedb=NULL;
	sosedb= new BOUND[maxbound];
	if (sosedb==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for sosedb constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}

    // Проход по всем внутренним КО принадлежащим расчётной области.
	integer i=0,j=0,jpos=0; // счётчик цикла.
	bool bfind; 
	integer G=-1; // текущая грань.
    doublereal x_c=0.0, y_c=0.0, z_c=0.0; // координаты центра грани
	integer iplane; // плоскость в которой лежит грань.


	for (i=0; i<maxelm; i++) {
		// проход по всем граням внутреннего КО в особом порядке.
		for (G=0; G<6; G++) {
			if (gran[G][i]>-1) {
				// граничная грань
				
				doublereal dS = 0.0;
                // узнать координаты центра грани и ориентацию в пространстве
				switch (G) {
					case ESIDE : x_c=pa[nvtx[1][i]-1].x;
						     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=YZ;
						     break;
					case WSIDE: x_c=pa[nvtx[0][i]-1].x;
						     y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=YZ;
						     break;
					case NSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						     y_c=pa[nvtx[2][i]-1].y;
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=XZ;
						     break;
					case SSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						     y_c=pa[nvtx[0][i]-1].y;
							 z_c=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[4][i] - 1].z - pa[nvtx[0][i] - 1].z);
							 iplane=XZ;
						     break;
					case TSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                             y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
							 z_c=pa[nvtx[4][i]-1].z;
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y);
                             iplane=XY;
						     break;
				    case BSIDE : x_c=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
                             y_c=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
					    	 z_c=pa[nvtx[0][i]-1].z;
							 dS = fabs(pa[nvtx[1][i] - 1].x - pa[nvtx[0][i] - 1].x)*fabs(pa[nvtx[3][i] - 1].y - pa[nvtx[0][i] - 1].y);
                             iplane=XY;
						     break;
				} // end case

				if (sosed[G][i]==0) {

                    // соседи в плоскости границы для граничных
                    // узлов. Это требуется для симметризации 
                    // портрета матрицы СЛАУ.
                    constr_boundary_sosed(G,i,sosedb, gran, sosed, sosedi, maxelm);


					// грань лежит на границе расчётной области.
					switch (G) {
						case ESIDE : sosedb[gran[G][i]].Norm=WSIDE;
							     sosedb[gran[G][i]].iII=sosedi[WSIDE][i].iNODE1;
							     break;
						case WSIDE: sosedb[gran[G][i]].Norm=ESIDE;
							sosedb[gran[G][i]].iII = sosedi[ESIDE][i].iNODE1;
							     break;
						case NSIDE : sosedb[gran[G][i]].Norm=SSIDE;
							sosedb[gran[G][i]].iII = sosedi[SSIDE][i].iNODE1;
							     break;
						case SSIDE : sosedb[gran[G][i]].Norm=NSIDE;
							sosedb[gran[G][i]].iII = sosedi[NSIDE][i].iNODE1;
							     break;
						case TSIDE : sosedb[gran[G][i]].Norm=BSIDE; 
							sosedb[gran[G][i]].iII = sosedi[BSIDE][i].iNODE1;
							     break;
						case BSIDE : sosedb[gran[G][i]].Norm=TSIDE; 
							sosedb[gran[G][i]].iII = sosedi[TSIDE][i].iNODE1;
							     break;
					} // end определеиие внутренней нормали
                    sosedb[gran[G][i]].iB=maxelm+gran[G][i];
					sosedb[gran[G][i]].iI=i;
					sosedb[gran[G][i]].dS = dS; // площадь грани.
					// координаты центра грани.
					sosedb[gran[G][i]].p_c.x = x_c;
					sosedb[gran[G][i]].p_c.y = y_c;
					sosedb[gran[G][i]].p_c.z = z_c;
					// Вычисляем emissivity:
					integer ibfound = -1;
					ibfound = whot_is_block[ptr[i]];
					if (ibfound >= 0) {
						// блок найден.
						// Определяем внутреннюю нормаль:
						switch (sosedb[gran[G][i]].Norm) {
						case WSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissE;
							break;
						case ESIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissW;
							break;
						case SSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissN;
							break;
						case NSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissS;
							break;
						case BSIDE : sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissT;
							break;
						case TSIDE: sosedb[gran[G][i]].emissivity = b[ibfound].radiation.emissB;
							break;
						}
					}
					else {
						printf("error: emissivity calculation block unfound. in constr_sosedb_temp.\n");
						system("PAUSE");
						exit(1);
					}
					

					bfind=false;
					for (j=0; j<lw; j++) {
						if (w[j].iPlane==iplane) {
							switch (iplane) {
								case XY : if ((x_c>w[j].g.xS) && (x_c<w[j].g.xE) && (y_c>w[j].g.yS) && (y_c<w[j].g.yE) && (fabs(z_c-w[j].g.zE)<admission)) 
										  {
											  //printf("z_c=%f, zE=%f\n",z_c, w[j].g.zE);
											  //system("PAUSE"); // debug
											  bfind=true; jpos=j;
										  }
										  break;
								case YZ : if ((z_c>w[j].g.zS) && (z_c<w[j].g.zE) && (y_c>w[j].g.yS) && (y_c<w[j].g.yE) && (fabs(x_c-w[j].g.xE)<admission)) { bfind=true; jpos=j; } break;
								case XZ : if ((x_c>w[j].g.xS) && (x_c<w[j].g.xE) && (z_c>w[j].g.zS) && (z_c<w[j].g.zE) && (fabs(y_c-w[j].g.yE)<admission)) { bfind=true; jpos=j; } break;
							}
						}
					}
					if (bfind) {
						// внешний источник тепла
						sosedb[gran[G][i]].MCB=ls+jpos;
#if doubleintprecision == 1
						//printf("MCB=%lld\n",sosedb[gran[G][i]].MCB);
#else
						//printf("MCB=%d\n",sosedb[gran[G][i]].MCB);
#endif
						
						// system("PAUSE"); // debug
					}
					else {
                         // граничное условие по умолчанию
					     sosedb[gran[G][i]].MCB=ls+lw;
						 //printf("default\n");
						 //system("PAUSE"); // debug
					}
				} 
			}
		}
	}

	
	// sosedb : maxbound -> maxp == maxelm + maxbound;

} // constr_sosedb_flow

// Свойства материала на границе твердотельной области.
void constr_prop_bound(doublereal **prop, doublereal** &prop_b, integer maxelm, integer maxbound,
					   integer **gran, integer **sosed, integer **nvtx, TOCHKA* pa,
					   BLOCK* b, integer lb) {
	// Выделение оперативной памяти:
	prop_b=NULL;
	prop_b=new doublereal*[6];
	if (prop_b==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for prop_b constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	integer i=0; // счётчик
	for (i=0; i<6; i++) prop_b[i]=NULL;
	for (i=0; i<6; i++) {
		prop_b[i]=new doublereal[maxbound];
		if (prop_b[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for prop_b[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for prop_b[%d] constr struct...\n", i);
#endif
		   
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}
	integer G; // текущая грань
	bool *bvisit=NULL;
	bvisit=new bool[maxbound];
	// Результат работы оператора new не требует проверки на null.
	//if (bvisit==NULL) {
	    // недостаточно памяти на данном оборудовании.
		//printf("Problem : not enough memory on your equipment for bvisit constr struct...\n");
		//printf("Please any key to exit...\n");
		//system("PAUSE");
		//system("pause");
		//exit(1);
	//}
	for (i=0; i<maxbound; i++) bvisit[i]=false; // признак посещаемости узлов
	TOCHKA p_c;
	integer ib;
	bool bi_fluid;


    // в цикле по всем внутренним КО принадлежащим расчётной области:
	for (i=0; i<maxelm; i++) {
		// проход по всем граням внутреннего КО в особом порядке.
		for (G=0; G<6; G++) {
			if (gran[G][i]>-1) {
				// граничная грань
				if (sosed[G][i]==0) {
					// грань лежит на границе расчётной области

					// В граничный узел сносятся свойства прилегающего КО:
					prop_b[RHO][gran[G][i]]=prop[RHO][i];
					prop_b[HEAT_CAPACITY][gran[G][i]]=prop[HEAT_CAPACITY][i];
					prop_b[LAM][gran[G][i]]=prop[LAM][i];
					prop_b[MULT_LAM_X][gran[G][i]] = prop[MULT_LAM_X][i];
					prop_b[MULT_LAM_Y][gran[G][i]] = prop[MULT_LAM_Y][i];
					prop_b[MULT_LAM_Z][gran[G][i]] = prop[MULT_LAM_Z][i];
					bvisit[gran[G][i]]=true;
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
						p_c.x=0.5*(pa[nvtx[0][i]-1].x+pa[nvtx[1][i]-1].x);
						p_c.y=0.5*(pa[nvtx[1][i]-1].y+pa[nvtx[3][i]-1].y);
						p_c.z=0.5*(pa[nvtx[0][i]-1].z+pa[nvtx[4][i]-1].z);

						bi_fluid=in_model_flow(p_c,ib,b,lb);

						if (bi_fluid) {
                            prop_b[RHO][gran[G][i]]=prop[RHO][sosed[G][i]-1];
					        prop_b[HEAT_CAPACITY][gran[G][i]]=prop[HEAT_CAPACITY][sosed[G][i]-1];
					        prop_b[LAM][gran[G][i]]=prop[LAM][sosed[G][i]-1];
							prop_b[MULT_LAM_X][gran[G][i]] = prop[MULT_LAM_X][sosed[G][i] - 1];
							prop_b[MULT_LAM_Y][gran[G][i]] = prop[MULT_LAM_Y][sosed[G][i] - 1];
							prop_b[MULT_LAM_Z][gran[G][i]] = prop[MULT_LAM_Z][sosed[G][i] - 1];
						}
						else {
                            prop_b[RHO][gran[G][i]]=prop[RHO][i];
					        prop_b[HEAT_CAPACITY][gran[G][i]]=prop[HEAT_CAPACITY][i];
					        prop_b[LAM][gran[G][i]]=prop[LAM][i];
							prop_b[MULT_LAM_X][gran[G][i]] = prop[MULT_LAM_X][i];
							prop_b[MULT_LAM_Y][gran[G][i]] = prop[MULT_LAM_Y][i];
							prop_b[MULT_LAM_Z][gran[G][i]] = prop[MULT_LAM_Z][i];
						}

					}

                    bvisit[gran[G][i]]=true;
				}
			}
		}
	}

	delete[] bvisit;
	bvisit=NULL;

} // constr_prop_bound 

// Свойства материала на границе жидкой области.
void constr_prop_bound_flow(doublereal **prop, doublereal** &prop_b, integer maxelm, integer maxbound,
					   integer **gran, integer **sosed) {
	// Выделение оперативной памяти:
	prop_b=NULL;
	prop_b=new doublereal*[3];
	if (prop_b==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for prop_b constr struct...\n");
		printf("Please any key to exit...\n");
	//	system("PAUSE");
		system("pause");
		exit(1);
	}
	integer i=0; // счётчик
	for (i=0; i<3; i++) prop_b[i]=NULL;
	for (i=0; i<3; i++) {
		prop_b[i]=new doublereal[maxbound];
		if (prop_b[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for prop_b[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for prop_b[%d] constr struct...\n", i);
#endif
		   
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}
	integer G; // текущая грань
	


    // в цикле по всем внутренним КО принадлежащим расчётной области:
	for (i=0; i<maxelm; i++) {
		// проход по всем граням внутреннего КО в особом порядке.
		for (G=0; G<6; G++) {
			if (gran[G][i]>-1) {
				// граничная грань
				if (sosed[G][i]==0) {
					// грань лежит на границе расчётной области

					// В граничный узел сносятся свойства прилегающего КО:
					prop_b[RHO][gran[G][i]]=prop[RHO][i];
					prop_b[MU][gran[G][i]]=prop[MU][i];
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
void determination_of_activity_flow(BOUND* sosedb, integer maxbound, 
									integer ls, integer lw, WALL* w, 
									bool &bactiv, bool &bpressure_fix, bool &bLR1free) {
	integer i=0; // счётчик цикла for
	bool bpfix=true; // по умолчанию давление нужно фиксировать в одной точке

	// 10.02.2017
	//bpfix = false; так делать нельзя т.к. решение развалится.
	
	for (i=0; i<maxbound; i++) {
		if ((sosedb[i].MCB>=ls) && (sosedb[i].MCB<(ls+lw))) {
			if (w[sosedb[i].MCB - ls].bpressure || w[sosedb[i].MCB - ls].bopening) {
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
    for (i=0; i<maxbound; i++) {
		if ((sosedb[i].MCB>=ls) && (sosedb[i].MCB<(ls+lw))) {
            if (!w[sosedb[i].MCB-ls].bpressure) {
				if ((fabs(w[sosedb[i].MCB-ls].Vx)>admission) || (fabs(w[sosedb[i].MCB-ls].Vy)>admission) || (fabs(w[sosedb[i].MCB-ls].Vz)>admission)) bnonzeromag=true;
			}
		}
	}
	if (!bactiv) bactiv=bnonzeromag;
	*/

    bpressure_fix=bpfix;

} // determination_of_activity_flow

// выделение оперативной памяти для задачи теплопроводности.
void allocation_memory_temp(doublereal* &potent, doublereal** &total_deformation, 
	equation3D* &sl, equation3D_bon* &slb, BOUND* sosedb, 
	integer maxelm, integer maxbound, integer ls, integer lw, WALL* w,
	doublereal operatingtemperature) {

    // выделение памяти под искомые полевые величины.
	potent=NULL;
    potent=new doublereal[maxelm+maxbound];
	if (potent==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for potent constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}

	total_deformation = NULL;
	total_deformation = new doublereal*[4];	
	if (total_deformation == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for total_deformation constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer i_1 = 0; i_1 < 4; i_1++) {
		total_deformation[i_1] = NULL;
	}
	for (integer i_1 = 0; i_1 < 4; i_1++) {
		total_deformation[i_1] = new doublereal[maxelm + maxbound];
		if (total_deformation[i_1] == NULL) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem : not enough memory on your equipment for total_deformation[%lld] constr struct...\n",i_1);
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}

    doublereal Tamb=0.0;
    bool bTamb=false;
       // Давление должно быть инициализировано с учётом граничных условий.
       if (lw>0) {
          Tamb=1.0e+20;
	      for (integer i=0; i<maxbound; i++) {
		      if (((sosedb[i].MCB<(ls+lw))&&(sosedb[i].MCB>=ls)) && (w[sosedb[i].MCB-ls].ifamily==1)) {
			      // на границе задана температура
			      Tamb=fmin(Tamb,w[sosedb[i].MCB-ls].Tamb);
                  bTamb=true;
		      }
	      }
	   }

	   // Внимание инициализация значением operatingtemperature !!!.
       if (!bTamb) Tamb= operatingtemperature; // 0.0;

	// обнуление
	integer j=0;
	for (j = 0; j < (maxelm + maxbound); j++) {
		potent[j] = Tamb; // ноль градусов цельсия начальная температура.
		for (integer i_1 = 0; i_1 < 4; i_1++) {
			total_deformation[i_1][j] = 0.0; // нулевая полная деформация
		}
	}

	sl=NULL;
	sl = new equation3D[maxelm]; // коэффициенты матрицы СЛАУ для внутренних КО.
	if (sl==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for slau temperature constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}

	slb=NULL;
	slb = new equation3D_bon[maxbound]; // коэффициенты матрицы СЛАУ для граничных КО
	if (slb==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for slau boundary temperature constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
} // allocation_memory_temp

bool is_EXTERNAL_FLOW(TOCHKA* pa1, integer** nvtx1, integer maxelm1,
	doublereal* &x47, doublereal* &y47, doublereal* &z47, integer** nvtx2, integer maxelm2,
	doublereal eps_minx, doublereal eps_miny, doublereal eps_minz,
	doublereal eps_maxx, doublereal eps_maxy, doublereal eps_maxz) {

	doublereal xmin1 = 1.0e60;
	doublereal ymin1 = 1.0e60;
	doublereal zmin1 = 1.0e60;
	doublereal xmax1 = 1.0e60;
	doublereal ymax1 = 1.0e60;
	doublereal zmax1 = 1.0e60;
	for (integer i = 0; i < maxelm1; i++) {
		for (integer j = 0; j < 8; j++) {
			if (pa1[nvtx1[j][i]-1].x < xmin1) xmin1 = pa1[nvtx1[j][i]-1].x;
			if (pa1[nvtx1[j][i]-1].y < ymin1) ymin1 = pa1[nvtx1[j][i]-1].y;
			if (pa1[nvtx1[j][i]-1].z < zmin1) zmin1 = pa1[nvtx1[j][i]-1].z;

			if (pa1[nvtx1[j][i]-1].x > xmax1) xmax1 = pa1[nvtx1[j][i]-1].x;
			if (pa1[nvtx1[j][i]-1].y > ymax1) ymax1 = pa1[nvtx1[j][i]-1].y;
			if (pa1[nvtx1[j][i]-1].z > zmax1) zmax1 = pa1[nvtx1[j][i]-1].z;
		}
	}

	doublereal xmin2 = 1.0e60;
	doublereal ymin2 = 1.0e60;
	doublereal zmin2 = 1.0e60;
	doublereal xmax2 = 1.0e60;
	doublereal ymax2 = 1.0e60;
	doublereal zmax2 = 1.0e60;
	for (integer i = 0; i < maxelm2; i++) {
		for (integer j = 0; j < 8; j++) {
			if (x47[nvtx2[j][i]] < xmin2) xmin2 = x47[nvtx2[j][i]];
			if (y47[nvtx2[j][i]] < ymin2) ymin2 = y47[nvtx2[j][i]];
			if (z47[nvtx2[j][i]] < zmin2) zmin2 = z47[nvtx2[j][i]];

			if (x47[nvtx2[j][i]] > xmax2) xmax2 = x47[nvtx2[j][i]];
			if (y47[nvtx2[j][i]] > ymax2) ymax2 = y47[nvtx2[j][i]];
			if (z47[nvtx2[j][i]] > zmax2) zmax2 = z47[nvtx2[j][i]];
		}
	}

	return ((fabs(xmin1-xmin2)<eps_minx)&&(fabs(ymin1 - ymin2)<eps_miny)&&(fabs(zmin1 - zmin2)<eps_minz)&&(fabs(xmax1 - xmax2)<eps_maxx)&&(fabs(ymax1 - ymax2)<eps_maxy)&&(fabs(zmax1 - zmax2)<eps_maxz));
}//is_EXTERNAL_FLOW


  // возвращает скорректированный массовый поток.
  // скоректированный массовый поток mf ВЫЧИСЛЯЕТСЯ на основе использования
  // скорректированной скорости или просто скорости на основе простейшей интерполляции.
  // Это нужно для отдельного решения уравнения конвекции-диффузии где задана пользовательская скорость,
  // никакой поправки Рхи-Чоу просто интерполляция.
  // 26.03.2017 15.09.2018 Декларация для использования уже здесь. Реализация в файле pamendment3.c.
void return_calc_correct_mass_flux_only_interpolation(integer iP, doublereal** potent, TOCHKA* pa, doublereal** prop, doublereal** prop_b,
	integer** nvtx, ALICE_PARTITION** sosedi, integer maxelm,
	doublereal* &mfcurrentretune, BOUND* &sosedb, integer &ls, integer &lw);

// выделение оперативной памяти для задачи гидродинамики.
void allocation_memory_flow(doublereal** &potent, equation3D** &sl, equation3D_bon** &slb,
	BOUND* sosedb, integer maxelm, integer maxbound,
	doublereal* &alpha, integer ls, integer lw, WALL* w,
	doublereal dgx, doublereal dgy, doublereal dgz, integer** nvtx,
	TOCHKA* pa, doublereal** prop, ALICE_PARTITION** sosedi,
	doublereal eps_minx, doublereal eps_miny, doublereal eps_minz,
	doublereal eps_maxx, doublereal eps_maxy, doublereal eps_maxz) {

	// выделение памяти под искомые полевые величины.
	potent=NULL;
    potent=new doublereal*[33];
	if (potent==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for potent constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	for (integer i=0; i<33; i++) potent[i]=NULL;
	for (integer i=0; i<33; i++) {
		potent[i]=new doublereal[maxelm+maxbound];
		if (potent[i]==NULL) {
	        // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for potent[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for potent[%d] constr struct...\n", i);
#endif
		    
		    printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
		    exit(1);
	    }
	}

	// обнуление
    for (integer i=0; i<33; i++) for (integer j=0; j<(maxelm+maxbound); j++) potent[i][j]=0.0;

	// Паскалевская составляющая давления при естественной конвекции 
	// не вычисляется из уравнений Навье-Стокса а просто добавляется при 
	// инициализации, хотя она влияет на скорости давления в Навье-Стоксе.
	TOCHKA p;
	doublereal minX1 = 1.0e20;
	doublereal minY1 = 1.0e20;
	doublereal minZ1 = 1.0e20;
	minX1 = 0.0;
	minY1 = 0.0;
	minZ1 = 0.0;
	for(integer iP = 0; iP < maxelm; iP++) {
		center_cord3D(iP, nvtx, pa, p,100);
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
		volume3D(iP, nvtx, pa, dx, dy, dz);
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
	for (integer iP = 0; iP < maxelm; iP++) {
		center_cord3D(iP, nvtx, pa, p,100);
		// вычисление размеров текущего контрольного объёма:
		doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
		volume3D(iP, nvtx, pa, dx, dy, dz);
		potent[PRESS][iP] += prop[RHO][iP] * (dgx*(p.x - minX1) + dgy*(p.y - minY1) + dgz*(p.z - minZ1));
		integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
		integer iE2, iN2, iT2, iW2, iS2, iB2;
		integer iE3, iN3, iT3, iW3, iS3, iB3;
		integer iE4, iN4, iT4, iW4, iS4, iB4;

		iE1 = sosedi[ESIDE][iP].iNODE1; iE2 = sosedi[ESIDE][iP].iNODE2; iE3 = sosedi[ESIDE][iP].iNODE3; iE4 = sosedi[ESIDE][iP].iNODE4;
		iN1 = sosedi[NSIDE][iP].iNODE1; iN2 = sosedi[NSIDE][iP].iNODE2; iN3 = sosedi[NSIDE][iP].iNODE3; iN4 = sosedi[NSIDE][iP].iNODE4;
		iT1 = sosedi[TSIDE][iP].iNODE1; iT2 = sosedi[TSIDE][iP].iNODE2; iT3 = sosedi[TSIDE][iP].iNODE3; iT4 = sosedi[TSIDE][iP].iNODE4;
		iW1 = sosedi[WSIDE][iP].iNODE1; iW2 = sosedi[WSIDE][iP].iNODE2; iW3 = sosedi[WSIDE][iP].iNODE3; iW4 = sosedi[WSIDE][iP].iNODE4;
		iS1 = sosedi[SSIDE][iP].iNODE1; iS2 = sosedi[SSIDE][iP].iNODE2; iS3 = sosedi[SSIDE][iP].iNODE3; iS4 = sosedi[SSIDE][iP].iNODE4;
		iB1 = sosedi[BSIDE][iP].iNODE1; iB2 = sosedi[BSIDE][iP].iNODE2; iB3 = sosedi[BSIDE][iP].iNODE3; iB4 = sosedi[BSIDE][iP].iNODE4;

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
	for (integer iP = 0; iP < maxelm+maxbound; iP++) {
		//potent[PRESS][iP] *= -1.0;
		// ANSYS Icepak не учитывает гидростатический перепад давления.
		// Для проверки этой гипотезы здесь произведено зануление.
		potent[PRESS][iP] = 0.0;
	}


	/*
	// Давление должно быть инициализировано с учётом граничных условий.
	for (i=0; i<maxbound; i++) {
		if (((sosedb[i].MCB<(ls+lw))&&(sosedb[i].MCB>=ls)) && (w[sosedb[i].MCB-ls].bpressure)) {
			// на границе задано давление
			potent[PRESS][sosedb[i].iB]=w[sosedb[i].MCB-ls].P;
		}
	}
	*/


	// Загрузка распределения начальной скорости.
	errno_t err_inicialization_data=0;
	FILE* fp_inicialization_data;	
#ifdef MINGW_COMPILLER
	fp_inicialization_data=fopen64("load.txt", "r");
#else
    err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif

	if (maxelm > 0) {
		// Только если есть жидкие ячейки.
		if (err_inicialization_data != 0) {
			// 04.04.2019
		    // Если мы решаем гидродинамическую задачу то считывать файл скоростей ненадо.

			// открытие неудачно или файл отсутствует.



			// Инициализация компонент скорости во внутренности расчётной области.
			// 26 марта 2017.
			for (integer i = 0; i < maxelm; i++) {
				potent[VX][i] = starting_speed_Vx;
				potent[VY][i] = starting_speed_Vy;
				potent[VZ][i] = starting_speed_Vz;
				// Скоректированнное поле скорости должно удовлетворять уравнению неразрывности.
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
			// При этом т.к. в тепловой модели сетка изменилась то осуществляется переинтерполляция с 
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
			// nvtx

			integer maxnode47 = 0;
			integer maxelm47 = 0;
			integer din47 = 0;
#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
			fscanf(fp_inicialization_data, "%lld", &din47);
			maxnode47 = din47;
			fscanf(fp_inicialization_data, "%lld", &din47);
#else
			fscanf(fp_inicialization_data, "%d", &din47);
			maxnode47 = din47;
			fscanf(fp_inicialization_data, "%d", &din47);
#endif	
#else
#if doubleintprecision == 1
			fscanf_s(fp_inicialization_data, "%lld", &din47);
			maxnode47 = din47;
			fscanf_s(fp_inicialization_data, "%lld", &din47);
#else
			fscanf_s(fp_inicialization_data, "%d", &din47);
			maxnode47 = din47;
			fscanf_s(fp_inicialization_data, "%d", &din47);
#endif	
#endif
			maxelm47 = din47;

			float fin47 = 0.0;
			doublereal *x47 = NULL, *y47 = NULL, *z47 = NULL, *Vx47 = NULL, *Vy47 = NULL, *Vz47 = NULL;
			integer** nvtx47 = NULL;
			// Везде нумерация с нуля.
			x47 = new doublereal[maxnode47];
			y47 = new doublereal[maxnode47];
			z47 = new doublereal[maxnode47];
			Vx47 = new doublereal[maxnode47];
			Vy47 = new doublereal[maxnode47];
			Vz47 = new doublereal[maxnode47];
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
				Vx47[i_47] = fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf(fp_inicialization_data, "%f", &fin47);
				Vy47[i_47] = fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf(fp_inicialization_data, "%f", &fin47);
				Vz47[i_47] = fin47;
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
				Vx47[i_47] = fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf_s(fp_inicialization_data, "%f", &fin47);
				Vy47[i_47] = fin47;
			}
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				fscanf_s(fp_inicialization_data, "%f", &fin47);
				Vz47[i_47] = fin47;
			}
#endif
			doublereal Speed_min = 1.0e60, Speed_max=-1.0e60;
			for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
				doublereal mag_speed = sqrt((Vx47[i_47])*(Vx47[i_47])+ (Vy47[i_47])*(Vy47[i_47])+ (Vz47[i_47])*(Vz47[i_47]));
				if (mag_speed < Speed_min) {
					Speed_min = mag_speed;
				}
				if (mag_speed > Speed_max) {
					Speed_max = mag_speed;
				}
			}
			printf("LOAD: Minimum speed = %e, Maximum speed = %e \n", Speed_min, Speed_max);

			nvtx47 = new integer*[8];
			for (integer i_47 = 0; i_47 < 8; i_47++) {
				nvtx47[i_47] = new integer[maxelm47];
			}

			for (integer j_47 = 0; j_47 < maxelm47; j_47++) {
				for (integer i_47 = 0; i_47 < 8; i_47++) {
#ifdef MINGW_COMPILLER
#if doubleintprecision == 1
					fscanf(fp_inicialization_data, "%lld", &din47);
#else
					fscanf(fp_inicialization_data, "%d", &din47);
#endif
#else
#if doubleintprecision == 1
					fscanf_s(fp_inicialization_data, "%lld", &din47);
#else
					fscanf_s(fp_inicialization_data, "%d", &din47);
#endif
#endif

					nvtx47[i_47][j_47] = din47;
				}
			}
			fclose(fp_inicialization_data);

			printf("load done...\n");
			printf("interpolation start...\n");

			// Инициализация компонент скорости во внутренности расчётной области.
			// 26 марта 2017.

			// Метод линейного порядка.
			doublereal min_x = 1e60;
			doublereal min_y = 1e60;
			doublereal min_z = 1e60;
			doublereal max_x = -1e60;
			doublereal max_y = -1e60;
			doublereal max_z = -1e60;

			for (integer i_48 = 0; i_48 < maxelm47; i_48++) {
				for (integer i_47 = 0; i_47 < 8; i_47++) {
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
			doublereal avgx_SH = avgx_separator;
			{
				doublereal delta = 1.0e60;
				for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
					if (fabs(avgx_separator - x47[i_47]) < delta) {
						delta = fabs(avgx_SH - x47[i_47]);
						avgx_SH = x47[i_47];
					}
				}
			}
			avgx_separator = avgx_SH;

			// Ищем сеточный узел максимально близкий к среднему значению.
			doublereal avgy_separator = 0.5*(min_y + max_y);
			doublereal avgy_SH = avgy_separator;
			{
				doublereal delta = 1.0e60;
				for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
					if (fabs(avgy_separator - y47[i_47]) < delta) {
						delta = fabs(avgy_SH - y47[i_47]);
						avgy_SH = y47[i_47];
					}
				}
			}
			avgy_separator = avgy_SH;

			// Ищем сеточный узел максимально близкий к среднему значению.
			doublereal avgz_separator = 0.5*(min_z + max_z);
			doublereal avgz_SH = avgz_separator;
			{
				doublereal delta = 1.0e60;
				for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
					if (fabs(avgz_separator - z47[i_47]) < delta) {
						delta = fabs(avgz_SH - z47[i_47]);
						avgz_SH = z47[i_47];
					}
				}
			}
			avgz_separator = avgz_SH;

			min_x = 1.05*fabs(max_x - min_x);
			if (min_x < 1.0e-40) {
				min_x = 1.05*fabs(max_x);
			}
			min_y = 1.05*fabs(max_y - min_y);
			if (min_y < 1.0e-40) {
				min_y = 1.05*fabs(max_y);
			}
			min_z = 1.05*fabs(max_z - min_z);
			if (min_z < 1.0e-40) {
				min_z = 1.05*fabs(max_z);
			}


			TOCHKA** pointerlist = new TOCHKA*[maxelm];
			doublereal** rthdsd_Gauss = new doublereal*[maxelm];
			for (integer i_47 = 0; i_47 < maxelm; i_47++) {
				pointerlist[i_47] = new TOCHKA[8];
				rthdsd_Gauss[i_47] = new doublereal[8];
			}


			
			integer i_11[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
			integer*** oct_load1 = NULL;
			doublereal avgx[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			doublereal avgy[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			doublereal avgz[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			integer** i_22 = NULL;
			i_22 = new integer*[8];
			for (integer i_1 = 0; i_1 < 8; i_1++) {
				i_22[i_1] = new integer[8];
				for (integer i_2 = 0; i_2 < 8; i_2++) {
					i_22[i_1][i_2] = 0;
				}
			}

			bool bext_flow=is_EXTERNAL_FLOW(pa, nvtx, maxelm,
				x47, y47, z47, nvtx47, maxelm47, eps_minx, eps_miny, eps_minz, eps_maxx, eps_maxy, eps_maxz );

			if (bext_flow) printf("EXTERNAL FLOW ACELERATOR is true..\n");

			// Были ли вычислены компоненты скорости на АЛИС ? 22,03,2019
			bool b_load_has_been_saved_on_ALICE = !bext_flow;

			if (!b_load_has_been_saved_on_ALICE) {

				printf("maxelm CFD=%lld \n", maxelm47);	
				integer** oct_load = NULL;
				oct_load = new integer*[8];
				if (oct_load == NULL) {
					printf("ERROR ALLOCATION MEMORY: array oct_load is NULL.");
					system("PAUSE");
				}
				for (integer i = 0; i < 8; i++) {
					oct_load[i] = NULL;
					oct_load[i] = new integer[maxelm47];
					if (oct_load[i] == NULL) {
						printf("ERROR ALLOCATION MEMORY: array oct_load[%lld] is NULL.", i);
						system("PAUSE");
					}
					for (integer i_1 = 0; i_1 < maxelm47; i_1++) {
						oct_load[i][i_1] = -1;
					}
				}
				
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
				
				oct_load1 = new integer**[8];
				if (oct_load1 == NULL) {
					printf("ERROR ALLOCATION MEMORY: array oct_load1 is NULL.");
					system("PAUSE");
				}
				for (integer i = 0; i < 8; i++) {
					oct_load1[i] = NULL;
					oct_load1[i] = new integer*[8];
					if (oct_load1 == NULL) {
						printf("ERROR ALLOCATION MEMORY: array oct_load1 is NULL.");
						system("PAUSE");
					}
					for (integer j = 0; j < 8; j++) {
						oct_load1[i][j] = new integer[i_11[i]];
						if (oct_load1[i][j] == NULL) {
							printf("ERROR ALLOCATION MEMORY: array oct_load1[%lld][%lld] is NULL.", i, j);
							system("PAUSE");
						}

						for (integer i_1 = 0; i_1 < i_11[i]; i_1++) {
							oct_load1[i][j][i_1] = -1;
						}
					}
				}

				printf("balance octree level 1: %e %e %e %e %e %e %e %e\n", (double)(i_11[0]) / (double)(maxelm47 + 1), (double)(i_11[1]) / (double)(maxelm47 + 1), (double)(i_11[2]) / (double)(maxelm47 + 1), (double)(i_11[3]) / (double)(maxelm47 + 1), (double)(i_11[4]) / (double)(maxelm47 + 1), (double)(i_11[5]) / (double)(maxelm47 + 1), (double)(i_11[6]) / (double)(maxelm47 + 1), (double)(i_11[7]) / (double)(maxelm47 + 1));


				
				for (integer i = 0; i < 8; i++) {

					// Разбиение максимально близкое 
					// к спомобу построения octree дерева
					switch (i) {
					case 0: avgx[i] = 0.25*(min_x1+max_x1);
						avgy[i] = 0.25*(min_y1 + max_y1);
						avgz[i] = 0.25*(min_z1 + max_z1);
						break;
					case 1: avgx[i] = 0.25*(min_x1 + max_x1);
						avgy[i] = 0.25*(min_y1 + max_y1);
						avgz[i] = 0.75*(min_z1 + max_z1);
						break;
					case 2: avgx[i] = 0.25*(min_x1 + max_x1);
						avgy[i] = 0.75*(min_y1 + max_y1);
						avgz[i] = 0.25*(min_z1 + max_z1);
						break;
					case 3: avgx[i] = 0.25*(min_x1 + max_x1);
						avgy[i] = 0.75*(min_y1 + max_y1);
						avgz[i] = 0.75*(min_z1 + max_z1);
						break;
					case 4: avgx[i] = 0.75*(min_x1 + max_x1);
						avgy[i] = 0.25*(min_y1 + max_y1);
						avgz[i] = 0.25*(min_z1 + max_z1);
						break;
					case 5: avgx[i] = 0.75*(min_x1 + max_x1);
						avgy[i] = 0.25*(min_y1 + max_y1);
						avgz[i] = 0.75*(min_z1 + max_z1);
						break;
					case 6: avgx[i] = 0.75*(min_x1 + max_x1);
						avgy[i] = 0.75*(min_y1 + max_y1);
						avgz[i] = 0.25*(min_z1 + max_z1);
						break;
					case 7: avgx[i] = 0.75*(min_x1 + max_x1);
						avgy[i] = 0.75*(min_y1 + max_y1);
						avgz[i] = 0.75*(min_z1 + max_z1);
						break;
					}

					/*
					// Новое разбиение.
					for (integer i_1 = 0; i_1 < i_11[i]; i_1++) {
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
						doublereal delta = 1.0e60;
						for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
							if (fabs(avgx_separator1 - x47[i_47]) < delta) {
								delta = fabs(avgx_SH1 - x47[i_47]);
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
						doublereal delta = 1.0e60;
						for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
							if (fabs(avgy_separator1 - y47[i_47]) < delta) {
								delta = fabs(avgy_SH1 - y47[i_47]);
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
						doublereal delta = 1.0e60;
						for (integer i_47 = 0; i_47 < maxnode47; i_47++) {
							if (fabs(avgz_separator1 - z47[i_47]) < delta) {
								delta = fabs(avgz_SH1 - z47[i_47]);
								avgz_SH1 = z47[i_47];
							}
						}
					}
					avgz_separator1 = avgz_SH1;
					avgz[i] = avgz_separator1; // Мы точно попали на сеточную линию.

				}

				for (integer i = 0; i < 8; i++) {
					delete[] oct_load[i];
					oct_load[i] = NULL;
				}
				delete[] oct_load;
				oct_load = NULL;

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

			


			for (integer i = 0; i < maxelm; i++) {
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
				if ((xc47 >= min_x1) && (xc47 <= max_x1) && (yc47 >= min_y1) && (yc47 <= max_y1) &&
					(zc47 >= min_z1) && (zc47 <= max_z1))
				{

					if (b_load_has_been_saved_on_ALICE) {
						
						// 22.03.2019 Для АЛИС сетки мы не можем разбить на три уровня octree дерева так как 
						// границы ячеек сетки не вдоль координатных линий а уступами, ступеньками сильно неравномерные
						// отклоняющиеся от линии разбиения тов меньшую то в большую сторону.
						for (integer i_47 = 0; i_47 < maxelm47; i_47++) {
							if ((xc47 >= x47[nvtx47[0][i_47]]) && (xc47 <= x47[nvtx47[1][i_47]]) && (yc47 >= y47[nvtx47[0][i_47]]) && (yc47 <= y47[nvtx47[3][i_47]]) &&
								(zc47 >= z47[nvtx47[0][i_47]]) && (zc47 <= z47[nvtx47[4][i_47]]))
							{
								ifound = i_47;
								bfound = true;
								break;
							}
						}
						
					}
					else {
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

				}




				// МНК
				// Три раза, для каждой компоненты скорости.
				doublereal starting_speed_Vx_47 = 0.0;
				doublereal starting_speed_Vy_47 = 0.0;
				doublereal starting_speed_Vz_47 = 0.0;

				if (bfound) {

					//TOCHKA p;
					//center_cord3D(i, nvtx, pa, p, 100);
					p.x = p.x + min_x;
					p.y = p.y + min_y;
					p.z = p.z + min_z;

					// VX

					for (integer j = 0; j <= 7; j++) {
						TOCHKA p1;
						p1.x = x47[nvtx47[j][ifound]];
						p1.y = y47[nvtx47[j][ifound]];
						p1.z = z47[nvtx47[j][ifound]];
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
						rthdsd_Gauss[i][j] = Vx47[nvtx47[j][ifound]];
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

					starting_speed_Vx_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));

					// VY

					for (integer j = 0; j <= 7; j++) {
						rthdsd_Gauss[i][j] = Vy47[nvtx47[j][ifound]];
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

					starting_speed_Vy_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));


					// VZ

					for (integer j = 0; j <= 7; j++) {
						rthdsd_Gauss[i][j] = Vz47[nvtx47[j][ifound]];
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

					starting_speed_Vz_47 = (koefmatr[0] + koefmatr[1] * (p.x) + koefmatr[2] * (p.y) + koefmatr[3] * (p.z));

					for (integer j = 0; j <= 3; j++) {
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

				potent[VX][i] = starting_speed_Vx_47;
				potent[VY][i] = starting_speed_Vy_47;
				potent[VZ][i] = starting_speed_Vz_47;
				// Скоректированнное поле скорости должно удовлетворять уравнению неразрывности.
				potent[VXCOR][i] = starting_speed_Vx_47;
				potent[VYCOR][i] = starting_speed_Vy_47;
				potent[VZCOR][i] = starting_speed_Vz_47;
			}


			doublereal Speed_min1 = 1.0e60, Speed_max1 = -1.0e60;
			for (integer i_47 = 0; i_47 < maxelm; i_47++) {
				doublereal mag_speed = sqrt((potent[VX][i_47])*(potent[VX][i_47]) + (potent[VY][i_47])*(potent[VY][i_47]) + (potent[VZ][i_47])*(potent[VZ][i_47]));
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
				potent[VX][i_47] *= Speed_max / Speed_max1;
				potent[VY][i_47] *= Speed_max / Speed_max1;
				potent[VZ][i_47] *= Speed_max / Speed_max1;
			}

			for (integer i_1 = 0; i_1 < 8; i_1++) {
				delete[] i_22[i_1];
			}
			delete[] i_22;

			if (oct_load1 != NULL) {
				for (integer i = 0; i < 8; i++) {
					for (integer j = 0; j < 8; j++) {
						delete[] oct_load1[i][j];
						oct_load1[i][j] = NULL;
					}
				}

				for (integer i = 0; i < 8; i++) {
					delete[] oct_load1[i];
					oct_load1[i] = NULL;
				}
				delete[] oct_load1;
				oct_load1 = NULL;
			}

			for (integer i = 0; i < maxelm; i++) {
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
			for (integer i_47 = 0; i_47 < 8; i_47++) {
				delete[] nvtx47[i_47];
			}
			delete[] nvtx47;

			printf("done.\n");
			//system("PAUSE");
		}
	}

	// Лучше начинать с нулевого поля скорости.
	// на границе должны быть выполнены граничные условия 8 мая 2013г (revised 2 апреля 2019г).
	// Скорость тоже должна быть инициализирована с учётом граничных условий Дирихле.
	for (integer i = 0; i < maxbound; i++) {		
		
		if ((sosedb[i].MCB < (ls + lw)) && (sosedb[i].MCB >= ls)) {

			// Добавил УСЛОВИЕ СИММЕТРИИ 02.04.2019.

			if ((w[sosedb[i].MCB - ls].bpressure) || (w[sosedb[i].MCB - ls].bopening)) {
				// Выходная граница потока на которой стоит условие Неймана.
				// Присваиваем скорость из ближайшей внутренней точки расчётной области.
				potent[VX][sosedb[i].iB] = potent[VX][sosedb[i].iI];
				potent[VY][sosedb[i].iB] = potent[VY][sosedb[i].iI];
				potent[VZ][sosedb[i].iB] = potent[VZ][sosedb[i].iI];
				// Скоректированнное поле скорости должно удовлетворять уравнению неразрывности.
				potent[VXCOR][sosedb[i].iB] = potent[VX][sosedb[i].iI];
				potent[VYCOR][sosedb[i].iB] = potent[VY][sosedb[i].iI];
				potent[VZCOR][sosedb[i].iB] = potent[VZ][sosedb[i].iI];
				//printf("Vx=%e, Vy=%e, Vz=%e\n",w[sosedb[i].MCB-ls].Vx, w[sosedb[i].MCB-ls].Vy, w[sosedb[i].MCB-ls].Vz);
				//system("PAUSE"); // debug
			}
			else if (w[sosedb[i].MCB - ls].bsymmetry) {
				// Тангенсальные компоненты скорости условие Неймана,
				// нормальная компонента скорости равна нулю.

				switch (w[sosedb[i].MCB - ls].iPlane) {
				case XY :
					potent[VX][sosedb[i].iB] = potent[VX][sosedb[i].iI];
					potent[VY][sosedb[i].iB] = potent[VY][sosedb[i].iI];
					potent[VXCOR][sosedb[i].iB] = potent[VX][sosedb[i].iI];
					potent[VYCOR][sosedb[i].iB] = potent[VY][sosedb[i].iI];
					potent[VZ][sosedb[i].iB] = 0.0;
					potent[VZCOR][sosedb[i].iB] = 0.0;
					break;
				case XZ:
					potent[VY][sosedb[i].iB] = 0.0;
					potent[VYCOR][sosedb[i].iB] = 0.0;
					potent[VX][sosedb[i].iB] = potent[VX][sosedb[i].iI];
					potent[VZ][sosedb[i].iB] = potent[VZ][sosedb[i].iI];
					potent[VXCOR][sosedb[i].iB] = potent[VX][sosedb[i].iI];
					potent[VZCOR][sosedb[i].iB] = potent[VZ][sosedb[i].iI];
					break;
				case YZ: 
					potent[VX][sosedb[i].iB] = 0.0;
					potent[VXCOR][sosedb[i].iB] = 0.0;
					potent[VY][sosedb[i].iB] = potent[VY][sosedb[i].iI];
					potent[VZ][sosedb[i].iB] = potent[VZ][sosedb[i].iI];
					potent[VYCOR][sosedb[i].iB] = potent[VY][sosedb[i].iI];
					potent[VZCOR][sosedb[i].iB] = potent[VZ][sosedb[i].iI];
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
				potent[VX][sosedb[i].iB] = w[sosedb[i].MCB - ls].Vx;
				potent[VY][sosedb[i].iB] = w[sosedb[i].MCB - ls].Vy;
				potent[VZ][sosedb[i].iB] = w[sosedb[i].MCB - ls].Vz;
				// Скоректированнное поле скорости должно удовлетворять уравнению неразрывности.
				potent[VXCOR][sosedb[i].iB] = w[sosedb[i].MCB - ls].Vx;
				potent[VYCOR][sosedb[i].iB] = w[sosedb[i].MCB - ls].Vy;
				potent[VZCOR][sosedb[i].iB] = w[sosedb[i].MCB - ls].Vz;
				//printf("Vx=%e, Vy=%e, Vz=%e\n",w[sosedb[i].MCB-ls].Vx, w[sosedb[i].MCB-ls].Vy, w[sosedb[i].MCB-ls].Vz);
				//system("PAUSE"); // debug
			}
		}
		else {
			// Абсолютно твёрдая неподвижная стенка на которой выполнено условие прилипания,
			// скорость равна нулю.
			// 02.04.2019
			potent[VX][sosedb[i].iB] = 0.0;
			potent[VY][sosedb[i].iB] = 0.0;
			potent[VZ][sosedb[i].iB] = 0.0;
			// Скоректированнное поле скорости должно удовлетворять уравнению неразрывности.
			potent[VXCOR][sosedb[i].iB] = 0.0;
			potent[VYCOR][sosedb[i].iB] = 0.0;
			potent[VZCOR][sosedb[i].iB] = 0.0;
		}
		
	}
	//*/
	 
	alpha=NULL;
	alpha = new doublereal[4];
	if (alpha==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for alpha constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	// Патанкар Fluent Гаврилов Андрей
    alpha[VX]=0.7; // 0.5 0.7 0.8
    alpha[VY]=0.7;
    alpha[VZ]=0.7;
    alpha[PRESS]=0.3; //0.8 0.3 (0.2, но вообще говоря он использует SIMPLEC по умолчанию).

   // коэффициенты матрицы СЛАУ для внутренних КО.
   sl=NULL;
   sl=new equation3D*[5];
   if (sl==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for slau flow constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
   
   for (integer i=0; i<5; i++) {
	   sl[i]=NULL;
   }

   for (integer i=0; i<5; i++) {
	   switch (i) {
		   case VX : sl[VX]=new equation3D[maxelm]; 
			         if (sl[VX]==NULL) {
	                     // недостаточно памяти на данном оборудовании.
		                 printf("Problem : not enough memory on your equipment for slau[VX] constr struct...\n");
		                 printf("Please any key to exit...\n");
						// system("PAUSE");
						 system("pause");
		                 exit(1);
	                 }
			         break;
		   case VY : sl[VY]=new equation3D[maxelm];
			         if (sl[VY]==NULL) {
	                     // недостаточно памяти на данном оборудовании.
		                 printf("Problem : not enough memory on your equipment for slau[VY] constr struct...\n");
		                 printf("Please any key to exit...\n");
						 //system("PAUSE");
						 system("pause");
		                 exit(1);
	                 } 
			         break;
		   case VZ : sl[VZ]=new equation3D[maxelm]; 
			         if (sl[VZ]==NULL) {
	                     // недостаточно памяти на данном оборудовании.
		                 printf("Problem : not enough memory on your equipment for slau[VZ] constr struct...\n");
		                 printf("Please any key to exit...\n");
						 //system("PAUSE");
						 system("pause");
		                 exit(1);
	                 }
			         break;
		   case PRESS : sl[PRESS]=NULL; break;
		   case PAM : sl[PAM]=new equation3D[maxelm];
			          if (sl[PAM]==NULL) {
	                     // недостаточно памяти на данном оборудовании.
		                 printf("Problem : not enough memory on your equipment for slau[PAM] constr struct...\n");
		                 printf("Please any key to exit...\n");
						 //system("PAUSE");
						 system("pause");
		                 exit(1);
	                  }
			          break;
	   }
   }

   // коэффициенты матрицы СЛАУ для граничных КО
   slb=NULL;
   slb = new equation3D_bon*[5];
   if (slb==NULL) {
	    // недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for slau boundary constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
   }

   for (integer i=0; i<5; i++) {
	   slb[i]=NULL;
   }

   for (integer i=0; i<5; i++) {
	   switch (i) {
		   case VX : slb[VX]=new equation3D_bon[maxbound];
			         if (slb[VX]==NULL) {
	                     // недостаточно памяти на данном оборудовании.
		                 printf("Problem : not enough memory on your equipment for slau_bon[VX] constr struct...\n");
		                 printf("Please any key to exit...\n");
						 //system("PAUSE");
						 system("pause");
		                 exit(1);
	                 } 
			         break;
		   case VY : slb[VY]=new equation3D_bon[maxbound]; 
			         if (slb[VY]==NULL) {
	                     // недостаточно памяти на данном оборудовании.
		                 printf("Problem : not enough memory on your equipment for slau_bon[VY] constr struct...\n");
		                 printf("Please any key to exit...\n");
						 //system("PAUSE");
						 system("pause");
		                 exit(1);
	                 }
			         break;
		   case VZ : slb[VZ]=new equation3D_bon[maxbound]; 
			         if (slb[VZ]==NULL) {
	                     // недостаточно памяти на данном оборудовании.
		                 printf("Problem : not enough memory on your equipment for slau_bon[VZ] constr struct...\n");
		                 printf("Please any key to exit...\n");
						 //system("PAUSE");
						 system("pause");
		                 exit(1);
	                 }  
			         break;
		   case PRESS : slb[PRESS]=NULL; break;
		   case PAM : slb[PAM]=new equation3D_bon[maxbound]; 
			         if (slb[PAM]==NULL) {
	                     // недостаточно памяти на данном оборудовании.
		                 printf("Problem : not enough memory on your equipment for slau_bon[PAM] constr struct...\n");
		                 printf("Please any key to exit...\n");
						 //system("PAUSE");
						 system("pause");
		                 exit(1);
	                 }
			         break;
	   }
   }

   
				
  
} // allocation_memory_flow


  // выделение оперативной памяти для задачи гидродинамики.
void allocation_memory_flow_2(
	doublereal** &diag_coef,  integer maxelm, integer maxbound,
	 doublereal* &SInvariantStrainRateTensor,
	doublereal** &mf) {

	// диагональные коэффициенты компонент скорости 
	// необходимо хранить для вычисления поправки Рхи-Чоу.
	diag_coef = NULL;
	diag_coef = new doublereal*[3];
	if (diag_coef == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for diag_coef constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}

	for (integer i = 0; i<3; i++) {
		diag_coef[i] = NULL;
		diag_coef[i] = new doublereal[maxelm + maxbound];
		if (diag_coef[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for diag_coef[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for diag_coef[%d] constr struct...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// инициализация
	for (integer i = 0; i<3; i++) for (integer j = 0; j<(maxelm + maxbound); j++) diag_coef[i][j] = 1.0;

	// Выделение памяти под S инвариант тензора скоростей-деформаций.
	SInvariantStrainRateTensor = NULL;
	SInvariantStrainRateTensor = new doublereal[maxelm + maxbound];
	if (SInvariantStrainRateTensor == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for SInvariantStrainRateTensor constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}
	// инициализация S инварианта strain-rate тензора:
	for (integer i = 0; i<(maxelm + maxbound); i++) {
		SInvariantStrainRateTensor[i] = 0.0; // инициализация
	}

	mf = NULL;
	mf = new doublereal*[maxelm];
	if (mf == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for mf constr struct...\n");
		printf("Please any key to exit...\n");
		//system("PAUSE");
		system("pause");
		exit(1);
	}

	for (integer i = 0; i<maxelm; i++) {
		mf[i] = NULL;
		mf[i] = new doublereal[6];
		if (mf[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for mf[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for mf[%d] constr struct...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	// инициализация :
	for (integer i = 0; i<maxelm; i++) {
		for (integer j = 0; j<6; j++) {
			mf[i][j] = 0.0;
		}
	}

}

// создаёт связи между контрольными объёмами для графической 
// визуализации.
// универсальна: подходит и для температуры и для течения
void constr_nvtxcell(integer* evt, BOUND* sosedb, integer maxbound, integer maxelm, bool bextendedprint, integer** &nvtxcell, integer &ncell,
	integer inx, integer iny, integer inz, TOCKA_INT* &tck_int_list) {

	// Если bextendedprinteger = true то мы имеем дело с расширенной печатью включая граничные значения узлов.

	// для графической визуализации
	// nvtxcell - связи для контрольных объёмов.
	// ncell - количество связей для контрольных объёмов.
	ncell = 0;
	// ncell > maxelm.
	integer i1, i2, i3, i4, i5, i6, i7, i8;
	//for (i = 0; i<(inx - 1); i++) for (j = 0; j<(iny - 1); j++) for (k = 0; k<(inz - 1); k++) {
	// Сокращает число просмотров на больших моделях.
	for (integer iscan = 0; iscan < maxelm; iscan++) {

		integer i = tck_int_list[iscan].i;
		integer j = tck_int_list[iscan].j;
		integer k = tck_int_list[iscan].k;

		if ((i < (inx - 1)) && (j < (iny - 1)) && (k < (inz - 1)))
		{

			// Подсчёт количества связей между внутренними КО в цикле
			i1 = i + j * inx + k * inx*iny;
		i2 = (i + 1) + j * inx + k * inx*iny;
		i3 = (i + 1) + (j + 1)*inx + k * inx*iny;
		i4 = i + (j + 1)*inx + k * inx*iny;
		i5 = i + j * inx + (k + 1)*inx*iny;
		i6 = (i + 1) + j * inx + (k + 1)*inx*iny;
		i7 = (i + 1) + (j + 1)*inx + (k + 1)*inx*iny;
		i8 = i + (j + 1)*inx + (k + 1)*inx*iny;
		/*i1 = i + j*inx + (k+1)*inx*iny;
		i2 = (i + 1) + j*inx + (k+1)*inx*iny;
		i3 = (i + 1) + (j + 1)*inx + (k+1)*inx*iny;
		i4 = i + (j + 1)*inx + (k+1)*inx*iny;
		i5 = i + j*inx + (k + 2)*inx*iny;
		i6 = (i + 1) + j*inx + (k + 2)*inx*iny;
		i7 = (i + 1) + (j + 1)*inx + (k + 2)*inx*iny;
		i8 = i + (j + 1)*inx + (k + 2)*inx*iny;
		*/
		if ((evt[i1] > 0) && (evt[i2] > 0) && (evt[i3] > 0) && (evt[i4] > 0) && (evt[i5] > 0) && (evt[i6] > 0) && (evt[i7] > 0) && (evt[i8] > 0)) {
			ncell++;

			if (bextendedprint) {

				// сканируем каждую внутреннюю нормаль в поисках граничных 
				// плоских бесконечно-тонких объёмов на которых вытавлено граничное условие.
				integer inorm = ESIDE;
				bool bfound1 = false, bfound2 = false, bfound3 = false, bfound4 = false;

				// поиск.
				for (integer isel = 0; isel < maxbound; isel++) {
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i1] - 1)) {
						bfound1 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i4] - 1)) {
						bfound2 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i5] - 1)) {
						bfound3 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i8] - 1)) {
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

				inorm = WSIDE;
				// поиск.
				for (integer isel = 0; isel < maxbound; isel++) {
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i2] - 1)) {
						bfound1 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i3] - 1)) {
						bfound2 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i6] - 1)) {
						bfound3 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i7] - 1)) {
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

				inorm = NSIDE;
				// поиск.
				for (integer isel = 0; isel < maxbound; isel++) {
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i1] - 1)) {
						bfound1 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i2] - 1)) {
						bfound2 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i5] - 1)) {
						bfound3 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i6] - 1)) {
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
				inorm = SSIDE;
				// поиск.
				for (integer isel = 0; isel < maxbound; isel++) {
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i3] - 1)) {
						bfound1 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i4] - 1)) {
						bfound2 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i7] - 1)) {
						bfound3 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i8] - 1)) {
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
				inorm = TSIDE;
				// поиск.
				for (integer isel = 0; isel < maxbound; isel++) {
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i1] - 1)) {
						bfound1 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i2] - 1)) {
						bfound2 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i3] - 1)) {
						bfound3 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i4] - 1)) {
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
				inorm = BSIDE;
				// поиск.
				for (integer isel = 0; isel < maxbound; isel++) {
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i5] - 1)) {
						bfound1 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i6] - 1)) {
						bfound2 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i7] - 1)) {
						bfound3 = true;
					}
					if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i8] - 1)) {
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
	// Выделение ОП 
	nvtxcell = NULL;
	nvtxcell = new integer*[8];
	if (nvtxcell == NULL) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem : not enough memory on your equipment for nvtxcell constr struct...\n");
		printf("Please any key to exit...\n");
		// system("PAUSE");
		system("pause");
		exit(1);
	}
	integer l = 0;
	for (integer i = 0; i<8; i++) {
		nvtxcell[i] = NULL;
		nvtxcell[i] = new integer[ncell];
		if (nvtxcell[i] == NULL) {
			// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for nvtxcell[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for nvtxcell[%d] constr struct...\n", i);
#endif
			
			printf("Please any key to exit...\n");
			//system("PAUSE");
			system("pause");
			exit(1);
		}
	}
	//for (i = 0; i<(inx - 1); i++) for (j = 0; j<(iny - 1); j++) for (k = 0; k<(inz - 1); k++) {
	// Сокращает число просмотров на больших моделях.
	for (integer iscan = 0; iscan < maxelm; iscan++) {

		integer i = tck_int_list[iscan].i;
		integer j = tck_int_list[iscan].j;
		integer k = tck_int_list[iscan].k;

		if ((i < (inx - 1)) && (j < (iny - 1)) && (k < (inz - 1)))
		{


			// Установка восьмиточечных связей для связывания КО.
			i1 = i + j * inx + k * inx*iny;
			i2 = (i + 1) + j * inx + k * inx*iny;
			i3 = (i + 1) + (j + 1)*inx + k * inx*iny;
			i4 = i + (j + 1)*inx + k * inx*iny;
			i5 = i + j * inx + (k + 1)*inx*iny;
			i6 = (i + 1) + j * inx + (k + 1)*inx*iny;
			i7 = (i + 1) + (j + 1)*inx + (k + 1)*inx*iny;
			i8 = i + (j + 1)*inx + (k + 1)*inx*iny;

			/*i1 = i + j*inx + (k + 1)*inx*iny;
			i2 = (i + 1) + j*inx + (k + 1)*inx*iny;
			i3 = (i + 1) + (j + 1)*inx + (k + 1)*inx*iny;
			i4 = i + (j + 1)*inx + (k + 1)*inx*iny;
			i5 = i + j*inx + (k + 2)*inx*iny;
			i6 = (i + 1) + j*inx + (k + 2)*inx*iny;
			i7 = (i + 1) + (j + 1)*inx + (k + 2)*inx*iny;
			i8 = i + (j + 1)*inx + (k + 2)*inx*iny;
			*/
			if ((evt[i1] > 0) && (evt[i2] > 0) && (evt[i3] > 0) && (evt[i4] > 0) && (evt[i5] > 0) && (evt[i6] > 0) && (evt[i7] > 0) && (evt[i8] > 0)) {
				nvtxcell[0][l] = evt[i1]; nvtxcell[1][l] = evt[i2]; nvtxcell[2][l] = evt[i3]; nvtxcell[3][l] = evt[i4];
				nvtxcell[4][l] = evt[i5]; nvtxcell[5][l] = evt[i6]; nvtxcell[6][l] = evt[i7]; nvtxcell[7][l] = evt[i8];
				l++;

				if (bextendedprint) {

					// сканируем каждую внутреннюю нормаль в поисках граничных 
					// плоских бесконечно-тонких объёмов на которых вытавлено граничное условие.
					integer inorm = ESIDE;
					bool bfound1 = false, bfound2 = false, bfound3 = false, bfound4 = false;
					integer im1 = -1, im2 = -1, im3 = -1, im4 = -1;


					// поиск.
					for (integer isel = 0; isel < maxbound; isel++) {
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i1] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
							//im1=sosedb[isel].iB+1;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i4] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
							//im2=sosedb[isel].iB+1;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i5] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
							//im3=sosedb[isel].iB+1;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i8] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
							//im4=sosedb[isel].iB+1;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						nvtxcell[0][l] = im1; nvtxcell[1][l] = evt[i2]; nvtxcell[2][l] = evt[i3]; nvtxcell[3][l] = im2;
						nvtxcell[4][l] = im3; nvtxcell[5][l] = evt[i6]; nvtxcell[6][l] = evt[i7]; nvtxcell[7][l] = im4;
						l++;

					}
					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;

					inorm = WSIDE;
					// поиск.
					for (integer isel = 0; isel < maxbound; isel++) {
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i2] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
							//im1=sosedb[isel].iB+1;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i3] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
							//im2=sosedb[isel].iB+1;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i6] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i7] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						nvtxcell[0][l] = evt[i1]; nvtxcell[1][l] = im1; nvtxcell[2][l] = im2; nvtxcell[3][l] = evt[i4];
						nvtxcell[4][l] = evt[i5]; nvtxcell[5][l] = im3; nvtxcell[6][l] = im4; nvtxcell[7][l] = evt[i8];
						l++;
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;

					inorm = NSIDE;
					// поиск.
					for (integer isel = 0; isel < maxbound; isel++) {
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i1] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i2] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i5] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i6] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						nvtxcell[0][l] = im1; nvtxcell[1][l] = im2; nvtxcell[2][l] = evt[i3]; nvtxcell[3][l] = evt[i4];
						nvtxcell[4][l] = im3; nvtxcell[5][l] = im4; nvtxcell[6][l] = evt[i7]; nvtxcell[7][l] = evt[i8];
						l++;
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;
					inorm = SSIDE;
					// поиск.
					for (integer isel = 0; isel < maxbound; isel++) {
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i3] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i4] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i7] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i8] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						nvtxcell[0][l] = evt[i1]; nvtxcell[1][l] = evt[i2]; nvtxcell[2][l] = im1; nvtxcell[3][l] = im2;
						nvtxcell[4][l] = evt[i5]; nvtxcell[5][l] = evt[i6]; nvtxcell[6][l] = im3; nvtxcell[7][l] = im4;
						l++;
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;
					inorm = TSIDE;
					// поиск.
					for (integer isel = 0; isel < maxbound; isel++) {
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i1] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i2] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i3] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i4] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						nvtxcell[0][l] = im1; nvtxcell[1][l] = im2; nvtxcell[2][l] = im3; nvtxcell[3][l] = im4;
						nvtxcell[4][l] = evt[i5]; nvtxcell[5][l] = evt[i6]; nvtxcell[6][l] = evt[i7]; nvtxcell[7][l] = evt[i8];
						l++;
					}

					bfound1 = false;
					bfound2 = false;
					bfound3 = false;
					bfound4 = false;
					inorm = BSIDE;
					// поиск.
					for (integer isel = 0; isel < maxbound; isel++) {
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i5] - 1)) {
							bfound1 = true;
							im1 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i6] - 1)) {
							bfound2 = true;
							im2 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i7] - 1)) {
							bfound3 = true;
							im3 = isel + 1 + maxelm;
						}
						if ((sosedb[isel].Norm == inorm) && (sosedb[isel].iI == evt[i8] - 1)) {
							bfound4 = true;
							im4 = isel + 1 + maxelm;
						}
					}

					if (bfound1&&bfound2&&bfound3&&bfound4) {
						nvtxcell[0][l] = evt[i1]; nvtxcell[1][l] = evt[i2]; nvtxcell[2][l] = evt[i3]; nvtxcell[3][l] = evt[i4];
						nvtxcell[4][l] = im1; nvtxcell[5][l] = im2; nvtxcell[6][l] = im3; nvtxcell[7][l] = im4;
						l++;
					}

				} // bextendedprint

			}
		}
	}
} // constr_nvtxcell

// Вычисление расстояния до ближайшей стенки
void calcdistwall(FLOW &f, integer ls, integer lw, WALL* w) {
	// Вычисление расстояния до ближайшей стенки требуется во многих моделях турбулентности.
    f.rdistWall=NULL;
	f.rdistWall=new doublereal[f.maxelm+f.maxbound];
	if (f.rdistWall==NULL) {
	      // недостаточно памяти на данном оборудовании.
		  printf("Problem : not enough memory on your equipment for rdistWall constr struct...\n");
		  printf("Please any key to exit...\n");
		  //system("PAUSE");
		  system("pause");
		  exit(1);
	}
	for (integer i=0; i<(f.maxelm+f.maxbound); i++) f.rdistWall[i]=1e30; // присваиваем очень большое вещественное число.

	bool *bwall=NULL;
	bwall=new bool[f.maxbound]; 
	if (bwall==NULL) {
	      // недостаточно памяти на данном оборудовании.
		  printf("Problem : not enough memory on your equipment for bwall constr struct...\n");
		  printf("Please any key to exit...\n");
		  //system("PAUSE");
		  system("pause");
		  exit(1);
	}
	for (integer i=0; i<f.maxbound; i++) bwall[i]=false; // граница не принадлежит твёрдой стенке.
	doublereal **rposition=NULL;
	rposition=new doublereal*[3];
	if (rposition==NULL) {
	      // недостаточно памяти на данном оборудовании.
		  printf("Problem : not enough memory on your equipment for rposition constr struct...\n");
		  printf("Please any key to exit...\n");
		  //system("PAUSE");
		  system("pause");
		  exit(1);
	}
	for (integer i=0; i<3; i++) {
		rposition[i]=NULL;
		rposition[i]=new doublereal[f.maxbound];
		if (rposition[i]==NULL) {
	      // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
			printf("Problem : not enough memory on your equipment for rposition[%lld] constr struct...\n", i);
#else
			printf("Problem : not enough memory on your equipment for rposition[%d] constr struct...\n", i);
#endif
		  
		  printf("Please any key to exit...\n");
		  //system("PAUSE");
		  system("pause");
		  exit(1);
	    }
	}
	const integer X=0;
	const integer Y=1;
	const integer Z=2;
	for (integer i=0; i<3; i++) {
		for (integer j=0; j<(f.maxbound); j++) {
			rposition[i][j]=0.0; // инициализация.
		}
	}

	for (integer i=0; i<f.maxelm; i++) {
		// цикл по всем внутренним контрольным объёмам
		for (integer G=0; G<6; G++) {

			
			if (f.sosedi[G][i].iNODE1>=f.maxelm) {
				// это граничный узел
				integer inumber=f.sosedi[G][i].iNODE1-f.maxelm; // номер граничного узла
				//TOCHKA p; // координаты центра КО.
				//center_cord3D(i, f.nvtx, f.pa, p,100); // вычисление координат центра КО.
				// вычисление размеров текущего контрольного объёма:
	            //doublereal dx=0.0, dy=0.0, dz=0.0;// объём текущего контроольного объёма
	           // volume3D(i, f.nvtx, f.pa, dx, dy, dz);

				if ((f.sosedb[inumber].MCB==(ls+lw))||(f.sosedb[inumber].MCB<ls)) {
					bwall[inumber]=true;					
				}

				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && (!w[f.sosedb[inumber].MCB - ls].bsymmetry) && (!w[f.sosedb[inumber].MCB - ls].bopening) && (!w[f.sosedb[inumber].MCB - ls].bpressure)) {
					doublereal Magnitude=fabs(w[f.sosedb[inumber].MCB-ls].Vx)+fabs(w[f.sosedb[inumber].MCB-ls].Vy)+fabs(w[f.sosedb[inumber].MCB-ls].Vz);
					if (Magnitude<admission) {
						// Стенка на которой задана нулевая скорость
						// Следовательно это твёрдая стенка
						bwall[inumber]=true;
					} // Magnitude
				} // WALL

				/*
				switch (G) {
					case ESIDE : rposition[X][inumber]=p.x+0.5*dx;
						     rposition[Y][inumber]=p.y;
							 rposition[Z][inumber]=p.z;
						     break;
					case WSIDE: rposition[X][inumber]=p.x-0.5*dx;
						     rposition[Y][inumber]=p.y;
							 rposition[Z][inumber]=p.z;
						     break;
					case NSIDE : rposition[X][inumber]=p.x;
						     rposition[Y][inumber]=p.y+0.5*dy;
							 rposition[Z][inumber]=p.z;
						     break;
					case SSIDE : rposition[X][inumber]=p.x;
						     rposition[Y][inumber]=p.y-0.5*dy;
							 rposition[Z][inumber]=p.z;
						     break;
					case TSIDE : rposition[X][inumber]=p.x;
						     rposition[Y][inumber]=p.y;
							 rposition[Z][inumber]=p.z+0.5*dz;
						     break;
					case BSIDE : rposition[X][inumber]=p.x;
						     rposition[Y][inumber]=p.y;
							 rposition[Z][inumber]=p.z-0.5*dz;
						     break;
					}
					*/
					rposition[X][inumber]= f.sosedb[inumber].p_c.x;
					rposition[Y][inumber]= f.sosedb[inumber].p_c.y;
					rposition[Z][inumber]= f.sosedb[inumber].p_c.z;

					
			}
		
			if (f.sosedi[G][i].iNODE2 >= f.maxelm) {
				// это граничный узел
				integer inumber = f.sosedi[G][i].iNODE2 - f.maxelm; // номер граничного узла
				//TOCHKA p; // координаты центра КО.
				//center_cord3D(i, f.nvtx, f.pa, p, 100); // вычисление координат центра КО.
														// вычисление размеров текущего контрольного объёма:
				//doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
				//volume3D(i, f.nvtx, f.pa, dx, dy, dz);

				if ((f.sosedb[inumber].MCB == (ls + lw)) || (f.sosedb[inumber].MCB<ls)) {
					bwall[inumber] = true;
				}

				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && (!w[f.sosedb[inumber].MCB - ls].bsymmetry) && (!w[f.sosedb[inumber].MCB - ls].bopening) && (!w[f.sosedb[inumber].MCB - ls].bpressure)) {
					doublereal Magnitude = fabs(w[f.sosedb[inumber].MCB - ls].Vx) + fabs(w[f.sosedb[inumber].MCB - ls].Vy) + fabs(w[f.sosedb[inumber].MCB - ls].Vz);
					if (Magnitude<admission) {
						// Стенка на которой задана нулевая скорость
						// Следовательно это твёрдая стенка
						bwall[inumber] = true;
					} // Magnitude
				} // WALL

				rposition[X][inumber] = f.sosedb[inumber].p_c.x;
				rposition[Y][inumber] = f.sosedb[inumber].p_c.y;
				rposition[Z][inumber] = f.sosedb[inumber].p_c.z;
				//printf("%e %e %e\n", rposition[X][inumber], rposition[Y][inumber], rposition[Z][inumber]);
				//getchar();

			}
		
			if (f.sosedi[G][i].iNODE3 >= f.maxelm) {
				// это граничный узел
				integer inumber = f.sosedi[G][i].iNODE3 - f.maxelm; // номер граничного узла
																	//TOCHKA p; // координаты центра КО.
																	//center_cord3D(i, f.nvtx, f.pa, p, 100); // вычисление координат центра КО.
																	// вычисление размеров текущего контрольного объёма:
																	//doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
																	//volume3D(i, f.nvtx, f.pa, dx, dy, dz);

				if ((f.sosedb[inumber].MCB == (ls + lw)) || (f.sosedb[inumber].MCB<ls)) {
					bwall[inumber] = true;
				}

				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && (!w[f.sosedb[inumber].MCB - ls].bsymmetry) && (!w[f.sosedb[inumber].MCB - ls].bopening) && (!w[f.sosedb[inumber].MCB - ls].bpressure)) {
					doublereal Magnitude = fabs(w[f.sosedb[inumber].MCB - ls].Vx) + fabs(w[f.sosedb[inumber].MCB - ls].Vy) + fabs(w[f.sosedb[inumber].MCB - ls].Vz);
					if (Magnitude<admission) {
						// Стенка на которой задана нулевая скорость
						// Следовательно это твёрдая стенка
						bwall[inumber] = true;
					} // Magnitude
				} // WALL

				rposition[X][inumber] = f.sosedb[inumber].p_c.x;
				rposition[Y][inumber] = f.sosedb[inumber].p_c.y;
				rposition[Z][inumber] = f.sosedb[inumber].p_c.z;

				//printf("%e %e %e\n", rposition[X][inumber], rposition[Y][inumber], rposition[Z][inumber]);
				//getchar();

			}

			if (f.sosedi[G][i].iNODE4 >= f.maxelm) {
				// это граничный узел
				integer inumber = f.sosedi[G][i].iNODE4 - f.maxelm; // номер граничного узла
																	//TOCHKA p; // координаты центра КО.
																	//center_cord3D(i, f.nvtx, f.pa, p, 100); // вычисление координат центра КО.
																	// вычисление размеров текущего контрольного объёма:
																	//doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
																	//volume3D(i, f.nvtx, f.pa, dx, dy, dz);

				if ((f.sosedb[inumber].MCB == (ls + lw)) || (f.sosedb[inumber].MCB<ls)) {
					bwall[inumber] = true;
				}

				if ((f.sosedb[inumber].MCB<(ls + lw)) && (f.sosedb[inumber].MCB >= ls) && (!w[f.sosedb[inumber].MCB - ls].bsymmetry) && (!w[f.sosedb[inumber].MCB - ls].bopening) && (!w[f.sosedb[inumber].MCB - ls].bpressure)) {
					doublereal Magnitude = fabs(w[f.sosedb[inumber].MCB - ls].Vx) + fabs(w[f.sosedb[inumber].MCB - ls].Vy) + fabs(w[f.sosedb[inumber].MCB - ls].Vz);
					if (Magnitude<admission) {
						// Стенка на которой задана нулевая скорость
						// Следовательно это твёрдая стенка
						bwall[inumber] = true;
					} // Magnitude
				} // WALL

				rposition[X][inumber] = f.sosedb[inumber].p_c.x;
				rposition[Y][inumber] = f.sosedb[inumber].p_c.y;
				rposition[Z][inumber] = f.sosedb[inumber].p_c.z;

				//printf("%e %e %e\n", rposition[X][inumber], rposition[Y][inumber], rposition[Z][inumber]);
				//getchar();

			}
		
		}
	}

	// Непосредственное вычисление расстояния:
	// Только для внутренних КО.
	for (integer i=0; i<f.maxelm; i++) {
		TOCHKA p; // координаты центра КО.
		center_cord3D(i, f.nvtx, f.pa, p,100); // вычисление координат центра КО.
		for (integer j=0; j<f.maxbound; j++) {
			if (bwall[j]) {
				doublereal dist=sqrt((p.x-rposition[X][j])*(p.x-rposition[X][j])+(p.y-rposition[Y][j])*(p.y-rposition[Y][j])+(p.z-rposition[Z][j])*(p.z-rposition[Z][j]));
				if (dist<f.rdistWall[i]) f.rdistWall[i]=dist;
			}
		}
	}

	// Непосредственное вычисление расстояния для граничных КО:
	// принадлежащих входной границе, выходной границе или границе с условиями симметрии.
	for (integer i=f.maxelm; i<f.maxelm+f.maxbound; i++) {
		integer inumber=i-f.maxelm; // номер граничного узла
		if (bwall[inumber]) {
			// расстояние равно нулю
			f.rdistWall[i]=0.0;
		}
		else {
            for (integer j=0; j<f.maxbound; j++) {
		 	    if (bwall[j]) {
                    doublereal dist=sqrt((rposition[X][inumber]-rposition[X][j])*(rposition[X][inumber]-rposition[X][j])+(rposition[Y][inumber]-rposition[Y][j])*(rposition[Y][inumber]-rposition[Y][j])+(rposition[Z][inumber]-rposition[Z][j])*(rposition[Z][inumber]-rposition[Z][j]));
				    if (dist<f.rdistWall[i]) f.rdistWall[i]=dist;
			    }
			}
		}
	}

	// Вычисление толщины пограничного слоя в формуле Эскудиера
	f.rdistWallmax=-1e20; // очень маленькое значение
	for (integer i=0; i<(f.maxelm+f.maxbound); i++) {
		if (f.rdistWall[i]>f.rdistWallmax) f.rdistWallmax=f.rdistWall[i];
	}


	// Освобождение оперативной памяти.
	delete[] bwall;
	bwall=NULL;
	for (integer i=0; i<3; i++) {
		delete[] rposition[i];
		rposition[i]=NULL;
	}
	delete[] rposition;
	rposition=NULL;
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
		printf("f[0].maxelm=%lld\n", f[0].maxelm);
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
				for (integer i=0; i<flow_interior; i++) {
				    f[i].rdistWall=NULL; // инициализация
					f[i].iflowregime=LAMINAR; // LAMINAR
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
			// до каждой програмной жидкой зоны (минимальное расстояние).
			if (flow_interior==1) {
				// всё впорядке т.к. имеем только одну жидкую зону.
				
				if (eqin.fluidinfo[0].iflowregime==0) {
					// Ламинарный режим
					f[0].iflowregime=LAMINAR; // LAMINAR
				}
				else {
					if (eqin.fluidinfo[0].iturbmodel==0) {
                       f[0].iflowregime=ZEROEQMOD; // ZEROEQMOD
					}
					if (eqin.fluidinfo[0].iturbmodel==1) {
                       f[0].iflowregime=SMAGORINSKY; // SMAGORINSKY
					}
					if (eqin.fluidinfo[0].iturbmodel==2) {
                       f[0].iflowregime=RNG_LES; // RNG_LES
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

				// Вычисление вспомогательных величин для моделей турбулентности :
			    bool bdist=false; // вычислять ли расстояние до стенки.
			    switch (f[0].iflowregime) {
			      case LAMINAR : bdist=false; break;  // Ламинарное течение.
			      case ZEROEQMOD : bdist=true; break; // Zero Equation Turbulence model (RANS)
				  case SMAGORINSKY : bdist=true; break; // Smagorinsky model (LES) (содержит как опцию модель Германо 1991 года).
				  case RNG_LES : bdist=false; break; // RNG_LES Renormalization Group Theory (LES)
			    }
				f[0].rdistWall=NULL; // инициализация
			    if (bdist) {
				   // Вычисление расстояния до ближайшей стенки
				   printf("please wait! program calculate distance to the wall...\n");
				   if (!bCFXcalcdistancetotheWALL) {
					   // простой переборный алгоритм.
					   calcdistwall(f[0], ls, lw, w);
				   }
				   else {
					   // способ используемый в ANSYS CFX.
					   calcdistwallCFX(f[0], ls, lw, w);
				   }
			    }
			}
			else {
				// имеем две жидкие зоны или больше.
			    doublereal** distzone=NULL;
				// distzone[i][j] - кратчайшее расстояние от опорной точки i-ого жидкого отсека пользователя до j-ой распознаной программой жидкой зоны.
				// выделение памяти оперативной памяти
				distzone=new doublereal*[flow_interior];
				if (distzone==NULL) {
	                // недостаточно памяти на данном оборудовании.
		            printf("Problem : not enough memory on your equipment for distzone constr struct...\n");
		            printf("Please any key to exit...\n");
					//system("PAUSE");
					system("pause");
		            exit(1);
	            }
				for (i=0; i<flow_interior; i++) {
					distzone[i]=NULL;
					distzone[i]=new doublereal[flow_interior];
					if (distzone[i]==NULL) {
	                     // недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
						printf("Problem : not enough memory on your equipment for distzone[%lld] constr struct...\n", i);
#else
						printf("Problem : not enough memory on your equipment for distzone[%d] constr struct...\n", i);
#endif
		                
		                 printf("Please any key to exit...\n");
						 //system("PAUSE");
						 system("pause");
		                 exit(1);
	                }
				}
				// инициализация:
				for (i=0; i<flow_interior; i++) for (integer j=0; j<flow_interior; j++) distzone[i][j]=1e30; // присваиваем максимально большое число.
				// Вычисление расстояний:
				for (i=0; i<flow_interior; i++) {
					for (integer j=0; j<flow_interior; j++) {
						// только до внутренних контрольных объёмов
						TOCHKA p; // координаты центра КО.
						for (integer iP=0; iP<f[j].maxelm; iP++) {
							center_cord3D(iP, f[j].nvtx, f[j].pa, p,100); // вычисление координат центра КО.
							distzone[i][j]=fmin(distzone[i][j],sqrt((eqin.fluidinfo[i].xc-p.x)*(eqin.fluidinfo[i].xc-p.x)+(eqin.fluidinfo[i].yc-p.y)*(eqin.fluidinfo[i].yc-p.y)+(eqin.fluidinfo[i].zc-p.z)*(eqin.fluidinfo[i].zc-p.z)));
						}
					}
				}
				integer *key=NULL;
				key=new integer[flow_interior]; // ключ для сортировки
				if (key==NULL) {
	                // недостаточно памяти на данном оборудовании.
		            printf("Problem : not enough memory on your equipment for key constr struct...\n");
		            printf("Please any key to exit...\n");
					//system("PAUSE");
					system("pause");
		            exit(1);
	            }
				// инициализация.
				for (i=0; i<flow_interior; i++) key[i]=flow_interior; // присваиваем максимальное значение.
				for (i=0; i<flow_interior; i++) {
					// поиск минимума в i-ой строке distzone
					doublereal rmin=1e30;
					integer imin=flow_interior;
					for (integer j=0; j<flow_interior; j++) {
						if (distzone[i][j]<rmin) {
							rmin=distzone[i][j];
							imin=j;
						}
					}
					// минимальное расстояние в i-ой строке находится в ячейке imin
					key[i]=imin;
				}

				bool berror=false;
				for (i=0; i<flow_interior; i++) if (key[i]>=flow_interior) {
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
					for (i=0; i<numberOfPairs; i++) {
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
				for (i=0; i<(flow_interior-1); i++) {
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
				for (i=0; i<flow_interior; i++) {
					delete[] distzone[i];
					distzone[i]=NULL;
				}
				delete[] distzone;
				distzone=NULL;
				delete[] key;
				key=NULL;

				for (i=0; i<flow_interior; i++) {
					if (eqin.fluidinfo[i].iflowregime==0) {
					    // Ламинарный режим
					    f[i].iflowregime=LAMINAR; // LAMINAR
				    }
				    else {
					    if (eqin.fluidinfo[i].iturbmodel==0) {
                            f[i].iflowregime=ZEROEQMOD; // ZEROEQMOD
					    }
						if (eqin.fluidinfo[i].iturbmodel==1) {
                            f[i].iflowregime=SMAGORINSKY; // SMAGORINSKY	
					    }
                        if (eqin.fluidinfo[i].iturbmodel==2) {
                            f[i].iflowregime=RNG_LES; // RNG_LES	
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

				    // Вычисление вспомогательных величин для моделей турбулентности :
			        bool bdist=false; // вычислять ли расстояние до стенки.
			        switch (f[i].iflowregime) {
			           case LAMINAR : bdist=false; break;
			           case ZEROEQMOD : bdist=true; break; // Zero Equation Model (RANS)
					   case SMAGORINSKY : bdist=true; break; // Smagorinsky model (LES) (содержит как опцию модель Германо 1991 года).
					   case RNG_LES : bdist=false; break; // RNG_LES Renormalization Group Theory (LES)
			        }
					f[i].rdistWall=NULL; // инициализация
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
void constr_link_on_surface_for_radiation_model(integer maxelm, integer* &whot_is_block, ALICE_PARTITION** sosedi,
	                                            integer** nvtx, TOCHKA* pa,  BLOCK* &b, integer lb)
{
	for (integer i = 0; i < lb; i++) {
		if (b[i].radiation.binternalRadiation) {
			integer itek = i;

			MY_PAIR* list = new MY_PAIR[maxelm];
			integer ic = 0;

			

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek) 
				{
					integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
					iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iW >= maxelm) {
						// граничный узел.
						list[ic].node1 = iP;
						list[ic].node21 = iW;
						list[ic].node22 = -1;
						list[ic].node23 = -1;
						list[ic].node24 = -1;
						//list[ic].dS = dy*dz; // Площадь ячейки.
						list[ic].dS1 = dy*dz;
						list[ic].dS2 = 0.0;
						list[ic].dS3 = 0.0;
						list[ic].dS4 = 0.0;
						ic++;
					}
					else
					{
						// Внутренняя грань, но соседи принадлежат разным блокам.
						if (itek != whot_is_block[iW]) {
							list[ic].node1 = iP;
							list[ic].node21 = iW;
							list[ic].node22 = -1;
							list[ic].node23 = -1;
							list[ic].node24 = -1;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							list[ic].dS1 = dy*dz;
							list[ic].dS2 = 0.0;
							list[ic].dS3 = 0.0;
							list[ic].dS4 = 0.0;
							ic++;
						}
					}

				}
			}

			b[itek].radiation.nodelistW = new MY_PAIR[ic];
			b[itek].radiation.nodelistWsize = ic;
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistW[j] = list[j];

			delete[] list;
			list = NULL;

			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
					integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
					iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iE >= maxelm) {
						list[ic].node1 = iP;
						list[ic].node21 = iE;
						list[ic].node22 = -1;
						list[ic].node23 = -1;
						list[ic].node24 = -1;
						//list[ic].dS = dy*dz; // Площадь ячейки.
						list[ic].dS1 = dy*dz;
						list[ic].dS2 = 0.0;
						list[ic].dS3 = 0.0;
						list[ic].dS4 = 0.0;
						ic++;
					}
					else
					{
						if (itek != whot_is_block[iE]) {
							list[ic].node1 = iP;
							list[ic].node21 = iE;
							list[ic].node22 = -1;
							list[ic].node23 = -1;
							list[ic].node24 = -1;
							//list[ic].dS = dy*dz; // Площадь ячейки.
							list[ic].dS1 = dy*dz;
							list[ic].dS2 = 0.0;
							list[ic].dS3 = 0.0;
							list[ic].dS4 = 0.0;
							ic++;
						}
					}

				}
			}

			b[itek].radiation.nodelistE = new MY_PAIR[ic];
			b[itek].radiation.nodelistEsize = ic;
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistE[j] = list[j];

			delete[] list;
			list = NULL;

			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
					integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
					iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iS >= maxelm) {
						list[ic].node1 = iP;
						list[ic].node21 = iS;
						list[ic].node22 = -1;
						list[ic].node23 = -1;
						list[ic].node24 = -1;
						//list[ic].dS = dx*dz; // Площадь ячейки.
						list[ic].dS1 = dx*dz;
						list[ic].dS2 = 0.0;
						list[ic].dS3 = 0.0;
						list[ic].dS4 = 0.0;
						ic++;
					}
					else
					{
						if (itek != whot_is_block[iS]) {
							list[ic].node1 = iP;
							list[ic].node21 = iS;
							list[ic].node22 = -1;
							list[ic].node23 = -1;
							list[ic].node24 = -1;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							list[ic].dS1 = dx*dz;
							list[ic].dS2 = 0.0;
							list[ic].dS3 = 0.0;
							list[ic].dS4 = 0.0;
							ic++;
						}
					}

				}
			}

			b[itek].radiation.nodelistS = new MY_PAIR[ic];
			b[itek].radiation.nodelistSsize = ic;
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistS[j] = list[j];

			delete[] list;
			list = NULL;

			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
					integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
					iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iN >= maxelm) {
						list[ic].node1 = iP;
						list[ic].node21 = iN;
						list[ic].node22 = -1;
						list[ic].node23 = -1;
						list[ic].node24 = -1;
						//list[ic].dS = dx*dz; // Площадь ячейки.
						list[ic].dS1 = dx*dz;
						list[ic].dS2 = 0.0;
						list[ic].dS3 = 0.0;
						list[ic].dS4 = 0.0;
						ic++;
					}
					else
					{
						if (itek != whot_is_block[iN]) {
							list[ic].node1 = iP;
							list[ic].node21 = iN;
							list[ic].node22 = -1;
							list[ic].node23 = -1;
							list[ic].node24 = -1;
							//list[ic].dS = dx*dz; // Площадь ячейки.
							list[ic].dS1 = dx*dz;
							list[ic].dS2 = 0.0;
							list[ic].dS3 = 0.0;
							list[ic].dS4 = 0.0;
							ic++;
						}
					}

				}
			}

			b[itek].radiation.nodelistN = new MY_PAIR[ic];
			b[itek].radiation.nodelistNsize = ic;
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistN[j] = list[j];

			delete[] list;
			list = NULL;


			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
					integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
					iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iB >= maxelm) {
						list[ic].node1 = iP;
						list[ic].node21 = iB;
						list[ic].node22 = -1;
						list[ic].node23 = -1;
						list[ic].node24 = -1;
						//list[ic].dS = dx*dy; // Площадь ячейки.
						list[ic].dS1 = dx*dy;
						list[ic].dS2 = 0.0;
						list[ic].dS3 = 0.0;
						list[ic].dS4 = 0.0;
						ic++;
					}
					else
					{
						if (itek != whot_is_block[iB]) {
							list[ic].node1 = iP;
							list[ic].node21 = iB;
							list[ic].node22 = -1;
							list[ic].node23 = -1;
							list[ic].node24 = -1;
							//list[ic].dS = dx*dy; // Площадь ячейки.
							list[ic].dS1 = dx*dy;
							list[ic].dS2 = 0.0;
							list[ic].dS3 = 0.0;
							list[ic].dS4 = 0.0;
							ic++;
						}
					}

				}
			}

			b[itek].radiation.nodelistB = new MY_PAIR[ic];
			b[itek].radiation.nodelistBsize = ic;
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistB[j] = list[j];

			delete[] list;
			list = NULL;

			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
					integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
					iE = sosedi[ESIDE][iP].iNODE1; iN = sosedi[NSIDE][iP].iNODE1; iT = sosedi[TSIDE][iP].iNODE1; iW = sosedi[WSIDE][iP].iNODE1; iS = sosedi[SSIDE][iP].iNODE1; iB = sosedi[BSIDE][iP].iNODE1;


					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
					volume3D(iP, nvtx, pa, dx, dy, dz);

					if (iT >= maxelm) {
						list[ic].node1 = iP;
						list[ic].node21 = iT;
						list[ic].node22 = -1;
						list[ic].node23 = -1;
						list[ic].node24 = -1;
						//list[ic].dS = dx*dy; // Площадь ячейки.
						list[ic].dS1 = dx*dy;
						list[ic].dS2 = 0.0;
						list[ic].dS3 = 0.0;
						list[ic].dS4 = 0.0;
						ic++;
					}
					else
					{
						if (itek != whot_is_block[iT]) {
							list[ic].node1 = iP;
							list[ic].node21 = iT;
							list[ic].node22 = -1;
							list[ic].node23 = -1;
							list[ic].node24 = -1;
							//list[ic].dS = dx*dy; // Площадь ячейки.
							list[ic].dS1 = dx*dy;
							list[ic].dS2 = 0.0;
							list[ic].dS3 = 0.0;
							list[ic].dS4 = 0.0;
							ic++;
						}
					}

				}
			}

			b[itek].radiation.nodelistT = new MY_PAIR[ic];
			b[itek].radiation.nodelistTsize = ic;
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistT[j] = list[j];

			delete[] list;
			list = NULL;


		}
	}
} // constr_link_on_surface_for_radiation_model


// Начало написания 14 сентября 2016. 15_06.
// Продолжение написания 20 сентября 2016. 14_45.
// Правка 12 июня 2019.
// Также необходимо менять код в модуле my_solver_v0_03.cpp,
// и в модуле my_elmatr_quad.cpp
// Создаёт необходимые связи для модели излучения внутри блока.
void constr_link_on_surface_for_radiation_model_alice(integer maxelm, BOUND* &sosedb, integer* &whot_is_block, ALICE_PARTITION** sosedi,
	integer** nvtx, TOCHKA* pa, BLOCK* &b, integer lb)
{
	for (integer i = 0; i < lb; i++) {
		if (b[i].radiation.binternalRadiation) {
			integer itek = i;

			MY_PAIR* list = new MY_PAIR[maxelm];
			integer ic = 0;



			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
					integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					integer iE2, iN2, iT2, iW2, iS2, iB2;
					integer iE3, iN3, iT3, iW3, iS3, iB3;
					integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					iE1 = sosedi[ESIDE][iP].iNODE1; iN1 = sosedi[NSIDE][iP].iNODE1; iT1 = sosedi[TSIDE][iP].iNODE1; 
					iW1 = sosedi[WSIDE][iP].iNODE1; iS1 = sosedi[SSIDE][iP].iNODE1; iB1 = sosedi[BSIDE][iP].iNODE1;

					iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
					iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;

					iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
					iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;

					iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
					iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;

					// вычисление размеров текущего контрольного объёма:
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS1 = sosedb[iW1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							// Внутренняя грань, но соседи принадлежат разным блокам.
							if (itek != whot_is_block[iW1]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS2 = sosedb[iW2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							// Внутренняя грань, но соседи принадлежат разным блокам.
							if (itek != whot_is_block[iW2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS3 = sosedb[iW3 - maxelm].dS; // Площадь ячейки.

							binc = true;
						}
						else
						{
							// Внутренняя грань, но соседи принадлежат разным блокам.
							if (itek != whot_is_block[iW3]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS4 = sosedb[iW4 - maxelm].dS; // Площадь ячейки.

							binc = true;
						}
						else
						{
							// Внутренняя грань, но соседи принадлежат разным блокам.
							if (itek != whot_is_block[iW4]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistW[j] = list[j];

			delete[] list;
			list = NULL;

			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
				

					integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					integer iE2, iN2, iT2, iW2, iS2, iB2;
					integer iE3, iN3, iT3, iW3, iS3, iB3;
					integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					iE1 = sosedi[ESIDE][iP].iNODE1; iN1 = sosedi[NSIDE][iP].iNODE1; iT1 = sosedi[TSIDE][iP].iNODE1;
					iW1 = sosedi[WSIDE][iP].iNODE1; iS1 = sosedi[SSIDE][iP].iNODE1; iB1 = sosedi[BSIDE][iP].iNODE1;

					iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
					iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;

					iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
					iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;

					iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
					iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;


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
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS1 = sosedb[iE1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iE1]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS2 = sosedb[iE2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iE2]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS3 = sosedb[iE3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iE3]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS4 = sosedb[iE4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iE4]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistE[j] = list[j];

			delete[] list;
			list = NULL;

			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
					
					integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					integer iE2, iN2, iT2, iW2, iS2, iB2;
					integer iE3, iN3, iT3, iW3, iS3, iB3;
					integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					iE1 = sosedi[ESIDE][iP].iNODE1; iN1 = sosedi[NSIDE][iP].iNODE1; iT1 = sosedi[TSIDE][iP].iNODE1;
					iW1 = sosedi[WSIDE][iP].iNODE1; iS1 = sosedi[SSIDE][iP].iNODE1; iB1 = sosedi[BSIDE][iP].iNODE1;

					iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
					iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;

					iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
					iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;

					iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
					iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;


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
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS1 = sosedb[iS1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iS1]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS2 = sosedb[iS2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iS2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS3 = sosedb[iS3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iS3]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS4 = sosedb[iS4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iS4]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistS[j] = list[j];

			delete[] list;
			list = NULL;

			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
					
					integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					integer iE2, iN2, iT2, iW2, iS2, iB2;
					integer iE3, iN3, iT3, iW3, iS3, iB3;
					integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					iE1 = sosedi[ESIDE][iP].iNODE1; iN1 = sosedi[NSIDE][iP].iNODE1; iT1 = sosedi[TSIDE][iP].iNODE1;
					iW1 = sosedi[WSIDE][iP].iNODE1; iS1 = sosedi[SSIDE][iP].iNODE1; iB1 = sosedi[BSIDE][iP].iNODE1;

					iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
					iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;

					iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
					iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;

					iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
					iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;


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
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS1 = sosedb[iN1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iN1]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS2 = sosedb[iN2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iN2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS3 = sosedb[iN3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iN3]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS4 = sosedb[iN4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iN4]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistN[j] = list[j];

			delete[] list;
			list = NULL;


			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
				
					integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					integer iE2, iN2, iT2, iW2, iS2, iB2;
					integer iE3, iN3, iT3, iW3, iS3, iB3;
					integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					iE1 = sosedi[ESIDE][iP].iNODE1; iN1 = sosedi[NSIDE][iP].iNODE1; iT1 = sosedi[TSIDE][iP].iNODE1;
					iW1 = sosedi[WSIDE][iP].iNODE1; iS1 = sosedi[SSIDE][iP].iNODE1; iB1 = sosedi[BSIDE][iP].iNODE1;

					iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
					iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;

					iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
					iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;

					iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
					iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;


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
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS1 = sosedb[iB1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iB1]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS2 = sosedb[iB2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iB2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS3 = sosedb[iB3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iB3]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS4 = sosedb[iB4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iB4]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistB[j] = list[j];

			delete[] list;
			list = NULL;

			list = new MY_PAIR[maxelm];
			ic = 0;

			for (integer iP = 0; iP < maxelm; iP++) {
				integer ib = whot_is_block[iP];
				if (ib == itek)
				{
				
					integer iE1, iN1, iT1, iW1, iS1, iB1; // номера соседних контрольных объёмов
					integer iE2, iN2, iT2, iW2, iS2, iB2;
					integer iE3, iN3, iT3, iW3, iS3, iB3;
					integer iE4, iN4, iT4, iW4, iS4, iB4;

					// -1 если узел не существует.
					// 0 .. maxelm-1 - строго внутренний узел.
					// maxelm .. maxelm + maxbound-1 - граничный узел.

					iE1 = sosedi[ESIDE][iP].iNODE1; iN1 = sosedi[NSIDE][iP].iNODE1; iT1 = sosedi[TSIDE][iP].iNODE1;
					iW1 = sosedi[WSIDE][iP].iNODE1; iS1 = sosedi[SSIDE][iP].iNODE1; iB1 = sosedi[BSIDE][iP].iNODE1;

					iE2 = sosedi[ESIDE][iP].iNODE2; iN2 = sosedi[NSIDE][iP].iNODE2; iT2 = sosedi[TSIDE][iP].iNODE2;
					iW2 = sosedi[WSIDE][iP].iNODE2; iS2 = sosedi[SSIDE][iP].iNODE2; iB2 = sosedi[BSIDE][iP].iNODE2;

					iE3 = sosedi[ESIDE][iP].iNODE3; iN3 = sosedi[NSIDE][iP].iNODE3; iT3 = sosedi[TSIDE][iP].iNODE3;
					iW3 = sosedi[WSIDE][iP].iNODE3; iS3 = sosedi[SSIDE][iP].iNODE3; iB3 = sosedi[BSIDE][iP].iNODE3;

					iE4 = sosedi[ESIDE][iP].iNODE4; iN4 = sosedi[NSIDE][iP].iNODE4; iT4 = sosedi[TSIDE][iP].iNODE4;
					iW4 = sosedi[WSIDE][iP].iNODE4; iS4 = sosedi[SSIDE][iP].iNODE4; iB4 = sosedi[BSIDE][iP].iNODE4;


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
					doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS1 = sosedb[iT1 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iT1]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS2 = sosedb[iT2 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iT2]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS3 = sosedb[iT3 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iT3]) {

								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
							list[ic].dS4 = sosedb[iT4 - maxelm].dS; // Площадь ячейки.
							binc = true;
						}
						else
						{
							if (itek != whot_is_block[iT4]) {
								// вычисление размеров текущего контрольного объёма:
								doublereal dx_loc = 0.0, dy_loc = 0.0, dz_loc = 0.0;// объём текущего контроольного объёма
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
			for (integer j = 0; j < ic; j++) b[itek].radiation.nodelistT[j] = list[j];

			delete[] list;
			list = NULL;


		}
	}
} // constr_link_on_surface_for_radiation_model_alice


  // освобожение памяти уровня 1
  // для уравнений гидродинамики.
void free_level1_flow(FLOW* &fglobal, integer &flow_interior) {
	integer iflow = 0;
	for (integer iflow = 0; iflow<flow_interior; iflow++) {
		integer j = 0;
		if (fglobal != NULL) {

#if doubleintprecision == 1
			printf("delete flow %lld icolor_diferent_fluid_domain\n", iflow);
#else
			printf("delete flow %d icolor_diferent_fluid_domain\n", iflow);
#endif

			if (fglobal[iflow].whot_is_block != NULL) {
				delete[] fglobal[iflow].whot_is_block;
				fglobal[iflow].whot_is_block = NULL;
			}
			
			if (fglobal[iflow].icolor_different_fluid_domain != NULL) {
				delete[] fglobal[iflow].icolor_different_fluid_domain;
				fglobal[iflow].icolor_different_fluid_domain = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld sosedb\n", iflow);
#else
			printf("delete flow %d sosedb\n", iflow);
#endif
			
			if (fglobal[iflow].sosedb != NULL) {
				delete[] fglobal[iflow].sosedb;
				fglobal[iflow].sosedb = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld rdistWall\n", iflow);
#else
			printf("delete flow %d rdistWall\n", iflow);
#endif
			
			if (fglobal[iflow].rdistWall != NULL) { // -1N
				delete[] fglobal[iflow].rdistWall;
				fglobal[iflow].rdistWall = NULL;
			}


#if doubleintprecision == 1
			printf("delete flow %lld sosedi\n", iflow);
#else
			printf("delete flow %d sosedi\n", iflow);
#endif
			
			if (fglobal[iflow].sosedi != NULL) {
				for (integer i = 0; i<12; i++) {
					if (fglobal[iflow].sosedi[i] != NULL) {
						delete[] fglobal[iflow].sosedi[i]; // -12N
					}
				}
			}
			if (fglobal[iflow].sosedi != NULL) {
				delete[] fglobal[iflow].sosedi;
				fglobal[iflow].sosedi = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld nvtx\n", iflow);
#else
			printf("delete flow %d nvtx\n", iflow);
#endif
			
			for (integer i = 0; i<8; i++) { // -8N
#if doubleintprecision == 1
											//printf("%lld\n",i); system("PAUSE"); // debug
#else
											//printf("%d\n",i); system("PAUSE"); // debug
#endif
										
				if (fglobal[iflow].nvtx != NULL) {
					if (fglobal[iflow].nvtx[i] != NULL) {
						delete[] fglobal[iflow].nvtx[i];
					}
				}

			}
			if (fglobal[iflow].nvtx != NULL) {
				delete[] fglobal[iflow].nvtx;
				fglobal[iflow].nvtx = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld prop_b\n", iflow);
#else
			printf("delete flow %d prop_b\n", iflow);
#endif
			
			if (fglobal[iflow].prop_b != NULL) {
				for (integer i = 0; i<3; i++) {
					if (fglobal[iflow].prop_b[i] != NULL) {
						delete[] fglobal[iflow].prop_b[i]; // -3N
					}
				}
			}
			if (fglobal[iflow].prop_b != NULL) {
				delete[] fglobal[iflow].prop_b;
				fglobal[iflow].prop_b = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld prop\n", iflow);
#else
			printf("delete flow %d prop\n", iflow);
#endif
			
			if (fglobal[iflow].prop != NULL) {
				for (integer i = 0; i<3; i++) {
					if (fglobal[iflow].prop[i] != NULL) {
						delete[] fglobal[iflow].prop[i]; // -3N
					}
				}
			}
			if (fglobal[iflow].prop != NULL) {
				delete[] fglobal[iflow].prop;
				fglobal[iflow].prop = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld mf\n", iflow);
#else
			printf("delete flow %d mf\n", iflow);
#endif
			
			if (fglobal[iflow].mf != NULL) {
				for (integer i = 0; i<fglobal[iflow].maxelm; i++) {
					if (fglobal[iflow].mf[i] != NULL) {
						delete[] fglobal[iflow].mf[i]; // -6N
						fglobal[iflow].mf[i] = NULL;
					}
				}
			}
			if (fglobal[iflow].mf != NULL) {
				delete[] fglobal[iflow].mf;
				fglobal[iflow].mf = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld SInvariantStrainRateTensor\n", iflow);
#else
			printf("delete flow %d SInvariantStrainRateTensor\n", iflow);
#endif
			
			if (fglobal[iflow].SInvariantStrainRateTensor != NULL) { // -1N
				delete[] fglobal[iflow].SInvariantStrainRateTensor;
				fglobal[iflow].SInvariantStrainRateTensor = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld diag_coef\n", iflow);
#else
			printf("delete flow %d diag_coef\n", iflow);
#endif
			
			if (fglobal[iflow].diag_coef != NULL) {
				for (integer i = 0; i<3; i++) {
					if (fglobal[iflow].diag_coef[i] != NULL) {
						delete[] fglobal[iflow].diag_coef[i]; // -3N
						fglobal[iflow].diag_coef[i] = NULL;
					}
				}
			}
			if (fglobal[iflow].diag_coef != NULL) {
				delete[] fglobal[iflow].diag_coef;
				fglobal[iflow].diag_coef = NULL;
			}


#if doubleintprecision == 1
			printf("delete flow %lld pa\n", iflow);
#else
			printf("delete flow %d pa\n", iflow);
#endif
			
			if (fglobal[iflow].pa != NULL) { // -3N
				delete[] fglobal[iflow].pa;
				fglobal[iflow].pa = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld ptr\n", iflow);
#else
			printf("delete flow %d ptr\n", iflow);
#endif
			
			if (fglobal[iflow].ptr != NULL) { // -N
				delete[] fglobal[iflow].ptr;
				fglobal[iflow].ptr = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld iN\n", iflow);
#else
			printf("delete flow %d iN\n", iflow);
#endif
			
			if (fglobal[iflow].iN != NULL) {
				for (integer i = 0; i<3; i++) {
					if (fglobal[iflow].iN[i] != NULL) {
						delete fglobal[iflow].iN[i]; // -max3(iWE, iSN, iBT)
						fglobal[iflow].iN[i] = NULL;
					}
				}
			}
			if (fglobal[iflow].iN != NULL) {
				delete fglobal[iflow].iN;
				fglobal[iflow].iN = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld id\n", iflow);
#else
			printf("delete flow %d id\n", iflow);
#endif
			
			if (fglobal[iflow].id != NULL) {
				for (integer i = 0; i<3; i++) {
					if (fglobal[iflow].id[i] != NULL) {
						integer n1 = 0;
						switch (i) {
						case 0: n1 = fglobal[iflow].iWE;  break;
						case 1: n1 = fglobal[iflow].iSN;  break;
						case 2: n1 = fglobal[iflow].iBT;  break;
						}
						for (integer j = 0; j<n1; j++) {
							if (fglobal[iflow].id[i][j] != NULL) {
								delete fglobal[iflow].id[i][j];
								fglobal[iflow].id[i][j] = NULL;
							}
						}
					}
				}
			}
			if (fglobal[iflow].id != NULL) {
				for (integer i = 0; i<3; i++) {
					if (fglobal[iflow].id[i] != NULL) {
						delete fglobal[iflow].id[i]; // 
						fglobal[iflow].id[i] = NULL;
					}
				}
			}
			if (fglobal[iflow].id != NULL) {
				delete fglobal[iflow].id;
				fglobal[iflow].id = NULL;
			}

		}
	}
} // free_level1_flow


  // освобожение памяти уровня 2
  // для уравнений гидродинамики.
void free_level2_flow(FLOW* &fglobal, integer &flow_interior) {
	integer iflow = 0;
	for (integer iflow = 0; iflow<flow_interior; iflow++) {
		integer j = 0;
		if (fglobal != NULL) {

#if doubleintprecision == 1
			printf("delete flow%lld slau\n", iflow);
#else
			printf("delete flow%d slau\n", iflow);
#endif
			
			if (fglobal[iflow].slau != NULL) {
				for (j = 0; j<5; j++) {
					if (fglobal[iflow].slau[j] != NULL) {
						delete[] fglobal[iflow].slau[j];
						fglobal[iflow].slau[j] = NULL;
					}
				}
				delete[] fglobal[iflow].slau;
				fglobal[iflow].slau = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow%lld slau_bon\n", iflow);
#else
			printf("delete flow%d slau_bon\n", iflow);
#endif
			
			if (fglobal[iflow].slau_bon != NULL) {
				for (j = 0; j<5; j++) {
					if (fglobal[iflow].slau_bon[j] != NULL) {
						delete[] fglobal[iflow].slau_bon[j];
						fglobal[iflow].slau_bon[j] = NULL;
					}
				}
				delete[] fglobal[iflow].slau_bon;
				fglobal[iflow].slau_bon = NULL;
			}

#if doubleintprecision == 1
			printf("delete fluid %lld potent\n", iflow);
#else
			printf("delete fluid %d potent\n", iflow);
#endif
		
			if (fglobal[iflow].potent != NULL) { // -26N
				for (integer i = 0; i<33; i++) {
					if (fglobal[iflow].potent[i] != NULL) {
						delete[] fglobal[iflow].potent[i];
						fglobal[iflow].potent[i] = NULL;
					}
				}
				delete[] fglobal[iflow].potent;
				fglobal[iflow].potent = NULL;
			}

#if doubleintprecision == 1
			printf("delete flow %lld alpha\n", iflow);
#else
			printf("delete flow %d alpha\n", iflow);
#endif
			
			if (fglobal[iflow].alpha != NULL) { // -4
				delete[] fglobal[iflow].alpha;
				fglobal[iflow].alpha = NULL;
			}

		}

	}
	//flow_interior=0;
} // free_level2_flow


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
void load_TEMPER_and_FLOW(TEMPER &t, FLOW* &f, integer &inx, integer &iny, integer &inz,
				 doublereal* &xpos, doublereal* &ypos, doublereal* &zpos, integer &flow_interior,
				 BLOCK* b, integer lb, integer lw, WALL* w, SOURCE* s, integer ls, 
	             integer lu, UNION* &my_union, doublereal temp_ref,
				 TPROP* matlist,  bool bextendedprint,
				 doublereal dgx, doublereal dgy, doublereal dgz, bool bALICEflag,
	             bool breconstruct, integer &iunion_id_p1) {

	// eqin - информация о решаемом наборе уравнений.

	integer icount_part = 1; // 19 счётчик хода выполнения

	integer ipolygon_count_limit = 0;

	integer iPnodes_count_shadow = 0, inx_shadow = 0, iny_shadow = 0, inz_shadow = 0;

	//инициализация.
	inx_shadow = inx;
	iny_shadow = iny;
	inz_shadow = inz;
	iPnodes_count_shadow = iPnodes_count_shadow_memo;

	RESTART_CONSTRUCT : // возвращаемся к перестроению всех структур данных:

	ipolygon_count_limit++;

#if doubleintprecision == 1
	printf("inx=%lld, iny=%lld, inz=%lld. Nodes=%lld\n", inx + 1, iny + 1, inz + 1, (inx + 1)*(iny + 1)*(inz + 1));
#else
	printf("inx=%d, iny=%d, inz=%d. Nodes=%d\n", inx + 1, iny + 1, inz + 1, (inx + 1)*(iny + 1)*(inz + 1));
#endif
    
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

	// ТЕПЛОПРОВОДНОСТЬ :
	t.ilevel_alice = NULL;

    integer *evt_t=NULL; // глобальная нумерация контрольных объёмов
	TOCKA_INT* tck_int_list = NULL;

	if (!bALICEflag) {
		enumerate_volume(evt_t, t.maxelm, TEMPERATURE, xpos, ypos, zpos, t.whot_is_block, inx, iny, inz, b, lb, lu, my_union, iunion_id_p1, tck_int_list,w,lw,s,ls);
		t.ilevel_alice = new integer[t.maxelm];
		// При конформной сетке всюду один и тот-же уровень - первый.
		for (integer i_3 = 0; i_3 < t.maxelm; i_3++) {
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
		// Осуществляет уникальную нумерацию листьев принадлежащих модели. Замена constr_sosed.
		// Оставляет пустым sosed=NULL.
		// Эта функция за раз вычисляет части 1-6. part 1..6.
		constr_nodes_nvtx_prop_alice(oc_global, inx, iny, inz, t.maxelm, xpos, ypos, zpos, TEMPERATURE, b, lb, t.whot_is_block, t.pa, t.maxnod, t.nvtx,
			t.prop, t.Sc, t.ipower_time_depend, matlist, t.ilevel_alice);

		if (b_on_adaptive_local_refinement_mesh) {
			// Инвариант корректности АЛИС сетки.
			printf("the invariant correctness 8 november 2016...\n");
			printf("ANES_ALICE_CORRECT in constr_struct_alice.cpp file\n");
			ANES_ALICE_CORRECT(t.maxnod, t.pa, t.maxelm, t.nvtx);
		}
	}
	//printf("part 1\n");
#if doubleintprecision == 1
	printf("part %lld enumerate_volume\n", icount_part++); //22
#else
	printf("part %d enumerate_volume\n", icount_part++); //22
#endif
	 
	integer *ent_t = NULL; // глобальная нумерация узлов.
	if (!bALICEflag) {
		// pa[0..maxnod-1]
		constr_nodes(t.pa, t.maxnod, ent_t, TEMPERATURE, t.whot_is_block, evt_t, inx, iny, inz, xpos, ypos, zpos, b, lb, tck_int_list, t.maxelm);
	}
   
	//printf("part 2\n");    
#if doubleintprecision == 1
	printf("part %lld constr_nodes\n", icount_part++); //22
#else
	printf("part %d constr_nodes\n", icount_part++); //22
#endif

	// находит соседей только среди внутренних КО 
	// для каждого внутреннего контрольного
	// объёма или 0 если соседа нет.
	// sosed[0..11][0..maxelm-1]
	integer **sosed = NULL;
	if (!bALICEflag) {
		constr_sosed(evt_t, ent_t, sosed, t.maxelm, inx, iny, inz, tck_int_list);
	}
	//printf("part 3\n");
#if doubleintprecision == 1
	printf("part %lld constr_neighbours\n", icount_part++); //22
#else
	printf("part %d constr_neighbours\n", icount_part++); //22
#endif

    // для каждого контрольного объёма принадлежащему
    // расчётной области определяет номера его вершин.
	// nvtx[0..7][0..maxelm-1]
	// Нумерация начинается с единицы.
	if (!bALICEflag) {
		// maxelm передается уже посчитанным заранее, его не надо вычислять заново.
		constr_nvtx(evt_t, ent_t, t.nvtx, t.maxelm, inx, iny, inz, tck_int_list);
	}
	//printf("part 4\n");
#if doubleintprecision == 1
	printf("part %lld constr_nvtx\n", icount_part++); //22
#else
	printf("part %d constr_nvtx\n", icount_part++); //22
#endif
	

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
	if (!bALICEflag) {
		constr_prop(evt_t, t.whot_is_block, ent_t, t.prop, t.maxelm, TEMPERATURE, b, lb, inx, iny, inz, t.Sc, t.ipower_time_depend, xpos, ypos, zpos, matlist, tck_int_list);
	}
	//printf("part 6\n");
#if doubleintprecision == 1
	printf("part %lld constr_prop\n", icount_part++); //22
#else
	printf("part %d constr_prop\n", icount_part++); //22
#endif
	

	/*
	*  sosedi, sosedb ( i - internal, b - boundary).
	*/

	// Алгоритм уникальной нумерации граней.
    // В результате каждая граничная грань,
    // принадлежащая расчётной области имеет
    // уникальный номер.
	integer **gran_t=NULL;
	if (!bALICEflag) {
		enumerate_gran_temp(gran_t, t.maxelm, t.nvtx, t.maxbound, sosed, t.pa, s, ls);
	}
	else {
		// calculate_max_bound_temp закодирована 10 сентября 2016. в 10_34.
		// gran_t==NULL не используется.
		calculate_max_bound_temp(oc_global, t.maxbound, (inx + 1)*(iny + 1)*(inz + 1), b, lb, s, ls);
#if doubleintprecision == 1
		printf("maxbound=%lld\n", t.maxbound);
#else
		printf("maxbound=%d\n", t.maxbound);
#endif
		
		// constr_sosedi_alice закодирована 14 сентября 2016. в 12_36.
		// заполнение prop_b. Замена constr_prop_bound. закодирована 14 сентября 2016. в 14_24.
		// nvtxcell остаётся тождественный NULL. Визуализация в техплот на адаптивной локально-измельчённой сетке
		// будет производится как описано в документации ANES решателя.
		// Заполнение t.sosedb, t.binternalsource. Замена constr_sosedb_temp закодирована 20 сентября 2016. в 11_24.
		constr_sosedi_prop_b_alice(oc_global, t.sosedb, t.binternalsource, t.sosedi,
			t.prop, t.prop_b, t.maxbound, t.maxelm, (inx + 1)*(iny + 1)*(inz + 1), b, lb, s, ls, w, lw, t.whot_is_block, t.nvtx, t.pa);
	}
	t.maxp = t.maxelm + t.maxbound;
#if doubleintprecision == 1
	printf("maxelm=%lld maxbound=%lld maxp=%lld\n", t.maxelm, t.maxbound, t.maxp);
#else
	printf("maxelm=%d maxbound=%d maxp=%d\n", t.maxelm, t.maxbound, t.maxp);
#endif
	
	//printf("part 7\n");
#if doubleintprecision == 1
	printf("part %lld enumerate_gran_temp\n", icount_part++); 
#else
	printf("part %d enumerate_gran_temp\n", icount_part++); 
#endif
	

    // Вычисление соседей для каждого внутреннего узла. 
    // Причём множество соседей выбирается и среди  
    // внутренних КО и среди граничных КО.
	if (!bALICEflag) {
		constr_sosedi(t.sosedi, t.maxelm, gran_t, sosed);
	}
	//printf("part 8\n");
#if doubleintprecision == 1
	printf("part %lld constr_sosedi\n", icount_part++); //8
#else
	printf("part %d constr_sosedi\n", icount_part++); //8
#endif
	

	// Заполнение информации о граничных узлах:
	if (!bALICEflag) {
        // Внимание !!! не заполнено поле dS.
		constr_sosedb_temp(t.sosedb, t.whot_is_block, t.binternalsource, t.maxelm, t.maxbound, gran_t, sosed, t.sosedi, t.nvtx, t.pa, b, lb, lw, w, s, ls);
	}
	//printf("part 9\n");
#if doubleintprecision == 1
	printf("part %lld constr_sosedb_temp\n", icount_part++); //9
#else
	printf("part %d constr_sosedb_temp\n", icount_part++); //9
#endif

	// Свойства материала на границе твердотельной области.
	if (!bALICEflag) {
		constr_prop_bound(t.prop, t.prop_b, t.maxelm, t.maxbound, gran_t, sosed, t.nvtx, t.pa, b, lb);
	}
	//printf("part 10\n");
#if doubleintprecision == 1
	printf("part %lld constr_prop_bound\n", icount_part++); //10
#else
	printf("part %d constr_prop_bound\n", icount_part++); //10
#endif

	// Нужно не забыть выделить память под potent.
	//t.alpha=0.9; // параметр нижней релаксации для уравнения теплопроводности
	// На прошлой неделе 12-16 сентября 2016, было выяснено что для стабильности сходимости алгоритма
	// должно стоять значение именно 1.0, а нижняя релаксация обеспечивается с помщью контстанты 
	// bHORF == 0.25 как рекомендовано в Theory Guide ANSYS Fluent. 
	t.alpha = 1.0; // ВНИМАНИЕ !!! только значение t.alpha = 1.0 стабильно. 

	
	// выделение оперативной памяти для задачи теплопроводности.
    // Инициализация также производится. Производится автоматическое определение минимальной температуры на существующих стенках WALL.
	// Данная функция универсальна  подходит также и для AЛИС сетки.
	if (!breconstruct) {
		allocation_memory_temp(t.potent, t.total_deformation, t.slau, t.slau_bon, t.sosedb, t.maxelm, t.maxbound, ls, lw, w, temp_ref);
		//printf("part 11\n");
#if doubleintprecision == 1
		printf("part %lld allocation_memory_temp\n", icount_part++); //11
#else
		printf("part %d allocation_memory_temp\n", icount_part++); //11
#endif
	}


	

	// Создаёт необходимые связи для модели излучения внутри блока.
	if (!breconstruct) {
		if (!bALICEflag) {
			constr_link_on_surface_for_radiation_model(t.maxelm, t.whot_is_block, t.sosedi, t.nvtx, t.pa, b, lb);
		}
		else {
			// sosedb - введено для хранения площади грани граничного узла на АЛИС сетке.
			constr_link_on_surface_for_radiation_model_alice(t.maxelm, t.sosedb, t.whot_is_block, t.sosedi, t.nvtx, t.pa, b, lb);
		}
	}
	// создаёт связи между контрольными объёмами для графической 
    // визуализации.

	//printf("part 11\n");
#if doubleintprecision == 1
	printf("part %lld constr_link_on_surface_for_radiation_model\n", icount_part++); //12
#else
	printf("part %d constr_link_on_surface_for_radiation_model\n", icount_part++); //12
#endif
	
	if (!bALICEflag) {
		if (bextendedprint) {
			printf("extended print. please wait...\n");
		}
		constr_nvtxcell(evt_t, t.sosedb, t.maxbound, t.maxelm, bextendedprint, t.nvtxcell, t.ncell, inx, iny, inz, tck_int_list);
	}

	//if (tck_int_list != NULL) {
		delete[] tck_int_list;
		tck_int_list = NULL;
	//}

#if doubleintprecision == 1
	/*
	for (integer i = 0; i <= inx; i++) {
		printf("xpos[%lld]=%e\n", i, xpos[i]);
	}
	// индексация начинается с нуля и заканчивается значением inx
	for (integer i = 0; i <= iny; i++) {
		printf("ypos[%lld]=%e\n", i, ypos[i]);
	}
	// индексация начинается с нуля и заканчивается значением inx
	for (integer i = 0; i <= inz; i++) {
		printf("zpos[%lld]=%e\n", i, zpos[i]);
	}
	system("PAUSE");
	*/
#else
	/*
	for (integer i = 0; i <= inx; i++) {
		printf("xpos[%d]=%e\n", i, xpos[i]);
	}
	// индексация начинается с нуля и заканчивается значением inx
	for (integer i = 0; i <= iny; i++) {
		printf("ypos[%d]=%e\n", i, ypos[i]);
	}
	// индексация начинается с нуля и заканчивается значением inx
	for (integer i = 0; i <= inz; i++) {
		printf("zpos[%d]=%e\n", i, zpos[i]);
	}
	system("PAUSE");
	*/
#endif
	

	//printf("part 12\n");
#if doubleintprecision == 1
		printf("part %lld constr_nvtxcell\n", icount_part++); //12
#else
		printf("part %d constr_nvtxcell\n", icount_part++); //12
#endif

	// При использовании АЛИС сетки визуализация производится в соответствии с документацией 
	// ANES солвера, поэтому добавление сеточных линий совсем не требуется.
	if( (!bALICEflag)) {

		// 15 августа. В случае нехватки узлов для качественной визуализации на основе 
		// Coarce Mesh сетки произойдёт автоматическое добавление недостающих линий и возврат.
		// Проверка nvtxcell
		integer iPnodes_count = 0;
		bool* found_hash_table = new bool[t.maxelm + 1];
		for (integer i = 1; i <= t.maxelm; i++) {
			found_hash_table[i] = false; // инициализация.
			integer ib = t.whot_is_block[i-1];
			if (b[ib].g.itypegeom == POLYGON) {
				// 30.08.2017. на полигонах мы не используем этот метод.
				found_hash_table[i] = true; // инициализация.
			}
		}

		for (integer i_1 = 0; i_1 < t.ncell; i_1++) {
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
		for (integer i = 1; i <= t.maxelm; i++) {
			if (!found_hash_table[i]) {
				// узел i не найден в nvtxcell
				iPnodes_count++;
			}
		}

#if doubleintprecision == 1
		printf("problem visualisation nodes count iPnodes_count=%lld\n", iPnodes_count);
#else
		printf("problem visualisation nodes count iPnodes_count=%d\n", iPnodes_count);
#endif
		
		bool bpolygon = false;
		for (integer i83 = 1; i83 < lb; i83++) {
			if (b[i83].g.itypegeom == POLYGON) {
				bpolygon = true;
			}
		}

		//system("PAUSE");
		if ((iPnodes_count > 0)&&(ipolygon_count_limit<iGLOBAL_RESTART_LIMIT)/*&&(!bpolygon)*/) {

			for (integer i = 1; i <= t.maxelm; i++) {
				if (!found_hash_table[i]) {
					integer iP = i - 1;
					// iP - номер центрального контрольного объёма
					// iP внутренний КО 0..maxelm-1
					integer iE, iN, iT, iW, iS, iB; // номера соседних контрольных объёмов
					iE = t.sosedi[ESIDE][iP].iNODE1; iN = t.sosedi[NSIDE][iP].iNODE1; iT = t.sosedi[TSIDE][iP].iNODE1; iW = t.sosedi[WSIDE][iP].iNODE1; iS = t.sosedi[SSIDE][iP].iNODE1; iB = t.sosedi[BSIDE][iP].iNODE1;

					// вычисление размеров текущего контрольного объёма:
					//doublereal dx = 0.0, dy = 0.0, dz = 0.0;// объём текущего контроольного объёма
					//volume3D(iP, t.nvtx, t.pa, dx, dy, dz);
					TOCHKA pointP;
					center_cord3D(iP, t.nvtx, t.pa, pointP,100);


					if ((iE >= t.maxelm) && (iW >= t.maxelm)) {
						// бавить в xpos
						addboundary(xpos, inx, pointP.x, YZ, b, lb, w, lw, s, ls);
					}
					if ((iN >= t.maxelm) && (iS >= t.maxelm)) {
						// бавить в ypos
						addboundary(ypos, iny, pointP.y, XZ, b, lb, w, lw, s, ls);
					}
					if ((iT >= t.maxelm) && (iB >= t.maxelm)) {
						// бавить в zpos
						addboundary(zpos, inz, pointP.z, XY, b, lb, w, lw, s, ls);
					}

				}
			}

			// упорядочивание по возрастанию
			Sort_method<doublereal>(xpos,inx);
			Sort_method<doublereal>(ypos,iny);
			Sort_method<doublereal>(zpos,inz);

			// Нужно освободить оперативную память из под всех структур данных:
			free_level1_temp(t);
			if (!breconstruct) {
				free_level2_temp(t); // освобождение памяти из под матриц.
			}

			delete[] found_hash_table;

			// 22.01.2017
			if (sosed != NULL) {
				for (integer i_97 = 0; i_97 < 12; i_97++) {
					if (sosed[i_97] != NULL) {
						delete[] sosed[i_97];
						sosed[i_97] = NULL;
					}
				}
				delete[] sosed;
			}
			sosed = NULL;

			if (gran_t != NULL) {
				for (integer i_97 = 0; i_97 < 6; i_97++) {
					if (gran_t[i_97] != NULL) {
						delete[] gran_t[i_97];
						gran_t[i_97] = NULL;
					}
				}
				delete[] gran_t;
			}
			gran_t = NULL;

			//if (ent_t != NULL) {
				delete[] ent_t;
			//}
			ent_t = NULL;

			//if (evt_t != NULL) {
				delete[] evt_t;
				evt_t = NULL;
			//}

			iPnodes_count_shadow_memo = iPnodes_count;

			printf("progon number is ipolygon_count_limit=%lld\n", ipolygon_count_limit);
			printf("iPnodes_count = %lld=%lld, x= %lld=%lld, y= %lld=%lld, z= %lld=%lld\n", iPnodes_count_shadow, iPnodes_count, inx_shadow,inx, iny_shadow,iny, inz_shadow, inz);
			if (abs(iPnodes_count_shadow - iPnodes_count) + abs(inx_shadow - inx) + abs(iny_shadow - iny) + abs(inz_shadow - inz) == 0) {
				bglobal_restart_06_10_2018_stagnation[ipolygon_count_limit] = true;
			}
			else {
				bglobal_restart_06_10_2018_stagnation[ipolygon_count_limit] = false;
			}
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
			
			iPnodes_count_shadow = iPnodes_count;
			// горячий рестарт рестарт.
			goto RESTART_CONSTRUCT;

		}
		delete[] found_hash_table;
	}


	
   
    // ГИДРОДИНАМИКА:
	integer i=0;

	integer maxelm_global_flow=0; // максимальное суммарное количество внутренних КО принадлежащих расчётной области.
	integer *evt_f=NULL; // глобальная нумерация контрольных объёмов
	integer *whot_is_block_fl=NULL; // побочная инормация.

	//printf("part 13\n");
#if doubleintprecision == 1
	printf("part %lld FLUID DYNAMIC start construct\n", icount_part++); //13
#else
	printf("part %d FLUID DYNAMIC start construct\n", icount_part++); //13
#endif

	

	TOCKA_INT* tck_int_list_flow = NULL;
	
	if (!bALICEflag) {
		enumerate_volume(evt_f, maxelm_global_flow, HYDRODINAMIC, xpos, ypos, zpos, whot_is_block_fl, inx, iny, inz, b, lb, lu, my_union, iunion_id_p1, tck_int_list_flow,w,lw,s,ls);
	}
	else {
		calculate_max_elm(oc_global, maxelm_global_flow, HYDRODINAMIC, b, lb, false);
		printf("maxelm calculate succseful.\n");
		// Это потребуется для алгоритма заливки.
		init_evt_f_alice(evt_f, HYDRODINAMIC, xpos, ypos, zpos, inx, iny, inz, b, lb, tck_int_list_flow,w,lw,s,ls);
	}
	
	//if (whot_is_block_fl != NULL) {
		// информация найдена преждевременно и мы её уничтожим.
		delete[] whot_is_block_fl;
	//}

#if doubleintprecision == 1
	printf("maxelm_global_flow=%lld\n", maxelm_global_flow);
#else
	printf("maxelm_global_flow=%d\n", maxelm_global_flow);
#endif
	
	//printf("part 14\n");
#if doubleintprecision == 1
	printf("part %lld enumerate_volume flow\n", icount_part++); //22
#else
	printf("part %d enumerate_volume flow\n", icount_part++); //22
#endif
	t.rootWE=NULL; t.rootSN=NULL; t.rootBT=NULL;

	
	// 22 сентября 2016 была выделена первая часть, она является
	// универсальной и подходит для АЛИС сетки тоже.
	// Части 11, 12, 12.5, 13.
	// Заполнение evt_f, domain_id, evt_f2.
	integer** evt_f2 = NULL;
	integer* domain_id = NULL;
	constr_ptr_temp_part1(flow_interior, evt_f, evt_f2, domain_id, inx, iny, inz, icount_part);

	

	   // Установка связей между
	   // гидродинамикой и теплообменом.
       // Части   14, 15, 16, 17, 18
	   // После 22 сентября 2016 эта функция несёт нагрузку второй части,
	   // первая часть универсальна и подходит и для АЛИС сетки.
	   if (!bALICEflag) {
		   constr_ptr_temp_part2(flow_interior, evt_f, evt_f2, domain_id, inx, iny, inz, tck_int_list_flow, maxelm_global_flow, icount_part);

		   constr_ptr_temp(flow_interior, f, t.maxelm, t.ptr, evt_t, evt_f, evt_f2, domain_id, inx, iny, inz, breconstruct, icount_part);

		   if (evt_f != NULL) {
			   delete[] evt_f;
		   }
		   evt_f = NULL;
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
	   domain_id = NULL;
		
		// Экспорт в графический визуализатор.
		if (!bALICEflag) {
			printf("please wait...\n");
			// записывает первую и третью части экспортируемого в tecplot360
			// файла на диск. Это делается для того чтобы выгрузить в дальнейшем 
			// 8*t.maxelm*sizeof(doublereal) из оперативной памяти компьютера.
			// Эти же файлы можно будет использовать и для печати гидродинамики.
			// Т.к. печать планируется производить только для внутренних КО.
			if (iunion_id_p1 == 0) {
				exporttecplotxy360T_3D_part1and3(t, t.maxelm, t.maxbound, bextendedprint, t.ncell, t.nvtx, t.nvtxcell, t.pa, t.sosedb, 1, t.ptr);
			}
			else {
				exporttecplotxy360T_3D_part1and3(my_union[iunion_id_p1-1].t, my_union[iunion_id_p1 - 1].t.maxelm, 
					my_union[iunion_id_p1 - 1].t.maxbound, bextendedprint, my_union[iunion_id_p1 - 1].t.ncell,
					my_union[iunion_id_p1 - 1].t.nvtx, my_union[iunion_id_p1 - 1].t.nvtxcell,
					my_union[iunion_id_p1 - 1].t.pa, my_union[iunion_id_p1 - 1].t.sosedb, 1, my_union[iunion_id_p1 - 1].t.ptr);
			}
		}
		if (t.nvtxcell != NULL) {
			for (integer i74 = 0; i74<8; i74++) { // -8N
				if (t.nvtxcell[i74] != NULL) {
					delete[] t.nvtxcell[i74];
					t.nvtxcell[i74] = NULL;
				}
			}
			delete[] t.nvtxcell;
		}
		t.nvtxcell = NULL; // если указатель равен NULL то это предохраняет от повторного удаления памяти.


		if (sosed != NULL) {
			for (i = 0; i < 12; i++) {
				if (sosed[i] != NULL) {
					delete[] sosed[i];
					sosed[i] = NULL;
				}
			}
			delete[] sosed;
		}
		sosed = NULL;

		// освобождение оперативной памяти TEMPER
		if (ent_t != NULL) {
			delete[] ent_t;
		}
		ent_t = NULL;
		//if (evt_t != NULL) {
			//delete evt_t;
		//}
		//evt_t = NULL;
		if (sosed != NULL) {
			for (i = 0; i < 12; i++) {
				if (sosed[i] != NULL) {
					delete[] sosed[i];
					sosed[i] = NULL;
				}
			}
			delete[] sosed;
		}
		sosed = NULL;

		if (gran_t != NULL) {
			for (i = 0; i < 6; i++) {
				if (gran_t[i] != NULL) {
					delete[] gran_t[i];
					gran_t[i] = NULL;
				}
			}
			delete[] gran_t;
		}
		gran_t = NULL;

	if (maxelm_global_flow > 0) {
			// есть зоны FLUID

		// освобождение оперативной памяти TEMPER
		if (ent_t != NULL) {
			delete[] ent_t;
		}
		ent_t=NULL;
		
		if (sosed != NULL) {
			for (i = 0; i < 12; i++) {
				if (sosed[i] != NULL) {
					delete[] sosed[i];
					sosed[i] = NULL;
				}
			}
			delete[] sosed;
		}
		sosed=NULL;

		if (gran_t != NULL) {
			for (i = 0; i < 6; i++) {
				if (gran_t[i] != NULL) {
					delete[] gran_t[i];
					gran_t[i] = NULL;
				}
			}
			delete[] gran_t;
		}
		gran_t=NULL;

		integer j=0; // счётчики
        integer *ent_f=NULL; // глобальная нумерация узлов.
		integer **sosed_f=NULL; // соседи для внутренних КО только среди внутренних КО.
		integer **gran_f=NULL; // уникальная нумерация граней

		

        // Проход по всем жидким зонам:
		for (i=0; i<flow_interior; i++) {
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
					f[i].whot_is_block, t.ptr, t.maxelm);

				// 22 сентября 2016 вычисление функции цвета для АЛИС сетки.
				// Данная функция использует поле  evt_f2[MASKDOMAINFLUIDCOLOR] и только.
				// Вычисление icolor_different_fluid_domain.
				constr_ptr_temp_part2(flow_interior, evt_f, evt_f2, domain_id, inx, iny, inz, tck_int_list_flow, maxelm_global_flow, icount_part);
				if (evt_f != NULL) {
					delete evt_f;
				}
				evt_f = NULL;
				constr_icolor_different_fluid_domain_alice(f[i].maxelm, f[i].icolor_different_fluid_domain,
					inx, iny, inz, evt_f2, xpos, ypos, zpos, f[i].pa, f[i].nvtx, oc_global);
			}

			// integer *ent_f; - глобальная нумерация узлов.
			if (!bALICEflag) {
				constr_nodes_flow(f[i].pa, f[i].maxnod, ent_f, evt_f2, i, inx, iny, inz, xpos, ypos, zpos);
			}

#if doubleintprecision == 1
			printf("part %lld constr_nodes_flow\n", icount_part++); //19
#else
			printf("part %d constr_nodes_flow\n", icount_part++); //19
#endif
	       

            // для каждого контрольного объёма принадлежащему
            // расчётной области определяет номера его вершин.
			if (!bALICEflag) {
				constr_nvtx_flow(evt_f2, f[i].icolor_different_fluid_domain, i, ent_f, f[i].nvtx, f[i].maxelm, inx, iny, inz);
			}
#if doubleintprecision == 1
			printf("part %lld constr_nvtx_flow\n", icount_part++); //20
#else
			printf("part %d constr_nvtx_flow\n", icount_part++); //20
#endif
	        

	        // находит соседей для каждого внутреннего контрольного
            // объёма или 0 если соседа нет.
			if (!bALICEflag) {
				// Это просто пропускаем т.к. информация о sosed_f при сетках АЛИС не требуется, она хранится в octree дереве.
				constr_sosed_flow(evt_f2, i, sosed_f, f[i].maxelm, inx, iny, inz);
			}
#if doubleintprecision == 1
			printf("part %lld constr_neighbour_flow\n", icount_part++); //21
#else
			printf("part %d constr_neighbour_flow\n", icount_part++); //21
#endif
	        

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
#if doubleintprecision == 1
			printf("%lld\n", f[i].maxelm);
#else
			printf("%d\n", f[i].maxelm);
#endif
			 
			if (bwait) {
				//system("PAUSE");
				system("pause");
			}
			if (!bALICEflag) {
				constr_prop_flow(evt_t, t.whot_is_block, evt_f2, i, f[i].prop, f[i].maxelm, b, lb, inx, iny, inz, xpos, ypos, zpos, matlist);
			}
			if (evt_t != NULL) {
				delete[] evt_t;
			}
			evt_t = NULL;
#if doubleintprecision == 1
			printf("part %lld constr_properties_flow\n", icount_part++); //22
#else
			printf("part %d constr_properties_flow\n", icount_part++); //22
#endif
	        

			
	
			/*
	         *  sosedi, sosedb ( i - internal, b - boundary).
	         */

            // Алгоритм уникальной нумерации граней.
            // В результате каждая граничная грань,
            // принадлежащая расчётной области имеет
            // уникальный номер.
			if (!bALICEflag) {
				enumerate_gran_flow(gran_f, f[i].maxelm, f[i].nvtx, f[i].maxbound, sosed_f, f[i].pa);
			}
			else {
				calculate_max_bound_flow(oc_global, f[i].maxbound, (inx + 1)*(iny + 1)*(inz + 1), b, lb, s, ls);
#if doubleintprecision == 1
				printf("maxbound=%lld\n", f[i].maxbound);
#else
				printf("maxbound=%d\n", f[i].maxbound);
#endif
				
			}
	        f[i].maxp = f[i].maxelm + f[i].maxbound;
#if doubleintprecision == 1
			printf("maxelm=%lld maxbound=%lld maxp=%lld\n", f[i].maxelm, f[i].maxbound, f[i].maxp);
			printf("part %lld enumerate_gran_flow\n", icount_part++);
#else
			printf("maxelm=%d maxbound=%d maxp=%d\n", f[i].maxelm, f[i].maxbound, f[i].maxp);
			printf("part %d enumerate_gran_flow\n", icount_part++);
#endif
			 

			


			if (bALICEflag) {
				// constr_sosedi_prop_b_flow_alice закодирована  22 сентября 2016 на основе
				// constr_sosedi_prop_b_alice которая была закодирована 14 сентября 2016. в 12_36.
				// заполнение prop_b. Замена constr_prop_bound_flow. закодирована 14 сентября 2016. в 14_24.
				// nvtxcell остаётся тождественный NULL. Визуализация в техплот на адаптивной локально-измельчённой сетке
				// будет производится как описано в документации ANES решателя.
				// Заполнение f[i].sosedb Замена constr_sosedb_flow закодирована 20 сентября 2016. в 11_24.
				// Работоспособна 24 сентября 2016.
				constr_sosedi_prop_b_flow_alice(oc_global, f[i].sosedb,
					f[i].sosedi, f[i].prop, f[i].prop_b,
					f[i].maxbound, f[i].maxelm, (inx + 1)*(iny + 1)*(inz + 1),
					b, lb, s, ls, w, lw, f[i].whot_is_block, f[i].nvtx, f[i].pa, matlist);

				/*
				// 26.09.2016 Добавок для АЛИС сетки.
				// Убрано 26.03.2019 Как  безсмысленная проверка.
				integer iE2, iN2, iT2, iW2, iS2, iB2; // номера соседних контрольных объёмов
				integer iE3, iN3, iT3, iW3, iS3, iB3; // номера соседних контрольных объёмов
				integer iE4, iN4, iT4, iW4, iS4, iB4; // номера соседних контрольных объёмов
				integer iP = 13; // !!!! внимание почему именно 13. Проверка безсмысленна.

				// -1 если не используется и [0..maxelm+maxbound-1] если используется.

				iE2 = f[i].sosedi[ESIDE][iP].iNODE2; iN2 = f[i].sosedi[NSIDE][iP].iNODE2; iT2 = f[i].sosedi[TSIDE][iP].iNODE2;
				iW2 = f[i].sosedi[WSIDE][iP].iNODE2; iS2 = f[i].sosedi[SSIDE][iP].iNODE2; iB2 = f[i].sosedi[BSIDE][iP].iNODE2;
				iE3 = f[i].sosedi[ESIDE][iP].iNODE3; iN3 = f[i].sosedi[NSIDE][iP].iNODE3; iT3 = f[i].sosedi[TSIDE][iP].iNODE3;
				iW3 = f[i].sosedi[WSIDE][iP].iNODE3; iS3 = f[i].sosedi[SSIDE][iP].iNODE3; iB3 = f[i].sosedi[BSIDE][iP].iNODE3;
				iE4 = f[i].sosedi[ESIDE][iP].iNODE4; iN4 = f[i].sosedi[NSIDE][iP].iNODE4; iT4 = f[i].sosedi[TSIDE][iP].iNODE4;
				iW4 = f[i].sosedi[WSIDE][iP].iNODE4; iS4 = f[i].sosedi[SSIDE][iP].iNODE4; iB4 = f[i].sosedi[BSIDE][iP].iNODE4;

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
				constr_sosedi(f[i].sosedi, f[i].maxelm, gran_f, sosed_f);
			}
#if doubleintprecision == 1
			printf("part %lld constr_neighbour\n", icount_part++); //24
#else
			printf("part %d constr_neighbour\n", icount_part++); //24
#endif
	       

	        // Заполнение информации о граничных узлах:
			if (!bALICEflag) {
				constr_sosedb_flow(f[i].sosedb, t.whot_is_block, f[i].ptr, f[i].maxelm, f[i].maxbound, gran_f, sosed_f, f[i].sosedi, f[i].nvtx, f[i].pa, lw, w, ls, b);
			}
#if doubleintprecision == 1
			printf("part %lld constr_neighbour_db_flow\n", icount_part++); //25
#else
			printf("part %d constr_neighbour_db_flow\n", icount_part++); //25
#endif
			


	        // Свойства жидкого материала на границе жидкой области.
			if (!bALICEflag) {
				constr_prop_bound_flow(f[i].prop, f[i].prop_b, f[i].maxelm, f[i].maxbound, gran_f, sosed_f);
			}
#if doubleintprecision == 1
			printf("part %lld constr_prop_bound_flow\n", icount_part++); //26
#else
			printf("part %d constr_prop_bound_flow\n", icount_part++); //26
#endif
			

			
			// 24 сентября 2016 (данная функция универсальна и подходит также и для АЛИС сетки).
			// Нужно ли фиксировать давление в одной точке.
            // По умолчанию используется последняя граничная точка.
			determination_of_activity_flow(f[i].sosedb, f[i].maxbound, ls, lw, w, f[i].bactive, f[i].bPressureFix, f[i].bLR1free);
#if doubleintprecision == 1
			printf("part %lld determination_of_activity_flow\n", icount_part++); //27
#else
			printf("part %d determination_of_activity_flow\n", icount_part++); //27
#endif
            

			
			


	        // Нужно не забыть выделить память под potent.
	        // Также необходимо установить параметры релаксации.

	        // выделение оперативной памяти для задачи гидродинамики.
			// инициализирует правильным образом вектор potent.
			// 24 сентября 2016 (теперь этот метод работает и на АЛИС сетке).
			if (!breconstruct) {
				// 31.03.2019 Передаём внутрь также значения допусков
				// для функции is_EXTERNAL_FLOW
				allocation_memory_flow( f[i].potent, f[i].slau, f[i].slau_bon,  f[i].sosedb,
					f[i].maxelm, f[i].maxbound, f[i].alpha, ls, lw, w,  dgx, dgy, dgz, f[i].nvtx,
					f[i].pa, f[i].prop, f[i].sosedi, fabs(xpos[1]-xpos[0]), fabs(ypos[1] - ypos[0]), fabs(zpos[1] - zpos[0]),
					fabs(xpos[inx] - xpos[inx-1]), fabs(ypos[iny] - ypos[iny-1]), fabs(zpos[inz] - zpos[inz-1]));
			}
			allocation_memory_flow_2(
				f[i].diag_coef, f[i].maxelm, f[i].maxbound,
				f[i].SInvariantStrainRateTensor,
				f[i].mf);

			{
				// Загрузка распределения начальной скорости.
				errno_t err_inicialization_data = 0;
				FILE* fp_inicialization_data;
#ifdef MINGW_COMPILLER
				fp_inicialization_data=fopen64("load.txt", "r");
#else
				err_inicialization_data = fopen_s(&fp_inicialization_data, "load.txt", "r");
#endif
				if ((f[i].maxelm>0)&&((err_inicialization_data == 0) || ((starting_speed_Vx*starting_speed_Vx + starting_speed_Vy * starting_speed_Vy + starting_speed_Vz * starting_speed_Vz > 1.0e-30) && (steady_or_unsteady_global_determinant != 3)))) {
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

						bool bsimplelinearinterpol = true; // выполняется простая линейная интерполяция скорости на грань.

						integer iflow = 0;
						return_calc_correct_mass_flux_only_interpolation(i32,
							f[i].potent,
							f[i].pa,
							f[i].prop,
							f[i].prop_b,
							f[i].nvtx,
							f[i].sosedi,
							f[i].maxelm,
							f[i].mf[i32],
							f[i].sosedb,
							ls,lw);

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

#if doubleintprecision == 1
			printf("part %lld allocation_memory_flow_2.\n", icount_part++); //28
#else
			printf("part %d allocation_memory_flow_2\n", icount_part++); //28
#endif
			//system("PAUSE");
			 if (0) {
		        xyplot( f, flow_interior, t);
		        printf("after allocate memory. OK.\n");
				system("pause");
	            //system("PAUSE"); // debug 
	         }


			// Освобождение оперативной памяти для каждой зоны FLUID.
			 if (ent_f != NULL) {
				 delete[] ent_f;
			 }
			ent_f=NULL;
			if (sosed_f != NULL) {
				for (j = 0; j < 12; j++) {
					if (sosed_f[j] != NULL) {
						delete[] sosed_f[j];
						sosed_f[j] = NULL;
					}
				}
				delete[] sosed_f;
			}
			sosed_f=NULL;

			if (gran_f != NULL) {
				for (j = 0; j < 6; j++) {
					if (gran_f[j] != NULL) {
						delete[] gran_f[j];
						gran_f[j] = NULL;
					}
				}
				delete[] gran_f;
			}			
			gran_f=NULL;

			// инициализация для полилинейного метода:
			// TODO

			f[i].OpTemp=temp_ref; // Operating Temperature

#if doubleintprecision == 1
			printf("end FLUID %lld\n", i);
#else
			printf("end FLUID %d\n", i);
#endif
			
			//system("PAUSE");

		} // конец цикла по всем FLUID зонам.


		// Вычисление расстояния до стенки не предназначено для 
		// АЛИС сетки т.к. там содержится сборка матрицы не учитывающая  
		// возможного наличия четырёх соседей.
		if (!bALICEflag) {
			printf("share information about the turbulence model\n");
			constr_fluid_equation(f, flow_interior, ls, lw, w);
			printf("end fluid\n");
		}
		else {
			// АЛИС сетка
			if (eqin.fluidinfo[0].iflowregime == 0) {
				// Ламинарный режим
				printf("share information about the turbulence model\n");
				constr_fluid_equation(f, flow_interior, ls, lw, w);
				printf("end fluid\n");
			}
			else {
				// Турбулентный режим.
				printf("share information about the turbulence model\n");
				constr_fluid_equation(f, flow_interior, ls, lw, w); // 13.04.2019
				printf("end fluid\n");
			}
		}
		if (evt_f2 != NULL) {
			for (i = 0; i < 3; i++) {
				if (evt_f2[i] != NULL) {
					delete[] evt_f2[i];
					evt_f2[i] = NULL;
				}
			}
			delete[] evt_f2;
		}
		evt_f2=NULL;

	} else {
       
	   if (!breconstruct) {
		   if (f != NULL) {
			   // Освобождение памяти для структуры flow.
			   free_level1_flow(f, flow_interior);
			   free_level2_flow(f, flow_interior); // освобождение памяти из под матриц.
			   delete[] f;
			   f = NULL;
		   }
		   f = new FLOW[1];
		   f[0].ptr = NULL;
	   }
	   flow_interior = 0; // нет зон FLUID
	   // Связи гидродинамики с теплообменом нет, т.к. 
	   // нет гидродинамической составляющей.
	   if (!breconstruct) {
		   // Если память занята то мы её освободим.
		   if (t.ptr != NULL) {
			   for (i = 0; i<2; i++) {
				   if (t.ptr[i] != NULL) {
					   delete[] t.ptr[i];
					   t.ptr[i] = NULL;
				   }
			   }
			   delete[] t.ptr;
			   t.ptr = NULL;
		   }
		   t.ptr = NULL;
	   }
	   // обнуляющая инициализация ссылок.
	   f[0].nvtx=NULL;
	   f[0].maxnod=0;
	   f[0].maxelm=0;
	   f[0].pa=NULL;
	   f[0].sosedi=NULL;
	   f[0].maxbound=0;
	   f[0].maxp=0;
	   f[0].sosedb=NULL;
	   f[0].icolor_different_fluid_domain = NULL;
	   if (f[0].ptr != NULL) {
		   delete[] f[0].ptr;
	   }
	   f[0].ptr=NULL;
	   f[0].potent=NULL;
	   f[0].prop=NULL;
	   f[0].prop_b=NULL;
	   f[0].alpha=NULL;
	   if (!breconstruct) {
		   f[0].slau = NULL;
		   f[0].slau_bon = NULL;
	   }
	   f[0].mf=NULL;
	   f[0].rdistWall=NULL;
	   f[0].SInvariantStrainRateTensor=NULL;
	   f[0].diag_coef=NULL;
	   f[0].bactive=false;
	   f[0].OpTemp=temp_ref;
	   f[0].bLR1free=true;
	   f[0].bPressureFix=true;
	   if (!breconstruct) {
		   f[0].iN = NULL;
		   f[0].id = NULL;
	   }
	   f[0].iflowregime=LAMINAR;

	   // Вся выделенная память должна быть освобождена.
	   if (evt_f2 != NULL) {
		   for (i = 0; i < 3; i++) {
			   if (evt_f2[i] != NULL) {
				   delete[] evt_f2[i];
				   evt_f2[i] = NULL;
			   }
		   }
		   delete[] evt_f2;
	   }
	   evt_f2 = NULL;
	}

	if (tck_int_list_flow != NULL) {
		delete[] tck_int_list_flow;
		tck_int_list_flow = NULL;
	}

	if (evt_t != NULL) {
		delete[] evt_t;
	}
	evt_t = NULL;
	
	if (maxelm_global_flow == 0) {
        // освобождение оперативной памяти TEMPER
		if (ent_t != NULL) {
			delete[] ent_t;
		}
		ent_t=NULL;
		if (evt_t != NULL) {
			delete[] evt_t;
		}
		evt_t=NULL;
		if (sosed != NULL) {
			for (i = 0; i < 6; i++) {
				if (sosed[i] != NULL) {
					delete[] sosed[i];
					sosed[i] = NULL;
				}
			}
			delete[] sosed;
		}
		sosed=NULL;
		if (gran_t != NULL) {
			for (i = 0; i < 6; i++) {
				if (gran_t[i] != NULL) {
					delete[] gran_t[i];
					gran_t[i] = NULL;
				}
			}
			delete[] gran_t;
		}
		gran_t=NULL;

		// освобождение оперативной памяти FLOW
		if (evt_f != NULL) {
			delete[] evt_f;
		}
		evt_f=NULL;
	}

	if (0) {
		// Старая версия xyplot на АЛИС сетке использовать невозможно и нужно прописать интерполляцию.
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
		for (i=0; i<t.maxnod; i++) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld %e %e %e\n",i+1,t.pa[i].x,t.pa[i].y,t.pa[i].z);
		#else
			fprintf(fp,"%d %e %e %e\n",i+1,t.pa[i].x,t.pa[i].y,t.pa[i].z);
		#endif
			
		}
        fprintf(fp,"\n");
		// вторая секция - контрольные объёмы.
		for (i=0; i<t.maxelm; i++) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld %lld %lld %lld %lld",i+1,t.nvtx[0][i],t.nvtx[1][i],t.nvtx[2][i],t.nvtx[3][i]);
			fprintf(fp," %lld %lld %lld %lld", t.nvtx[4][i],t.nvtx[5][i],t.nvtx[6][i],t.nvtx[7][i]);
			fprintf(fp," %lld %lld %lld", t.sosedi[ESIDE][i].iNODE1, t.sosedi[NSIDE][i].iNODE1, t.sosedi[TSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", t.sosedi[WSIDE][i].iNODE1, t.sosedi[SSIDE][i].iNODE1, t.sosedi[BSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", t.sosedi[EE][i].iNODE1, t.sosedi[NN][i].iNODE1, t.sosedi[TTSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", t.sosedi[WW][i].iNODE1, t.sosedi[SS][i].iNODE1, t.sosedi[BB][i].iNODE1);
		#else
			fprintf(fp,"%d %d %d %d %d",i+1,t.nvtx[0][i],t.nvtx[1][i],t.nvtx[2][i],t.nvtx[3][i]);
			fprintf(fp," %d %d %d %d", t.nvtx[4][i],t.nvtx[5][i],t.nvtx[6][i],t.nvtx[7][i]);
			fprintf(fp," %d %d %d", t.sosedi[ESIDE][i].iNODE1, t.sosedi[NSIDE][i].iNODE1, t.sosedi[TSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", t.sosedi[WSIDE][i].iNODE1, t.sosedi[SSIDE][i].iNODE1, t.sosedi[BSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", t.sosedi[EE][i].iNODE1, t.sosedi[NN][i].iNODE1, t.sosedi[TTSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", t.sosedi[WW][i].iNODE1, t.sosedi[SS][i].iNODE1, t.sosedi[BB][i].iNODE1);
		#endif
			
			fprintf(fp," %e %e %e", t.prop[RHO][i], t.prop[HEAT_CAPACITY][i], t.prop[LAM][i]);
			fprintf(fp," %e %e %e", t.prop[MULT_LAM_X][i], t.prop[MULT_LAM_Y][i], t.prop[MULT_LAM_Z][i]);
		}
        fprintf(fp,"\n");
        // четвёртая секция - связей контрольных объёмов для графической визуализации.
		for (i=0; i<t.ncell; i++) {
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
		for (i=0; i<f.maxnod; i++) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld %e %e %e\n",i+1,f.pa[i].x,f.pa[i].y,f.pa[i].z);
		#else
			fprintf(fp,"%d %e %e %e\n",i+1,f.pa[i].x,f.pa[i].y,f.pa[i].z);
		#endif
			
		}
		fprintf(fp,"\n");
		// вторая секция - контрольные объёмы.
		for (i=0; i<f.maxelm; i++) {
		#if doubleintprecision == 1
			fprintf(fp,"%lld %lld %lld %lld %lld",i+1,f.nvtx[0][i],f.nvtx[1][i],f.nvtx[2][i],f.nvtx[3][i]);
			fprintf(fp," %lld %lld %lld %lld",f.nvtx[4][i],f.nvtx[5][i],f.nvtx[6][i],f.nvtx[7][i]);
			fprintf(fp," %lld %lld %lld", f.sosedi[ESIDE][i].iNODE1, f.sosedi[NSIDE][i].iNODE1, f.sosedi[TSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", f.sosedi[WSIDE][i].iNODE1, f.sosedi[SSIDE][i].iNODE1, f.sosedi[BSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", f.sosedi[EE][i].iNODE1, f.sosedi[NN][i].iNODE1, f.sosedi[TTSIDE][i].iNODE1);
			fprintf(fp," %lld %lld %lld", f.sosedi[WW][i].iNODE1, f.sosedi[SS][i].iNODE1, f.sosedi[BB][i].iNODE1);
		#else
			fprintf(fp,"%d %d %d %d %d",i+1,f.nvtx[0][i],f.nvtx[1][i],f.nvtx[2][i],f.nvtx[3][i]);
			fprintf(fp," %d %d %d %d",f.nvtx[4][i],f.nvtx[5][i],f.nvtx[6][i],f.nvtx[7][i]);
			fprintf(fp," %d %d %d", f.sosedi[ESIDE][i].iNODE1, f.sosedi[NSIDE][i].iNODE1, f.sosedi[TSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", f.sosedi[WSIDE][i].iNODE1, f.sosedi[SSIDE][i].iNODE1, f.sosedi[BSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", f.sosedi[EE][i].iNODE1, f.sosedi[NN][i].iNODE1, f.sosedi[TTSIDE][i].iNODE1);
			fprintf(fp," %d %d %d", f.sosedi[WW][i].iNODE1, f.sosedi[SS][i].iNODE1, f.sosedi[BB][i].iNODE1);
		#endif
			
			fprintf(fp," %e %e %e\n", f.prop[RHO][i], f.prop[MU][i], f.prop[BETA_T][i]);
		}
        fprintf(fp,"\n");
		
        fclose(fp); // закрытие файла
	}
} // meshin_nodes_flow
*/


#endif